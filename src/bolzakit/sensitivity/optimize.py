"""Repeated-eigenvalue-aware steepest ascent for lambda_k(p) and for gaps lambda_{k+1}(p) - lambda_k(p).

At a multiple eigenvalue lambda_k is not differentiable, only directionally: f'(p; d) = lambda_min(A(d)) with
A(d) = sum_l d_l A_l the splitting matrices of the active cluster (Seyranian-Lund-Olhoff 1994; Overton 1988/1992).
The steepest-ascent direction solves max_{d feasible, |d| <= 1} lambda_min(A(d)), whose dual is
    phi* = min_{Y in S, mu >= 0} || P_eq (g(Y) + sum_{l in B} mu_l s_l e_l) ||,   g(Y)_l = <Y, A_l>,
over the spectraplex S = {Y >= 0, tr Y = 1} and bound multipliers mu (spec 5.8). phi* = 0 is the SLO/KKT first-order
condition. With eps-active (near-tied) clusters the model is the reduced model of spec 5.2, so "stationary" means
eps-stationary with eps = mult_rtol * s.

NUMERICAL: computed for the discrete pencil (K, M) as given; discretisation error relative to any continuum operator
is not controlled. Imports only numpy/scipy and bolzakit.sensitivity.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Callable

import numpy as np
import scipy.linalg as sla
from scipy.optimize import least_squares

from bolzakit.sensitivity import core as _core
from bolzakit.sensitivity.core import DISCRETE_CAVEAT, ritz_upper_bounds, splitting_matrices

__all__ = ["ParametricPencil", "LowestEigenvalue", "Gap", "AscentDirection", "ascent_direction", "IterRecord",
           "EigenOptResult", "maximize_eigenvalue"]


# ----------------------------------------------------------------------------------------------- types
@dataclass(frozen=True)
class ParametricPencil:
    K: Callable              # p -> K(p), dense or scipy.sparse
    M: Callable              # p -> M(p)
    dK: Callable             # p -> [dK/dp_l for l in range(d)]
    dM: Callable | None = None   # None: M independent of p


@dataclass(frozen=True)
class LowestEigenvalue:      # maximise lambda_k (k = 0: lambda_min)
    k: int = 0


@dataclass(frozen=True)
class Gap:                   # maximise lambda_{k+1} - lambda_k
    k: int


@dataclass(frozen=True)
class AscentDirection:
    d: np.ndarray                  # unit feasible direction (zeros if phi_upper certifies stationarity)
    phi: float                     # lambda_min(A(d)) (primal); the steepest-ascent rate
    phi_upper: float               # ||P g(Y)|| (dual); phi <= phi* <= phi_upper
    duality_gap: float             # phi_upper - phi
    Y_upper: np.ndarray            # dual matrix for the maximised (upper) cluster's minimum
    Y_lower: np.ndarray | None     # Gap only
    mu: np.ndarray                 # (d,) bound multipliers; 0 for free variables
    rates: np.ndarray              # eigvalsh(A_upper(d))
    rates_lower: np.ndarray | None # eigvalsh(A_lower(d)); Gap only
    n_iter: int                    # FISTA iterations
    stop: str                      # "gap", "stationary", "no_progress" or "max_iter" (contract: ascent_direction)


@dataclass(frozen=True)
class IterRecord:
    it: int
    value: float
    multiplicity: int
    phi: float
    phi_upper: float
    tau: float
    accepted: bool
    screened: bool
    rates: np.ndarray
    refined: bool = False    # an eps-refinement round (tau = 0): mult_rtol was shrunk at this iterate


@dataclass(frozen=True)
class EigenOptResult:
    p: np.ndarray
    value: float
    eigvals: np.ndarray
    multiplicity: int
    phi: float
    phi_upper: float
    Y_upper: np.ndarray
    Y_lower: np.ndarray | None
    mu: np.ndarray
    status: str
    n_iter: int
    n_eigensolves: int
    n_screened: int
    history: tuple
    label: str = "NUMERICAL"
    caveat: str = DISCRETE_CAVEAT
    mult_rtol: float = float("nan")    # the clustering tolerance in force at the end (after any eps-refinement)
    eps: float = float("nan")          # largest eigenvalue gap inside the active clusters at p (see the driver)


# ----------------------------------------------------------------------------------------------- dual subproblem
def _sym_basis(m: int) -> np.ndarray:
    """Frobenius-orthonormal basis of symmetric m x m matrices, shape (m(m+1)/2, m, m)."""
    out = []
    for a in range(m):
        for b in range(a, m):
            E = np.zeros((m, m))
            if a == b:
                E[a, a] = 1.0
            else:
                E[a, b] = E[b, a] = 1 / np.sqrt(2)
            out.append(E)
    return np.array(out)


def _proj_simplex(x: np.ndarray) -> np.ndarray:
    """Euclidean projection onto {x >= 0, sum x = 1}."""
    u = np.sort(x)[::-1]
    css = np.cumsum(u) - 1.0
    r = np.nonzero(u - css / np.arange(1, len(x) + 1) > 0)[0][-1]
    return np.maximum(x - css[r] / (r + 1), 0.0)


def _proj_spectraplex(Y: np.ndarray) -> np.ndarray:
    w, V = np.linalg.eigh((Y + Y.T) / 2)
    return (V * _proj_simplex(w)) @ V.T


# Rank tolerance of the face projection. The rows of N[fixed] have norm <= 1 (N has orthonormal columns), so a pivoted QR
# diagonal |R_kk| <= _FACE_RANK_RTOL * max(1, |R_00|) means the corresponding direction moves the fixed coordinates
# by at most ~1e-12 per unit length. Treating it as free costs at most ~1e-12 of span-N error once _feasible_direction
# zeroes the fixed coordinates, two orders below the 1e-10 certificate tolerance checked by ascent_direction.
_FACE_RANK_RTOL = 1e-12
# _feasible_direction treats a projected residual of norm <= _NOISE_RTOL * |r| as zero: below that it is cancellation
# noise of the projection (the face holds no ascent direction), and normalising it would invent a "unit" direction.
_NOISE_RTOL = 1e3 * np.finfo(float).eps
# ... and rejects a direction whose fixed coordinates, zeroed at the end, exceed _LEAK_RTOL of its norm: zeroing them
# moves d off span N by at most that much, which must stay below the 1e-10 certificate tolerance.
_LEAK_RTOL = 1e-11


def _project_face(N: np.ndarray, r: np.ndarray, fixed: np.ndarray) -> np.ndarray:
    """Orthogonal projection of r onto {d in span N : d_l = 0 for l in fixed} (N orthonormal columns). Exactly zero
    when the fixed coordinates pin all of span N (numerical rank of N[fixed] = q, tolerance _FACE_RANK_RTOL)."""
    if N.shape[0] == N.shape[1]:                 # span N = R^d: just zero the fixed coordinates
        d = r.copy()
        d[fixed] = 0.0
        return d
    if not fixed.any():
        return N @ (N.T @ r)
    E = N[fixed]                                  # rows of N for the fixed coordinates
    # z = N^T r minus its component in range(E^T) (rank-revealing QR), i.e. the projection onto null(E)
    Q, Rq, _ = sla.qr(E.T, mode="economic", pivoting=True)
    diag = np.abs(np.diag(Rq))
    rank = int(np.sum(diag > _FACE_RANK_RTOL * max(1.0, diag[0] if diag.size else 1.0)))
    if rank >= N.shape[1]:                        # null(E) = {0}: the face is {0}; z - Q Q^T z would be pure noise
        return np.zeros(N.shape[0])
    Q = Q[:, :rank]
    z = N.T @ r
    return N @ (z - Q @ (Q.T @ z))


def _feasible_direction(N, r, signs, mu):
    """A feasible unit-or-zero direction from the dual residual r: project r onto the face of the cone that fixes
    every at-bound variable with a positive multiplier or an outward component, adding violators until feasible.
    Returns 0 when the projection is at the noise level of r (_NOISE_RTOL) or when zeroing the fixed coordinates
    would move d off span N by more than _LEAK_RTOL relative (both: no reliable direction on this face)."""
    fixed = (signs != 0) & (mu > 0)
    for _ in range(len(r) + 1):
        d = _project_face(N, r, fixed)
        bad = (signs != 0) & ~fixed & (signs * d < 0)
        if not bad.any():
            break
        fixed |= bad
    nd = np.linalg.norm(d)
    if nd <= _NOISE_RTOL * np.linalg.norm(r) or np.linalg.norm(d[fixed]) > _LEAK_RTOL * nd:
        return np.zeros_like(d)
    d[fixed] = 0.0
    return d / np.linalg.norm(d)


_CERT_SPAN_TOL = 1e-10   # a returned d must lie in span N to this (absolute, d is unit or zero)
_CERT_RTOL = 1e-12       # phi recomputation and weak duality phi <= phi_upper, relative to gscale
_POLISH_RELS = (0.0, 1e-8, 1e-4, -1.0)   # polish faces: current range, two shrunken ranges, whole space
_ROTATE_AFTER = 3     # rotate the range after this many checks (every 10 iterations) without halving the gap
_ROTATE_NFEV = 15     # Gauss-Newton evaluations per rotation


def ascent_direction(A_upper, A_lower=None, *, null_basis=None, bound_signs=None, tol=1e-10,
                     max_iter=10_000, stat_tol: float = 0.0, f_value: float = 0.0) -> AscentDirection:
    """Steepest-ascent direction for lambda_min(A_upper(d)) [- lambda_max(A_lower(d))] over feasible unit d.

    A_upper, A_lower: splitting-matrix stacks (d, m, m). null_basis: orthonormal basis (d, q) of null(A_eq), None
    for R^d. bound_signs: (d,) in {-1, 0, +1}: +1 at a lower bound, -1 at an upper bound, 0 free. Stops at the first
    of (a) phi_upper - phi <= tol * max(1, ||G||_2); (b) phi_upper <= stat_tol * max(1, |f_value|); (c) max_iter;
    (e) no progress (two consecutive restarts: f stalled at rounding level).

    Contract on the returned stop (always phi <= phi* <= phi_upper, d feasible, unit or zero):
    - "gap": phi_upper - phi <= tol * gscale with phi > 0; d is a certified unit ascent direction of rate phi.
    - "stationary": d = 0 and phi = 0, and either phi_upper <= stat_tol * max(1, |f_value|) (rule b) or
      phi_upper <= tol * gscale (rule a reached with phi <= 0): phi* = 0 to that accuracy (a first-order KKT point
      of the reduced model).
    - "no_progress" / "max_iter": the best certified point seen; the gap may exceed tol * gscale. d may be 0 with
      phi = 0 (no feasible direction with a positive certified rate was found); then phi_upper bounds phi*.
    Every returned result passes a final recomputation of these properties (d in span N to 1e-10, sign feasible,
    phi = lambda_min(...) and phi <= phi_upper to 1e-12 * gscale, Y and mu dual feasible); a result that fails is
    replaced by d = 0, phi = 0 with the stop relabelled "stationary" or "no_progress" by the size of phi_upper.
    Callers (the Task 7/9 driver) take a step only for phi > 0, treat "stationary" as convergence, and treat
    "no_progress"/"max_iter" with d = 0 as stationary to accuracy phi_upper.

    Cost: typically 10-30 iterations and a few ms (d = 64, m <= 4, with bounds and equalities: median ~10-15 ms).
    Known tail, for the driver's subproblem budget: small stacks with bounds at a near-degenerate optimum (e.g.
    d = 5, m = 4, several bounds active, phi* ~ 0) can need ~1000 iterations (~0.3 s); a census of 200 such
    instances had median 10 iterations and maximum 1060."""
    Au = np.asarray(A_upper, dtype=float)
    if Au.ndim != 3 or Au.shape[1] != Au.shape[2]:
        raise ValueError("A_upper must have shape (d, m, m)")
    dim, mu_ = Au.shape[0], Au.shape[1]
    Al = None if A_lower is None else np.asarray(A_lower, dtype=float)
    if Al is not None and (Al.ndim != 3 or Al.shape[0] != dim or Al.shape[1] != Al.shape[2]):
        raise ValueError("A_lower must have shape (d, m', m') with the same d as A_upper")
    N = np.eye(dim) if null_basis is None else np.asarray(null_basis, dtype=float)
    if N.ndim != 2 or N.shape[0] != dim:
        raise ValueError("null_basis must have shape (d, q)")
    if np.linalg.norm(N.T @ N - np.eye(N.shape[1])) > 1e-10:
        raise ValueError("null_basis must have orthonormal columns")
    signs = np.zeros(dim) if bound_signs is None else np.asarray(bound_signs, dtype=float)
    if signs.shape != (dim,) or not np.all(np.isin(signs, (-1, 0, 1))):
        raise ValueError("bound_signs must be a (d,) array with entries in {-1, 0, +1}")
    B = np.nonzero(signs)[0]
    P = N @ N.T

    # Linear map T: x = [svec(Y_up), svec(Y_lo), mu_B] -> P (g_up(Y_up) - g_lo(Y_lo) + sum mu_l s_l e_l)
    Eu = _sym_basis(mu_)
    cols = [np.einsum("lab,kab->lk", Au, Eu)]
    nu = Eu.shape[0]
    if Al is not None:
        El = _sym_basis(Al.shape[1])
        cols.append(-np.einsum("lab,kab->lk", Al, El))
        nl = El.shape[0]
    else:
        El, nl = None, 0
    S = np.zeros((dim, len(B)))
    S[B, np.arange(len(B))] = signs[B]
    cols.append(S)
    T = P @ np.hstack(cols)
    L = max(np.linalg.norm(T, 2) ** 2, 1e-300)
    G = np.hstack([Au.reshape(dim, -1)] + ([Al.reshape(dim, -1)] if Al is not None else []))
    gscale = max(1.0, np.linalg.norm(G, 2))
    stat_thresh = stat_tol * max(1.0, abs(f_value))

    def unpack(x):
        Yu = np.einsum("k,kab->ab", x[:nu], Eu)
        Yl = np.einsum("k,kab->ab", x[nu:nu + nl], El) if El is not None else None
        return Yu, Yl, x[nu + nl:]

    def pack(Yu, Yl, m):
        parts = [np.einsum("ab,kab->k", Yu, Eu)]
        if El is not None:
            parts.append(np.einsum("ab,kab->k", Yl, El))
        parts.append(m)
        return np.concatenate(parts)

    def project(x):
        Yu, Yl, m = unpack(x)
        return pack(_proj_spectraplex(Yu), None if Yl is None else _proj_spectraplex(Yl), np.maximum(m, 0.0))

    def primal_ok(d):
        """d is feasible for the cone: zero, or unit, in span N to _CERT_SPAN_TOL and sign feasible on the bounds."""
        nd = np.linalg.norm(d)
        return nd == 0.0 or (abs(nd - 1.0) <= 1e-12 and np.linalg.norm(d - P @ d) <= _CERT_SPAN_TOL
                             and bool(np.all(signs * d >= 0)))

    def rates_of(d):
        ru = np.linalg.eigvalsh(np.einsum("l,lab->ab", d, Au))
        rl = None if Al is None else np.linalg.eigvalsh(np.einsum("l,lab->ab", d, Al))
        return ru, rl, float(ru[0] - (rl[-1] if rl is not None else 0.0))

    def certificate_ok(d, phi, Yu, Yl, mu_full, phi_upper):
        """Recompute the certificate from the fields to be returned: d feasible (primal_ok), phi its rate, Y on the
        spectraplex, mu >= 0 on the bounds only, phi_upper = |P(g(Y) + mu s)|, and weak duality phi <= phi_upper."""
        ctol = _CERT_RTOL * gscale
        if not primal_ok(d) or abs(rates_of(d)[2] - phi) > ctol:
            return False
        for Y in (Yu, Yl):
            if Y is not None and (np.linalg.eigvalsh(Y)[0] < -1e-12 or abs(np.trace(Y) - 1.0) > 1e-12):
                return False
        if np.any(mu_full < 0) or np.any(mu_full[signs == 0] != 0):
            return False
        gv = np.einsum("lab,ab->l", Au, Yu) + mu_full * signs
        if Al is not None:
            gv -= np.einsum("lab,ab->l", Al, Yl)
        return abs(float(np.linalg.norm(P @ gv)) - phi_upper) <= ctol and phi <= phi_upper + ctol

    def primal(r, m_full):
        """Feasible d from the residual r and its rate phi. d = 0 is always feasible with rate 0, so a unit d with
        phi < 0 (r is rounding noise near a stationary point) is replaced by d = 0, phi = 0."""
        d = _feasible_direction(N, r, signs, m_full)
        if not primal_ok(d):
            d = np.zeros(dim)
        ru, rl, phi = rates_of(d)
        if phi < 0.0:
            d, phi, ru = np.zeros(dim), 0.0, np.zeros(mu_)
            rl = None if Al is None else np.zeros(Al.shape[1])
        return d, phi, ru, rl

    def polish(x, rel):
        """Exact solve on a face near the current iterate: Y = W Z W^T with W the eigenvectors of each Y block whose
        eigenvalue exceeds rel * (largest), the multipliers above rel * (largest) free, the others zero; minimise
        ||T x||^2 subject to tr Z = 1 per block (a KKT system). rel < 0 takes the whole space. Returns the polished
        point if it is feasible, else None. This removes the sqrt(eps) accuracy floor of the first-order iteration,
        so the duality-gap test (a) can be met; several rel values let the face grow or shrink."""
        Yu, Yl, m = unpack(x)
        faces, col_blocks, cons, z0 = [], [], [], []
        for Y, Ab, sgn in ((Yu, Au, 1.0), (Yl, Al, -1.0)):
            if Y is None:
                continue
            w, V = np.linalg.eigh(Y)
            W = V[:, w > max(rel * w[-1], 1e-14)] if rel >= 0 else V
            Fb = _sym_basis(W.shape[1])
            Wb = np.einsum("ia,kab,jb->kij", W, Fb, W)            # basis of {W Z W^T}
            col_blocks.append(sgn * np.einsum("lab,kab->lk", Ab, Wb))
            faces.append((Fb, W))
            cons.append(np.array([np.trace(F) for F in Fb]))
            z0.append(np.einsum("ia,ij,jb,kab->k", W, Y, W, Fb))   # current Y in the face coordinates
        act = np.nonzero(m > (rel * m.max() if rel >= 0 and m.size else -1.0))[0]
        TY = P @ np.hstack(col_blocks)
        nY = TY.shape[1]
        C = np.zeros((len(cons), nY))
        off = 0
        for b, c in enumerate(cons):
            C[b, off:off + len(c)] = c
            off += len(c)
        zY0 = np.concatenate(z0)
        for _ in range(len(act) + 1):         # active-set drop: free multipliers that come out negative are removed
            Tp = np.hstack([TY, P @ S[:, act]])
            nz = Tp.shape[1]
            Cp = np.hstack([C, np.zeros((len(cons), len(act)))])
            KKT = np.block([[Tp.T @ Tp, Cp.T], [Cp, np.zeros((len(cons), len(cons)))]])
            # solve for the step from the current iterate zc: the KKT null space is null([Tp; C]) x {0}, so the
            # minimum-norm step gives the minimiser nearest the iterate, which stays PSD when the iterate lies inside
            # a non-unique optimal face (e.g. phi* = 0 with a family of dual certificates)
            zc = np.concatenate([zY0, m[act]])
            rhs = np.concatenate([-(Tp.T @ (Tp @ zc)), 1.0 - Cp @ zc])
            z = zc + np.linalg.lstsq(KKT, rhs, rcond=None)[0][:nz]
            neg = z[nY:] < 0
            if not neg.any():
                break
            act = act[~neg]
        off, Ys = 0, []
        for Fb, W in faces:
            Ys.append(W @ np.einsum("k,kab->ab", z[off:off + len(Fb)], Fb) @ W.T)
            off += len(Fb)
        mB = np.zeros(len(B))
        mB[act] = z[off:]
        xz = pack(Ys[0], Ys[1] if Al is not None else None, mB)
        return project(xz) if feasible(xz) else None   # project: removes rounding-level infeasibility

    def feasible(x):
        Yu, Yl, m = unpack(x)
        return (np.linalg.eigvalsh(Yu)[0] >= -1e-12 and (Yl is None or np.linalg.eigvalsh(Yl)[0] >= -1e-12)
                and bool(np.all(m >= 0)))

    PAu = np.einsum("kl,lab->kab", P, Au)           # projected stacks: residual = sum_b sgn_b <PA_b, Y_b> + P S mu
    PAl = None if Al is None else np.einsum("kl,lab->kab", P, Al)

    def rotate(x):
        """Range-rotating refinement (Burer-Monteiro): Y = R R^T / |R|_F^2 per block and mu_l = s_l^2 for the positive
        multipliers (the others stay 0), so every point is feasible; minimise |T x|^2 over (R, s) by a few
        Gauss-Newton (trust-region) steps from the current iterate, with an analytic Jacobian. The face polish keeps
        the range of Y fixed; this moves it (jointly with the multipliers), which the first-order iteration does only
        slowly when the optimal Y is rank-deficient (a sharp optimum where several small eigenvalues of A(d) cross)."""
        Yu, Yl, m = unpack(x)
        blocks = [(Yu, PAu, 1.0)] + ([] if Al is None else [(Yl, PAl, -1.0)])
        act = np.nonzero(m > 0)[0]
        TB = T[:, nu + nl:]
        TBa = TB[:, act]
        nR = sum(Y.shape[0] ** 2 for Y, _, _ in blocks)

        def Rs(v):
            out, off = [], 0
            for Y, _, _ in blocks:
                k = Y.shape[0]
                out.append(v[off:off + k * k].reshape(k, k))
                off += k * k
            return out

        def res(v):
            out = TBa @ v[nR:] ** 2
            for (_, PA, sg), R in zip(blocks, Rs(v)):
                out += sg * np.einsum("kab,ab->k", PA, R @ R.T) / max(np.sum(R * R), 1e-300)
            return out

        def jac(v):                     # d<PA_l, R R^T / c>/dR = (2/c) (PA_l R - <PA_l, Y> R);  d(s^2)/ds = 2 s
            J = []
            for (_, PA, sg), R in zip(blocks, Rs(v)):
                c = max(np.sum(R * R), 1e-300)
                g = np.einsum("kab,ab->k", PA, R @ R.T) / c
                J.append((sg * 2 / c) * (np.einsum("kab,bc->kac", PA, R) - g[:, None, None] * R).reshape(dim, -1))
            J.append(TBa * (2 * v[nR:]))
            return np.hstack(J)

        v0 = []
        for Y, _, _ in blocks:
            w, V = np.linalg.eigh(Y)
            v0.append((V * np.sqrt(np.maximum(w, 0.0) + 1e-12)).ravel())
        v0.append(np.sqrt(m[act]))
        v = least_squares(res, np.concatenate(v0), jac=jac, method="trf", xtol=1e-15, ftol=1e-15, gtol=1e-15,
                          max_nfev=_ROTATE_NFEV).x
        Ys = [R @ R.T / max(np.sum(R * R), 1e-300) for R in Rs(v)]
        mm = np.zeros_like(m)
        mm[act] = v[nR:] ** 2
        return pack(Ys[0], Ys[1] if Al is not None else None, mm)

    best = None   # (gap, x, d, phi, rates, rates_lower) with the smallest certified gap so far

    def certify(cand):
        nonlocal best
        rc = T @ cand
        m_full = np.zeros(dim)
        m_full[B] = unpack(cand)[2]
        d, phi, ru, rl = primal(rc, m_full)
        m_cut = np.where(m_full > tol * gscale, m_full, 0.0)
        if np.any(m_cut != m_full):           # negligible multipliers must not pin coordinates: try them free too
            alt = primal(rc, m_cut)
            if alt[1] > phi:
                d, phi, ru, rl = alt
        gap = max(float(np.linalg.norm(rc)) - phi, 0.0)
        if best is None or gap < best[0]:
            best = (gap, cand, d, phi, ru, rl)
        return 0.5 * float(rc @ rc)

    x0 = pack(np.eye(mu_) / mu_, None if Al is None else np.eye(Al.shape[1]) / Al.shape[1], np.zeros(len(B)))
    x = x0
    y, t = x.copy(), 1.0
    fprev = 0.5 * float(np.sum((T @ x) ** 2))
    stop, it, restarted, stalled = "max_iter", 0, False, False
    stale, prev_gap, f_at_rot = 0, np.inf, np.inf

    def adopt(cand):
        """Certify cand; if it lowers the dual objective, move the iterate there and restart the momentum."""
        nonlocal x, y, t, fprev, restarted, stalled
        if cand is None:
            return
        fc = certify(cand)
        if fc < fprev:
            x, y, t, fprev, restarted, stalled = cand, cand.copy(), 1.0, fc, False, False

    def polish_all():
        for rel in _POLISH_RELS:
            adopt(polish(x, rel))

    for it in range(1, max_iter + 1):
        x_new = project(y - (T.T @ (T @ y)) / L)
        r = T @ x_new
        fval = 0.5 * float(r @ r)
        if fval > fprev:                      # adaptive restart (O'Donoghue-Candes)
            # a restart right after a restart: even a plain projected-gradient step from x failed to decrease f,
            # so the iteration is at rounding level
            stalled, restarted = restarted, True
            t, y = 1.0, x.copy()
        else:
            stalled, restarted = False, False
            t_new = (1 + np.sqrt(1 + 4 * t * t)) / 2
            y = x_new + ((t - 1) / t_new) * (x_new - x)
            x, t, fprev = x_new, t_new, fval
        if np.sqrt(2 * fprev) <= stat_thresh:
            stop = "stationary"
            break
        if it == 1 or it % 10 == 0 or restarted:
            certify(x)
            polish_all()
            # rotate only when the fixed-range polish has stalled: several checks without halving the certified gap,
            # or a stall of the iteration itself, and only if f moved since the last rotation
            if best[0] > tol * gscale and (stale >= _ROTATE_AFTER or stalled) and fprev < f_at_rot:
                f_at_rot, stale = fprev, 0
                adopt(rotate(x))
                polish_all()
            if best[0] <= tol * gscale:
                stop = "gap"
                break
            stale = stale + 1 if best[0] > 0.5 * prev_gap else 0
            prev_gap = best[0]
        if stalled:
            stop = "no_progress"
            break
    if stop == "stationary":
        r = T @ x
        Yu, Yl, mB = unpack(x)
        d, phi, ru = np.zeros(dim), 0.0, np.zeros(mu_)
        rl = None if Al is None else np.zeros(Al.shape[1])
    else:
        certify(x)
        _, x, d, phi, ru, rl = best
        r = T @ x
        Yu, Yl, mB = unpack(x)
        if stop == "gap" and phi <= 0.0:      # converged with no ascent direction: phi* = 0 to tolerance
            stop = "stationary"
    phi_upper = float(np.linalg.norm(r))
    mu_full = np.zeros(dim)
    mu_full[B] = mB
    # Final soundness guard: never return an uncertified result. If the primal side fails, fall back to d = 0 (always
    # feasible, rate 0) with the same dual point; if the dual side fails too, to the initial dual point Y = I/m, mu = 0.
    if not certificate_ok(d, phi, Yu, Yl, mu_full, phi_upper):
        d, phi = np.zeros(dim), 0.0
        ru, rl = np.zeros(mu_), None if Al is None else np.zeros(Al.shape[1])
        if not certificate_ok(d, phi, Yu, Yl, mu_full, phi_upper):
            Yu, Yl, mB = unpack(x0)
            mu_full = np.zeros(dim)
            phi_upper = float(np.linalg.norm(T @ x0))
        if phi_upper <= max(tol * gscale, stat_thresh):
            stop = "stationary"
        elif stop != "max_iter":
            stop = "no_progress"
    return AscentDirection(d=d, phi=phi, phi_upper=phi_upper, duality_gap=max(phi_upper - phi, 0.0), Y_upper=Yu,
                           Y_lower=Yl, mu=mu_full, rates=ru, rates_lower=rl, n_iter=it, stop=stop)


# ----------------------------------------------------------------------------------------------- driver helpers
def _window(K, M, sparse, k_low: int, k_up: int, mult_rtol: float, validate: bool):
    """Eigenpairs through the top of the single-linkage cluster based at k_up (plus one), and that cluster's top j.
    Returns w, V, j, s. The scale s = max(|lambda_0|, |lambda_{k_low+1}|, 1e-300) is fixed before the growth."""
    n = K.shape[0]
    if validate:
        _core._validate(K, M, sparse)
    count = min(n, k_up + 4)
    while True:
        w, V = _core._eigs(K, M, count, sparse)
        s = max(abs(w[0]), abs(w[k_low + 1]) if k_low + 1 < len(w) else 0.0, 1e-300)
        j = k_up
        while j + 1 < len(w) and w[j + 1] - w[j] <= mult_rtol * s:
            j += 1
        if j + 1 < len(w) or count == n:
            return w, V, j, s
        count = min(n, count + 4)


def _as_list(X, d):
    X = list(X)
    if len(X) != d:
        raise ValueError(f"pencil.dK(p) / pencil.dM(p) must return {d} matrices")
    return X


def _normalised_ascent(A_up, A_lo, f, c, **kw) -> AscentDirection:
    """ascent_direction on the stacks divided by the eigenvalue scale c > 0, with f_value = f / c, mapped back to the
    original units: d and Y are scale-free; phi, phi_upper, the duality gap, mu and the rates scale by c.

    c is the eigenvalue scale s of the pencil (not ||A||: the splitting matrices vanish at a smooth stationary point,
    and dividing by them would blow the residual back up to O(1)). The stacks then hold O(||dK/dp|| / s) data
    whatever the units of K and M, so the subproblem's absolute floors act relatively: in the original units its
    stopping rules read phi_upper <= stat_tol * max(c, |f|) and phi_upper - phi <= tol * max(c, ||A||_2)."""
    ad = ascent_direction(A_up / c, None if A_lo is None else A_lo / c, f_value=f / c, **kw)
    return replace(ad, phi=ad.phi * c, phi_upper=ad.phi_upper * c, duality_gap=ad.duality_gap * c, mu=ad.mu * c,
                   rates=ad.rates * c, rates_lower=None if ad.rates_lower is None else ad.rates_lower * c)


def _stack_norm(A) -> float:
    """||[vec A_1, ..., vec A_d]||_2 >= max over unit d of ||A(d)||_2: a bound on the cluster eigenvalues' speed."""
    return float(np.linalg.norm(A.reshape(A.shape[0], -1), 2)) if A.size else 0.0


def _internal_gap(w, clusters) -> float:
    """Largest consecutive eigenvalue gap inside the given (i, j) clusters (0 if every cluster is a singleton)."""
    g = 0.0
    for i, j in clusters:
        if j > i:
            g = max(g, float(np.max(np.diff(w[i:j + 1]))))
    return g


def maximize_eigenvalue(pencil: ParametricPencil, p0, objective, *, lower=None, upper=None, A_eq=None, b_eq=None,
                        step0: float | None = None, max_iter: int = 100, mult_rtol: float = 1e-3,
                        stat_tol: float = 1e-8, armijo: float = 1e-4, screen: bool = True,
                        refine: bool = True, mult_rtol_min: float = 1e-11, refine_max_iter: int = 100,
                        callback: Callable | None = None, _debug_screen: bool = False,
                        _direction: Callable | None = None, _max_eigensolves: int | None = None,
                        _subproblem_max_iter: int = 2_000) -> EigenOptResult:
    """Maximise lambda_k(p) (LowestEigenvalue(k)) or lambda_{k+1}(p) - lambda_k(p) (Gap(k)) subject to
    lower <= p <= upper and A_eq p = b_eq, by multiple-eigenvalue steepest ascent with an Armijo trust region
    (spec 7). Raises ValueError for an infeasible p0 or an out-of-range objective.k, and TypeError for an unknown
    objective. Statuses:
    - "stationary": phi* = 0 to the subproblem's accuracy (its "stationary" stop, or its rule (b) met by the
      certified bound phi_upper): a stationary point of the reduced model that ties the eigenvalues of each active
      cluster (internal gaps <= result.eps; see eps-refinement below).
      It is the only status that certifies first-order stationarity, and only for that eps-reduced model.
    - "subproblem_not_converged": the subproblem stopped "no_progress" or "max_iter" with no certified ascent
      direction (d = 0, phi = 0); phi_upper bounds phi* but is NOT certified to stat_tol. (A "no_progress" /
      "max_iter" result with phi > 0 still carries a certified ascent direction of rate phi and is stepped along.)
    - "step_too_small": the trust radius fell below 1e-10 * step0 without a strictly improving Armijo step, with no
      refinement left to do (singleton clusters, or internal gaps <= mult_rtol_min * s). At a smooth maximum
      (multiplicity 1) this is the normal end: f(p + tau d) - f = O(tau^2) falls below ulp(f) while phi_upper is
      still above stat_tol * s, so p is optimal to rounding. At multiplicity >= 2 it certifies nothing: the tied
      model (internal gap result.eps) stalled; p is the best iterate found, typically accurate to O(eps).
    - "max_iter": max_iter outer iterations (max_iter + refine_max_iter once eps-refinement starts), or
      _max_eigensolves eigensolves.
    - "lower_coalescence" (LowestEigenvalue(k), k > 0): lambda_k coalesced with lambda_{k-1}.
    - "gap_closed" (Gap): lambda_{k+1} - lambda_k <= mult_rtol * s.
    - "direction_nonpositive": only with the _direction hook, which returned d = 0 or a rate <= 0 (no certificate).

    A trial step is accepted only if it passes the Armijo test AND strictly increases f: once armijo * tau * phi is
    below ulp(f) the Armijo target rounds to f, and equal-value steps would otherwise cycle until max_iter.

    eps-refinement (refine=True): the tied model of a cluster whose eigenvalues are not exactly tied (largest
    internal gap eps > mult_rtol_min * s) is only an eps-model. It is refined when it stalls: at a "stationary" point
    of the eps-model, and when the trust radius has collapsed to where no step can close the internal spread
    (Delta * sigma < eps, sigma = a bound on how fast that spread changes per unit step, or Delta < 1e-10 * step0).
    A round sets mult_rtol <- max(min(eps / (2 s), mult_rtol / 10), mult_rtol_min), which splits the cluster at its
    largest gap, re-evaluates, and is recorded as an IterRecord with refined=True (and passed to the callback); it
    respects _max_eigensolves. Meanwhile a neighbour whose gap could close within the trust radius (gap <= Delta *
    speed, speeds from the eigenvector derivatives) rejoins its cluster, up to the user's mult_rtol * s, so long
    steps do not zigzag across a crossing; a refinement at a stationary point also caps Delta below the split gap.
    Once refinement starts, the outer budget grows to max_iter + refine_max_iter. The result reports the final
    mult_rtol and eps (the largest gap inside the active clusters at p: the reduced model treats eigenvalues that
    far apart as tied; 0 for simple or exactly repeated eigenvalues). Every threshold is relative to s or to
    speeds of the same units, so the iteration is invariant under K -> a K, M -> b M (bit-identical for powers of
    two; otherwise up to the run's own sensitivity to rounding, which for an unconverged or stalled run can be as
    large as its distance from the optimum).
    refine=False keeps the fixed mult_rtol clustering (and no trust-region links).

    An accepted step with poor model agreement (gain < tau * phi / 4) halves the trust radius, which damps the
    fixed-radius zigzag of steepest ascent on an ill-conditioned smooth maximum (it remains first order: such
    instances can still end at "max_iter", close to but short of the maximum).

    The splitting stacks are divided by the eigenvalue scale s before each subproblem (_normalised_ascent), so the
    iteration is invariant under K -> a K, M -> b M, and "stationary" means phi_upper <= stat_tol * max(s, |f|).
    The subproblem budget is capped by _subproblem_max_iter (default 2000, about twice the longest tail in the
    ascent_direction census, i.e. well under a second for small stacks).

    Private keywords (tests and examples only): _debug_screen eigensolves every screened trial and asserts the
    Armijo test would have failed; _direction(sm_upper, null_basis, bound_signs) -> (d, rate) replaces the
    steepest-ascent direction (the plate example's naive baseline); _max_eigensolves caps eigensolves (status
    "max_iter"); _subproblem_max_iter is ascent_direction's max_iter."""
    if not isinstance(objective, (Gap, LowestEigenvalue)):
        raise TypeError("objective must be LowestEigenvalue or Gap")
    p = np.array(p0, dtype=float).ravel()
    dim = p.size
    lo = np.full(dim, -np.inf) if lower is None else np.broadcast_to(np.asarray(lower, float), (dim,)).copy()
    hi = np.full(dim, np.inf) if upper is None else np.broadcast_to(np.asarray(upper, float), (dim,)).copy()
    if np.any(p < lo) or np.any(p > hi):
        raise ValueError("p0 violates the bounds")
    N = None
    if A_eq is not None:
        A_eq = np.atleast_2d(np.asarray(A_eq, float))
        b_eq = np.zeros(A_eq.shape[0]) if b_eq is None else np.atleast_1d(np.asarray(b_eq, float))
        if np.linalg.norm(A_eq @ p - b_eq) > 1e-10 * max(1.0, np.linalg.norm(b_eq)):
            raise ValueError("p0 violates A_eq p = b_eq")
        N = sla.null_space(A_eq)
    is_gap = isinstance(objective, Gap)
    k = int(objective.k)
    K0, M0 = pencil.K(p), pencil.M(p)
    n = K0.shape[0]
    if not 0 <= k < (n - 1 if is_gap else n):
        raise ValueError(f"objective.k = {k} is out of range for an n = {n} pencil "
                         f"({'Gap needs 0 <= k < n - 1' if is_gap else 'LowestEigenvalue needs 0 <= k < n'})")
    span = hi - lo
    if step0 is None:
        step0 = (0.05 * float(np.max(span)) if np.all(np.isfinite(span))
                 else 0.05 * max(1.0, float(np.max(np.abs(p)))))
    D0 = float(step0)
    Delta, Dmax = D0, 100 * D0
    hist: list = []
    counters = {"eig": 0, "screened": 0}
    rtol = float(mult_rtol)          # current clustering tolerance (shrunk by the eps-refinement)

    def evaluate(q, KM=None):
        Kq, Mq = (pencil.K(q), pencil.M(q)) if KM is None else KM
        K, M, sparse = _core._prepare(Kq, Mq)
        counters["eig"] += 1
        w, V, j, s = _window(K, M, sparse, k, k + 1 if is_gap else k, rtol, validate=True)
        return K, M, sparse, w, V, j, s

    def value(w):
        return float(w[k + 1] - w[k]) if is_gap else float(w[k])

    def log(rec):
        hist.append(rec)
        if callback:
            callback(rec)

    def finish(status, st, it, ad):
        w, j = st[3], st[5]
        mult = j - (k + 1 if is_gap else k) + 1
        return EigenOptResult(p=p.copy(), value=value(w), eigvals=w[:j + 2].copy(), multiplicity=mult,
                              phi=ad.phi if ad else float("nan"), phi_upper=ad.phi_upper if ad else float("nan"),
                              Y_upper=ad.Y_upper if ad else np.zeros((0, 0)), Y_lower=ad.Y_lower if ad else None,
                              mu=ad.mu if ad else np.zeros(dim), status=status, n_iter=it,
                              n_eigensolves=counters["eig"], n_screened=counters["screened"], history=tuple(hist),
                              mult_rtol=rtol, eps=eps_now)

    def refine_round(st, it, ad, mult, gabs):
        """One eps-refinement round at the current iterate: shrink rtol below the largest internal gap gabs, extend
        the budget, record it. Returns a finished result if the eigensolve cap forbids the re-evaluation."""
        nonlocal rtol, limit
        if _max_eigensolves is not None and counters["eig"] >= _max_eigensolves:
            return finish("max_iter", st, it, ad)
        limit = max_iter + refine_max_iter
        rtol = max(min(0.5 * gabs / st[6], 0.1 * rtol), mult_rtol_min)
        log(IterRecord(it, value(st[3]), mult, ad.phi, ad.phi_upper, 0.0, False, False, ad.rates, refined=True))
        return None

    st = evaluate(p, (K0, M0))
    ad = None
    eps_now = float("nan")
    limit = max_iter
    it = 0
    while it < limit:
        K, M, sparse, w, V, j, s = st
        f = value(w)
        # active clusters and status conditions
        if is_gap:
            if w[k + 1] - w[k] <= rtol * s:
                return finish("gap_closed", st, it, ad)
            i = k
            while i > 0 and w[i] - w[i - 1] <= rtol * s:
                i -= 1
            up, low = (k + 1, j), (i, k)
        else:
            if k > 0 and w[k] - w[k - 1] <= rtol * s:
                return finish("lower_coalescence", st, it, ad)
            up, low = (k, j), None
        dKs = _as_list(pencil.dK(p), dim)
        dMs = None if pencil.dM is None else _as_list(pencil.dM(p), dim)
        sm_up = splitting_matrices(K, M, dKs, dMs, up, eig=(w, V), check_gap=False, validate=False)
        sm_lo = (splitting_matrices(K, M, dKs, dMs, low, eig=(w, V), check_gap=False, validate=False)
                 if low else None)
        # Trust-region links (eps-refinement only): a neighbour whose gap to a cluster could close within the trust
        # radius joins it (gap <= Delta * (its speed + the cluster's speed)), but never beyond the user's mult_rtol * s.
        # After the refinement has shrunk rtol, this re-forms the cluster while steps are long enough to cross it;
        # otherwise single-eigenvalue steps zigzag across the crossing. Speeds bound |d lambda / d tau|, unit d.
        links = []                                   # (gap, speed) of every trust-region link added
        if refine:
            def speed(a):
                v = V[:, a]
                g = np.array([float(v @ np.asarray(dKs[l] @ v).ravel()) - (0.0 if dMs is None else
                              w[a] * float(v @ np.asarray(dMs[l] @ v).ravel())) for l in range(dim)])
                return float(np.linalg.norm(g if N is None else N.T @ g))
            n = K.shape[0]
            rho = _stack_norm(sm_up.mats)
            j0 = j
            while j + 1 < n:
                if j + 2 > len(w) and len(w) < n:      # splitting needs eigenpair j + 2 too
                    w, V = _core._eigs(K, M, min(n, len(w) + 4), sparse)
                    counters["eig"] += 1
                    st = (K, M, sparse, w, V, j, s)
                sp_ = speed(j + 1) + rho
                if w[j + 1] - w[j] > min(Delta * sp_, mult_rtol * s):
                    break
                links.append((w[j + 1] - w[j], sp_))
                j += 1
            if j > j0:                                 # the loop kept len(w) >= min(j + 2, n)
                st = (K, M, sparse, w, V, j, s)
                up = (up[0], j)
                sm_up = splitting_matrices(K, M, dKs, dMs, up, eig=(w, V), check_gap=False, validate=False)
            if low:
                i0, rho = low[0], _stack_norm(sm_lo.mats)
                i = i0
                while i > 0:
                    sp_ = speed(i - 1) + rho
                    if w[i] - w[i - 1] > min(Delta * sp_, mult_rtol * s):
                        break
                    links.append((w[i] - w[i - 1], sp_))
                    i -= 1
                if i < i0:
                    low = (i, k)
                    sm_lo = splitting_matrices(K, M, dKs, dMs, low, eig=(w, V), check_gap=False, validate=False)
        eps_now = _internal_gap(w, [up] + ([low] if low else []))
        # speed bound for the spread inside the active clusters: two eigenvalues of a cluster separate or close at
        # most at 2 ||A(d)|| <= 2 ||stack||_2 per unit step, a linked neighbour at its link speed (scale-free)
        sigma = max([2 * _stack_norm(sm_up.mats)] + ([2 * _stack_norm(sm_lo.mats)] if low else [])
                    + [sp_ for _, sp_ in links])
        # at-bound test: within 1e-12 (relative) or within a negligible 1e-6 fraction of the trust radius, so a
        # variable a rounding distance from its bound never limits the step to ~1e-12
        tolb = np.maximum(1e-12 * np.where(np.isfinite(span), np.maximum(1.0, span), 1.0), 1e-6 * Delta)
        signs = np.where(p - lo <= tolb, 1.0, np.where(hi - p <= tolb, -1.0, 0.0))
        if _direction is None:
            c = s if s > 1e-300 else 1.0
            ad = _normalised_ascent(sm_up.mats, None if sm_lo is None else sm_lo.mats, f, c, null_basis=N,
                                    bound_signs=signs, max_iter=_subproblem_max_iter, stat_tol=stat_tol)
            # Stationary: the subproblem's "stationary" label, or its rule (b) met by the certified dual bound
            # (phi* <= phi_upper <= stat_tol * max(1, |f| / c) in the normalised units). The latter covers a "gap"
            # stop whose tiny certified phi > 0 the subproblem reported first (it checks rule (b) on its FISTA
            # iterate, not on the polished point); without it the driver zigzags at rate ~1e-9.
            if ad.stop == "stationary" or ad.phi_upper <= stat_tol * max(c, abs(f)):
                if not refine or eps_now <= mult_rtol_min * s:
                    return finish("stationary", st, it, ad)
                # eps-refinement: stationary only for the eps-reduced model; split the cluster and go on (shrink the
                # trust radius too, so that no trust-region link re-forms across the largest gap)
                if links:
                    Delta = min(Delta, 0.5 * eps_now / max(sp_ for _, sp_ in links))
                done = refine_round(st, it, ad, up[1] - up[0] + 1, eps_now)
                if done is not None:
                    return done
                st = evaluate(p, (K, M))
                it += 1
                continue
            # "no_progress"/"max_iter" with d = 0 certify nothing beyond phi* <= phi_upper: not converged. With
            # phi > 0 they still carry a certified ascent direction of rate phi, which the step below uses.
            if ad.phi <= 0:
                return finish("subproblem_not_converged", st, it, ad)
            d, phi = ad.d, ad.phi
        else:
            d, phi = _direction(sm_up, N, signs)
            d = np.asarray(d, dtype=float)
            phi = float(phi)
            ad = AscentDirection(d=d, phi=phi, phi_upper=phi, duality_gap=0.0, Y_upper=np.zeros((0, 0)),
                                 Y_lower=None, mu=np.zeros(dim), rates=np.linalg.eigvalsh(
                                     np.einsum("l,lab->ab", d, sm_up.mats)), rates_lower=None, n_iter=0, stop="gap")
            if not np.any(d) or phi <= 0:
                return finish("direction_nonpositive", st, it, ad)
        # maximal feasible step along d
        with np.errstate(divide="ignore", invalid="ignore"):
            steps = np.where(d > 1e-15, (hi - p) / d, np.where(d < -1e-15, (lo - p) / d, np.inf))
        tau_b = float(np.min(steps)) if steps.size else np.inf
        mult = up[1] - up[0] + 1
        accepted = refined = False
        while not accepted:
            # eps-limited stall: no step within the radius can close the spread inside the clusters, so the tied
            # model is invalid at this scale (it stalls at a small phi > stat_tol): refine instead of shrinking on
            if refine and eps_now > mult_rtol_min * s and (Delta * sigma < eps_now or Delta < 1e-10 * D0):
                done = refine_round(st, it, ad, mult, eps_now)
                if done is not None:
                    return done
                refined = True
                break
            if Delta < 1e-10 * D0:
                return finish("step_too_small", st, it, ad)
            if _max_eigensolves is not None and counters["eig"] >= _max_eigensolves:
                return finish("max_iter", st, it, ad)
            tau = min(Delta, tau_b)
            q = np.clip(p + tau * d, lo, hi)
            target = f + armijo * tau * phi
            Kq = Mq = None
            if screen and not is_gap:
                Kq, Mq = pencil.K(q), pencil.M(q)
                theta = ritz_upper_bounds(Kq, Mq, V[:, :j + 1])[k]
                # rounding margin (scale-free: relative to the eigenvalue scale s); lambda_k <= theta_k
                if theta + 1e-12 * max(abs(f), s) < target:
                    counters["screened"] += 1
                    if _debug_screen:
                        Kc, Mc, sq = _core._prepare(Kq, Mq)
                        wq = _core._eigs(Kc, Mc, k + 1, sq)[0]
                        assert wq[k] < target, "screen rejected a step that passes the Armijo test"
                    log(IterRecord(it, f, mult, phi, ad.phi_upper, tau, False, True, ad.rates))
                    Delta = tau / 2
                    continue
            st_q = evaluate(q, None if Kq is None else (Kq, Mq))
            fq = value(st_q[3])
            accepted = fq >= target and fq > f       # strict: a rounding-level target equals f (see docstring)
            log(IterRecord(it, f, mult, phi, ad.phi_upper, tau, accepted, False, ad.rates))
            if accepted:
                # spec 7 step 8; a bound-limited step (tau < Delta) never shrinks Delta
                if fq - f >= 0.75 * tau * phi:
                    Delta = min(max(2 * tau, Delta), Dmax)
                elif fq - f < 0.25 * tau * phi:          # poor agreement: shrink (avoids a zigzag at fixed tau)
                    Delta = min(Delta, tau / 2)
                p, st = q, st_q
            else:
                Delta = tau / 2
        if refined:
            st = evaluate(p, (K, M))
        it += 1
    return finish("max_iter", st, it, ad)
