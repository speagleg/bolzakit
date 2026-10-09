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

from dataclasses import dataclass
from typing import Callable

import numpy as np
import scipy.linalg as sla
from scipy.optimize import least_squares

from bolzakit.sensitivity import core as _core
from bolzakit.sensitivity.core import DISCRETE_CAVEAT, ritz_upper_bounds, splitting_matrices

__all__ = ["ParametricPencil", "LowestEigenvalue", "Gap", "AscentDirection", "ascent_direction", "IterRecord",
           "EigenOptResult"]


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
