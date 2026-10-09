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
    stop: str                      # "gap", "stationary", "max_iter" or "no_progress" (see ascent_direction)


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


def _project_face(N: np.ndarray, r: np.ndarray, fixed: np.ndarray) -> np.ndarray:
    """Orthogonal projection of r onto {d in span N : d_l = 0 for l in fixed} (N orthonormal columns)."""
    if not fixed.any():
        return N @ (N.T @ r)
    E = N[fixed]                                  # rows of N for the fixed coordinates
    _, sv, Vt = np.linalg.svd(E, full_matrices=True)
    rank = int(np.sum(sv > 1e-12 * max(1.0, sv[0] if sv.size else 1.0)))
    Z = Vt[rank:].T                               # null(E N) in the coordinates of N
    NZ = N @ Z
    return NZ @ (NZ.T @ r)


def _feasible_direction(N, r, signs, mu):
    """A feasible unit-or-zero direction from the dual residual r: project r onto the face of the cone that fixes
    every at-bound variable with a positive multiplier or an outward component, adding violators until feasible."""
    fixed = (signs != 0) & (mu > 0)
    for _ in range(len(r) + 1):
        d = _project_face(N, r, fixed)
        bad = (signs != 0) & ~fixed & (signs * d < 0)
        if not bad.any():
            break
        fixed |= bad
    d[fixed] = 0.0
    nd = np.linalg.norm(d)
    return d / nd if nd > 0 else d


_POLISH_RELS = (0.0, 1e-8, 1e-4, -1.0)   # polish faces: current range, two shrunken ranges, whole space
_ROTATE_EVERY = 20                        # iterations between range-rotating refinements (and at every restart)


def ascent_direction(A_upper, A_lower=None, *, null_basis=None, bound_signs=None, tol=1e-10,
                     max_iter=10_000, stat_tol: float = 0.0, f_value: float = 0.0) -> AscentDirection:
    """Steepest-ascent direction for lambda_min(A_upper(d)) [- lambda_max(A_lower(d))] over feasible unit d.

    A_upper, A_lower: splitting-matrix stacks (d, m, m). null_basis: orthonormal basis (d, q) of null(A_eq), None
    for R^d. bound_signs: (d,) in {-1, 0, +1}: +1 at a lower bound, -1 at an upper bound, 0 free. Stops at the first
    of (a) phi_upper - phi <= tol * max(1, ||G||_2); (b) phi_upper <= stat_tol * max(1, |f_value|); (c) max_iter;
    (e) no progress (two consecutive restarts: f stalled at rounding level). stop is "gap" for (a), "stationary" for
    (b) or for (a) with phi <= 0 (d = 0), "max_iter" for (c) and "no_progress" for (e); for (c) and (e) the returned
    point is the best certified one seen (smallest phi_upper - phi), so phi <= phi* <= phi_upper still holds."""
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

    def primal(r, m_full):
        """Feasible d from the residual r and its rate phi. d = 0 is always feasible with rate 0, so a unit d with
        phi < 0 (r is rounding noise near a stationary point) is replaced by d = 0, phi = 0."""
        d = _feasible_direction(N, r, signs, m_full)
        ru = np.linalg.eigvalsh(np.einsum("l,lab->ab", d, Au))
        rl = None if Al is None else np.linalg.eigvalsh(np.einsum("l,lab->ab", d, Al))
        phi = float(ru[0] - (rl[-1] if rl is not None else 0.0))
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
        Tp = P @ np.hstack(col_blocks + [S[:, act]])
        nz = Tp.shape[1]
        C = np.zeros((len(cons), nz))
        off = 0
        for b, c in enumerate(cons):
            C[b, off:off + len(c)] = c
            off += len(c)
        KKT = np.block([[Tp.T @ Tp, C.T], [C, np.zeros((len(cons), len(cons)))]])
        rhs = np.concatenate([np.zeros(nz), np.ones(len(cons))])
        z = np.linalg.lstsq(KKT, rhs, rcond=None)[0][:nz]
        # the minimisers form z + null([Tp; C]); take the one nearest the current iterate, which stays PSD when the
        # iterate lies inside a non-unique optimal face (e.g. phi* = 0 with a family of dual certificates)
        z0 = np.concatenate(z0 + [m[act]])
        Q = np.vstack([Tp, C])
        _, sv, Vt = np.linalg.svd(Q, full_matrices=True)
        Nn = Vt[int(np.sum(sv > 1e-12 * max(1.0, sv[0] if sv.size else 1.0))):].T
        z = z + Nn @ (Nn.T @ (z0 - z))
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

    def rotate(x):
        """Range-rotating refinement (Burer-Monteiro): Y = R R^T / |R|_F^2 per block and mu = s^2, so every point is
        feasible; minimise |T x|^2 over (R, s) by least squares from the current iterate. The face polish keeps the
        range of Y fixed; this moves it, which the first-order iteration does only slowly when the optimal Y is
        rank-deficient (a sharp optimum where several small eigenvalues of A(d) cross)."""
        Yu, Yl, m = unpack(x)
        blocks = [Y for Y in (Yu, Yl) if Y is not None]
        sizes = [Y.shape[0] for Y in blocks]

        def split(v):
            out, off = [], 0
            for k in sizes:
                R = v[off:off + k * k].reshape(k, k)
                out.append(R @ R.T / max(np.sum(R * R), 1e-300))
                off += k * k
            return out, v[off:] ** 2

        def res(v):
            Ys, mm = split(v)
            return T @ pack(Ys[0], Ys[1] if Al is not None else None, mm)

        v0 = []
        for Y in blocks:
            w, V = np.linalg.eigh(Y)
            v0.append((V * np.sqrt(np.maximum(w, 0.0) + 1e-12)).ravel())
        v0 = np.concatenate(v0 + [np.sqrt(np.maximum(m, 0.0))])
        v = least_squares(res, v0, method="trf", xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=200 * len(v0)).x
        Ys, mm = split(v)
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

    x = pack(np.eye(mu_) / mu_, None if Al is None else np.eye(Al.shape[1]) / Al.shape[1], np.zeros(len(B)))
    y, t = x.copy(), 1.0
    fprev = 0.5 * float(np.sum((T @ x) ** 2))
    stop, it, restarted = "max_iter", 0, False
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
            for rnd in range(2):
                if rnd == 1:                  # the fixed-range polish did not close the gap: rotate, then re-polish
                    if best[0] <= tol * gscale or not (it % _ROTATE_EVERY == 0 or restarted):
                        break
                    cand = rotate(x)
                    fc = certify(cand)
                    if fc < fprev:
                        x, y, t, fprev, restarted, stalled = cand, cand.copy(), 1.0, fc, False, False
                for rel in _POLISH_RELS:
                    cand = polish(x, rel)
                    if cand is None:
                        continue
                    fc = certify(cand)
                    if fc < fprev:            # adopt the better polished point and restart the momentum there
                        x, y, t, fprev, restarted, stalled = cand, cand.copy(), 1.0, fc, False, False
            if best[0] <= tol * gscale:
                stop = "gap"
                break
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
    return AscentDirection(d=d, phi=phi, phi_upper=phi_upper, duality_gap=max(phi_upper - phi, 0.0), Y_upper=Yu,
                           Y_lower=Yl, mu=mu_full, rates=ru, rates_lower=rl, n_iter=it, stop=stop)
