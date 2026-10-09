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
    stop: str                      # "gap" (a), "stationary" (b) or "max_iter" (c)


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


def ascent_direction(A_upper, A_lower=None, *, null_basis=None, bound_signs=None, tol=1e-10,
                     max_iter=10_000, stat_tol: float = 0.0, f_value: float = 0.0) -> AscentDirection:
    """Steepest-ascent direction for lambda_min(A_upper(d)) [- lambda_max(A_lower(d))] over feasible unit d.

    A_upper, A_lower: splitting-matrix stacks (d, m, m). null_basis: orthonormal basis (d, q) of null(A_eq), None
    for R^d. bound_signs: (d,) in {-1, 0, +1}: +1 at a lower bound, -1 at an upper bound, 0 free. Stops at the first
    of (a) phi_upper - phi <= tol * max(1, ||G||_2); (b) phi_upper <= stat_tol * max(1, |f_value|); (c) max_iter."""
    Au = np.asarray(A_upper, dtype=float)
    if Au.ndim != 3 or Au.shape[1] != Au.shape[2]:
        raise ValueError("A_upper must have shape (d, m, m)")
    dim, mu_ = Au.shape[0], Au.shape[1]
    Al = None if A_lower is None else np.asarray(A_lower, dtype=float)
    if Al is not None and (Al.ndim != 3 or Al.shape[0] != dim or Al.shape[1] != Al.shape[2]):
        raise ValueError("A_lower must have shape (d, m', m') with the same d as A_upper")
    N = np.eye(dim) if null_basis is None else np.asarray(null_basis, dtype=float)
    if N.shape[0] != dim:
        raise ValueError("null_basis must have d rows")
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
        d = _feasible_direction(N, r, signs, m_full)
        Ad = np.einsum("l,lab->ab", d, Au)
        ru = np.linalg.eigvalsh(Ad)
        rl = None if Al is None else np.linalg.eigvalsh(np.einsum("l,lab->ab", d, Al))
        phi = float(ru[0] - (rl[-1] if rl is not None else 0.0))
        return d, phi, ru, rl

    def polish(x):
        """Exact solve on the face of the current iterate: Y = W Z W^T on the range of each Y block, the positive
        multipliers free, the others zero; minimise ||T x||^2 subject to tr Z = 1 per block (a KKT system).
        Returns the polished point if it is feasible, else None. This removes the sqrt(eps) accuracy floor of the
        first-order iteration, so the duality-gap test (a) can be met."""
        Yu, Yl, m = unpack(x)
        faces, col_blocks, cons = [], [], []
        for Y, Ab, sgn in ((Yu, Au, 1.0), (Yl, Al, -1.0)):
            if Y is None:
                continue
            w, V = np.linalg.eigh(Y)
            W = V[:, w > 1e-14]
            Fb = _sym_basis(W.shape[1])
            Wb = np.einsum("ia,kab,jb->kij", W, Fb, W)            # basis of {W Z W^T}
            col_blocks.append(sgn * np.einsum("lab,kab->lk", Ab, Wb))
            faces.append((Fb, W))
            cons.append(np.array([np.trace(F) for F in Fb]))
        act = np.nonzero(m > 0)[0]
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
        off, Ys = 0, []
        for Fb, W in faces:
            Z = np.einsum("k,kab->ab", z[off:off + len(Fb)], Fb)
            if np.linalg.eigvalsh(Z)[0] < -1e-12:
                return None
            Ys.append(W @ Z @ W.T)
            off += len(Fb)
        if np.any(z[off:] < 0):
            return None
        mB = np.zeros(len(B))
        mB[act] = z[off:]
        return pack(Ys[0], Ys[1] if Al is not None else None, mB)

    x = pack(np.eye(mu_) / mu_, None if Al is None else np.eye(Al.shape[1]) / Al.shape[1], np.zeros(len(B)))
    y, t, fprev = x.copy(), 1.0, np.inf
    stop, it = "max_iter", 0
    for it in range(1, max_iter + 1):
        x_new = project(y - (T.T @ (T @ y)) / L)
        r = T @ x_new
        fval = 0.5 * float(r @ r)
        if fval > fprev:                      # adaptive restart (O'Donoghue-Candes)
            t, y = 1.0, x.copy()
            continue
        t_new = (1 + np.sqrt(1 + 4 * t * t)) / 2
        y = x_new + ((t - 1) / t_new) * (x_new - x)
        x, t, fprev = x_new, t_new, fval
        phi_upper = float(np.sqrt(2 * fval))
        if phi_upper <= stat_thresh:
            stop = "stationary"
            break
        if it % 10 == 0 or it == 1:
            for cand in (x, polish(x)):
                if cand is None:
                    continue
                rc = T @ cand
                m_full = np.zeros(dim)
                m_full[B] = unpack(cand)[2]
                _, phi, _, _ = primal(rc, m_full)
                if float(np.linalg.norm(rc)) - phi <= tol * gscale:
                    x, stop = cand, "gap"
                    break
            if stop == "gap":
                break
    r = T @ x
    phi_upper = float(np.linalg.norm(r))
    Yu, Yl, mB = unpack(x)
    mu_full = np.zeros(dim)
    mu_full[B] = mB
    if stop == "stationary":
        d = np.zeros(dim)
        phi, ru = 0.0, np.zeros(mu_)
        rl = None if Al is None else np.zeros(Al.shape[1])
    else:
        d, phi, ru, rl = primal(r, mu_full)
    return AscentDirection(d=d, phi=phi, phi_upper=phi_upper, duality_gap=max(phi_upper - phi, 0.0), Y_upper=Yu,
                           Y_lower=Yl, mu=mu_full, rates=ru, rates_lower=rl, n_iter=it, stop=stop)
