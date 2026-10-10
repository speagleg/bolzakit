"""Frequency trimming of a split degenerate mode pair (an inverse problem; NUMERICAL, discrete problem).

A nominally degenerate pair (i, i + 1) of K x = lambda M x is split by imperfections. Candidate trims l = 0..L-1 (a
point mass at a site, removed material, a stiffness patch) change the pencil linearly in their amounts t_l:
    K(t) = K + sum_l t_l dK_l,    M(t) = M + sum_l t_l dM_l.
In the pair's M-orthonormal eigenbasis U the reduced 2x2 model is D + A(dt), D = diag(-g/2, g/2) (g the split),
A(dt) = sum_l dt_l A_l with the splitting matrices A_l = U^T (dK_l - lam_bar dM_l) U. The pair is tied to first order
iff (D + A)_00 = (D + A)_11 and (D + A)_01 = 0, two linear equations G dt = (g, 0) with
G = [A_00 - A_11; A_01]. The first step is the linear program
    min sum_l c_l |t_l + dt_l|  subject to  G dt = (g, 0),  lower <= t + dt <= upper
(with the default objective="sum" a vertex uses at most two trims). Later steps are Newton corrections on the trims the LP selected (least squares on
that support, exact re-solve each time); the LP is re-solved only if a correction leaves the bounds or flips the sign
of a trim. Convergence is quadratic until the split reaches the rounding floor of the eigensolver.

Objective: the trim magnitude sum_l c_l |t_l| (c >= 0), so removal (lower < 0, e.g. ablation) costs like addition. The
LP is posed over t = p - q with p, q >= 0 and p in [max(lo,0), max(hi,0)], q in [max(-hi,0), max(-lo,0)]; for
lower >= 0 it is exactly min c^T t. The objective is bounded below by 0, so the LP is never unbounded (negative costs
are rejected with a ValueError). With objective="sum" a vertex still has at most two active sites; objective="peak"
minimises max_l c_l |t_l| instead (an epigraph LP) and spreads the trim over equivalent sites.

Backtracking: a trial iterate whose split is not below the best so far is rejected and the step halved, up to
`backtrack` times (TrimRecord.step "halve"). Each trial costs one eigensolve and counts towards n_eigensolves, n_iter
and max_iter. A ClusterSplitError raised at a trial iterate propagates (it is not treated as a rejected step).

Resolution floor: TrimResult.rel_floor = eps * lambda_max / |lam_bar|, with lambda_max the largest eigenvalue of the
pencil (K, M) at t = t0 (a dense eigh of the top eigenvalue for small n, otherwise a loose-tolerance ARPACK
``eigsh(which="LA")`` run, a few Lanczos steps; if ARPACK fails, a dense eigh for n <= 4000, else the upper bound
||K||_inf / gamma(M) with gamma(M) > 0 the Gershgorin lower bound on lambda_min(M), else nan) and lam_bar the best
pair's mean. It is a conservative estimate: an upper estimate of what the eigensolver can resolve, which can be loose
(on the spoked ring of examples/gyro_trim.py splits fall well below it). The API returns the raw split_after /
rel_split_after and keeps status "converged" whenever the split met tol, even below rel_floor; a split below
rel_floor is not resolved by the eigensolver, so reports must quote it as "<= rel_floor", never as a measured value.

Statuses: "converged" (split <= tol * |lam_bar|), "stalled" (two consecutive iterations without halving the best
split; or, when the iteration budget ran out with tol below the floor (rel_floor >= tol), the best split is within a
factor 2 of rel_floor: the resolution floor), "max_iter" (budget exhausted otherwise), "infeasible" (the LP has no
solution: no trim pattern within the bounds ties the pair to first order; t is the best iterate found),
"lp_failed" (HiGHS stopped without a solution for another reason, linprog status 1 iteration limit or 4 numerical
difficulties; t is the best iterate found). The pair must stay isolated (the core gap rule; ClusterSplitError
otherwise, with a message that names the pair's own spread when the pair is not nearly degenerate).

NUMERICAL: computed for the discrete pencil (K, M) as given; discretisation error relative to any continuum operator
is not controlled. Imports only numpy/scipy and bolzakit.sensitivity.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
import scipy.linalg as sla
import scipy.sparse.linalg as spla
from scipy.optimize import linprog

from bolzakit.sensitivity.core import DISCRETE_CAVEAT, smallest_eigenpairs, splitting_matrices
from bolzakit.sensitivity.errors import ClusterSplitError

__all__ = ["TrimRecord", "TrimResult", "trim_degenerate_pair"]


@dataclass(frozen=True)
class TrimRecord:
    it: int
    t_cost: float          # sum c |t| at this iterate
    split: float           # lambda_{i+1} - lambda_i
    rel_split: float       # split / |lam_bar|
    lam_bar: float
    step: str              # how the NEXT iterate was obtained: "lp", "newton", "halve" (rejected trial; step halved) or "" (last record)


@dataclass(frozen=True)
class TrimResult:
    t: np.ndarray                  # trim amounts of the best iterate (lowest split)
    cost: float                    # sum_l c_l |t_l|
    split_before: float
    split_after: float
    rel_split_before: float
    rel_split_after: float
    eigvals_before: np.ndarray     # the pair at t = t0
    eigvals_after: np.ndarray      # the pair at t
    rel_floor: float               # eigensolver resolution estimate eps * lambda_max / |lam_bar| (see module docstring)
    status: str                    # "converged", "stalled", "max_iter", "infeasible" or "lp_failed"
    n_iter: int
    n_eigensolves: int
    history: tuple
    label: str = "NUMERICAL"
    caveat: str = DISCRETE_CAVEAT


def _bounds(L, lower, upper, cost):
    """Broadcast and validate (lower, upper, cost) for L trims."""
    for name, v in (("lower", lower), ("upper", upper), ("cost", cost)):
        if v is not None and np.ndim(v) > 0 and np.shape(v) != (L,):
            raise ValueError(f"{name} must be a scalar or have length {L} (one entry per trim), got shape "
                             f"{np.shape(v)}")
    lo = np.broadcast_to(np.asarray(lower, float), (L,)).copy()
    hi = np.full(L, np.inf) if upper is None else np.broadcast_to(np.asarray(upper, float), (L,)).copy()
    if np.any(lo > hi):
        raise ValueError("lower must not exceed upper")
    c = np.ones(L) if cost is None else np.asarray(cost, float)
    c = np.broadcast_to(c, (L,)).copy()
    if np.any(c < 0):
        raise ValueError("cost must be non-negative (the objective is the trim magnitude sum c |t|)")
    return lo, hi, c


def _lp(G, b, lo, hi, c, objective="sum"):
    """min sum c|x| ("sum") or min max c|x| ("peak") subject to G x = b, lo <= x <= hi, posed over x = p - q with
    p, q >= 0. Returns (x, "") on success, (None, "infeasible") for linprog status 2, (None, "lp_failed") otherwise.
    Status 3 (unbounded) cannot occur: c >= 0 bounds the objective below by 0."""
    L = len(c)
    bounds = list(zip(np.maximum(lo, 0.0), np.maximum(hi, 0.0))) + list(zip(np.maximum(-hi, 0.0),
                                                                             np.maximum(-lo, 0.0)))
    if objective == "sum":
        res = linprog(np.concatenate([c, c]), A_eq=np.hstack([G, -G]), b_eq=b, bounds=bounds, method="highs")
    else:                              # z >= c_l (p_l + q_l); a small sum term makes p_l q_l = 0 at the optimum
        A_ub = np.hstack([np.diag(c), np.diag(c), -np.ones((L, 1))])
        res = linprog(np.concatenate([1e-6 / L * c, 1e-6 / L * c, [1.0]]), A_ub=A_ub, b_ub=np.zeros(L),
                      A_eq=np.hstack([G, -G, np.zeros((G.shape[0], 1))]), b_eq=b, bounds=bounds + [(0, None)],
                      method="highs")
    if res.status == 0:
        return np.asarray(res.x[:L] - res.x[L:2 * L], dtype=float), ""
    return None, ("infeasible" if res.status == 2 else "lp_failed")


def _check_objective(objective):
    if objective not in ("sum", "peak"):
        raise ValueError('objective must be "sum" or "peak"')


def _zero_like(K):
    return sp.csc_matrix(K.shape) if sp.issparse(K) else np.zeros(K.shape)


def _lambda_max(K, M):
    """Largest eigenvalue of K x = lambda M x (cheap estimate; used only for the resolution floor)."""
    n = K.shape[0]
    if n <= 400 or not (sp.issparse(K) or sp.issparse(M)):
        return _lambda_max_dense(K, M, n)
    try:
        return float(spla.eigsh(sp.csc_matrix(K), k=1, M=sp.csc_matrix(M), which="LA", tol=1e-3,
                                return_eigenvectors=False)[0])
    except Exception:                  # ARPACK did not converge: never fall back on a lower bound (floor too small)
        pass
    if n <= 4000:
        return _lambda_max_dense(K, M, n)
    Kc, Mc = sp.csr_matrix(K), sp.csr_matrix(M)
    k_norm = float(abs(Kc).sum(axis=1).max())                    # lambda_max(K) <= ||K||_inf
    d = Mc.diagonal()
    gamma = float(np.min(d - (np.asarray(abs(Mc).sum(axis=1)).ravel() - np.abs(d))))   # Gershgorin: <= lambda_min(M)
    return k_norm / gamma if gamma > 0 else float("nan")


def _lambda_max_dense(K, M, n):
    Kd = K.toarray() if sp.issparse(K) else np.asarray(K)
    Md = M.toarray() if sp.issparse(M) else np.asarray(M)
    return float(sla.eigh(Kd, Md, eigvals_only=True, subset_by_index=[n - 1, n - 1])[0])


def _not_nearly_degenerate(w, i, j, gap_ratio_min, err):
    """Re-raise the core gap failure with the right cause when the pair's own split, not a close third eigenvalue, is
    what fails the rule (g_out >= gap_ratio_min |lam_bar| but g_out < 10 g_in)."""
    lb, g_in = float(w[i:j + 1].mean()), float(w[j] - w[i])
    g_out = min(float(w[i] - w[i - 1]) if i > 0 else np.inf, float(w[j + 1] - w[j]) if j + 1 < len(w) else np.inf)
    if g_out >= gap_ratio_min * abs(lb):
        ratio = g_out / max(gap_ratio_min * abs(lb), 10 * g_in)
        raise ClusterSplitError(
            f"pair {(i, j)} is not nearly degenerate: its split {g_in:.3g} is more than a tenth of the distance "
            f"{g_out:.3g} to the nearest other eigenvalue (gap ratio {ratio:.3g} < 1); trim_degenerate_pair needs an "
            "isolated, nearly degenerate pair") from err


def _strictly_inside(t, lo):
    """t_l > lo_l by more than a rounding margin (always true for lo_l = -inf)."""
    lo_f = np.where(np.isfinite(lo), lo, 0.0)
    return ~np.isfinite(lo) | (t > lo_f + 1e-14 * np.maximum(1.0, np.abs(lo_f)))


def trim_degenerate_pair(K, M, dKs, dMs, pair, *, lower=0.0, upper=None, cost=None, tol: float = 1e-10,
                         max_iter: int = 20, gap_ratio_min: float = 1e-3, objective: str = "sum",
                         backtrack: int = 4) -> TrimResult:
    """Smallest-cost trim t (lower <= t <= upper) that ties the split pair (i, i + 1) of K x = lambda M x.

    ``dKs`` / ``dMs``: length-L sequences of the per-unit changes of K and M for each candidate trim (dense or
    sparse); None means that matrix is unchanged by every trim (not both None). ``cost`` (default ones, >= 0) weights
    the objective sum_l c_l |t_l| (the trim magnitude), e.g. 1 for every gram of added or removed mass. Returns a
    TrimResult; see the module docstring (rel_floor is a conservative, possibly loose, estimate).

    ``objective="peak"`` minimises max_l c_l |t_l| instead of the sum (e.g. the largest tuning voltage squared).
    ``backtrack``: when a step increases the split while the previous split is well above the resolution floor
    (> 4 eps lambda_max), the step is halved, up to ``backtrack`` times in a row (history step "halve"); 0 disables
    it. Strongly nonlinear trims (a first-order effect that nearly cancels, e.g. depth ablation of a ring) need it."""
    _check_objective(objective)
    i, j = (int(pair[0]), int(pair[1]))
    if j != i + 1 or i < 0:
        raise ValueError("pair must be (i, i + 1) with i >= 0")
    if dKs is None and dMs is None:
        raise ValueError("dKs and dMs cannot both be None")
    L = len(dKs) if dKs is not None else len(dMs)
    dKs = [_zero_like(K)] * L if dKs is None else list(dKs)
    dMs_l = None if dMs is None else list(dMs)
    if len(dKs) != L or (dMs_l is not None and len(dMs_l) != L):
        raise ValueError("dKs and dMs must have the same length")
    if L == 0:
        raise ValueError("dKs / dMs must contain at least one candidate trim")
    lo, hi, c = _bounds(L, lower, upper, cost)

    def pencil(t):
        Kt, Mt = K, M
        for l in np.flatnonzero(t):
            Kt = Kt + t[l] * dKs[l]
            if dMs_l is not None:
                Mt = Mt + t[l] * dMs_l[l]
        return Kt, Mt

    def solve_lp(t, G, r):
        return _lp(G, r + G @ t, lo, hi, c, objective)

    lam_max = _lambda_max(K, M)
    floor_abs = float(np.finfo(float).eps * lam_max)

    t = np.clip(np.zeros(L), lo, hi)
    hist: list = []
    n_eig = 0
    best = None                        # (split, t, pair eigenvalues)
    first = None
    no_gain = 0
    status = "max_iter"
    t_prev, g_prev, n_half = None, None, 0
    for it in range(max_iter + 1):
        Kt, Mt = pencil(t)
        w, V = smallest_eigenpairs(Kt, Mt, min(j + 2, Kt.shape[0]))
        n_eig += 1
        try:
            sm = splitting_matrices(Kt, Mt, dKs, dMs_l, (i, j), eig=(w, V), gap_ratio_min=gap_ratio_min,
                                    validate=False)
        except ClusterSplitError as e:
            _not_nearly_degenerate(w, i, j, gap_ratio_min, e)
            raise
        g, lb = float(w[j] - w[i]), float(sm.lam_bar)
        rec = dict(it=it, t_cost=float(c @ np.abs(t)), split=g, rel_split=g / abs(lb), lam_bar=lb)
        if (backtrack and t_prev is not None and g > g_prev and g_prev > 4.0 * floor_abs and n_half < backtrack
                and it < max_iter):
            n_half += 1                # reject the iterate: halve the step (not counted as an iteration without gain)
            hist.append(TrimRecord(step="halve", **rec))
            t = t_prev + 0.5 * (t - t_prev)
            continue
        n_half = 0
        if first is None:
            first = (g, w[i:j + 1].copy())
        if best is None or g < best[0]:
            no_gain = 0 if best is None or g <= 0.5 * best[0] else no_gain + 1
            best = (g, t.copy(), w[i:j + 1].copy())
        else:
            no_gain += 1
        if g <= tol * abs(lb):
            status = "converged"
        elif no_gain >= 2:
            status = "stalled"
        elif it == max_iter:
            status = "max_iter"
        else:
            A = sm.mats
            G = np.vstack([A[:, 0, 0] - A[:, 1, 1], A[:, 0, 1]])
            r = np.array([g, 0.0])
            free = (np.abs(t) > 1e-14) & _strictly_inside(t, lo) & _strictly_inside(-t, -hi)
            step, t_new = "lp", None
            if it > 0 and np.any(free):
                S = np.flatnonzero(free)
                dt = np.zeros(L)
                dt[S] = np.linalg.lstsq(G[:, S], r, rcond=None)[0]
                tn = t + dt
                if np.all(tn >= lo) and np.all(tn <= hi) and np.all(tn[S] * t[S] > 0):   # in bounds, signs kept
                    step, t_new = "newton", t + dt
            if t_new is None:
                t_new, why = solve_lp(t, G, r)
                if t_new is None:
                    status = why
            if t_new is not None:
                hist.append(TrimRecord(step=step, **rec))
                t_prev, g_prev = t.copy(), g
                t = t_new
                continue
        hist.append(TrimRecord(step="", **rec))
        break
    g, t_best, pair_w = best
    lb_best = float(pair_w.mean())
    rel_floor = float(floor_abs / abs(lb_best))
    if status == "max_iter" and rel_floor >= tol and g / abs(lb_best) <= 2.0 * rel_floor:
        status = "stalled"             # out of iterations at the resolution floor (tol is below the floor)
    return TrimResult(t=t_best, cost=float(c @ np.abs(t_best)), split_before=first[0], split_after=g,
                      rel_split_before=first[0] / abs(float(first[1].mean())), rel_split_after=g / abs(lb_best),
                      eigvals_before=first[1], eigvals_after=pair_w, rel_floor=rel_floor, status=status,
                      n_iter=len(hist) - 1, n_eigensolves=n_eig, history=tuple(hist))
