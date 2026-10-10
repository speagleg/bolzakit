"""Frequency trimming of a split degenerate mode pair (an inverse problem; NUMERICAL, discrete problem).

A nominally degenerate pair (i, i + 1) of K x = lambda M x is split by imperfections. Candidate trims l = 0..L-1 (a
point mass at a site, removed material, a stiffness patch) change the pencil linearly in their amounts t_l:
    K(t) = K + sum_l t_l dK_l,    M(t) = M + sum_l t_l dM_l.
In the pair's M-orthonormal eigenbasis U the reduced 2x2 model is D + A(dt), D = diag(-g/2, g/2) (g the split),
A(dt) = sum_l dt_l A_l with the splitting matrices A_l = U^T (dK_l - lam_bar dM_l) U. The pair is tied to first order
iff (D + A)_00 = (D + A)_11 and (D + A)_01 = 0, two linear equations G dt = (g, 0) with
G = [A_00 - A_11; A_01]. The first step is the linear program
    min c^T (t + dt)  subject to  G dt = (g, 0),  lower <= t + dt <= upper
(a vertex uses at most two trims). Later steps are Newton corrections on the trims the LP selected (least squares on
that support, exact re-solve each time); the LP is re-solved only if a correction leaves the bounds. Convergence is
quadratic until the split reaches the rounding floor of the eigensolver (about eps * lambda_max / lambda_pair).

Statuses: "converged" (split <= tol * |lam_bar|), "stalled" (two consecutive iterations without halving the best
split: the rounding floor), "max_iter", "infeasible" (no trim pattern within the bounds ties the pair to first order;
t is the best iterate found). The pair must stay isolated (the core gap rule, ClusterSplitError otherwise).

NUMERICAL: computed for the discrete pencil (K, M) as given; discretisation error relative to any continuum operator
is not controlled. Imports only numpy/scipy and bolzakit.sensitivity.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
from scipy.optimize import linprog

from bolzakit.sensitivity.core import DISCRETE_CAVEAT, smallest_eigenpairs, splitting_matrices

__all__ = ["TrimRecord", "TrimResult", "trim_degenerate_pair"]


@dataclass(frozen=True)
class TrimRecord:
    it: int
    t_cost: float          # c^T t at this iterate
    split: float           # lambda_{i+1} - lambda_i
    rel_split: float       # split / |lam_bar|
    lam_bar: float
    step: str              # how the NEXT iterate was obtained: "lp", "newton" or "" (last record)


@dataclass(frozen=True)
class TrimResult:
    t: np.ndarray                  # trim amounts of the best iterate (lowest split)
    cost: float                    # c^T t
    split_before: float
    split_after: float
    rel_split_before: float
    rel_split_after: float
    eigvals_before: np.ndarray     # the pair at t = t0
    eigvals_after: np.ndarray      # the pair at t
    status: str                    # "converged", "stalled", "max_iter" or "infeasible"
    n_iter: int
    n_eigensolves: int
    history: tuple
    label: str = "NUMERICAL"
    caveat: str = DISCRETE_CAVEAT


def _zero_like(K):
    return sp.csc_matrix(K.shape) if sp.issparse(K) else np.zeros(K.shape)


def _strictly_inside(t, lo):
    """t_l > lo_l by more than a rounding margin (always true for lo_l = -inf)."""
    lo_f = np.where(np.isfinite(lo), lo, 0.0)
    return ~np.isfinite(lo) | (t > lo_f + 1e-14 * np.maximum(1.0, np.abs(lo_f)))


def trim_degenerate_pair(K, M, dKs, dMs, pair, *, lower=0.0, upper=None, cost=None, tol: float = 1e-10,
                         max_iter: int = 20, gap_ratio_min: float = 1e-3) -> TrimResult:
    """Smallest-cost trim t (lower <= t <= upper) that ties the split pair (i, i + 1) of K x = lambda M x.

    ``dKs`` / ``dMs``: length-L sequences of the per-unit changes of K and M for each candidate trim (dense or
    sparse); None means that matrix is unchanged by every trim (not both None). ``cost`` (default ones) weights the
    LP objective c^T t, e.g. 1 for every gram of added mass. Returns a TrimResult; see the module docstring."""
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
    lo = np.broadcast_to(np.asarray(lower, float), (L,)).copy()
    hi = np.full(L, np.inf) if upper is None else np.broadcast_to(np.asarray(upper, float), (L,)).copy()
    if np.any(lo > hi):
        raise ValueError("lower must not exceed upper")
    c = np.ones(L) if cost is None else np.broadcast_to(np.asarray(cost, float), (L,)).copy()

    def pencil(t):
        Kt, Mt = K, M
        for l in np.flatnonzero(t):
            Kt = Kt + t[l] * dKs[l]
            if dMs_l is not None:
                Mt = Mt + t[l] * dMs_l[l]
        return Kt, Mt

    def solve_lp(t, G, r):
        res = linprog(c, A_eq=G, b_eq=r + G @ t, bounds=list(zip(lo, hi)),
                      method="highs")
        return None if res.status != 0 else np.asarray(res.x, dtype=float)

    t = np.clip(np.zeros(L), lo, hi)
    hist: list = []
    n_eig = 0
    best = None                        # (split, t, pair eigenvalues)
    first = None
    no_gain = 0
    status = "max_iter"
    for it in range(max_iter + 1):
        Kt, Mt = pencil(t)
        w, V = smallest_eigenpairs(Kt, Mt, min(j + 2, Kt.shape[0]))
        n_eig += 1
        sm = splitting_matrices(Kt, Mt, dKs, dMs_l, (i, j), eig=(w, V), gap_ratio_min=gap_ratio_min, validate=False)
        g, lb = float(w[j] - w[i]), float(sm.lam_bar)
        if first is None:
            first = (g, w[i:j + 1].copy())
        if best is None or g < best[0]:
            no_gain = 0 if best is None or g <= 0.5 * best[0] else no_gain + 1
            best = (g, t.copy(), w[i:j + 1].copy())
        else:
            no_gain += 1
        rec = dict(it=it, t_cost=float(c @ t), split=g, rel_split=g / abs(lb), lam_bar=lb)
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
            free = _strictly_inside(t, lo) & _strictly_inside(-t, -hi)
            step, t_new = "lp", None
            if it > 0 and np.any(free):
                S = np.flatnonzero(free)
                dt = np.zeros(L)
                dt[S] = np.linalg.lstsq(G[:, S], r, rcond=None)[0]
                if np.all(t + dt >= lo) and np.all(t + dt <= hi):
                    step, t_new = "newton", t + dt
            if t_new is None:
                t_new = solve_lp(t, G, r)
                if t_new is None:
                    status = "infeasible"
            if t_new is not None:
                hist.append(TrimRecord(step=step, **rec))
                t = t_new
                continue
        hist.append(TrimRecord(step="", **rec))
        break
    g, t_best, pair_w = best
    lb_best = float(pair_w.mean())
    return TrimResult(t=t_best, cost=float(c @ t_best), split_before=first[0], split_after=g,
                      rel_split_before=first[0] / abs(float(first[1].mean())), rel_split_after=g / abs(lb_best),
                      eigvals_before=first[1], eigvals_after=pair_w, status=status, n_iter=len(hist) - 1,
                      n_eigensolves=n_eig, history=tuple(hist))
