"""B12: tau(P_-(h)) for Green's-function elements of the group von Neumann algebra L(Gamma_Bolza).

Definition of h (stated precisely).  For a base point z0 in the Poincare disk and a radius R,

    h_{R,z0} = sum_{gamma in Gamma, gamma != 1, d(z0, gamma z0) <= R}  G1(d(z0, gamma z0)) u_gamma,
    G1(d) = -(1/(2 pi)) ln tanh(d/2)      (free Green's function of -Delta on H^2, -Delta G1 = delta).

h is self-adjoint (G1(d(z0,gamma^-1 z0)) = G1(d(z0,gamma z0))), has positive coefficients and tau(h) = 0.
mu_h is its spectral measure w.r.t. the canonical trace tau(x) = <x delta_e, delta_e>, and
tau(P_-(h)) = mu_h((-inf, 0)).

Moments (exact group combinatorics):  v_j = h^j delta_e (finitely supported, keyed by the EXACT group element),
    m_{2j} = <v_j, v_j>,   m_{2j+1} = <v_j, v_{j+1}>.
NO finite word-ball sections are used: these are moments of the true L(Gamma) element.

Error model for float moments (all coefficients c_gamma > 0, so every sum is of positive terms): with unit
roundoff u = 2^-53, a sum of N positive terms has relative error <= (N-1)u; a product adds u.  Over j
convolution steps with |supp h| = S terms per entry, then a final inner product over N_j terms:
    rel.err(m_k) <= delta_k := 4 * (k (S + 2) + N + 10) * u      (factor 4 = safety, also covers c_gamma).
For the systole shell h = c * sum_k u_{g_k} the moments are c^k * (integer closed-walk counts), computed
EXACTLY with integers.

Bounds.
 (1) Christoffel bound (elementary): for any polynomial P with P(x0) = 1 and deg P <= n-1,
     mu({x0}) <= int |P|^2 d mu; minimising gives mu({x0}) <= lambda_n(x0) := 1 / (v^T H^{-1} v),
     H = (m_{i+j})_{i,j<n}, v = (x0^i).  Exact with Fractions when moments are integers.
 (2) Polynomial certificates (Chebyshev-Markov-Stieltjes type, self-certifying): with M >= ||h||
     (rigorous: M = ||h||_{l1} = sum c_gamma), find polynomials P >= 1_{[-M,0]} and Q <= 1_{[-M,0)} on
     [-M, M] (LP on a grid), then VERIFY the inequalities on all of [-M, M] by a Lipschitz/Markov bound
     |P'| <= sum k^2 |p_k| / M on a fine grid, shift P (resp. Q) by the certified violation, and conclude
         int Q d mu  <=  mu((-inf,0))  <=  mu((-inf,0])  <=  int P d mu,
     with int P d mu = sum p_k tau(T_k(h/M)) evaluated with the moment error budget.  No external theorem is
     needed; optimal certificates reproduce the Chebyshev-Markov-Stieltjes bounds.
 (3) KPM (Jackson kernel) point estimate of mu((-inf,0)) -- NOT a bound.
"""
from __future__ import annotations

import math
import time
from fractions import Fraction

import mpmath as mp
import numpy as np

from ..exact import group as G

U = 2.0 ** -53


def G1(d):
    return -np.log(np.tanh(np.asarray(d, dtype=float) / 2)) / (2 * math.pi)


# ---------------------------------------------------------------- the element h

def build_h(R: float, z0: complex = 0j):
    """Return (rows int64 [S,8], coeffs float [S], dists float [S])."""
    if z0 == 0:
        X = G.ball_by_displacement(R)
        ch = G.rows_cosh_disp0(X)
        d = np.arccosh(np.maximum(ch, 1.0))
    else:
        r0 = 2 * math.atanh(abs(z0))
        X = G.ball_by_displacement(R + 2 * r0)
        d = np.array([_dist_z0(G.elt(r), z0) for r in X])
    keep = (d > 1e-9) & (d <= R * (1 + 1e-12))
    X, d = X[keep], d[keep]
    return X, G1(d), d


def _dist_z0(g, z0):
    w = G.act(g, z0)
    return math.acosh(1 + 2 * abs(z0 - w) ** 2 / ((1 - abs(z0) ** 2) * (1 - abs(w) ** 2)))


def systole_shell():
    """h_sys = c * (u_{g0}+...+u_{g7}),  c = G1(l_sys)."""
    X = np.array([G.row(g) for g in G.GENS])
    return X, np.full(8, float(G1(G.L_SYS))), np.full(8, G.L_SYS)


# ---------------------------------------------------------------- convolution

def _void(X):
    X = np.ascontiguousarray(X)
    return X.view(np.dtype((np.void, X.dtype.itemsize * X.shape[1]))).ravel()


def _aggregate(keys, vals):
    v = _void(keys)
    uniq, inv = np.unique(v, return_inverse=True)
    out_vals = np.bincount(inv.ravel(), weights=vals) if vals.dtype.kind == "f" else _int_bincount(inv.ravel(), vals, len(uniq))
    out_keys = keys[np.unique(inv.ravel(), return_index=True)[1]]
    return out_keys, out_vals


def _int_bincount(inv, vals, n):
    out = np.zeros(n, dtype=np.int64)
    np.add.at(out, inv, vals)
    return out


def apply_h(hX, hc, keys, vals):
    """(h * v) where h = sum c_gamma u_gamma (left multiplication): supp -> gamma * eta."""
    Ls = [G.left_mult_matrix(G.elt(r)) for r in hX]
    allk = np.concatenate([keys @ L for L in Ls])
    allv = np.concatenate([vals * c for c in hc])
    return _aggregate(allk, allv)


def _dot_exact(a, b):
    """Exact dot product of int64 arrays (int64 if provably no overflow, else Python ints)."""
    if a.dtype.kind == "i" and b.dtype.kind == "i" and len(a):
        bound = int(np.abs(a).max()) * int(np.abs(b).max()) * len(a)
        if bound < 2 ** 62:
            return int(np.dot(a, b))
    return int(np.dot(a.astype(object), b.astype(object)))


def inner(k1, v1, k2, v2):
    a, b = _void(k1), _void(k2)
    common, i1, i2 = np.intersect1d(a, b, assume_unique=True, return_indices=True)
    if v1.dtype.kind == "i" and v2.dtype.kind == "i":
        return _dot_exact(v1[i1], v2[i2])
    return (v1[i1] * v2[i2]).sum()


_HASHV = np.array([1000003, 998244353, 1000000007, 754974721, 167772161, 469762049, 2013265921, 1811939329], dtype=np.int64)


def _final_step_partitioned(hX, hc, keys, vals, integer, max_rows):
    """<v_J, v_{J+1}> and ||v_{J+1}||^2 without storing v_{J+1}: output keys are split into P classes by a
    fixed linear hash; each class is aggregated separately (an element's class does not depend on how it
    was produced, so per-class aggregation is exact)."""
    Ls = [G.left_mult_matrix(G.elt(r)) for r in hX]
    total = len(keys) * len(hX)
    P = max(1, int(math.ceil(total / max_rows)))
    hk = (keys @ _HASHV)
    m_odd = 0 if integer else 0.0
    m_even = 0 if integer else 0.0
    nsupp = 0
    for part in range(P):
        ks, vs = [], []
        for L, c in zip(Ls, hc):
            K2 = keys @ L
            sel = ((K2 @ _HASHV) % P) == part
            ks.append(K2[sel])
            vs.append(vals[sel] * c)
        k2, v2 = _aggregate(np.concatenate(ks), np.concatenate(vs))
        nsupp += len(k2)
        selJ = (hk % P) == part
        if integer:
            m_odd += int(inner(keys[selJ], vals[selJ], k2, v2))
            m_even += _dot_exact(v2, v2)
        else:
            m_odd += float(inner(keys[selJ], vals[selJ], k2, v2))
            m_even += float(np.dot(v2, v2))
    return m_odd, m_even, nsupp, P


def moments(hX, hc, J, integer=False, verbose=False, max_rows=1.5e7):
    """m_0..m_{2J}.  v_1..v_{J-1} are stored; the last step is hash-partitioned (memory ~ max_rows rows).
    integer=True: hc must be ints, computation exact (Python ints for the inner products)."""
    keys = np.array([G.row(G.IDENTITY)])
    vals = np.array([1], dtype=np.int64) if integer else np.array([1.0])
    hc = np.asarray(hc, dtype=np.int64 if integer else float)
    m = [int(1) if integer else 1.0]
    sizes = [1]
    times = [0.0]
    for j in range(J):
        t = time.time()
        if j == J - 1:
            m_odd, m_even, nsupp, P = _final_step_partitioned(hX, hc, keys, vals, integer, max_rows)
            m += [m_odd, m_even]
            sizes.append(nsupp)
            times.append(time.time() - t)
            if verbose:
                print(f"j={j+1} |supp v_j|={nsupp} (partitioned P={P}) time={times[-1]:.1f}s", flush=True)
            break
        k2, v2 = apply_h(hX, hc, keys, vals)
        if integer:
            m_odd = int(inner(keys, vals, k2, v2))
            m_even = _dot_exact(v2, v2)
        else:
            m_odd = float(inner(keys, vals, k2, v2))
            m_even = float(np.dot(v2, v2))
        m += [m_odd, m_even]
        keys, vals = k2, v2
        sizes.append(len(keys))
        times.append(time.time() - t)
        if verbose:
            print(f"j={j+1} |supp v_j|={len(keys)} maxcoef={G.max_abs_coeff(keys)} time={times[-1]:.1f}s", flush=True)
        if G.max_abs_coeff(keys) > G.INT64_SAFE:
            raise OverflowError("int64 key coefficients too large")
    return m, sizes, times


def moment_rel_errors(m, S, sizes):
    """delta_k for float moments (see module docstring)."""
    out = []
    for k in range(len(m)):
        j = (k + 1) // 2
        N = sizes[min(j, len(sizes) - 1)]
        out.append(4 * (k * (S + 2) + N + 10) * U)
    return out


# ---------------------------------------------------------------- bounds

def christoffel_bound(m, x0=0, n=None):
    """Exact lambda_n(x0) = 1/(v^T H^{-1} v) using Fractions; m integers or Fractions."""
    if n is None:
        n = (len(m) - 1) // 2 + 1
    x0 = Fraction(x0)
    H = [[Fraction(m[i + j]) for j in range(n)] for i in range(n)]
    v = [x0 ** i for i in range(n)]
    # solve H y = v by exact Gaussian elimination
    A = [row[:] + [v[i]] for i, row in enumerate(H)]
    for c in range(n):
        p = next(r for r in range(c, n) if A[r][c] != 0)
        A[c], A[p] = A[p], A[c]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c] / A[c][c]
                A[r] = [a - f * b for a, b in zip(A[r], A[c])]
    y = [A[i][n] / A[i][i] for i in range(n)]
    q = sum(vi * yi for vi, yi in zip(v, y))
    return 1 / q


def cheb_moments(m, M, rel_err, dps=50):
    """mu_k = tau(T_k(h/M)) and absolute error bounds, from monomial moments m_j with relative errors."""
    with mp.workdps(dps):
        K = len(m) - 1
        # T_k coefficients
        T = [[mp.mpf(0)] * (K + 1) for _ in range(K + 1)]
        T[0][0] = mp.mpf(1)
        if K >= 1:
            T[1][1] = mp.mpf(1)
        for k in range(2, K + 1):
            for j in range(K + 1):
                T[k][j] = (2 * T[k - 1][j - 1] if j >= 1 else 0) - T[k - 2][j]
        Mm = mp.mpf(M)
        scaled = [mp.mpf(m[j]) / Mm ** j for j in range(K + 1)]
        mu = [mp.fsum(T[k][j] * scaled[j] for j in range(K + 1)) for k in range(K + 1)]
        err = [mp.fsum(abs(T[k][j]) * abs(scaled[j]) * rel_err[j] for j in range(K + 1)) for k in range(K + 1)]
        # + conversion of mu_k to float64 (relative 2^-53) and residual mp rounding
        return [float(x) for x in mu], [float(e) * 1.01 + abs(float(x)) * 2.0 ** -50 + 1e-40 for e, x in zip(err, mu)]


def _cheb_eval(p, x):
    """Clenshaw for sum p_k T_k(x), vectorised."""
    b1 = np.zeros_like(x)
    b2 = np.zeros_like(x)
    for c in p[:0:-1]:
        b1, b2 = 2 * x * b1 - b2 + c, b1
    return x * b1 - b2 + p[0]


def _lp_certificate(mu, upper, n_grid=3000, rounds=6):
    """LP on a grid + cutting planes: after each solve, the local minima of the violation on a fine grid are
    added as constraints (the final validity is NOT taken from the LP; `certify` checks it rigorously)."""
    from scipy.optimize import linprog
    K = len(mu) - 1
    th = np.linspace(0, np.pi, n_grid)
    xs = np.unique(np.concatenate([np.cos(th), np.linspace(-0.02, 0.02, 401), [0.0]]))
    fine = np.linspace(-1, 1, 400001)
    p = None
    for _ in range(rounds):
        V = np.cos(np.outer(np.arccos(np.clip(xs, -1, 1)), np.arange(K + 1)))
        target = (xs <= 0).astype(float) if upper else (xs < 0).astype(float)
        if upper:  # minimise sum p mu  s.t. V p >= target
            res = linprog(c=mu, A_ub=-V, b_ub=-target, bounds=[(-1e7, 1e7)] * (K + 1), method="highs")
        else:      # maximise sum q mu  s.t. V q <= target
            res = linprog(c=-np.array(mu), A_ub=V, b_ub=target, bounds=[(-1e7, 1e7)] * (K + 1), method="highs")
        if not res.success:
            raise RuntimeError(res.message)
        p = res.x
        v = _cheb_eval(p, fine)
        tgt = (fine <= 0).astype(float) if upper else (fine < 0).astype(float)
        viol = (tgt - v) if upper else (v - tgt)
        bad = fine[viol > 1e-9]
        if len(bad) == 0:
            break
        xs = np.unique(np.concatenate([xs, bad[:: max(1, len(bad) // 2000)]]))
    return p


def certify(p, upper, n_check=None, target_cell_err=1e-7):
    """Certified shift eps >= 0 such that P + eps >= 1_{[-1,0]} (upper) or Q - eps <= 1_{[-1,0)} (lower)
    on all of [-1, 1].  On each cell [x_i, x_{i+1}] of a uniform grid (n_check points on each of [-1,0] and
    [0,1]) the polynomial differs from its linear interpolant by at most h^2/8 * max|P''|, and
    |P''| <= sum_k |p_k| k^2 (k^2 - 1)/3 on [-1,1] (|T_k''| <= T_k''(1) = k^2(k^2-1)/3).  Float evaluation
    error of the Clenshaw sum is bounded crudely by 16 (K+2) u * max(sum|p_k|, max|P(x_i)|)."""
    K = len(p) - 1
    d2 = float(sum(abs(c) * k * k * (k * k - 1) / 3 for k, c in enumerate(p)))
    if n_check is None:  # choose the grid so that h^2/8 * d2 <= target_cell_err (capped)
        n_check = int(min(4e7, max(2e6, math.sqrt(d2 / (8 * target_cell_err)) + 2)))
    worst = 0.0
    for a, b, tgt in [(-1.0, 0.0, 1.0), (0.0, 1.0, 0.0)]:
        xs = np.linspace(a, b, n_check)
        h = (b - a) / (n_check - 1)
        vals = _cheb_eval(np.asarray(p, float), xs)
        evalerr = 16 * (K + 2) * U * max(float(sum(abs(c) for c in p)), float(np.abs(vals).max()))
        # min over a cell >= min of endpoint values - h^2/8 * d2 ; max <= max endpoints + h^2/8 d2
        if upper:
            viol = np.max(tgt - np.minimum(vals[:-1], vals[1:]))
        else:
            viol = np.max(np.maximum(vals[:-1], vals[1:]) - tgt)
        worst = max(worst, viol + h * h / 8 * d2 + evalerr)
    return max(0.0, worst)


def rigorous_bounds(mu, mu_err, n_grid=3000):
    """Return (lower for mu((-inf,0)), upper for mu((-inf,0])) for the measure of h/M (same as for h)."""
    out = {}
    for upper in (True, False):
        p = _lp_certificate(mu, upper, n_grid)
        eps = certify(p, upper)
        val = float(np.dot(p, mu))
        err = float(np.dot(np.abs(p), mu_err))
        if upper:
            out["upper_le0"] = min(1.0, val + eps * mu[0] + err)   # mu_0 = 1
        else:
            out["lower_lt0"] = max(0.0, val - eps * mu[0] - err)
        out["eps_" + ("upper" if upper else "lower")] = eps
        out["moment_err_" + ("upper" if upper else "lower")] = err
    return out


def kpm_estimate(mu):
    """Jackson-damped Chebyshev estimate of mu([-1, 0)) (NOT a bound)."""
    N = len(mu)
    g = [((N - k + 1) * math.cos(math.pi * k / (N + 1)) + math.sin(math.pi * k / (N + 1)) / math.tan(math.pi / (N + 1))) / (N + 1)
         for k in range(N)]
    s = g[0] * mu[0] / 2
    for k in range(1, N):
        s -= (2 / math.pi) * g[k] * mu[k] * math.sin(k * math.pi / 2) / k
    return s


from .. import data as _data  # noqa: E402
from ..provenance import provenance  # noqa: E402


@provenance("CERTIFIED", "exact-group moments τ(h^k) + self-certifying polynomial bounds (item 02, run_tau.py); "
            "recompute with `python scripts/regen_tau.py` (~30 min, ~8 GB RAM)")
def tau_bounds(case="systole_shell") -> dict:
    """Two-sided bounds on τ(P₋(h)), the von Neumann trace of the spectral projection of h onto (−∞, 0) (item 02 B12).
    case: 'systole_shell' (h = c·Σ_{k=0}^{7} u_{g_k}) or 0..3 (G₁-weighted ball operators; see details['R'], ['z0'])."""
    import json as _json
    d = _data.load("tau")
    if case == "systole_shell":
        e = d["systole_shell"]
    elif isinstance(case, int) and not isinstance(case, bool) and 0 <= case < len(d["general"]):
        e = d["general"][case]
    else:
        raise ValueError(f"case must be 'systole_shell' or an int in 0..{len(d['general']) - 1}")
    iv = e["tau_P_minus_rigorous_interval"]
    lo, hi = _json.loads(iv) if isinstance(iv, str) else iv   # stored as a JSON list in the bundled file; str accepted too
    return {"lower": float(lo), "upper": float(hi), "h": e["h"], "details": e}
