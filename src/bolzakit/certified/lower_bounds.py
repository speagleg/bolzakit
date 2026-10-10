"""Certified lower bound for lambda_1(Bolza) by Lehmann-Goerisch, P2 trial functions and RT2 fluxes.

A port of the Bolza lambda1 research code (the original research code: lb_core.py, lb_bolza.py; reviewed in
review_lb_1.md, ACCEPT).  Changes from SRC: the precision (128 bits) and the series order (K = 14) are fixed and set
internally, never at import; the precision-dependent tables are built lazily; the P2 shapes and polynomial helpers
are reused from ``forms`` (identical code); every SRC ``assert`` is an explicit fail-closed check; non-finite element
data and a non-negligible Green residual make the result uncertified.  Also changed (all on the untrusted float side
or in diagnostics; the certified step is unchanged):
  * the trial eigenvectors come from ``eigsh(k=5)`` (ntrial + 2) with a fixed start vector
    ``v0 = RandomState(0).rand(np2)`` (SRC: ``k=6`` and ARPACK's random start), so the trial vectors, hence the
    certificate, are reproducible;
  * SRC's per-element cross-check of K, M against cert_rr_p2.element_p2 (a diagnostic, ``crosscheck_vs_cert_rr_p2_max``)
    is dropped: the P2 shapes are the shared ``forms`` code, and the independent checker re-derives the LG step;
  * the resulting enclosure differs from SRC's published one by about 2.6e-16 (other trial vectors in the degenerate
    lambda_1 eigenspace); both are valid lower bounds;
  * ``data["nubar"]`` is stored as a plain float string (``repr(float(x))``), not a numpy repr.

Method (LB_PLAN.md §1(c), Liu-Plum arXiv:2512.23182 Thm 4.1/4.3, Behnke-Mertins-Plum-Wieners 2000 Thm 2.1).
Shifted problem on D = H^1(B): M = q + m, N = m, X = L^2(TB) x L^2(B), T u = (grad u, u), b_G = L^2 product, so the
shifted eigenvalues are lambda~ = lambda + 1: lambda~_1 = 1 (constants), lambda~_{2,3,4} = lambda_{1,2,3} + 1,
lambda~_5 = lambda_4 + 1.  n = 3 trial functions v_i (P2, float coefficients used as exact data), w_i =
(sigma_i, v_i + div sigma_i) with sigma_i in RT2 (H(div)-conforming on the glued surface), so b_G(w_i, T phi) =
N(v_i, phi) exactly.  A0 = M(v_i, v_j), A1 = N(v_i, v_j), A2 = b_G(w_i, w_j); A = A0 - rho A1,
B = A0 - 2 rho A1 + rho^2 A2 (certified positive definite).  If rho <= lambda~_5 (m = 4), and nu_3 < nubar < 0 is
certified (nubar B - A positive definite), then lambda_1 + 1 = lambda~_2 >= rho - rho/(1 - nubar).

HONESTY NOTE (review_lb_1 row 4c): the computation cannot detect a wrong rho.  rho = 7.5 violates the hypothesis
(lambda~_5 ~ 6.354) and still returns a plausible sub-true value.  The validity of rho = 5 rests entirely on the gap
input lambda_4 >= 4.0 (the trace-formula count certificates lp_A / lp_B), which is recorded as a dependency.

Per triangle T the geometry enters only through (Q, M, D) in the reference simplex coordinates (a, b):
   sqrt(g) = D Q^{-3/2},  sqrt(g) g^{-1} = adj(M) / (D sqrt Q),  g / sqrt(g) = M / (D sqrt Q).
Hyperbolic triangle (hyperboloid-barycentric map): Q = B(W,W), M = L L^T - Q G.  Flat triangle: Q = 1, M = J^T J,
D = |det J|.  RT2: j = sqrt(g) sigma, j in P2^2 + x P2hom (dim 15); DOFs: outward flux density per unit edge
parameter at s in {0, 1/2, 1} (s = barycentric weight of the edge endpoint with the smaller global label) and the
6 interior moments int j_c q, q in {1, a, b}.  Element integrals int_simplex poly Q^{-s} by the binomial series in
eps = Q/Qc - 1 with a Lagrange remainder (as in forms.element_p2).
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from math import factorial

import numpy as np
from flint import arb, arb_mat

from ..errors import CertificationFailed, OutOfValidityRange
from . import forms
from ._arb import ball_str, certainly_gt, certainly_lt, json_float, prec, rational_lower
from ._hyp import Bm, det3
from .certificate import Certificate, Claim, Dependency, make_producer, trial_provenance
from .cluster import GapInput
from .forms import binom_neg, mons, padd, peval, pmul, pscale, psup_centered
from .jets import J
from .mesh import Mesh, bolza_vertices_hyp, positions_jet

PREC = 128        # SRC LB_PREC default
K_SERIES = 14     # SRC K_SERIES default (fixed; no environment override)
DEGMU = 8
GREEN_TOL = arb("1e-20")   # |Green residual| bound: conforming sigma gives ~1e-36 (review_lb_1 row 2a)

ipmul, ipadd, ider = forms._ipmul, forms._ipadd, forms._ider


def ieval(p, a, b):
    return sum((c * a ** i * b ** j for (i, j), c in p.items()), Fraction(0))


def wq(i, j):
    return Fraction(factorial(i) * factorial(j), factorial(i + j + 2))


def iint(p):
    return sum((c * wq(i, j) for (i, j), c in p.items()), Fraction(0))


# ------------------------------------------------------------------ reference P2 basis (forms.SHAPES == SRC P2S)
LAM = forms._LAM
VERTS = [(1, 0), (0, 1), (0, 0)]          # local vertex r <-> barycentric lambda_r = 1
LEDGES = ((0, 1), (1, 2), (2, 0))
P2S = forms.SHAPES
P2G = forms.GRADS


# ------------------------------------------------------------------ reference RT2 basis (exact rationals)
def _rt2_basis():
    import sympy as sp
    raw = []
    for m in mons(2):
        raw.append(({m: 1}, {}))
    for m in mons(2):
        raw.append(({}, {m: 1}))
    for (i, j) in [(2, 0), (1, 1), (0, 2)]:
        raw.append(({(i + 1, j): 1}, {(i, j + 1): 1}))
    nodes = [Fraction(0), Fraction(1, 2), Fraction(1)]
    dofs = []
    for (i, j) in LEDGES:
        Vi, Vj = VERTS[i], VERTS[j]
        d = (Vj[0] - Vi[0], Vj[1] - Vi[1])
        nrm = (d[1], -d[0])                      # outward normal (CCW triangle), length |edge|
        for s in nodes:
            x = (s * Vi[0] + (1 - s) * Vj[0], s * Vi[1] + (1 - s) * Vj[1])
            dofs.append(("e", nrm, x))
    for c in range(2):
        for q in ({(0, 0): 1}, {(1, 0): 1}, {(0, 1): 1}):
            dofs.append(("i", c, q))

    def apply(dof, v):
        if dof[0] == "e":
            _, nrm, x = dof
            return nrm[0] * ieval(v[0], *x) + nrm[1] * ieval(v[1], *x)
        _, c, q = dof
        return iint(ipmul(v[c], q))
    Dm = sp.Matrix([[sp.Rational(apply(d, v).numerator, apply(d, v).denominator) for v in raw] for d in dofs])
    C = Dm.inv()
    basis = []
    for l in range(15):
        comp = [{}, {}]
        for r in range(15):
            cr = Fraction(int(C[r, l].p), int(C[r, l].q))
            if cr != 0:
                comp[0] = ipadd(comp[0], {k: v * cr for k, v in raw[r][0].items()})
                comp[1] = ipadd(comp[1], {k: v * cr for k, v in raw[r][1].items()})
        basis.append(comp)
    # self-check: DOF_k(phi_l) = delta_kl
    for k, d in enumerate(dofs):
        for l in range(15):
            if apply(d, basis[l]) != (1 if k == l else 0):
                raise RuntimeError("RT2 basis self-check failed")
    return basis


# ------------------------------------------------------------------ coefficient tables (integer, common denominator)
def _table(entries, cols):
    """entries: list (rows) of lists of (col-block, poly) ; returns (integer arb_mat, denominator)"""
    colidx = {c: i for i, c in enumerate(cols)}
    den = 1
    for row in entries:
        for blk, p in row:
            for v in p.values():
                den = den * Fraction(v).denominator // math.gcd(den, Fraction(v).denominator)
    mat = [[0] * len(cols) for _ in entries]
    for r, row in enumerate(entries):
        for blk, p in row:
            for m, v in p.items():
                mat[r][colidx[(blk, m)]] += int(Fraction(v) * den)
    return arb_mat(mat), den


def _arbfrac_mat(F):
    den = 1
    for r in F:
        for v in r:
            den = den * v.denominator // math.gcd(den, v.denominator)
    return arb_mat([[int(v * den) for v in r] for r in F]) * (arb(1) / den)


class _Exact:
    pass


@lru_cache(maxsize=None)
def _exact():
    """Exact (precision-independent) tables: RT2 basis, its divergence, the integer coefficient tables and the
    rational cross matrices BX, GX (built once, lazily: sympy is only needed here)."""
    t = _Exact()
    t.RT = RT = _rt2_basis()
    t.RTDIV = RTDIV = [ipadd(ider(f[0], 0), ider(f[1], 1)) for f in RT]
    M2, M4, M6 = mons(2), mons(4), mons(6)
    t.M2, t.M4, t.M6 = M2, M4, M6
    # P2 stiffness: K_pq = sum_{c,d} int dphi_p/dc adjM_cd dphi_q/dd Q^{-1/2}/D ; adjM = [[M22,-M12],[-M12,M11]]
    rows = []
    for p in range(6):
        for q in range(6):
            row = []
            for c in range(2):
                for d in range(2):
                    pr = ipmul(P2G[p][c], P2G[q][d])
                    blk = {(0, 0): "A11", (1, 1): "A22"}.get((c, d), "A12")
                    row.append((blk, pr))
            rows.append(row)
    t.TAB_K, t.DEN_K = _table(rows, [(b, m) for b in ("A11", "A12", "A22") for m in M2])
    # P2 mass: D int phi_p phi_q Q^{-3/2}
    t.TAB_M, t.DEN_M = _table([[("R", ipmul(P2S[p], P2S[q]))] for p in range(6) for q in range(6)],
                              [("R", m) for m in M4])
    # RT sigma: int phi_l^T M phi_m Q^{-1/2}/D
    rows = []
    for l in range(15):
        for m in range(15):
            row = []
            for c in range(2):
                for d in range(2):
                    blk = {(0, 0): "M11", (1, 1): "M22"}.get((c, d), "M12")
                    row.append((blk, ipmul(RT[l][c], RT[m][d])))
            rows.append(row)
    t.TAB_S, t.DEN_S = _table(rows, [(b, m) for b in ("M11", "M12", "M22") for m in M6])
    # RT div-div: int div_l div_m Q^2 Q^{-1/2}/D
    t.TAB_DD, t.DEN_DD = _table([[("Z", ipmul(RTDIV[l], RTDIV[m]))] for l in range(15) for m in range(15)],
                                [("Z", m) for m in M4])
    # cross (exact): BX_lp = int div_l psi_p ; GX_lp = int phi_l . grad psi_p   (both polynomial)
    t.BX = [[iint(ipmul(RTDIV[l], P2S[p])) for p in range(6)] for l in range(15)]
    t.GX = [[iint(ipadd(ipmul(RT[l][0], P2G[p][0]), ipmul(RT[l][1], P2G[p][1]))) for p in range(6)]
            for l in range(15)]
    t.BX_F = np.array([[float(v) for v in r] for r in t.BX])
    t.GX_F = np.array([[float(v) for v in r] for r in t.GX])
    return t


class _Tables:
    pass


_TABLES = {}


def _tables():
    """Precision-dependent tables (moments, binomials, rational cross matrices) at the current flint precision."""
    import flint
    key = flint.ctx.prec
    if key not in _TABLES:
        e = _exact()
        t = _Tables()
        t.BX_A = _arbfrac_mat(e.BX)
        t.GX_A = _arbfrac_mat(e.GX)
        t.MONS_MU = mons(DEGMU)
        t.MONS_SER = mons(2 * K_SERIES)
        FACT = [arb(1)]
        for i in range(1, 80):
            FACT.append(FACT[-1] * i)
        t.WMAT = arb_mat([[FACT[m[0] + j[0]] * FACT[m[1] + j[1]] / FACT[m[0] + j[0] + m[1] + j[1] + 2]
                           for j in t.MONS_SER] for m in t.MONS_MU])
        t.WMU = [FACT[m[0]] * FACT[m[1]] / FACT[m[0] + m[1] + 2] for m in t.MONS_MU]
        t.S_LIST = [arb(1) / 2, arb(3) / 2]
        t.BINOM = arb_mat([[binom_neg(s, k) for s in t.S_LIST] for k in range(K_SERIES + 1)])
        t.IDX_MU = {m: i for i, m in enumerate(t.MONS_MU)}
        _TABLES[key] = t
    return _TABLES[key]


# ------------------------------------------------------------------ moments of Q^{-s}
def moments(Q):
    """returns mu[si][m] enclosing int_simplex m Q^{-s}, s = 1/2, 3/2, m of degree <= 8."""
    tb = _tables()
    third = arb(1) / 3
    Qc = peval(Q, third, third)
    eps = padd(pscale(Q, 1 / Qc), {(0, 0): arb(1)}, -1)
    emax = psup_centered(eps)
    if not certainly_lt(emax, 0.25):
        raise OutOfValidityRange("triangle too large for the series (emax >= 0.25)")
    epows = [{(0, 0): arb(1)}]
    for _ in range(K_SERIES):
        epows.append(pmul(epows[-1], eps))
    E = arb_mat([[ep.get(j, arb(0)) for ep in epows] for j in tb.MONS_SER])
    Z = tb.WMAT * (E * tb.BINOM)
    out = []
    for si, s in enumerate(tb.S_LIST):
        if emax == 0:
            rr = 0.0
        else:
            rem = abs(binom_neg(s, K_SERIES + 1)) * emax ** (K_SERIES + 1) * (1 - emax) ** (-s - K_SERIES - 1)
            rr = float(rem.upper()) * (1 + 1e-12)
        sc = Qc ** (-s)
        col = []
        for r, m in enumerate(tb.MONS_MU):
            v = Z[r, si]
            if rr:
                v = v + arb(0, float((tb.WMU[r] * rr).upper()) * (1 + 1e-12))
            col.append(v * sc)
        out.append(col)
    return out


def wvec(mu, P, monos, scale):
    """vector over monos of  scale * int mono * P * Q^{-s}  (mu = moment list for that s)."""
    idx = _tables().IDX_MU
    return [scale * sum((c * mu[idx[(k[0] + m[0], k[1] + m[1])]] for k, c in P.items()), arb(0)) for m in monos]


def element(Q, Mpol, D):
    """Q: poly dict (arb), Mpol = (M11, M12, M22) poly dicts (arb), D: arb > 0.
    Returns arb_mats: K (6x6), Mm (6x6), S = Ssigma + Sdd (15x15)."""
    e = _exact()
    M2, M4, M6 = e.M2, e.M4, e.M6
    mu = moments(Q)
    M11, M12, M22 = Mpol
    iD = 1 / D
    wk = wvec(mu[0], M22, M2, iD) + wvec(mu[0], pscale(M12, -1), M2, iD) + wvec(mu[0], M11, M2, iD)
    K = arb_mat(6, 6)
    kv = e.TAB_K * arb_mat([[x] for x in wk])
    wm = wvec(mu[1], {(0, 0): arb(1)}, M4, D)
    mv = e.TAB_M * arb_mat([[x] for x in wm])
    ws = wvec(mu[0], M11, M6, iD) + wvec(mu[0], M12, M6, iD) + wvec(mu[0], M22, M6, iD)
    sv = e.TAB_S * arb_mat([[x] for x in ws])
    Q2 = pmul(Q, Q)
    wd = wvec(mu[0], Q2, M4, iD)
    dv = e.TAB_DD * arb_mat([[x] for x in wd])
    Mm = arb_mat(6, 6)
    for p in range(6):
        for q in range(6):
            K[p, q] = kv[6 * p + q, 0] / e.DEN_K
            Mm[p, q] = mv[6 * p + q, 0] / e.DEN_M
    S = arb_mat(15, 15)
    for l in range(15):
        for m in range(15):
            S[l, m] = sv[15 * l + m, 0] / e.DEN_S + dv[15 * l + m, 0] / e.DEN_DD
    return K, Mm, S


def element_finite(el) -> bool:
    return all(Mx[i, j].is_finite() for Mx in el for i in range(Mx.nrows()) for j in range(Mx.ncols()))


# ------------------------------------------------------------------ topology: P2 and RT2 dof maps
class Topo:
    """tri_labels: (Nt,3) int array of global vertex labels (local vertex r <-> barycentric lambda_r)."""

    def __init__(self, tri_labels, nverts):
        ekey = {}
        owner = {}
        p2loc, rtloc, rtsign = [], [], []
        Nt = len(tri_labels)
        for t, L in enumerate(tri_labels):
            L = [int(x) for x in L]
            if len(set(L)) != 3:
                raise ValueError("degenerate triangle labels")
            eids, rl, rs = [], [], []
            for r, (u, v) in enumerate(LEDGES):
                k = (min(L[u], L[v]), max(L[u], L[v]))
                if k not in ekey:
                    ekey[k] = len(ekey); owner[k] = [t]
                else:
                    owner[k].append(t)
                e = ekey[k]
                eids.append(nverts + e)
                sg = 1 if owner[k][0] == t else -1
                # local s = lambda_u ; global s = lambda of the endpoint with the smaller label
                flip = (L[u] != k[0])
                for tl in range(3):
                    tg = 2 - tl if flip else tl
                    rl.append(3 * e + tg); rs.append(sg)
            p2loc.append(L + eids)
            rtloc.append(rl + [None] * 6); rtsign.append(rs + [1] * 6)
        self.Ne = len(ekey)
        if not all(len(v) == 2 and v[0] != v[1] for v in owner.values()):
            raise ValueError("edge not shared by 2 distinct triangles")
        self.ekey, self.owner = ekey, owner
        self.chi = nverts - self.Ne + Nt
        base = 3 * self.Ne
        for t in range(Nt):
            for i in range(6):
                rtloc[t][9 + i] = base + 6 * t + i
        self.p2loc = np.array(p2loc); self.rtloc = np.array(rtloc); self.rtsign = np.array(rtsign, float)
        self.np2 = nverts + self.Ne
        self.nrt = base + 6 * Nt
        self.Nt = Nt


# ------------------------------------------------------------------ assembly + certificate
def assemble_float(topo, elems):
    import scipy.sparse as sps
    e = _exact()
    Nt = topo.Nt
    r6 = np.repeat(topo.p2loc, 6, axis=1).ravel(); c6 = np.tile(topo.p2loc, (1, 6)).ravel()
    r15 = np.repeat(topo.rtloc, 15, axis=1).ravel(); c15 = np.tile(topo.rtloc, (1, 15)).ravel()
    sg15 = (topo.rtsign[:, :, None] * topo.rtsign[:, None, :]).reshape(Nt, 225)
    kv = np.empty((Nt, 36)); mv = np.empty((Nt, 36)); sv = np.empty((Nt, 225))
    for t, (K, Mm, S) in enumerate(elems):
        kv[t] = [float(K[p, q].mid()) for p in range(6) for q in range(6)]
        mv[t] = [float(Mm[p, q].mid()) for p in range(6) for q in range(6)]
        sv[t] = [float(S[l, m].mid()) for l in range(15) for m in range(15)]
    Kg = sps.csr_matrix((kv.ravel(), (r6, c6)), shape=(topo.np2, topo.np2))
    Mg = sps.csr_matrix((mv.ravel(), (r6, c6)), shape=(topo.np2, topo.np2))
    Sg = sps.csr_matrix(((sv * sg15).ravel(), (r15, c15)), shape=(topo.nrt, topo.nrt))
    rb = np.repeat(topo.rtloc, 6, axis=1).ravel(); cb = np.tile(topo.p2loc, (1, 15)).ravel()
    bv = (topo.rtsign[:, :, None] * e.BX_F[None, :, :]).reshape(Nt, 90)
    gv = (topo.rtsign[:, :, None] * e.GX_F[None, :, :]).reshape(Nt, 90)
    Bg = sps.csr_matrix((bv.ravel(), (rb, cb)), shape=(topo.nrt, topo.np2))
    Gg = sps.csr_matrix((gv.ravel(), (rb, cb)), shape=(topo.nrt, topo.np2))
    return Kg, Mg, Sg, Bg, Gg


def solve_fluxes(Sg, Bg, V):
    """c_i = argmin ||sigma||^2 + ||v_i + div sigma||^2  ->  S c = -B v."""
    import scipy.sparse.linalg as sla
    lu = sla.splu(Sg.tocsc())
    return np.column_stack([-lu.solve(Bg @ V[:, i]) for i in range(V.shape[1])])


def arb_forms(topo, elems, V, C):
    """exact (Arb) A0 = V^T(K+M)V, A1 = V^T M V, A2 = b_G(w_i,w_j) with w_i = (sigma_i, v_i + div sigma_i),
    and the Green-formula residual R = sum_T C^T (GX + BX) V (should be ~0: int sigma.grad v + div sigma v = 0)."""
    tb = _tables()
    n = V.shape[1]
    A0 = arb_mat(n, n); A1 = arb_mat(n, n); A2 = arb_mat(n, n); R = arb_mat(n, n)
    GB = tb.GX_A + tb.BX_A
    for t, (K, Mm, S) in enumerate(elems):
        Vl = arb_mat(V[topo.p2loc[t]].tolist())
        Cl = arb_mat((topo.rtsign[t][:, None] * C[topo.rtloc[t]]).tolist())
        VtM = Vl.transpose() * Mm
        A1 += VtM * Vl
        A0 += Vl.transpose() * (K * Vl) + VtM * Vl
        X = Cl.transpose() * (tb.BX_A * Vl)
        A2 += Cl.transpose() * (S * Cl) + X + X.transpose() + VtM * Vl
        R += Cl.transpose() * (GB * Vl)
    return A0, A1, A2, R


def is_pd(A):
    """Sylvester's criterion: every leading principal minor certainly > 0 (NaN never passes)."""
    m = A.nrows()
    for k in range(1, m + 1):
        d = arb_mat([[A[i, j] for j in range(k)] for i in range(k)]).det()
        if not certainly_gt(d, 0):
            return False
    return True


def lehmann(A0, A1, A2, rho):
    """Lehmann-Goerisch (Liu-Plum arXiv:2512.23182 Thm 4.1/4.3): with rho <= lambda~_{m+1}, A = A0 - rho A1,
    B = A0 - 2 rho A1 + rho^2 A2 (must be PD), nu_1<=..<=nu_n eigenvalues of A z = nu B z (all negative):
    lambda~_{m+1-k} >= rho - rho/(1 - nu_k).  Certifies k = n through an upper bound nubar >= nu_n (nubar B - A PD,
    nubar < 0).  Fail closed: returns dict with ok=False (and a reason) instead of SRC's asserts."""
    import scipy.linalg as sl
    rho = arb(rho)
    n = A0.nrows()
    A = A0 - rho * A1
    B = A0 - 2 * rho * A1 + rho * rho * A2
    out = dict(ok=False, nu_float=None, nubar_max=None, lower_bound_k_n=None, float_bounds_all_k=None, reason="")
    if not is_pd(B):
        out["reason"] = "B not positive definite"
        return out
    Af = np.array([[float(A[i, j].mid()) for j in range(n)] for i in range(n)])
    Bf = np.array([[float(B[i, j].mid()) for j in range(n)] for i in range(n)])
    nus = sl.eigh(Af, Bf, eigvals_only=True)
    out["nu_float"] = nus.tolist()
    out["float_bounds_all_k"] = [float(rho.mid()) - float(rho.mid()) / (1 - x) for x in nus]
    # certify the largest: nubar with nubar*B - A PD  => all nu < nubar
    top = nus[-1]
    step = 1e-12
    while True:
        nub = top + abs(top) * step
        if is_pd(arb(nub) * B - A):
            break
        step *= 2
        if not step < 1e-2:
            out["reason"] = "could not certify nu_max"
            return out
    lo, hi = top - abs(top) * 1e-9, nub
    for _ in range(40):
        mid = (lo + hi) / 2
        if is_pd(arb(mid) * B - A):
            hi = mid
        else:
            lo = mid
    nub = arb(hi)
    out["nubar_max"] = hi
    if not certainly_lt(nub, 0):
        out["reason"] = "nu_max not negative"
        return out
    lb = rho - rho / (1 - nub)
    out["lower_bound_k_n"] = lb
    out["ok"] = lb.is_finite()
    return out


# ------------------------------------------------------------------ Bolza driver (SRC lb_bolza.py)
def geometry(P0, P1, P2_):
    W1 = [P0[i] - P2_[i] for i in range(3)]; W2 = [P1[i] - P2_[i] for i in range(3)]
    G11, G12, G22 = Bm(W1, W1), Bm(W1, W2), Bm(W2, W2)
    b1, b2 = Bm(W1, P2_), Bm(W2, P2_)
    Q = {(0, 0): Bm(P2_, P2_), (1, 0): 2 * b1, (0, 1): 2 * b2, (2, 0): G11, (1, 1): 2 * G12, (0, 2): G22}
    L1 = {(0, 0): b1, (1, 0): G11, (0, 1): G12}
    L2 = {(0, 0): b2, (1, 0): G12, (0, 1): G22}
    M11 = padd(pmul(L1, L1), pscale(Q, G11), -1)
    M12 = padd(pmul(L1, L2), pscale(Q, G12), -1)
    M22 = padd(pmul(L2, L2), pscale(Q, G22), -1)
    det = det3(P0, P1, P2_)
    if certainly_gt(det, 0):
        D = det
    elif certainly_lt(det, 0):
        D = -det
    else:
        raise OutOfValidityRange("degenerate element: orientation not certain")
    return Q, (M11, M12, M22), D


def lg_bound(topo, elems, ntrial, rho, v0):
    """The Lehmann-Goerisch step on assembled element data: float P2 eigenvectors (untrusted trial functions; the
    constant mode U[:, 0] is dropped), RT2 flux solve, exact Arb forms, certified nubar.  Returns
    (lehmann dict, A0, A1, A2, R, V, float eigenvalues)."""
    import scipy.sparse.linalg as sla
    Kg, Mg, Sg, Bg, Gg = assemble_float(topo, elems)
    w, U = sla.eigsh(Kg, k=ntrial + 2, M=Mg, sigma=-0.5, v0=v0)
    o = np.argsort(w); w, U = w[o], U[:, o]
    V = U[:, 1:1 + ntrial].copy()
    Cf = solve_fluxes(Sg, Bg, V)
    A0, A1, A2, R = arb_forms(topo, elems, V, Cf)
    res = lehmann(A0, A1, A2, rho)
    return res, A0, A1, A2, R, V, w


def green_ok(R) -> bool:
    """Green residual certainly negligible (conforming sigma: int sigma.grad v + div sigma v = 0 exactly)."""
    return all(certainly_lt(abs(R[i, j]), GREEN_TOL) for i in range(R.nrows()) for j in range(R.ncols()))


def _float_down(q: Fraction) -> float:
    f = float(q)
    while Fraction(f) > q:
        f = math.nextafter(f, -math.inf)
    return f


def _arb_float_down(x) -> float:
    """float certainly <= the exact arb point x."""
    f = float(x)
    while not arb(f) <= x:
        f = math.nextafter(f, -math.inf)
    return f


def _attested_lower(ball: str):
    """Largest-ish float <= the lower end of a stored ball string '[mid +/- rad]' (exact Fraction arithmetic, as
    the checker reads it); None for a non-finite ball."""
    s = ball.strip()
    try:
        if s.startswith("["):
            mid, rad = s[1:-1].split("+/-")
            q = Fraction(mid.strip()) - Fraction(rad.strip())
        else:
            q = Fraction(s)
    except (ValueError, ZeroDivisionError, OverflowError):
        return None
    return _float_down(q)


@dataclass
class LowerBoundResult:
    lam1_lower: float | None             # certificate-attested: rounded DOWN from the stored enclosure; None if not
                                         # certified
    lam1_lower_inprocess: float | None   # in-memory ball's lower end, rounded down (fidelity pins)
    certified: bool
    certificate: Certificate


def _run(n, rho, gap, sym):
    mesh = Mesh(n)
    tri_labels = mesh.lab[mesh.T]
    topo = Topo(tri_labels, mesh.Nf)
    if topo.chi != -2:
        raise ValueError(f"Euler characteristic {topo.chi} != -2")
    Xj = positions_jet(mesh, [[J(c) for c in V] for V in bolza_vertices_hyp()])
    X = [[c.v for c in x] for x in Xj]
    per = len(mesh.T) // 8
    elems = [None] * len(mesh.T)
    for ti, (i0, i1, i2) in enumerate(mesh.T):
        sec, loc = divmod(ti, per)
        if sym and sec >= 2:
            elems[ti] = elems[(sec % 2) * per + loc]; continue
        elems[ti] = element(*geometry(X[i0], X[i1], X[i2]))
    rho_q = Fraction(rho)
    inputs = {"surface": "bolza", "n": n, "sym": sym, "trial_space": "P2", "flux_space": "RT2",
              "gap": {"lower": gap.lower, "source": gap.source}}
    deps = [Dependency("D15", "Lehmann-Goerisch theorem (shifted form): for rho <= lambda~_{m+1}, B positive definite "
                              "and nu_k < 0, lambda~_{m+1-k} >= rho - rho/(1 - nu_k)",
                       "Liu-Plum arXiv:2512.23182 Thm 4.1/4.3 (A1-A5, rho > 0); Behnke-Mertins-Plum-Wieners 2000 "
                       "Thm 2.1 (via Vejchodsky arXiv:1801.02561 §2)"),
            Dependency("gap", f"lambda_4 >= {gap.lower}, hence rho = {rho_q} <= lambda~_5 = lambda_4 + 1 (the only "
                              "check of rho: the computation cannot detect a wrong rho)", gap.source)]
    data = {"rho": str(rho_q), "shift": 1, "k": 3}
    claim_lhs = "lambda_1"
    finite = all(element_finite(el) for el in elems)
    if not finite:
        data.update({"A0": None, "A1": None, "A2": None, "lam1_lower_enclosure": "nan",
                     "diagnostics": {"error": "non-finite element data"}})
        cert = Certificate(kind="lg_lower", inputs=inputs, claims=[Claim(claim_lhs, ">=", "nan", False)],
                           data=data, dependencies=deps, producer=make_producer())
        return None, cert
    v0 = np.random.RandomState(0).rand(topo.np2)          # fixed start: reproducible trial vectors
    res, A0, A1, A2, R, V, w = lg_bound(topo, elems, 3, arb(rho_q.numerator) / rho_q.denominator, v0)
    inputs["trial"] = trial_provenance([V[:, j] for j in range(V.shape[1])], "RandomState(0).rand(np2)")
    gok = green_ok(R)
    ok = res["ok"] and gok
    lb = res["lower_bound_k_n"] - 1 if res["lower_bound_k_n"] is not None else arb("nan")
    rhs = rational_lower(lb) if (ok and lb.is_finite()) else "nan"
    ok = ok and rhs != "nan"
    data.update({
        "A0": [[ball_str(A0[i, j]) for j in range(3)] for i in range(3)],
        "A1": [[ball_str(A1[i, j]) for j in range(3)] for i in range(3)],
        "A2": [[ball_str(A2[i, j]) for j in range(3)] for i in range(3)],
        "nubar": repr(float(res["nubar_max"])) if res["nubar_max"] is not None else None,   # plain float string
        "lam1_lower_enclosure": ball_str(lb),
        "diagnostics": {
            "reason": res["reason"] or ("" if gok else "Green residual not negligible"),
            "float_fe_eigs": [json_float(x) for x in w],
            "nu_float": [json_float(x) for x in (res["nu_float"] or [])],
            "float_lower_bounds_lambda_123": [json_float(x - 1) for x in (res["float_bounds_all_k"] or [])][::-1],
            "green_residual_max": json_float(max(abs(R[i, j]).upper() for i in range(3) for j in range(3))),
            "n_elements": len(mesh.T), "np2": topo.np2, "nrt": topo.nrt,
        },
    })
    cert = Certificate(kind="lg_lower", inputs=inputs, claims=[Claim(claim_lhs, ">=", rhs, bool(ok))],
                       data=data, dependencies=deps, producer=make_producer())
    return (_arb_float_down(lb.lower()) if ok else None), cert


def certify_lam1_lower(n=48, rho=5, gap: GapInput = None, sym=False, strict=False) -> LowerBoundResult:
    """Certified lower bound lambda_1(Bolza) >= lam1_lower by Lehmann-Goerisch (module docstring).

    gap: required lambda_4 lower bound (GapInput; e.g. 4.0 from the trace-formula count certificates lp_A / lp_B).
    The hypothesis rho <= lambda~_5 = lambda_4 + 1 is checked against it (ValueError otherwise) and recorded as a
    dependency.  It is the only check of rho: the computation itself cannot detect a wrong rho (see the module
    docstring).  sym=True reuses sectors 0-1 for sectors 2-7 (SRC SYM=1; the published n = 48 run used sym=False).
    lam1_lower is certificate-attested (rounded down from the stored enclosure the checker re-derives).
    Precision is set internally (128 bits).  ~4 min for n = 48."""
    if gap is None:
        raise ValueError("gap input (lambda_4 lower bound) is required")
    rho_q = Fraction(rho)
    if not rho_q > 0:
        raise ValueError("rho must be > 0 (Liu-Plum A5)")
    if Fraction(gap.lower) + 1 < rho_q:
        raise ValueError(f"Lehmann-Goerisch hypothesis rho <= lambda~_5 = lambda_4 + 1 fails: rho = {rho}, "
                         f"lambda_4 >= {gap.lower}")
    with prec(PREC):
        inproc, cert = _run(n, rho, gap, sym)
    ok = cert.certified
    if strict and not ok:
        raise CertificationFailed("could not certify the Lehmann-Goerisch lower bound for lambda_1")
    attested = _attested_lower(cert.data["lam1_lower_enclosure"]) if ok else None
    return LowerBoundResult(lam1_lower=attested, lam1_lower_inprocess=inproc if ok else None, certified=ok,
                            certificate=cert)
