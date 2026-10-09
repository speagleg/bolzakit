"""P1 and P2 (hyperboloid-barycentric) element forms in Arb ball arithmetic.

Ported from the Bolza lambda1 research code (cert_rr.py): exact element integrals of the hyperboloid-barycentric P1
trial space on a geodesic triangle, with a Lagrange-remainder-controlled series in Q^{-s}; assembly into float
matrices (untrusted, for the eigenvector guess) and exact Arb projections onto a float trial space.
Callers must set the working precision (``_arb.prec``); this module never changes global state.
"""
import math
import time
from math import comb

import numpy as np
import scipy.sparse as sp
from flint import arb, arb_mat

from ..errors import OutOfValidityRange
from ._arb import certainly_gt, certainly_lt
from ._hyp import Bm, det3
from .jets import J

K_SERIES = 3      # P1 series order
K_SERIES_P2 = 6   # P2 series order (fixed; no environment override)


# ----- polynomial helpers: dict {(i,j): arb} in variables (a, b) -----
def pmul(p, q):
    r = {}
    for (i, j), c in p.items():
        for (k, l), d in q.items():
            key = (i + k, j + l)
            r[key] = r[key] + c * d if key in r else c * d
    return r


def padd(p, q, s=1):
    r = dict(p)
    for k, v in q.items():
        r[k] = r[k] + s * v if k in r else s * v
    return r


def pscale(p, s):
    return {k: v * s for k, v in p.items()}


FACT = [math.factorial(_i) for _i in range(40)]  # exact Python ints (no import-time precision dependence)


def pint(p):
    """exact integral over the simplex {a,b >= 0, a+b <= 1}."""
    return sum((c * FACT[i] * FACT[j] / FACT[i + j + 2] for (i, j), c in p.items()), arb(0))


def psup_simplex(p):
    """crude upper bound of |p| on the simplex: a, b in [0,1]."""
    return sum((abs(c) for c in p.values()), arb(0))


def psup_centered(p):
    """upper bound of |p| on the simplex using centered variables a = 1/3 + u, b = 1/3 + w, |u|,|w| <= 2/3."""
    # re-expand around the centroid
    third = arb(1) / 3
    sh = {}
    for (i, j), c in p.items():
        # (u + 1/3)^i (w + 1/3)^j
        for k in range(i + 1):
            for l in range(j + 1):
                key = (k, l)
                v = c * comb(i, k) * comb(j, l) * third ** (i - k) * third ** (j - l)
                sh[key] = sh[key] + v if key in sh else v
    tt = arb(2) / 3
    return sum((abs(c) * tt ** (k + l) for (k, l), c in sh.items()), arb(0))


def binom_neg(s, k):
    """binom(-s, k) for arb s."""
    r = arb(1)
    for i in range(k):
        r = r * (-s - i) / (i + 1)
    return r


def element_p1(P0, P1, P2):
    """3x3 stiffness and mass matrices (arb) of the hyperboloid-barycentric P1 element on triangle (P0,P1,P2)."""
    W1 = [P0[i] - P2[i] for i in range(3)]
    W2 = [P1[i] - P2[i] for i in range(3)]
    G11, G12, G22 = Bm(W1, W1), Bm(W1, W2), Bm(W2, W2)
    b1, b2 = Bm(W1, P2), Bm(W2, P2)
    Q = {(0, 0): Bm(P2, P2), (1, 0): 2 * b1, (0, 1): 2 * b2, (2, 0): G11, (1, 1): 2 * G12, (0, 2): G22}
    L1 = {(0, 0): b1, (1, 0): G11, (0, 1): G12}
    L2 = {(0, 0): b2, (1, 0): G12, (0, 1): G22}
    D = abs(det3(P0, P1, P2))
    # N = adj(L L^T - Q G) = [[L2^2 - Q G22, -L1 L2 + Q G12], [ ., L1^2 - Q G11]]
    N11 = padd(pmul(L2, L2), pscale(Q, G22), -1)
    N12 = padd(pscale(Q, G12), pmul(L1, L2), -1)
    N22 = padd(pmul(L1, L1), pscale(Q, G11), -1)
    third = arb(1) / 3
    Qc = sum((c * third ** (i + j) for (i, j), c in Q.items()), arb(0))
    eps = pscale(Q, 1 / Qc); eps = padd(eps, {(0, 0): arb(1)}, -1)
    emax = psup_centered(eps)
    if not certainly_lt(emax, 0.5):
        raise OutOfValidityRange("element too large for the series bound (emax >= 0.5)")
    epows = [{(0, 0): arb(1)}]
    for k in range(K_SERIES):
        epows.append(pmul(epows[-1], eps))

    def integ(p, s):
        tot = arb(0)
        for k in range(K_SERIES + 1):
            tot += binom_neg(s, k) * pint(pmul(p, epows[k]))
        rem = abs(binom_neg(s, K_SERIES + 1)) * emax ** (K_SERIES + 1) * (1 - emax) ** (-s - K_SERIES - 1)
        err = psup_simplex(p) * rem / 2
        return (tot + arb(0, float(err.upper()) * (1 + 1e-9) + 1e-300)) * Qc ** (-s)
    half, three_half = arb(1) / 2, arb(3) / 2
    # gradients of phi0 = a, phi1 = b, phi2 = 1-a-b in (a,b)
    grads = [(1, 0), (0, 1), (-1, -1)]
    Nint = {}
    Nint[(0, 0)] = integ(N11, half); Nint[(0, 1)] = integ(N12, half); Nint[(1, 1)] = integ(N22, half)
    Nm = [[Nint[(0, 0)], Nint[(0, 1)]], [Nint[(0, 1)], Nint[(1, 1)]]]
    Ke = [[None] * 3 for _ in range(3)]
    for p in range(3):
        for q in range(p, 3):
            gp, gq = grads[p], grads[q]
            v = sum((gp[i] * gq[j] * Nm[i][j] for i in range(2) for j in range(2)), arb(0)) / D
            Ke[p][q] = Ke[q][p] = v
    phis = [{(1, 0): arb(1)}, {(0, 1): arb(1)}, {(0, 0): arb(1), (1, 0): arb(-1), (0, 1): arb(-1)}]
    Me = [[None] * 3 for _ in range(3)]
    for p in range(3):
        for q in range(p, 3):
            Me[p][q] = Me[q][p] = D * integ(pmul(phis[p], phis[q]), three_half)
    return Ke, Me


def assemble_p1(mesh, X):
    t0 = time.time()
    T = mesh.T; lab = mesh.lab
    elems = []
    rows, cols, kv, mv = [], [], [], []
    for t, (i0, i1, i2) in enumerate(T):
        Ke, Me = element_p1(X[i0], X[i1], X[i2])
        elems.append((Ke, Me))
        L = [lab[i0], lab[i1], lab[i2]]
        for p in range(3):
            for q in range(3):
                rows.append(L[p]); cols.append(L[q])
                kv.append(float(Ke[p][q].mid())); mv.append(float(Me[p][q].mid()))
    Kf = sp.csr_matrix((kv, (rows, cols)), shape=(mesh.Nf, mesh.Nf))
    Mf = sp.csr_matrix((mv, (rows, cols)), shape=(mesh.Nf, mesh.Nf))
    return elems, Kf, Mf, time.time() - t0


def project_p1(mesh, elems, Xvecs):
    """exact (arb) projections  Kp[a][b] = sum_T x_a|T^T K_T x_b|T  for float (dyadic) vectors Xvecs[a]."""
    m = len(Xvecs)
    Kp = [[arb(0)] * m for _ in range(m)]; Mp = [[arb(0)] * m for _ in range(m)]
    lab = mesh.lab
    xa = [[arb(float(v)) for v in x] for x in Xvecs]  # exact conversion of doubles
    for t, (i0, i1, i2) in enumerate(mesh.T):
        Ke, Me = elems[t]
        L = (lab[i0], lab[i1], lab[i2])
        loc = [[xa[a][L[p]] for p in range(3)] for a in range(m)]
        for a in range(m):
            Ka = [sum((Ke[p][q] * loc[a][q] for q in range(3)), arb(0)) for p in range(3)]
            Ma = [sum((Me[p][q] * loc[a][q] for q in range(3)), arb(0)) for p in range(3)]
            for b in range(a, m):
                kab = sum((loc[b][p] * Ka[p] for p in range(3)), arb(0))
                mab = sum((loc[b][p] * Ma[p] for p in range(3)), arb(0))
                Kp[a][b] += kab; Mp[a][b] += mab
    for a in range(m):
        for b in range(a):
            Kp[a][b] = Kp[b][a]; Mp[a][b] = Mp[b][a]
    return Kp, Mp



# ======================================================================= P2 (quadratic) elements
# ------------------------------------------------------------------ P2 dof map
class P2Space:
    def __init__(self, mesh):
        self.mesh = mesh
        lab = mesh.lab
        ekey = {}
        loc = []
        cnt = {}
        for (i0, i1, i2) in mesh.T:
            L = (int(lab[i0]), int(lab[i1]), int(lab[i2]))
            ids = []
            for (u, v) in ((0, 1), (1, 2), (2, 0)):
                k = (min(L[u], L[v]), max(L[u], L[v]))
                if k not in ekey:
                    ekey[k] = len(ekey)
                cnt[k] = cnt.get(k, 0) + 1
                ids.append(mesh.Nf + ekey[k])
            loc.append(list(L) + ids)
        self.loc = np.array(loc)
        self.Ne = len(ekey)
        self.ndof = mesh.Nf + self.Ne
        chi = mesh.Nf - self.Ne + len(mesh.T)
        # structural mesh checks: explicit raises (a bare assert is stripped under python -O)
        if chi != -2:
            raise ValueError(f"Euler characteristic {chi}")
        if not all(v == 2 for v in cnt.values()):
            raise ValueError("edge not shared by exactly 2 triangles")
        if any(L[0] == L[1] or L[1] == L[2] or L[0] == L[2] for L in loc):
            raise ValueError("degenerate triangle labels")


def mons(d):
    return [(i, t - i) for t in range(d + 1) for i in range(t, -1, -1)]


_LAM = [{(1, 0): 1}, {(0, 1): 1}, {(0, 0): 1, (1, 0): -1, (0, 1): -1}]


def _ipmul(p, q):
    r = {}
    for (i, j), c in p.items():
        for (k, l), d in q.items():
            r[(i + k, j + l)] = r.get((i + k, j + l), 0) + c * d
    return {k: v for k, v in r.items() if v != 0}


def _ipadd(p, q, s=1):
    r = dict(p)
    for k, v in q.items():
        r[k] = r.get(k, 0) + s * v
    return {k: v for k, v in r.items() if v != 0}


SHAPES = []
for _i in range(3):
    SHAPES.append(_ipadd(_ipmul({(0, 0): 2}, _ipmul(_LAM[_i], _LAM[_i])), _LAM[_i], -1))
for (_u, _v) in ((0, 1), (1, 2), (2, 0)):
    SHAPES.append(_ipmul({(0, 0): 4}, _ipmul(_LAM[_u], _LAM[_v])))


def _ider(p, var):
    r = {}
    for (i, j), c in p.items():
        if var == 0 and i > 0: r[(i - 1, j)] = r.get((i - 1, j), 0) + c * i
        if var == 1 and j > 0: r[(i, j - 1)] = r.get((i, j - 1), 0) + c * j
    return r


GRADS = [[_ider(s, 0), _ider(s, 1)] for s in SHAPES]
KL = {(0, 0): 0, (0, 1): 1, (1, 0): 1, (1, 1): 2}
# stiffness table: K_pq = sum over (kl, m) of coef * C_kl(m), C_kl(m) = int m a_kl
KTAB = {}
MTAB = {}
for _p in range(6):
    for _q in range(_p, 6):
        t = {}
        for k in range(2):
            for l in range(2):
                for m1, c1 in GRADS[_p][k].items():
                    for m2, c2 in GRADS[_q][l].items():
                        key = (KL[(k, l)], (m1[0] + m2[0], m1[1] + m2[1]))
                        t[key] = t.get(key, 0) + c1 * c2
        KTAB[(_p, _q)] = [(kk, m, c) for (kk, m), c in t.items() if c != 0]
        MTAB[(_p, _q)] = list(_ipmul(SHAPES[_p], SHAPES[_q]).items())


# ----- fixed moment tables (precision dependent, hence built lazily at the current working precision) -----
DEG_MAX = 8
MONS8 = mons(DEG_MAX)
MONS_SER = mons(2 * K_SERIES_P2)
_TABLES = {}


class _Tables:
    pass


def wmono(i, j):
    return arb(FACT[i]) * FACT[j] / FACT[i + j + 2]


def _tables():
    """WMAT, W8, S_LIST, BINOM, BOXES at the current flint precision (cached per precision)."""
    import flint
    key = flint.ctx.prec
    if key not in _TABLES:
        t = _Tables()
        t.WMAT = arb_mat([[wmono(m[0] + j[0], m[1] + j[1]) for j in MONS_SER] for m in MONS8])
        t.W8 = {m: wmono(*m) for m in MONS8}
        t.S_LIST = [arb(1) / 2, arb(3) / 2, arb(5) / 2, arb(7) / 2]
        t.BINOM = arb_mat([[binom_neg(s, k) for s in t.S_LIST] for k in range(K_SERIES_P2 + 1)])
        t.BOXES = [(arb("0.25", "0.25"), arb("0.25", "0.25")), (arb("0.75", "0.25"), arb("0.25", "0.25")),
                   (arb("0.25", "0.25"), arb("0.75", "0.25"))]
        _TABLES[key] = t
    return _TABLES[key]


def peval(p, a, b):
    return sum((c * a ** i * b ** j for (i, j), c in p.items()), arb(0))


_INF = arb("inf")


def element_p2(P0, P1, P2, order=2):
    """P2 element on the triangle with jet vertices P0, P1, P2.
    Returns K, M: lists [K0, K1, K2] of 6x6 arb lists (Taylor coefficients in p; K1 = q', 2*K2 = q''), and
    betas = (beta1, beta1r, beta2, beta2r) (upper bounds, arb) when order == 2; all +inf if a relative bound cannot
    be certified (never a pass downstream).  Raises OutOfValidityRange for an element whose orientation or series
    bound cannot be certified."""
    tb = _tables()
    W1 = [P0[i] - P2[i] for i in range(3)]; W2 = [P1[i] - P2[i] for i in range(3)]
    G11, G12, G22 = Bm(W1, W1), Bm(W1, W2), Bm(W2, W2)
    b1, b2 = Bm(W1, P2), Bm(W2, P2)
    Qj = {(0, 0): Bm(P2, P2), (1, 0): 2 * b1, (0, 1): 2 * b2, (2, 0): G11, (1, 1): 2 * G12, (0, 2): G22}
    L1 = {(0, 0): b1, (1, 0): G11, (0, 1): G12}
    L2 = {(0, 0): b2, (1, 0): G12, (0, 1): G22}
    det = det3(P0, P1, P2)
    if certainly_gt(det.v, 0):
        Dj = det
    elif certainly_lt(det.v, 0):
        Dj = -det
    else:
        raise OutOfValidityRange("degenerate element: orientation not certain")
    Nj = [padd(pmul(L2, L2), pscale(Qj, G22), -1), padd(pscale(Qj, G12), pmul(L1, L2), -1),
          padd(pmul(L1, L1), pscale(Qj, G11), -1)]           # N11, N12, N22
    comp = lambda P, a: {k: getattr(c, a) for k, c in P.items()}
    Q0, Q1, Q2 = comp(Qj, "v"), comp(Qj, "d"), comp(Qj, "e")
    N0 = [comp(P, "v") for P in Nj]
    third = arb(1) / 3
    Qc = peval(Q0, third, third)
    eps = padd(pscale(Q0, 1 / Qc), {(0, 0): arb(1)}, -1)
    emax = psup_centered(eps)
    if not certainly_lt(emax, 0.25):
        raise OutOfValidityRange("element too large for the series bound (emax >= 0.25)")
    epows = [{(0, 0): arb(1)}]
    for _ in range(K_SERIES_P2):
        epows.append(pmul(epows[-1], eps))
    E = arb_mat([[ep.get(j, arb(0)) for ep in epows] for j in MONS_SER])
    Z = tb.WMAT * E * tb.BINOM                                     # (45 x 4): series moments for s = 1/2..7/2
    mu = []
    for si, s in enumerate(tb.S_LIST):
        rem = abs(binom_neg(s, K_SERIES_P2 + 1)) * emax ** (K_SERIES_P2 + 1) * (1 - emax) ** (-s - K_SERIES_P2 - 1)
        rr = float(rem.upper()) * (1 + 1e-12)
        sc = Qc ** (-s)
        mu.append({m: (Z[r, si] + arb(0, float((tb.W8[m] * rr).upper()) * (1 + 1e-12))) * sc
                   for r, m in enumerate(MONS8)})

    def I(si, P, m):
        mm = mu[si]
        return sum((c * mm[(k[0] + m[0], k[1] + m[1])] for k, c in P.items()), arb(0))
    M2m = mons(2); M4m = mons(4)
    half, th, fh = arb(1) / 2, arb(3) / 2, arb(5) / 2
    if order == 0:
        Cv = [{m: I(0, N0[kk], m) / Dj.v for m in M2m} for kk in range(3)]
        Rv = {m: Dj.v * I(1, {(0, 0): arb(1)}, m) for m in M4m}
        K0 = [[None] * 6 for _ in range(6)]; M0 = [[None] * 6 for _ in range(6)]
        for (p, q), tab in KTAB.items():
            K0[p][q] = K0[q][p] = sum((c * Cv[kk][m] for kk, m, c in tab), arb(0))
            M0[p][q] = M0[q][p] = sum((c * Rv[m] for m, c in MTAB[(p, q)]), arb(0))
        return [K0], [M0], None
    N1 = [comp(P, "d") for P in Nj]; N2 = [comp(P, "e") for P in Nj]
    Q1sq = pmul(Q1, Q1)
    C = []
    for kk in range(3):
        n0q1 = pmul(N0[kk], Q1); n1q1 = pmul(N1[kk], Q1); n0q2 = pmul(N0[kk], Q2); n0q1q1 = pmul(N0[kk], Q1sq)
        Ck = {}
        for m in M2m:
            s0 = I(0, N0[kk], m)
            s1 = I(0, N1[kk], m) - half * I(1, n0q1, m)
            s2 = I(0, N2[kk], m) - half * I(1, n1q1, m) - half * I(1, n0q2, m) + arb(3) / 8 * I(2, n0q1q1, m)
            Ck[m] = J(s0, s1, s2) / Dj
        C.append(Ck)
    R = {}
    one = {(0, 0): arb(1)}
    for m in M4m:
        j0 = I(1, one, m)
        j1 = -th * I(2, Q1, m)
        j2 = -th * I(2, Q2, m) + arb(15) / 8 * I(3, Q1sq, m)
        R[m] = Dj * J(j0, j1, j2)
    Ks = [[[None] * 6 for _ in range(6)] for _ in range(3)]
    Ms = [[[None] * 6 for _ in range(6)] for _ in range(3)]
    for (p, q), tab in KTAB.items():
        kv = J(0)
        for kk, m, c in tab:
            kv = kv + C[kk][m] * c
        mv = J(0)
        for m, c in MTAB[(p, q)]:
            mv = mv + R[m] * c
        for o, a in enumerate(("v", "d", "e")):
            Ks[o][p][q] = Ks[o][q][p] = getattr(kv, a)
            Ms[o][p][q] = Ms[o][q][p] = getattr(mv, a)
    # ---- relative derivative bounds over boxes covering the simplex
    b1m = b1r = b2m = b2r = arb(0)
    for (A_, B_) in tb.BOXES:
        Qb = J(peval(Q0, A_, B_), peval(Q1, A_, B_), peval(Q2, A_, B_))
        Nb = [J(peval(N0[k], A_, B_), peval(N1[k], A_, B_), peval(N2[k], A_, B_)) for k in range(3)]
        f = 1 / (Dj * Qb.sqrt())
        rho = Dj / (Qb * Qb.sqrt())
        n11, n12, n22 = Nb[0].v, Nb[1].v, Nb[2].v
        dN = n11 * n22 - n12 * n12
        if not (certainly_gt(dN, 0) and certainly_gt(n11, 0)):
            # N not certainly positive definite on the box (SRC asserted): no relative bound, the element is
            # marked non-finite so that no certificate built on it can pass (fail closed, survives python -O)
            return Ks, Ms, (_INF, _INF, _INF, _INF)
        Ninv = [[n22 / dN, -n12 / dN], [-n12 / dN, n11 / dN]]
        Nm1 = [[Nb[0].d, Nb[1].d], [Nb[1].d, Nb[2].d]]
        Nm2 = [[Nb[0].e, Nb[1].e], [Nb[1].e, Nb[2].e]]
        X1 = [[sum((Ninv[i][k] * Nm1[k][j] for k in range(2)), arb(0)) for j in range(2)] for i in range(2)]
        X2 = [[sum((Ninv[i][k] * Nm2[k][j] for k in range(2)), arb(0)) for j in range(2)] for i in range(2)]
        r1 = f.d / f.v; r2 = f.e / f.v
        A1 = [[X1[i][j] + (r1 if i == j else 0) for j in range(2)] for i in range(2)]
        A2 = [[2 * (X2[i][j] + r1 * X1[i][j] + (r2 if i == j else 0)) for j in range(2)] for i in range(2)]
        fro = lambda A: arb(sum((A[i][j] ** 2 for i in range(2) for j in range(2)), arb(0)).abs_upper()).sqrt()
        v1, v2 = fro(A1).upper(), fro(A2).upper()
        w1, w2 = abs(rho.d / rho.v).upper(), abs(2 * rho.e / rho.v).upper()
        if not all(x.is_finite() for x in (v1, v2, w1, w2)):
            return Ks, Ms, (_INF, _INF, _INF, _INF)      # non-finite beta bound: fail closed (SRC asserted)
        b1m = v1 if v1 > b1m else b1m; b2m = v2 if v2 > b2m else b2m
        b1r = w1 if w1 > b1r else b1r; b2r = w2 if w2 > b2r else b2r
    return Ks, Ms, (b1m, b1r, b2m, b2r)


def compute_elements(mesh, X, order=2, sym=True, log=None):
    """element data for all triangles; with sym=True only sectors 0 and 1 are computed (pi/2-rotation
    invariance of the vertex configuration and of the mesh construction, CERT_RESULTS Sec. 2.4)."""
    T = mesh.T; per = len(T) // 8
    out = [None] * len(T)
    t0 = time.time()
    for ti, (i0, i1, i2) in enumerate(T):
        sec, loc = divmod(ti, per)
        if sym and sec >= 2:
            out[ti] = out[(sec % 2) * per + loc]
            continue
        out[ti] = element_p2(X[i0], X[i1], X[i2], order)
        if log and ti % 2000 == 0:
            log(f"  element {ti}/{len(T)} ({time.time() - t0:.0f}s)")
    return out


def assemble_float(space, elems, order_idx, which):
    """float (midpoint) sparse global matrix of Taylor coefficient order_idx; which = 0 (K) or 1 (M)."""
    loc = space.loc
    rows = np.repeat(loc, 6, axis=1).ravel(); cols = np.tile(loc, (1, 6)).ravel()
    vals = np.empty(len(loc) * 36)
    for t, e in enumerate(elems):
        Mx = e[which][order_idx]
        vals[36 * t:36 * t + 36] = [float(Mx[p][q].mid()) for p in range(6) for q in range(6)]
    return sp.csr_matrix((vals, (rows, cols)), shape=(space.ndof, space.ndof))


def to_arbmats(elems, nforms_K, nforms_M):
    """convert element matrices to arb_mat objects (shared between symmetric copies)."""
    cache = {}
    out = []
    for e in elems:
        k = id(e)
        if k not in cache:
            cache[k] = ([arb_mat(e[0][o]) for o in range(nforms_K)], [arb_mat(e[1][o]) for o in range(nforms_M)])
        out.append(cache[k])
    return out



def gram(space, amats, vecs, form_idx_K, form_idx_M):
    """exact (arb) Gram matrices  V^T K V, V^T M V  summed over elements, for float vectors vecs (list)."""
    m = len(vecs)
    Vf = np.array(vecs).T
    GK = arb_mat(m, m); GM = arb_mat(m, m)
    for t, (Km, Mm) in enumerate(amats):
        Vl = arb_mat(Vf[space.loc[t]].tolist())
        GK += Vl.transpose() * (Km[form_idx_K] * Vl)
        GM += Vl.transpose() * (Mm[form_idx_M] * Vl)
    return GK, GM

