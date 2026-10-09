"""P1 (hyperboloid-barycentric) element forms in Arb ball arithmetic.

Ported from the Bolza lambda1 research code (cert_rr.py): exact element integrals of the hyperboloid-barycentric P1
trial space on a geodesic triangle, with a Lagrange-remainder-controlled series in Q^{-s}; assembly into float
matrices (untrusted, for the eigenvector guess) and exact Arb projections onto a float trial space.
Callers must set the working precision (``_arb.prec``); this module never changes global state.
"""
import math
import time

import scipy.sparse as sp
from math import comb
from flint import arb

from ._hyp import Bm, det3

K_SERIES = 3


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
    assert emax < 0.5, "triangle too large for the series"
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

