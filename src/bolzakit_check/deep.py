"""--deep: NON-rigorous consistency check of sampled rr_upper element matrices.

Each sampled element is a geodesic triangle with hyperboloid vertices P0, P1, P2 (B(x, y) = x0 y0 - x1 y1 - x2 y2,
B(P, P) = 1).  Hyperboloid-barycentric chart: phi(a, b) = w / sqrt(B(w, w)), w = a P0 + b P1 + (1 - a - b) P2 on the
simplex a, b >= 0, a + b <= 1.  The pulled-back hyperbolic metric is g_ij = -B(d_i phi, d_j phi), and
    K_pq = int grad s_p^T g^-1 grad s_q sqrt(det g) da db,     M_pq = int s_p s_q sqrt(det g) da db
for the Lagrange shapes s_p in the barycentric coordinates (l0, l1, l2) = (a, b, 1 - a - b):
    P1: l0, l1, l2;   P2: 2 l_i^2 - l_i (i = 0, 1, 2), then 4 l0 l1, 4 l1 l2, 4 l2 l0.
Here the integrals are evaluated by tensor Gauss-Legendre quadrature (collapsed square, mpmath, 30 digits) and
compared with the stored balls to 1e-8 relative.  This is a consistency check only, not a proof.
"""
from __future__ import annotations

import mpmath as mp

from ._cert import CheckError, ball

TOL = 1e-8
NODES = 16


def _gl(n):
    """Gauss-Legendre nodes / weights on [0, 1] (mpmath, current precision)."""
    out = []
    for k in range(1, n + 1):
        x = mp.cos(mp.pi * (k - mp.mpf(1) / 4) / (n + mp.mpf(1) / 2))
        for _ in range(100):
            p0, p1 = mp.mpf(1), x
            for j in range(2, n + 1):
                p0, p1 = p1, ((2 * j - 1) * x * p1 - (j - 1) * p0) / j
            dp = n * (x * p1 - p0) / (x * x - 1)
            dx = p1 / dp
            x -= dx
            if abs(dx) < mp.mpf(10) ** (-mp.mp.dps - 5):
                break
        out.append(((1 - x) / 2, 1 / ((1 - x * x) * dp * dp)))   # weight 2/((1-x^2) P'^2), halved for [0, 1]
    return out


def _shapes(deg, a, b):
    l = (a, b, 1 - a - b)
    dl = ((1, 0), (0, 1), (-1, -1))
    if deg == 1:
        return list(l), [list(d) for d in dl]
    s, g = [], []
    for i in range(3):
        s.append(2 * l[i] ** 2 - l[i])
        g.append([(4 * l[i] - 1) * dl[i][c] for c in range(2)])
    for i, j in ((0, 1), (1, 2), (2, 0)):
        s.append(4 * l[i] * l[j])
        g.append([4 * (dl[i][c] * l[j] + l[i] * dl[j][c]) for c in range(2)])
    return s, g


def element(P, deg):
    """(K, M) of one element by quadrature (lists of mpf)."""
    B = lambda x, y: x[0] * y[0] - x[1] * y[1] - x[2] * y[2]  # noqa: E731
    nb = 3 if deg == 1 else 6
    K = [[mp.mpf(0)] * nb for _ in range(nb)]
    M = [[mp.mpf(0)] * nb for _ in range(nb)]
    rule = _gl(NODES)
    for xi, wx in rule:
        for et, we in rule:
            a, b = xi, (1 - xi) * et                 # collapsed square -> simplex, Jacobian (1 - xi)
            wt = wx * we * (1 - xi)
            w = [a * P[0][i] + b * P[1][i] + (1 - a - b) * P[2][i] for i in range(3)]
            Q = B(w, w)
            sQ = mp.sqrt(Q)
            dw = [[P[0][i] - P[2][i] for i in range(3)], [P[1][i] - P[2][i] for i in range(3)]]
            # d_i phi = d_i w / sqrt(Q) - w B(d_i w, w) / Q^{3/2}
            dphi = [[dw[c][i] / sQ - w[i] * B(dw[c], w) / (Q * sQ) for i in range(3)] for c in range(2)]
            g = [[-B(dphi[c], dphi[d]) for d in range(2)] for c in range(2)]
            det = g[0][0] * g[1][1] - g[0][1] * g[1][0]
            vol = mp.sqrt(det)
            gi = [[g[1][1] / det, -g[0][1] / det], [-g[1][0] / det, g[0][0] / det]]
            s, gr = _shapes(deg, a, b)
            for p in range(nb):
                hp = [gi[c][0] * gr[p][0] + gi[c][1] * gr[p][1] for c in range(2)]
                for q in range(p, nb):
                    K[p][q] += wt * vol * (gr[q][0] * hp[0] + gr[q][1] * hp[1])
                    M[p][q] += wt * vol * s[p] * s[q]
    for p in range(nb):
        for q in range(p):
            K[p][q], M[p][q] = K[q][p], M[q][p]
    return K, M


def _close(stored, comp):
    scale = max(abs(mp.mpf(ball(x).mid().str(30, radius=False))) for r in stored for x in r)
    for i, r in enumerate(stored):
        for j, x in enumerate(r):
            b = ball(x)
            mid = mp.mpf(b.mid().str(30, radius=False))
            rad = mp.mpf(b.rad().str(5, radius=False))
            if abs(comp[i][j] - mid) > TOL * scale + rad:
                return False
    return True


def check_elements(data, degree):
    """(ok, message) for data["element_sample"] of an rr_upper certificate."""
    es = data.get("element_sample")
    if not es or not es.get("items"):
        return True, "(deep: consistency only; no element_sample)"
    with mp.workdps(30):
        for it in es["items"]:
            P = [[mp.mpf(ball(c).mid().str(30, radius=False)) for c in V] for V in it["vertices"]]
            if len(P) != 3:
                raise CheckError("element_sample vertices must be 3 points")
            K, M = element(P, degree)
            if not (_close(it["K"], K) and _close(it["M"], M)):
                return False, f"(deep: consistency only) element {it.get('index')} disagrees with quadrature"
    return True, f"(deep: consistency only; {len(es['items'])} elements agree to {TOL:g})"
