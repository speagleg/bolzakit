"""Conforming fan mesh of the glued Bolza octagon and certified (Arb / jet / dual) node positions.

Ported from the Bolza lambda1 research code: octagon.Mesh (topology and float positions; the P1-FEM assembly is not
ported), cert_rr.{bolza_vertices_hyp,positions_arb}, cert_rr_p2.positions_jet and
cert_first_order.{vertices_dual,positions_dual}.  pairing_residual is new (Arb version of octagon.check_pairings'
side_map test).
"""
import numpy as np
from flint import arb, arb_mat

from ._hyp import Bm, hyp_interp, mink, tohyp, tok, hdist
from .jets import J, D


def _hyp_interp_float(A, B, s):
    """Points at fraction s (array) of hyperbolic arclength from A to B (float hyperboloid points)."""
    s = np.atleast_1d(s)[:, None]
    L = hdist(A, B)
    if L < 1e-15:
        return np.repeat(A[None], len(s), 0)
    return (np.sinh((1 - s) * L) * A + np.sinh(s * L) * B) / np.sinh(L)


class Mesh:
    """Fixed-topology fan mesh of an octagon; positions recomputed from vertices."""

    def __init__(self, n):
        self.n = n
        idx = {}
        cnt = [0]

        def new(key):
            idx[key] = cnt[0]; cnt[0] += 1
            return idx[key]
        new(("C",))
        for k in range(8):
            for i in range(1, n + 1):
                new(("S", k, i))           # spoke k, i/n from centre to V_k
        for k in range(8):
            for i in range(2, n + 1):
                for j in range(1, i):
                    new(("T", k, i, j))    # triangle k interior / side nodes (row n = side k)
        self.N = cnt[0]
        self.idx = idx

        def node(k, i, j):
            if i == 0:
                return idx[("C",)]
            if j == 0:
                return idx[("S", k, i)]
            if j == i:
                return idx[("S", (k + 1) % 8, i)]
            return idx[("T", k, i, j)]
        self.node = node
        tris = []
        for k in range(8):
            for i in range(n):
                for j in range(i + 1):
                    tris.append((node(k, i, j), node(k, i + 1, j), node(k, i + 1, j + 1)))
                    if j < i:
                        tris.append((node(k, i, j), node(k, i + 1, j + 1), node(k, i, j + 1)))
        self.T = np.array(tris)
        # identification: side k+4 node j <-> side k node n-j
        parent = list(range(self.N))

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]; x = parent[x]
            return x
        for k in range(4):
            for j in range(n + 1):
                a = node(k + 4, n, j)
                b = node(k, n, n - j)
                parent[find(a)] = find(b)
        roots = sorted(set(find(i) for i in range(self.N)))
        rid = {r: i for i, r in enumerate(roots)}
        self.lab = np.array([rid[find(i)] for i in range(self.N)])
        self.Nf = len(roots)

    def positions(self, V, center=None):
        n = self.n
        H = tohyp(V)
        if center is None:
            c = H.sum(0); c = c / np.sqrt(mink(c, c))
        else:
            c = tohyp(center)
        X = np.zeros((self.N, 3))
        X[self.idx[("C",)]] = c
        spokes = []
        for k in range(8):
            sp_k = _hyp_interp_float(c, H[k], np.arange(n + 1) / n)
            sp_k[-1] = H[k]; sp_k[0] = c
            spokes.append(sp_k)
            for i in range(1, n + 1):
                X[self.idx[("S", k, i)]] = sp_k[i]
        for k in range(8):
            for i in range(2, n + 1):
                P, Q = spokes[k][i], spokes[(k + 1) % 8][i]
                pts = _hyp_interp_float(P, Q, np.arange(1, i) / i)
                for j in range(1, i):
                    X[self.idx[("T", k, i, j)]] = pts[j - 1]
        return tok(X)

# ------------------------------------------------------------------ certified node positions
def bolza_vertices_hyp():
    R = (3 + 2 * arb(2).sqrt()).acosh()
    ch, sh = R.cosh(), R.sinh()
    out = []
    for k in range(8):
        th = arb.pi() * k / 4
        out.append([ch, sh * th.cos(), sh * th.sin()])
    return out


def positions_arb(mesh, Vh, center=None):
    """Arb hyperboloid positions of all mesh nodes, same exact formulas as octagon.Mesh.positions."""
    n = mesh.n
    if center is None:
        c = [sum((V[i] for V in Vh), arb(0)) for i in range(3)]
        nc = Bm(c, c).sqrt()
        c = [x / nc for x in c]
    else:
        c = center
    X = [None] * mesh.N
    X[mesh.idx[("C",)]] = c
    spokes = []
    for k in range(8):
        sp_k = [c] + [hyp_interp(c, Vh[k], arb(i) / n) for i in range(1, n)] + [Vh[k]]
        spokes.append(sp_k)
        for i in range(1, n + 1):
            X[mesh.idx[("S", k, i)]] = sp_k[i]
    for k in range(8):
        for i in range(2, n + 1):
            P, Q = spokes[k][i], spokes[(k + 1) % 8][i]
            for j in range(1, i):
                X[mesh.idx[("T", k, i, j)]] = hyp_interp(P, Q, arb(j) / i)
    return X


def positions_jet(mesh, Vh):
    """jets of the hyperboloid positions of all mesh nodes; same exact formulas as octagon.Mesh.positions."""
    n = mesh.n
    c = [sum((V[i] for V in Vh), J(0)) for i in range(3)]
    nc = Bm(c, c).sqrt(); c = [x / nc for x in c]
    X = [None] * mesh.N
    X[mesh.idx[("C",)]] = c
    spokes = []
    for k in range(8):
        sp_k = [c] + [hyp_interp(c, Vh[k], arb(i) / n) for i in range(1, n)] + [Vh[k]]
        spokes.append(sp_k)
        for i in range(1, n + 1):
            X[mesh.idx[("S", k, i)]] = sp_k[i]
    for k in range(8):
        for i in range(2, n + 1):
            P, Q = spokes[k][i], spokes[(k + 1) % 8][i]
            for j in range(1, i):
                X[mesh.idx[("T", k, i, j)]] = hyp_interp(P, Q, arb(j) / i)
    return X


def vertices_dual():
    """Bolza vertices with derivative along delta: R_k = R + t (-1)^k (sigma'(0) = 0)."""
    R = (3 + 2 * arb(2).sqrt()).acosh()
    out = []
    for k in range(8):
        th = arb.pi() * k / 4
        Rk = D(R, arb((-1) ** k))
        ch, sh = Rk.cosh(), Rk.sinh()
        out.append([ch, sh * th.cos(), sh * th.sin()])
    return out


def positions_dual(mesh, Vh):
    n = mesh.n
    c = [sum((V[i] for V in Vh), D(0)) for i in range(3)]
    nc = Bm(c, c).sqrt(); c = [x / nc for x in c]
    X = [None] * mesh.N
    X[mesh.idx[("C",)]] = c
    spokes = []
    for k in range(8):
        sp_k = [c] + [hyp_interp(c, Vh[k], D(arb(i) / n)) for i in range(1, n)] + [Vh[k]]
        spokes.append(sp_k)
        for i in range(1, n + 1):
            X[mesh.idx[("S", k, i)]] = sp_k[i]
    for k in range(8):
        for i in range(2, n + 1):
            P, Q = spokes[k][i], spokes[(k + 1) % 8][i]
            for j in range(1, i):
                X[mesh.idx[("T", k, i, j)]] = hyp_interp(P, Q, D(arb(j) / i))
    return X




# ------------------------------------------------------------------ side pairings in Arb
def _frame_arb(P, Q):
    """Arb version of octagon.frame: Lorentz-orthonormal frame [P, unit tangent toward Q, normal], det > 0."""
    c = Bm(P, Q)
    e = [Q[i] - c * P[i] for i in range(3)]
    ne = (-Bm(e, e)).sqrt()
    e = [x / ne for x in e]
    # Lorentz-orthogonal complement of span{P, e}: n = G (P x e), G = diag(1,-1,-1)
    cr = [P[1] * e[2] - P[2] * e[1], P[2] * e[0] - P[0] * e[2], P[0] * e[1] - P[1] * e[0]]
    n = [cr[0], -cr[1], -cr[2]]
    nn = (-Bm(n, n)).sqrt()
    n = [x / nn for x in n]
    d = P[0] * (e[1] * n[2] - e[2] * n[1]) - e[0] * (P[1] * n[2] - P[2] * n[1]) + n[0] * (P[1] * e[2] - P[2] * e[1])
    if d.upper() < 0:
        n = [-x for x in n]
    elif not d.lower() > 0:
        raise ValueError("frame orientation is not certain at this precision")
    return arb_mat([[P[r], e[r], n[r]] for r in range(3)])


def pairings_arb(Vh):
    """Arb version of octagon.pairings: g_k (k = 0..3) maps side k+4 onto side k, V_{k+4} -> V_{k+1}, V_{k+5} -> V_k."""
    out = []
    for k in range(4):
        F1 = _frame_arb(Vh[k + 4], Vh[(k + 5) % 8])
        F2 = _frame_arb(Vh[k + 1], Vh[k])
        out.append(F2 * F1.inv())
    return out


def pairing_residual(mesh, X):
    """Upper bound (arb) of the max over paired side nodes (side k+4 node j <-> side k node n-j) of the hyperbolic
    distance between the node of side k+4 mapped by g_k and its partner.  X: Arb hyperboloid positions."""
    n = mesh.n
    Vh = [X[mesh.idx[("S", k, n)]] for k in range(8)]  # spoke tips are the vertices
    g = pairings_arb(Vh)
    worst = arb(0)
    for k in range(4):
        for j in range(n + 1):
            a = X[_side_node(mesh, k + 4, j)]
            b = X[_side_node(mesh, k, n - j)]
            ga = [sum((g[k][r, s] * a[s] for s in range(3)), arb(0)) for r in range(3)]
            u = [ga[r] - b[r] for r in range(3)]
            w = (-Bm(u, u)).upper()  # -B(u,u) = 4 sinh^2(d/2) on the hyperboloid
            d = arb(0) if w <= 0 else 2 * (w.sqrt() / 2).asinh()
            d = arb(d.upper())
            if d > worst:
                worst = d
    return worst


def _side_node(mesh, k, j):
    return mesh.node(k % 8, mesh.n, j)
