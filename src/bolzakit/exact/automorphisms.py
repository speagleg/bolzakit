"""Automorphisms of the Bolza surface: N(Gamma)/Gamma with N(Gamma) = Delta(2,3,8).

Exact part.  Elements of PSU(1,1) normalising Gamma are represented projectively as pairs (A, B) in
Z[zeta8]^2 standing for the matrix [[A, mu B], [mu conj B, conj A]] with positive determinant
d = |A|^2 - mu^2|B|^2 in Z[sqrt2] (scalar multiples act identically on the disk).  Conjugation
  n gamma n^{-1} = n gamma adj(n) / d,   adj(A,B) = (conj A, -B),
is computed exactly; divisibility by d in Z[zeta8] and membership of the result in Gamma are checked
exactly (membership by lookup in an exact enumeration of a displacement ball).

  rho  : z -> zeta8 z  (rotation by pi/4 about the octagon centre 0).  It is not of the (A,B) shape with
         entries in Q(zeta8) (it needs e^{i pi/8}), so its conjugation action is coded directly:
         rho (A,B) rho^{-1} = (A, zeta8 B).
  sigma: rotation by pi about Q = point at distance r_in/2 from 0 in direction 0 (half-way to the midpoint
         of side 0); (A,B) = (2i, -i(2 - sqrt2)), det 4(2 - sqrt2), trace 0, sigma(0) = midpoint M_0.
         (Q is the right-angle vertex of a (2,3,8) triangle with its pi/8 angle at 0:
          cosh d(0,Q) = cos(pi/3)/sin(pi/8) = 1/(2 sin(pi/8)), and cosh(2 d(0,Q)) = 1 + sqrt2 = cosh r_in.)
  kappa: z -> conj(z) (orientation reversing), kappa (A,B) kappa = (conj A, conj B).

Float part.  The quotient N/Gamma is computed with Mobius maps normalised so that n(p) lies in F for a
generic base point p; cosets are keyed by the rounded reduced image of p.
"""
from __future__ import annotations

import cmath
import itertools
import math

import numpy as np

from . import group as G
from . import zeta8 as Z
from ..geometry import octagon as Geo

SIGMA = ((0, 0, 2, 0), Z.neg(Z.mul(Z.I, Z.sub((2, 0, 0, 0), Z.SQRT2))))


def pdet(n):
    return G.det(n)


def conj_by(n, g):
    """n g n^{-1} exactly; raises if not integral."""
    A, B = n
    adj = (Z.conj(A), Z.neg(B))
    M = G.mul(G.mul(n, g), adj)
    d = pdet(n)  # (p, q) in Z[sqrt2]
    # divide by d = p + q sqrt2: multiply by Galois conjugate p - q sqrt2, divide by integer norm p^2 - 2 q^2
    dbar = Z.sub((d[0], 0, 0, 0), Z.scal(d[1], Z.SQRT2))
    N = d[0] * d[0] - 2 * d[1] * d[1]
    out = []
    for X in M:
        Y = Z.mul(X, dbar)
        if any(c % N for c in Y):
            raise ValueError("conjugate not integral")
        out.append(tuple(c // N for c in Y))
    return tuple(out)


def conj_rho(g, k=1):
    A, B = g
    return (A, Z.mul(Z.zeta_pow(k), B))


def conj_kappa(g):
    return (Z.conj(g[0]), Z.conj(g[1]))


def normalises_exactly(R=8.0):  # images have displacement <= 4.75
    """Check that rho, sigma, kappa map each generator g_k into Gamma (exact membership in the
    displacement ball of radius R).  Returns dict of word-images."""
    X = G.ball_by_displacement(R)
    keys = {tuple(int(v) for v in r) for r in X}
    res = {}
    for name, f in [("rho", conj_rho), ("sigma", lambda g: conj_by(SIGMA, g)), ("kappa", conj_kappa)]:
        imgs = [f(g) for g in G.GENS]
        # Gamma here is the lift in SU(1,1); conjugation may change the lift by a sign character, so
        # membership is tested in PSU(1,1): h or -h.  `signs` records which.
        signs = []
        for h in imgs:
            kp = tuple(list(h[0]) + list(h[1]))
            km = tuple(-v for v in kp)
            signs.append(1 if kp in keys else (-1 if km in keys else 0))
        ok = all(s_ != 0 for s_ in signs)
        dets = [G.det(h) for h in imgs]
        res[name] = dict(in_Gamma=ok, signs=signs, dets_one=all(d == (1, 0) for d in dets), images=imgs)
    return res


# ---------------------------------------------------------------- float quotient N/Gamma

def _mat_sigma():
    A, B = SIGMA
    a = Z.to_complex(A)
    b = G.MU * Z.to_complex(B)
    M = np.array([[a, b], [np.conj(b), np.conj(a)]], dtype=complex)
    return M / cmath.sqrt(np.linalg.det(M))


M_RHO = np.array([[cmath.exp(1j * math.pi / 8), 0], [0, cmath.exp(-1j * math.pi / 8)]], dtype=complex)
M_SIGMA = _mat_sigma()
P_GENERIC = 0.0731 + 0.0417j


def _normalise(M):
    """Left-multiply M by an element of Gamma so that M(P_GENERIC) lies in F."""
    z = Geo.mobius(M, P_GENERIC)
    _, _, R = Geo.reduce_point(z)
    M = R @ M
    return M / cmath.sqrt(np.linalg.det(M))


def _key(M, nd=7):
    z = Geo.mobius(M, P_GENERIC)
    return (round(z.real, nd), round(z.imag, nd))


def quotient_group(gens=(M_RHO, M_SIGMA)):
    """Enumerate N/Gamma for N = <Gamma, gens>.  Returns (elements as normalised matrices, mult table)."""
    I = np.eye(2, dtype=complex)
    elems = [_normalise(I)]
    index = {_key(elems[0]): 0}
    frontier = [0]
    while frontier:
        nf = []
        for i in frontier:
            for g in gens:
                M = _normalise(elems[i] @ g)
                k = _key(M)
                if k not in index:
                    index[k] = len(elems)
                    elems.append(M)
                    nf.append(index[k])
        frontier = nf
        if len(elems) > 10_000:
            raise RuntimeError("quotient too large -- generators do not normalise Gamma?")
    n = len(elems)
    table = np.zeros((n, n), dtype=int)
    for i, j in itertools.product(range(n), range(n)):
        table[i, j] = index[_key(_normalise(elems[i] @ elems[j]))]
    return elems, table


def element_orders(table):
    n = len(table)
    orders = []
    for i in range(n):
        k, x = 1, i
        while x != 0:
            x = table[x, i]
            k += 1
        orders.append(k)
    return orders


def center(table):
    n = len(table)
    return [i for i in range(n) if all(table[i, j] == table[j, i] for j in range(n))]


# ---------------------------------------------------------------- GL(2,3)

def gl23():
    els = []
    for a, b, c, d in itertools.product(range(3), repeat=4):
        if (a * d - b * c) % 3:
            els.append((a, b, c, d))
    return els


def _m3(x, y):
    a, b, c, d = x
    e, f, g, h = y
    return ((a * e + b * g) % 3, (a * f + b * h) % 3, (c * e + d * g) % 3, (c * f + d * h) % 3)


def find_isomorphism_to_gl23(table, gen_idx):
    """Search images of the generators in GL(2,3) such that the induced map (defined along a BFS spanning
    tree) is a well-defined bijective homomorphism.  Returns dict elem_index -> GL(2,3) element, or None."""
    n = len(table)
    G3 = gl23()
    # spanning tree words
    parent = {0: None}
    order = [0]
    for i in order:
        for gi, g in enumerate(gen_idx):
            j = table[i, g]
            if j not in parent:
                parent[j] = (i, gi)
                order.append(j)
    for imgs in itertools.product(G3, repeat=len(gen_idx)):
        phi = {0: (1, 0, 0, 1)}
        for j in order[1:]:
            i, gi = parent[j]
            phi[j] = _m3(phi[i], imgs[gi])
        if len(set(phi.values())) != n:
            continue
        if all(phi[table[i, j]] == _m3(phi[i], phi[j]) for i in range(n) for j in range(n)):
            return phi
    return None


# ---------------------------------------------------------------- special points

def weierstrass_points_from_involution(R=9.0):
    """Fixed points on Sigma of the hyperelliptic involution J = rho^4 (z -> -z).  A point z is fixed iff
    -gamma(z) = z for some gamma, i.e. z is the fixed point of the order-2 elliptic z -> -gamma(z), whose
    trace is -2 Im A (so we need Im A = 0, decided exactly).  Returns the distinct reduced fixed points."""
    X = G.ball_by_displacement(R)
    pts = []
    for r in X:
        g = G.elt(r)
        A = g[0]
        if Z.conj(A) != A:
            continue
        M = G.to_matrix(g)
        Mj = np.array([[-M[0, 0], -M[0, 1]], [M[1, 0], M[1, 1]]])  # z -> -gamma(z)
        a, b, c, d = Mj[0, 0], Mj[0, 1], Mj[1, 0], Mj[1, 1]
        if abs(c) < 1e-14:
            z = 0j if abs(b) < 1e-14 else None
        else:
            disc = cmath.sqrt((d - a) ** 2 + 4 * b * c)
            z1 = (a - d + disc) / (2 * c)
            z2 = (a - d - disc) / (2 * c)
            z = z1 if abs(z1) < 1 else z2
        if z is None or abs(z) >= 1:
            continue
        zr, _, _ = Geo.reduce_point(z)
        pts.append(zr)
    return dedupe_mod_gamma(pts)


def same_point_mod_gamma(z, w, tol=1e-7):
    """True if w = gamma z for some gamma (checked against images of z under a small ball)."""
    X = _SMALL_BALL()
    for r in X:
        if abs(G.act(G.elt(r), z) - w) < tol:
            return True
    return False


_cache = {}


def _SMALL_BALL():
    if "b" not in _cache:
        _cache["b"] = G.ball_by_displacement(2 * G.R_C + 0.5)
    return _cache["b"]


def dedupe_mod_gamma(pts):
    out = []
    for z in pts:
        if not any(same_point_mod_gamma(z, w) for w in out):
            out.append(z)
    return out


def stabiliser(elems, table, z):
    """Indices of cosets n with n(z) = z mod Gamma."""
    return [i for i, M in enumerate(elems) if same_point_mod_gamma(Geo.mobius(M, z), z)]


def is_cyclic(table, idx):
    orders = element_orders(table)
    return any(orders[i] == len(idx) for i in idx)


def orbit_on_surface(elems, z):
    """Distinct points of Sigma in the Aut+ orbit of [z] (reduced representatives)."""
    pts = [Geo.reduce_point(Geo.mobius(M, z))[0] for M in elems]
    return dedupe_mod_gamma(pts)


# Cone points of Delta(2,3,8) relative to the octagon centre (orbifold points of Sigma/Aut+):
#   order 8: z = 0;   order 2: Q = tanh(r_in/4) (on the real axis, d(0,Q) = r_in/2);
#   order 3: P with cosh d(0,P) = cot(pi/8)/sqrt3, in direction pi/8 (angle between 0Q and 0P is pi/8).
Q_ORDER2 = math.tanh(G.R_IN / 4)
P_ORDER3 = math.tanh(math.acosh((1 + G.SQ2) / math.sqrt(3)) / 2) * cmath.exp(1j * math.pi / 8)


# ---------------------------------------------------------------- bolzakit API

import functools  # noqa: E402
from ..provenance import provenance  # noqa: E402


class Automorphisms:
    """Aut+(Bolza) = N/Gamma, with N = <Gamma, rho, sigma> (the (2,3,8) triangle group); |Aut+| = 48 = GL(2,3)."""

    def __init__(self):
        self.elements, self.table = quotient_group()
        self.order = len(self.elements)
        self.full_order = 2 * self.order          # with the reflection kappa: z -> conj(z)
        self._phi = None
        self._W = None

    def element_order_counts(self) -> dict[int, int]:
        out: dict[int, int] = {}
        for o in element_orders(self.table):
            out[int(o)] = out.get(int(o), 0) + 1
        return dict(sorted(out.items()))

    def isomorphism_to_GL23(self) -> dict[int, tuple]:
        """index -> 2x2 matrix over F_3 as a tuple of tuples ((a, b), (c, d))."""
        if self._phi is None:
            idx = {_key(M): i for i, M in enumerate(self.elements)}
            gi = [idx[_key(_normalise(M_RHO))], idx[_key(_normalise(M_SIGMA))]]
            phi = find_isomorphism_to_gl23(self.table, gi)
            if phi is None:
                raise RuntimeError("no isomorphism to GL(2,3) found")
            self._phi = {int(i): ((a, b), (c, d)) for i, (a, b, c, d) in phi.items()}
        return dict(self._phi)

    def weierstrass_points(self) -> list[complex]:
        if self._W is None:
            self._W = weierstrass_points_from_involution()
        return list(self._W)

    def stabiliser(self, z: complex) -> list[int]:
        return stabiliser(self.elements, self.table, z)

    def chart_to(self, p: complex) -> np.ndarray:
        """A normaliser element M (2x2 matrix) with M(p) = 0 mod Gamma (M includes the Gamma-reduction)."""
        for M in self.elements:
            w = Geo.mobius(M, p)
            if same_point_mod_gamma(w, 0j):
                _, _, R = Geo.reduce_point(w)
                return R @ M
        raise ValueError(f"{p} is not in the Aut+ orbit of the octagon centre")


@functools.lru_cache(maxsize=1)
def _automorphisms() -> Automorphisms:
    return Automorphisms()


@provenance("PROVED", "generators of the normaliser conjugate Γ into Γ exactly; quotient enumerated (item 02); "
            "isomorphism to GL(2,3) by explicit map")
def automorphisms() -> Automorphisms:
    return _automorphisms()
