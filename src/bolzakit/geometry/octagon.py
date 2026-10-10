"""Floating-point hyperbolic geometry for the Bolza octagon (Poincare disk unless stated).

Everything here is float64 geometry used for *reduction*, *lengths* and *integer-valued counts that are
rounded*; the group elements themselves are exact (group.py).
"""
from __future__ import annotations

import cmath
import math

import numpy as np

from ..exact import group as G
from ..provenance import provenance

ZETA16 = [cmath.exp(1j * (2 * k + 1) * math.pi / 8) for k in range(8)]
VERTICES = [G.VERTEX_RADIUS * w for w in ZETA16]                       # angles (2k+1) pi/8
MIDPOINTS = [G.MIDPOINT_RADIUS * cmath.exp(1j * k * math.pi / 4) for k in range(8)]  # angles k pi/4
GEN_MATS = [G.to_matrix(g) for g in G.GENS]
GEN_INV_MATS = [G.to_matrix(G.inv(g)) for g in G.GENS]


@provenance("NUMERICAL", "closed-form hyperbolic geometry in float64 (original research code, geometry.py)", "float64 rounding")
def mobius(M, z):
    return (M[0, 0] * z + M[0, 1]) / (M[1, 0] * z + M[1, 1])


@provenance("NUMERICAL", "closed-form hyperbolic geometry in float64 (original research code, geometry.py)", "float64 rounding")
def dist(z, w):
    """Hyperbolic distance in the Poincare disk (curvature -1)."""
    num = abs(z - w) ** 2
    den = (1 - abs(z) ** 2) * (1 - abs(w) ** 2)
    return math.acosh(1 + 2 * num / den)


@provenance("NUMERICAL", "closed-form hyperbolic geometry in float64 (original research code, geometry.py)", "float64 rounding")
def cosh_dist(z, w):
    return 1 + 2 * abs(z - w) ** 2 / ((1 - abs(z) ** 2) * (1 - abs(w) ** 2))


@provenance("NUMERICAL", "closed-form hyperbolic geometry in float64 (original research code, geometry.py)", "float64 rounding")
def poincare_to_klein(z):
    return 2 * z / (1 + abs(z) ** 2)


@provenance("NUMERICAL", "closed-form hyperbolic geometry in float64 (original research code, geometry.py)", "float64 rounding")
def klein_to_poincare(k):
    r2 = abs(k) ** 2
    return k / (1 + math.sqrt(max(0.0, 1 - r2)))


GEN_IMAGES_OF_0 = [mobius(M, 0) for M in GEN_MATS]    # g_k(0), the 8 neighbour centres


@provenance("NUMERICAL", "closed-form hyperbolic geometry in float64 (original research code, geometry.py)", "float64 rounding")
def in_F(z, tol=1e-12):
    """z in the closed Dirichlet octagon: d(z,0) <= d(z, g_k 0) for all k (with tolerance)."""
    c0 = cosh_dist(z, 0)
    return all(c0 <= cosh_dist(z, w) * (1 + tol) for w in GEN_IMAGES_OF_0)


@provenance("NUMERICAL", "closed-form hyperbolic geometry in float64 (original research code, geometry.py)", "float64 rounding")
def reduce_point(z, max_steps=10_000):
    """Return (z', word, M): z' = M(z) lies in the octagon F, M is the matrix of the total map, and `word` lists the
    generator indices used (applied left-to-right as z -> g_{w1}^{...}).  Implemented as: while some g_k 0 is closer
    to z than 0, replace z by g_k^{-1} z.  Raises ValueError unless z lies in the open unit disk (|z| < 1)."""
    if not abs(z) < 1:
        raise ValueError(f"reduce_point needs a point of the open unit disk, got |z| = {abs(z)!r}")
    M = np.eye(2, dtype=complex)
    word = []
    for _ in range(max_steps):
        c0 = cosh_dist(z, 0)
        best, kbest = c0, -1
        for k, w in enumerate(GEN_IMAGES_OF_0):
            c = cosh_dist(z, w)
            if c < best * (1 - 1e-13):
                best, kbest = c, k
        if kbest < 0:
            return z, word, M
        z = mobius(GEN_INV_MATS[kbest], z)
        M = GEN_INV_MATS[kbest] @ M
        word.append(G.inv_index(kbest))
    raise RuntimeError("reduction did not terminate")


# ---------------------------------------------------------------- axes and clipping

def fixed_points(M):
    """Boundary fixed points (repelling, attracting) of a hyperbolic SU(1,1)-type matrix acting on the disk."""
    a, b, c, d = M[0, 0], M[0, 1], M[1, 0], M[1, 1]
    # c z^2 + (d - a) z - b = 0
    disc = cmath.sqrt((d - a) ** 2 + 4 * b * c)
    z1 = (a - d + disc) / (2 * c)
    z2 = (a - d - disc) / (2 * c)
    # attracting fixed point: |derivative| < 1 ; derivative = 1/(c z + d)^2
    if abs(c * z1 + d) > abs(c * z2 + d):
        return z2, z1
    return z1, z2


# Klein-model octagon: vertices are Klein images of the Poincare vertices.
KLEIN_VERTS = np.array([poincare_to_klein(v) for v in VERTICES])
_KV = np.stack([KLEIN_VERTS.real, KLEIN_VERTS.imag], axis=1)
_EDGE_N = []  # inward normals and offsets: n . x >= c for x in F
for k in range(8):
    p, q = _KV[k], _KV[(k + 1) % 8]
    e = q - p
    n = np.array([-e[1], e[0]])
    if n @ (-p) < 0:
        n = -n
    _EDGE_N.append((n, n @ p))


def clip_chord(p, q):
    """Intersect the Klein chord p->q (points on the unit circle, as complex) with the closed octagon.
    Returns (t0, t1) parameters on p + t(q-p) or None."""
    P = np.array([p.real, p.imag])
    D = np.array([(q - p).real, (q - p).imag])
    t0, t1 = 0.0, 1.0
    for n, c in _EDGE_N:
        num = c - n @ P
        den = n @ D
        if abs(den) < 1e-15:
            if num > 1e-13:
                return None
            continue
        t = num / den
        if den > 0:
            t0 = max(t0, t)
        else:
            t1 = min(t1, t)
        if t0 > t1:
            return None
    return t0, t1


def klein_dist(x, y):
    """Distance between two points of the Klein disk (complex)."""
    num = 1 - (x.real * y.real + x.imag * y.imag)
    den = math.sqrt((1 - abs(x) ** 2) * (1 - abs(y) ** 2))
    return math.acosh(max(1.0, num / den))


def axis_length_in_F(M, edge_tol=1e-9):
    """(length of axis(M) inside the closed octagon F, weight) where weight = 1/2 if the clipped segment
    lies on a side of F (it is then shared with the neighbouring tile), else 1.  Returns 0 if disjoint."""
    zr, za = fixed_points(M)
    p, q = zr, za  # boundary points are their own Klein images
    r = clip_chord(p, q)
    if r is None:
        return 0.0
    t0, t1 = r
    if t1 - t0 < 1e-14:
        return 0.0
    x = p + t0 * (q - p)
    y = p + t1 * (q - p)
    L = klein_dist(x, y)
    w = 1.0
    for n, c in _EDGE_N:
        if abs(n @ np.array([x.real, x.imag]) - c) < edge_tol and abs(n @ np.array([y.real, y.imag]) - c) < edge_tol:
            w = 0.5
            break
    return L * w


@provenance("NUMERICAL", "closed-form hyperbolic geometry in float64 (original research code, geometry.py)", "float64 rounding")
def disk_to_upper(z: complex) -> complex:
    """Cayley map from the Poincaré disk to the upper half-plane, w = i(1+z)/(1-z)."""
    return 1j * (1 + z) / (1 - z)


@provenance("NUMERICAL", "closed-form hyperbolic geometry in float64 (original research code, geometry.py)", "float64 rounding")
def upper_to_disk(w: complex) -> complex:
    """Inverse Cayley map, z = (w - i)/(w + i)."""
    return (w - 1j) / (w + 1j)
