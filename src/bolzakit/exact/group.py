"""Exact model of the Bolza Fuchsian group Gamma in SU(1,1).

Geometry (README C1).  Regular hyperbolic octagon F centred at 0 in the Poincare disk with
interior angles pi/4:
  * circumradius R_c:  cosh R_c = cot^2(pi/8) = 3 + 2 sqrt2;  vertex Poincare radius tanh(R_c/2) = 2^(-1/4);
  * inradius r_in:      cosh r_in = cot(pi/8) = 1 + sqrt2;   side-midpoint Poincare radius sqrt(sqrt2 - 1);
  * CONVENTION: side midpoints at angles k*pi/4, vertices at angles (2k+1)*pi/8 (k = 0..7).
Side k+4 is paired with side k by the hyperbolic translation g_k along the axis through 0 at
angle k*pi/4 with translation length l_sys = 2 r_in = 2 arccosh(1+sqrt2):

    g_k = [[ cosh(l/2),              sinh(l/2) e^{i k pi/4} ],
           [ sinh(l/2) e^{-i k pi/4}, cosh(l/2)             ]],
    cosh(l/2) = 1 + sqrt2,  sinh(l/2) = sqrt(2 + 2 sqrt2) =: mu.

Exact representation.  mu is irrational over Q(zeta8), but mu^2 = 2 + 2 sqrt2 is in Z[zeta8], so
every product of generators has the form

    gamma = [[ A, mu*B ], [ mu*conj(B), conj(A) ]],   A, B in Z[zeta8],  |A|^2 - mu^2 |B|^2 = 1,

and we store gamma as the 8-tuple of ints (A, B).  Products:
    (A1,B1)(A2,B2) = (A1 A2 + mu^2 B1 conj(B2),  A1 B2 + B1 conj(A2)),     (A,B)^{-1} = (conj A, -B).
Traces tr(gamma) = 2 Re A lie in Z[sqrt2]; half-traces Re A are proved/checked to lie in Z[sqrt2].
NOTE: the invariant quaternion algebra of Gamma is a division algebra over Q(sqrt2) (README C1), so Gamma
cannot be conjugated into SL(2, Q(sqrt2)); some quadratic extension (here Q(zeta8, mu)) is unavoidable.
(-I is not in Gamma: checked in tests by enumeration.)

Generators: GENS[k] = g_k for k = 0..7 with g_{k+4} = g_k^{-1}.
Relation (convention of checks/bgrp.py):  g0 g1^-1 g2 g3^-1 g0^-1 g1 g2^-1 g3 = I, i.e. with
g_k^{-1} = g_{k+4}:  RELATOR = (0, 5, 2, 7, 4, 1, 6, 3).
"""
from __future__ import annotations

import math

import numpy as np

from . import zeta8 as Z

# ---------------------------------------------------------------- constants (floats, for geometry)
SQ2 = math.sqrt(2.0)
COSH_HALF_LSYS = 1 + SQ2
L_SYS = 2 * math.acosh(1 + SQ2)               # 3.0571418390...
COSH_RC = 3 + 2 * SQ2                         # circumradius
R_C = math.acosh(COSH_RC)                     # 2.4485...
R_IN = math.acosh(1 + SQ2)                    # 1.5286...
VERTEX_RADIUS = 2 ** -0.25                    # Poincare radius of the vertices
MIDPOINT_RADIUS = math.sqrt(SQ2 - 1)          # Poincare radius of side midpoints
MU = math.sqrt(2 + 2 * SQ2)
AREA = 4 * math.pi
GENUS = 2

Elt = tuple  # (A, B) with A, B 4-tuples

IDENTITY = (Z.ONE, Z.ZERO)
MINUS_IDENTITY = (Z.neg(Z.ONE), Z.ZERO)
A_GEN = (1, 1, 0, -1)  # 1 + sqrt2


def gen(k: int) -> Elt:
    """g_k, k mod 8; g_{k+4} = g_k^{-1}."""
    k %= 8
    if k < 4:
        return (A_GEN, Z.zeta_pow(k))
    return (A_GEN, Z.neg(Z.zeta_pow(k - 4)))


GENS = [gen(k) for k in range(8)]
RELATOR = (0, 5, 2, 7, 4, 1, 6, 3)


def inv_index(k: int) -> int:
    return (k + 4) % 8


def mul(g: Elt, h: Elt) -> Elt:
    A1, B1 = g
    A2, B2 = h
    A = Z.add(Z.mul(A1, A2), Z.mul(Z.MU2, Z.mul(B1, Z.conj(B2))))
    B = Z.add(Z.mul(A1, B2), Z.mul(B1, Z.conj(A2)))
    return (A, B)


def inv(g: Elt) -> Elt:
    A, B = g
    return (Z.conj(A), Z.neg(B))


def word(ws) -> Elt:
    g = IDENTITY
    for k in ws:
        g = mul(g, GENS[k])
    return g


def det(g: Elt):
    """|A|^2 - mu^2 |B|^2 as (p, q) in Z[sqrt2]; equals (1, 0) for group elements."""
    A, B = g
    a = Z.mul(A, Z.conj(A))
    b = Z.mul(Z.MU2, Z.mul(B, Z.conj(B)))
    return Z.real_to_zsqrt2(Z.sub(a, b))


def trace(g: Elt):
    """tr(gamma) = 2 Re A in Z[sqrt2] as (p, q)."""
    return Z.re2(g[0])


def half_trace(g: Elt):
    """Re A as (p, q) with p + q sqrt2; raises if not in Z[sqrt2]."""
    p, q = trace(g)
    if p % 2 or q % 2:
        raise ValueError(f"half-trace not in Z[sqrt2]: trace = {p}+{q}sqrt2")
    return (p // 2, q // 2)


def abs_half_trace(g: Elt):
    p, q = half_trace(g)
    return (p, q) if Z.zs_sign(p, q) >= 0 else (-p, -q)


def classify(g: Elt) -> str:
    """'identity' | 'elliptic' | 'parabolic' | 'hyperbolic', decided exactly from |tr| vs 2."""
    if g == IDENTITY or g == MINUS_IDENTITY:
        return "identity"
    p, q = trace(g)
    if Z.zs_sign(p, q) < 0:
        p, q = -p, -q
    s = Z.zs_cmp((p, q), (2, 0))
    return {-1: "elliptic", 0: "parabolic", 1: "hyperbolic"}[s]


def rows_classify_nonhyperbolic(X: np.ndarray) -> np.ndarray:
    """Indices of rows that are NOT hyperbolic (|tr| <= 2), decided exactly: float pre-filter with a
    wide margin, exact Z[sqrt2] comparison for every row within the margin."""
    P, Q = rows_trace(X)
    t = np.abs(P.astype(float) + Q.astype(float) * _S2)
    cand = np.nonzero(t < 2 + 1e-6)[0]
    out = []
    for i in cand:
        if classify(elt(X[i])) != "hyperbolic":
            out.append(int(i))
    return np.array(out, dtype=np.int64)


def cosh_disp0(g: Elt):
    """cosh d(0, gamma 0) = 2|A|^2 - 1 in Z[sqrt2] (exact)."""
    p, q = Z.abs2(g[0])
    return (2 * p - 1, 2 * q)


# ---------------------------------------------------------------- float conversions

def to_matrix(g: Elt) -> np.ndarray:
    A, B = g
    a = Z.to_complex(A)
    b = MU * Z.to_complex(B)
    return np.array([[a, b], [np.conj(b), np.conj(a)]], dtype=complex)


def act(g: Elt, z: complex) -> complex:
    M = to_matrix(g)
    return (M[0, 0] * z + M[0, 1]) / (M[1, 0] * z + M[1, 1])


def translation_length(g: Elt) -> float:
    t = abs(Z.zs_float(trace(g))) / 2
    return 2 * math.acosh(t) if t > 1 else 0.0


def disp0(g: Elt) -> float:
    return math.acosh(max(1.0, Z.zs_float(cosh_disp0(g))))


# ---------------------------------------------------------------- vectorised enumeration
# Elements as int64 rows of length 8: (A0..A3, B0..B3).  Right/left multiplication by a fixed element is
# Z-linear in these coordinates (conj is linear), hence an 8x8 integer matrix.

def right_mult_matrix(h: Elt) -> np.ndarray:
    """R with row(g*h) = row(g) @ R."""
    A2, B2 = h
    # A = A1*A2 + B1*(mu^2 conj B2);  B = A1*B2 + B1*conj(A2)
    R = np.zeros((8, 8), dtype=np.int64)
    R[:4, :4] = Z.mat_right(A2)
    R[4:, :4] = Z.mat_right(Z.mul(Z.MU2, Z.conj(B2)))
    R[:4, 4:] = Z.mat_right(B2)
    R[4:, 4:] = Z.mat_right(Z.conj(A2))
    return R


def left_mult_matrix(h: Elt) -> np.ndarray:
    """L with row(h*g) = row(g) @ L."""
    A1, B1 = h
    # A = A1*A2 + mu^2 B1 conj(B2);  B = A1*B2 + B1*conj(A2)
    L = np.zeros((8, 8), dtype=np.int64)
    L[:4, :4] = Z.mat_right(A1)                                 # A2 -> A2*A1
    L[4:, :4] = Z.CONJ_MAT @ Z.mat_right(Z.mul(Z.MU2, B1))     # B2 -> conj(B2)*mu^2 B1
    L[4:, 4:] = Z.mat_right(A1)                                 # B2 -> B2*A1
    L[:4, 4:] = Z.CONJ_MAT @ Z.mat_right(B1)                   # A2 -> conj(A2)*B1
    return L


def row(g: Elt) -> np.ndarray:
    return np.array(list(g[0]) + list(g[1]), dtype=np.int64)


def elt(r) -> Elt:
    r = [int(x) for x in r]
    return (tuple(r[:4]), tuple(r[4:]))


R_GENS = [right_mult_matrix(g) for g in GENS]

_S2 = math.sqrt(2.0)


def rows_abs2A(X: np.ndarray):
    """|A|^2 = p + q sqrt2 for each row, returned as int arrays (p, q).  Uses A*conj(A):
    with A = a0 + a1 z + a2 z^2 + a3 z^3, |A|^2 = (a0^2+a1^2+a2^2+a3^2) + (a0 a1 + a1 a2 + a2 a3 - a3 a0) sqrt2."""
    a0, a1, a2, a3 = (X[:, i] for i in range(4))
    p = a0 * a0 + a1 * a1 + a2 * a2 + a3 * a3
    q = a0 * a1 + a1 * a2 + a2 * a3 - a3 * a0
    return p, q


def rows_cosh_disp0(X: np.ndarray) -> np.ndarray:
    p, q = rows_abs2A(X)
    return (2 * p - 1).astype(float) + 2 * q.astype(float) * _S2


def rows_trace(X: np.ndarray):
    """trace = 2 Re A = 2 a0 + (a1 - a3) sqrt2 as int arrays (P, Q)."""
    return 2 * X[:, 0], X[:, 1] - X[:, 3]


def _unique_rows(X: np.ndarray) -> np.ndarray:
    if len(X) == 0:
        return X
    Xc = np.ascontiguousarray(X)
    v = Xc.view(np.dtype((np.void, Xc.dtype.itemsize * Xc.shape[1])))
    _, idx = np.unique(v, return_index=True)
    return Xc[np.sort(idx)]


def _setdiff_rows(X: np.ndarray, Y: np.ndarray) -> np.ndarray:
    """Rows of X not in Y."""
    if len(X) == 0 or len(Y) == 0:
        return X
    Xc = np.ascontiguousarray(X)
    Yc = np.ascontiguousarray(Y)
    dt = np.dtype((np.void, Xc.dtype.itemsize * Xc.shape[1]))
    m = np.isin(Xc.view(dt).ravel(), Yc.view(dt).ravel())
    return Xc[~m]


def spheres_by_word_length(L: int):
    """Exact spheres S_0..S_L of the Cayley graph w.r.t. the 8 generators (as int64 row arrays).
    Uses bipartiteness (relator has even length) : S_{n+1} = (S_n * gens) minus S_{n-1}.  The S_n are
    disjoint as long as the parity map is well defined, which the tests verify independently."""
    S = [np.array([row(IDENTITY)])]
    prev = np.zeros((0, 8), dtype=np.int64)
    for n in range(L):
        cand = np.concatenate([S[-1] @ R for R in R_GENS])
        cand = _unique_rows(cand)
        cand = _setdiff_rows(cand, prev)
        cand = _setdiff_rows(cand, S[-1])  # empty if bipartite; kept for safety
        prev = S[-1]
        S.append(cand)
    return S


def ball_by_displacement(Rmax: float, tol: float = 1e-9) -> np.ndarray:
    """All gamma with d(0, gamma 0) <= Rmax (exact int rows), including the identity.

    Completeness: F is the Dirichlet domain of Gamma centred at 0 (regular octagon, side k is the bisector
    of 0 and g_k 0).  For gamma != 1, 0 is not in gamma F, so 0 violates one of the 8 half-plane conditions
    defining gamma F, i.e. d(0, gamma g_k 0) < d(0, gamma 0) for some k.  Hence every element of the ball
    is reached from 1 by right multiplications through elements of the ball, so BFS with cutoff Rmax
    is complete.  The cutoff uses floats with a tolerance; boundary cases (equality) are kept."""
    ch = math.cosh(Rmax) * (1 + tol)
    seen = np.array([row(IDENTITY)])
    front = seen
    while len(front):
        cand = np.concatenate([front @ R for R in R_GENS])
        cand = cand[rows_cosh_disp0(cand) <= ch]
        cand = _unique_rows(cand)
        cand = _setdiff_rows(cand, seen)
        seen = np.concatenate([seen, cand])
        front = cand
    return seen


def max_abs_coeff(X: np.ndarray) -> int:
    return int(np.abs(X).max()) if len(X) else 0


INT64_SAFE = 2 ** 28  # products of two coordinates (and short sums) stay far below 2^63


from ..provenance import provenance  # noqa: E402


@provenance("EXACT", "2·arccosh(1+√2): half-trace of the generators is 1+√2 (exact, item 02)")
def systole() -> float:
    """Length of the shortest closed geodesic of the Bolza surface, 2·arccosh(1+√2) = 3.05714183896..."""
    return L_SYS


@provenance("EXACT", "exact enumeration by word length in Z[zeta8] (item 02)")
def elements(max_word_length: int) -> list:
    """All distinct group elements of word length <= max_word_length, as Elt objects."""
    return [elt(r) for S in spheres_by_word_length(max_word_length) for r in S]
