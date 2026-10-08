"""Exact arithmetic in Z[zeta8] and Z[sqrt2].

An element of Z[zeta8] (zeta = exp(i*pi/4), zeta^4 = -1) is a 4-tuple of Python ints
(a0, a1, a2, a3) meaning a0 + a1*zeta + a2*zeta^2 + a3*zeta^3.  Python ints are unbounded,
so every operation here is exact.

An element of Z[sqrt2] is a pair (p, q) meaning p + q*sqrt2.

Vectorised helpers (int64 numpy, used for large enumerations) express multiplication by a
FIXED element y as a 4x4 integer matrix acting on row vectors: x*y == x @ mat(y).
"""
from __future__ import annotations

import math
from fractions import Fraction

import numpy as np

Z8 = tuple  # alias for readability

ZERO = (0, 0, 0, 0)
ONE = (1, 0, 0, 0)
ZETA = (0, 1, 0, 0)
I = (0, 0, 1, 0)                   # zeta^2 = i
SQRT2 = (0, 1, 0, -1)              # zeta + zeta^-1 = zeta - zeta^3
MU2 = (2, 2, 0, -2)                # mu^2 = 2 + 2 sqrt2  (mu = sqrt(2+2sqrt2) = 2 sinh(l_sys/2)/sqrt... see group.py)


def add(x, y):
    return (x[0] + y[0], x[1] + y[1], x[2] + y[2], x[3] + y[3])


def sub(x, y):
    return (x[0] - y[0], x[1] - y[1], x[2] - y[2], x[3] - y[3])


def neg(x):
    return (-x[0], -x[1], -x[2], -x[3])


def mul(x, y):
    """Negacyclic convolution (zeta^4 = -1)."""
    a0, a1, a2, a3 = x
    b0, b1, b2, b3 = y
    return (
        a0 * b0 - a1 * b3 - a2 * b2 - a3 * b1,
        a0 * b1 + a1 * b0 - a2 * b3 - a3 * b2,
        a0 * b2 + a1 * b1 + a2 * b0 - a3 * b3,
        a0 * b3 + a1 * b2 + a2 * b1 + a3 * b0,
    )


def conj(x):
    """Complex conjugation: zeta^k -> zeta^-k = -zeta^(4-k)."""
    return (x[0], -x[3], -x[2], -x[1])


def scal(n, x):
    return (n * x[0], n * x[1], n * x[2], n * x[3])


def zeta_pow(k):
    k %= 8
    s = 1 if k < 4 else -1
    v = [0, 0, 0, 0]
    v[k % 4] = s
    return tuple(v)


def to_complex(x):
    z = complex(math.cos(math.pi / 4), math.sin(math.pi / 4))
    return x[0] + x[1] * z + x[2] * z * z + x[3] * z ** 3


def is_real(x):
    """x == conj(x)."""
    return x == conj(x)


def real_to_zsqrt2(x):
    """For a REAL x in Z[zeta8] return (p, q) with x = p + q*sqrt2 (p, q may be half-integers
    in general; returned as Fractions only if needed)."""
    assert is_real(x), x
    # real => a2 = 0 and a3 = -a1, x = a0 + a1*(zeta - zeta^3) = a0 + a1*sqrt2
    return (x[0], x[1])


def re2(x):
    """2*Re(x) = x + conj(x) as an element (p, q) of Z[sqrt2] (always integral)."""
    s = add(x, conj(x))
    return real_to_zsqrt2(s)


def abs2(x):
    """|x|^2 = x*conj(x) in Z[sqrt2] as (p, q)."""
    return real_to_zsqrt2(mul(x, conj(x)))


# ---------------------------------------------------------------- Z[sqrt2] (also Q(sqrt2))

def zs_sign(p, q):
    """Exact sign of p + q*sqrt2 for integers or Fractions p, q."""
    if p >= 0 and q >= 0:
        return 0 if (p == 0 and q == 0) else 1
    if p <= 0 and q <= 0:
        return -1
    # opposite signs: compare p^2 with 2 q^2
    d = p * p - 2 * q * q
    if d == 0:
        return 0  # impossible for nonzero rationals (sqrt2 irrational) but keep exact
    return (1 if p > 0 else -1) if d > 0 else (1 if q > 0 else -1)


def zs_cmp(x, y):
    """Compare x=(p,q), y=(p',q') in Q(sqrt2) exactly: -1, 0, 1."""
    return zs_sign(x[0] - y[0], x[1] - y[1])


def zs_float(x):
    return float(x[0]) + float(x[1]) * math.sqrt(2)


def zs_str(x):
    p, q = x
    if isinstance(p, Fraction) or isinstance(q, Fraction):
        return f"{p}+{q}*sqrt2"
    return f"{p}{'+' if q >= 0 else '-'}{abs(q)}sqrt2"


# ---------------------------------------------------------------- vectorised (int64)

def mat_right(y):
    """4x4 int matrix M with (x*y) == x @ M for row vectors x (coefficient tuples)."""
    rows = []
    for k in range(4):
        e = [0, 0, 0, 0]
        e[k] = 1
        rows.append(mul(tuple(e), y))
    return np.array(rows, dtype=np.int64)


CONJ_MAT = np.array([[1, 0, 0, 0], [0, 0, 0, -1], [0, 0, -1, 0], [0, -1, 0, 0]], dtype=np.int64)
"""conj(x) == x @ CONJ_MAT."""
