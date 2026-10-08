"""Γ(p) ⊂ SL(2,Z): the golden geodesic (axis of Q = [[1,1],[1,0]]) is a systole of Γ(p)\\H iff p is a prime Lucas
number L_{2^j} (3, 7, 47, 2207, ...).  Item 05; Lean: MathPursuits.Item05_GoldenGeodesic.golden_systolic_iff.
This concerns congruence covers of the modular surface, not the Bolza group."""
from __future__ import annotations
import math
from ..provenance import provenance


def lucas(n: int) -> int:
    a, b = 2, 1
    for _ in range(n):
        a, b = b, a + b
    return a


def _is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    f = 3
    while f * f <= n:
        if n % f == 0:
            return False
        f += 2
    return True


def _n_of_p(p: int, cap: int):
    a, b, k = 0, 1, 0
    while k < cap:
        a, b = b, (a + b) % p
        k += 1
        if a == 0:
            n = k if k % 2 == 0 else 2 * k
            return n if n <= cap else None
    return None


@provenance("PROVED", "item 05 proof + Lean theorem golden_systolic_iff (no sorry)")
def golden_geodesic_is_systole(p: int) -> bool:
    if p == 2 or not _is_prime(p):
        return False
    phi = (1 + 5 ** 0.5) / 2
    n = _n_of_p(p, int(2 * math.log(p, phi)) + 4)
    return n is not None and lucas(n) == p * p - 2
