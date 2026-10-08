"""Γ(p) ⊂ SL(2,Z): the golden geodesic (axis of Q = [[1,1],[1,0]]) is a systole of Γ(p)\\H iff p is a prime Lucas
number L_{2^j} (3, 7, 47, 2207, ...).  Item 05; Lean: MathPursuits.Item05_GoldenGeodesic.golden_systolic_iff.
This concerns congruence covers of the modular surface, not the Bolza group."""
from __future__ import annotations
import math
from ..provenance import provenance


@provenance("EXACT", "integer recurrence L_0=2, L_1=1, L_{n+1}=L_n+L_{n-1}")
def lucas(n: int) -> int:
    a, b = 2, 1
    for _ in range(n):
        a, b = b, a + b
    return a


def _is_prime(n: int) -> bool:
    """Trial division by odd numbers: O(sqrt(n)) time. Fine for the small primes used here, not for large n."""
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
    """True iff the golden geodesic is a systole of Gamma(p)\\H, for p prime.

    Theorem (item 05): let k be the rank of apparition of p (least k >= 1 with p | F_k) and
    n(p) = k if k is even, else 2k. Then the geodesic is a systole iff L_{n(p)} = p^2 - 2,
    equivalently p is a prime Lucas number L_{2^j} (3, 7, 47, 2207, ...). Non-primes and p = 2 give False.
    """
    if p == 2 or not _is_prime(p):
        return False
    phi = (1 + 5 ** 0.5) / 2
    n = _n_of_p(p, int(2 * math.log(p, phi)) + 4)
    return n is not None and lucas(n) == p * p - 2
