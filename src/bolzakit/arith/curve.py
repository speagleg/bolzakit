"""Arithmetic of the Bolza curve y² = x⁵ − x over Q (item 15c). Needs the optional extra: pip install bolzakit[pari].
The Jacobian is Q-simple (L_3 = x⁴−2x²+9 is irreducible over Q) and splits over Q(√2) as E′² with CM by Z[√−2]."""
from __future__ import annotations
from ..provenance import provenance, require


def _pari():
    return require("cypari2", "pari").Pari()


@provenance("EXACT", "PARI hyperellcharpoly")
def euler_factor(p: int) -> list[int]:
    pari = _pari()
    P = pari(f"hyperellcharpoly(Mod(1,{p})*(x^5-x))")
    return [int(c) for c in pari.Vec(P)]


def _check_good_prime(p) -> int:
    if isinstance(p, bool) or not isinstance(p, int) or p <= 0:
        raise ValueError(f"prime must be a positive int, got {p!r}")
    if p == 1 or 256 % p == 0:
        raise ValueError(f"p={p}: bad reduction (or not prime); poldisc(x^5-x) = -2^8, so only odd primes are good")
    return p


@provenance("PROVED", "item 15c: an abelian surface whose Frobenius char. poly at a good prime is irreducible over Q "
            "is Q-simple; L_3 = x^4-2x^2+9 is irreducible (PARI polisirreducible)")
def jacobian_is_Q_simple(primes=(3,)) -> bool:
    """True if for some prime in `primes` the Frobenius char. poly of y^2=x^5-x is irreducible over Q.

    3 is a prime of good reduction (disc(x^5-x) = -2^8); p=2 and other primes dividing 256 (bad
    reduction) raise ValueError, as do non-positive or non-int p."""
    primes = [_check_good_prime(p) for p in primes]
    pari = _pari()
    return any(bool(pari(f"polisirreducible(hyperellcharpoly(Mod(1,{p})*(x^5-x)))")) for p in primes)
