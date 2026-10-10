"""Arithmetic of the Bolza curve y² = x⁵ − x over Q (see docs/arith.md). Needs the optional extra: pip install bolzakit[pari].
The Jacobian is Q-simple (L_3 = x⁴−2x²+9 is irreducible over Q) and splits over Q(√2) as E′² with CM by Z[√−2]."""
from __future__ import annotations
from ..provenance import provenance, require


def _pari():
    return require("cypari2", "pari").Pari()


@provenance("EXACT", "PARI hyperellcharpoly")
def euler_factor(p: int) -> list[int]:
    p = _check_good_prime(p)
    pari = _pari()
    P = pari(f"hyperellcharpoly(Mod(1,{p})*(x^5-x))")
    return [int(c) for c in pari.Vec(P)]


def _check_good_prime(p) -> int:
    from .golden import _is_prime
    if isinstance(p, bool) or not isinstance(p, int) or p <= 0:
        raise ValueError(f"p must be a positive int, got {p!r}")
    if not _is_prime(p):
        raise ValueError(f"p={p} is not prime")
    if p == 2:
        raise ValueError("p=2 is a prime of bad reduction: poldisc(x^5-x) = -2^8")
    return p


@provenance("PROVED", "proof in docs/arith.md: an abelian surface whose Frobenius char. poly at a good prime is irreducible over Q "
            "is Q-simple; L_3 = x^4-2x^2+9 is irreducible (PARI polisirreducible). True is proved; None means "
            "inconclusive at the given primes")
def jacobian_is_Q_simple(primes=(3,)) -> bool | None:
    """True if for some prime in `primes` the Frobenius char. poly of y^2=x^5-x is irreducible over Q; else None.

    The criterion is one-sided: irreducibility of L_p at one good prime proves Q-simplicity, but a reducible L_p
    proves nothing (e.g. L_5, L_7, L_13, L_17, L_23 are reducible although the Jacobian is Q-simple). So the result
    is True (proved) or None (inconclusive at the given primes); it is never False.

    3 is a prime of good reduction (disc(x^5-x) = -2^8); p=2 (bad reduction) raises ValueError, as do
    non-prime, non-positive or non-int p."""
    primes = [_check_good_prime(p) for p in primes]
    pari = _pari()
    if any(bool(pari(f"polisirreducible(hyperellcharpoly(Mod(1,{p})*(x^5-x)))")) for p in primes):
        return True
    return None
