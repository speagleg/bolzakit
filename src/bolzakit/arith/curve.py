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


@provenance("PROVED", "item 15c: an abelian surface whose Frobenius char. poly at a good prime is irreducible over Q "
            "is Q-simple; L_3 = x^4-2x^2+9 is irreducible (PARI polisirreducible)")
def jacobian_is_Q_simple(primes=(3,)) -> bool:
    pari = _pari()
    return any(bool(pari(f"polisirreducible(hyperellcharpoly(Mod(1,{p})*(x^5-x)))")) for p in primes)
