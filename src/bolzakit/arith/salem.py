"""The Bolza quartic x⁴−4x³−2x²−4x+1 = N_{Q(√2)/Q}(x²−(2+2√2)x+1); its largest root e^{ℓ_sys/2} is a Salem number."""
import mpmath as mp
from ..provenance import provenance


@provenance("EXACT", "item 15b (PARI: reciprocal, irreducible, norm from Q(√2))")
def bolza_quartic() -> list[int]:
    return [1, -4, -2, -4, 1]


@provenance("PROVED", "item 15b: largest root = 1+√2+√(2+2√2) = e^{ℓ_sys/2}; two roots on |x|=1")
def salem_number(dps: int = 50):
    with mp.workdps(dps):
        return 1 + mp.sqrt(2) + mp.sqrt(2 + 2 * mp.sqrt(2))
