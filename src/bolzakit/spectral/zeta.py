from fractions import Fraction
from ..provenance import provenance


@provenance("PROVED", "heat-trace constant term χ/6 = -1/3 (identity term of the Selberg trace formula) minus "
            "dim ker Δ = 1")
def zeta_zero() -> Fraction:
    """Spectral zeta function at 0: ζ_Δ(0) = χ/6 − 1 = −4/3."""
    return Fraction(-4, 3)
