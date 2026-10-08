"""Arithmetic companions: Γ(p) golden geodesic, Bolza quartic / Salem number, curve arithmetic (PARI extra)."""
from .golden import golden_geodesic_is_systole, lucas  # noqa: F401
from .salem import bolza_quartic, salem_number  # noqa: F401
from . import curve  # noqa: F401  (functions raise OptionalDependencyMissing without cypari2)
