"""Regular N-rings of unit vortices centred at a Weierstrass point (Bolza) or at the origin (flat / H² models).

flat / sinh / tanh: energy-momentum Hessian in mpmath by central differences (item 03 ringmodes.py), step 10^(-dps/3).
bolza: H_Sigma = H_tanh + S_ring (item 03 dynamics.py): the tanh part in mpmath as above, the float64 smooth part
S_ring by central differences with step h = 2e-3 in the cos(m theta_k) radial direction, plus the Omega correction
-dOmega j''(s0), dOmega = (dS/ds)/j'(s0). Units: C1 (c_norm = 1/2pi, as item 03).
"""
from __future__ import annotations
import mpmath as mp
from ..errors import OutOfValidityRange
from ..provenance import provenance
from . import ringmodes as RM

_MODELS = {"flat": RM.H_flat, "sinh": RM.H_sinh, "tanh": RM.H_tanh}
BOLZA_H = 2e-3            # item 03 dynamics.analyse default step for the float64 S part


class Ring:
    def __init__(self, N: int, xi, model: str, dps: int):
        if model not in (*_MODELS, "bolza"):
            raise ValueError(f"model must be one of {(*_MODELS, 'bolza')}")
        if N < 2:
            raise ValueError("N >= 2")
        self.N, self.xi, self.model, self.dps = N, mp.mpf(xi), model, dps

    def _check_bolza(self):
        from ..green import VALIDITY_RADIUS
        if mp.sqrt(self.xi) > VALIDITY_RADIUS:
            raise OutOfValidityRange(f"ring radius sqrt(xi)={float(mp.sqrt(self.xi)):.3f} exceeds {VALIDITY_RADIUS}")

    def modes(self) -> dict:
        """{m: lambda_m} for m = 1..N-1, C1 units (radial cos(m theta_k) Rayleigh quotient of Hess(H - Omega J)).
        For model="bolza" this is item 03's `lam_sigma` (cos copy), computed for m = 1..N//2; lambda_{N-m} = lambda_m
        since cos((N-m) theta_k) = cos(m theta_k). The sin copy differs at O(xi^4) (C8-vs-C_N mismatch; 1e-5 at
        xi = 0.01) and is not returned."""
        if self.model == "bolza":
            self._check_bolza()                        # before any expensive work
            from ..green import green
            from . import dynamics as D
            core = green()._core
            half = D.lam_sigma(core, float(self.xi), range(1, self.N // 2 + 1), h=BOLZA_H, N=self.N, dps=self.dps)
            return {m: float(half[min(m, self.N - m)]) for m in range(1, self.N)}
        with mp.workdps(self.dps):
            cn = RM.c_norm(self.N)
            r = RM.analyse(_MODELS[self.model], self.N, self.xi, flat=(self.model == "flat"))
            return {m: float(r["alpha"][m] / cn) for m in range(1, self.N)}


@provenance("NUMERICAL", "energy–momentum Hessian by central differences: mpmath for flat/sinh/tanh (item 03 "
            "ringmodes.py); bolza = tanh (mpmath) + float64 S_ring smooth part, step 2e-3 (item 03 dynamics.py)",
            "flat/sinh/tanh O(h^2), h = 10^(-dps/3); bolza ~3e-6 at xi = 0.01 vs the λ3 series (Green tail + O(h^2))")
def ring(N: int, xi, model: str = "bolza", dps: int = 30) -> Ring:
    return Ring(N, xi, model, dps)


@provenance("PROVED", "item 03 note.tex: C_k momentum-map structure + closed-form β3; valid for N=7 in the "
            "Boatto–Koiller model with the degree-7 truncation of the regular part")
def lambda3_series() -> dict:
    from ..green import constants
    b3 = constants()["beta3"]["value"]
    return {"coefficients": [mp.mpf(0), mp.mpf(12), mp.mpf(-4), mp.mpf(-61) / 2 + 126 * mp.pi * b3],
            "c3_closed_form": "-61/2 + 126*pi*beta3",
            "valid_for": "N=7, Boatto–Koiller model, degree-7 truncation"}
