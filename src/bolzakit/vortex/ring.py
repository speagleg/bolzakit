"""Regular N-rings of unit vortices centred at a Weierstrass point (Bolza) or at the origin (flat / H² models).

flat / sinh / tanh: energy-momentum Hessian in mpmath by central differences (item 03 ringmodes.py), step 10^(-dps/3).
bolza: H_Sigma = H_tanh + S_ring (item 03 dynamics.py): the tanh part in mpmath as above, the float64 smooth part
S_ring by central differences with step h = 2e-3 in the cos(m theta_k) radial direction, plus the Omega correction
-dOmega j''(s0), dOmega = (dS/ds)/j'(s0); the same for the sin(m theta_k) direction, and modes() returns the mean
of the two copies (C7-symmetrised). Units: C1 (c_norm = 1/2pi, as item 03).
"""
from __future__ import annotations
import mpmath as mp
from ..errors import OutOfValidityRange
from ..provenance import provenance
from . import ringmodes as RM

_MODELS = {"flat": RM.H_flat, "sinh": RM.H_sinh, "tanh": RM.H_tanh}
BOLZA_H = 2e-3            # item 03 dynamics.analyse default step for the float64 S part


def _real_xi(xi):
    """xi as a real mpf; ValueError unless 0 < xi < 1 (complex input allowed only with zero imaginary part)."""
    try:
        x = mp.mpmathify(xi)
    except (TypeError, ValueError) as e:
        raise ValueError(f"xi must be a real number in (0, 1), got {xi!r}") from e
    if isinstance(x, mp.mpc):
        if x.imag != 0:
            raise ValueError(f"xi must be real, got {xi!r}")
        x = x.real
    if not (0 < x < 1):
        raise ValueError(f"xi = r^2 must satisfy 0 < xi < 1, got {xi!r}")
    return x


class Ring:
    """Regular N-ring, xi = r^2 (Euclidean disk radius squared). Validation happens here, before any expensive work."""

    def __init__(self, N: int, xi, model: str, dps: int, T: float = 2.6, D_ball: float = 15.0):
        if model not in (*_MODELS, "bolza"):
            raise ValueError(f"model must be one of {(*_MODELS, 'bolza')}")
        if N < 2:
            raise ValueError("N >= 2")
        self.N, self.xi, self.model, self.dps = N, _real_xi(xi), model, dps
        self.T, self.D_ball = float(T), float(D_ball)
        self._bolza = None
        if model == "bolza":
            from ..green import VALIDITY_RADIUS
            from ..green.surface import check_T_D
            self.T, self.D_ball = check_T_D(self.T, self.D_ball)
            if mp.sqrt(self.xi) > VALIDITY_RADIUS:
                raise OutOfValidityRange(
                    f"ring radius sqrt(xi)={float(mp.sqrt(self.xi)):.3f} exceeds {VALIDITY_RADIUS}")

    def _bolza_pairs(self) -> dict:
        """{m: (cos copy, sin copy)} for m = 1..N//2 (item 03 dynamics.py), cached."""
        if self._bolza is None:
            from ..green import green
            from . import dynamics as D
            core = green(self.T, self.D_ball)._core
            self._bolza = D.lam_sigma(core, float(self.xi), range(1, self.N // 2 + 1), h=BOLZA_H, N=self.N,
                                      dps=self.dps)
        return self._bolza

    def modes(self) -> dict:
        """{m: lambda_m} for m = 1..N-1, C1 units (radial-mode Rayleigh quotient of Hess(H - Omega J)).

        flat/sinh/tanh: the cos(m theta_k) radial mode (cos and sin copies coincide by the exact SO(2) symmetry).
        bolza: the MEAN of item 03's cos copy `lam_sigma` and sin copy `lam_sigma_sin`, i.e. the C_N-symmetrised
        eigenvalue, which is what the degree-7 truncated model and `lambda3_series` describe. The two copies differ
        by `split(m)`, an O(xi^4) symmetry breaking between the surface's C8 (at the Weierstrass point) and the
        ring's C7: about 1000 xi^4 for m = 3. lambda_{N-m} = lambda_m."""
        if self.model == "bolza":
            P = self._bolza_pairs()
            return {m: float((P[k][0] + P[k][1]) / 2) for m in range(1, self.N) for k in [min(m, self.N - m)]}
        with mp.workdps(self.dps):
            cn = RM.c_norm(self.N)
            r = RM.analyse(_MODELS[self.model], self.N, self.xi, flat=(self.model == "flat"))
            return {m: float(r["alpha"][m] / cn) for m in range(1, self.N)}

    def split(self, m: int) -> float:
        """cos-copy minus sin-copy radial eigenvalue of mode m, C1 units. 0.0 for flat/sinh/tanh (rotation invariant
        Hamiltonians). For bolza: O(xi^4) C8-vs-C7 symmetry breaking, about 1000 xi^4 for m = 3 at the Bolza
        Weierstrass point. split(N - m) = split(m) (sin((N-m) theta_k) = -sin(m theta_k); Rayleigh quotients agree)."""
        if not 1 <= m <= self.N - 1:
            raise ValueError(f"m must be in 1..{self.N - 1}")
        if self.model != "bolza":
            return 0.0
        c, s = self._bolza_pairs()[min(m, self.N - m)]
        return float(c - s)


@provenance("NUMERICAL", "energy–momentum Hessian by central differences: mpmath for flat/sinh/tanh (item 03 "
            "ringmodes.py); bolza = tanh (mpmath) + float64 S_ring smooth part, step 2e-3 (item 03 dynamics.py), "
            "C7-symmetrised (mean of cos and sin copies)",
            "flat/sinh/tanh O(h^2), h = 10^(-dps/3). bolza: residual vs the λ3 series is dominated by the O(ξ^4) "
            "C8-vs-C7 cos/sin splitting (≈1000 ξ^4 for m=3, removed by the symmetrisation) plus an O(ξ) "
            "heat-kernel-tail offset at T=2.6 (≈ -1.5e-6 at ξ=0.01); the O(h^2) differencing error is negligible")
def ring(N: int, xi, model: str = "bolza", dps: int = 30, *, T: float = 2.6, D_ball: float = 15.0) -> Ring:
    return Ring(N, xi, model, dps, T=T, D_ball=D_ball)


@provenance("PROVED", "item 03 note.tex: C_k momentum-map structure + closed-form β3; valid for N=7 in the "
            "Boatto–Koiller model with the degree-7 truncation of the regular part")
def lambda3_series() -> dict:
    from ..green import constants
    b3 = constants()["beta3"]["value"]
    return {"coefficients": [mp.mpf(0), mp.mpf(12), mp.mpf(-4), mp.mpf(-61) / 2 + 126 * mp.pi * b3],
            "c3_closed_form": "-61/2 + 126*pi*beta3",
            "valid_for": "N=7, Boatto–Koiller model, degree-7 truncation"}
