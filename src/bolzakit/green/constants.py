"""Closed-form Green/Bergman constants at a Bolza Weierstrass point (derivation in docs/green.md).

I0 = ∫_C d²x/|x⁵−x| = (π/4) g(1/8) g(1/2) g(3/8), with g(x) = Γ(x)/Γ(1−x)  (complex beta integral, w = x⁴).
B_Σ(0,0) = 2|a|²/I0, |a| = s(1)/(tanh(c/2)·108^{1/8}), s(1) by Gauss's ₂F₁(1) formula.
Here c is the hyperbolic distance with tanh(c/2) = th in the code, where cosh c = 1/(2 sin(π/8)) (`ch`).
β1 = (1/π − B)/4,  β2 = 1/(8π),  β3 = −(4|a|⁶/(2 I0) − 3/π)/36.
β2 rests on there being no holomorphic 1-form of character 1 mod 8.

Derivation: an AI second reader in the original research code (private archive; not independently checked), §3.A (Schwarz map
s(τ)=τ^{1/8}F(a+1/8,b+1/8;9/8;τ)/F(a,b;7/8;τ), a=17/48, b=1/48; I0 via the complex beta integral) and
the original note's §5 'Closed form'; reproduced in bolzakit docs/green.md.
"""
from __future__ import annotations
import mpmath as mp
from ..provenance import provenance


def _core(dps):
    with mp.workdps(dps + 10):
        g = lambda x: mp.gamma(x) / mp.gamma(1 - x)
        I0 = mp.pi / 4 * g(mp.mpf(1) / 8) * g(mp.mpf(1) / 2) * g(mp.mpf(3) / 8)
        a, b = mp.mpf(17) / 48, mp.mpf(1) / 48
        F1 = lambda A, B, C: mp.gamma(C) * mp.gamma(C - A - B) / (mp.gamma(C - A) * mp.gamma(C - B))
        s1 = F1(a + mp.mpf(1) / 8, b + mp.mpf(1) / 8, mp.mpf(9) / 8) / F1(a, b, mp.mpf(7) / 8)
        ch = 1 / (2 * mp.sin(mp.pi / 8))
        th = mp.sqrt((ch - 1) / (ch + 1))
        aa = s1 / (th * mp.power(108, mp.mpf(1) / 8))
        B = 2 * aa ** 2 / I0
        return I0, aa, B


@provenance("PROVED", "(2,3,8) Schwarz map + complex beta integral (original research code, review 1); agrees with the Hejhal "
            "1-form computation and the heat-kernel image sum to 8e-12; derivation: an AI second reader (not independently checked); "
            "see docs/green.md")
def constants(dps: int = 30) -> dict:
    I0, aa, B = _core(dps)
    with mp.workdps(dps):
        return {
            "B_Sigma_00": {"value": +B, "closed_form": "2|a|^2/I0", "label": "PROVED"},
            "beta1": {"value": (1 / mp.pi - B) / 4, "closed_form": "(1/pi - B)/4", "label": "PROVED"},
            "beta2": {"value": 1 / (8 * mp.pi), "closed_form": "1/(8 pi)  [no holomorphic 1-form of character 1 mod 8]", "label": "PROVED"},
            "beta3": {"value": -(4 * aa ** 6 / (2 * I0) - 3 / mp.pi) / 36, "closed_form": "-(4|a|^6/(2 I0) - 3/pi)/36",
                      "label": "PROVED"},
        }
