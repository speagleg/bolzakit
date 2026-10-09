"""Deform the Bolza surface and certify upper bounds on its first three Laplace eigenvalues (CERTIFIED, about 15 s).

Three directions in Teichmueller space, each at t = 0.02 (deform() is validated for |t| <= 0.05 only):
  "F_A"  equilateral octagon with alternating angles;
  "F_B"  equiangular octagon with alternating sides;
  a unit 16-vector of Klein vertex coordinates, projected onto the constraint manifold (certified Krawczyk step).
For each deformed surface, certify_spectrum (P2 Rayleigh-Ritz, n = 12) returns certified UPPER bounds on
lambda_1..lambda_3. If the bound on lambda_1 is below the certified lower bound lambda_1(Bolza) >= 3.8388757010
(Lehmann-Goerisch, proofs/bolza-local-max/certificates/lg_lower_n48.json, conditional on lambda_4 >= 4.0), the
deformed surface certainly has a smaller lambda_1 than Bolza at that t. This is a check at one point; it is not the
local-maximum proof, which needs the second-order (Ky Fan) argument.
"""
import numpy as np

from bolzakit.certified import Surface, certify_spectrum
from bolzakit.errors import OutOfValidityRange

LAM1_BOLZA_LOWER = 3.8388757010   # certified (lg_lower_n48), conditional on lambda_4 >= 4.0


def main():
    bolza = Surface.bolza()
    v = np.random.default_rng(0).standard_normal(16)
    v /= np.linalg.norm(v)
    print(f"{'surface':24s} {'lambda_1 <=':>14s} {'lambda_2 <=':>14s} {'lambda_3 <=':>14s}  below Bolza's lambda_1?")
    for direction in ("F_A", "F_B", v):
        surface = bolza.deform(direction, 0.02)
        r = certify_spectrum(surface, k=3, degree=2, n=12)
        assert r.certified
        below = r.upper[0] < LAM1_BOLZA_LOWER
        print(f"{surface.label:24s} {r.upper[0]:14.9f} {r.upper[1]:14.9f} {r.upper[2]:14.9f}  {below}")
    try:
        bolza.deform("F_A", 0.06)
    except OutOfValidityRange as e:
        print(f"t = 0.06: OutOfValidityRange ({e})")
    else:
        raise AssertionError("deform(t=0.06) should raise OutOfValidityRange")


if __name__ == "__main__":
    main()
