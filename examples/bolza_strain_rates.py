"""Splitting rates of the triple eigenvalue lambda_1 = lambda_2 = lambda_3 of the Bolza surface (about 1 min).

At Bolza, lambda_1 has multiplicity 3. Under a deformation along a direction v it splits at first order with rates
equal to the eigenvalues of the 3 x 3 first-order matrix A(v) (Lemmas 2-4 of proofs/bolza-local-max/PROOF.md).
tr A = 0 (the cluster sum Lambda = lambda_1 + lambda_2 + lambda_3 is stationary at Bolza), so the rates sum to 0.
  F_A: A != 0, the triple splits linearly (Lemma 5).
  F_B: tangent to ker A (Lemma 7), so the rates vanish; the splitting is second order, and the Ky Fan bound
       D2 Lambda(v, v) <= kyfan_upper < 0 is certified (Lemma 8; P2 at n = 24 here, n = 96 in the proof package).
  a random unit vector: generic rates. Vector directions certify weakly at second order (minimum-norm second-order
       jet, docs/certified.md), so only the rates are printed for it.
The rates are float (NUMERICAL) eigenvalues of the discrete A (P2, n = 24), per unit t; only the Ky Fan bound is
CERTIFIED.
"""
import numpy as np

from bolzakit.certified import GapInput, Lam1Input, certify_cluster

# Explicit inputs: the certified lambda_1 ball (Lehmann-Goerisch lower end, P2 Rayleigh-Ritz upper end) and the
# lambda_4 >= 4.0 gap from the trace-formula counts. The certificates record them as dependencies.
LAM1 = Lam1Input("3.8388757010", "3.8388872830", "lg_lower_n48 (lower) and rr_p2_n96 (upper)", dep_id="D15")
GAP = GapInput("4.0", "trace-formula counts lp_A_3.3 and lp_B_3.3_4.0")


def main():
    v = np.random.default_rng(0).standard_normal(16)
    v /= np.linalg.norm(v)
    for name, direction in (("F_A", "F_A"), ("F_B", "F_B"), ("random vector", v)):
        r = certify_cluster(direction, degree=2, n=24, lam1=LAM1, gap=GAP)
        rates = ", ".join(f"{x:+.6f}" for x in r.A_float)
        print(f"{name:14s} splitting rates (eig A): {rates}   sum {sum(r.A_float):+.1e}")
    r = certify_cluster("F_B", degree=2, n=24, lam1=LAM1, gap=GAP, second_order=True)
    print(f"F_B Ky Fan: D2 Lambda(v, v) <= {r.kyfan_upper:.6f}  certified={r.certified}")
    assert r.certified and r.kyfan_upper < 0


if __name__ == "__main__":
    main()
