/-
# Certificate arithmetic (ARITHMETIC-CHECKED)

**Scope, stated plainly.**  These theorems check only the *final rational arithmetic* that turns
the printed Arb enclosures of the certificate runs into the inequalities used by PROOF.md.  The
numbers below are copied from the `.out` / `.json` files in `work/` (endpoints taken outward,
i.e. "midpoint + radius" for an upper bound and "midpoint − radius" for a lower bound).
Nothing here re-runs or verifies the Arb/python-flint computations themselves (trace-formula LP,
FEM element integrals, Sturm counts, gap/transfer bounds, determinants): those remain CERTIFIED
by the Python scripts (D13), not by Lean.
-/
import Mathlib

namespace Bolza.Cert

/-! ### Lemma 1(a): trace-formula LP (`lp_count.out`) -/

/-- run A: bound ∈ [0.912612352675290375 ± 9.95e-20], target < 1 ⇒ Spec ∩ (0, 3.3] = ∅. -/
theorem lpA_upper_lt_one : (912612352675290375 : ℚ) / 10 ^ 18 + 995 / 10 ^ 22 < 1 := by norm_num

/-- run B: bound ∈ [3.86966111976973759 ± 1.65e-18], target < 4 ⇒ #[3.3, 4.0] ≤ 3. -/
theorem lpB_upper_lt_four : (386966111976973759 : ℚ) / 10 ^ 17 + 165 / 10 ^ 20 < 4 := by norm_num

/-- `lpB_upper_lt_four` gives a count ≤ 3 only because the count is an integer:
an integer `N ≤ 3.8697 < 4` is `≤ 3`. -/
theorem count_le_three (N : ℕ) (h : (N : ℚ) ≤ 386966111976973759 / 10 ^ 17 + 165 / 10 ^ 20) :
    N ≤ 3 := by
  have h4 : (N : ℚ) < 4 := lt_of_le_of_lt h lpB_upper_lt_four
  have : N < 4 := by exact_mod_cast h4
  omega

/-! ### Lemma 1(b),(c): Rayleigh–Ritz bounds and the λ₁ interval -/

/-- P1, n = 192: certified μ₃ ≤ 3.839184804702046 ≤ 3.8391848048 (rounded outward). -/
theorem p1_mu3 : (3839184804702046 : ℚ) / 10 ^ 15 ≤ 38391848048 / 10 ^ 10 := by norm_num

/-- gap: λ₄ > 4.0 and λ₃ ≤ 3.8391848048 give a gap ≥ 0.1608. -/
theorem gap_ge : (4 : ℚ) - 38391848048 / 10 ^ 10 ≥ 1608 / 10 ^ 4 := by norm_num

/-- the gap used by F1 (Lemma 2, r₀ = 0.16) -/
theorem gap_ge_r0 : (4 : ℚ) - 38391848048 / 10 ^ 10 > 16 / 100 := by norm_num

/-- P2, n = 96: certified μ₃ ≤ 3.8388872829261946 ≤ 3.8388872830. -/
theorem p2_mu3 : (38388872829261946 : ℚ) / 10 ^ 16 ≤ 38388872830 / 10 ^ 10 := by norm_num

/-- SU13 interval [3.83888625884, 3.83888825884] (midpoint 3.83888725884 ± 1e-6) intersected with
the P2 upper bound is the nonempty interval [3.83888625884, 3.8388872830] of Lemma 1(c). -/
theorem lambda1_interval :
    (383888625884 : ℚ) / 10 ^ 11 = 383888725884 / 10 ^ 11 - 1 / 10 ^ 6 ∧
    (383888625884 : ℚ) / 10 ^ 11 < 38388872830 / 10 ^ 10 ∧
    (38388872830 : ℚ) / 10 ^ 10 < 383888825884 / 10 ^ 11 := by norm_num

/-- Note G1: the scripts' λ-ball [3.8388862587999, 3.8388882588001] contains the Lemma 1(c)
interval [3.83888625884, 3.8388872830], so every certificate run with that ball is valid. -/
theorem script_ball_covers :
    (38388862587999 : ℚ) / 10 ^ 13 ≤ 383888625884 / 10 ^ 11 ∧
    (38388872830 : ℚ) / 10 ^ 10 ≤ 38388882588001 / 10 ^ 13 := by norm_num

/-- KMP's rigorous bound λ₁ < 3.8388976481 for genus-2 orbifolds lies ≈ 1.04e-5 above the
certified upper end 3.8388872830 (consistency, not used). -/
theorem kmp_consistent :
    (38388872830 : ℚ) / 10 ^ 10 < 38388976481 / 10 ^ 10 ∧
    (38388976481 : ℚ) / 10 ^ 10 - 38388872830 / 10 ^ 10 = 103651 / 10 ^ 10 := by norm_num

/-! ### D1: SU13 slack -/

/-- The sign certificate needs only λ₁ ≥ 3.8388029653 (bisection, review_sign_1); the SU13 lower
end 3.83888625884 exceeds it by 8.329e-5, i.e. ≈ 83.3 × SU13's 1e-6.  (PROOF.md quotes
"≈ 8.4e-5 / 84×" in one place, which is measured from the *upper* end 3.8388872830; from the
guaranteed lower end the slack is 83.3×.  Either way ≥ 80×, as stated in §Dependencies.) -/
theorem su13_slack :
    (383888625884 : ℚ) / 10 ^ 11 - 38388029653 / 10 ^ 10 = 8329354 / 10 ^ 11 ∧
    (383888625884 : ℚ) / 10 ^ 11 - 38388029653 / 10 ^ 10 ≥ 80 * (1 / 10 ^ 6) ∧
    (38388872830 : ℚ) / 10 ^ 10 - 38388029653 / 10 ^ 10 ≥ 84 * (1 / 10 ^ 6) := by norm_num

/-- Lemma 5 needs only λ₁ ≥ 3.83788, satisfied by the SU13 lower end. -/
theorem lemma5_threshold : (383788 : ℚ) / 10 ^ 5 ≤ 383888625884 / 10 ^ 11 := by norm_num

/-! ### Lemma 5: A ≠ 0 (`cert_first_order.out`) -/

/-- n = 192, joint certified upper end b(e,e)/m(u) ≤ −3.0875694043934345 < 0. -/
theorem lemma5_joint_192 : (-30875694043934345 : ℚ) / 10 ^ 16 < 0 := by norm_num

/-- n = 192, re-assembled from the printed component balls:
b(e,e)/m ≤ b(u,u)/m + |transfer|/m ≤ (−4.806192 + 1.96e-7) + (1.71 + 8.63e-3) = −3.087561804 < 0.
(Weaker than the joint value, as it must be; still negative.) -/
theorem lemma5_assembled_192 :
    ((-4806192 : ℚ) / 10 ^ 6 + 196 / 10 ^ 9) + (171 / 100 + 863 / 10 ^ 5) = -3087561804 / 10 ^ 9 ∧
    (-3087561804 : ℚ) / 10 ^ 9 < 0 := by norm_num

/-- n = 96 likewise: (−4.806689 + 1.39e-7) + (3.78 + 7.16e-3) = −1.019528861 < 0. -/
theorem lemma5_assembled_96 :
    ((-4806689 : ℚ) / 10 ^ 6 + 139 / 10 ^ 9) + (378 / 100 + 716 / 10 ^ 5) < 0 := by norm_num

/-! ### Lemma 8: Ky Fan sign bound (`cert_sign_48_96.json`, n = 96) -/

/-- The certified joint upper end −12.850476214196528 ≤ −12.8504 < 0 (outward rounding). -/
theorem sign_joint : (-12850476214196528 : ℚ) / 10 ^ 15 ≤ -128504 / 10 ^ 4 ∧
    (-128504 : ℚ) / 10 ^ 4 < 0 := by norm_num

/-- The hypothesis constant used by `BolzaInputs.g_sign` (Lehmann–Goerisch certificate, primary):
−9.3411 < 0, and the SU13 alternative −12.8504 is at least as strong, so it also discharges it. -/
theorem sign_su13_implies_lg : (-128504 : ℚ) / 10 ^ 4 ≤ -93411 / 10 ^ 4 ∧
    (-93411 : ℚ) / 10 ^ 4 < 0 := by norm_num

/-- **Finding (consistent with PROOF.md Lemma 8 step 5, "one-at-a-time" rows):** the sign bound
can NOT be assembled by adding the printed one-at-a-time error rows to the central value.
Central upper end −14.2742 + 2.23e-5 plus the six rows
(A 1.5016520626076613e-5, M0 4.6164514537494483e-4, M1 1.5021804037118613e-5,
 X1 0.839368051438381, X2 0.5700597204815345, c1 7.367181513551202e-5)
gives about −12.8642, which is *below* the joint certified upper end −12.8505: the rows sum to
≈1.4100 while the joint error is ≈1.4237.  So the row sum is not an upper bound; only the joint
Arb ball (`sign_joint`) certifies the sign. -/
theorem sign_rows_not_additive :
    let rows : ℚ := 15016520626076613 / 10 ^ 21 + 46164514537494483 / 10 ^ 20 +
      15021804037118613 / 10 ^ 21 + 839368051438381 / 10 ^ 15 + 5700597204815345 / 10 ^ 16 +
      7367181513551202 / 10 ^ 20
    let central_up : ℚ := -142742 / 10 ^ 4 + 223 / 10 ^ 7
    central_up + rows < -12850476214196528 / 10 ^ 15 ∧
    (-12850476214196528 : ℚ) / 10 ^ 15 - central_up > rows := by
  norm_num

/-! ### Lemma 4 / 9: determinants (`cert_rank.out`, `cert_dpi_rank.out`) -/

/-- d(angle sum)/dR = [−5.572971195984104378981991909236 ± 5.18e-31] excludes 0. -/
theorem dangle_ne_zero :
    (-5572971195984104378981991909236 : ℚ) / 10 ^ 30 + 518 / 10 ^ 33 < 0 := by norm_num

/-- rank J = 5: 5×5 minor det = [46151927.03406868578682884302258 ± 2.89e-24] excludes 0. -/
theorem rankJ_minor : (4615192703406868578682884302258 : ℚ) / 10 ^ 23 - 289 / 10 ^ 26 > 0 := by
  norm_num

/-- G3: best 6×6 minor det = [7872428136359854080.0 ± 3e-6] excludes 0 (rank dπ = 6). -/
theorem dpi_minor : (7872428136359854080 : ℚ) - 3 / 10 ^ 6 > 0 := by norm_num

/-- G3 cross-check: 11×11 minor det = [−1.3827971011809151404e24 ± 2.74e4] excludes 0. -/
theorem dpi_cross : (-13827971011809151404 : ℚ) * 10 ^ 5 + 274 * 10 ^ 2 < 0 := by norm_num

/-! ### Exact algebraic identities at the regular octagon (PROVED, not just arithmetic) -/

/-- Lemma 7 setup: at p = 0, C = cosh R = 3 + 2√2 and x = y = tan(π/8) = √2 − 1 satisfy the
F_B surface condition C²xy − C(x+y) − 1 = 0 (the identification tan(π/8) = √2 − 1 is standard
and not re-proved here). -/
theorem FB_surface_condition :
    let C := 3 + 2 * Real.sqrt 2
    let x := Real.sqrt 2 - 1
    C ^ 2 * x * x - C * (x + x) - 1 = 0 := by
  intro C x
  have h : Real.sqrt 2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  simp only [C, x]
  ring_nf
  rw [show Real.sqrt 2 ^ 4 = (Real.sqrt 2 ^ 2) ^ 2 by ring, show Real.sqrt 2 ^ 3 = Real.sqrt 2 ^ 2 * Real.sqrt 2 by ring, h]
  ring

/-- Side pairing trace: with cosh(L/2) = 1 + √2, tr g_k = 1 + 2 cosh L = 11 + 8√2, and this lies
in the certified ball [22.3137084989848 ± 3.97e-14] of `cert_rank.out`. -/
theorem pairing_trace :
    1 + 2 * (2 * (1 + Real.sqrt 2) ^ 2 - 1) = 11 + 8 * Real.sqrt 2 ∧
    (223137084989848 : ℝ) / 10 ^ 13 - 397 / 10 ^ 16 < 11 + 8 * Real.sqrt 2 ∧
    11 + 8 * Real.sqrt 2 < (223137084989848 : ℝ) / 10 ^ 13 + 397 / 10 ^ 16 := by
  have h : Real.sqrt 2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  have hs0 : 0 ≤ Real.sqrt 2 := Real.sqrt_nonneg 2
  refine ⟨by ring_nf; rw [h]; ring, ?_, ?_⟩
  · -- √2 > 1.4142135623730950488 (the margin is only ≈ 8e-17)
    have : (14142135623730950488 : ℝ) / 10 ^ 19 < Real.sqrt 2 := by
      rw [Real.lt_sqrt (by norm_num)]; norm_num
    linarith
  · have : Real.sqrt 2 < (141421356237310 : ℝ) / 10 ^ 14 := by
      rw [Real.sqrt_lt' (by norm_num)]; norm_num
    linarith

end Bolza.Cert
