/-
Copyright: item 23 (Bolza λ₁ extremality), Lean skeleton.

# Lemma 10 of work/PROOF.md: the local-max criterion (PROVED IN LEAN, no axioms beyond Mathlib's)

Setting.  `V` is any finite-dimensional real normed space (in the application `V = ℝ⁶`, the
tangent space of the chart `S` of Lemma 9, identified with TT tensors via σ).
`M : V → Mat₃` is the 3×3 cluster matrix of Lemma 2 (F1).  We assume the second-order
expansion

  M(v) = λ*·I + A(v) + B(v,v) + R(v),  ‖R(v)‖ = o(‖v‖²),

with `A` linear, every `A v` symmetric and traceless, and the trace of the quadratic part
negative definite on `K = ker A`:  tr B(k,k) ≤ -κ‖k‖².

Conclusion.  For every small `v ≠ 0` there is a vector `x ≠ 0` with
`xᵀ M(v) x < λ* · xᵀx`, i.e. `λ_min(M v) < λ*` (see `exists_eigenvalue_lt` for the
eigenvalue form when `M v` is symmetric).

Compared with PROOF.md Lemma 10:
* the hypothesis "A ≠ 0" is NOT needed (if `A = 0` the trace argument covers everything);
* symmetry of `M v` is only needed for the eigenvalue reformulation;
* the case split is |w| ≥ C‖v‖² vs |w| < C‖v‖² (PROOF.md uses C|k|², equivalent);
* `W` is an arbitrary algebraic complement of `K` (no inner product needed);
* in case (ii) we use `tr A(v) = 0` directly (PROOF.md writes `tr A(w)`, same thing since
  `A k = 0`).
The constant `c_A > 0` of PROOF.md step 1 is obtained here from the explicit elementary bound
`traceless_symm_exists_neg` (a traceless symmetric 3×3 matrix `S` has a test vector with entries
in {0,±1} and `xᵀSx ≤ -‖S‖/4`, sup-norm) together with the antilipschitz bound of the
injective map `A|_W`.
-/
import Mathlib

open Matrix Filter Topology

namespace Bolza

/-- 3×3 real matrices. -/
abbrev Mat3 := Matrix (Fin 3) (Fin 3) ℝ

/-- 3×3 real arrays with the (global) sup norm `‖S‖ = max |S i j|`.  We use the plain function
type rather than `Matrix` so that Mathlib's canonical normed-space instances apply. -/
abbrev Sq3 := Fin 3 → Fin 3 → ℝ

/-- the quadratic form `xᵀ S x` -/
def qf (S : Sq3) (x : Fin 3 → ℝ) : ℝ := ∑ i, ∑ j, x i * S i j * x j

/-- trace -/
def tr3 (S : Sq3) : ℝ := S 0 0 + S 1 1 + S 2 2

/-- symmetry -/
def IsSym3 (S : Sq3) : Prop := ∀ i j, S i j = S j i

/-- identity -/
def one3 : Sq3 := fun i j => if i = j then 1 else 0

lemma qf_eq_matrix (S : Sq3) (x : Fin 3 → ℝ) : qf S x = x ⬝ᵥ (Matrix.of S *ᵥ x) := by
  simp [qf, dotProduct, mulVec, Finset.mul_sum, mul_assoc]

lemma qf_add (S T : Sq3) (x : Fin 3 → ℝ) : qf (S + T) x = qf S x + qf T x := by
  simp [qf, mul_add, add_mul, Finset.sum_add_distrib]

lemma qf_sub (S T : Sq3) (x : Fin 3 → ℝ) : qf (S - T) x = qf S x - qf T x := by
  simp [qf, mul_sub, sub_mul, Finset.sum_sub_distrib]

lemma qf_smul_one (c : ℝ) (x : Fin 3 → ℝ) : qf (c • one3) x = c * (x ⬝ᵥ x) := by
  simp [qf, one3, dotProduct, Fin.sum_univ_three]; ring

lemma tr3_add (S T : Sq3) : tr3 (S + T) = tr3 S + tr3 T := by simp [tr3]; ring
lemma tr3_sub (S T : Sq3) : tr3 (S - T) = tr3 S - tr3 T := by simp [tr3]; ring

lemma entry_le_norm (S : Sq3) (i j : Fin 3) : |S i j| ≤ ‖S‖ := by
  have h1 : ‖S i j‖ ≤ ‖S i‖ := norm_le_pi_norm (S i) j
  have h2 : ‖S i‖ ≤ ‖S‖ := norm_le_pi_norm S i
  simpa [Real.norm_eq_abs] using h1.trans h2

lemma abs_qf_le (S : Sq3) (x : Fin 3 → ℝ) (hx : ∀ i, |x i| ≤ 1) : |qf S x| ≤ 9 * ‖S‖ := by
  have h1 : ∀ i j, |x i * S i j * x j| ≤ ‖S‖ := by
    intro i j
    rw [abs_mul, abs_mul]
    have := entry_le_norm S i j
    have := hx i; have := hx j
    have := abs_nonneg (x i); have := abs_nonneg (x j); have := abs_nonneg (S i j)
    calc |x i| * |S i j| * |x j| ≤ 1 * |S i j| * 1 := by gcongr
      _ ≤ ‖S‖ := by linarith
  calc |qf S x| ≤ ∑ i, |∑ j, x i * S i j * x j| := Finset.abs_sum_le_sum_abs _ _
    _ ≤ ∑ i, ∑ j, |x i * S i j * x j| :=
        Finset.sum_le_sum fun i _ => Finset.abs_sum_le_sum_abs _ _
    _ ≤ ∑ _i : Fin 3, ∑ _j : Fin 3, ‖S‖ :=
        Finset.sum_le_sum fun i _ => Finset.sum_le_sum fun j _ => h1 i j
    _ = 9 * ‖S‖ := by simp; ring

lemma abs_tr3_le (S : Sq3) : |tr3 S| ≤ 3 * ‖S‖ := by
  unfold tr3
  have := entry_le_norm S 0 0
  have := entry_le_norm S 1 1
  have := entry_le_norm S 2 2
  have := abs_add_three (S 0 0) (S 1 1) (S 2 2)
  linarith

/-- Elementary core of `c_A > 0` (PROOF.md Lemma 10 step 1): a traceless symmetric 3×3 matrix
has a test vector with entries of absolute value ≤ 1 on which its quadratic form is
`≤ -‖S‖/4` (sup norm). -/
lemma traceless_symm_exists_neg (S : Sq3) (hS : IsSym3 S) (htr : tr3 S = 0) :
    ∃ x : Fin 3 → ℝ, (∀ i, |x i| ≤ 1) ∧ qf S x ≤ -‖S‖ / 4 := by
  by_contra hcon
  push Not at hcon
  have h10 : S 1 0 = S 0 1 := hS 1 0
  have h20 : S 2 0 = S 0 2 := hS 2 0
  have h21 : S 2 1 = S 1 2 := hS 2 1
  unfold tr3 at htr
  have test : ∀ a b c : ℝ, |a| ≤ 1 → |b| ≤ 1 → |c| ≤ 1 →
      -‖S‖ / 4 < qf S ![a, b, c] := by
    intro a b c ha hb hc
    apply hcon
    intro i; fin_cases i <;> simpa
  have q : ∀ a b c : ℝ, qf S ![a, b, c] =
      a * (S 0 0 * a + S 0 1 * b + S 0 2 * c) + b * (S 1 0 * a + S 1 1 * b + S 1 2 * c)
        + c * (S 2 0 * a + S 2 1 * b + S 2 2 * c) := by
    intro a b c
    simp [qf, Fin.sum_univ_three]; ring
  have t1 := test 1 0 0 (by norm_num) (by norm_num) (by norm_num)
  have t2 := test 0 1 0 (by norm_num) (by norm_num) (by norm_num)
  have t3 := test 0 0 1 (by norm_num) (by norm_num) (by norm_num)
  have t4 := test 1 1 0 (by norm_num) (by norm_num) (by norm_num)
  have t5 := test 1 (-1) 0 (by norm_num) (by norm_num) (by norm_num)
  have t6 := test 1 0 1 (by norm_num) (by norm_num) (by norm_num)
  have t7 := test 1 0 (-1) (by norm_num) (by norm_num) (by norm_num)
  have t8 := test 0 1 1 (by norm_num) (by norm_num) (by norm_num)
  have t9 := test 0 1 (-1) (by norm_num) (by norm_num) (by norm_num)
  simp only [q, h10, h20, h21] at t1 t2 t3 t4 t5 t6 t7 t8 t9
  ring_nf at t1 t2 t3 t4 t5 t6 t7 t8 t9
  set t := ‖S‖ with ht
  have t0 : 0 ≤ t := norm_nonneg S
  rcases t0.lt_or_eq with tpos | tzero
  · have : ‖S‖ < t := by
      refine (pi_norm_lt_iff tpos).2 fun i => (pi_norm_lt_iff tpos).2 fun j => ?_
      rw [Real.norm_eq_abs]
      fin_cases i <;> fin_cases j <;> simp only [Fin.zero_eta, Fin.mk_one, Fin.reduceFinMk, h10, h20, h21] <;>
        (rw [abs_lt]; constructor <;> linarith)
    exact lt_irrefl _ this
  · rw [← tzero] at t1 t2 t3
    linarith

/-! ## Two real-arithmetic lemmas (the inequalities of Lemma 10, cases (i) and (ii)) -/

lemma case1_arith {qa qb qr a C β ε r w : ℝ} (ha : 0 < a) (hr : 0 < r) (hβ : 0 ≤ β)
    (hε1 : ε ≤ 1) (hC : a * C = 36 * (β + 1) + a) (hw : C * r ^ 2 ≤ w)
    (hqa : qa ≤ -(a * w) / 4) (hqb : qb ≤ 9 * (β * r ^ 2)) (hqr : qr ≤ 9 * (ε * r ^ 2)) :
    qa + qb + qr < 0 := by
  have hr2 : 0 < r ^ 2 := by positivity
  have h1 : a * (C * r ^ 2) ≤ a * w := mul_le_mul_of_nonneg_left hw ha.le
  have h2 : a * (C * r ^ 2) = (36 * (β + 1) + a) * r ^ 2 := by rw [← mul_assoc, hC]
  have h3 : ε * r ^ 2 ≤ 1 * r ^ 2 := mul_le_mul_of_nonneg_right hε1 hr2.le
  have h4 : 0 < a * r ^ 2 := mul_pos ha hr2
  nlinarith

lemma case2_arith {tk t2 t3 t4 t5 κ C β ε r kn wn D : ℝ} (hκ : 0 < κ) (hr : 0 < r)
    (hr1 : r ≤ 1) (hCr : C * r ≤ 1) (hC : 0 < C) (hβ : 0 ≤ β) (hεκ : ε ≤ κ / 12)
    (hwn : wn ≤ C * r ^ 2) (hwn0 : 0 ≤ wn) (hkn : r - C * r ^ 2 ≤ kn) (hkn2 : kn ≤ 2 * r)
    (htk : tk ≤ -κ * kn ^ 2) (ht2 : t2 ≤ 3 * (β * kn * wn)) (ht3 : t3 ≤ 3 * (β * wn * kn))
    (ht4 : t4 ≤ 3 * (β * wn * wn)) (ht5 : t5 ≤ 3 * (ε * r ^ 2))
    (hD : D = 2 * κ * C + 12 * β * C + 3 * β * C ^ 2) (hDr : D * r ≤ κ / 4) :
    tk + t2 + t3 + t4 + t5 < 0 := by
  have hr2 : 0 < r ^ 2 := by positivity
  have hr3 : 0 < r ^ 3 := by positivity
  have hCr2 : C * r ^ 2 ≤ r := by
    have := mul_le_mul_of_nonneg_right hCr hr.le
    calc C * r ^ 2 = C * r * r := by ring
      _ ≤ 1 * r := this
      _ = r := one_mul r
  have hpos : 0 ≤ r - C * r ^ 2 := by linarith
  have hk2 : (r - C * r ^ 2) ^ 2 ≤ kn ^ 2 := pow_le_pow_left₀ hpos hkn 2
  have hsq : r ^ 2 - 2 * C * r ^ 3 ≤ (r - C * r ^ 2) ^ 2 := by
    have hx : (r - C * r ^ 2) ^ 2 = r ^ 2 - 2 * C * r ^ 3 + (C * r ^ 2) ^ 2 := by ring
    rw [hx]; linarith [sq_nonneg (C * r ^ 2)]
  have e1 : tk ≤ -κ * r ^ 2 + 2 * κ * C * r ^ 3 := by
    have hx := mul_le_mul_of_nonneg_left (hsq.trans hk2) hκ.le
    have hy : κ * (r ^ 2 - 2 * C * r ^ 3) = κ * r ^ 2 - 2 * κ * C * r ^ 3 := by ring
    linarith
  have hkw : kn * wn ≤ 2 * C * r ^ 3 := by
    have hx := mul_le_mul hkn2 hwn hwn0 (by positivity)
    have hy : 2 * r * (C * r ^ 2) = 2 * C * r ^ 3 := by ring
    linarith
  have hww : wn * wn ≤ C ^ 2 * r ^ 4 := by
    have hx := mul_le_mul hwn hwn hwn0 (by positivity)
    have hy : C * r ^ 2 * (C * r ^ 2) = C ^ 2 * r ^ 4 := by ring
    linarith
  have e2 : t2 + t3 ≤ 12 * β * C * r ^ 3 := by
    have hx := mul_le_mul_of_nonneg_left hkw hβ
    have hy1 : β * kn * wn = β * (kn * wn) := by ring
    have hy2 : β * wn * kn = β * (kn * wn) := by ring
    have hy3 : β * (2 * C * r ^ 3) = 2 * β * C * r ^ 3 := by ring
    linarith
  have hr4 : r ^ 4 ≤ r ^ 3 := by
    have hx := mul_le_mul_of_nonneg_left hr1 hr3.le
    have hy : r ^ 4 = r ^ 3 * r := by ring
    linarith
  have e4 : t4 ≤ 3 * β * C ^ 2 * r ^ 3 := by
    have hx1 := mul_le_mul_of_nonneg_left hww hβ
    have hx2 := mul_le_mul_of_nonneg_left hr4 (by positivity : (0:ℝ) ≤ β * C ^ 2)
    have hy1 : β * (wn * wn) = β * wn * wn := by ring
    have hy2 : β * (C ^ 2 * r ^ 4) = β * C ^ 2 * r ^ 4 := by ring
    have hy3 : β * C ^ 2 * (r ^ 4) = β * C ^ 2 * r ^ 4 := by ring
    linarith
  have e5 : t5 ≤ κ / 4 * r ^ 2 := by
    have hx := mul_le_mul_of_nonneg_right hεκ hr2.le
    linarith
  have hDr3 : D * r ^ 3 ≤ κ / 4 * r ^ 2 := by
    have hx : D * r ^ 3 = (D * r) * r ^ 2 := by ring
    rw [hx]; exact mul_le_mul_of_nonneg_right hDr hr2.le
  have hDexp : D * r ^ 3 = 2 * κ * C * r ^ 3 + 12 * β * C * r ^ 3 + 3 * β * C ^ 2 * r ^ 3 := by
    rw [hD]; ring
  have hκr : 0 < κ * r ^ 2 := mul_pos hκ hr2
  linarith

/-! ## The criterion -/

section criterion

variable {V : Type*} [NormedAddCommGroup V] [NormedSpace ℝ V] [FiniteDimensional ℝ V]

/-- **PROOF.md Lemma 10 (local-max criterion), expansion form.  PROVED IN LEAN.**

`M v = λ•1 + A v + B v v + o(‖v‖²)`, every `A v` symmetric with `tr (A v) = 0`, and
`tr (B k k) ≤ -κ‖k‖²` on `ker A` (κ > 0).  Then for all small `v ≠ 0` there is `x ≠ 0` with
`qf (M v) x = xᵀ M(v) x < λ·xᵀx`. -/
theorem local_max_criterion (M : V → Sq3) (lam : ℝ) (A : V →L[ℝ] Sq3)
    (B : V →L[ℝ] V →L[ℝ] Sq3)
    (hA_symm : ∀ v, IsSym3 (A v)) (hA_tr : ∀ v, tr3 (A v) = 0)
    (hrem : ∀ ε > 0, ∀ᶠ v in 𝓝 (0 : V), ‖M v - (lam • one3 + A v + B v v)‖ ≤ ε * ‖v‖ ^ 2)
    (κ : ℝ) (hκ : 0 < κ)
    (hneg : ∀ k ∈ LinearMap.ker (A : V →ₗ[ℝ] Sq3), tr3 (B k k) ≤ -κ * ‖k‖ ^ 2) :
    ∀ᶠ v in 𝓝 (0 : V), v ≠ 0 → ∃ x : Fin 3 → ℝ, x ≠ 0 ∧ qf (M v) x < lam * (x ⬝ᵥ x) := by
  obtain ⟨W, hKW⟩ := Submodule.exists_isCompl (LinearMap.ker (A : V →ₗ[ℝ] Sq3))
  -- lower bound for A on the complement W (this is c_A > 0 of PROOF.md)
  have hinj : Function.Injective ((A : V →ₗ[ℝ] Sq3).comp W.subtype) := by
    rw [← LinearMap.ker_eq_bot, LinearMap.ker_eq_bot']
    intro w hw
    have hwK : (w : V) ∈ LinearMap.ker (A : V →ₗ[ℝ] Sq3) := by simpa using hw
    have : (w : V) ∈ LinearMap.ker (A : V →ₗ[ℝ] Sq3) ⊓ W := ⟨hwK, w.2⟩
    rw [hKW.inf_eq_bot, Submodule.mem_bot] at this
    exact Subtype.ext this
  obtain ⟨L, hL0, hL⟩ := (LinearMap.injective_iff_antilipschitz _).1 hinj
  have hLpos : (0 : ℝ) < L := by exact_mod_cast hL0
  have hAW : ∀ w ∈ W, ‖w‖ ≤ L * ‖A w‖ := by
    intro w hw
    have := hL.le_mul_dist (⟨w, hw⟩ : W) 0
    simpa [dist_eq_norm] using this
  obtain ⟨a, ha, hapos⟩ : ∃ a : ℝ, a = 1 / L ∧ 0 < a := ⟨1 / L, rfl, by positivity⟩
  obtain ⟨β, hβ, hβ0⟩ : ∃ β : ℝ, β = ‖B‖ ∧ 0 ≤ β := ⟨‖B‖, rfl, norm_nonneg B⟩
  have hB : ∀ u w : V, ‖B u w‖ ≤ β * ‖u‖ * ‖w‖ := fun u w => by rw [hβ]; exact B.le_opNorm₂ u w
  obtain ⟨C, hC, hCpos⟩ : ∃ C : ℝ, C = 36 * (β + 1) / a + 1 ∧ 0 < C :=
    ⟨_, rfl, by positivity⟩
  have haC : a * C = 36 * (β + 1) + a := by rw [hC]; field_simp
  obtain ⟨D, hD, hD0⟩ : ∃ D : ℝ, D = 2 * κ * C + 12 * β * C + 3 * β * C ^ 2 ∧ 0 ≤ D :=
    ⟨_, rfl, by positivity⟩
  obtain ⟨ε, hε, hεpos⟩ : ∃ ε : ℝ, ε = min 1 (κ / 12) ∧ 0 < ε := ⟨_, rfl, by positivity⟩
  have hεle1 : ε ≤ 1 := hε ▸ min_le_left _ _
  have hεleκ : ε ≤ κ / 12 := hε ▸ min_le_right _ _
  obtain ⟨ρ, hρ, hρpos⟩ : ∃ ρ : ℝ, ρ = min (min 1 (1 / (C + 1))) (κ / (4 * (D + 1))) ∧ 0 < ρ :=
    ⟨_, rfl, by positivity⟩
  filter_upwards [hrem ε hεpos, Metric.ball_mem_nhds (0 : V) hρpos] with v hv hball hv0
  have hrpos : 0 < ‖v‖ := norm_pos_iff.2 hv0
  have hrρ : ‖v‖ < ρ := by simpa using hball
  have hr1 : ‖v‖ < 1 := lt_of_lt_of_le hrρ (hρ ▸ (min_le_left _ _).trans (min_le_left _ _))
  have hrC : ‖v‖ < 1 / (C + 1) :=
    lt_of_lt_of_le hrρ (hρ ▸ (min_le_left _ _).trans (min_le_right _ _))
  have hrD : ‖v‖ < κ / (4 * (D + 1)) := lt_of_lt_of_le hrρ (hρ ▸ min_le_right _ _)
  have hCr : C * ‖v‖ ≤ 1 := by
    rw [lt_div_iff₀ (by positivity)] at hrC; nlinarith
  have hDr : D * ‖v‖ ≤ κ / 4 := by
    rw [lt_div_iff₀ (by positivity)] at hrD; nlinarith
  -- decomposition v = k + w
  have hmem : v ∈ LinearMap.ker (A : V →ₗ[ℝ] Sq3) ⊔ W := by
    rw [hKW.sup_eq_top]; exact Submodule.mem_top
  obtain ⟨k, hk, w, hw, hvkw⟩ := Submodule.mem_sup.1 hmem
  have hAk : A k = 0 := by simpa using hk
  have hMv : M v - lam • one3 = A v + B v v + (M v - (lam • one3 + A v + B v v)) := by abel
  have hRle := hv
  generalize M v - (lam • one3 + A v + B v v) = R at hMv hRle
  by_cases hcase : C * ‖v‖ ^ 2 ≤ ‖w‖
  · -- case (i): far from K, first-order argument
    have hAv : A v = A w := by rw [← hvkw, map_add, hAk, zero_add]
    obtain ⟨x, hx1, hxneg⟩ := traceless_symm_exists_neg (A w) (hA_symm w) (hA_tr w)
    have hAwlow : a * ‖w‖ ≤ ‖A w‖ := by
      have := hAW w hw
      rw [ha, div_mul_eq_mul_div, one_mul, div_le_iff₀ hLpos]; linarith
    have hBvv : ‖B v v‖ ≤ β * ‖v‖ ^ 2 := by
      have := hB v v; nlinarith [norm_nonneg v]
    have e1 : qf (B v v) x ≤ 9 * (β * ‖v‖ ^ 2) := by
      have := (abs_le.1 (abs_qf_le (B v v) x hx1)).2; nlinarith
    have e2 : qf R x ≤ 9 * (ε * ‖v‖ ^ 2) := by
      have := (abs_le.1 (abs_qf_le R x hx1)).2; nlinarith
    have e3 : qf (A w) x ≤ -(a * ‖w‖) / 4 := by linarith
    have key : qf (M v - lam • one3) x < 0 := by
      rw [hMv, qf_add, qf_add, hAv]
      exact case1_arith hapos hrpos hβ0 hεle1 haC hcase e3 e1 e2
    refine ⟨x, ?_, ?_⟩
    · rintro rfl; simp [qf] at key
    · rw [qf_sub, qf_smul_one] at key; linarith
  · -- case (ii): near K, trace argument (uses tr A ≡ 0 essentially)
    push Not at hcase
    have hwle : ‖w‖ ≤ C * ‖v‖ ^ 2 := hcase.le
    have hkv : k = v - w := by rw [← hvkw]; abel
    have hklow : ‖v‖ - C * ‖v‖ ^ 2 ≤ ‖k‖ := by
      have := norm_sub_norm_le v w
      rw [hkv]; linarith
    have hkup : ‖k‖ ≤ 2 * ‖v‖ := by
      have := norm_sub_le v w
      have : C * ‖v‖ ^ 2 ≤ ‖v‖ := by
        have hx := mul_le_mul_of_nonneg_right hCr hrpos.le
        have hy : C * ‖v‖ ^ 2 = C * ‖v‖ * ‖v‖ := by ring
        linarith
      rw [hkv]; linarith
    have hBexp : B v v = B k k + B k w + B w k + B w w := by
      rw [← hvkw]; simp only [map_add, _root_.add_apply]; abel
    have t2 : tr3 (B k w) ≤ 3 * (β * ‖k‖ * ‖w‖) :=
      (le_abs_self _).trans ((abs_tr3_le _).trans (by linarith [hB k w]))
    have t3 : tr3 (B w k) ≤ 3 * (β * ‖w‖ * ‖k‖) :=
      (le_abs_self _).trans ((abs_tr3_le _).trans (by linarith [hB w k]))
    have t4 : tr3 (B w w) ≤ 3 * (β * ‖w‖ * ‖w‖) :=
      (le_abs_self _).trans ((abs_tr3_le _).trans (by linarith [hB w w]))
    have t5 : tr3 R ≤ 3 * (ε * ‖v‖ ^ 2) :=
      (le_abs_self _).trans ((abs_tr3_le _).trans (by linarith))
    have htr : tr3 (M v - lam • one3) < 0 := by
      rw [hMv, tr3_add, tr3_add, hA_tr, hBexp, tr3_add, tr3_add, tr3_add, zero_add]
      exact case2_arith hκ hrpos hr1.le hCr hCpos hβ0 hεleκ hwle (norm_nonneg w) hklow hkup
        (hneg k hk) t2 t3 t4 t5 hD hDr
    -- some diagonal entry of M v - λ is < 0
    unfold tr3 at htr
    have : ∃ i : Fin 3, (M v - lam • one3) i i < 0 := by
      by_contra h; push Not at h
      linarith [h 0, h 1, h 2]
    obtain ⟨i, hi⟩ := this
    refine ⟨Pi.single i 1, ?_, ?_⟩
    · intro h; have := congrFun h i; simp at this
    · have e1 : qf (M v) (Pi.single i 1) = M v i i := by
        simp [qf, Pi.single_apply]
      have e2 : (Pi.single i (1 : ℝ) : Fin 3 → ℝ) ⬝ᵥ (Pi.single i 1) = 1 := by
        simp [dotProduct_single]
      rw [e1, e2]
      simp [one3] at hi
      linarith

end criterion

/-- If `xᵀ S x < c·xᵀx` for some `x` and `S` is symmetric, some eigenvalue of `S` is `< c`.
(So the conclusion of `local_max_criterion` is `λ_min (M v) < λ*`.) -/
lemma exists_eigenvalue_lt (S0 : Sq3) (hS : (Matrix.of S0).IsHermitian) (c : ℝ)
    (h : ∃ x : Fin 3 → ℝ, qf S0 x < c * (x ⬝ᵥ x)) :
    ∃ i, hS.eigenvalues i < c := by
  by_contra hcon
  push Not at hcon
  obtain ⟨x, hx⟩ := h
  rw [qf_eq_matrix] at hx
  set S : Mat3 := Matrix.of S0
  have halg : algebraMap ℝ Mat3 c = c • (1 : Mat3) := Algebra.algebraMap_eq_smul_one c
  have hH : (S - c • (1 : Mat3)).IsHermitian := by
    refine hS.sub ?_
    exact (isHermitian_one).smul (by simp [IsSelfAdjoint])
  have hpsd : (S - c • (1 : Mat3)).PosSemidef := by
    rw [hH.posSemidef_iff_eigenvalues_nonneg]
    intro i
    have hmem := hH.eigenvalues_mem_spectrum_real i
    set μ := hH.eigenvalues i
    have : μ + c ∈ spectrum ℝ S := by
      rw [spectrum.mem_iff] at hmem ⊢
      convert hmem using 2
      rw [← halg, map_add]; abel
    rw [hS.spectrum_real_eq_range_eigenvalues] at this
    obtain ⟨j, hj⟩ := this
    have := hcon j
    show 0 ≤ μ
    linarith
  have h2 := hpsd.dotProduct_mulVec_nonneg x
  simp [sub_mulVec, dotProduct_sub, smul_mulVec, dotProduct_smul] at h2
  linarith

end Bolza
