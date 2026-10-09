/-
# Ky Fan step of Lemma 8: (Λ∘x)″(0) ≤ g″(0)  (PROVED IN LEAN)

PROOF.md Lemma 8 step 1: Ky Fan gives Λ(x(p)) ≤ g(p) with equality at p = 0 (g = trace of the
compressed pencil on a 4-dimensional trial space).  The analytic consequence used is
"a function bounded above by g, with equality at 0, has second derivative ≤ g″(0)".
That is `second_deriv_le_of_le` below.

The Ky Fan inequality itself (Λ ≤ tr_V(q|m) for 4-dimensional V ⊂ H¹, an infinite-dimensional
min-max statement, D9) is NOT formalized; it enters only as a hypothesis.  We do prove the
finite-dimensional *bathtub* lemma that is the combinatorial heart of finite-dimensional Ky Fan
(`bathtub`): for ascending μ and weights 0 ≤ d ≤ 1 with Σ d = k, Σ μᵢdᵢ ≥ Σ_{i<k} μᵢ.
-/
import Mathlib

open Filter Topology Set Matrix

namespace Bolza

/-- If `h ≥ 0` near 0, `h 0 = 0`, `h` is differentiable near 0 and `h'` is differentiable at 0,
then `h''(0) ≥ 0`. -/
theorem second_deriv_nonneg_of_min {h : ℝ → ℝ} (hpos : ∀ᶠ t in 𝓝 0, 0 ≤ h t) (h0 : h 0 = 0)
    (hd : ∀ᶠ t in 𝓝 0, DifferentiableAt ℝ h t) (hd2 : DifferentiableAt ℝ (deriv h) 0) :
    0 ≤ deriv (deriv h) 0 := by
  have hmin : IsLocalMin h 0 := by
    filter_upwards [hpos] with t ht
    rw [h0]; exact ht
  have h'0 : deriv h 0 = 0 := hmin.deriv_eq_zero
  by_contra hneg
  push Not at hneg
  set c := deriv (deriv h) 0
  -- slope of h' at 0 tends to c < 0 from the right
  have hslope := hd2.hasDerivAt.tendsto_slope_zero_right
  have hev1 : ∀ᶠ t in 𝓝[>] (0 : ℝ), t⁻¹ • (deriv h (0 + t) - deriv h 0) < 0 :=
    hslope (Iio_mem_nhds hneg)
  have hev2 : ∀ᶠ t in 𝓝[>] (0 : ℝ), DifferentiableAt ℝ h t ∧ 0 ≤ h t :=
    nhdsWithin_le_nhds (hd.and hpos)
  obtain ⟨ε, hε, hεsub⟩ := mem_nhdsGT_iff_exists_Ioo_subset.1 (hev1.and hev2)
  have hε0 : 0 < ε := hε
  set t := ε / 2
  have ht : t ∈ Ioo 0 ε := ⟨by positivity, by simp only [t]; linarith⟩
  have hcont : ContinuousOn h (Icc 0 t) := by
    intro s hs
    rcases eq_or_lt_of_le hs.1 with hs0 | hs0
    · subst hs0
      exact (hd.self_of_nhds).continuousAt.continuousWithinAt
    · exact (hεsub ⟨hs0, lt_of_le_of_lt hs.2 ht.2⟩).2.1.continuousAt.continuousWithinAt
  have hderiv : ∀ s ∈ Ioo 0 t, HasDerivAt h (deriv h s) s := fun s hs =>
    (hεsub ⟨hs.1, hs.2.trans ht.2⟩).2.1.hasDerivAt
  obtain ⟨ξ, hξ, hξeq⟩ := exists_hasDerivAt_eq_slope h (deriv h) ht.1 hcont hderiv
  have hξneg : deriv h ξ < 0 := by
    have := (hεsub ⟨hξ.1, hξ.2.trans ht.2⟩).1
    simp only [zero_add, h'0, sub_zero, smul_eq_mul] at this
    have hξpos : 0 < ξ := hξ.1
    by_contra hc; push Not at hc
    have : 0 ≤ ξ⁻¹ * deriv h ξ := mul_nonneg (inv_nonneg.2 hξpos.le) hc
    linarith
  have htpos := (hεsub ht).2.2
  rw [h0, sub_zero, sub_zero] at hξeq
  have hq : h t / t < 0 := hξeq ▸ hξneg
  have hnn : 0 ≤ h t / t := div_nonneg htpos ht.1.le
  linarith

/-- **Ky Fan → (Λ∘x)″ ≤ g″ (PROOF.md Lemma 8 step 1).**  If `f ≤ g` near 0 with equality at 0,
both differentiable near 0 with derivatives differentiable at 0, then `f''(0) ≤ g''(0)`. -/
theorem second_deriv_le_of_le {f g : ℝ → ℝ} (hle : ∀ᶠ t in 𝓝 0, f t ≤ g t) (h0 : f 0 = g 0)
    (hf : ∀ᶠ t in 𝓝 0, DifferentiableAt ℝ f t) (hg : ∀ᶠ t in 𝓝 0, DifferentiableAt ℝ g t)
    (hf2 : DifferentiableAt ℝ (deriv f) 0) (hg2 : DifferentiableAt ℝ (deriv g) 0) :
    deriv (deriv f) 0 ≤ deriv (deriv g) 0 := by
  set h := fun t => g t - f t
  have hderiv_eq : deriv h =ᶠ[𝓝 0] fun t => deriv g t - deriv f t := by
    filter_upwards [hf, hg] with t hft hgt
    exact (hgt.hasDerivAt.sub hft.hasDerivAt).deriv
  have hd2 : DifferentiableAt ℝ (deriv h) 0 :=
    (hg2.sub hf2).congr_of_eventuallyEq hderiv_eq
  have hdd : deriv (deriv h) 0 = deriv (deriv g) 0 - deriv (deriv f) 0 := by
    rw [hderiv_eq.deriv_eq]
    exact deriv_sub hg2 hf2
  have := second_deriv_nonneg_of_min (h := h)
    (by filter_upwards [hle] with t ht; simp [h]; linarith)
    (by simp [h, h0])
    (by filter_upwards [hf, hg] with t hft hgt; exact hgt.sub hft)
    hd2
  linarith

/-- **Bathtub lemma** (combinatorial core of finite-dimensional Ky Fan): if `μ` is ascending on
`Fin n`, `0 ≤ dᵢ ≤ 1` and `Σ dᵢ = k`, then `Σ_{i<k} μᵢ ≤ Σ μᵢ dᵢ`. -/
theorem bathtub {n : ℕ} (μ d : Fin n → ℝ) (hμ : Monotone μ) (hd0 : ∀ i, 0 ≤ d i)
    (hd1 : ∀ i, d i ≤ 1) (k : ℕ) (hk : ∑ i, d i = k) :
    ∑ i ∈ Finset.univ.filter (fun i : Fin n => (i : ℕ) < k), μ i ≤ ∑ i, μ i * d i := by
  classical
  set L := Finset.univ.filter (fun i : Fin n => (i : ℕ) < k)
  -- a threshold value m with μ i ≤ m on L and μ i ≥ m off L
  by_cases hkn : k < n
  · set m := μ ⟨k, hkn⟩
    have hL : ∀ i ∈ L, μ i ≤ m := by
      intro i hi
      simp only [L, Finset.mem_filter] at hi
      exact hμ (show i ≤ ⟨k, hkn⟩ from Fin.le_iff_val_le_val.2 hi.2.le)
    have hLc : ∀ i ∉ L, m ≤ μ i := by
      intro i hi
      simp only [L, Finset.mem_filter, Finset.mem_univ, true_and, not_lt] at hi
      exact hμ (Fin.le_iff_val_le_val.2 hi)
    have hcard : (L.card : ℝ) = k := by
      have : L.card = k := by
        simp only [L]
        rw [Fin.card_filter_val_lt]
        omega
      exact_mod_cast this
    -- Σ μ d - Σ_L μ = Σ_L μ (d - 1) + Σ_{Lᶜ} μ d ≥ m (Σ_L (d-1) + Σ_{Lᶜ} d) = m (Σ d - k) = 0
    have key : ∀ i, μ i * d i - (if i ∈ L then μ i else 0) ≥
        m * d i - (if i ∈ L then m else 0) := by
      intro i
      by_cases hi : i ∈ L
      · simp only [hi, ite_true]
        have := hL i hi
        nlinarith [hd1 i]
      · simp only [hi, ite_false, sub_zero]
        exact mul_le_mul_of_nonneg_right (hLc i hi) (hd0 i)
    have hsum := Finset.sum_le_sum (fun i (_ : i ∈ Finset.univ) => key i)
    simp only [Finset.sum_sub_distrib, ← Finset.mul_sum, Finset.sum_ite_mem, Finset.univ_inter,
      Finset.sum_const, nsmul_eq_mul] at hsum
    rw [hk, hcard] at hsum
    linarith
  · -- k ≥ n: then Σ d = k ≤ ... forces all d = 1 when k = n; L = univ
    push Not at hkn
    have hL : L = Finset.univ := by
      ext i; simp [L]; omega
    have hsumd : ∑ i, d i ≤ n := by
      calc ∑ i, d i ≤ ∑ _i : Fin n, (1 : ℝ) := Finset.sum_le_sum fun i _ => hd1 i
        _ = n := by simp
    have hkn' : (k : ℝ) = n := by
      have : (n : ℝ) ≤ k := by exact_mod_cast hkn
      linarith [hk ▸ hsumd]
    -- all d i = 1
    have hall : ∀ i, d i = 1 := by
      intro i
      by_contra hne
      have hlt : d i < 1 := lt_of_le_of_ne (hd1 i) hne
      have : ∑ j, d j < ∑ _j : Fin n, (1 : ℝ) :=
        Finset.sum_lt_sum (fun j _ => hd1 j) ⟨i, Finset.mem_univ _, hlt⟩
      simp at this
      linarith
    rw [hL]
    simp [hall]

/-- Bathtub lemma, unsorted form: some `k`-element index set `T` has `Σ_T μ ≤ Σ μᵢdᵢ`. -/
theorem bathtub_exists {n : ℕ} (μ d : Fin n → ℝ) (hd0 : ∀ i, 0 ≤ d i) (hd1 : ∀ i, d i ≤ 1)
    (k : ℕ) (hk : ∑ i, d i = k) :
    ∃ T : Finset (Fin n), T.card = k ∧ ∑ i ∈ T, μ i ≤ ∑ i, μ i * d i := by
  classical
  set σ := Tuple.sort μ
  have hmono : Monotone (μ ∘ σ) := Tuple.monotone_sort μ
  have hk' : ∑ i, (d ∘ σ) i = k := by rw [← hk]; exact Equiv.sum_comp σ d
  have hb := bathtub (μ ∘ σ) (d ∘ σ) hmono (fun i => hd0 _) (fun i => hd1 _) k hk'
  have hkn : k ≤ n := by
    have : (k : ℝ) ≤ n := by
      rw [← hk]
      calc ∑ i, d i ≤ ∑ _i : Fin n, (1 : ℝ) := Finset.sum_le_sum fun i _ => hd1 i
        _ = n := by simp
    exact_mod_cast this
  refine ⟨(Finset.univ.filter (fun i : Fin n => (i : ℕ) < k)).map σ.toEmbedding, ?_, ?_⟩
  · rw [Finset.card_map, Fin.card_filter_val_lt]; omega
  · rw [Finset.sum_map]
    exact le_of_le_of_eq hb (Equiv.sum_comp σ (fun i => μ i * d i))

/-- **Finite-dimensional Ky Fan inequality (PROVED IN LEAN).**  For a real symmetric `n×n`
matrix `S` and `V` with orthonormal columns (`Vᵀ V = 1`, `n×k`), the trace of the compression
`Vᵀ S V` is at least the sum of some `k` eigenvalues of `S` (hence at least the sum of the `k`
smallest). -/
theorem ky_fan {n k : ℕ} (S : Matrix (Fin n) (Fin n) ℝ) (hS : S.IsHermitian)
    (V : Matrix (Fin n) (Fin k) ℝ) (hV : Vᵀ * V = 1) :
    ∃ T : Finset (Fin n), T.card = k ∧
      ∑ i ∈ T, hS.eigenvalues i ≤ Matrix.trace (Vᵀ * S * V) := by
  classical
  set μ := hS.eigenvalues
  set W : Matrix (Fin n) (Fin n) ℝ := (hS.eigenvectorUnitary : Matrix (Fin n) (Fin n) ℝ)
  have hstar : star W = Wᵀ := by
    rw [Matrix.star_eq_conjTranspose, Matrix.conjTranspose_eq_transpose_of_trivial]
  have hWW : W * Wᵀ = 1 := by
    rw [← hstar]; exact Matrix.mem_unitaryGroup_iff.1 hS.eigenvectorUnitary.2
  have hspec : S = W * Matrix.diagonal μ * Wᵀ := by
    conv_lhs => rw [hS.spectral_theorem]
    rw [Unitary.conjStarAlgAut_apply, ← hstar]
    congr 2
  set Y := Wᵀ * V
  have hYY : Yᵀ * Y = 1 := by
    simp only [Y, Matrix.transpose_mul, Matrix.transpose_transpose]
    rw [Matrix.mul_assoc, ← Matrix.mul_assoc W, hWW, Matrix.one_mul, hV]
  set d : Fin n → ℝ := fun j => ∑ l, Y j l * Y j l
  have htr : Matrix.trace (Vᵀ * S * V) = ∑ j, μ j * d j := by
    have : Vᵀ * S * V = Yᵀ * Matrix.diagonal μ * Y := by
      rw [hspec]; simp only [Y, Matrix.transpose_mul, Matrix.transpose_transpose, Matrix.mul_assoc]
    rw [this]
    simp only [Matrix.trace, Matrix.diag, Matrix.mul_apply, Matrix.diagonal_apply,
      Matrix.transpose_apply, d]
    simp only [mul_ite, mul_zero, Finset.sum_ite_eq', Finset.mem_univ, ite_true]
    rw [Finset.sum_comm]
    refine Finset.sum_congr rfl fun j _ => ?_
    rw [Finset.mul_sum]
    refine Finset.sum_congr rfl fun l _ => ?_
    ring
  have hd0 : ∀ j, 0 ≤ d j := fun j => Finset.sum_nonneg fun l _ => mul_self_nonneg _
  have hsum : ∑ j, d j = k := by
    have : ∑ j, d j = Matrix.trace (Yᵀ * Y) := by
      simp only [d, Matrix.trace, Matrix.diag, Matrix.mul_apply, Matrix.transpose_apply]
      rw [Finset.sum_comm]
    rw [this, hYY, Matrix.trace_one, Fintype.card_fin]
  have hd1 : ∀ j, d j ≤ 1 := by
    intro j
    set P := Y * Yᵀ
    have hPP : P * P = P := by
      simp only [P]; rw [Matrix.mul_assoc, ← Matrix.mul_assoc Yᵀ, hYY, Matrix.one_mul]
    have hPjj : P j j = d j := by simp [P, Matrix.mul_apply, d]
    have hPsym : ∀ a b, P a b = P b a := by
      intro a b; simp only [P, Matrix.mul_apply, Matrix.transpose_apply]
      exact Finset.sum_congr rfl fun l _ => mul_comm _ _
    have h2 : P j j = ∑ m, P j m * P j m := by
      have := congrFun (congrFun hPP j) j
      rw [Matrix.mul_apply] at this
      rw [← this]
      exact Finset.sum_congr rfl fun m _ => by rw [hPsym m j]
    have h3 : P j j * P j j ≤ ∑ m, P j m * P j m :=
      Finset.single_le_sum (f := fun m => P j m * P j m) (fun m _ => mul_self_nonneg _)
        (Finset.mem_univ j)
    rw [hPjj] at h2 h3
    nlinarith [hd0 j]
  obtain ⟨T, hT, hTle⟩ := bathtub_exists μ d hd0 hd1 k hsum
  exact ⟨T, hT, htr ▸ hTle⟩

end Bolza
