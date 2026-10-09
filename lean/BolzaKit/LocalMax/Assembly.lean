/-
# Assembly (PROOF.md "Theorem: assembly"), conditional on explicitly named inputs

Everything that is *not* finite-dimensional algebra/analysis is packaged as a field of the
structure `BolzaInputs`; each field's docstring cites the PROOF.md lemma / dependency it stands
for.  The theorem `bolza_strict_local_max` is then PROVED IN LEAN from those fields, using
* Lemma 10 in C² form (`local_max_criterion_C2`, Taylor + Weyl/trace argument),
* the real Schur-type lemmas of `RepReal.lean` (Lemmas 4, 6(iv)–(vi), 7),
* nothing else.

Coordinates.  `TB = ℝ³ × ℝ³` are the block coordinates of T_B ≅ TT (the chart S of Lemma 9
transported by the isomorphism σ of Lemma 9d): first factor T₋ = K₀, second factor T₊ = W₀.
`RepTheory.lean` checks (exactly, by kernel evaluation) that the change of basis `P` conjugates
the cert_chars.py matrices of r, t, κ, a to the block matrices
(Xr,Xr), (XmT,XpT), (XmK,XpK), (−1,+1) used below.

The conclusion "λ₁(X) < λ₁(B) for X ≠ B near B" is exactly PROOF.md §0 Theorem.
-/
import BolzaKit.LocalMax.Taylor
import BolzaKit.LocalMax.RepReal
import BolzaKit.LocalMax.KyFan

open Matrix Filter Topology

namespace Bolza

/-- block coordinates of T_B: K₀ × W₀ -/
abbrev TB := (Fin 3 → ℝ) × (Fin 3 → ℝ)

/-- block-diagonal action -/
def act (X Y : Mat3) (v : TB) : TB := (X *ᵥ v.1, Y *ᵥ v.2)

/-- The external inputs.  Each field is a hypothesis of the main theorem; none is proved in Lean. -/
structure BolzaInputs where
  /-- Teichmüller space T₂ (D12: a topological 6-manifold). -/
  T2 : Type
  /-- its topology -/
  topT2 : TopologicalSpace T2
  /-- the Bolza surface -/
  B : T2
  /-- λ₁ on T₂ -/
  lam1 : T2 → ℝ
  /-- the chart π|_S of Lemma 9, written in TT block coordinates (σ, Lemma 9d) -/
  chart : TB → T2
  /-- Lemma 9(c): the chart sends 0 to B -/
  chart_zero : chart 0 = B
  /-- Lemma 9(c): π|_S is a homeomorphism onto an open neighbourhood of B (invariance of domain,
  certified rank dπ = 6 of CERT_G3); we only use that images of neighbourhoods of 0 are
  neighbourhoods of B. -/
  chart_nhds : ∀ s ∈ 𝓝 (0 : TB), chart '' s ∈ @nhds T2 topT2 B
  /-- Lemma 2 (F1) on the slice S: the 3×3 cluster matrix M_S (in TT block coordinates) -/
  M : TB → Sq3
  /-- λ* = λ₁(B) -/
  lamStar : ℝ
  /-- Lemma 1(c): λ₁(B) = λ* (the multiplicity-3 cluster value) -/
  lam1_B : lam1 B = lamStar
  /-- Lemma 2 (F1): M is real-analytic on S, in particular C² at 0 -/
  M_C2 : ContDiffAt ℝ 2 M 0
  /-- Lemma 2 (F1) item 3: M is symmetric -/
  M_sym : ∀ v, IsSym3 (M v)
  /-- Lemma 2 (F1) item 4: M(0) = λ* I -/
  M_zero : M 0 = lamStar • one3
  /-- Lemma 2 (F1) items 1, 3 + Lemma 1 (isolation, λ₀ = 0 < λ₁): near B, λ₁(π x) is the smallest
  eigenvalue of M(x), hence ≤ every Rayleigh quotient of M(x) (min-max, D9). -/
  lam1_le_rayleigh : ∀ᶠ v in 𝓝 (0 : TB), ∀ x : Fin 3 → ℝ,
    lam1 (chart v) * (x ⬝ᵥ x) ≤ qf (M v) x
  /-- ρ_E(r) (E = λ₁-eigenspace, m₀-orthonormal basis) -/
  ρr : Mat3
  /-- ρ_E(t) -/
  ρt : Mat3
  /-- ρ_E(κ) -/
  ρk : Mat3
  /-- ρ_E(a), a the antipodal isometry -/
  ρa : Mat3
  /-- ρ_E(J), J the hyperelliptic involution -/
  ρJ : Mat3
  /-- isometries act orthogonally on (E, m₀) (only needed for r) -/
  ρr_orth : ρrᵀ * ρr = 1
  /-- Lemma 3(b) (Lemma N): G-equivariance A(g·σ) = ρ_E(g) A(σ) ρ_E(g)ᵀ, for g = r -/
  equiv_r : ∀ v, Matrix.of (fderiv ℝ M 0 (act Xr Xr v)) = ρr * Matrix.of (fderiv ℝ M 0 v) * ρrᵀ
  /-- Lemma 3(b) for g = t -/
  equiv_t : ∀ v, Matrix.of (fderiv ℝ M 0 (act XmT XpT v)) = ρt * Matrix.of (fderiv ℝ M 0 v) * ρtᵀ
  /-- Lemma 3(b) for g = κ (= s_mid) -/
  equiv_k : ∀ v, Matrix.of (fderiv ℝ M 0 (act XmK XpK v)) = ρk * Matrix.of (fderiv ℝ M 0 v) * ρkᵀ
  /-- Lemma 3(b) for g = a (block matrix (−1,+1), checked in `RepTheory.conj_a`) -/
  equiv_a : ∀ v, Matrix.of (fderiv ℝ M 0 (act (-1) 1 v)) = ρa * Matrix.of (fderiv ℝ M 0 v) * ρaᵀ
  /-- Lemma 1(c) (via D4 = FBP Lemmas 7.2–7.3): E is a real-irreducible G-module -/
  E_irred : Irreducible3 [ρr, ρt, ρk]
  /-- group relation in G (D5; exact in cert_chars.py): J² = 1 -/
  J_sq : ρJ * ρJ = 1
  /-- group relation in G (D5; exact in cert_chars.py): J is central -/
  J_central : ∀ g ∈ [ρr, ρt, ρk], ρJ * g = g * ρJ
  /-- group relation in G (Lemma 6(iii); class {a, a³ = aJ}, Cook class 4a): a² = J -/
  a_sq : ρa * ρa = ρJ
  /-- group relation in G (Lemma 6(iii)): a is central modulo J -/
  a_central_modJ : ∀ g ∈ [ρr, ρt, ρk], g * ρa = ρa * g ∨ g * ρa = ρa * ρJ * g
  /-- Lemma 5 (CERTIFIED, cert_first_order.py): A ≠ 0 -/
  A_ne_zero : fderiv ℝ M 0 ≠ 0
  /-- Lemma 11 (OPEN-1 closure): D²Λ = tr D²M is invariant under r -/
  hess_inv_r : ∀ v w, tr3 (fderiv ℝ (fderiv ℝ M) 0 (act Xr Xr v) (act Xr Xr w)) =
    tr3 (fderiv ℝ (fderiv ℝ M) 0 v w)
  /-- Lemma 11: invariance under t -/
  hess_inv_t : ∀ v w, tr3 (fderiv ℝ (fderiv ℝ M) 0 (act XmT XpT v) (act XmT XpT w)) =
    tr3 (fderiv ℝ (fderiv ℝ M) 0 v w)
  /-- Lemma 11: invariance under κ -/
  hess_inv_k : ∀ v w, tr3 (fderiv ℝ (fderiv ℝ M) 0 (act XmK XpK v) (act XmK XpK w)) =
    tr3 (fderiv ℝ (fderiv ℝ M) 0 v w)
  /-- Lemma 7 + Lemma 9(d)3: w_B = σ(dp v_B), the TT tangent of the family F_B -/
  wB : TB
  /-- Lemma 7 (hyperbolic trigonometry dℓ(g₀)/dp = −√(2+2√2) ≠ 0): w_B ≠ 0 -/
  wB_ne : wB ≠ 0
  /-- Lemma 3 relabelling line: r maps F_B(p) to F_B(−p), so r·w_B = −w_B -/
  wB_r : act Xr Xr wB = -wB
  /-- Lemma 3 relabelling line: F_B is s_mid-invariant (s_mid = κ in cert_chars.py; any
  s_mid-type reflection gives the same line, `Rep.smid_line_in_Kminus`) -/
  wB_smid : act XmK XpK wB = wB
  /-- Lemma 9(d)1 + Lemma 7: the F_B family read in the chart S, γ = p ∘ x_{F_B} : ℝ → S,
  in TT block coordinates; γ(0) = 0 -/
  γ : ℝ → TB
  /-- γ(0) = x_B -/
  γ_zero : γ 0 = 0
  /-- γ is analytic (F_B is an analytic family, p analytic): differentiable near 0 -/
  γ_diff : ∀ᶠ t in 𝓝 0, DifferentiableAt ℝ γ t
  /-- … and twice differentiable at 0 -/
  γ_diff2 : DifferentiableAt ℝ (deriv γ) 0
  /-- Lemma 9(d)3: σ(dp v_B) = w_B, i.e. γ'(0) = w_B -/
  γ_tangent : deriv γ 0 = wB
  /-- Lemma 8 step 1–3: the Ky Fan comparison function g(p) = tr(G_m⁻¹G_q) of the P2 trial space -/
  g : ℝ → ℝ
  /-- g is analytic: differentiable near 0 -/
  g_diff : ∀ᶠ t in 𝓝 0, DifferentiableAt ℝ g t
  /-- … and twice differentiable at 0 -/
  g_diff2 : DifferentiableAt ℝ (deriv g) 0
  /-- Lemma 8 step 1 (Ky Fan / min-max, D9, with Λ = tr M on S by F1): Λ(γ(p)) ≤ g(p) near 0 -/
  kyfan_le : ∀ᶠ t in 𝓝 0, tr3 (M (γ t)) ≤ g t
  /-- Lemma 8 step 1: equality at p = 0 (exact eigenfunctions e_i ∈ E in the trial space) -/
  kyfan_eq : tr3 (M (γ 0)) = g 0
  /-- Lemma 8 (CERTIFIED): g″(0) ≤ −9.3411 (Lehmann–Goerisch route, primary certificate).  The
  SU13 route gives the stronger g″(0) ≤ −12.850476… ≤ −12.8504 (`Cert.sign_joint`), which implies
  this hypothesis (`Cert.sign_su13_implies_lg`); any negative constant suffices for the argument. -/
  g_sign : deriv (deriv g) 0 ≤ -(93411 / 10000 : ℝ)

/-! ### small helpers -/

/-- trace as a linear map -/
def tr3L : Sq3 →ₗ[ℝ] ℝ where
  toFun := tr3
  map_add' S T := tr3_add S T
  map_smul' c S := by simp [tr3]; ring

lemma tr3_smul (c : ℝ) (S : Sq3) : tr3 (c • S) = c • tr3 S := by
  simp [tr3]; ring

lemma tr3_eq_trace (S : Sq3) : tr3 S = Matrix.trace (Matrix.of S) := by
  simp [tr3, Matrix.trace, Fin.sum_univ_three]

lemma tr3_conj_orth (S : Sq3) (Q : Mat3) (hQ : Qᵀ * Q = 1) (S' : Sq3)
    (h : Matrix.of S' = Q * Matrix.of S * Qᵀ) : tr3 S' = tr3 S := by
  rw [tr3_eq_trace, tr3_eq_trace, h, Matrix.trace_mul_comm, ← Matrix.mul_assoc, hQ,
    Matrix.one_mul]

lemma eq_zero_of_conj_zero (S S' : Sq3) (Q : Mat3) (h : Matrix.of S' = Q * Matrix.of S * Qᵀ)
    (hS : S = 0) : S' = 0 := by
  have : Matrix.of S' = 0 := by rw [h, hS]; simp
  funext i j; exact congrFun (congrFun this i) j

lemma dot_self_pos {x : Fin 3 → ℝ} (hx : x ≠ 0) : 0 < x ⬝ᵥ x := by
  simp only [dotProduct, Fin.sum_univ_three]
  by_contra h
  push Not at h
  apply hx
  have h0 : x 0 = 0 := by nlinarith [sq_nonneg (x 0), sq_nonneg (x 1), sq_nonneg (x 2)]
  have h1 : x 1 = 0 := by nlinarith [sq_nonneg (x 0), sq_nonneg (x 1), sq_nonneg (x 2)]
  have h2 : x 2 = 0 := by nlinarith [sq_nonneg (x 0), sq_nonneg (x 1), sq_nonneg (x 2)]
  ext i; fin_cases i <;> simp [h0, h1, h2]

lemma norm_sq_le_dot (y : Fin 3 → ℝ) : ‖y‖ ^ 2 ≤ y ⬝ᵥ y := by
  have hs : 0 ≤ y ⬝ᵥ y := by
    simp only [dotProduct, Fin.sum_univ_three]
    nlinarith [mul_self_nonneg (y 0), mul_self_nonneg (y 1), mul_self_nonneg (y 2)]
  have : ‖y‖ ≤ Real.sqrt (y ⬝ᵥ y) := by
    refine (pi_norm_le_iff_of_nonneg (Real.sqrt_nonneg _)).2 fun i => ?_
    rw [Real.norm_eq_abs]
    apply Real.abs_le_sqrt
    simp only [dotProduct, Fin.sum_univ_three]
    fin_cases i <;> simp <;> nlinarith [sq_nonneg (y 0), sq_nonneg (y 1), sq_nonneg (y 2)]
  calc ‖y‖ ^ 2 ≤ (Real.sqrt (y ⬝ᵥ y)) ^ 2 := by gcongr
    _ = y ⬝ᵥ y := Real.sq_sqrt hs

lemma act_inl (X Y : Mat3) (y : Fin 3 → ℝ) : act X Y (y, 0) = (X *ᵥ y, 0) := by
  simp [act]

lemma act_inr (X Y : Mat3) (z : Fin 3 → ℝ) : act X Y (0, z) = (0, Y *ᵥ z) := by
  simp [act]

/-- trace as a continuous linear map -/
noncomputable def tr3CLM : Sq3 →L[ℝ] ℝ := LinearMap.toContinuousLinearMap tr3L

lemma tr3CLM_apply (S : Sq3) : tr3CLM S = tr3 S := rfl

/-- **Second-order chain rule at a critical point (PROOF.md Lemma 8 "Consequence").  PROVED.**
If `M` is C² at 0 with `tr dM(0) = 0` and `γ(0) = 0` is twice differentiable at 0, then
`(tr M ∘ γ)″(0) = tr D²M(0)(γ′(0), γ′(0))`. -/
lemma second_deriv_trace_comp (M : TB → Sq3) (hM : ContDiffAt ℝ 2 M 0) (γ : ℝ → TB)
    (hγ0 : γ 0 = 0) (hγd : ∀ᶠ t in 𝓝 0, DifferentiableAt ℝ γ t)
    (hγ2 : DifferentiableAt ℝ (deriv γ) 0) (htr : ∀ v, tr3 (fderiv ℝ M 0 v) = 0) :
    (∀ᶠ t in 𝓝 0, DifferentiableAt ℝ (fun t => tr3 (M (γ t))) t) ∧
    DifferentiableAt ℝ (deriv (fun t => tr3 (M (γ t)))) 0 ∧
    deriv (deriv (fun t => tr3 (M (γ t)))) 0 =
      tr3 (fderiv ℝ (fderiv ℝ M) 0 (deriv γ 0) (deriv γ 0)) := by
  have hMd : ∀ᶠ y in 𝓝 (0 : TB), HasFDerivAt M (fderiv ℝ M y) y := by
    filter_upwards [hM.eventually (by simp)] with y hy
    exact (hy.differentiableAt (by norm_num)).hasFDerivAt
  have hγc : Filter.Tendsto γ (𝓝 0) (𝓝 0) := by
    have := (hγd.self_of_nhds).continuousAt.tendsto
    rwa [hγ0] at this
  have hMγ : ∀ᶠ t in 𝓝 (0 : ℝ), HasFDerivAt M (fderiv ℝ M (γ t)) (γ t) := hγc.eventually hMd
  have hφ : ∀ᶠ t in 𝓝 (0 : ℝ), HasDerivAt (fun t => tr3 (M (γ t)))
      (tr3 (fderiv ℝ M (γ t) (deriv γ t))) t := by
    filter_upwards [hMγ, hγd] with t h1 h2
    have := tr3CLM.hasFDerivAt.comp_hasDerivAt t (h1.comp_hasDerivAt t h2.hasDerivAt)
    simpa [tr3CLM_apply, Function.comp_def] using this
  have hderiv_eq : deriv (fun t => tr3 (M (γ t))) =ᶠ[𝓝 0]
      fun t => tr3 (fderiv ℝ M (γ t) (deriv γ t)) := by
    filter_upwards [hφ] with t ht; exact ht.deriv
  -- derivative of ψ(t) = tr (dM(γ t) (γ' t)) at 0
  have hf1 : ContDiffAt ℝ 1 (fderiv ℝ M) 0 := hM.fderiv_right (m := 1) (by norm_num)
  have hD2 : HasFDerivAt (fderiv ℝ M) (fderiv ℝ (fderiv ℝ M) 0) (γ 0) := by
    rw [hγ0]; exact (hf1.differentiableAt (by norm_num)).hasFDerivAt
  have hc : HasDerivAt (fun t => fderiv ℝ M (γ t)) (fderiv ℝ (fderiv ℝ M) 0 (deriv γ 0)) 0 :=
    hD2.comp_hasDerivAt 0 (hγd.self_of_nhds.hasDerivAt)
  have hψ := tr3CLM.hasFDerivAt.comp_hasDerivAt 0 (hc.clm_apply hγ2.hasDerivAt)
  have hval : tr3CLM (fderiv ℝ (fderiv ℝ M) 0 (deriv γ 0) (deriv γ 0) +
      fderiv ℝ M (γ 0) (deriv (deriv γ) 0)) = tr3 (fderiv ℝ (fderiv ℝ M) 0 (deriv γ 0) (deriv γ 0)) := by
    rw [map_add, tr3CLM_apply, tr3CLM_apply, hγ0, htr, add_zero]
  rw [hval] at hψ
  have hψ' : HasDerivAt (fun t => tr3 (fderiv ℝ M (γ t) (deriv γ t)))
      (tr3 (fderiv ℝ (fderiv ℝ M) 0 (deriv γ 0) (deriv γ 0))) 0 := by
    simpa [tr3CLM_apply, Function.comp_def] using hψ
  refine ⟨?_, ?_, ?_⟩
  · filter_upwards [hφ] with t ht; exact ht.differentiableAt
  · exact hψ'.differentiableAt.congr_of_eventuallyEq hderiv_eq
  · rw [hderiv_eq.deriv_eq]; exact hψ'.deriv

section steps

variable (I : BolzaInputs)

/-- **Lemma 4 (tr A ≡ 0), PROVED IN LEAN** from r-equivariance (Lemma 3b) and the exact fact that
r has no invariant vector on T_B (`no_invariant_functional`; ⟨χ_T, 1⟩ = 0). -/
theorem trace_A_zero : ∀ v, tr3 (fderiv ℝ I.M 0 v) = 0 := by
  set A := fderiv ℝ I.M 0
  have hinv : ∀ v, tr3 (A (act Xr Xr v)) = tr3 (A v) := fun v =>
    tr3_conj_orth _ _ I.ρr_orth _ (I.equiv_r v)
  set φ : TB →ₗ[ℝ] ℝ := tr3L.comp (A : TB →ₗ[ℝ] Sq3)
  have φ1 : φ.comp (LinearMap.inl ℝ _ _) = 0 := by
    apply no_invariant_functional
    intro y
    have := hinv (y, 0)
    rw [act_inl] at this
    simpa [φ, tr3L] using this
  have φ2 : φ.comp (LinearMap.inr ℝ _ _) = 0 := by
    apply no_invariant_functional
    intro z
    have := hinv (0, z)
    rw [act_inr] at this
    simpa [φ, tr3L] using this
  intro v
  have hv : v = (v.1, 0) + (0, v.2) := by simp
  have e1 := LinearMap.congr_fun φ1 v.1
  have e2 := LinearMap.congr_fun φ2 v.2
  simp only [LinearMap.comp_apply, LinearMap.inl_apply, LinearMap.inr_apply,
    LinearMap.zero_apply] at e1 e2
  have : φ v = 0 := by rw [hv, map_add, e1, e2, add_zero]
  simpa [φ, tr3L] using this

/-- **Lemma 8 "Consequence", PROVED IN LEAN:** from Ky Fan (Λ∘γ ≤ g, equality at 0), the
certified g″(0) ≤ −9.3411 (Lehmann–Goerisch; the SU13 value −12.8504 is an alternative, any
negative constant suffices) and dΛ(0) = 0 (Lemma 4), D²Λ(w_B, w_B) ≤ −9.3411. -/
theorem wB_sign : tr3 (fderiv ℝ (fderiv ℝ I.M) 0 I.wB I.wB) ≤ -(93411 / 10000 : ℝ) := by
  obtain ⟨hd, hd2, hval⟩ := second_deriv_trace_comp I.M I.M_C2 I.γ I.γ_zero I.γ_diff I.γ_diff2
    (trace_A_zero I)
  have := second_deriv_le_of_le I.kyfan_le I.kyfan_eq hd I.g_diff hd2 I.g_diff2
  rw [hval, I.γ_tangent] at this
  exact this.trans I.g_sign

/-- **Lemma 6(iv): K₀ ⊂ ker A, PROVED IN LEAN** (antipodal argument). -/
theorem A_vanishes_K0 : ∀ y, fderiv ℝ I.M 0 (y, 0) = 0 := by
  intro y
  rcases rhoE_a_scalar [I.ρr, I.ρt, I.ρk] I.E_irred I.ρa I.ρJ I.J_sq I.J_central I.a_sq
    I.a_central_modJ with h | h
  all_goals
    have e := I.equiv_a (y, 0)
    rw [h] at e
    simp only [transpose_one, transpose_neg, Matrix.one_mul, Matrix.mul_one, neg_mul, mul_neg,
      neg_neg] at e
    have hact : act (-1) 1 (y, 0) = -(y, 0) := by
      simp [act, neg_mulVec, one_mulVec]
    rw [hact, map_neg] at e
    have e' : -(fderiv ℝ I.M 0 (y, 0)) = fderiv ℝ I.M 0 (y, 0) := by
      funext i j; exact congrFun (congrFun e i) j
    funext i j
    have := congrFun (congrFun e' i) j
    simp only [Pi.neg_apply] at this
    simp only [Pi.zero_apply]
    linarith

/-- **Lemma 6(v): ker A ⊂ K₀, PROVED IN LEAN** from W₀ irreducible, A|K₀ = 0 and A ≠ 0. -/
theorem kerA_sub_K0 : ∀ k ∈ LinearMap.ker (fderiv ℝ I.M 0 : TB →ₗ[ℝ] Sq3), k.2 = 0 := by
  set A := fderiv ℝ I.M 0
  have hK0 := A_vanishes_K0 I
  set U : Submodule ℝ (Fin 3 → ℝ) := LinearMap.ker ((A : TB →ₗ[ℝ] Sq3).comp (LinearMap.inr ℝ _ _))
  have hmem : ∀ z, z ∈ U ↔ A (0, z) = 0 := by intro z; simp [U]
  have hU := Wplus_irreducible U
    (by
      intro u hu
      rw [hmem] at hu ⊢
      have e := I.equiv_t (0, u)
      rw [act_inr, hu] at e
      exact eq_zero_of_conj_zero 0 _ _ e rfl)
    (by
      intro u hu
      rw [hmem] at hu ⊢
      have e := I.equiv_k (0, u)
      rw [act_inr, hu] at e
      exact eq_zero_of_conj_zero 0 _ _ e rfl)
  rcases hU with hbot | htop
  · intro k hk
    have hk' : A k = 0 := hk
    have hsplit : k = (k.1, 0) + (0, k.2) := by simp
    rw [hsplit, map_add, hK0, zero_add] at hk'
    have : k.2 ∈ U := (hmem _).2 hk'
    rw [hbot, Submodule.mem_bot] at this
    exact this
  · exfalso
    apply I.A_ne_zero
    apply ContinuousLinearMap.ext
    intro v
    have hz : A (0, v.2) = 0 := (hmem _).1 (by rw [htop]; exact Submodule.mem_top)
    have hsplit : v = (v.1, 0) + (0, v.2) := by simp
    rw [hsplit, map_add, hK0, hz, add_zero]
    rfl

/-- the restriction of D²Λ = tr D²M to K₀, as a bilinear form on ℝ³ -/
noncomputable def hessK0 : (Fin 3 → ℝ) →ₗ[ℝ] (Fin 3 → ℝ) →ₗ[ℝ] ℝ :=
  LinearMap.mk₂ ℝ (fun y y' => tr3 (fderiv ℝ (fderiv ℝ I.M) 0 (y, 0) (y', 0)))
    (by
      intro y₁ y₂ y'
      rw [show ((y₁ + y₂, 0) : TB) = (y₁, 0) + (y₂, 0) by simp, map_add,
        ContinuousLinearMap.add_apply, tr3_add])
    (by
      intro c y y'
      rw [show ((c • y, 0) : TB) = c • (y, 0) by simp, map_smul,
        ContinuousLinearMap.smul_apply, tr3_smul, smul_eq_mul])
    (by
      intro y y₁ y₂
      rw [show ((y₁ + y₂, 0) : TB) = (y₁, 0) + (y₂, 0) by simp, map_add, tr3_add])
    (by
      intro c y y'
      rw [show ((c • y', 0) : TB) = c • (y', 0) by simp, map_smul, tr3_smul, smul_eq_mul])

lemma hessK0_inv (X Y : Mat3)
    (h : ∀ v w, tr3 (fderiv ℝ (fderiv ℝ I.M) 0 (act X Y v) (act X Y w)) =
      tr3 (fderiv ℝ (fderiv ℝ I.M) 0 v w)) :
    Xᵀ * LinearMap.toMatrix₂' ℝ (hessK0 I) * X = LinearMap.toMatrix₂' ℝ (hessK0 I) := by
  have : (hessK0 I).compl₁₂ (Matrix.toLin' X) (Matrix.toLin' X) = hessK0 I := by
    apply LinearMap.ext₂
    intro y y'
    have := h (y, 0) (y', 0)
    rw [act_inl, act_inl] at this
    simpa [hessK0, Matrix.toLin'_apply] using this
  have e := LinearMap.toMatrix₂'_compl₁₂ (hessK0 I) (Matrix.toLin' X) (Matrix.toLin' X)
  rw [this, LinearMap.toMatrix'_toLin'] at e
  exact e.symm

/-- **Lemma 6(vi) + Lemma 11 + Lemma 8 ⇒ D²Λ|_{K₀} = c·g_inv with c < 0, PROVED IN LEAN.** -/
theorem hess_K0_neg : ∃ c < 0, ∀ y, tr3 (fderiv ℝ (fderiv ℝ I.M) 0 (y, 0) (y, 0)) = c * (y ⬝ᵥ y) := by
  set S := LinearMap.toMatrix₂' ℝ (hessK0 I)
  have hS : S = S 0 0 • (1 : Mat3) :=
    invariant_form_Kminus S (hessK0_inv I _ _ I.hess_inv_r) (hessK0_inv I _ _ I.hess_inv_t)
      (hessK0_inv I _ _ I.hess_inv_k)
  have hform : ∀ y, tr3 (fderiv ℝ (fderiv ℝ I.M) 0 (y, 0) (y, 0)) = S 0 0 * (y ⬝ᵥ y) := by
    intro y
    have e1 : hessK0 I y y = tr3 (fderiv ℝ (fderiv ℝ I.M) 0 (y, 0) (y, 0)) := by simp [hessK0]
    have e2 : hessK0 I = Matrix.toLinearMap₂' ℝ S := by simp [S]
    rw [← e1, e2, Matrix.toLinearMap₂'_apply', hS]
    simp [smul_mulVec, dotProduct_smul]
  refine ⟨S 0 0, ?_, hform⟩
  -- the certified direction w_B lies on the line of (e₁, 0)
  obtain ⟨hy, hz⟩ := line_r_smid I.wB.1 I.wB.2
    (by have := congrArg Prod.fst I.wB_r; simpa [act] using this)
    (by have := congrArg Prod.snd I.wB_r; simpa [act] using this)
    (by have := congrArg Prod.fst I.wB_smid; simpa [act] using this)
    (by have := congrArg Prod.snd I.wB_smid; simpa [act] using this)
  have hw : I.wB = (I.wB.1, 0) := by ext <;> simp [hz]
  have hsign := wB_sign I
  rw [hw, hform] at hsign
  have hpos : 0 < I.wB.1 ⬝ᵥ I.wB.1 := by
    apply dot_self_pos
    intro h0
    apply I.wB_ne
    rw [hw, h0]; rfl
  by_contra hc
  push Not at hc
  have := mul_nonneg hc hpos.le
  linarith

end steps

/-- **Main theorem (PROOF.md §0), conditional on `BolzaInputs`.  PROVED IN LEAN.**
The Bolza surface is a strict local maximum of λ₁ on T₂. -/
theorem bolza_strict_local_max (I : BolzaInputs) :
    ∀ᶠ X in @nhds I.T2 I.topT2 I.B, X ≠ I.B → I.lam1 X < I.lam1 I.B := by
  let _ := I.topT2
  obtain ⟨c, hc, hform⟩ := hess_K0_neg I
  have hneg : ∀ k ∈ LinearMap.ker (fderiv ℝ I.M 0 : TB →ₗ[ℝ] Sq3),
      tr3 (fderiv ℝ (fderiv ℝ I.M) 0 k k) ≤ -(-c) * ‖k‖ ^ 2 := by
    intro k hk
    have h2 := kerA_sub_K0 I k hk
    have hk : k = (k.1, 0) := by ext <;> simp [h2]
    rw [hk, hform]
    have hn : ‖((k.1, 0) : TB)‖ = ‖k.1‖ := by simp [Prod.norm_def]
    rw [hn, neg_neg]
    have := norm_sq_le_dot k.1
    nlinarith
  have hcrit := local_max_criterion_C2 I.M I.lamStar I.M_C2 I.M_zero I.M_sym (trace_A_zero I)
    (-c) (by linarith) hneg
  have hev : ∀ᶠ v in 𝓝 (0 : TB), v ≠ 0 → I.lam1 (I.chart v) < I.lamStar := by
    filter_upwards [hcrit, I.lam1_le_rayleigh] with v h1 h2 hv0
    obtain ⟨x, hx0, hx⟩ := h1 hv0
    have hpos := dot_self_pos hx0
    have := h2 x
    by_contra hge
    push Not at hge
    have := mul_le_mul_of_nonneg_right hge hpos.le
    linarith
  have himg := I.chart_nhds _ hev
  filter_upwards [himg] with X hX hXB
  obtain ⟨v, hv, rfl⟩ := hX
  rw [I.lam1_B]
  apply hv
  rintro rfl
  exact hXB I.chart_zero

end Bolza
