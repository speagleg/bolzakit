/-
# Second-order Taylor expansion with uniform Peano remainder, and the C² form of Lemma 10

PROVED IN LEAN.  Mathlib has Taylor's theorem in one variable; the uniform multivariate
second-order Peano remainder needed by PROOF.md Lemma 10 is proved here from
* the mean value inequality on a ball,
* differentiability of `fderiv f` at 0 (from C²),
* symmetry of the second derivative (`ContDiffAt.isSymmSndFDerivAt`, Schwarz).
-/
import BolzaKit.LocalMax.LocalMax

open Filter Topology

namespace Bolza

variable {E F : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
  [NormedAddCommGroup F] [NormedSpace ℝ F]

/-- **Uniform second-order Taylor expansion at 0 for a C² map.**
`f v = f 0 + f'(0) v + ½ f''(0)(v,v) + o(‖v‖²)`. -/
theorem taylor2_littleO (f : E → F) (hf : ContDiffAt ℝ 2 f 0) :
    ∀ ε > 0, ∀ᶠ v in 𝓝 (0 : E),
      ‖f v - (f 0 + fderiv ℝ f 0 v + (1 / 2 : ℝ) • fderiv ℝ (fderiv ℝ f) 0 v v)‖ ≤ ε * ‖v‖ ^ 2 := by
  intro ε hε
  set D := fderiv ℝ f 0 with hD
  set D2 := fderiv ℝ (fderiv ℝ f) 0 with hD2
  -- f is differentiable near 0
  have hev : ∀ᶠ y in 𝓝 (0 : E), HasFDerivAt f (fderiv ℝ f y) y := by
    filter_upwards [hf.eventually (by simp)] with y hy
    exact (hy.differentiableAt (by norm_num)).hasFDerivAt
  -- fderiv f is differentiable at 0 with derivative D2
  have hf1 : ContDiffAt ℝ 1 (fderiv ℝ f) 0 := hf.fderiv_right (m := 1) (by norm_num)
  have hD2d : HasFDerivAt (fderiv ℝ f) D2 0 := (hf1.differentiableAt (by norm_num)).hasFDerivAt
  -- symmetry of the second derivative
  have hsymm : ∀ v w, D2 v w = D2 w v := by
    have := hf.isSymmSndFDerivAt (by simp [minSmoothness])
    intro v w; exact this v w
  -- little-o estimate for fderiv f
  have hlo := (hasFDerivAt_iff_isLittleO_nhds_zero.1 hD2d).def hε
  simp only [zero_add] at hlo
  obtain ⟨δ, hδ, hball⟩ := Metric.eventually_nhds_iff_ball.1 (hev.and hlo)
  filter_upwards [Metric.ball_mem_nhds (0 : E) hδ] with v hv
  set r := ‖v‖
  have hsub : Metric.closedBall (0 : E) r ⊆ Metric.ball 0 δ := by
    intro y hy
    simp only [Metric.mem_closedBall, dist_zero_right] at hy
    simp only [Metric.mem_ball, dist_zero_right] at hv ⊢
    linarith
  set g : E → F := fun y => f y - f 0 - D y - (1 / 2 : ℝ) • D2 y y with hg
  set g' : E → E →L[ℝ] F := fun y => fderiv ℝ f y - D - D2 y with hg'
  have hderiv : ∀ y ∈ Metric.closedBall (0 : E) r, HasFDerivWithinAt g (g' y) (Metric.closedBall 0 r) y := by
    intro y hy
    have hyb := hball y (hsub hy)
    have h1 : HasFDerivAt f (fderiv ℝ f y) y := hyb.1
    have h2 : HasFDerivAt (fun y => D2 y y) ((D2 y).comp (ContinuousLinearMap.id ℝ E) + D2.flip y) y :=
      (D2.hasFDerivAt).clm_apply (hasFDerivAt_id y)
    have h3 := ((h1.sub_const (f 0)).sub D.hasFDerivAt).sub (h2.const_smul (1 / 2 : ℝ))
    refine (h3.congr_fderiv ?_).hasFDerivWithinAt
    ext w
    simp only [hg', ContinuousLinearMap.sub_apply, ContinuousLinearMap.smul_apply,
      ContinuousLinearMap.add_apply, ContinuousLinearMap.comp_apply, ContinuousLinearMap.id_apply,
      ContinuousLinearMap.flip_apply]
    rw [hsymm w y]
    module
  have hbound : ∀ y ∈ Metric.closedBall (0 : E) r, ‖g' y‖ ≤ ε * r := by
    intro y hy
    have hyb := (hball y (hsub hy)).2
    simp only [Metric.mem_closedBall, dist_zero_right] at hy
    have : ‖g' y‖ ≤ ε * ‖y‖ := by
      simpa [hg'] using hyb
    exact this.trans (by gcongr)
  have hmv := (convex_closedBall (0 : E) r).norm_image_sub_le_of_norm_hasFDerivWithin_le
    hderiv hbound (Metric.mem_closedBall_self (norm_nonneg v))
    (show v ∈ Metric.closedBall (0 : E) r by simp [r])
  have hg0 : g 0 = 0 := by simp [hg]
  rw [hg0, sub_zero, sub_zero] at hmv
  have : f v - (f 0 + D v + (1 / 2 : ℝ) • D2 v v) = g v := by simp only [hg]; abel
  rw [this]
  calc ‖g v‖ ≤ ε * r * ‖v‖ := hmv
    _ = ε * ‖v‖ ^ 2 := by ring

/-- `fderiv` of a symmetric-matrix-valued map at a point is symmetric. -/
lemma fderiv_isSym3 {V : Type*} [NormedAddCommGroup V] [NormedSpace ℝ V]
    (M : V → Sq3) (hM : DifferentiableAt ℝ M 0) (hsym : ∀ v, IsSym3 (M v)) (v : V) :
    IsSym3 (fderiv ℝ M 0 v) := by
  intro i j
  -- evaluation at (i,j) and (j,i) are continuous linear maps
  set eij : Sq3 →L[ℝ] ℝ := (ContinuousLinearMap.proj (R := ℝ) (φ := fun _ : Fin 3 => ℝ) j).comp
    (ContinuousLinearMap.proj (R := ℝ) (φ := fun _ : Fin 3 => Fin 3 → ℝ) i)
  set eji : Sq3 →L[ℝ] ℝ := (ContinuousLinearMap.proj (R := ℝ) (φ := fun _ : Fin 3 => ℝ) i).comp
    (ContinuousLinearMap.proj (R := ℝ) (φ := fun _ : Fin 3 => Fin 3 → ℝ) j)
  have h1 : fderiv ℝ (fun x => eij (M x)) 0 = eij.comp (fderiv ℝ M 0) :=
    (eij.hasFDerivAt.comp 0 hM.hasFDerivAt).fderiv
  have h2 : fderiv ℝ (fun x => eji (M x)) 0 = eji.comp (fderiv ℝ M 0) :=
    (eji.hasFDerivAt.comp 0 hM.hasFDerivAt).fderiv
  have heq : (fun x => eij (M x)) = (fun x => eji (M x)) := by
    funext x; simp only [eij, eji, ContinuousLinearMap.comp_apply, ContinuousLinearMap.proj_apply]
    exact hsym x i j
  have := congrArg (fun L => L v) (h1.symm.trans (heq ▸ h2))
  simpa [eij, eji] using this

/-- **PROOF.md Lemma 10, C² form.  PROVED IN LEAN.**
`M : V → Sym(3)` C² at 0, `M 0 = λ•1`, `A := dM(0)` traceless, and the Hessian of `tr M`
negative definite on `ker A` in the quantitative form `tr D²M(0)(k,k) ≤ -κ‖k‖²`.
Then `λ_min M(v) < λ` (Rayleigh form) for all small `v ≠ 0`.
(Negative definiteness and the quantitative form are equivalent on a finite-dimensional
subspace; the assembly produces the quantitative form directly.) -/
theorem local_max_criterion_C2 {V : Type*} [NormedAddCommGroup V] [NormedSpace ℝ V]
    [FiniteDimensional ℝ V]
    (M : V → Sq3) (lam : ℝ) (hM : ContDiffAt ℝ 2 M 0) (hM0 : M 0 = lam • one3)
    (hsym : ∀ v, IsSym3 (M v)) (htr : ∀ v, tr3 (fderiv ℝ M 0 v) = 0)
    (κ : ℝ) (hκ : 0 < κ)
    (hneg : ∀ k ∈ LinearMap.ker (fderiv ℝ M 0 : V →ₗ[ℝ] Sq3),
      tr3 (fderiv ℝ (fderiv ℝ M) 0 k k) ≤ -κ * ‖k‖ ^ 2) :
    ∀ᶠ v in 𝓝 (0 : V), v ≠ 0 → ∃ x : Fin 3 → ℝ, x ≠ 0 ∧ qf (M v) x < lam * (x ⬝ᵥ x) := by
  have hdiff : DifferentiableAt ℝ M 0 := hM.differentiableAt (by norm_num)
  refine local_max_criterion M lam (fderiv ℝ M 0) ((1 / 2 : ℝ) • fderiv ℝ (fderiv ℝ M) 0)
    (fderiv_isSym3 M hdiff hsym) htr ?_ (κ / 2) (by positivity) ?_
  · intro ε hε
    filter_upwards [taylor2_littleO M hM ε hε] with v hv
    rw [hM0] at hv
    simpa using hv
  · intro k hk
    have := hneg k hk
    have e : tr3 (((1 / 2 : ℝ) • fderiv ℝ (fderiv ℝ M) 0) k k) =
        (1 / 2) * tr3 (fderiv ℝ (fderiv ℝ M) 0 k k) := by
      simp [tr3]; ring
    rw [e]; nlinarith

end Bolza
