/-
# Non-vacuity check: a toy model of `BolzaInputs` (PROVED IN LEAN)

A structure with ~40 hypothesis fields could be accidentally contradictory, in which case
`bolza_strict_local_max` would be vacuous.  Here we build an explicit (toy, NOT geometric) instance
of `BolzaInputs`, so the hypotheses are jointly satisfiable:

* T₂ := TB, chart := id, B := 0;
* E := the K₀-type representation (ρ_r, ρ_t, ρ_κ) = (Xr, XmT, XmK), ρ_a = −1, ρ_J = 1;
* A(y, z) = [[0, −z₀, z₁], [−z₀, 0, z₂], [z₁, z₂, 0]] — the (unique up to scale) equivariant map
  W₀ → Sym(E), found numerically and verified here;
* M(v) = λ·1 + A(v) − 3‖v‖₂²·1, so D²(tr M)(v, w) = −18 ⟨v, w⟩;
* w_B = (e₁, 0), γ(t) = t·w_B, g = tr M ∘ γ (so g″(0) = −18 ≤ −9.3411);
* λ₁(v) := λ − 9‖M(v) − λ·1‖ (a lower bound for all Rayleigh quotients).

This says nothing about the Bolza surface; it only shows that the interface is consistent.
-/
import BolzaKit.LocalMax.Assembly

open Matrix Filter Topology

namespace Bolza.Model

/-- the intertwiner W₀ → Sym(E) -/
def Aw (z : Fin 3 → ℝ) : Sq3 := ![![0, -z 0, z 1], ![-z 0, 0, z 2], ![z 1, z 2, 0]]

lemma Aw_add (z z' : Fin 3 → ℝ) : Aw (z + z') = Aw z + Aw z' := by
  funext i j; fin_cases i <;> fin_cases j <;> simp [Aw] <;> ring

lemma Aw_smul (c : ℝ) (z : Fin 3 → ℝ) : Aw (c • z) = c • Aw z := by
  funext i j; fin_cases i <;> fin_cases j <;> simp [Aw]

/-- A as a linear map on TB (depends only on the W₀ component) -/
def AlinL : TB →ₗ[ℝ] Sq3 where
  toFun v := Aw v.2
  map_add' v w := by simp [Aw_add]
  map_smul' c v := by simp [Aw_smul]

noncomputable def Alin : TB →L[ℝ] Sq3 := LinearMap.toContinuousLinearMap AlinL

/-- dot product on TB -/
def dotTB (v w : TB) : ℝ := v.1 ⬝ᵥ w.1 + v.2 ⬝ᵥ w.2

def BqL : TB →ₗ[ℝ] TB →ₗ[ℝ] Sq3 :=
  LinearMap.mk₂ ℝ (fun v w => (-3 * dotTB v w) • one3)
    (by intro v v' w; simp [dotTB, add_dotProduct]; module)
    (by intro c v w; simp [dotTB, smul_dotProduct, smul_smul]; module)
    (by intro v w w'; simp [dotTB, dotProduct_add]; module)
    (by intro c v w; simp [dotTB, dotProduct_smul, smul_smul]; module)

noncomputable def Bq : TB →L[ℝ] TB →L[ℝ] Sq3 :=
  LinearMap.toContinuousLinearMap
    ((LinearMap.toContinuousLinearMap : (TB →ₗ[ℝ] Sq3) ≃ₗ[ℝ] (TB →L[ℝ] Sq3)).toLinearMap ∘ₗ BqL)

lemma Bq_apply (v w : TB) : Bq v w = (-3 * dotTB v w) • one3 := rfl

lemma Alin_apply (v : TB) : Alin v = Aw v.2 := rfl

def lamS : ℝ := 4

noncomputable def Mm (v : TB) : Sq3 := lamS • one3 + Alin v + Bq v v

lemma Mm_hasFDerivAt (y : TB) : HasFDerivAt Mm (Alin + (Bq + Bq.flip) y) y := by
  have h := ((hasFDerivAt_const (lamS • one3) y).add Alin.hasFDerivAt).add
    (Bq.hasFDerivAt.clm_apply (hasFDerivAt_id y))
  have h' : HasFDerivAt Mm (0 + Alin + ((Bq y).comp (ContinuousLinearMap.id ℝ TB) + Bq.flip y)) y :=
    h
  refine h'.congr_fderiv ?_
  apply ContinuousLinearMap.ext
  intro w
  simp

lemma fderiv_Mm : fderiv ℝ Mm = fun y => Alin + (Bq + Bq.flip) y := by
  funext y; exact (Mm_hasFDerivAt y).fderiv

lemma fderiv_Mm_zero : fderiv ℝ Mm 0 = Alin := by
  rw [fderiv_Mm]; simp

lemma fderiv2_Mm : fderiv ℝ (fderiv ℝ Mm) 0 = Bq + Bq.flip := by
  rw [fderiv_Mm]
  exact ((hasFDerivAt_const Alin 0).add (Bq + Bq.flip).hasFDerivAt).fderiv.trans (by simp)

lemma D2_apply (v w : TB) : fderiv ℝ (fderiv ℝ Mm) 0 v w = (-6 * dotTB v w) • one3 := by
  rw [fderiv2_Mm]
  simp only [ContinuousLinearMap.add_apply, ContinuousLinearMap.flip_apply, Bq_apply]
  have : dotTB w v = dotTB v w := by simp [dotTB, dotProduct_comm]
  rw [this, ← add_smul]; ring_nf

lemma tr3_one3 : tr3 one3 = 3 := by simp [tr3, one3]; norm_num

lemma Mm_contDiff : ContDiff ℝ 2 Mm := by
  have hB : ContDiff ℝ 2 (fun v : TB => Bq v v) :=
    Bq.contDiff.clm_apply contDiff_id
  have h1 : ContDiff ℝ 2 (fun _ : TB => lamS • one3) := contDiff_const
  have h2 : ContDiff ℝ 2 (fun v : TB => Alin v) := Alin.contDiff
  unfold Mm
  exact (h1.add h2).add hB

lemma orth_dot (X : Mat3) (hX : Xᵀ * X = 1) (a b : Fin 3 → ℝ) : (X *ᵥ a) ⬝ᵥ (X *ᵥ b) = a ⬝ᵥ b := by
  calc (X *ᵥ a) ⬝ᵥ (X *ᵥ b) = ((a ᵥ* Xᵀ) ᵥ* X) ⬝ᵥ b := by
        rw [dotProduct_mulVec, vecMul_transpose]
    _ = a ⬝ᵥ b := by rw [vecMul_vecMul, hX, vecMul_one]

lemma orth_Xr : Xrᵀ * Xr = 1 := by
  ext i j; fin_cases i <;> fin_cases j <;> simp [Xr, Matrix.mul_apply, Fin.sum_univ_three]
lemma orth_XmT : XmTᵀ * XmT = 1 := by
  ext i j; fin_cases i <;> fin_cases j <;> simp [XmT, Matrix.mul_apply, Fin.sum_univ_three]
lemma orth_XmK : XmKᵀ * XmK = 1 := by
  ext i j; fin_cases i <;> fin_cases j <;> simp [XmK, Matrix.mul_apply, Fin.sum_univ_three]
lemma orth_XpT : XpTᵀ * XpT = 1 := by
  ext i j; fin_cases i <;> fin_cases j <;> simp [XpT, Matrix.mul_apply, Fin.sum_univ_three]
lemma orth_XpK : XpKᵀ * XpK = 1 := by
  ext i j; fin_cases i <;> fin_cases j <;> simp [XpK, Matrix.mul_apply, Fin.sum_univ_three]

lemma dot_act (X Y : Mat3) (hX : Xᵀ * X = 1) (hY : Yᵀ * Y = 1) (v w : TB) :
    dotTB (act X Y v) (act X Y w) = dotTB v w := by
  simp [dotTB, act, orth_dot X hX, orth_dot Y hY]

lemma hess_inv_model (X Y : Mat3) (hX : Xᵀ * X = 1) (hY : Yᵀ * Y = 1) (v w : TB) :
    tr3 (fderiv ℝ (fderiv ℝ Mm) 0 (act X Y v) (act X Y w)) = tr3 (fderiv ℝ (fderiv ℝ Mm) 0 v w) := by
  rw [D2_apply, D2_apply, dot_act X Y hX hY]

lemma equiv_model (Xm Xp : Mat3) (h : ∀ z, Matrix.of (Aw (Xp *ᵥ z)) = Xm * Matrix.of (Aw z) * Xmᵀ)
    (v : TB) : Matrix.of (fderiv ℝ Mm 0 (act Xm Xp v)) = Xm * Matrix.of (fderiv ℝ Mm 0 v) * Xmᵀ := by
  rw [fderiv_Mm_zero, Alin_apply, Alin_apply]
  exact h v.2

lemma intertw_r (z : Fin 3 → ℝ) : Matrix.of (Aw (Xr *ᵥ z)) = Xr * Matrix.of (Aw z) * Xrᵀ := by
  ext i j; fin_cases i <;> fin_cases j <;>
    simp [Aw, Xr, Matrix.mul_apply, Fin.sum_univ_three, mulVec, dotProduct]
lemma intertw_t (z : Fin 3 → ℝ) : Matrix.of (Aw (XpT *ᵥ z)) = XmT * Matrix.of (Aw z) * XmTᵀ := by
  ext i j; fin_cases i <;> fin_cases j <;>
    simp [Aw, XpT, XmT, Matrix.mul_apply, Fin.sum_univ_three, mulVec, dotProduct]
lemma intertw_k (z : Fin 3 → ℝ) : Matrix.of (Aw (XpK *ᵥ z)) = XmK * Matrix.of (Aw z) * XmKᵀ := by
  ext i j; fin_cases i <;> fin_cases j <;>
    simp [Aw, XpK, XmK, Matrix.mul_apply, Fin.sum_univ_three, mulVec, dotProduct]

lemma XmT_mulVec (u : Fin 3 → ℝ) : XmT *ᵥ u = ![u 2, u 0, u 1] := by
  ext i; fin_cases i <;> simp [XmT, mulVec, dotProduct, Fin.sum_univ_three]
lemma XmK_mulVec (u : Fin 3 → ℝ) : XmK *ᵥ u = ![u 0, u 1, -u 2] := by
  ext i; fin_cases i <;> simp [XmK, mulVec, dotProduct, Fin.sum_univ_three]
lemma Xr_mulVec (u : Fin 3 → ℝ) : Xr *ᵥ u = ![-u 2, -u 1, u 0] := by
  ext i; fin_cases i <;> simp [Xr, mulVec, dotProduct, Fin.sum_univ_three]

/-- the K₀-type representation is irreducible (model's E) -/
theorem Kminus_irreducible : Irreducible3 [Xr, XmT, XmK] := by
  intro U hU
  have hT : ∀ u ∈ U, XmT *ᵥ u ∈ U := hU XmT (by simp)
  have hK : ∀ u ∈ U, XmK *ᵥ u ∈ U := hU XmK (by simp)
  have hR : ∀ u ∈ U, Xr *ᵥ u ∈ U := hU Xr (by simp)
  by_cases h : ∃ u ∈ U, u ≠ 0
  · right
    obtain ⟨u, hu, hu0⟩ := h
    have proj : ∀ w ∈ U, (w 2) • (![0, 0, 1] : Fin 3 → ℝ) ∈ U := by
      intro w hw
      have : (w 2) • (![0, 0, 1] : Fin 3 → ℝ) = (1 / 2 : ℝ) • (w - XmK *ᵥ w) := by
        rw [XmK_mulVec]; ext i; fin_cases i <;> simp [vecHead, vecTail] <;> ring
      rw [this]; exact U.smul_mem _ (U.sub_mem hw (hK w hw))
    have h1 := proj u hu
    have h2 := proj _ (hT u hu)
    have h3 := proj _ (hT _ (hT u hu))
    rw [XmT_mulVec] at h2
    rw [XmT_mulVec, XmT_mulVec] at h3
    simp only [Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two,
      Matrix.head_cons, Matrix.tail_cons] at h2 h3
    have unit : ∀ c : ℝ, c ≠ 0 → c • (![0, 0, 1] : Fin 3 → ℝ) ∈ U → (![0, 0, 1] : Fin 3 → ℝ) ∈ U := by
      intro c hc hcU
      have := U.smul_mem c⁻¹ hcU
      rwa [smul_smul, inv_mul_cancel₀ hc, one_smul] at this
    have e2 : (![0, 0, 1] : Fin 3 → ℝ) ∈ U := by
      by_cases a2 : u 2 = 0
      · by_cases a1 : u 1 = 0
        · have a0 : u 0 ≠ 0 := by
            intro a0; apply hu0; ext i; fin_cases i <;> simp [a0, a1, a2]
          exact unit _ a0 h3
        · exact unit _ a1 h2
      · exact unit _ a2 h1
    have e0 : (![1, 0, 0] : Fin 3 → ℝ) ∈ U := by
      have := U.neg_mem (hR _ e2)
      convert this using 1; rw [Xr_mulVec]; ext i; fin_cases i <;> simp
    have e1 : (![0, 1, 0] : Fin 3 → ℝ) ∈ U := by
      have := hT _ e0
      convert this using 1; rw [XmT_mulVec]; ext i; fin_cases i <;> simp
    rw [eq_top_iff]
    intro v _
    have : v = v 0 • (![1, 0, 0] : Fin 3 → ℝ) + v 1 • ![0, 1, 0] + v 2 • ![0, 0, 1] := by
      ext i; fin_cases i <;> simp
    rw [this]
    exact U.add_mem (U.add_mem (U.smul_mem _ e0) (U.smul_mem _ e1)) (U.smul_mem _ e2)
  · left
    push Not at h
    exact (Submodule.eq_bot_iff U).2 h

lemma abs_qf_le_dot (S : Sq3) (x : Fin 3 → ℝ) : |qf S x| ≤ 9 * ‖S‖ * (x ⬝ᵥ x) := by
  have hdot : x ⬝ᵥ x = x 0 * x 0 + x 1 * x 1 + x 2 * x 2 := by
    simp [dotProduct, Fin.sum_univ_three]
  have hterm : ∀ i j, |x i * S i j * x j| ≤ ‖S‖ * (x ⬝ᵥ x) := by
    intro i j
    have hS := entry_le_norm S i j
    have hxx : |x i| * |x j| ≤ x ⬝ᵥ x := by
      rw [hdot]
      have h1 : |x i| * |x j| ≤ (x i * x i + x j * x j) / 2 := by
        nlinarith [sq_nonneg (|x i| - |x j|), sq_abs (x i), sq_abs (x j), abs_mul_abs_self (x i),
          abs_mul_abs_self (x j)]
      fin_cases i <;> fin_cases j <;> simp at h1 ⊢ <;>
        nlinarith [mul_self_nonneg (x 0), mul_self_nonneg (x 1), mul_self_nonneg (x 2)]
    rw [abs_mul, abs_mul]
    calc |x i| * |S i j| * |x j| = |S i j| * (|x i| * |x j|) := by ring
      _ ≤ ‖S‖ * (x ⬝ᵥ x) := mul_le_mul hS hxx (by positivity) (norm_nonneg _)
  calc |qf S x| ≤ ∑ i, |∑ j, x i * S i j * x j| := Finset.abs_sum_le_sum_abs _ _
    _ ≤ ∑ i, ∑ j, |x i * S i j * x j| :=
        Finset.sum_le_sum fun i _ => Finset.abs_sum_le_sum_abs _ _
    _ ≤ ∑ _i : Fin 3, ∑ _j : Fin 3, ‖S‖ * (x ⬝ᵥ x) :=
        Finset.sum_le_sum fun i _ => Finset.sum_le_sum fun j _ => hterm i j
    _ = 9 * ‖S‖ * (x ⬝ᵥ x) := by simp; ring

lemma Mm_apply (v : TB) (i j : Fin 3) :
    Mm v i j = lamS * one3 i j + Aw v.2 i j + (-3 * dotTB v v) * one3 i j := by
  simp only [Mm, Pi.add_apply, Pi.smul_apply, smul_eq_mul, Alin_apply, Bq_apply]

def wBm : TB := (![0, 1, 0], 0)

noncomputable def γm (t : ℝ) : TB := t • wBm

noncomputable def gm (t : ℝ) : ℝ := tr3 (Mm (γm t))

lemma gm_eq : gm = fun t => 12 - 9 * t ^ 2 := by
  funext t
  simp only [gm, tr3, Mm_apply]
  simp [γm, Aw, wBm, dotTB, dotProduct, Fin.sum_univ_three, one3, lamS]
  ring

lemma deriv_gm : deriv gm = fun t => -18 * t := by
  rw [gm_eq]; funext t
  have := ((hasDerivAt_pow 2 t).const_mul 9).const_sub 12
  rw [this.deriv]; push_cast; ring

lemma hasDerivAt_γm (t : ℝ) : HasDerivAt γm wBm t := by
  have h := (hasDerivAt_id t).smul_const wBm
  simp only [id, one_smul] at h
  exact h

lemma deriv_γm : deriv γm = fun _ => wBm := by
  funext t; exact (hasDerivAt_γm t).deriv

lemma Mm_sym (v : TB) : IsSym3 (Mm v) := by
  intro i j
  rw [Mm_apply, Mm_apply]
  fin_cases i <;> fin_cases j <;> simp [Aw, one3]

/-- **The toy model: `BolzaInputs` is satisfiable.** -/
noncomputable def toyInputs : BolzaInputs where
  T2 := TB
  topT2 := inferInstance
  B := 0
  lam1 v := lamS - 9 * ‖Mm v - lamS • one3‖
  chart := id
  chart_zero := rfl
  chart_nhds s hs := by simpa using hs
  M := Mm
  lamStar := lamS
  lam1_B := by simp [Mm, Bq_apply, Alin_apply, Aw]
  M_C2 := Mm_contDiff.contDiffAt
  M_sym := Mm_sym
  M_zero := by
    have : Aw 0 = 0 := by funext i j; fin_cases i <;> fin_cases j <;> simp [Aw]
    simp [Mm, Bq_apply, Alin_apply, this, dotTB]
  lam1_le_rayleigh := Filter.Eventually.of_forall fun v x => by
    have h := abs_qf_le_dot (Mm v - lamS • one3) x
    have e : qf (Mm v) x = lamS * (x ⬝ᵥ x) + qf (Mm v - lamS • one3) x := by
      rw [qf_sub, qf_smul_one]; ring
    have hx : 0 ≤ x ⬝ᵥ x := by
      simp only [dotProduct, Fin.sum_univ_three]
      nlinarith [mul_self_nonneg (x 0), mul_self_nonneg (x 1), mul_self_nonneg (x 2)]
    show (lamS - 9 * ‖Mm v - lamS • one3‖) * (x ⬝ᵥ x) ≤ qf (Mm v) x
    rw [e]
    have := (abs_le.1 h).1
    nlinarith
  ρr := Xr
  ρt := XmT
  ρk := XmK
  ρa := -1
  ρJ := 1
  ρr_orth := orth_Xr
  equiv_r := equiv_model Xr Xr intertw_r
  equiv_t := equiv_model XmT XpT intertw_t
  equiv_k := equiv_model XmK XpK intertw_k
  equiv_a v := by
    rw [fderiv_Mm_zero, Alin_apply, Alin_apply]
    simp [act]
  E_irred := Kminus_irreducible
  J_sq := by simp
  J_central := by intro g _; simp
  a_sq := by simp
  a_central_modJ := by intro g _; left; simp
  A_ne_zero := by
    rw [fderiv_Mm_zero]
    intro h
    have := congrFun (congrFun (congrArg (fun L => L ((0 : Fin 3 → ℝ), ![1, 0, 0])) h) 0) 1
    simp [Alin_apply, Aw] at this
  hess_inv_r := hess_inv_model Xr Xr orth_Xr orth_Xr
  hess_inv_t := hess_inv_model XmT XpT orth_XmT orth_XpT
  hess_inv_k := hess_inv_model XmK XpK orth_XmK orth_XpK
  wB := wBm
  wB_ne := by
    intro h
    have := congrFun (congrArg Prod.fst h) 1
    simp [wBm] at this
  wB_r := by
    ext i <;> fin_cases i <;> simp [act, wBm, Xr, mulVec, dotProduct, Fin.sum_univ_three]
  wB_smid := by
    ext i <;> fin_cases i <;> simp [act, wBm, XmK, XpK, mulVec, dotProduct, Fin.sum_univ_three]
  γ := γm
  γ_zero := by simp [γm]
  γ_diff := Filter.Eventually.of_forall fun t => (hasDerivAt_γm t).differentiableAt
  γ_diff2 := by rw [deriv_γm]; exact differentiableAt_const _
  γ_tangent := by rw [deriv_γm]
  g := gm
  g_diff := Filter.Eventually.of_forall fun t => by
    rw [gm_eq]; fun_prop
  g_diff2 := by rw [deriv_gm]; fun_prop
  kyfan_le := Filter.Eventually.of_forall fun t => le_refl _
  kyfan_eq := rfl
  g_sign := by
    rw [deriv_gm]
    have := ((hasDerivAt_id (0 : ℝ)).const_mul (-18))
    rw [show (fun t : ℝ => -18 * t) = fun t => -18 * id t from rfl, this.deriv]
    norm_num

/-- the conclusion of the main theorem in the toy model (sanity check of the interface) -/
example : ∀ᶠ X in @nhds toyInputs.T2 toyInputs.topT2 toyInputs.B,
    X ≠ toyInputs.B → toyInputs.lam1 X < toyInputs.lam1 toyInputs.B :=
  bolza_strict_local_max toyInputs

end Bolza.Model
