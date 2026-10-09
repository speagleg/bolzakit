/-
# Schur-type facts for the block representation, over ℝ (PROVED IN LEAN)

Block coordinates (see `RepTheory.lean`): T_B ≅ ℝ³ × ℝ³ = T₋ × T₊ with T₋ = K₀ (−1-eigenspace of
the antipodal isometry a) and T₊ = W₀.  The generators act by signed permutation blocks
(Lean-checked to be the `P`-conjugates of the exact cert_chars.py matrices in `RepTheory.lean`):

  r : (Xr, Xr),   t : (XmT, XpT),   κ = s_mid : (XmK, XpK),   a : (−1, +1).

Proved here, for real matrices:
* `invariant_form_Kminus` (PROOF.md Lemma 6(vi)): every bilinear form on K₀ invariant under
  r, t, κ is a scalar multiple of the standard inner product (g_inv = I in this basis).
* `Wplus_irreducible` (used in PROOF.md Lemma 6(v)): T₊ = W₀ has no invariant subspaces
  other than 0 and itself (only t and κ are needed).
* `no_invariant_functional` (PROOF.md Lemma 4, proof 2 / Schur): a linear functional on ℝ³
  invariant under Xr vanishes (applies to both blocks).
* `involution_schur` (PROOF.md Lemma 6(iv)): an involution commuting with an irreducible family
  of 3×3 real matrices is ±1.
* `rhoE_a_scalar` (PROOF.md Lemma 6(iv)): with J central, J² = 1, a² = J, a central mod J,
  and E irreducible: ρ_E(a) = ±1.
* `line_r_smid` (PROOF.md Lemma 7): {v : r v = −v, s_mid v = v} = span((e₁, 0)) ⊂ K₀.
-/
import BolzaKit.LocalMax.LocalMax
import BolzaKit.LocalMax.RepTheory

open Matrix

namespace Bolza

/-- r on either block -/
def Xr : Mat3 := !![0, 0, -1; 0, -1, 0; 1, 0, 0]
/-- t on T₋ -/
def XmT : Mat3 := !![0, 0, 1; 1, 0, 0; 0, 1, 0]
/-- κ = s_mid on T₋ -/
def XmK : Mat3 := !![1, 0, 0; 0, 1, 0; 0, 0, -1]
/-- t on T₊ -/
def XpT : Mat3 := !![0, -1, 0; 0, 0, 1; -1, 0, 0]
/-- κ = s_mid on T₊ -/
def XpK : Mat3 := !![1, 0, 0; 0, -1, 0; 0, 0, -1]

/-- real matrix of a row-major integer list -/
def ofList (l : List Int) : Mat3 := Matrix.of fun i j => (Rep.getM 3 l i j : ℝ)

/-- The real matrices above are exactly the integer blocks used in `RepTheory`. -/
theorem blocks_match :
    ofList Rep.gR.m = Xr ∧ ofList Rep.gR.p = Xr ∧ ofList Rep.gT.m = XmT ∧
    ofList Rep.gT.p = XpT ∧ ofList Rep.gK.m = XmK ∧ ofList Rep.gK.p = XpK := by
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_⟩ <;>
  · ext i j
    fin_cases i <;> fin_cases j <;>
      simp [ofList, Rep.getM, Rep.gR, Rep.gT, Rep.gK, Xr, XmT, XmK, XpT, XpK]

/-! ### Invariant bilinear forms on K₀ -/

/-- **Lemma 6(vi).** A bilinear form `S` on K₀ ≅ ℝ³ with `gᵀ S g = S` for g = r, t, κ is `c • 1`. -/
theorem invariant_form_Kminus (S : Mat3) (hr : Xrᵀ * S * Xr = S) (ht : XmTᵀ * S * XmT = S)
    (hk : XmKᵀ * S * XmK = S) : S = S 0 0 • (1 : Mat3) := by
  have er := fun i j => congrFun (congrFun hr i) j
  have et := fun i j => congrFun (congrFun ht i) j
  have ek := fun i j => congrFun (congrFun hk i) j
  simp only [Matrix.mul_apply, Fin.sum_univ_three, Matrix.transpose_apply, Xr, XmT, XmK,
    Matrix.of_apply, Matrix.cons_val', Matrix.cons_val_zero, Matrix.cons_val_one,
    Matrix.cons_val_two, Matrix.empty_val', Matrix.cons_val_fin_one, Matrix.head_cons,
    Matrix.head_fin_const] at er et ek
  have r01 := er 0 1; have r02 := er 0 2; have r12 := er 1 2; have r00 := er 0 0
  have r10 := er 1 0; have r20 := er 2 0; have r21 := er 2 1
  have t00 := et 0 0; have t11 := et 1 1; have t01 := et 0 1; have t10 := et 1 0
  have t02 := et 0 2; have t20 := et 2 0; have t12 := et 1 2; have t21 := et 2 1
  have k02 := ek 0 2; have k12 := ek 1 2; have k20 := ek 2 0; have k21 := ek 2 1
  norm_num at r01 r02 r12 r00 r10 r20 r21 t00 t11 t01 t10 t02 t20 t12 t21 k02 k12 k20 k21
  ext i j
  fin_cases i <;> fin_cases j <;> simp <;> linarith

/-! ### Irreducibility of W₀ = T₊ -/

lemma XpT_mulVec (u : Fin 3 → ℝ) : XpT *ᵥ u = ![-u 1, u 2, -u 0] := by
  ext i; fin_cases i <;> simp [XpT, mulVec, dotProduct, Fin.sum_univ_three]

lemma XpK_mulVec (u : Fin 3 → ℝ) : XpK *ᵥ u = ![u 0, -u 1, -u 2] := by
  ext i; fin_cases i <;> simp [XpK, mulVec, dotProduct, Fin.sum_univ_three]

/-- **W₀ is irreducible** (invariance under t and κ suffices). -/
theorem Wplus_irreducible (U : Submodule ℝ (Fin 3 → ℝ))
    (hT : ∀ u ∈ U, XpT *ᵥ u ∈ U) (hK : ∀ u ∈ U, XpK *ᵥ u ∈ U) : U = ⊥ ∨ U = ⊤ := by
  by_cases h : ∃ u ∈ U, u ≠ 0
  · right
    obtain ⟨u, hu, hu0⟩ := h
    have proj : ∀ w ∈ U, (w 0) • (![1, 0, 0] : Fin 3 → ℝ) ∈ U := by
      intro w hw
      have : (w 0) • (![1, 0, 0] : Fin 3 → ℝ) = (1 / 2 : ℝ) • (w + XpK *ᵥ w) := by
        rw [XpK_mulVec]; ext i; fin_cases i <;> simp [vecHead, vecTail] <;> ring
      rw [this]; exact U.smul_mem _ (U.add_mem hw (hK w hw))
    have h1 := proj u hu
    have h2 := proj _ (hT u hu)
    have h3 := proj _ (hT _ (hT u hu))
    rw [XpT_mulVec] at h2
    rw [XpT_mulVec, XpT_mulVec] at h3
    simp only [Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two,
      Matrix.head_cons, Matrix.tail_cons, neg_neg] at h2 h3
    have unit : ∀ c : ℝ, c ≠ 0 → c • (![1, 0, 0] : Fin 3 → ℝ) ∈ U → (![1, 0, 0] : Fin 3 → ℝ) ∈ U := by
      intro c hc hcU
      have := U.smul_mem c⁻¹ hcU
      rwa [smul_smul, inv_mul_cancel₀ hc, one_smul] at this
    have e0 : (![1, 0, 0] : Fin 3 → ℝ) ∈ U := by
      by_cases a0 : u 0 = 0
      · by_cases a1 : u 1 = 0
        · have a2 : u 2 ≠ 0 := by
            intro a2; apply hu0; ext i; fin_cases i <;> simp [a0, a1, a2]
          exact unit _ (neg_ne_zero.2 a2) h3
        · exact unit _ (neg_ne_zero.2 a1) h2
      · exact unit _ a0 h1
    have e2 : (![0, 0, 1] : Fin 3 → ℝ) ∈ U := by
      have := U.neg_mem (hT _ e0)
      convert this using 1; rw [XpT_mulVec]; ext i; fin_cases i <;> simp
    have e1 : (![0, 1, 0] : Fin 3 → ℝ) ∈ U := by
      have := hT _ e2
      convert this using 1; rw [XpT_mulVec]; ext i; fin_cases i <;> simp
    rw [eq_top_iff]
    intro v _
    have : v = v 0 • (![1, 0, 0] : Fin 3 → ℝ) + v 1 • ![0, 1, 0] + v 2 • ![0, 0, 1] := by
      ext i; fin_cases i <;> simp
    rw [this]
    exact U.add_mem (U.add_mem (U.smul_mem _ e0) (U.smul_mem _ e1)) (U.smul_mem _ e2)
  · left
    push Not at h
    exact (Submodule.eq_bot_iff U).2 h

/-! ### No invariant linear functional (Lemma 4, proof 2) -/

/-- A linear functional on ℝ³ invariant under `Xr` is zero; Xr is the action of r on both
blocks, so T_B has no r-invariant functional, in particular tr A ≡ 0 once tr A is r-invariant. -/
theorem no_invariant_functional (φ : (Fin 3 → ℝ) →ₗ[ℝ] ℝ) (h : ∀ v, φ (Xr *ᵥ v) = φ v) :
    φ = 0 := by
  have hv : ∀ v : Fin 3 → ℝ, Xr *ᵥ v = ![-v 2, -v 1, v 0] := by
    intro v; ext i; fin_cases i <;> simp [Xr, mulVec, dotProduct, Fin.sum_univ_three]
  have b : ∀ v : Fin 3 → ℝ, φ v = v 0 * φ ![1, 0, 0] + v 1 * φ ![0, 1, 0] + v 2 * φ ![0, 0, 1] := by
    intro v
    have : v = v 0 • (![1, 0, 0] : Fin 3 → ℝ) + v 1 • ![0, 1, 0] + v 2 • ![0, 0, 1] := by
      ext i; fin_cases i <;> simp
    nth_rewrite 1 [this]
    simp only [map_add, map_smul, smul_eq_mul]
  have h0 := h ![1, 0, 0]
  have h1 := h ![0, 1, 0]
  have h2 := h ![0, 0, 1]
  rw [hv, b, b ![1, 0, 0]] at h0
  rw [hv, b, b ![0, 1, 0]] at h1
  rw [hv, b, b ![0, 0, 1]] at h2
  simp at h0 h1 h2
  have z1 : φ ![0, 1, 0] = 0 := by linarith
  have z0 : φ ![1, 0, 0] = 0 := by linarith
  have z2 : φ ![0, 0, 1] = 0 := by linarith
  apply LinearMap.ext
  intro v
  rw [b v, z0, z1, z2]
  simp

/-! ### Real Schur lemma for involutions (Lemma 6(iv)) -/

/-- Irreducibility of a family of 3×3 real matrices. -/
def Irreducible3 (gs : List Mat3) : Prop :=
  ∀ U : Submodule ℝ (Fin 3 → ℝ), (∀ g ∈ gs, ∀ u ∈ U, g *ᵥ u ∈ U) → U = ⊥ ∨ U = ⊤

lemma mat_eq_of_mulVec {X Y : Mat3} (h : ∀ v, X *ᵥ v = Y *ᵥ v) : X = Y := by
  have : Matrix.toLin' X = Matrix.toLin' Y := LinearMap.ext fun v => by
    simp [Matrix.toLin'_apply, h v]
  exact Matrix.toLin'.injective this

/-- **Involution Schur lemma.** `Z² = 1` commuting with an irreducible family ⇒ `Z = ±1`. -/
theorem involution_schur (gs : List Mat3) (hirr : Irreducible3 gs) (Z : Mat3) (hZ2 : Z * Z = 1)
    (hcomm : ∀ g ∈ gs, Z * g = g * Z) : Z = 1 ∨ Z = -1 := by
  set U := LinearMap.ker (Matrix.toLin' (Z - 1)) with hU
  have hmem : ∀ u, u ∈ U ↔ Z *ᵥ u = u := by
    intro u
    simp [hU, Matrix.toLin'_apply, sub_mulVec, sub_eq_zero]
  have hinv : ∀ g ∈ gs, ∀ u ∈ U, g *ᵥ u ∈ U := by
    intro g hg u hu
    rw [hmem] at hu ⊢
    rw [mulVec_mulVec, hcomm g hg, ← mulVec_mulVec, hu]
  rcases hirr U hinv with hbot | htop
  · right
    apply mat_eq_of_mulVec
    intro v
    -- (Z + 1) v ∈ U, hence 0
    have h1 : (Z + 1) *ᵥ v ∈ U := by
      rw [hmem, mulVec_mulVec, mul_add, hZ2, mul_one, add_comm]
    rw [hbot, Submodule.mem_bot, add_mulVec, one_mulVec] at h1
    rw [neg_mulVec, one_mulVec]
    exact eq_neg_of_add_eq_zero_left h1
  · left
    apply mat_eq_of_mulVec
    intro v
    have : v ∈ U := by rw [htop]; exact Submodule.mem_top
    rw [hmem] at this
    rw [this, one_mulVec]

/-- **Lemma 6(iv), E side.**  ρ_E(J) is a central involution, ρ_E(a)² = ρ_E(J), and every
generator g satisfies g a g⁻¹ ∈ {a, aJ}.  If E is irreducible then ρ_E(a) = ±1. -/
theorem rhoE_a_scalar (gs : List Mat3) (hirr : Irreducible3 gs) (ρa ρJ : Mat3)
    (hJ2 : ρJ * ρJ = 1) (hJc : ∀ g ∈ gs, ρJ * g = g * ρJ) (ha2 : ρa * ρa = ρJ)
    (hac : ∀ g ∈ gs, g * ρa = ρa * g ∨ g * ρa = ρa * ρJ * g) : ρa = 1 ∨ ρa = -1 := by
  rcases involution_schur gs hirr ρJ hJ2 hJc with hJ | hJ
  · subst hJ
    refine involution_schur gs hirr ρa ha2 fun g hg => ?_
    rcases hac g hg with h | h
    · exact h.symm
    · rw [mul_one] at h; exact h.symm
  · exfalso
    have hdet := congrArg Matrix.det ha2
    rw [Matrix.det_mul, hJ, Matrix.det_neg, Matrix.det_one] at hdet
    simp at hdet
    nlinarith [sq_nonneg ρa.det]

/-! ### Lemma 7: the F_B line -/

/-- In block coordinates, `{v = (y,z) : r v = −v, κ v = v}` is the line of `(e₁, 0)`. -/
theorem line_r_smid (y z : Fin 3 → ℝ) (hr1 : Xr *ᵥ y = -y) (hr2 : Xr *ᵥ z = -z)
    (hk1 : XmK *ᵥ y = y) (hk2 : XpK *ᵥ z = z) : y = ![0, y 1, 0] ∧ z = 0 := by
  have a1 := congrFun hr1 0; have a2 := congrFun hr1 2
  have b1 := congrFun hr2 0; have b2 := congrFun hr2 2
  have c := congrFun hk1 2
  have d1 := congrFun hk2 1; have d2 := congrFun hk2 2
  simp [Xr, XmK, XpK, mulVec, dotProduct, Fin.sum_univ_three] at a1 a2 b1 b2 c d1 d2
  constructor
  · ext i; fin_cases i <;> simp <;> linarith
  · ext i; fin_cases i <;> simp <;> linarith

end Bolza
