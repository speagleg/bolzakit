# Second reading of the Lemma N proof (review_theory_1.md §3)

Reader: independent agent, 2026-10-08. I did not write Lemma N, F1_PROOF or review_theory_1.

Sources re-downloaded and read with pdftotext, then grepped:
- Mazzoli arXiv:2101.09380 (§1 and Remark 3.3);
- El Soufi–Ilias arXiv:math/0701777 (Lemma 2.1).

Tromba's book, Fischer–Marsden 1975 and Earle–Eells were **not read**.

## Verdict: **CORRECT WITH FIXES**
- The statement of Lemma N is correct. Steps 1, 2, 3 and 5 are correct as written, up to cosmetic points.
- Step 4 has one real gap (N1). It assumes a C¹ family k_t of smooth hyperbolic metrics with [k_t] = [h_t], taken
  from "any Teichmüller slice". That needs t ↦ [h_t] to be C¹ in the slice's smooth structure, i.e. smoothness of π
  into T₂. This is the very kind of fact G3 is about, and it is neither proved nor cited.
- Fix 1 below removes it with an explicit smoothing, so that no slice, no lift of conjugation and no principal-bundle
  fact are needed.
- With Fix 1, the only external inputs are:
  - elliptic regularity;
  - the Fischer–Marsden/Tromba splitting at the single smooth metric h₀ (now QUOTED, Fix 2);
  - Riemannian centre of mass (standard).

## Step-by-step

| Step | Claim | Verdict | Notes |
|---|---|---|---|
| 1 | b_ḣ(e,e′) = −∫⟨ḣ,T(e,e′)⟩dA, T = Sym(de⊗de′) − ½(⟨de,de′⟩ − λ₀ee′)h₀ | **CORRECT** | Re-derived: d(h^{ij}) = −h^{ik}ḣ_{kl}h^{lj}, d√h = ½(tr ḣ)√h, and tr ḣ = ⟨ḣ,h₀⟩. **Independent confirmation (QUOTED, re-typed from pdftotext, El Soufi–Ilias Lemma 2.1):** "Q_h(u) := ∫_M u P_{k,h}u v_g = −∫_M (du ⊗ du + ¼ ∆_g u² g, h) v_g", for analytic deformations with h = d/dt g_t. Since Δ(u²) = 2uΔu − 2\|du\|² = 2λu² − 2\|du\|², we get du⊗du + ¼Δ(u²)g = T(u,u). So the two formulas agree. Cosmetic: "reference-simplex coordinates, in which h₀ is smooth" is inaccurate. h₀ is only cellwise analytic in those coordinates. The intended setting (Step 4) is h₀ = the Bolza metric on B and h_t = Φ_t*h_{B_t} with Φ₀ = id. The formula is tensorial and can be evaluated cellwise, so this changes nothing. Agreement with F1's b_v: q_{h_t} equals F1's q_t by (H3), and t ↦ h_t is C¹ into L^∞ (cellwise analytic, jointly in (t,p)) ✓ |
| 2 | div T = 0 for e, e′ ∈ E | **CORRECT** | Re-derived in index form: 2∇^iS_ij = (tr∇²e)e′_j + (tr∇²e′)e_j + ∇²e′(∇e,·)_j + ∇²e(∇e′,·)_j, and div(φh₀) = dφ. The Hessian terms cancel against d⟨de,de′⟩, and the λ₀ terms cancel because tr∇²e = −λ₀e. This is the standard "stress–energy tensor is divergence-free on-shell" |
| 3 | b_{L_Xh₀} = 2∫⟨X, div T⟩ = 0 for Lipschitz X | **CORRECT** | L_Xh₀ = 2Sym∇X♭ ∈ L^∞ a.e. by Rademacher. Integration by parts of a W^{1,∞} field against a smooth tensor on a closed manifold is valid, and T is symmetric. Because b is linear in ḣ ∈ L^∞, no flow and no H¹-differentiability of e∘ψ_s⁻¹ is needed. This correctly removes the defect (i) of the original sketch |
| 4 | ḣ = σ + L_{Y}h₀ + L_Xh₀ with σ TT, Y smooth, X Lipschitz | **GAP (N1), and minor points N2–N4**; fixed below | See below |
| 5 | A(g·σ) = ρ_E(g)A(σ)ρ_E(g)ᵀ for all g ∈ Isom(B), both orientations | **CORRECT** | Re-derived: b_{(g⁻¹)*σ}(e,e′) = −∫⟨σ, g*T(e,e′)⟩ = −∫⟨σ, T(e∘g, e′∘g)⟩ = b_σ(ρ(g)⁻¹e, ρ(g)⁻¹e′). With ρ_E orthogonal in the m₀-orthonormal basis, this gives ρ(g)A(σ)ρ(g)ᵀ. Orientation plays no role, since dA is used as a measure. The convention (g⁻¹)* versus g* (CERT_P6 uses pull-back, an anti-homomorphism) does not affect the real characters, because χ(g⁻¹) = χ(g). Minor (N5): "it is also the derivative of the mapping-class action that the octagon relabellings induce" is asserted, not shown. See Fix 3 |

### Step 4 in detail
**N1 (gap): existence of k_t.**
- The step needs a C¹ curve k_t of smooth hyperbolic metrics with [k_t] = [h_t] as marked surfaces. Taking
  k_t := slice(s(t)) requires the coordinate s(t) of π(x(t)) in the slice to be C¹ in t.
- That is: π: C → T₂ must be C¹ for the slice's smooth structure. Equivalently, the holonomy map from Tromba's T₂ to
  the Fuchsian character variety must be a diffeomorphism (Weil/Goldman/Earle–Eells type), and π is smooth into the
  character variety.
- This is true and standard, but it is not cited. It is also not automatic from "Tromba's slice exists".
- My G3 certificate alone does not supply it either. It shows that the τ̄_w are topological coordinates, but it does
  not show that the slice's τ̄_w-Jacobian is invertible.

**N2 (minor): the PSL(2,ℝ)-lift.** "Conjugation acts freely and properly on Fuchsian representations, so the quotient
map is a principal bundle" is more than is needed. A direct argument suffices:
- c_t is unique, because the centraliser of a non-elementary group is trivial.
- c_t is C¹ by the implicit function theorem applied to c ↦ (cρ_t(a)c⁻¹, cρ_t(b)c⁻¹), for two non-commuting
  hyperbolic generators a, b. This map is an injective immersion of PSL(2,ℝ).

This point disappears under Fix 1.

**N3 (minor): regularity.** "φ_t jointly C¹ on each closed cell" is not quite enough to differentiate φ_t*k_t,
which involves dφ_t. One needs ∂_t∇_pφ to exist and be continuous on each closed cell.
- It does hold here: D_t is jointly analytic per cell, c_t is C¹ in a finite-dimensional Lie group, and D′_t is
  smooth.
- It should be stated.
- "Piecewise C¹ on finitely many closed polygonal cells + continuous ⇒ Lipschitz" ✓.

**N4 (citation): the splitting k̇ = σ + L_Yh₀.** The Mazzoli sentence quoted in CERT_P6 §1 only *characterises*
S²_tt as T_hH ∩ (T_h(Diff₀·h))^⊥. It does not by itself give the direct-sum splitting with smooth Y.

The splitting is QUOTED from Mazzoli arXiv:2101.09380, Remark 3.3, citing Fischer–Marsden [FM75]: "Given g a hyperbolic
metric, every symmetric tensor S admits an orthogonal decomposition of the following form: S = Stt + LX g +
((−∆g f + f ) g + g∇2 f ), where: • Stt is transverse traceless with respect to g; • Stt + LX g is tangent to the
space of Riemannian metrics with constant Gaussian curvature equal to −1. In other words, … Stt + LX g ∈ ker dKg ;
• (−∆g f + f ) g + g∇2 f lies in the L2 -orthogonal of ker dKg."

Applied to S = k̇ ∈ ker dK_{h₀}: the third component lies in ker dK ∩ (ker dK)^⊥ = 0, so k̇ = σ + L_Yh₀.
- Smoothness of Y (indeed of σ, which is Re of a holomorphic quadratic differential) follows from elliptic
  regularity. The quote does not state it.
- Cite it this way, rather than "Fischer–Tromba decomposition (quoted in CERT_P6 §1)".

### Fix 1 (replaces the slice in Step 4; PROVED here)
Let D_t: H² → H² be the developing map of h_t := Φ_t*h_{B_t}, with Φ₀ = id. It is ρ_t-equivariant, cellwise
analytic jointly in (t,p), and D₀ = id after normalisation. Fix ε > 0 small, and a smooth radial kernel
η_ε(d_{h₀}(p,y)) supported in an ε-ball. Define the Riemannian centre of mass

  Ψ̃_t(p) := argmin_{q∈H²} ∫ d(q, D_t(y))² η_ε(d(p,y)) dA(y).

Its properties:
- **Well defined.** The functional is strictly convex on H² and proper.
- **Equivariant.** Ψ̃_t(γp) = ρ_t(γ)Ψ̃_t(p), because the kernel is deck-invariant and D_t is equivariant.
- **Smooth.** It is jointly C^∞ in (t,p), by the IFT on the gradient equation: p-derivatives fall on the smooth
  kernel, and t-derivatives fall on the cellwise analytic D_t.
- **Ψ̃₀ = id.** D₀ is an isometry and the kernel is symmetric, so the centre of mass of the pushed-forward measure is
  the image of the centre.
- **A diffeomorphism for small |t|.** dΨ̃_t = I + O(t) uniformly, so Ψ̃_t is an equivariant local diffeomorphism and
  proper, hence a diffeomorphism.

Set Ψ_t: B → B_t (descended), k_t := Ψ_t*h_{B_t} and φ_t := Ψ_t⁻¹∘Φ_t. Then:
- k_t is a jointly smooth family of smooth hyperbolic metrics with k₀ = h₀;
- φ_t: (B,h_t) → (B,k_t) is an isometry, isotopic to the identity (it commutes with the deck group, so it is the
  identity on π₁);
- φ₀ = id, φ_t is jointly analytic×smooth on each closed cell and continuous, so X := ∂_tφ_t|₀ is Lipschitz;
- ḣ = k̇ + L_Xh₀ a.e.

Steps 1, 3 and the splitting of N4 then give A_oct(v) = b_σ with σ the TT part of k̇. σ(v) is independent of the
choice of ε and kernel: two choices of Ψ_t differ by Diff₀, so k̇ changes by a Lie derivative.

This is exactly the original F1_PROOF sketch's Ψ_t, now constructed. The reviewer's Step 3 makes the Lipschitz flow
harmless. **The map v ↦ σ(v) is what "dπ" means in Lemma N(a)**. It is onto TT by the G3 certificate
(CERT_G3.md §3(d)3).

### Fix 3 (N5: relabellings)
For an octagon symmetry γ (r, s_mid, …) acting linearly on C, the fan mesh is γ-equivariant (CERT_P6 §4.1:
centre = normalised vertex sum, and hyp_interp is symmetric and isometry-equivariant). So Φ_{γx} = Γ∘Φ_x∘γ_B⁻¹,
with Γ the hyperbolic isometry of the disk and γ_B the induced isometry of B. Hence:
- h_{γx(t)} = (γ_B⁻¹)*h_{x(t)};
- ḣ_{γv} = (γ_B⁻¹)*ḣ_v;
- σ(γv) = γ_B·σ(v), since the TT projection is natural under isometries of h₀.

This is the identity "r·w_B = −w_B" used in CERT_P6 §2.5. One line in the write-up suffices.

## Summary of required edits to review_theory_1 §3 / F1_PROOF §4
1. Step 4: replace "any Teichmüller slice … C¹ lifts exist" by Fix 1. Alternatively, cite explicitly that π is smooth
   into T₂ (holonomy ↔ Tromba T₂ diffeomorphism), which is a citation the project has not read.
2. Cite the splitting as Mazzoli Rem. 3.3 / Fischer–Marsden (quoted above), plus elliptic regularity for Y.
3. Step 1: drop "reference-simplex coordinates, in which h₀ is smooth"; take h₀ = Bolza metric and Φ₀ = id.
4. Optionally cite El Soufi–Ilias Lemma 2.1 as an independent source for the first-variation formula.
5. Add Fix 3 (one line) for the relabelling equivariance used in CERT_P6 §2.5.

With these, **Lemma N is PROVED**. Its residual inputs are elliptic regularity, the Fischer–Marsden splitting
(quoted), and the existence of the Riemannian centre of mass in H² (standard; not quoted).
