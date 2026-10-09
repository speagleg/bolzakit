# Review (theory 1): CERT_P6.md + cert_chars.py, F1_PROOF.md

Reviewer: independent adversarial second reader (agent), 2026-10-08.
Scope: the representation theory (CERT_P6, cert_chars.py), F_B ⊂ K, the certified character of E, Lemma N, and
F1_PROOF (including its quoted sources).
Not in scope: CERT_RESULTS §1–2 (LP, Rayleigh–Ritz) and §4 (A ≠ 0), beyond how they are used here.

Independent evidence produced (files in `work/scratch_RT/`):
- `indep_chars.py` / `indep_chars.out`. A new model, sharing no code with cert_chars.py. It works with
  Ḡ = Isom(B)/⟨J⟩ ≅ O_h, the 48 (anti-)Möbius maps permuting W = {0, ∞, ±1, ±i}.
  - Holomorphic g: pull back φ = p(x)dx²/y² directly, g*φ = (det M)²λ⁻²(cx+d)²p(Mx)·dx²/y², exact in sympy.
  - Antiholomorphic g: φ ↦ conj(g*φ).
  - No Sym² formula is used, and only λ² enters. This also shows that J acts trivially on Q(B).
- `fb_trig.py` / `fb_trig.out`: an mpmath check of the F_B length derivative.
- Sources were re-downloaded from arXiv into a fresh directory, and the quotes were grepped from the pdftotext
  output.

**Overall verdict.** The theory is sound.
- No mathematical error was found in any claim that the main theorem uses.
- Lemma N is correct. Its "sketch" status is a real but closable gap. A complete proof is given in §3 below.
- There is one further gap, which is easy to close: F1 is verified only for the 11-dimensional octagon chart. The
  local-max Lemma on T₂ also needs F1 on a neighbourhood in T₂, which in turn needs either a certified rank 6 for
  dπ or a Teichmüller slice. See G3.

## 1. Table of claims

| # | Claim (source) | Verdict | Evidence |
|---|---|---|---|
| 1a | G = Isom(B) has order 96. Generators r, t, κ of y² = x⁵−x; G⁺ ≅ GL(2,3) (CERT_P6 §1) | CONFIRMED | Independent construction: exactly 24 Möbius and 24 anti-Möbius maps permute W (= O_h, order 48). Each lifts twice (±λ), giving 96. The curve condition (ax+b)(cx+d)[(ax+b)⁴−(cx+d)⁴] = λ²(x⁵−x) was re-derived by hand. The pull-back matrix P = (det M/λ)[[d,b],[c,a]] (columns = images) was re-derived by hand and is correct |
| 1b | T_B ≅ Re Q(B) as a real G-module; G⁺ acts C-linearly, reversing elements act C-antilinearly by φ ↦ conj(g*φ) (§1) | CONFIRMED | The Mazzoli quote was checked verbatim (arXiv:2101.09380 p. 400–402 of the text: "Transverse traceless tensors are exactly those 2-tensors that can be written as Re φ…"; "It is shown in [Tro92] that, for every hyperbolic metric h, S²_tt(Σ, h) coincides with…"). For antiholomorphic w = g(z), g*(φ dw²) = φ(g)(∂̄g)²dz̄², whose real part is Re of the holomorphic differential conj(φ∘g)·conj(∂̄g)²dz². So the action is antilinear, and its real trace is 0 (J⁻¹TJ = −T) |
| 1c | χ_T: r, r³ ↦ −2; r² ↦ −2; J ↦ 6; t ↦ 0; every reversing g ↦ 0 (§0 item 1) | CONFIRMED (two independent methods) | (i) Direct pull-back (indep_chars.out): r −2, r² −2, t 0, u 2, s_mid 0, s_vert 0, antipodal 0; ⟨χ_T,χ_T⟩ = 2, ⟨χ_T,1⟩ = 0. (ii) Orbifold dimension count (consistency check that uses the F3-type fact dim T_B^H = dim Teich(B/H) = 2(k−3) for a genus-0 quotient with k cone points), via (1/\|H\|)Σχ_T: ⟨J⟩: k = 6 ⇒ 6 ✓; ⟨r²⟩ ∋ J: k = 4 (0, ∞, {±1}, {±i}) ⇒ 2 = (12+2χ(r²))/4 ⇒ χ(r²) = −2; ⟨r⟩: k = 3 ⇒ 0 = (12+2χ(r²)+4χ(r))/8 ⇒ χ(r) = −2; ⟨t,J⟩: k = 4 (2 W-orbits + 2 face centres) ⇒ 2 = (12+4χ(t))/6 ⇒ χ(t) = 0; ⟨u,J⟩: k = 5 ⇒ χ(u) = 2. A reflection's fixed locus is the 3-dimensional real Teichmüller space ⇒ χ(s) = 0 |
| 1d | Octagon r ↔ order-8 element; s_mid ↔ the reflection with srs = r⁻¹ that fixes 4 Weierstrass points; s_vert fixes 2 (§1) | CONFIRMED, with a cosmetic wording error | Geometric check: s_mid is θ ↦ π/4−θ, which sends midpoint class k to −k. It fixes classes M₀ (on the axis) and M₂ (M₂ ↦ M₆ ~ M₂; *not on the axis*), plus the centre and the vertex: 4 points. The phrase "two midpoint classes on its axis" is wrong: only one class lies on the axis. The count is right. s_vert sends k ↦ −k−1 and fixes no class: 2 points. In O_h terms s_mid = σ_h (coordinate plane, x ↦ x̄) and s_vert = σ_d (x ↦ i x̄). All pairs (C4, σ_h ∋ axis) are conjugate in O_h, so nothing depends on which representatives are chosen |
| 2a | T_B = K₀ ⊕ W₀, non-isomorphic and absolutely irreducible; T_B^G = 0; K₀ = χ₇, W₀ = χ₁₀ (§2.1–2.2) | CONFIRMED | A new, more transparent description: the antipodal isometry a (lift of x ↦ −1/x̄, Cook class 4a) is central mod J. On T_B, a² = J = id, so T_B = T₊ ⊕ T₋ (its ±1-eigenspaces), each of rank 3 (exact). ⟨χ±,χ±⟩ = 1 and ⟨χ₊,χ₋⟩ = 0. In O_h language W₀ = T₊ = T₂g and K₀ = T₋ = T₂u. Values: T₋ has s_mid +1, s_vert −1, a −3; T₊ has s_mid −1, s_vert +1, a +3. These match Cook X.7 and X.10 (re-read from arXiv:2108.11825 App. A; the rows quoted in CERT_P6 §4.3 are verbatim) |
| 2b | The four 3-dim irreducibles are X₁, X₂, O, Oε; every 3-dim real irreducible E has Sym₀²E ≅ D ⊕ W₀, with no K₀ and no trivial summand (§2.2–2.3) | CONFIRMED | The Clifford/Gallagher argument is correct (odd dimension ⇒ real type, irreducible on G⁺; two extensions each). Independently, with E = vec ⊗ {1, det, sgn_S₄, det·sgn_S₄} (T₁u, T₁g, T₂u, T₂g): ⟨E,E⟩ = 1 and Sym₀²: norm 2, mult(T₊) = 1, mult(T₋) = 0, mult(1) = 0 in all four cases. This is the textbook O_h fact Sym₀²T₁ = Sym₀²T₂ = E_g ⊕ T₂g. Note: on O_h, orientation = det, so twisting by "ε" alone gives only the two T₁'s; the S₄-sign twist is needed for the T₂'s (cert_chars handles this correctly via its 4 linear characters) |
| 2c | Schur: tr A ≡ 0, A\|K₀ = 0, A ≠ 0 ⇒ ker A = K₀, rank 3 (§2.4) | CONFIRMED (given Lemma N) | Hom_G(K₀⊕W₀, 1⊕D⊕W₀) = End_G(W₀) = ℝ. **Shorter proof of K₀ ⊂ ker A (suggested for the paper):** E is real-irreducible of odd dimension, so the central J acts by ±1. If J = −1, then a² = −1 on E would be a complex structure on ℝ³, which is impossible. So J = 1 on E, and a acts by a scalar ±1. Then A(a·v) = ρ(a)A(v)ρ(a)⁻¹ = A(v), so for v ∈ K₀ = {a·v = −v} we get A(v) = −A(v) = 0. This needs no character table |
| 3a | The line {r = −1, s_mid = +1} in T_B is 1-dim and lies in K₀ (§2.5) | CONFIRMED | indep_chars: nullspace dimension 1; a·v = −v, so v ∈ T₋ = K₀. Restriction to C₄v: mult(B₁) = 1 in T₂u and 0 in T₂g. \|Stab_G(v)\| = 16, and u: x ↦ i/x fixes v |
| 3b | w_B := dπ(dF_B/dp) lies on that line (§2.5) | CONFIRMED | F_B(p) (equiangular, sides alternating) is invariant under the π/2 rotation and the bisector reflection of side 0, and r swaps p ↔ −p. Relabellings commute with π |
| 3c | w_B ≠ 0 by hyperbolic trigonometry (§2.5) | CONFIRMED | All formulas re-derived: right-triangle laws cosh c = cot A cot B and cos A = cosh a sin B; dθ′ = −dθ at π/8. Closed form: dℓ(g₀)/dp = −cos(π/4)/(sin²(π/8) sinh r_in) = **−sinh r_in = −√(2+2√2)** = −2.19736822693562 (fb_trig.out, mpmath with 30 digits, implicit solve). g₀ is the translation along M₄OM₀ by 2h_a, by the central and bisector symmetry; at p = 0, ℓ = 2arccosh(1+√2) ✓ |
| 3d | d²Λ/dp²(F_B) = c·g_inv(w_B,w_B), with no rank-6 condition needed (§2.5) | CONFIRMED | Chain rule plus dΛ(B) = 0. Λ is G-invariant, so D²Λ is G-invariant; K₀ is absolutely irreducible, so D²Λ\|K₀ = c·g_inv |
| 4 | tr(g\|E) ∈ [1 ± 0.0352] for g = r, s_mid, s_vert (n = 96), so E = χ₈ (§4) | CONFIRMED (rigorous, conditional on the cited inputs) | Bound re-derived: m(e_i, g e_j) = H − N, because the cross terms vanish (g preserves E and (E⊕1)^⊥). ‖H‖ ≤ 1 (compression of a unitary). ‖N‖₁ = ‖XᵀUX‖₁ ≤ ‖X‖²_HS = tr W. W ≤ (K̃ − λI)/(λ₄ − λ), from q(w) ≥ λ₄m(w) on (E⊕1)^⊥; the worst case over λ ∈ [λ_lo, λ_hi] gives δ. ‖(I−W)⁻¹W‖₁ ≤ δ/(1−δ). Total ≤ 2δ/(1−δ) ✓. Using G⁻¹ instead of orthonormalising is a congruence, so it is fine. **Exact g-invariance of the FE space** (needed for g·u = permuted vector) was checked in the formulas, not just numerically: the centre is the normalised vertex sum, `hyp_interp(A,B,s) = hyp_interp(B,A,1−s)` is isometry-equivariant, the chord nodes T(k,i,j) ↦ T(−k−1,i,i−j) under s_vert, and barycentric P1 on the hyperboloid is linear-isometry equivariant. Inputs: dim E = 3, λ₁ ∈ SU13 interval, λ₄ ≥ 4.0 (CERT_RESULTS rows 1b, 3). The main argument does not need item 4 (§2.3) |
| 5 | Lemma N (F1_PROOF §4) | **GAP (closable); statement CORRECT** | The sketch has the right idea, but: (i) it differentiates s ↦ e∘ψ_s⁻¹ in H¹ for Lipschitz flows without justification; (ii) the existence of a C¹ family Ψ_t and the C¹/Lipschitz regularity of ψ_t = Ψ_t⁻¹Φ_t are asserted, not shown; (iii) the TT decomposition of a smooth tangent is quoted. A complete proof that avoids (i) is in §3 below |
| 6a | F1 Theorem, Steps 1–6 (F1_PROOF §2) | CONFIRMED (minor remarks) | Each step checked. Resolvent expansion; uniform Neumann inversion on compact Γ; P² = P for complex t (resolvent identity T(z)−T(w) = (w−z)T(z)T(w) with T = RM, or the identity theorem on ℝ^d); E4 constant rank; second contour giving exactly λ₁..λ₃ (needs no spectrum in (0, λ₁), certified); M′ = q′ − λ₀m′. Remarks: Step 1's "norm ≤ 1 in the H¹ norm of h₀" holds only for the reference-atlas norm (otherwise a constant), which is harmless. (H1) for the octagon: D = \|det\| is analytic because det ≠ 0 near x_B; the arclength interpolation sinh(s·L)/sinh L is analytic in cosh L ✓ |
| 6b | Octagon family satisfies (H1)–(H3) with k = ω; hence Λ is analytic on the chart C | CONFIRMED | As above |
| 6c | Use of F1 for the local-max Lemma on T₂ | **GAP (G3)** | F1 is verified on the 11-dim C. The Lemma needs M(v) for v in a neighbourhood of B in T₂, i.e. a local analytic section of π: C → T₂. That requires rank dπ(x_B) = 6, which is only numerical (singular values 6.4e2 vs ≤ 1e-5, review_1). Fix: certify the 6×11 trace-coordinate Jacobian rank in Arb (trivial margin), or apply F1 to an analytic Teichmüller slice of smooth metrics (H1 is immediate there) and transport A via Lemma N. CERT_RESULTS §4 does not need this: A_oct(δ) ≠ 0 already implies dπδ ≠ 0 by Lemma N(a) |
| 6d | Ky Fan reduction from min-max (§5) | CONFIRMED | λ_j ≤ μ_j via V_j; equality at eigenvectors. The P5 inequality D²Λ ≤ g″(0) is correct because g − Λ ≥ 0 is C², vanishes at 0, and so has zero derivative and nonnegative second derivative there |
| 6e | Cai arXiv:2608.02459 Prop. 2.5 quotes | CONFIRMED verbatim (v2) | "For every t0 ∈ R^N, there exists a complex neighborhood U ⊂ C^N of t0 such that … is a holomorphic family of closed sectorial forms with common domain; in Kato's terminology, it is a family of type (a)." and "the corresponding Riesz projection is operator-norm holomorphic in z, and therefore real analytic for real t; its rank is locally constant." Caveat: F1_PROOF's "that is the same structure as ours" overstates it. Cai's family is a *fixed* Dirichlet form E(m_z u, m_z̄ v) with a multiplicative weight, conjugated by a unitary to a fixed L²(µ₀). It is not a variable-matrix-coefficient pencil. Cite it as an analogue ("the same argument"), not as covering our case. His frame uses H_Γ, B_Γ and the characteristic polynomial, not the G^{−1/2}-orthonormalised M |
| 6f | Duran arXiv:2507.16895 §4.1 quotes | CONFIRMED verbatim (v3, ref [31] = Kato 1995; in v2 the reference is [27]) | Add "v3" to the citation |
| 6g | Moslehian arXiv:1108.1467 Ky Fan quote | CONFIRMED (verbatim up to the omitted lead-in "They assert that"; attribution [21, 1950] = PNAS 36, 31–35; eigenvalues in decreasing order ✓) | |
| 7 | Cook Thm 3.9 and App. A rows X.7–X.10 (CERT_P6 §4.3) | CONFIRMED verbatim | Re-read from arXiv:2108.11825, pdftotext -layout. My O_h identification matches on all classes: 4a = antipodal lifts (X.8 −3 = T₁u, X.10 +3 = T₂g), 2a = σ_d = s_vert, 2d = σ_h = s_mid |
| 8 | §2.6 (λ₂ ≡ λ₃ along F_B) | NOT REVIEWED in depth (not needed) | The exact parts (\|Stab\| = 16, u ∈ Stab) were reproduced |

## 2. Remarks for the write-up
- **Cleaner presentation.** The whole of §2 can be stated through the antipodal isometry a, the order-4 orientation-
  reversing lift of x ↦ −1/x̄:
  - K = (−1)-eigenspace of a on T_B, and W = (+1)-eigenspace;
  - K ⊂ ker A for every odd-dimensional irreducible E (argument in row 2c);
  - ker A = K once A ≠ 0, because W is absolutely irreducible;
  - F_B ⊂ K because a·w_B = −w_B.

  Only "W is absolutely irreducible" needs a character computation (norm 1, exact).
- Fix the wording in CERT_P6 §1 (s_mid's fixed midpoint classes: one on the axis, one fixed via the gluing).
- Cite Duran as v3.

## 3. Lemma N: complete proof (REVIEWER'S PROOF, not in the original files)
Notation:
- h₀ is the hyperbolic metric of B and dA its area.
- Δ ≥ 0, and E is the λ₀-eigenspace.
- For an L^∞ symmetric 2-tensor ḣ, b_ḣ := d/ds|₀ (q − λ₀m)_{h₀+sḣ} restricted to E × E. Here q_h(f,g) =
  ∫ h^{ij}∂_if∂_jg √det h, and m_h(f,g) = ∫ fg √det h, written in any smooth atlas of B. For the octagon these are
  the reference-simplex coordinates, in which h₀ is smooth.

**Step 1 (first-variation formula).** From d(h^{ij}√h) = (−h^{ik}ḣ_{kl}h^{lj} + ½(tr ḣ)h^{ij})√h and
d√h = ½(tr ḣ)√h:

  b_ḣ(e,e′) = −∫_B ⟨ḣ, T(e,e′)⟩ dA,  T(e,e′) := Sym(de⊗de′) − ½(⟨de,de′⟩ − λ₀ee′)h₀.

This is linear and continuous in ḣ ∈ L^∞. It agrees with F1's b_v = ∂_v q − λ₀∂_v m, because F1's coefficients
are exactly a = √h h⁻¹ and ρ = √h in that atlas.

**Step 2 (T is divergence-free on E).** With ∇² the Hessian (symmetric) and tr∇²e = −λ₀e:

  2 div T = (tr∇²e)de′ + (tr∇²e′)de + ∇²e(∇e′,·) + ∇²e′(∇e,·) − d⟨de,de′⟩ + λ₀(e de′ + e′ de) = 0,

since d⟨de,de′⟩ = ∇²e(∇e′,·) + ∇²e′(∇e,·). (e and e′ are smooth by elliptic regularity.)

**Step 3 (Lie derivatives drop out).** For X ∈ W^{1,∞} (Lipschitz) vector fields, L_Xh₀ = 2Sym∇X♭ ∈ L^∞, and

  b_{L_Xh₀} = −2∫⟨∇X♭, T⟩ = 2∫⟨X, div T⟩ = 0.

Integration by parts is valid for X ∈ W^{1,∞} and smooth T on a closed manifold. This replaces the flow
argument, which needed H¹-differentiability of e∘ψ_s⁻¹.

**Step 4 (octagon tangent = TT + Lipschitz Lie derivative).** Let h_t := Φ_t*h_{B_t}. This is a piecewise-analytic
L^∞ family on B with h₀ the Bolza metric, and ḣ := ∂_th_t|₀ ∈ L^∞. Let k_t be a C¹ family of smooth hyperbolic
metrics on B, with k₀ = h₀ and [k_t] = [h_t] in T₂ as marked surfaces. Any real-analytic Teichmüller slice works,
e.g. Tromba's [Tro92]; this is the one standard input.
- Let φ_t: (B,h_t) → (B,k_t) be the isometry isotopic to the identity. It is unique, since genus ≥ 2.
- Lift to H². With developing maps D_t (explicit: the mesh map on the octagon, jointly piecewise-analytic in (t,p))
  and D′_t (smooth in (t,p)), with holonomies ρ_t and ρ′_t = c_tρ_tc_t⁻¹:
  φ̃_t = D′_t⁻¹ ∘ c_t ∘ D_t.
- c_t ∈ PSL(2,ℝ) can be chosen C¹ in t. Conjugation acts freely and properly on Fuchsian representations, so the
  quotient map is a principal bundle and C¹ lifts exist.
- So φ_t is jointly C¹ on each closed mesh cell and continuous across cells. Hence X := ∂_tφ_t|₀ is piecewise
  smooth and continuous, i.e. Lipschitz, and φ₀ = id.
- Differentiating h_t = φ_t*k_t cellwise gives ḣ = k̇ + L_Xh₀ a.e.
- By the Fischer–Tromba decomposition (quoted in CERT_P6 §1 via Mazzoli), k̇ = σ + L_Yh₀ with σ = Re φ TT and Y
  smooth, where σ is the TT representative of [k̇] = dπ(v).
- By Steps 1 and 3, A_oct(v) = [b_σ(e_i,e_j)] =: A(dπ v).

So A is well defined on T_B ≅ {TT tensors}, and it equals F1's dM(0).

**Step 5 (equivariance).** For an isometry g of h₀ (either orientation), with g acting by (g⁻¹)* on TT tensors and
ρ_E(g)e = e∘g⁻¹: T(e∘g⁻¹, e′∘g⁻¹) = (g⁻¹)*T(e,e′), and ∫⟨(g⁻¹)*σ, (g⁻¹)*T⟩dA = ∫⟨σ,T⟩dA, because dA is
g-invariant as a measure. Hence A(g·σ) = ρ_E(g)A(σ)ρ_E(g)ᵀ.
- This is the action of CERT_P6 §1: (g⁻¹)*Re φ = Re((g⁻¹)*φ) for holomorphic g, and Re(conj((g⁻¹)*φ)) for
  antiholomorphic g.
- It is also the derivative of the mapping-class action that the octagon relabellings induce at the fixed point B.
∎

Residual citations, all standard: elliptic regularity; the existence of an analytic (or C¹) slice of hyperbolic
metrics over T₂ together with the TT decomposition (Tromba); and smooth local sections for PSL(2,ℝ)-conjugation on
Fuchsian representations. With these cited, Lemma N is PROVED.
