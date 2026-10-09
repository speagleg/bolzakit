# PROOF: the Bolza surface is a strict local maximum of λ₁ on genus-2 Teichmüller space

Integrated 2026-10-08 from the item-23 working files. This document **supersedes** RESULTS.md, CERT_RESULTS.md,
CERT_SIGN.md, CERT_P6.md, F1_PROOF.md and CERT_G3.md for the final statement. Those files remain the detailed
records, and every number below is copied from them, from their output files, or from a review. Nothing was
recomputed for this document.

**Status of this document.** Every component lemma has been second-read; Lemma 9 (CERT_G3) and the OPEN-1 closure were second-read in review_proof_1 (ACCEPT-WITH-FIXES). The theorem is PROVED (computer-assisted), conditional on the external inputs of §Dependencies.

**Update 2026-10-09 (SU13 removed).** The λ₁ lower bound no longer comes from Strohmaier–Uski (SU13). It is now
CERTIFIED in-house by Lehmann–Goerisch, **λ₁ ≥ 3.8388757010** (Lemma 1(d); LB_RESULTS.md; second-read in
review_lb_1: ACCEPT, fixes F1–F3 applied). The certified λ₁ ball is now **[3.8388757010, 3.8388872830]**, and
Lemmas 5 and 8 were re-certified with it (−3.0383476·m(u) and −9.341154560). SU13 (D1) is kept only as an optional
cross-check; the older SU13-based numbers are kept as alternative rows. The weakest inputs are now FBP's quoted
Lemmas 7.2–7.3 and §7.4 geodesic data, and the Lehmann–Goerisch theorem as quoted from secondary sources (D15).

Labels:
- **EXACT**: exact arithmetic, either integers or Q(ζ₂₄) with rational coefficients.
- **CERTIFIED**: Arb ball arithmetic (python-flint 0.9.0) or exact rationals. Every truncation is bounded, and
  floats only *propose* candidates.
- **PROVED**: a written argument.
- **QUOTED**: a published statement we read and quote, not re-verified.

---

## 0. Theorem

**Theorem.** Let B be the Bolza surface and T₂ the Teichmüller space of closed hyperbolic surfaces of genus 2.
Then B is a strict local maximum of λ₁ on T₂: there is a neighbourhood N of B in T₂ such that λ₁(X) < λ₁(B) for
every X ∈ N∖{B}.

More precisely:
- λ₁(B) has multiplicity 3, and λ₁(B) ∈ [3.8388757010, 3.8388872830] (Lemma 1(c),(d); with the optional SU13
  input D1 this narrows to [3.83888625884, 3.8388872830]).
- The first-order cluster map A: T_B → Sym₀(E) has rank 3. Its kernel K is 3-dimensional and is an absolutely
  irreducible Isom(B)-module.
- Off K, λ₁ − λ₁(B) ≤ −c_A|w|/2 to first order (Lemma 10 case (i)). Along K, the cluster sum Λ = λ₁+λ₂+λ₃ satisfies D²Λ|_K = c·g_inv with c < 0.
- Along the equiangular family F_B: (Λ∘x)″(0) ≤ −9.3411 with the SU13-free λ ball (the central value is −14.2742;
  with the optional SU13 ball the certified bound is −12.8504).

**Conditional on** the external inputs of §Dependencies, namely:
- (a) **FBP Lemmas 7.2–7.3** (no 1- or 2-dimensional irreducible representation in the λ₁-eigenspace): D4. Used
  only for multiplicity 3 / irreducibility of E (Lemma 1(c)), hence in Lemmas 5, 6 and 8. Not used by the λ₁
  lower bound.
- (b) **FBP §7.4 geodesic data** (systole, multiplicity lower bounds 12, 12, 24, primitivity): D3. Used in the LP
  of Lemma 1(a), hence in λ₄ > 4.0 and in ρ = 5 of Lemma 1(d).
- (c) **The Selberg trace formula** in FBP's normalisation (FBP-K (3.1), from Buser) and FBP Lemma 3.2 with its
  proof (Lemma A): D2. Used in Lemma 1(a).
- (d) **The Lehmann–Goerisch theorem** (Liu–Plum arXiv:2512.23182 Thm 4.1/4.3 with A1–A5; second source
  Behnke–Mertins–Plum–Wieners 2000 Thm 2.1 as reproduced in Vejchodský arXiv:1801.02561 §2): D15. Used in Lemma
  1(d).
- (e) The other quoted or standard facts D5–D13: FBP §7.1 group data; T_B ≅ Re Q(B) (Tromba via Mazzoli,
  Daskalopoulos–Wentworth); the Fischer–Marsden splitting (via Mazzoli); elliptic regularity, Rellich
  compactness, the spectral theorem, Neumann series, the Riemannian centre of mass; min-max and Ky Fan; Kato's
  one-parameter theory (remark only); the Poincaré polygon theorem; T₂ as a topological 6-manifold, invariance of
  domain and continuity of the mapping-class action; and the trusted software (python-flint/Arb, sympy).
- **Not required:** SU13 (D1). It is an optional cross-check only.

## 1. Setup

**The octagon coordinate space C.**
- x = (V₀, …, V₇) ∈ ℝ¹⁶ are the Klein-model coordinates of the vertices of a convex hyperbolic octagon.
- Side k joins V_k and V_{k+1}. The Bolza gluing pairs side k+4 with side k (k = 0…3) by the orientation-preserving
  isometry g_k ∈ SO⁺(2,1), with V_{k+4} ↦ V_{k+1} and V_{k+5} ↦ V_k.
- C is the set of x satisfying the 5 constraints L_{k+4} = L_k (k = 0…3) and Σθ_k = 2π. The Poincaré polygon
  theorem then gives a closed hyperbolic surface B_x of genus 2, with the single vertex cycle 0→5→2→7→4→1→6→3→0.
- x_B is the regular octagon: sides 2 arccosh(1+√2) = 3.0571418389620, angles π/4.
- By Lemma 9a (certified rank J = 5) and the analytic IFT, C is an analytic 11-manifold near x_B.

**The map π: C → T₂.**
- B_x is marked by the mesh homeomorphism Φ_x: B → B_x, defined below.
- The marking is the isomorphism Γ_{x_B} → Γ_x, g_k(x_B) ↦ g_k(x). It is constant in x up to isotopy.
- π is continuous: it factors through the analytic holonomy map x ↦ (g_k(x)) and the conjugation quotient.

**The mesh map and the forms.**
- The octagon is triangulated by a fixed fan of geodesic triangles. Side nodes sit at exact hyperbolic-arclength
  fractions j/n, matched j ↔ n−j by the pairings. Interior nodes are given by explicit, isometry-equivariant
  formulas: the centre is the normalised vertex sum, and the hyperbolic interpolation is symmetric.
- On a triangle with hyperboloid vertices P₀, P₁, P₂, set X(a,b) = W/√B(W,W), with W = aP₀ + bP₁ + (1−a−b)P₂ and
  B = diag(1,−1,−1). Φ_x is the cellwise map built from this. It is a piecewise-analytic bi-Lipschitz
  homeomorphism, compatible with the gluings.
- The pulled-back forms on H¹(B) are
  - q_x(f) = Σ_T ∫(∇f)ᵀ a_T ∇f, with a_T = adj(LLᵀ − QG)/(D√Q);
  - m_x(f) = Σ_T ∫ ρ_T f², with ρ_T = D·Q^{−3/2};
  - here Q = B(W,W) ≥ 1, L_i = B(W_i, W), G_ij = B(W_i, W_j), W₁ = P₀ − P₂, W₂ = P₁ − P₂ and D = |det[P₀P₁P₂]|
    (CERT_RESULTS §2.2, §4.2).
- The pencil (q_x, m_x) has the Laplace spectrum of B_x.

**The objects of the proof.**
- λ₀ = 0 < λ₁ ≤ λ₂ ≤ … is the Laplace spectrum, and Λ := λ₁ + λ₂ + λ₃.
- E is the λ₁-eigenspace of B, and G = Isom(B), of order 96. G⁺ = Isom⁺(B) ≅ GL(2,3), and J = r⁴ is the
  hyperelliptic involution.
- T_B := T_B T₂ ≅ {transverse-traceless (TT) tensors on (B, h₀)} = {Re φ : φ ∈ Q(B)}. This is a real G-module,
  with G⁺ acting C-linearly by pull-back and reversing elements acting C-antilinearly (Lemma 6, QUOTED
  identification).
- Octagon elements of G: r is the rotation by π/4 (relabelling k → k+1); s_mid is the reflection in the axis
  through two opposite side midpoints; s_vert is the reflection in a diagonal through two opposite vertices.
- t has order 3; u: x ↦ i/x; a is the antipodal isometry, the order-4 orientation-reversing lift of x ↦ −1/x̄.
  These are defined in the model y² = x⁵ − x.
- A(v) ∈ Sym(E) is the first-order cluster matrix (Lemma 2), K := ker A, and D²Λ is the Hessian of Λ at B.

---

## 2. The logical chain

### Lemma 1 (isolation, multiplicity and the λ₁ enclosure). CERTIFIED, citing the trace formula, FBP geodesic data, FBP Lemmas 7.2–7.3 and the Lehmann–Goerisch theorem (SU13 optional)

**Statement.**
- (a) Spec(B) ∩ (0, 3.3] = ∅, and #(Spec(B) ∩ [3.3, 4.0]) ≤ 3, counted with multiplicity.
- (b) λ₃(B) ≤ 3.8391848048 (P1). Hence λ₁, λ₂, λ₃ ∈ (3.3, 3.8391849] and λ₄ > 4.0, a gap ≥ 0.1608. This part uses
  neither SU13 nor FBP's multiplicity theorem.
- (c) λ₁ = λ₂ = λ₃, and dim E = 3. E is a real-irreducible G-module. λ₁ ≤ 3.8388872830 (P2).
- (d) **λ₁ ≥ 3.8388757010** (CERTIFIED, Lehmann–Goerisch; uses only (a) and D15, neither SU13 nor FBP Lemmas
  7.2–7.3). Hence **λ₁ ∈ [3.8388757010, 3.8388872830]**. (Optional, with SU13 = D1: λ₁ ∈ [3.83888625884,
  3.8388872830].)

**Proof of (a).** Use the Selberg trace formula in FBP's normalisation, FBP-K (3.1), from Buser [Bus10, p.253]. The
test functions are f(x) = u(x²/s²)e^{−x²/(2s²)}, where u is an odd-degree (N = 21) Hermite combination with rational
coefficients (CERT_RESULTS §1.2).
- Sign conditions: f ≤ 0 on [sys, ∞) and F ≥ 0 beyond the window. Both are checked by exact Sturm root counts:
  0 roots in each case.
- c is a lower bound for F on the window, from 4000 Arb balls.
- The integral term I uses exact Gaussian moments, plus `acb.integral` on [0, 12] with an explicit tail bound.
- The geodesic term G_P uses lower bounds 12/12/24 on the multiplicities of ℓ₁ = 2 arccosh(1+√2),
  ℓ₂ = 2 arccosh(3+2√2) and ℓ₃ = 2 arccosh(5+3√2).
- Run A uses Lemma A (FBP Lemma 3.2 without the hypothesis a ≤ λ₁; the proof is FBP's, unchanged). Run B uses FBP
  Lemma 3.2 verbatim, with a = 3.3 ≤ λ₁ from run A.

| run | N, s | window | c ≥ | certified bound | target |
|---|---|---|---|---|---|
| A | 21, 1 | [0, 3.3] | 0.9816424759937855 | [0.912612352675290375 ± 9.95e-20] | < 1 ✓ |
| B | 21, 13/10 | [3.3, 4.0] | 0.9964283521820437 | [3.86966111976973759 ± 1.65e-18] | < 4 ✓ |

Script: `lp_count.py` → `lp_count.out`, `lp_A_3.3.json`, `lp_B_3.3_4.0.json`.

Sanity windows give 3.4803 ≥ 3 for [3.3, 3.9] and 7.5466 ≥ 7 for [3.3, 5.4]. X = 4.05 fails (bound 4.0536,
review_cert_1).

**Proof of (b).** The hyperboloid-barycentric P1 space (§1) is exactly H¹-conforming on the glued surface:
- an edge trace depends only on its endpoint values;
- isometries act linearly on W;
- the space contains the constants.

Element integrals ∫_simplex p·Q^{−s} are enclosed by a truncated binomial series with a Lagrange remainder and exact
monomial integration. Min-max (QUOTED: Cook Thm 1.22, from Chavel) on V = span{1, x₁, x₂, x₃} gives λ₃ ≤ μ₃. μ₃ < t is
certified by Sylvester's criterion for tM − K in Arb.
- n = 192: certified_mu3_upper = 3.839184804702046, area/4π = [1.00000000000 ± 2.23e-13].
- Script: `cert_rr.py 48 96 192` → `cert_rr_bolza.out`.

**Proof of (c).**
- FBP Lemmas 7.2–7.3 with j = 1 (QUOTED): the λ₁-eigenspace contains no 1- or 2-dimensional irreducible real
  representation of Isom(B). So m₁ ≥ 3, and with (b), λ₁ = λ₂ = λ₃ and dim E = 3.
- Any decomposition of E would contain a summand of dimension 1 or 2, so E is irreducible.
- The P2 Rayleigh–Ritz bound (Lemma 8, `cert_rr_p2.py`, n = 96: certified_mu3_upper = 3.8388872829261946) gives
  λ₁ ≤ λ₃ ≤ 3.8388872830.
- *Optional cross-check (D1, no longer needed).* SU13 §5 (QUOTED) locates an eigenvalue in [3.83888625884,
  3.83888825884], with midpoint 3.83888725884. By (a) and (b), the only eigenvalue in (0, 4.0] is λ₁, so this one is
  λ₁. It is consistent with (d).

**Proof of (d) (LB_RESULTS.md §0–§1; plan LB_PLAN.md route (c); second-read in review_lb_1: ACCEPT).**
- Lehmann–Goerisch (QUOTED, D15) for D = H¹(B), M = q + m, N = m, i.e. for the shifted eigenvalues λ̃ = λ + 1:
  λ̃₁ = 1 (constants), λ̃₂,₃,₄ = λ₁,₂,₃ + 1, λ̃₅ = λ₄ + 1. By (a), λ₄ > 4.0, so **ρ = 5 ≤ λ̃₅**; ρ > 0 (A5).
  Take m = 4, n = 3. Only (a) is used here; (a) counts with multiplicity and uses neither SU13 nor FBP 7.2–7.3
  (review_lb_1 rows 1e, 4c).
- X = L²(TB) × L²(B), T u = (∇u, u), b_G the L² product (A3). Trial v₁, v₂, v₃: conforming P2 functions on the
  n = 48 mesh (untrusted float eigenvectors, used as exact data). w_i = (σ_i, v_i + div σ_i) with σ_i in an
  RT₂ space in the reference charts; equal-and-opposite edge flux DOFs make σ_i ∈ H(div, B) on the glued surface
  (the flux ι_σdA is intrinsic and the edge parameter s is preserved by the pairings), so
  b_G(w_i, Tφ) = N(v_i, φ) for all φ ∈ H¹(B) exactly (A4), for every σ_i.
- A0, A1, A2 summed over 18432 elements in Arb (element integrals ∫poly·Q^{−s}, K = 14 series terms plus
  Lagrange remainder). Certified in Arb: B = A0 − 2ρA1 + ρ²A2 ≻ 0 and ν̄B − A ≻ 0 with
  ν̄ = −30.031942622361008 < 0. Liu–Plum Thm 4.3 with k = 3: λ₁ + 1 = λ̃₂ ≥ ρ − ρ/(1 − ν̄), i.e.
  **λ₁ ≥ 4 − 5/(1 − ν̄) = [3.8388757010527243542306205221870856407 ± 6.06e-38]**, headline value
  3.8388757010527243 (rounded down; the json field `certified_lambda1_lower` = 3.8388757010527246 is
  rounded to nearest and must not be quoted, review_lb_1 F2). The k = 3 bound alone implies λ₁,₂,₃ ≥ L.
- Script: `lb_bolza.py 48` → `lb_bolza_48.out`, `lb_bolza_48.json` (211 s). n = 24 gives 3.8387028952 (below the
  needed threshold; not used).
- Checks: flat-torus known answer (`lb_torus.out`, N = 32: 39.4783342462 ≤ 4π² = 39.4784176044, O(h⁴));
  exact global Green residual ≤ 7.4e-36; non-conforming sign flip makes B indefinite (negative control); P2
  element matrices agree with `cert_rr_p2` to 2e-24.
- **Independent re-checks (review_lb_1):** the final LG step redone from the stored A0, A1, A2 in mpmath
  interval arithmetic (same 3.83887570105272435423…); element integrals re-derived from a from-scratch metric
  to relative ≤ 1.2e-32; independent geometric gluing check (0 BAD edges); exact toy model showing that wrong
  indexing of the zero mode or an invalid ρ produces false bounds; n = 24 and torus re-runs identical.
- **ρ is not self-checking.** With ρ = 6.5 or 7.5 (> λ̃₅, hypothesis violated) the code still returns sub-true
  values (`lb_bolza_24_rho6.5.out`; review_lb_1 row 4c). So the LP certificate of (a) is load-bearing.
- (d) does **not** give multiplicity: the interval [3.8388757010, 3.8388872830] (width 1.16e-5) does not imply
  λ₁ = λ₂ = λ₃. (c) still needs FBP 7.2–7.3 (review_lb_1 row 6).

**Note (G1 / review_cert_1 row 3; applies to the SU13-based alternative runs only).**
- The SU13-free runs (`lb_cert_sign.py`, `lb_cert_first_order.py`) use instead the ball arb("3.83888149202636215",
  "5.791e-6") = [3.838875701026362, 3.838887283026362], which contains [3.8388757010527243, 3.8388872829261946]
  (exact Fraction check, review_lb_1 row 5a).
- Several scripts use the λ-ball 3.8388872588 ± 1.0000001e-6, i.e. [3.8388862587999, 3.8388882588001]. Its lower
  end lies 4.2e-11 *below* SU13's, which is conservative. Its upper end lies 4.2e-11 below SU13's upper end
  3.83888825884.
- The upper end is nevertheless covered: the certified P2 bound λ₁ ≤ 3.8388872830 lies inside the ball.
- So every certificate that assumes "λ₁ is in the ball" is valid as run. The reviewer re-ran the sign certificate
  with the exact SU13 lower end 3.83888625884 and got −12.850504 (n = 96).

### Lemma 2 (F1: regularity of the cluster). PROVED (F1_PROOF §0–2; second-read: review_theory_1 rows 6a, 6b CONFIRMED)

**Setting.** A pulled-back family is t ↦ (a_t, ρ_t) ∈ L^∞, t ∈ U ⊂ ℝ^d, satisfying:
- (H1) C^k or analytic regularity;
- (H2) ellipticity at 0;
- (H3) (q_t, m_t) is the Dirichlet/L² pair of a closed hyperbolic surface, pulled back by a bi-Lipschitz
  homeomorphism;
- (H4) isolation of the cluster at 0 (Lemma 1: r₀ = 0.16).

The octagon family on C satisfies (H1)–(H3) with k = ω:
- a_T and ρ_T are analytic in the vertex coordinates;
- at x_B, Q ≥ 1 and D > 0 uniformly over the finitely many elements;
- the hyperbolic-arclength interpolation sinh(sL)/sinh L is analytic in cosh L.

The same holds on any analytic sub-family, in particular on the slice S of Lemma 9.

**Statement.** Let Γ = {|z − λ_*| = r₀/2}. Near 0:
1. exactly λ₁(t), λ₂(t), λ₃(t) lie inside Γ;
2. the Riesz projection P(t) = (2πi)⁻¹∮_Γ(zM_t − Q_t)⁻¹M_t dz is analytic;
3. there is an analytic symmetric M(t) ∈ Sym(3) whose eigenvalues are λ₁, λ₂, λ₃, so **Λ = tr M is real-analytic**;
4. M(0) = λ_*I and dM(0)[v] = A(v) := [b_v(e_i, e_j)], with b_v = ∂_v q − λ_*∂_v m and (e_i) m₀-orthonormal in E;
5. M(v) = λ_*I + A(v) + Q₂(v) + R(v), with |R(v)| ≤ C|v|³ uniformly on a ball.

**Proof.** F1_PROOF §2, Steps 1–6. The ingredients are:
- Rellich compactness;
- the spectral theorem for the compact self-adjoint S_t = (Q_t + M_t)⁻¹M_t;
- Neumann-series analyticity of inversion, uniformly on the compact Γ;
- differentiation under the contour integral;
- constant rank of nearby idempotents;
- a second contour Γ′ ⊃ [−½, λ_* + r₀/2], whose projection has rank 4 = {0, λ₁, λ₂, λ₃}. This uses
  Spec ∩ (0, 3.3] = ∅;
- the frame f̃ = P(t)e·G^{−1/2} and M_ij = q_t(f̃_i, f̃_j).

M′ = q′ − λ_*m′ follows by differentiating the m_t-orthonormality. Several parameters are harmless: only the
*individual* eigenvalues can fail to be differentiable, while the trace cannot.

**Literature.** Kato, *Perturbation Theory for Linear Operators* Ch. VII: **not read** by us.
- Cai, arXiv:2608.02459 v2, Prop. 2.5 (QUOTED, read) runs the same Riesz-projection argument for t ∈ ℝ^N. **It is
  cited as an analogue only.** Cai's family is a fixed Dirichlet form with a multiplicative weight, conjugated to a
  fixed L²(µ₀). It is not our variable-coefficient pencil (review_theory_1 row 6e).
- Duran, arXiv:2507.16895 **v3**, §4.1 (QUOTED, second-hand Kato references, one parameter; [31] = Kato 1995).

### Lemma 3 (Lemma N: A is intrinsic and G-equivariant). PROVED

The proof below is review_theory_1 §3, with the fixes of review_lemmaN applied. It was second-read in
review_lemmaN: CORRECT WITH FIXES, and all fixes are incorporated here.

**Statement.** Let h₀ be the Bolza metric. For an L^∞ symmetric 2-tensor ḣ, set
b_ḣ(e, e′) := d/ds|₀ (q − λ_*m)_{h₀+sḣ}(e, e′) on E × E. Then:
- (a) for every v ∈ T_{x}C at x = x_B, A_oct(v) := [b_v(e_i,e_j)] equals A(σ(v)) := [b_{σ(v)}(e_i,e_j)], where
  σ(v) is a TT tensor depending only on v (constructed in Step 4). So A is well defined on TT ≅ T_B.
- (b) A(g·σ) = ρ_E(g)A(σ)ρ_E(g)ᵀ for every g ∈ G, of either orientation.

**Step 1 (first-variation formula).** With d(h^{ij}) = −h^{ik}ḣ_{kl}h^{lj} and d√h = ½(tr ḣ)√h,

  b_ḣ(e, e′) = −∫_B ⟨ḣ, T(e,e′)⟩ dA,  where T(e,e′) := Sym(de⊗de′) − ½(⟨de,de′⟩ − λ_*ee′)h₀.

- The formula is tensorial, so it can be evaluated cellwise. Here h₀ is the Bolza metric on B, and
  h_t := Φ_t*h_{B_t} with Φ₀ = id. (In the reference-simplex coordinates h₀ is only cellwise analytic; this is the
  corrected wording of review_lemmaN Step 1.)
- It agrees with F1's b_v, because F1's coefficients are a = √h h⁻¹ and ρ = √h, and t ↦ h_t is C¹ into L^∞.
- Independent confirmation, QUOTED from El Soufi–Ilias arXiv:math/0701777 Lemma 2.1: "Q_h(u) = −∫_M (du ⊗ du +
  ¼ ∆_g u² g, h) v_g". Since Δ(u²) = 2λu² − 2|du|², their integrand equals T(u,u).

**Step 2 (T is divergence-free on E).** The eigenfunctions e and e′ are smooth by elliptic regularity, and
tr∇²e = −λ_*e. Then

  2 div T = (tr∇²e)de′ + (tr∇²e′)de + ∇²e(∇e′,·) + ∇²e′(∇e,·) − d⟨de,de′⟩ + λ_*(e de′ + e′ de) = 0.

**Step 3 (Lie derivatives drop out).** Let X be a W^{1,∞} vector field. Then L_Xh₀ = 2Sym∇X♭ ∈ L^∞ (Rademacher),
and integrating by parts against the smooth T on the closed surface gives b_{L_Xh₀} = 2∫⟨X, div T⟩ = 0. b is linear
in ḣ ∈ L^∞, so no flow and no H¹-differentiability of e∘ψ_s⁻¹ is needed.

**Step 4 (octagon tangent = TT + Lie derivative; the equivariant smoothing of review_lemmaN Fix 1).**
- Let x(t) be a C¹ curve in C through x_B, h_t := Φ_t*h_{B_t}, and D_t: H² → H² the developing map of h_t. D_t is
  ρ_t-equivariant, cellwise analytic jointly in (t, p), and D₀ = id after normalisation.
- Fix ε > 0 and a smooth radial kernel η_ε(d_{h₀}(p, y)) supported in an ε-ball. Define the Riemannian centre of
  mass

    Ψ̃_t(p) := argmin_{q∈H²} ∫ d(q, D_t(y))² η_ε(d(p,y)) dA(y).

**Properties of Ψ̃_t.**
- *Well defined:* the functional is strictly convex on H² and proper.
- *Equivariant:* Ψ̃_t(γp) = ρ_t(γ)Ψ̃_t(p), because the kernel is deck-invariant and D_t is equivariant.
- *Smooth (the IFT argument):* Ψ̃_t(p) is the unique zero of the gradient of the functional, whose Hessian in q is
  positive definite. By the implicit function theorem it is jointly C^∞ in (t, p). The p-derivatives fall on the
  smooth kernel, and the t-derivatives fall on the cellwise-analytic D_t inside the integral.
- *Ψ̃₀ = id:* D₀ is an isometry and the kernel is symmetric.
- *A diffeomorphism for small |t|:* dΨ̃_t = I + O(t) uniformly, so Ψ̃_t is an equivariant proper local
  diffeomorphism, hence a diffeomorphism.

**The smooth comparison family.** Descend Ψ̃_t to Ψ_t: B → B_t, and set k_t := Ψ_t*h_{B_t} and
φ_t := Ψ_t⁻¹∘Φ_t. Then:
- k_t is a jointly smooth family of smooth hyperbolic metrics with k₀ = h₀;
- φ_t: (B, h_t) → (B, k_t) is an isometry, isotopic to the identity (it commutes with the deck group);
- *Mixed derivatives (review_lemmaN N3):* φ_t is jointly analytic×smooth on each closed cell and continuous across
  cells. The mixed derivative ∂_t∇_pφ_t exists and is continuous on each closed cell, which is what
  differentiating φ_t*k_t requires. Hence X := ∂_tφ_t|₀ is piecewise C¹ and continuous, i.e. Lipschitz, and
  **ḣ = k̇ + L_Xh₀ almost everywhere**.

**The splitting (Fischer–Marsden, via Mazzoli; review_lemmaN N4).** QUOTED from Mazzoli arXiv:2101.09380,
Remark 3.3, citing Fischer–Marsden [FM75]: "Given g a hyperbolic metric, every symmetric tensor S admits an
orthogonal decomposition of the following form: S = Stt + LX g + ((−∆g f + f) g + g∇² f), where: • Stt is transverse
traceless with respect to g; • Stt + LX g is tangent to the space of Riemannian metrics with constant Gaussian
curvature equal to −1. In other words, … Stt + LX g ∈ ker dKg; • (−∆g f + f) g + g∇² f lies in the L²-orthogonal of
ker dKg."
- Apply this to S = k̇ ∈ ker dK_{h₀}. The third component lies in ker dK ∩ (ker dK)^⊥ = 0, so k̇ = σ + L_Yh₀ with
  σ TT.
- Y is smooth by elliptic regularity (σ = Re φ is smooth). The quote does not state this.

**Conclusion.**
- Set σ(v) := σ, the TT part of k̇. By Steps 1 and 3, A_oct(v) = [b_{σ(v)}(e_i,e_j)] = A(σ(v)).
- σ(v) does not depend on ε or on the kernel: two choices of Ψ_t differ by Diff₀, so k̇ changes by a Lie
  derivative.
- v ↦ σ(v) is linear, and it is what "dπ" means throughout.

**Step 5 (equivariance).** Let g be an isometry of h₀, acting on TT by (g⁻¹)* and on E by ρ_E(g)e = e∘g⁻¹. Then
T(e∘g⁻¹, e′∘g⁻¹) = (g⁻¹)*T(e,e′). Since dA is g-invariant as a measure, A(g·σ) = ρ_E(g)A(σ)ρ_E(g)ᵀ.
- This is the action of Lemma 6: (g⁻¹)*Re φ = Re((g⁻¹)*φ) for holomorphic g, and Re(conj((g⁻¹)*φ)) for
  antiholomorphic g.
- The convention (g⁻¹)* versus g* does not affect real characters, because χ(g⁻¹) = χ(g).

**Relabelling line (review_lemmaN Fix 3).** Let γ be an octagon symmetry (r, s_mid, s_vert, …) acting linearly on
C. The fan mesh is γ-equivariant, so Φ_{γx} = Γ∘Φ_x∘γ_B⁻¹. Hence h_{γx(t)} = (γ_B⁻¹)*h_{x(t)} and
**σ(γv) = γ_B·σ(v)**, because the TT projection is natural under isometries of h₀.

**Residual inputs:**
- elliptic regularity;
- the Fischer–Marsden splitting at the single smooth metric h₀ (QUOTED via Mazzoli);
- existence and uniqueness of the Riemannian centre of mass in H² (standard, not quoted).

### Lemma 4 (tr A ≡ 0, i.e. dΛ(x_B) = 0). PROVED, two independent proofs (second-read: review_cert_1 row 3′ CONFIRMED; review_theory_1 row 2c CONFIRMED)

**Proof 1** (CERT_RESULTS §3; no Bochner, no orbifold Teichmüller theory).
- r: ℝ¹⁶ → ℝ¹⁶, (rx)_k = Rot_{π/4}(x_{k−1}), is linear and maps C to C. B_{rx} is isometric to B_x, so Λ∘r = Λ.
- Hence dΛ = dΛ∘P on T := T_{x_B}C, where P = (1/8)Σ_{j<8} r^j.
- (ℝ¹⁶)^r = {δ_k = Rot^k δ₀} is 2-dimensional: radial and tangential. The radial direction is not in T, since
  d(angle sum)/dR = **−5.5729711960** ≠ 0 (CERTIFIED, `cert_rank.py`; equal to the closed form from
  cot(θ/2) = cosh R·tan(π/8)). The tangential direction is the infinitesimal rotation, along which Λ is constant.
- So T^r = span(rotation) and dΛ = 0.

Inputs: Lemma 2 (Λ is C¹ on C); rank J = 5 (Lemma 9a). Independent reviewer check: tr(r | ker J) = −1.000, and
ker[J; r−I] is the rotation field.

**Proof 2** (Schur). tr A is a G-invariant linear functional on T_B (Lemma 3b), and ⟨χ_T, 1⟩ = 0 (EXACT, Lemma 6),
so tr A ≡ 0.

**Use.** tr A ≡ 0 is used twice in Lemma 10: for c_A > 0, and essentially in case (ii).

### Lemma 5 (A ≠ 0). CERTIFIED, citing Lemma 1(c),(d) (Lehmann–Goerisch lower bound, FBP 7.2–7.3) (second-read: review_cert_1 rows 4.1–4.6 CONFIRMED; SU13-free re-run: review_lb_1 rows 5a–5b CONFIRMED)

**Statement.** Let δ be the alternating-radius tangent: vertex radii R + t(−1)^k + σ(t), with σ′(0) = 0. This is
tangent to F_A up to scale. There is an exact eigenfunction e = P_E u ∈ E with

  **b_δ(e, e) ≤ −3.0383476242423058·m(u) < 0** (n = 192, SU13-free λ ball [3.8388757010, 3.8388872830]);
  at n = 96 the bound is −0.9989201873540878·m(u).

| λ₁ ball | n = 96 | n = 192 | script → output |
|---|---|---|---|
| **SU13-free** (Lemma 1(d) + P2) | −0.9989201873540878 | **−3.0383476242423058** | `lb_cert_first_order.py 96 192` → `lb_cert_first_order.out`, `lb_cert_first_order_96_192.json` |
| SU13-free, reviewer re-run, SYM=0 | — | −3.0383476316928864 | review_lb_1, `scratch_RLB/lb_cert_first_order_192.json` |
| SU13 (alternative, D1) | −1.0282409843057394 | −3.0875694043934345 | `SYM=1 cert_first_order.py 96 192` → `cert_first_order.out` |

(The SYM setting of the author's SU13-free run is not recorded in its output; the script default is SYM=0. The
reviewer's SYM=0 run certifies the same bound to 7e-9, so the n = 192 number does not depend on the shortcut;
review_lb_1 F3.)

So A_oct(δ) ≠ 0. By Lemma 3a, A(σ(δ)) ≠ 0, so A ≠ 0 on T_B (and σ(δ) ≠ 0).

**Proof.**
- u = x − (m(x,1)/m(1,1))·1 is exactly m-orthogonal to the constants; x is an untrusted float eigencombination.
  Put w = u − e, which is ⊥ E ⊕ 1 in both m and q.
- With μ = R(u), λ₄ ≥ 4.0 and λ₁ in the λ ball (SU13-free: [3.8388757010, 3.8388872830]), the gap argument gives
  m(w) ≤ (μ − λ₁)/(λ₄ − λ₁)·m(u) and q(w) ≤ (μ − λ₁)λ₄/(λ₄ − λ₁)·m(u).
- b(e,e) = b(u,u) − 2b(u,w) + b(w,w). The cross and quadratic terms are bounded by elementwise relative bounds,
  −β_T a_T ≤ ȧ_T ≤ β_T a_T and |ρ̇_T| ≤ β_{ρ,T}ρ_T, enclosed in Arb over a parameter box ⊇ simplex, with an
  energy-weighted β for the leading term.
- n = 192 values, SU13-free ball (`lb_cert_first_order_96_192.json`):
  - μ = [3.839086294670 ± 5.19e-13];
  - b(u,u)/m = [−4.806193 ± 8.81e-7];
  - m(w)/m = [0.0013 ± 6.49e-5] and q(w)/m = [0.005 ± 2.29e-4];
  - transfer error [1.7 ± 0.0679].
- n = 192 values, SU13 ball (alternative): b(u,u)/m = [−4.806192 ± 1.96e-7]; m(w)/m ≤ 0.0012 and q(w)/m ≤ 0.0049;
  transfer error [1.71 ± 8.63e-3].
- Scripts: `lb_cert_first_order.py 96 192` → `lb_cert_first_order.out` (SU13-free; differs from
  `cert_first_order.py` only in the λ ball (line 290) and the json name); `SYM=1 cert_first_order.py 96 192` →
  `cert_first_order.out` (SU13 alternative).
- The π/2 symmetry shortcut is justified by exact invariance. It was also checked against the full computation at
  n = 48, and by `sym_check_elements.py 16` (55 296 entry pairs, 0 non-overlapping).

**Corrected remark on β_max (review_cert_1 row 4.4′).** The printed β_max ≈ 33 (32.99 at n = 192) is **mostly
enclosure overestimation, not geometry**: the enclosure is taken over the box [0,1]² ⊇ simplex, with the Frobenius
norm and dependency blow-up. At the side-midpoint elements the sampled true sup is ≈ 3.3, against an enclosure of
33.4. Overestimation is safe.

**Robustness (review_cert_1 row 4.6).**
- Robustness in λ₁'s lower bound: evaluated monotonically at the endpoint (exc/(λ₄−λ) is decreasing in λ), the bound still certifies at λ₁ ≥ 3.83788 (≈ −0.0078, thin margin). Under the certificate's plain ball assembly, which treats μ−λ and λ₄−λ as independent, it does NOT certify there (+0.0103; bolzakit v0.2 T8 measurement, confirmed by its reviewer). The SU13-free bound 3.8388757010527… certifies comfortably (b(e,e) ≤ −3.0383·m).
- With λ₄ ≥ 3.9 it still gives −1.89.
- A mesh-only deformation (where A = 0 exactly) is correctly *not* certified (+1.409).

The identification of [b(e_i,e_j)] with the first-order cluster variation along δ needs only the one-parameter
case of Lemma 2.

### Lemma 6 (P6: the decomposition of T_B). EXACT + PROVED (second-read: review_theory_1 rows 1a–2c CONFIRMED, with an independent model)

**(i) The model (EXACT; identification QUOTED).**
- B ≅ {y² = x⁵ − x}. The model rests on QUOTED FBP §7.1: the S₄ quotient of Isom⁺(B), the Weierstrass points, and
  Isom(B) ≅ GL(2,3) ⋊ ℤ/2. It also uses the standard fact that a finite S₄ ⊂ PSL(2,ℂ) is conjugate to the
  octahedral group (not re-read).
- The generators r: (x,y) ↦ (ix, ζ₈y), t: x ↦ (x+i)/(x−i) and κ: (x,y) ↦ (x̄,ȳ) generate a group of order 96, of
  which 48 elements are holomorphic and form GL(2,3) (element orders 1:1, 2:13, 3:8, 4:6, 6:8, 8:12). Since
  |Isom(B)| = 96 (FBP), this is all of G.
- T_B ≅ {Re φ : φ ∈ Q(B)} as a real G-module. QUOTED: Mazzoli arXiv:2101.09380 §1, citing Tromba [Tro92]:
  "Transverse traceless tensors are exactly those 2-tensors that can be written as Re φ"; also
  Daskalopoulos–Wentworth arXiv:math/0604403 §2.
- Q(B) = Sym²H⁰(K) has basis (dx/y)², x(dx/y)², x²(dx/y)².
- Script: `cert_chars.py` part A → `cert_chars_A.out` / `.json` (EXACT).

**(ii) Characters on T_B (EXACT).**
- χ_T(r) = χ_T(r³) = −2, χ_T(r²) = −2, χ_T(J) = 6 and χ_T(t) = 0.
- χ_T(g) = 0 for every orientation-reversing g, because an antilinear map has real trace 0.
- Hence ⟨χ_T, χ_T⟩ = 2 and ⟨χ_T, 1⟩ = 0 (so T_B^G = 0).
- These values were reproduced by an independent pull-back model in review_theory_1 (`scratch_RT/indep_chars.out`).

**(iii) The antipodal splitting (review_theory_1 §2, the simpler argument).**
- a is the order-4 orientation-reversing lift of x ↦ −1/x̄. Its conjugacy class {a, a³ = aJ} has size 2 (Cook
  class 4a), so a is central modulo J, and a² = J.
- J acts trivially on T_B (χ_T(J) = 6 = dim), so a² = id on T_B. Hence T_B = T₊ ⊕ T₋, the ±1-eigenspaces of a, and
  both are G-invariant.
- Both have rank 3 (EXACT). ⟨χ₊,χ₊⟩ = ⟨χ₋,χ₋⟩ = 1 and ⟨χ₊,χ₋⟩ = 0 (EXACT).
- Define **K₀ := T₋** and **W₀ := T₊**. They are absolutely irreducible and non-isomorphic. In Cook's labelling,
  K₀ ≅ χ₇ and W₀ ≅ χ₁₀ (matched on all 13 classes against Coo18 App. A, QUOTED; the labels are not used
  logically).

**(iv) K₀ ⊂ ker A (PROVED; the antipodal argument is valid here because E is real-irreducible of odd dimension,
Lemma 1c).**
- E has odd real dimension, so End_G(E) = ℝ, and the central J acts on E by ±1.
- If J = −1, then ρ_E(a)² = −I on ℝ³, which is impossible (det² = −1). So J = 1 on E.
- Then ρ_E(a) commutes with ρ_E(G) (a is central mod J), so ρ_E(a) = ±I.
- By Lemma 3b, A(a·v) = ρ_E(a)A(v)ρ_E(a)⁻¹ = A(v). For v ∈ T₋ we have a·v = −v, so A(v) = −A(v) = 0.
- This needs no character table and no identification of E.

**(v) ker A = K₀ exactly, and rank A = 3 (PROVED).**
- A|W₀: W₀ → Sym(E) is G-equivariant (Lemma 3b) and W₀ is irreducible, so A|W₀ is 0 or injective.
- A ≠ 0 (Lemma 5) and A|K₀ = 0 force A|W₀ ≠ 0, hence injective.
- So **K = ker A = K₀** and rank A = 3.

Here the antipodal argument alone is *not* enough. It needs the exact fact that W₀ is irreducible (norm 1). This
is the only place where (iii)'s character computation is used.

**(vi) The Hessian on K is a scalar, given G-invariance (PROVED; see OPEN-1 for the G-invariance).** K₀ is
absolutely irreducible, so every G-invariant symmetric bilinear form on K₀ is c·g_inv for a fixed invariant inner
product g_inv.

**Consistency (not used).**
- The full Clifford/Gallagher computation (CERT_P6 §2.3) shows Sym₀²E ≅ D ⊕ W₀ for each of the four candidate
  3-dimensional irreducibles. This gives the same conclusions by Schur.
- The certified character of E is tr(g|E) ∈ [1 ± 0.0352] for g = r, s_mid, s_vert (n = 96;
  `cert_chars.py` part B → `cert_chars_B.out`). So E ≅ χ₈, the octahedral representation, in agreement with Cook
  Thm 3.9. This identification is **not used** by the main argument.
- Corrected wording (review_theory_1 row 1d). s_mid fixes 4 Weierstrass points: the centre, the vertex point, the
  midpoint class M₀ **on its axis**, and the class M₂, which is **fixed via the gluing** (M₂ ↦ M₆ ~ M₂) and does
  not lie on the axis. s_vert fixes 2: the centre and the vertex point.

### Lemma 7 (F_B ⊂ K₀). EXACT + PROVED (second-read: review_theory_1 rows 3a–3c CONFIRMED)

**The family.**
- F_B(p): all vertices at radius R(p); V_{2j} at angle jπ/2 and V_{2j+1} at φ + jπ/2, with p = φ − π/4.
- The surface condition, CERTIFIED as Arb jets (CERT_SIGN §1.5), is C²xy − C(x+y) − 1 = 0, with C = cosh R,
  x = tan(φ/2) and y = tan(π/4 − φ/2). At p = 0, C = 3 + 2√2.
- F_B(p) is invariant under r² and under s_mid (the bisector reflection of side 0), and r maps F_B(p) to F_B(−p).

**The line.**
- The exact null space {v ∈ T_B : r·v = −v, s_mid·v = v} is 1-dimensional (EXACT, `cert_chars_A.json`:
  dim_T_B_eps_s_line = 1).
- Its spanning vector satisfies a·v = −v, so **v ∈ T₋ = K₀**. This is EXACT, from cert_chars (v_B_line_lies_in X₁)
  and independently from review_theory_1 (`indep_chars.out`: "A v == -v … True").

**w_B lies on the line.** By the relabelling line of Lemma 3, w_B := σ(dF_B/dp) satisfies r·w_B = −w_B and
s_mid·w_B = w_B. So w_B ∈ ℝv. **The reflection type matters (LEAN.md):** only s_mid-type reflections (side-bisector) give a line inside K₀; the analogous line for an s_vert-type (vertex) reflection lies in W₀ (`Rep.smid_line_in_Kminus`).

**w_B ≠ 0** (PROVED by hyperbolic trigonometry; review_theory_1 re-derived it).
- The pairing g₀ of F_B(p) is the translation along the bisector of side 0 by 2h_a.
- dℓ(g₀)/dp = −cos(π/4)/(sin²(π/8) sinh r_in) = −√(2+2√2) = **−2.19736822693562**, and dℓ(g₁)/dp = +2.19736822693562
  (mpmath, `scratch_RT/fb_trig.out`).
- ℓ(g₀) − ℓ(g₁) is a function on T₂ with nonzero derivative along π∘F_B.
- Independently, the certificate of Lemma 8 itself forces dπ(v_B) ≠ 0 (review_sign_1 row 5b).

**Conclusion.** **w_B ∈ K₀ ∖ {0}.**

### Lemma 8 (the Ky Fan sign bound along F_B). CERTIFIED, citing Lemma 1 (λ₄ > 4.0, FBP 7.2–7.3, Lehmann–Goerisch lower bound of 1(d)) (second-read: review_sign_1, ACCEPT; SU13-free re-run: review_lb_1 rows 5a–5b CONFIRMED)

**Statement.** Let x(p) be F_B, with v_B = x′(0). Then

  (Λ∘x)″(0) ≤ g″(0) ≤ **−9.341154560446736** (n = 96, P2, SYM=1, SU13-free λ ball [3.8388757010, 3.8388872830]),

with central value g″(0) = [−14.274 ± 2.13e-4] (SU13-free run; [−14.2742 ± 2.23e-5] in the SU13 run).

| λ₁ ball | certified g″(0) ≤ (n = 96) | script → output |
|---|---|---|
| **SU13-free** (Lemma 1(d) + P2) | **−9.341154560446736** | `lb_cert_sign.py 96` → `lb_cert_sign.out`, `lb_cert_sign_96.json` |
| SU13-free, reviewer re-run (SYM=1) | −9.341154566965994 | review_lb_1, `scratch_RLB/lb_cert_sign_96.json` |
| SU13 (alternative, D1) | −12.850476214196528 | `cert_sign.py 48 96` → `cert_sign.out`, `cert_sign_48_96.json` |

The rest of this lemma (symmetry-shortcut checks, error budget, sensitivity) describes the SU13-based run unless
stated otherwise; the SU13-free run uses the same script with only the λ ball changed (`diff cert_sign.py
lb_cert_sign.py`: line 34 and the json name).

The symmetry shortcut is not needed:
- the reviewer's SYM=0 re-run at n = 96 gives −12.850476213730866;
- the author's SYM=0 runs certify at n = 24 (`cert_sign_nosym.out`, from the run `12 24`; n = 12 alone does not
  certify) and at n = 48, with −12.566084126010532 (`cert_sign_48_nosym.out`).

Other certified values (CERT_SIGN §0): −9.7138 at n = 24 and −12.5661 at n = 48.

**Proof.**
1. **Ky Fan (PROVED from min-max; F1_PROOF §5; review_theory_1 row 6d CONFIRMED).**
   - λ₀ ≡ 0, with the constants as eigenfunctions. For any 4-dimensional V ⊂ H¹, Λ(p) ≤ tr_V(q_p | m_p).
   - Take V(p) = span{1, e_i + pψ_i}, where e_i ∈ E are exact and the ψ_i are fixed P2 correctors (untrusted floats).
     Set g(p) = tr(G_m⁻¹G_q). Then g ≥ Λ∘x and g(0) = Λ(0).
   - g is analytic, and Λ∘x is analytic by Kato's one-parameter theory (type (a) holomorphic family, isolated
     cluster with gap ≥ 0.16). So h = g − Λ∘x ≥ 0 has a minimum at 0, and **(Λ∘x)″(0) ≤ g″(0)**.
   - **F1 (Lemma 2 on C) is not needed for this inequality itself (review_sign_1 G2).** Even without twice
     differentiability, the limsup of the second difference quotients of Λ∘x is ≤ g″(0). F1 and Lemma 4 are
     needed only to read the result as D²Λ(v_B, v_B) < 0.
   - Finite-dimensional Ky Fan is QUOTED from Moslehian arXiv:1108.1467, attributed to K. Fan, PNAS 36 (1950).
2. **Formula (re-derived by review_sign_1; exact pencil test agrees to relative 3e-39).**

     g″(0) = tr(N X) − 2 tr(N M₁ N A) + 2λ c₁ᵀ N c₁ / m(1,1),  N = M₀⁻¹,

   where:
   - X_ij = (q″−λm″)(e_i,e_j) + 2b(ψ_i,e_j) + 2b(e_i,ψ_j) + 2(q−λm)(ψ_i,ψ_j);
   - A_ij = b(e_i,e_j);
   - M₁,ij = m′(e_i,e_j) + m(ψ_i,e_j) + m(e_i,ψ_j);
   - c₁,i = m′(e_i,1) + m(ψ_i,1);
   - M₀,ij = m(e_i,e_j).
3. **P2 hyperboloid-barycentric elements (`cert_rr_p2.py`).**
   - The edge trace is a quadratic in the edge parameter s, fixed by the endpoints and the geodesic midpoint.
     Isometries preserve s, so the space is exactly conforming.
   - Element integrals ∫poly·Q^{−s} for s ∈ {½, 3/2, 5/2, 7/2} use K = 6 series terms plus a Lagrange remainder.
   - Order-2 Taylor jets in p over Arb.
   - Validation at Bolza (n = 96): certified λ₃ ≤ 3.8388872829261946; certified cluster sum ≤ 11.516661838577997;
     area/4π = [1 ± 1.19e-18].
4. **Transfer to the exact e_i (gap argument, as in Lemma 5).**
   - At n = 96, SU13-free ball (`lb_cert_sign_96.json`): √q(w_i) ≤ 0.016953498968304958, 0.016953498968304958 and
     0.016957233113081952; √m(w_i) ≤ 0.008476749477446675, 0.008476749477446675 and 0.00847861654983665.
   - At n = 96, SU13 ball (alternative): √q(w_i) ≤ 5.0299e-3, 5.0299e-3 and 5.0425e-3; √m(w_i) ≤ 2.515e-3,
     2.515e-3 and 2.521e-3.
   - μ_i − λ_low is almost entirely the width of the λ interval (discretisation part ≈ 2e-8): ≈ 1.02e-6 with
     SU13; with the SU13-free ball (width 1.16e-5) the distances are correspondingly larger, so **η is at the
     floor set by the λ₁ lower bound**. This is why the certified constant drops from −12.85 to −9.34.
   - Every entry of X, A, M₁, c₁, M₀ is its value at u plus bounded cross terms in w. Relative bounds
     −β¹a ≤ ȧ ≤ β¹a and −β²a ≤ ä ≤ β²a are enclosed over 3 boxes covering the simplex.
   - The 3×3 algebra is done in Arb (arb_mat inverse).
5. **Error budget at n = 96** (`cert_sign.out`, `cert_sign_48_96.json`).

   | term | value |
   |---|---|
   | T₂ = tr(q″−λm″)(u,u), float | 45.635701833848785 |
   | coupling (ψ terms), float | −59.9099 |
   | central g″(0) | [−14.2742 ± 2.23e-5] |
   | X1 (4·b(ψ,w) cross terms) | +0.839368 |
   | X2 (2·(q″−λm″)(u,w) and w–w terms) | +0.570060 |
   | M0 | +4.6e-4 |
   | A, M1, c1 | 1.5e-5, 1.5e-5, 7.4e-5 |
   | **certified upper bound** | **−12.850476214196528** |

   **The breakdown rows are "one-at-a-time" (review_sign_1 G3):** each is the upper end minus the central value
   with only that error group switched on. They do not add up exactly to the joint total (≈ 1.41 against the
   joint 1.42). The certified number is the joint Arb ball. Adding the rows to the central value would give ≈ −12.864, which is *below* the certified −12.8504762; so the rows must not be summed into a bound (checked in Lean: `Cert.sign_rows_not_additive`).

   **On the β values (review_cert_1 row 4.4′, applied here too).** The P2 maxima β¹_max = 13.2 and β²_max = 38.1
   are enclosure values. They should not be read as the true pointwise sup. review_sign_1 sampled the true ratios
   on ~150 elements (n = 12, 48) and found observed/enclosure ≤ 0.70 (β¹), 0.89 (β²), 0.99 (β¹ρ) and 0.97 (β²ρ).
   The leading terms use the energy-weighted rms values: 4.2/9.1 for u, and 4.6/4.0 for ψ.
6. **Sensitivity (review_sign_1 §2, `scratch_RS/bisect96.out`).** The certificate holds for every rigorous λ₁
   lower bound **≥ 3.8388029653**, and fails at 3.8388029297. That threshold is ≈ 8.4e-5 below λ₁, so SU13's claimed ±1e-6
   has ≈ 83× slack (exactly (3.83888625884 − 3.8388029653)/1e-6 ≈ 83.3, LEAN.md). If SU13's 1e-6 were relative (±3.84e-6), the bound would still certify (n = 48, radius 1e-5,
   gives −9.43). The SU13-free bound 3.8388757010 (Lemma 1(d)) clears the threshold by 7.3e-5 (LB_RESULTS §0).
7. **Negative controls (review_sign_1 §3), all of which fail as they should:** ψ = 0; λ-radius 1e-4; λ₄ = 3.84;
   a λ-ball excluding the true λ (raises).
8. **Independent float check (NUMERICAL).**
   - review_sign_1's own conforming P1 code, on a different mesh and with re-meshed finite differences, gives
     Λ″ ≈ −14.26 by Richardson at n = 64/128.
   - The 5-point finite difference of the P2 cluster gives −14.2739 (n = 24).

**Consequence.** By Lemma 4, dΛ_C(x_B) = 0, so (Λ∘x)″(0) = D²Λ_C(v_B, v_B) + dΛ_C(x″(0)) = D²Λ_C(v_B, v_B). Hence
**D²Λ_C(v_B, v_B) ≤ −9.3411 < 0** (SU13-free; −12.8504 with the optional SU13 ball).

### Lemma 9 (G3: local chart of T₂ at B). CERTIFIED + PROVED (CERT_G3; second-read: review_proof_1 §2, re-run reproduced)

**(a) rank J = 5** (CERTIFIED, `cert_rank.py` → `cert_rank.out`). The minor in columns (0,2,5,7,8) has
det = [46151927.03406868578682884302258 ± 2.89e-24]. Re-run by review_cert_1.

**(b) rank dπ(x_B) = 6** (CERTIFIED, `cert_dpi_rank.py` → `cert_dpi_rank.out`).
- The six SO(2,1)-trace functions τ_w are those of the words g₀g₁, g₀G₂, g₁g₂, g₁g₃, g₁G₃ and g₂g₃ (G_k = g_k⁻¹).
  At x_B their values are 340.705627484771 / 134.882250993909.
- Restricted to T_{x_B}C = ker J (an Arb basis), they give a 6×11 Jacobian JT. 392 of its 462 6×6 minors are
  certified nonzero.
- The best minor has tangent columns (1,3,5,7,8,10), i.e. free coordinates 3, 6, 10, 12, 13, 15, with
  det = **[7872428136359854080.0 ± 3e-6]**.
- Cross-check: the 11×11 minor of [Jc; Jf] has det = [−1.3827971011809151404e+24 ± 2.74e+4].
- Float cross-check (`scratch_G3/xcheck.out`): pairings agree to 2.4e-13 and the Jacobian to relative 3.3e-9.

**(c) Topology (PROVED).**
- τ_w = τ̄_w∘π, because tr_{SO(2,1)} is conjugation-invariant.
- S ⊂ C is the analytic 6-dimensional coordinate slice through x_B in the free coordinates {3,6,10,12,13,15}.
  Φ := (τ_w)_w is an analytic diffeomorphism from a neighbourhood of x_B in S onto an open U ⊂ ℝ⁶.
- π|_S is continuous and injective, so by invariance of domain it is a homeomorphism onto an open N ∋ B. T₂ is a
  topological 6-manifold (QUOTED, Felikson–Natanzon arXiv:1102.4861, Fenchel–Nielsen).
- So the τ̄_w are local coordinates of T₂ at B.

**(d) Consequences (PROVED).**
1. p := (Φ|_S)⁻¹∘Φ is an analytic retraction C → S with π∘p = π. Hence Λ_C = Λ_S∘p,
   D²Λ_C(w,w) = D²Λ_S(dp w, dp w) and dΛ_S(x_B) = 0.
2. For a curve x(t) in C, with k_t from Lemma 3 Step 4: d(τ_w∘x)(v) = dτ_w(σ(v)), where dτ_w is the derivative
   of the smooth Diff₀-invariant length function on smooth metrics.
   - rank JT = 6 then forces σ to be onto the 6-dimensional TT space, and an isomorphism on T_{x_B}S.
   - **(dτ_w)_w is injective on TT.**
   - (Smoothness of τ_w: τ_w(k) is the trace of the holonomy of the flat so(2,1) Cartan connection ω_LC + θ of k, which is defined and C¹ in k for every metric along a fixed loop, and is flat, hence free-homotopy-invariant, exactly when K = −1. The identity holds for every word w, not only the six.)
3. Consequently **σ∘dp = σ** on T_{x_B}C. This follows by applying (dτ_w)_w (injective) to
   dτ(σ(dp v)) = d(τ∘p)(v) = dτ(v) = dτ(σ(v)), since τ∘p = τ.
   - In particular σ(dp v_B) = w_B ∈ K₀ and dp v_B ∈ K_S := σ⁻¹(K₀) = ker(A∘σ|_{T S}).
   - (This one-line deduction is made explicit here; it uses only (d)2.)
4. Lemma 2 applies on S, an analytic sub-family of C. So in the chart S: M_S(v) = λ_*I + A(σv) + Q₂(v) + O(|v|³),
   and tr A ≡ 0.

### Lemma 10 (local-max criterion). PROVED (RESULTS §2.3; second-read: review_1 §1.3 "valid")

**Statement.** Let M(v) be a symmetric 3×3 matrix function on a neighbourhood of 0 in ℝ⁶ with
M(v) = λ_*I + A(v) + Q₂(v) + R(v), where A is linear, Q₂ is quadratic and |R| = o(|v|²) uniformly. Assume:
- tr A ≡ 0;
- tr Q₂(k) = ½D²Λ(k,k) ≤ −κ|k|² for all k ∈ K = ker A, with κ > 0.

Then λ_min M(v) < λ_* for all small v ≠ 0. The criterion is sufficient, not necessary; the sharper condition is
λ_min Q₂(k) < 0 on K.

**Proof.**
1. Let W be a complement of K, so A|_W is injective. Set c_A = min_{|w|=1}(−λ_min A(w)). Each A(w) is traceless
   and nonzero, so it has a negative eigenvalue; the minimum is over a compact set, so **c_A > 0 (first use of
   tr A ≡ 0)**.
2. Let C′ bound ‖Q₂(v)‖ ≤ C′|v|², and |R(v)| ≤ C′|v|² with R = o(|v|²), on a fixed ball. Write v = w + k.
3. **Case (i): |w| ≥ C|k|².** By Weyl's inequality,
   λ_min M(v) − λ_* ≤ −c_A|w| + 2C′(|w|² + |k|²) + o.
   With C > 8C′/c_A this is ≤ −c_A|w|/2 + O(|w|²) < 0.
4. **Case (ii): |w| < C|k|².** Then |v| ≍ |k|, and
   tr(M(v) − λ_*I) = **tr A(w)** + tr Q₂(k) + O(|k|³) + o(|k|²) ≤ −κ|k|² + o(|k|²) < 0,
   **using tr A(w) = 0 (second, essential use)**. So λ_min ≤ tr/3 + λ_* < λ_*. ∎

### Lemma 11 (G-invariance of the Hessian on K). PROVED (§OPEN-1 closure; second-read: review_proof_1 §1)

**Statement.** The quadratic form H(σ) := D²Λ_S(σ⁻¹σ, σ⁻¹σ) on TT (with σ: T_{x_B}S → TT the isomorphism of
Lemma 9d) is invariant under the natural action of all of G. With Lemma 6(vi): H|_{K₀} = c·g_inv.

### Theorem: assembly

1. Work in the chart S (Lemma 9). Λ_S is real-analytic (Lemma 2 on S), M_S has the expansion of Lemma 2 with first
   term A∘σ (Lemma 3a), and tr(A∘σ) ≡ 0 (Lemma 4).
2. ker(A∘σ) = σ⁻¹(K₀) =: K_S is 3-dimensional (Lemmas 5 and 6v; σ is an isomorphism, Lemma 9d).
3. dp v_B ∈ K_S ∖ {0} (Lemmas 7 and 9d3), and D²Λ_S(dp v_B, dp v_B) = D²Λ_C(v_B, v_B) ≤ −9.3411 (Lemmas 8 and 9d1; −12.8504 with the optional SU13 ball).
4. By Lemma 11, D²Λ_S|_{K_S} = c·(σ*g_inv). By step 3, c·g_inv(w_B, w_B) ≤ −9.3411, so c < 0 and
   D²Λ_S|_{K_S} is negative definite. Hence κ > 0 in Lemma 10.
5. Lemma 10 gives λ₁(π(x)) = λ_min M_S(x) < λ_* = λ₁(B) for x ∈ S near x_B, x ≠ x_B.
6. π|_S is a homeomorphism onto a neighbourhood N of B (Lemma 9c). So λ₁ < λ₁(B) on N ∖ {B}. ∎

---

## OPEN

### OPEN-1 (CLOSED): G-equivariance of the chart for the non-octagon isometries (blocks the step "one direction ⇒ all of K")

**The gap.** Assembly step 4 needs D²Λ|_K to be invariant under a group acting absolutely irreducibly on K₀:
- CERT_P6 row 6 and §2.5 assert this, and review_theory_1 row 3d confirms it, by "Λ is G-invariant, so D²Λ is
  G-invariant".
- But the only structure in which Λ has been shown to be C² is the chart S of Lemma 9, and the only map identifying
  T_{x_B}S with the G-module TT is σ.
- σ is shown to intertwine the actions **only for the octagon relabellings** (Lemma 3, relabelling line =
  review_lemmaN Fix 3).

**Why the octagon group is not enough.** By the exact table (CERT_P6 §3), K₀ ≅ χ₇ restricted to the octagon D₈
(16 elements: e, J, r^{±1}, r^{±3}, r^{±2}, 4 × s_mid-type, 4 × s_vert-type) has norm
(9 + 9 + 4·1 + 2·1 + 4·1 + 4·1)/16 = 2. So K₀|_{D₈} is reducible: ε_s ⊕ (a 2-dimensional irreducible), with ε_s of
multiplicity 1 (`cert_chars_A.json`: D8_mult_eps_s in X1 = 1).
- D₈-invariant quadratic forms on K₀ therefore form a 2-parameter family.
- The certified direction w_B spans exactly the ε_s line (Lemma 7). So D₈-invariance plus Lemma 8 controls only
  that line, not the 2-dimensional complement in K₀.
- An element outside D₈ (t, u or a) is needed, and its compatibility with the chart is proved in no reviewed
  document. review_lemmaN N1 flagged the closely related smooth-structure issue for k_t; Fix 1 resolved that, but
  not this one.

**Proposed closure (written here for the first time; second-read in review_proof_1: CORRECT).** For g ∈ G, let g act on T₂ by change of
marking, normalised so that g·[k] = [(g⁻¹)*k] for smooth hyperbolic metrics k on B. Then g fixes B and
Λ∘g = Λ (the surfaces are isometric).
1. *g is analytic in the chart S.*
   - For every word w′ in the g_k: τ̄_{w′}(g·X) = τ̄_{g⁻¹(w′)}(X), because (g⁻¹)*k-lengths of the class w′ are the
     k-lengths of the class g⁻¹(w′).
   - On S, τ̄_{w″}∘π|_S = τ_{w″}|_S is analytic for every word w″, and the chart is (Φ|_S)⁻¹∘τ̄.
   - So g_S := (π|_S)⁻¹∘g∘π|_S is analytic near x_B and fixes x_B.
   - Since dΛ_S(x_B) = 0, the Hessian is invariant: D²Λ_S(dg_S v, dg_S v) = D²Λ_S(v, v).
2. *σ intertwines:* σ(dg_S v) = g·σ(v).
   - Let x(t) be a curve in S and y(t) := g_S(x(t)). By Lemma 9d2 applied to y and to x:
     dτ_w(σ(y′)) = d/dt τ_w(y(t)) = d/dt τ_{g⁻¹w}(x(t)) = dτ_{g⁻¹w}(σ(x′)).
   - On smooth metrics, τ_{g⁻¹w}(k) = τ_w((g⁻¹)*k). Differentiating at the g-fixed h₀ gives
     dτ_{g⁻¹w}(σ) = dτ_w(g·σ).
   - Hence dτ_w(σ(y′)) = dτ_w(g·σ(x′)) for the six words w. Since (dτ_w)_w is injective on TT (Lemma 9d2),
     σ(y′) = g·σ(x′).
3. Steps 1 and 2 give Lemma 11 for all g ∈ G, in particular for t. With Lemma 6(vi), D²Λ_S|_{K_S} = c·σ*g_inv.

**Inputs of the proposed closure:**
- Lemma 9d2. This is the identity d(τ_w∘x)(v) = dτ_w(σ(v)) and the smoothness of τ_w on smooth metrics, both from
  CERT_G3, second-read in review_proof_1.
- The (trivial) naturality of geodesic length under isometries.
- The definition of the G-action on T₂.
- continuity of the extended mapping-class action on T₂ (standard), used so that g·π(x) ∈ N for x near x_B.

**Status:** this closure and CERT_G3 §3 were accepted in review_proof_1, so OPEN-1 is closed and the Theorem holds as stated, conditional on §Dependencies.

**Alternative closure:** cite that π: C → T₂ is real-analytic into the standard analytic structure of T₂, and that
the stabiliser of B in the extended mapping-class group acts analytically with derivative equal to the natural
action on Re Q(B) (Earle–Eells / Tromba / Weil–Goldman). Not read by us.

No other logical hole was found. The remaining items are **status** items, not holes: FBP Lemmas 7.2–7.3 and §7.4 data are quoted, the Lehmann–Goerisch theorem is quoted from secondary sources (Goerisch's papers not read), and Tromba's book and Kato were not read. (SU13's weak provenance no longer matters: D1 is optional since 2026-10-09.) They are listed in §Dependencies
and §Reviews.

---

## Dependencies (complete list of external inputs)

| # | Input | Used in | Status |
|---|---|---|---|
| D1 | *(OPTIONAL since 2026-10-09; not required)* SU13 (Strohmaier–Uski, CMP 317 (2013), arXiv:1110.2150) §5: an eigenvalue of B lies within 1e-6 of 3.83888725884 | Optional cross-check of Lemma 1(d), and the alternative (SU13-ball) rows of Lemmas 5 and 8 | **No longer required.** Replaced by the certified Lehmann–Goerisch bound λ₁ ≥ 3.8388757010 (Lemma 1(d), D15). Kept as a consistency check: SU13's interval [3.83888625884, 3.83888825884] is consistent with (overlaps) the certified [3.8388757010, 3.8388872830]. Provenance remains weak (unpublished code, one sentence: "The rigorous implementation that uses very simple bounds in its current form gives a rigorous error bound of 10−6", absolute/relative unstated) |
| D2 | Selberg trace formula in FBP's normalisation (FBP-K arXiv:2111.14699 (3.1), from Buser [Bus10, p.253]); FBP Lemma 3.2 and its proof (Lemma A) | Lemma 1a, hence λ₄ > 4.0 and ρ = 5 in Lemma 1d | QUOTED; checked verbatim by review_cert_1 |
| D3 | FBP §7.4 geodesic data: sys(B) = 2 arccosh(1+√2) (Jenni, Schmutz); unoriented multiplicities ≥ 12, 12, 24 for ℓ₁, ℓ₂, ℓ₃; primitivity since ℓ₃ < 2ℓ₁ | Lemma 1a, hence ρ = 5 in Lemma 1d | QUOTED; only lower bounds are used. **One of the weakest inputs** |
| D4 | FBP Lemmas 7.2 and 7.3 (no 1- or 2-dimensional irreducible representation in the λ₁-eigenspace) | Lemma 1c (multiplicity 3, E irreducible), hence Lemmas 5, 6, 8 (dim E = 3 in the gap transfer). **Not** used by Lemma 1d | QUOTED; their proofs use symmetry, nodal domains and FBP Prop. 4.4, with no eigenvalue input. **One of the weakest inputs** |
| D5 | FBP §7.1: Isom⁺(B)/⟨r⁴⟩ ≅ S₄, Isom(B) ≅ GL(2,3) ⋊ ℤ/2, \|Isom(B)\| = 96, the hyperelliptic involution is r⁴; plus the standard fact that a finite S₄ ⊂ PSL(2,ℂ) is octahedral | Lemma 6(i) (the model y² = x⁵ − x) | QUOTED / standard (not re-read) |
| D6 | T_B ≅ {Re φ : φ ∈ Q(B)} (Tromba [Tro92], via Mazzoli arXiv:2101.09380 §1; Daskalopoulos–Wentworth arXiv:math/0604403 §2) | Lemma 6(i), Lemma 3 | QUOTED; Tromba's book not read |
| D7 | Fischer–Marsden splitting (Mazzoli Remark 3.3, citing [FM75]) | Lemma 3 Step 4 | QUOTED; FM75 not read |
| D8 | Elliptic regularity; Rellich compactness; spectral theorem for compact self-adjoint operators; Neumann series; Riemannian centre of mass in H² | Lemmas 2, 3 | standard |
| D9 | Min-max (Cook 2018 Thm 1.22, from Chavel) and Ky Fan (Moslehian arXiv:1108.1467, attributed to Fan 1950; reduced to min-max in F1_PROOF §5) | Lemmas 1b, 8 | QUOTED + PROVED reduction |
| D10 | Kato one-parameter analytic perturbation (type (a)), only for the remark that (Λ∘x)″ ≤ g″ needs no F1; Kato not read (Duran v3 quotes it second-hand) | Lemma 8 remark | standard; not load-bearing |
| D11 | Poincaré polygon theorem (the constraints give a closed surface with the Bolza gluing) | Setup | standard |
| D12 | T₂ is a topological 6-manifold (Fenchel–Nielsen; Felikson–Natanzon arXiv:1102.4861 §2, QUOTED); invariance of domain; continuity of the quotient to T₂; continuity of the extended mapping-class group action on T₂ | Lemma 9c | QUOTED / standard |
| D13 | Software: python-flint 0.9.0 (Arb ball arithmetic, `acb.integral`), sympy 1.14.0 exact Sturm counts and exact Q(ζ₂₄) arithmetic. numpy/scipy are used only for untrusted proposals | all CERTIFIED steps | trusted |
| D15 | **Lehmann–Goerisch theorem**: Liu–Plum arXiv:2512.23182 Thm 4.1/4.3 with assumptions A1–A5 (§4.1; A5 includes ρ > 0), from Goerisch–Albrecht / Goerisch–Haunhorst; second source Behnke–Mertins–Plum–Wieners (2000) Thm 2.1, in shifted form, as reproduced in Vejchodský arXiv:1801.02561 §2 | Lemma 1d, hence Lemmas 5 and 8 | QUOTED from **secondary sources**; statement checked word for word against the Liu–Plum PDF and matched to BMPW's shifted form by review_lb_1 (rows 1a, 1b). Goerisch's original papers and the NPW book chapter were not read. **One of the weakest inputs** |
| D14 | (OPEN-1 alternative only) analyticity of π into T₂ and of the mapping-class action | Lemma 11 | not needed (OPEN-1 closed in the chart S) |

Not used, although recorded elsewhere: SU13 (D1; optional cross-check only since 2026-10-09); FBP Thm 7.1 (m₁ = 3; our Lemma 1 rederives it from D2–D4); Cook's character
table (only for labels); Bochner linearisation (F2) and (2,3,8)-orbifold rigidity (F3), which are bypassed by
Lemma 4; the certified character of E (E ≅ χ₈); CERT_P6 §2.6 (λ₂ ≡ λ₃ along F_B).

## Reviews: verdicts and how each GAP was resolved

| Review | Scope | Verdict | GAPs / fixes and resolution |
|---|---|---|---|
| `review_1.md` | RESULTS §2 (criterion, numerics, literature) | **ACCEPT WITH FIXES** | #1 (Lemma proof incomplete; Kato/Riesz citation; T_B^G = 0; both uses of tr A ≡ 0): Lemma written out (Lemma 10). tr A ≡ 0 now has two proofs (Lemma 4) without Bochner/(2,3,8). F1 proved (Lemma 2). #2 (coordinate-dependent Hessian eigenvalues): relabelled; replaced by the representation-theoretic reduction, now exact (Lemma 6). #3 (sorted-member Richardson guardrail): redone at n = 160/320 (RESULTS §5d); not used in the proof. #4 (Observation 2 wrong): retracted; CERT_P6 §2.6 gives a hidden-involution explanation (not needed). #5 (Str16 quote/wording): fixed in LITERATURE. #6 (docstring): fixed. #7 (single discretisation family): superseded by the certificates; review_sign_1 adds an independent P1 code. #8 (SU13 locates *an* eigenvalue): handled in Lemma 1c |
| `review_cert_1.md` | CERT_RESULTS §1–§4 (LP, P1 RR, tr A, A ≠ 0) | **CONFIRMED**; one minor GAP, two remarks, nothing WRONG | Row 3 GAP (SU13 midpoint 3.83888725884, not 3.8388872588; upper end 4.2e-11 short): stated correctly in Lemma 1c; the upper end is covered by the certified P2 bound 3.8388872830 (G1 note); §4 needs only λ₁ ≥ 3.83788 under monotone endpoint evaluation (not under plain ball assembly; see Lemma 5). Remark 4.4′ (β_max "thin strip"): corrected to enclosure overestimation, true sup ≈ 3.3 (Lemma 5). Remark on SU13 provenance: D1. NaN safety: CONFIRMED safe (a NaN cannot pass). F1 phrasing for piecewise-smooth pullbacks: Lemma 2 (H1)–(H3) |
| `review_sign_1.md` | CERT_SIGN (P5) | **ACCEPT** (no WRONG; three presentation GAPs) | G1 (λ-ball upper end): covered by P2 λ₃ ≤ 3.8388872830 (Lemma 1c note); re-run with the exact SU13 low end gives −12.850504. G2 (F1 overstated for Claim S): Lemma 8 step 1 states that F1 is not needed for (Λ∘x)″ ≤ g″ itself. G3 (error breakdown): rows labelled one-at-a-time (Lemma 8 step 5). Sensitivity answer, λ₁ ≥ 3.8388029653: D1 |
| `review_theory_1.md` | CERT_P6, cert_chars.py, F1_PROOF, Lemma N | **Sound**; no mathematical error; two GAPs | Lemma N GAP (sketch only): reviewer's complete proof (§3), itself second-read in review_lemmaN; incorporated as Lemma 3. G3 GAP (F1 only on the 11-dimensional chart; rank dπ = 6 numerical): closed by CERT_G3 (Lemma 9; second-read in review_proof_1). Minor: s_mid wording (fixed, Lemma 6 consistency), Cai cited as an analogue (Lemma 2), Duran v3 (Lemma 2). Suggestion adopted: the antipodal-isometry argument for K₀ ⊂ ker A (Lemma 6 iv) |
| `review_lemmaN.md` | Lemma N proof of review_theory_1 §3 | **CORRECT WITH FIXES** | N1 (slice existence assumes π smooth): Fix 1, equivariant centre-of-mass smoothing (Lemma 3 Step 4). N2 (PSL₂ lift): moot under Fix 1. N3 (mixed derivatives ∂_t∇_pφ): stated (Lemma 3 Step 4). N4 (splitting citation): Mazzoli Rem. 3.3 / Fischer–Marsden quoted, plus elliptic regularity for Y. N5 (relabellings): Fix 3, the relabelling line. Step 1 wording: corrected. Optional El Soufi–Ilias citation: added |
| review_proof_1 | CERT_G3, OPEN-1, whole document | ACCEPT-WITH-FIXES | fixes 1–9 above (applied) |
| `review_lb_1.md` | LB_PLAN, LB_RESULTS, lb_core/lb_bolza/lb_torus, lb_cert_sign, lb_cert_first_order (SU13-free lower bound, Lemma 1d, and the re-certified Lemmas 5 and 8) | **ACCEPT** (three cosmetic fixes, no number changed) | F1 (A5 needs ρ > 0; cite BMPW Thm 2.1 via Vejchodský as second source): added to LB_RESULTS §0, LB_PLAN §1(c), D15. F2 (json `certified_lambda1_lower` rounded to nearest, 2.5e-16 above the Arb endpoint): LB_RESULTS §0 and Lemma 1d quote the downward-rounded 3.8388757010527243 / the Arb string. F3 (Lemma 5 rows mislabelled SYM=1; script default SYM=0): label withdrawn in LB_RESULTS §4; reviewer's SYM=0 re-run −3.0383476316928864 recorded in Lemma 5. Scope note (row 6): Lemma 1d does not give multiplicity; FBP 7.2–7.3 still needed (stated in Lemma 1d and D4) |

Not reviewed and not needed: CERT_P6 §2.6.

Integrated-proof read: review_proof_1.
