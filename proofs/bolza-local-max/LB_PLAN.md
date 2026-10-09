# LB_PLAN: an SU13-free rigorous lower bound λ₁(Bolza) ≥ 3.83881

Written 2026-10-09, before any substantial coding. Goal: replace PROOF.md §Dependencies D1 (SU13) with an
in-house Arb certificate λ₁ ≥ L, with L ≥ 3.8388029653 required and L ≥ 3.83881 desired.

## 0. What is available without SU13
- (LP-A) Spec(B) ∩ (0, 3.3] = ∅ and (LP-B) #Spec(B) ∩ [3.3, 4.0] ≤ 3, hence **λ₄ > 4.0** (Lemma 1a; trace formula
  + FBP geodesic data; no SU13, no FBP multiplicity theorem).
- λ₃ ≤ 3.8388872830 (P2 conforming Rayleigh–Ritz, Lemma 8 / `cert_rr_p2.py`).
- The exactly conforming hyperboloid-barycentric P2 space and its element integrals ∫ poly · Q^{−s} on every
  geodesic triangle (`cert_rr_p2.py`), with the pulled-back coefficients a_T = adj(LLᵀ − QG)/(D√Q),
  ρ_T = D Q^{−3/2}.

The gap ρ_gap := 4.0 − 3.8388873 ≈ 0.161 is what a Temple/Lehmann-type bound runs on.

## 1. Routes considered

### (a) Symmetry reduction + guaranteed lower bounds on a planar polygon (Liu / CR / Carstensen–Gedicke)
- Reduction: E ≅ χ₈ (octahedral). For a reflection s with tr(s|E) = +1 (e.g. s_mid, s_vert, CERT_P6 row 7),
  E has a 1-dimensional (−1)-eigenspace, so the s-odd sector contains a λ₁-eigenfunction and no constants.
  Iterating with a reflection group gives a mixed D/N problem on a fundamental polygon (down to the (2,3,8)
  triangle of area π/24 for the full group Δ(2,3,8)).
- Liu 2015 (AMC 267, 341–355), abstract (via the journal record): "proposes a universal framework to give
  explicit lower and upper bounds for the eigenvalues of eigenvalue problems of self-adjoint differential
  operators. In the case of the Laplacian operator, by applying Crouzeix–Raviart finite elements, an efficient
  algorithm is developed…". The framework is abstract (a(·,·), b(·,·), projection P_h with explicit error
  constant C_h): λ_k ≥ λ_{k,h}/(1 + C_h² λ_{k,h}) (quoted in Liu–Plum arXiv:2512.23182, eq. (1.1)).
  Variable coefficients/weights are covered *in principle* (the theorem is abstract) but C_h must be re-derived:
  for −div(A∇u) = λwu the CR interpolation constant picks up sup/inf ratios of A and w per element. In the
  Klein model w = (1−|x|²)^{−3/2}, A = (1−|x|²)^{−1/2}(I − …) vary by O(1) across the polygon.
- Precision: Liu–Plum Remark 3.4 (quoted): "the lower eigenvalue bound obtained … only converges at rate
  O(h)". We need relative precision 2e-5 (7.7e-5/3.84). With C_h = O(h) and CR rate h², λ − λ_low ≈ C_h²λ² ⇒
  h ≲ 1e-3 in hyperbolic units on a polygon of diameter ~1: ~10⁶–10⁷ CR elements, each needing rigorous
  variable-coefficient quadrature. Feasible in principle, very expensive, and the C_h re-derivation for
  weights is new work.
- **Verdict: rejected as primary** (cost); the projective-gluing obstruction of CERT_PLAN §1(b) is indeed
  absent on a single polygon, but precision O(h) kills it at 1e-5.

### (b) Selberg trace formula with the exact length spectrum (BSV-style)
- BSV (Booker–Strömbergsson–Venkatesh, IMRN 2006) as described in Seymour-Howell arXiv:2201.08760 §1/§2.3
  (quoted): verification uses the trace formula "with Hecke operators … The inclusion of the Hecke operators
  allows us to use linear algebra to help pick out individual eigenvalues"; the enclosure comes from a
  Rayleigh quotient ε_i² = Q(c_i, H̃_i)/Q(c_i, H) with H̃(λ) = H(λ)(λ − λ̃)², "a weighted average of
  (λ − λ̃_i)²".
- Without Hecke operators (we have none usable for a cocompact quaternionic group in bolzakit), the plain trace
  formula localises an eigenvalue only to the resolution of the test function: f supported in [−L, L] ⇒
  f̂ of exponential type L (Paley–Wiener), so separating λ₁ from a window at distance δλ = 7.7e-5, i.e.
  δr = δλ/(2r) ≈ 2e-5, needs L ≳ 1/δr ~ 10⁴–10⁵. The number of closed geodesics up to L grows like e^L/L.
  FBP/our LP already fail beyond X ≈ 4.05 with three lengths.
- Complete enumeration up to moderate L (≈ 10–15) is certifiable with bolzakit's exact group, but it cannot
  reach the resolution needed.
- **Verdict: rejected** (infeasible resolution, not just cost).

### (c) Temple / Kato / Lehmann–Goerisch with a flux (H(div)) reconstruction. **CHOSEN**
- Temple: for mean-zero u ∈ D(Δ) with ρ = R(u) < 4 ≤ λ₄ and Spec ∩ (λ₁, 4) = ∅,
  λ₁ ≥ ρ − (‖Δu‖² − ρ²)/(4 − ρ). Needs u ∈ H² of the closed surface: conforming P2 is not.
- **Lehmann–Goerisch** replaces ‖Hu‖² by b_G(w,w) for any w with b_G(w, Tφ) = N(v, φ) for all φ
  ("Goerisch's XbT concept"). Quoted from Liu–Plum arXiv:2512.23182 §4.1:
  - Assumptions A1–A5 (A3: "Let X be a real vector space and b_G be a symmetric positive semidefinite bilinear
    form on X. Let T : D → X be a linear operator such that b_G(Tu,Tv) = M(u,v)"; A4: "Given v₁,…,v_n ∈ D, let
    w₁,…,w_n ∈ X satisfy b_G(w_i,Tv) = N(v_i,v) ∀v ∈ D"; A5: A0 = (M(v_i,v_j)), A1 = (N(v_i,v_j)),
    A2 = (b_G(w_i,w_j)), A := A0 − ρA1, B := A0 − 2ρA1 + ρ²A2, "Suppose B is positive definite. Let
    ν₁ ≤ … ≤ ν_n be the eigenvalues of Az = νBz").
  - Theorem 4.3: "Suppose ρ is a lower bound to the (m+1)th eigenvalue λ_{m+1} for some m ≥ n … then
    λ_{m+1−k} ≥ ρ − ρ/(1 − ν_k) (1 ≤ k ≤ q)" (q = number of negative ν).
  - *(Added after review_lb_1, fix F1.)* A5 also requires **ρ > 0** (review_lb_1 row 1a, from the PDF §4.1). Here
    ρ = 5, so it is met.
  - *Second source (added after review_lb_1, fix F1):* Behnke–Mertins–Plum–Wieners (2000), Thm 2.1, as
    reproduced in Vejchodský arXiv:1801.02561 §2, states the method directly in shifted form (shift γ). With
    γ = 1 their matrices are A0 − ρA1 and A0 − 2ρA1 + ρ²A2 with A2 = ‖σ‖² + ‖v + div σ‖², and their conclusion
    ([ρ − γ − ρ/(1 − μ_i), ρ − γ) contains ≥ i eigenvalues of the unshifted problem) is λ₁ ≥ 4 − 5/(1 − ν̄)
    (review_lb_1 row 1b). Goerisch–Haunhorst and the NPW book chapter were not read.
- **Our instantiation** (no new theory):
  - D = H¹(B), M = q + m (positive definite), N = m. Shifted eigenvalues λ̃ = λ + 1:
    λ̃₁ = 1 (constants), λ̃₂,₃,₄ = λ₁,₂,₃ + 1, λ̃₅ = λ₄ + 1 > 5. Take **ρ = 5, m = 4, n = 3**.
  - X = L²(TB) × L²(B), T u = (∇u, u), b_G = L² inner product. A3 holds.
  - For a trial v ∈ P2 and any σ ∈ H(div, B), put **w = (σ, v + div σ)**. Then for all φ ∈ H¹:
    b_G(w, Tφ) = ∫σ·∇φ + ∫(v + div σ)φ = ∫vφ = N(v, φ), by integration by parts on the closed surface.
    **A4 holds exactly for every σ**; no divergence constraint has to be met. This is why the shift is used.
  - b_G(w_i,w_j) = ⟨σ_i,σ_j⟩ + ⟨v_i + div σ_i, v_j + div σ_j⟩. Best σ ≈ ∇ũ with ũ = (Δ+1)⁻¹v.
  - k = 3 bounds λ̃₂ = λ₁ + 1. **No SU13 and no FBP multiplicity theorem** are needed: three trial functions
    give lower bounds for all of λ₁, λ₂, λ₃ at once.
- **σ space.** Raviart–Thomas RT₂ on every reference simplex (a,b), unknown j = √g σ (contravariant density):
  - div σ = (1/√g) ∂_i j^i, √g = ρ_T = D Q^{−3/2};
  - |σ|²dA = jᵀ (a_T)^{−1} j da db with (a_T)^{−1} = (LLᵀ − QG)/(D√Q) (since det(√g g⁻¹) = 1 in 2D and the
    Gram identity det(LLᵀ − QG) = Q D²; to be checked symbolically);
  - outward flux through an edge per unit edge parameter s is intrinsic (ι_σ dA); the edge point
    W = sP_A + (1−s)P_B is the same on both sides, including across the projective gluings (isometries are
    linear on W). So matching the RT edge DOFs (flux values at s ∈ {0, ½, 1}, opposite signs) makes σ exactly
    H(div)-conforming on the glued surface. Same mechanism as the P2 conformity argument.
  - All integrals are ∫ poly · Q^{−s} with s ∈ {½, 3/2} or polynomial: ‖σ‖² (deg 8, Q^{−1/2});
    ‖v + div σ‖² = ∫[D Q^{−3/2}v² + 2v ∂·j + (∂·j)² Q²·Q^{−1/2}/D].
- **Cost / sensitivity estimate.** With A1 normalised, the bound loses ≈ 30·(energy error of v) +
  ≈ 730·‖w − Tũ‖². P2 at n = 48 has eigenvalue error ≈ 3e-7 ⇒ ≈ 1e-5 loss; RT₂ flux error O(h³). Expected
  L ≈ 3.83887 at n = 48, comfortably ≥ 3.83881. A float SPD solve (~2e5 unknowns) plus Arb element integrals:
  minutes.
- Fully rigorous: the only floats are the trial coefficients (untrusted, they only enter as exact data).

### (c′) Lehmann–Goerisch with H(div) Whitney fluxes was named in CERT_PLAN §1(b) as a "research project";
that assessment was for *guaranteed lower bounds via CR + gluing*. The present use needs only a priori
λ₄ > 4 (already certified) and a conforming flux, which the reference-coordinate RT space provides directly.

## 2. Choice
- **Primary: route (c)**, Lehmann–Goerisch (n = 3, ρ = 5) with P2 trial functions and RT₂ fluxes on the existing
  mesh. Final inequality: λ₁ ≥ ρ − ρ/(1 − ν̄) − 1, where ν̄ < 0 is an Arb-certified upper bound for ν₃ (Sylvester
  test of ν̄B − A ≻ 0 and of B ≻ 0).
- **Fallback 1:** same with n = 1 (one trial function) plus FBP Lemmas 7.2–7.3 (λ₁ = λ₂ = λ₃), already used in
  PROOF Lemma 1c; or RT₃ / finer n if precision is short.
- **Fallback 2:** route (a) on the s-odd sector with Liu's projection bound (very expensive).

## 3. Validation (known answer first)
- Flat square torus ℝ²/ℤ² with the *same* code (Q ≡ 1, M = JᵀJ constant): λ₁ = 4π² (multiplicity 4),
  λ₅ = 8π². Use n = 4 trial functions, ρ = 1 + 78.9 ≤ 1 + 8π². The certified lower bound must be ≤ 4π² and
  close to it.
- Negative control on Bolza: ρ above the true λ₄ + 1 (≈ 6.354) must be able to produce a "bound" above the
  certified upper bound 3.8388872830 (shows the hypothesis ρ ≤ λ̃₅ is load-bearing).
- Consistency: L ≤ 3.8388872830 (P2 upper bound).

## 4. Downstream
If certified, copy cert_sign.py / cert_first_order.py to lb_cert_sign.py / lb_cert_first_order.py with the λ
ball replaced by [L, 3.8388872830] and re-run (n = 96 and 96/192 respectively).

## 5. Files (new only)
lb_core.py (geometry-agnostic Arb element integrals for P2 + RT₂), lb_torus.py (validation), lb_bolza.py
(certificate), lb_*.json/.out, LB_RESULTS.md, scratch_LB/.
