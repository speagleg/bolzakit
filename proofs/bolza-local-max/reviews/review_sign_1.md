# Second reading of CERT_SIGN.md (P5: certified g″(0) ≤ −12.85 along F_B)

Reader: independent adversarial reviewer (track RS), 2026-10-08. The files reviewed are `CERT_SIGN.md`,
`cert_rr_p2.py` and `cert_sign.py`. Context came from README, RESULTS §0 and §2.3, CERT_PLAN (c-ii) and
CERT_RESULTS §1–4.

I wrote only this file and files under `scratch_RS/`. No existing file was edited.

**Verdict: ACCEPT (Claim S as stated, conditional on the listed inputs).**
- No WRONG items.
- There are three small GAPs, all documentation or presentation; none changes the number.
- The certificate fails when it should under every perturbation tried.
- Sensitivity to SU13 (controller's question): the sign is still certified for any rigorous λ₁ lower bound
  ≥ 3.8388030, i.e. ≈ 8.4e-5 below λ₁. SU13's claimed precision is 1e-6, so there is roughly 80× slack.

## 1. Verdict table

| # | Claim | Verdict | Evidence |
|---|---|---|---|
| 1a | Λ(p) ≤ g(p), equality at 0, hence (Λ∘x)″(0) ≤ g″(0) | CONFIRMED | See note 1a below. |
| 1b | The trial space uses the exact E, not an approximation | CONFIRMED | e_i = P_E u_i is never computed. Every quantity of g″(0) at e is rewritten as the u-value plus cross terms in w = u − e, and those cross terms are bounded. The identities used, q(e,·) = λm(e,·) and m(e,1) = 0, hold exactly for exact e. |
| 1c | g is C² | CONFIRMED | a_T(p) and ρ_T(p) are real-analytic in p, uniformly on each closed triangle (Q ≥ 1, D > 0). The e_i and ψ_i are fixed H¹ functions. So the Gram entries are analytic and g = tr(G_m⁻¹G_q) is analytic near 0. |
| 1d | The §1.2 formula for g″(0) | CONFIRMED | Re-derived by hand (Schur complement; R = Q − λS; R′(0) = A; R″(0) = X + 2λc₁c₁ᵀ/m₁₁). Independently tested on a random 9×9 quadratic pencil, with K(p)1 = 0, an exact 3-dimensional eigenspace, a non-orthonormal e, random ψ, and A, M₁, c₁ all nonzero. In mpmath at 60 digits, FD = 1341.4543723164838 = formula (relative difference 3e-39) (`scratch_RS/check_formula.py`). Ky Fan g ≥ λ₀+…+λ₃ was also checked there. |
| 2a | Exact conformity of P2 across the projective side pairings | CONFIRMED | Argument checked (note 2a). Code: edge dofs keyed by glued label pairs; χ = −2 asserted; exactly 2 triangles per key (this also excludes multi-edges, which would give a count of 4). My P1 code uses its own gluing and a different mesh, and reproduces area = 4π and the cluster. |
| 2b | Element integrals rigorous (series + remainder) | CONFIRMED | See note 2b below. |
| 2c | Element p-jets | CONFIRMED | Each Taylor coefficient of Q^{−1/2}·N and D·Q^{−3/2} was re-derived. The J-class rules for sqrt, acosh, tan, sin, cos, sinh, cosh and division were checked. Independent quadrature against FD in p (own F_B vertices by root-finding, own node rule): K′ and M′ agree to ≤ 3e-8, K″ to ≤ 3e-8, and M″ to ≤ 2e-6 (relative to max entry; FD-limited). |
| 3a | Gap inequalities m(w) ≤ exc/(λ₄−λ), q(w) ≤ exc·λ₄/(λ₄−λ) | CONFIRMED | Spectral decomposition. u ⊥_m 1 exactly, since c_i is computed in Arb with exact Tm. exc = q(u) − λm(u) is an Arb enclosure over the whole λ-ball. Both bounds are monotone in λ₄, so the lower bound 4 is valid. |
| 3b | Inputs λ₄ > 4.0 and λ₁ ∈ SU13 | CONFIRMED as cited | λ₄ > 4 comes from the LP (CERT_RESULTS row 1+2), which was not re-checked here. λ = λ₁ = λ₂ = λ₃ uses FBP 7.2–7.3 and SU13. See G1 on the λ-ball. |
| 3c | The Rayleigh quotients are certified | CONFIRMED | q(u_i) and m(u_i) are exact Arb sums of Arb element matrices against exact float vectors. Only ub() (an upward-rounded float) is applied before the sqrt. |
| 3d | η = √q(w) ≤ 5.03e-3 at n = 96 | CONFIRMED (re-run) | Re-run gives 5.0299e-3, 5.0299e-3 and 5.0425e-3. |
| 4a | Every term of g″ is transferred | CONFIRMED | Expanded entry by entry (note 4a). No omitted cross term. |
| 4b | No float arithmetic in the bounds | CONFIRMED | All bounds are Arb. The only floats are ub(), which rounds up by 4e-16·\|x\| on top of float(upper), and the remainder radii in element_p2, which are upward-rounded by (1+1e-12). Both are safe. |
| 4c | β enclosures | CONFIRMED | See note 4c below. |
| 4d | The total transfer error is 1.42 | CONFIRMED (re-run) | n = 96: X1 +0.839, X2 +0.570, M0 +4.6e-4. Upper bound −12.850476. |
| 5a | F_B is a valid smooth path of glued hyperbolic surfaces through B | CONFIRMED | Independent mpmath root-finding of "angle at V₀ = π/4" from raw hyperboloid geometry. At p = 0, 0.1 and −0.05: all 8 angles = π/4 to 1e-39, and L_k = L_{k+4} to 1e-40. cosh R jet value and second derivative (16.48528137…) match FD; the first derivative is 0 in both (`check_fb.py`). |
| 5b | Nonzero tangent | CONFIRMED | In ℝ¹⁶ the tangent is not an so(2,1) orbit direction (residual 8.1 vs \|v\| = 11.5). Stronger: if dπ(v) = 0 in Teichmüller space, then (Λ∘x)″(0) = dΛ(π∘x)″ = 0 because dΛ = 0. That contradicts the certified ≤ −12.85, so the certificate itself proves dπ(v) ≠ 0 (given F1 on C). |
| 5c | Pullback = spectrum of x(p) | CONFIRMED | The mesh map is piecewise smooth, bi-Lipschitz and gluing-equivariant (isometries act linearly on W and preserve s, and the arclength node rules commute with isometries). The fan stays star-shaped for small p. |
| 6 | Re-runs and perturbations | CONFIRMED | §3 below. |
| G1 | λ-ball LAM | GAP (harmless) | See G1 below. |
| G2 | Dependency 4 (F1) overstated for Claim S | GAP (presentation) | (Λ∘x) along a 1-parameter analytic family of forms is analytic by Kato (type (a) holomorphic family; isolated cluster, gap ≥ 0.16). So the first bullet of Claim S does not need F1. F1 (C² on C or T₂) is needed only to read the result as D²Λ(v_B, v_B) and for P6. |
| G3 | §3.2 table | GAP (cosmetic) | The X1/X2/M0 breakdown is "each group switched on alone", so the groups do not add to the total exactly (0.839 + 0.570 + … ≈ 1.41 vs 1.42 from the joint ball). Not an error. |

**Note 1a (Λ ≤ g).**
- The argument: λ₀ ≡ 0. V(p) = span{1, e_i + pψ_i} ⊂ H¹ is 4-dimensional for small p. Ky Fan for the closed form
  pencil (q_p, m_p) on H¹(B) gives λ₀+…+λ₃ ≤ tr(G_m⁻¹G_q). g(0) = 3λ.
- h = g − Λ∘x ≥ 0 with h(0) = 0, so h″(0) ≥ 0 whenever Λ∘x is twice differentiable at 0.
- Even without that, the limsup of the second difference quotients of Λ∘x is at most g″(0).

**Note 2a (conformity).**
- The edge restriction is X(s) = (sP_i + (1−s)P_j)/√Q, and a P2 polynomial restricted to the edge is quadratic
  in s.
- The pairing isometry is linear on ℝ³ and B-orthogonal, so gX(s) is the same-s point of the paired edge.
- The P2 edge basis is symmetric under s ↔ 1−s.

**Note 2b (element integrals).**
- The Lagrange remainder bound |C(−s,K+1)|·ε_max^{K+1}(1−ε_max)^{−s−K−1} is correct for ξ ∈ [−ε_max, ε_max].
- ε_max is a sup over the square [−1/3, 2/3]² ⊇ the simplex, and ε_max < 0.25 is asserted.
- The remainder is integrated against monomials a^i b^j ≥ 0.
- Polynomial degree ≤ 8 is covered, and a missing key would raise rather than truncate.
- The metric formula g = (LLᵀ − QG)/Q² was re-derived: √det g = |D|Q^{−3/2} and g⁻¹√det g = adj(LLᵀ−QG)/(|D|√Q).
- Independent check: Klein-model metric with analytic Jacobian, 576-point Duffy–Gauss, own node positions. At
  n = 12 and n = 48, over 6 triangles each (centre, side strip, β-maximal strips, other sectors), K₀ and M₀ agree
  with the Arb midpoints to ≤ 1.3e-13 relative (`scratch_RS/check_elements.py`).

**Note 4a (transfer of each term).**
- m(e_i,e_j) = m(u_i,u_j) − m(w_i,w_j). The code's bound sm_i·sm_j is valid (tighter than its generic form).
- b(e,e): the cross terms with W1q/W1m, and the quadratic term with β₁max and β₁ρmax.
- (q″−λm″)(e,e): the same with β₂.
- 2b(ψ,e) is transferred both ways (r1).
- (q−λm)(ψ,ψ) has no e.
- m′(e,e), m(ψ,e) with √m(ψ)·√m(w), and c₁ with W1m[1]·sm.
- A, M₁ and c₁ enter the A-term and the c₁-term, and are all ball-propagated through arb_mat inv.
- The 3×3 algebra is pure Arb.

**Note 4c (β enclosures).**
- The 3 boxes [0,½]², [½,1]×[0,½] and [0,½]×[½,1] cover the simplex.
- a⁻¹ȧ is similar to the symmetric a^{−1/2}ȧa^{−1/2}, so Frobenius ≥ spectral radius gives −βa ≤ ȧ ≤ βa.
- The second-derivative matrix 2(X₂ + r₁X₁ + r₂I) = a⁻¹ä was re-derived.
- Pointwise sampling with independent geometry and FD in p, over ~150 elements at n = 12 and n = 48: the observed
  ratio to β is ≤ 0.70 (β₁), 0.89 (β₂), 0.99 (β₁ρ) and 0.97 (β₂ρ). No violation; ρ's bound is nearly tight
  (`check_beta.py`).

**G1 (the λ-ball LAM).**
- LAM = 3.8388872588 ± 1.0000001e-6, while SU13's midpoint is 3.83888725884.
- The **lower** end is 4.2e-11 *below* SU13's lower end, so it is conservative. The lower end is the one that
  matters (it enters exc).
- The **upper** end is 4.2e-11 *below* SU13's upper end, so LAM does not contain SU13's full interval.
- This is harmless, because the certified P2 RR (CERT_SIGN §2) gives λ₁ = λ₃ ≤ 3.8388872830 anyway. That should
  be stated, or LAM should be set to [SU13_low, 3.8388872830].
- A re-run with low = 3.83888625884 exactly gives −12.850504 (n = 96).

## 2. Weakest λ₁ lower bound that still certifies (controller's question)
`scratch_RS/bisect_lam.py` sets LAM = [low, 3.8388882588] and bisects on low. Element data are cached, and the
script is otherwise unchanged `cert_sign.main`.

| n | certifies for λ_low ≥ | fails at | λ₁(SU13 mid) − threshold | g″ upper at the threshold |
|---|---|---|---|---|
| 48 | 3.8388093904 | 3.8388093576 | 7.8e-5 | −0.002 |
| 96 | **3.8388029653** | 3.8388029297 | **8.4e-5** | −0.0005 |

Other points along the n = 96 curve:

| λ_low | g″(0) ≤ | √q(w) |
|---|---|---|
| SU13 low | −12.85 | 5.0e-3 |
| μ − 3.2e-5 | −5.90 | 2.8e-2 |
| μ − 7.5e-5 | −0.90 | 4.3e-2 |
| μ − 1.8e-4 | +7.5 | 6.6e-2 |

**Consequences.**
- The sign certificate needs a rigorous λ₁ lower bound accurate to ≈ 8e-5 absolute (2.2e-5 relative). SU13's
  1e-6 has ~84× slack.
- If SU13's 1e-6 is *relative*, then ±3.84e-6 absolute is still certified. The interpolated bound is about −9.5
  (compare n = 48, where radius 1e-5 gave −9.43).
- What remains is the *validity* of SU13's bound (an unpublished-code claim), not its precision. Any independent
  rigorous lower bound λ₁ ≥ 3.83881 (e.g. a Lehmann–Goerisch / Liu-type guaranteed lower bound, or an MPS
  a-posteriori bound) would replace SU13 entirely.
- The FE side contributes almost nothing to η: μ₉₆ − λ₁ ≈ 2e-8.

## 3. Re-runs and perturbations (work/.venv, nice, OMP_NUM_THREADS=2)

| run | result | expected |
|---|---|---|
| n = 48, SYM=1 | upper −12.566084, central −14.2742, √q(w) ≤ 5.69e-3 (30 s) | matches CERT_SIGN |
| n = 96, SYM=1 | upper **−12.850476**, X1 0.839, X2 0.570 (199 s) | matches |
| n = 96, SYM=0 (no symmetry shortcut) | upper −12.850476213730866, identical to SYM=1 (421 s) | matches; new evidence for the shortcut at n = 96 |
| n = 48, λ-radius 1e-5 | −9.43, certified | degrades as expected |
| n = 48, λ-radius 1e-4 | +2.13, **not certified** | must fail ✓ |
| n = 48, λ₄ = 3.85 | −7.63, certified | degrades ✓ |
| n = 48, λ₄ = 3.84 (gap 1e-3) | +9.26, **not certified** | must fail ✓ |
| n = 24, λ = 3.8389 ± 1e-8 (excludes the true λ, exc < 0) | NaN → `arb_mat.inv` raises; no certificate | must fail loudly ✓ |
| n = 48, ψ = 0 (no corrector) | g″ ≤ +46.39 (= T₂), **not certified** | must fail ✓ |

Independent float cross-check (`scratch_RS/indep_p1.py`, `indep_p1.out`) uses my own code throughout:
- F_B by root-finding;
- a different mesh, projective-uniform subdivision of each sector rather than octagon.Mesh's arclength fan;
- my own gluing;
- conforming P1;
- Klein metric;
- FD in p over re-meshed surfaces (no pullback).

| n | dofs | λ₁..₃ | Λ″ (5-point, h = 0.02 / 0.01) |
|---|---|---|---|
| 16 | 1 022 | 3.874, 3.884 ×2 | −11.630 / −11.633 |
| 32 | 4 094 | 3.849, 3.851 ×2 | −13.273 / −13.276 |
| 64 | 16 382 | 3.8419, 3.8422 ×2 | −13.953 / −13.956 |
| 128 | 65 534 | 3.83968, 3.83972 ×2 | −14.182 / −14.185 |

- The increments shrink by 2.4 and then 3.0, approaching h².
- Richardson (ratio 4, n = 64/128) gives **Λ″ ≈ −14.26**. This agrees with the certificate's central g″ = −14.274
  and its FD of the P2 cluster (−14.274).
- The ordering Λ″ ≤ g″ is consistent with Λ″ ≈ g″: the correctors are nearly optimal.

## 4. External inputs (complete list for Claim S)
1. λ₄(B) > 4.0: the trace-formula LP (CERT_RESULTS 1a/1b, plus row 2). It cites the Selberg trace formula in
   FBP/Buser normalisation, FBP §7.4 geodesic lengths and lower multiplicities (12, 12, 24), and the systole
   (Jenni/Schmutz via FBP).
2. λ₁ = λ₂ = λ₃ (single eigenvalue on E): FBP Lemmas 7.2–7.3.
3. A rigorous λ₁ lower bound ≥ 3.8388030. SU13's 3.8388862588 is used; see §2 for the slack. SU13's claim is a
   one-sentence statement about unpublished code.
4. Ky Fan trace minimum principle for closed form pencils on H¹ (standard).
5. Kato analytic perturbation for the isolated cluster sum along the 1-parameter family. This suffices for
   (Λ∘x)″(0) ≤ g″(0) (G2). F1 (C² on C) is needed only for the D²Λ(v_B, v_B) reading and for c < 0.
6. dΛ(x_B) = 0 on T_{x_B}C: CERT_RESULTS §3, with certified rank J = 5 plus F1. It is used only for the D²Λ form
   and the nonzero-tangent corollary.
7. Poincaré polygon theorem: L_{k+4} = L_k and angle sum 2π give a closed surface with the Bolza gluing pattern.
8. P6 (v_B ∈ K, K irreducible) for c < 0. This is not provided by P5.
9. Software: python-flint 0.9.0 / Arb correctness. numpy and scipy are used only for untrusted proposals.

## 5. Files (scratch_RS/)
- `perturb.py`, `perturb48.out`: perturbation runs.
- `bisect_lam.py`, `bisect48.out`, `bisect96.out`: λ-threshold.
- `check_formula.py`: exact pencil test of §1.2.
- `check_fb.py`: F_B geometry and jets.
- `check_elements.py`: independent element integrals and p-derivatives.
- `check_beta.py`: pointwise β sampling.
- `indep_p1.py`, `indep_p1.out`: independent FE Λ″.
- `run48.out`, `run96.out`, `nosym/run96_nosym.out`: re-runs. The json files written by those runs are
  `cert_sign_48.json`, `cert_sign_96.json` and `nosym/cert_sign_96_nosym.json`.
