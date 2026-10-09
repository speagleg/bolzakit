# Second-reader review of CERT_RESULTS.md §1–§4 (track R, 2026-10-08)

**Scope.** This review covers CERT_RESULTS §1 (trace-formula LP), §2 (conforming Rayleigh–Ritz), §3 (tr A ≡ 0) and
§4 (A ≠ 0). It does not cover CERT_SIGN, CERT_P6 or F1_PROOF.

**Method.**
- I re-ran every script on copies in `scratch_R/`, using `work/.venv`, `nice -n 10` and `OMP_NUM_THREADS=2`.
- I wrote independent checks (`scratch_R/indep_lp.py`, `check_rr.py`, `check_fo*.py`, `check_rank.py`, the
  negative controls and `nan_inject.py`).
- I downloaded FBP (arXiv:2111.14699) and SU13 (arXiv:1110.2150) and compared every quoted statement with the PDF
  text.

**Verdict.** §1–§4 are **CONFIRMED** as CERTIFIED computer-assisted statements, under the inputs listed at the end.
There is one minor GAP, which is numerically harmless and fixed in one line: row 3's interval endpoint and the λ₁
ball in the code. There are also two remarks: the SU13 input is weak in its provenance, and the β_max explanation is
inaccurate. No float leaks into a certified inequality, and I found nothing WRONG.

## Claim table

| # | Claim | Verdict | Evidence |
|---|---|---|---|
| 1.1 | Lemma A is FBP Lemma 3.2 without a ≤ λ₁, and the proof carries over | CONFIRMED | FBP Lemma 3.2 and (3.1) are quoted verbatim and correctly in LITERATURE §5.1 (checked against the PDF, `scratch_R/dl_fbp/fbp.txt` l.385–449). In FBP's proof, a ≤ λ₁ only serves to put every nonzero eigenvalue in the region F ≥ 0. Lemma A gets the same from F ≥ c on [0,a] and F ≥ 0 on [a,∞). The count in (0,a] then follows from S − F(0) ≥ c·#. |
| 1.2 | Fourier pair and F(λ) = s·v(s²(λ−¼))·e^{−s²(λ−¼)/2} for all λ ≥ 0 | CONFIRMED | Hermite functions are eigenfunctions of the unitary FT with eigenvalue (−i)^{2n} = (−1)^n, and scaling gives the factor s. Independent mpmath quadrature of (2π)^{−1/2}∫f e^{−iξx} with the stored rational c_n agrees to 15 digits at ξ = 0, 0.7, 1.9 and at ξ = i/2 (F(0) = 43.28502984 for A, −927.0839188 for B). `indep_lp.out` |
| 1.3a | f ≤ 0 on [sys, ∞) (u ≤ 0 on [y₀, ∞)) | CONFIRMED | Exact sympy Sturm count. I checked the semantics of `count_roots(y0, None)`: it counts over [y0, +∞), closed. Independent flint `fmpq_poly.complex_roots`: the only real root of u is y = 9.2242 < y₀ = 9.3461 (run A) and −3.65 (run B), and the leading coefficient of u is < 0. y₀ < ℓ₁²/s² is asserted in Arb. |
| 1.3b | F ≥ 0 beyond the window (v ≥ 0 on [s²(b−¼), ∞)) | CONFIRMED | Same method. The real roots of v are −0.79 (A) and 4.697 < 6.3375 (B), and the leading coefficient of v is > 0. |
| 1.3c | c = rigorous lower bound of F on the window | CONFIRMED | 4000 Arb balls. Coverage of each closed subinterval is asserted in ball arithmetic, and `lower()` is used. On a 2001-point grid, F's minimum is ≈ 1.0 (the LP constraint); the certified 0.98164 / 0.99643 are conservative. |
| 1.3d | Integral term I, including the tail | CONFIRMED | The exact Gaussian-moment part was checked by hand (∫ r(sr)^{2k}e^{−s²r²/2} dr = 2^k k!/s²). `acb.integral` on [0,12]: the integrand is meromorphic with poles at Im z = ±½, and Arb returns a non-finite result on balls containing a pole, so this is rigorous. The tail bound e^{−24π}·(4/s)·Σ\|v_k\|2^k k! is correct. Independent mpmath values: I = 51.153059786184 (A), −336.259613312622 (B), matching the Arb balls. |
| 1.3e | Geodesic inputs: sys = ℓ₁; unoriented multiplicities 12/12/24 of ℓ₁, ℓ₂, ℓ₃; primitive since ℓ₃ < 2ℓ₁ | CONFIRMED (cited) | Checked verbatim against FBP §7.4 (PDF l.1238–1262). FBP divide AS88's oriented counts by 2, so the counts are unoriented, which matches the (3.2) normalisation ℓf/sinh(ℓ/2) used in the code. Only lower bounds are used, and f ≤ 0 makes dropping terms safe. ℓ₃ = 5.828070775 and 2ℓ₁ = 6.114283678 were recomputed. |
| 1.3f | No float leaks into the certified inequality | CONFIRMED | Floats appear only in: the LP proposal; rationalising c_n (after which everything is exact); the ball radius in step (3) (coverage is asserted); and `float(tail.upper())·1.0000001` (outward, and the tail is ≈ 1e-33). PASS is the Arb comparison `bound < target`. |
| 1.4 | Run A bound 0.91261 < 1; run B bound 3.86966 < 4; 4.05 fails | CONFIRMED | Re-run reproduces 0.912612352675… and 3.869661119769… (`rerun_lpA.out`, `rerun_lpB.out`). X = 4.05 gives 4.0536 and FAILS (`rerun_lpB405.out`). The mpmath recomputation gives the same I + G − F(0) (0.895859, 3.855840). |
| C1 | Spec(B) ∩ (0,3.3] = ∅; #[3.3, 4.0] ≤ 3 | CONFIRMED | 1.1–1.4. Run B applies FBP Lemma 3.2 verbatim, with a = 3.3 ≤ λ₁ from run A. |
| 2.1 | The hyperboloid-barycentric P1 space is exactly H¹-conforming on the glued surface | CONFIRMED | **Argument:** the trace on an edge is determined by its endpoints, since X ∝ aP₀ + (1−a)P₂ fixes a. Isometries are linear in W. The mesh is a fan with no hanging nodes; I checked the up/down triangles in `Mesh`. **Check** (`check_rr_12.out`): I built the pairing isometries g_k independently from frames (mpmath). They satisfy det = 1, gᵀJg = J to 4e-26, g₀₀ > 0, and their translation length is exactly ℓ₁ = 2arccosh(1+√2). They map the Arb side node j of side k+4 onto side node n−j of side k (residual 2.6e-26), and the `Mesh.lab` labels agree. All 8 vertices form one class, and all 1152 triangles (n = 12) have det[P₀P₁P₂] > 0, so the mesh is a straight-line Klein-model triangulation with no fold-over. At n = 192 the area enclosure 4π·[1 ± 2.2e-13] supports this. |
| 2.2a | Exact element formulas √det g = D·Q^{−3/2}, g^{−1}√det g = adj(LLᵀ−QG)/(D√Q) | CONFIRMED | Re-derived: g = (LLᵀ − QG)/Q². The Gram identity Q·detG − Lᵀ adjG L = det[W₁W₂W]²·det B = D² holds, with det B = +1. Independent check: mpmath 2-D quadrature of the metric, obtained by numerically differentiating X(a,b) (no closed form), on 6 elements including sides and corners. All K and M entries lie inside the Arb enclosures, with worst relative midpoint difference 8.7e-9. |
| 2.2b | Truncated-series remainder enclosure | CONFIRMED | The Lagrange remainder \|C(−s,K+1)\|·ε_max^{K+1}(1−ε_max)^{−s−K−1} is valid for ε ≥ −ε_max, and the centred sup bound for ε_max is valid on the simplex. The ∫\|p\| ≤ sup\|p\|·½ step is fine. Test: enclosures with K_SERIES = 0, 1, 2, 3 all overlap, with radii 1.1e-1 → 1.0e-4. An under-estimated remainder would have shown up at K = 0. |
| 2.3 | Min-max: λ₃(B) ≤ μ_max(K_p, M_p) on V = span{1, x₁, x₂, x₃} | CONFIRMED | Indexing is λ₀ = 0, so 4 dimensions bound λ₃. Positive definiteness of tM_p − K_p, with K_p ⪰ 0, implies M_p is nonsingular (so dim V = 4) and every Ritz value is < t. Negative control: t = μ₃ − 1e-7 and t = μ₃ + 1e-9 are both rejected (`perturb_rr_48.out`). The constant column is included, so multiplicities are handled correctly. |
| 2.3b | λ₃ ≤ 3.8391848048 (n = 192) | CONFIRMED | The re-run reproduces `certified_mu3_upper = 3.839184804702046` exactly, with area/4π = [1 ± 2.2e-13] (`rerun_rr_192.out`). n = 48 also reproduces (3.84363818843). |
| 2.4 | SYM = 1 shortcut | CONFIRMED | The radii R + t(−1)^k, the centre (normalised vertex sum) and all node formulas are invariant under r² ∈ SO(2,1). Element data depend only on B-products and \|det\|. `Mesh.T` is built with an outer loop over sectors and identical local order. The full computation and the shortcut match (author check), and the 55,296-pair overlap check re-runs. |
| C2 | Isolation: λ₁, λ₂, λ₃ ∈ (3.3, 3.8391849], λ₄ > 4.0, gap ≥ 0.1608 | CONFIRMED | From C1 and 2.3b. The window [3.3, 4.0] is closed, so λ₄ > 4 is strict. This uses neither SU13 nor FBP's multiplicity theorem. |
| 3 | λ₁ = λ₂ = λ₃ ∈ [3.8388862588, 3.8388882588] | **GAP (minor)** | m₁ ≥ 3 from FBP Lemmas 7.2 and 7.3 (j = 1) is quoted correctly. The proofs in FBP are symmetry and nodal-domain arguments with no eigenvalue input, as claimed. **But:** SU13 give 3.83888725884… ± 1e-6, so the upper end is 3.83888825884, not 3.8388882588: the stated interval is 4.2e-11 too short on top. `cert_first_order.py` uses `arb("3.8388872588", "1.0000001e-6")`, whose extra radius of 1e-13 does not cover this 4.2e-11. **Fix:** use midpoint `3.83888725884` or radius `1.0001e-6`. The effect on §4 is about 1e-10, against a margin of 3.09. **Provenance:** SU13 §5 states the ±1e-6 bound in one sentence. The rigorous code is not published, and the sentence does not say whether the bound is absolute or relative (relative would mean ±3.84e-6). See the sensitivity row 4.6. |
| 3′ (§3) | tr A ≡ 0 via the r-action, without Bochner or (2,3,8) | CONFIRMED (logic), conditional on F1 | r is linear, r⁸ = I, r maps C to C (relabelling k → k+1 preserves the pattern k ↔ k+4), and Λ∘r = Λ. Hence dΛ = dΛ∘P with P the average over ⟨r⟩. (ℝ¹⁶)^r = {δ_k = Rot^k δ₀} is 2-dimensional. The radial direction is ∉ T because d(Σθ)/dR ≠ 0, and the tangential direction is the rotation, along which dΛ = 0. Independent check (`check_rank.out`): J has singular values 82.0, 82.0, 66.9, 34.5, 34.5 (rank 5); r preserves ker J; ker[J; r−I] is 1-dimensional, with \|cos\| = 1.000 against the rotation field; tr(r \| ker J) = −1.000, matching CERT §5's hand computation. Remark: F1 is needed in the form "Λ is C¹ on C for the mesh-map (piecewise-smooth, bi-Lipschitz) pullback family". The RESULTS §2.1 wording, "C² family of diffeomorphisms", should be adjusted to say this. |
| 5 | rank J = 5 certified; d(Σθ)/dR = −5.5729711960 | CONFIRMED | Re-run reproduces the minor det 46151927.03… and the derivative, which agrees with the closed form to 1e-30 (`rerun_rank.out`). The angle and Klein-lift formulas were checked. |
| 4.1 | x(t) ⊂ C, σ′(0) = 0, tangent δ | CONFIRMED | Arb dual derivatives along δ of all 5 constraints are [± 2e-25] … [± 6e-22] (`check_fo.out` (c)). Dropping σ = O(t²) leaves first derivatives unchanged. |
| 4.2 | Dual-number element derivatives dK, dM | CONFIRMED | Central finite differences of `cert_rr.element` at R ± h(−1)^k (h = 1e-6, K_SERIES = 10, n = 24) agree with the dual values to 2.4e-13 relative on 13 elements (`check_fo24.out`). The formula ∫Ṅ Q^{−s} − s∫N Q̇ Q^{−s−1} was checked. |
| 4.3 | Criterion b(e,e) < 0 for some e ∈ E ⇒ A(δ) ≠ 0 | CONFIRMED | This holds trivially, since the diagonal entry of a nonzero e is nonzero. Tracelessness (§3) is not needed for A(δ) ≠ 0, only for interpreting it as an element of Sym₀. Identifying [b(e_i,e_j)] with the first-order cluster variation along the 1-parameter curve needs only 1-parameter degenerate (Rellich–Kato) perturbation theory for a C¹ family of forms, which is weaker than F1. |
| 4.4 | Transfer u → e = P_E u; m(w), q(w) bounds; cross/quadratic bounds; β bounds | CONFIRMED | Re-derived: μm(u) = λ₁m(e) + q(w) and q(w) ≥ λ₄m(w) give both bounds. The monotonicity in λ₁ and λ₄ is the right way round, and Arb takes care of the λ₁ ball. \|xᵀȧy\| ≤ β\|x\|_a\|y\|_a follows from −βa ≤ ȧ ≤ βa, then Cauchy–Schwarz over the elements. β_T ≥ ρ(N⁻¹Ṅ) + \|Ḋ/D + Q̇/2Q\| on the box [0,1]² ⊇ simplex. Sampled true pointwise sups always lie below the enclosures (`check_fo24.out`, `beta_probe.out`). m(e) > 0. |
| 4.4′ | "β_max ≈ 33 is a thin-strip artefact of the mesh map" | inaccurate remark (harmless) | At the side-midpoint elements, the sampled true sup is ≈ 3.3, against an enclosure of 33.4 (`beta_probe.out`). The factor of 10 comes from evaluating over the box [0,1]² rather than the simplex, plus the Frobenius norm and dependency blow-up. Overestimation is safe. For P5 (track S), subdividing the simplex when computing β should cut the transfer error a lot. |
| 4.4″ | NaN safety (controller's question) | CONFIRMED | python-flint returns NaN for `arb(0, r)**2` and for the sqrt of a ball straddling 0. That is the NaN source in `fro`. In `main`, `ba if ba > beta_a` silently skips a NaN in the *reported* β_max, but the NaN is stored per element and enters Sa (and Sr), so err, b_ee_upper and the verdict become NaN, and `NaN < 0` is False. Injection test (`nan_inject.out`): forcing a NaN β_a or β_ρ on one element at n = 48 gives `transfer_err = nan`, `A_nonzero_certified = false`, while the printed β_max stays finite (33.39). In the phase-1 runs Sa and Sr are finite (β_rms printed), so no element had NaN at n = 96 or 192. NaN does appear at n = 8 for large elements, where the assertion `detN > 0` also fails. |
| 4.5 | b(e,e) ≤ −1.028 (n = 96) and ≤ −3.088 (n = 192), times m(u) | CONFIRMED | The n = 192 re-run reproduces b_ee_upper = −3.0875694 (`rerun_fo_192.out`). |
| 4.6 | Robustness / negative controls | CONFIRMED | (i) A mesh-map-only deformation (centre node moved, same surface, so A = 0 exactly) gives b(e,e) ≤ +1.409: not certified, correctly (`neg_control_fo2.out`). (ii) The pure rotation gives NaN, so not certified. (iii) Sensitivity at n = 192: the certificate survives any lower bound λ₁ ≥ 3.83788 (1.2e-3 below SU13), and λ₄ ≥ 3.9 still gives −1.89. So the endpoint issue in row 3 and the absolute/relative ambiguity are both irrelevant numerically. |

## Inputs the certified statements depend on
1. **Selberg trace formula** in FBP normalisation (3.1), from Buser [Bus10, p.253]. Used in C1 and everything after
   it.
2. **FBP §7.4 geodesic data**, citing Jenni, Schmutz, Aurich–Steiner and ABS91:
   - sys(B) = 2arccosh(1+√2);
   - at least 12 / 12 / 24 unoriented primitive geodesics of lengths ℓ₁, ℓ₂, ℓ₃.

   Used in C1.
3. **FBP Lemmas 7.2 and 7.3** (m₁ ≥ 3), which rest on FBP Prop. 4.4 and Courant. Used in row 3 and §4.
4. **SU13 §5**: an eigenvalue lies within 1e-6 of 3.83888725884. This is an unpublished rigorous computation, stated
   in one sentence, and it does not say whether the bound is absolute or relative. §4 actually needs only
   λ₁ ≥ 3.8379. Used in row 3 and §4.
5. **F1**:
   - for §3: Λ is C¹ on the constraint manifold C under the mesh-map pullback;
   - for §4, only the 1-parameter first-order degenerate perturbation formula.
6. **Software:** python-flint 0.9.0 (Arb ball arithmetic, `acb.integral`), sympy's exact Sturm root counting, and
   the correctness of `octagon.Mesh`'s combinatorics. I checked the combinatorics independently: labels, pairings,
   orientation.

No input is misquoted. The only numerical slip is the row-3 upper endpoint, 4.2e-11.

## Suggested edits (for the controller)
- Row 3: write [3.8388862588, 3.8388882589]. In `cert_first_order.py`, use `arb("3.83888725884", "1e-6")`. Also add
  a sentence saying that §4 needs only λ₁ ≥ 3.8379.
- §4.5: replace the "thin-strip artefact" sentence with "mostly enclosure overestimation (box ⊇ simplex); true sup
  ≈ 3.3".
- `cert_first_order.py`: bound the Frobenius norm by its upper endpoint, e.g.
  `arb(sum(abs(m).upper()**2 ...)).sqrt()`. Note that `x*x` and `abs(x)*abs(x)` still give a ball straddling 0, and
  its sqrt is NaN (checked). Also fail loudly if `ba` is not finite. Robustness only; results are unaffected.
- RESULTS §2.1 F1: phrase it for piecewise-smooth (Lipschitz) pullback families of forms.
