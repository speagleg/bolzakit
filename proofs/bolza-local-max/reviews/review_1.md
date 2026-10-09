# Review 1 — item 23 (Bolza: local maximality of λ₁)

Reviewer: independent second reader (Rule R), 2026-10-08. I did not write the item.
Scope: README C1–C4, work/LITERATURE.md, PLAN.md, RESULTS.md, octagon.py and the scripts.
I wrote nothing else and committed nothing. My scripts are in the session scratchpad (`r23/chk1–4.py`).
Python: bolzakit venv (numpy 2.5.3, scipy 1.18.1), run with `-I`, `nice -n 10`, `OMP_NUM_THREADS=2`.

## Verdict: **ACCEPT WITH FIXES**

- The local-max criterion (RESULTS §2 Lemma) is **valid** under its stated hypotheses. The fixes below make it
  rigorous; none of them changes the conclusion.
- The construction is valid.
- The numerics reproduce at finer meshes.
- The labels are honest. Nothing is over-claimed as PROVED or CERTIFIED.
- One item-level explanation (Observation 2's symmetry reason) is wrong.
- Several reported numbers are coordinate artefacts and should be labelled as such.
- A representation-theoretic reduction (§1.4 below) strengthens C1 and C3 considerably, and should be added.

## 1. The local-max criterion (crux)

### 1.1 Setting
- The λ₁-cluster is isolated: the next eigenvalue is ≈ 5.35.
- The total projection P(v) is given by a Riesz integral, so it is jointly C^ω, or at least C², in all 6 real
  parameters. This holds whenever the pulled-back quadratic forms depend C^ω/C² on v in L(H¹, H⁻¹), and the mass
  form likewise.
- With a v-dependent L²(g_v)-orthonormal frame of Ran P(v), the cluster eigenvalues are the eigenvalues of a
  symmetric 3×3 matrix M(v) = λ₀I + A(v) + Q₂(v) + O(|v|³), with a uniform remainder on a neighbourhood.
- **C² is enough.** The analyticity reference (Bonifacio's [15], unchecked) is not needed: any C² family of
  pull-back diffeomorphisms works. Kato Ch. VII/VIII, or the multi-parameter Riesz projection argument, should be
  cited explicitly (Rule C).
- The non-smoothness of individual sorted eigenvalues is irrelevant. The argument uses only M(v) and the Weyl
  inequality.

### 1.2 Is tr A ≡ 0? **Yes, and it is essential.**
- Λ = tr M is C² and invariant under G = Isom(B). Since T_B^G = 0, dΛ(B) = 0, i.e. tr A ≡ 0.
- To be rigorous, T_B^G = 0 needs two facts:
  - (i) Bochner linearisation of a finite group acting smoothly at a fixed point, so that the tangent space of
    Fix(G) is T_B^G;
  - (ii) Fix(G) ≅ Teich of the (2,3,8) orbifold, which has dimension 0.

  Both should be cited.
- Numerically, tr A ≤ 1.1e-4 against ‖A‖ ≈ 80.
- tr A ≡ 0 is used **twice**:
  - (a) to get c = min_{|w|=1} (−λ_min A(w)) > 0: a nonzero traceless symmetric matrix has a negative eigenvalue;
  - (b) in step (ii), where |w| < C|k|². There tr A(w) = O(|k|²) would compete with tr Q₂(k) ~ −κ|k|². Without
    exact tracelessness the argument fails.

  The write-up should say this explicitly.

### 1.3 Proof sketch: checked
Write v = w + k, with k ∈ K = ker A and w ∈ W, a complement of K.
- **Step (i): |w| ≥ C|k|².** By Weyl, λ_min(M) − λ₀ ≤ λ_min A(w) + ‖Q₂(v)‖ + O(|v|³) ≤ −c|w| + C′(|w|² + |k|²).
  Since |k|² ≤ |w|/C, choosing C > 2C′/c gives ≤ −c|w|/2 + C′|w|² < 0 for 0 < |w| small.
- **Step (ii): |w| < C|k|².** Then |v| ≍ |k| and tr Q₂(w+k) = tr Q₂(k) + O(|k||w| + |w|²) = tr Q₂(k) + O(|k|³). So
  tr(M − λ₀) ≤ −κ|k|² + O(|k|³) < 0, and λ_min ≤ tr/3 < 0.
- All O-terms are uniform: Taylor's theorem with a C³, or even C², remainder bound on a fixed ball.
- The mixed-direction worry ("small K-component plus tiny off-K component") is exactly case (ii). It is handled by
  the trace bound, without needing the eigenvalues of A(w) + Q₂(k).

**Is the Hessian of the SUM the right object?**
- Yes, for this sufficient criterion: ½D²Λ(k,k) = tr Q₂(k), which is basis-free.
- The sharper object is λ_min of the second-order reduced matrix on K. Along K, M(tk) = λ₀ + t²Q₂(k) + O(t³), so a
  strict maximum along K needs only λ_min Q₂(k) < 0.
- So the criterion is **sufficient, not necessary**. State this.

### 1.4 Representation theory: new; should be added, and it upgrades C1/C3
I computed the action of the octagon's D₈ (the rotation r by π/4 and a reflection) on T_B, modulo trivial
directions (`chk1.py`; the tangent space is invariant to 2e-9):

| element | trace on T_B | on K | on W |
|---|---|---|---|
| r | −2.000000 | −1.000000 | −1.000000 |
| r² | −2.000000 | −1.000000 | −1.000000 |
| r⁴ (hyperelliptic) | +6 | | |
| reflection | 0 | **−1** | **+1** |

- r maps to a 4-cycle in PGL(2,3) ≅ S₄, and r² to a double transposition. The hyperelliptic involution acts
  trivially. Among S₄-representations of dimension 6, the only one with χ(4-cycle) = −2 and χ((13)(24)) = −2 is
  **V ⊕ V**, where V is the standard 3-dimensional representation.
- The 3-dimensional irreps of GL(2,3) factor through S₄, so E (dim 3) is V or V⊗sgn. In either case
  Sym₀E ≅ 2 ⊕ V.
- So A ∈ Hom_G(V ⊕ V, 2 ⊕ V), and **rank A ∈ {0, 3}**. Since the numerical A ≠ 0 (s₁ ≈ 81.6), rank A = 3 and
  dim K = 3 is forced. C1 can then be stated as "proved modulo A ≠ 0", which needs only a crude enclosure.
- The reflection separates K from W (traces −1 and +1). So K is an irreducible representation of the full isometry
  group, and **D²Λ|_K = c·g_inv** for the invariant metric. Negative definiteness is therefore automatic once
  c ≠ 0. The whole hypothesis of the Lemma reduces to **the sign of one number**, which F_B (⊂ K) already gives.
- **Fix:** the eigenvalues −14.29, −14.29, −7.29 are artefacts of the Euclidean vertex-coordinate metric, which is
  only D₈-invariant. The same holds for:
  - "q between 7.3 and 24";
  - the "slowest −7.3ε²" direction;
  - the K-sphere max and min.

  Only the sign (and λ_min Q₂ in an invariant metric) is intrinsic. RESULTS §0/§5d and README C3 should say this.
- Consistency check: the full 6×6 D²Λ has signature (3−, 3+). This matches c_K < 0 and c_W > 0, and is consistent
  with the cluster mean increasing along F_A ⊂ W.
- Still to do: an explicit check of the 3-cycle character. My decomposition is deduced from D₈ characters plus the
  S₄ character table, not computed directly.

## 2. Construction validity: confirmed
- The trace-coordinate Jacobian on the 11-dimensional constraint tangent has rank 6 (sv 6.8e4 … 6.4e2 | ≤ 1.1e-5).
- Independent check with 40 random words of length 2–5:
  - sv 1.3e7 … 5.0e4, then a gap to ≤ 7.9e-3;
  - restricted to the 6-dimensional slice E alone, rank 6 (min sv 5.0e4).

  Since trace functions factor through Teich₂, the slice → Teich₂ map is a local diffeomorphism at B. There are no
  degenerate or double-counted directions.
- Random deformed surface (|v| = 0.01 in span E):
  - constraints 7e-15;
  - Gauss–Bonnet area / 4π = 1.0000000000000007;
  - pairings in SO⁺(2,1) to 5e-13;
  - single vertex cycle 0→5→2→7→4→1→6→3→0, product − I = 3e-11;
  - glued mesh χ = −2;
  - FEM mass / 4π = 0.99999997 (n = 80).
- Code reading:
  - side nodes are hyperbolic-arclength uniform and identified j ↔ n−j, so they match exactly under the isometric
    gluing;
  - the Klein-model stiffness (1−r²)^{-1/2}(I − xxᵀ) and density (1−r²)^{-3/2} are correct;
  - P1 traces on glued edges differ by a projective reparametrisation (O(h²) nonconformity). This is consistent,
    and the Bolza validation covers it.

## 3. Numerics: spot-checks reproduce
- **Hessian on K, finer mesh** (`chk3.py`, t = 0.002/0.001 Richardson, along eigenvectors of the n = 80 H_K):

  | direction | n = 80 (item) | n = 160 | n = 320 |
  |---|---|---|---|
  | dir 0 | −14.290 | −14.337 | |
  | dir 2 | −7.290 | −7.338 | −7.350 |

  The h² trend from the item's n = 20/40/80 extrapolates to −7.35, consistent. The 0.05 non-monotone shift in
  dir 0 is plausibly O(h²) misalignment of K picked up by the large W-block (2651). It is irrelevant to the sign.
- **F_B p = 0.02**, R2 from n = 80/160/320 with an off-centre fan: λ₁ = 3.8304523, λ₂ = λ₃ = 3.8416925. This is
  identical to the item. (Weak independence: shifting the fan centre changes values by ≤ 3e-5 and does not remove the
  2+1 mesh splitting.)
- **Mixed directions:** worst K direction (refine_k) at ε = 0.03, plus Cε²·w for random unit w ∈ W, C ∈ {0.1, 0.5,
  2}, 4 each, R(80,160).
  - λ₁ − λ₁(B) ∈ [−7.3e-2, −7.9e-3], against −6.67e-3 for the pure K direction.
  - No mixed direction beats the pure-K decrease, as the Lemma predicts.
- n = 20–80 stability: adequate for the **sign**, given the margin (≥ 7 against a discretisation drift of ≈ 0.05).
  With §1.4 only one sign is needed, and that is robust.
- **Mesh-splitting pitfall recurs (Low/Medium).**
  - At the refine_k worst direction the physical λ₁–λ₂ gap is 1.65e-5 (ε = 0.015) and 4.4e-4 (ε = 0.03). Both are
    below the n = 80 mesh splitting (≈ 8e-4).
  - So the sorted-member R(80,160) values behind "−7.24/−7.41/−7.33" violate the item's own guardrail (AGENT.md).
    k_sphere near its maximum is similar.
  - The possible bias is ≲ 3e-4, against a signal of 1.6e-3, so the sign is safe but the 2-digit /ε² values are not.
  - Fix: use the lower-pair mean (≥ λ₁), or n = 160/320, and soften the precision.

## 4. Literature: quotes confirmed from the PDFs (pdftotext)
- **KMP arXiv:2111.12716, Thm 1.1(2):** "Any hyperbolic orbifold of genus 2 must satisfy λ1 < 3.8388976481. The
  quoted bounds are rigorously valid with the decimal expansion interpreted as an exact rational number." Correct.
  Gap 1.0389e-5.
- **FBP arXiv:2111.14699, "Theorem 7.1. The Bolza surface satisfies m1(B) = 3."** Correct. Proved via Cook's lower
  bound plus LP (verify8.sage).
- **SU13 §5:** "rigorous error bound of 10−6" and "supported by our calculations in Teichmüller space (the
  corresponding analysis will be discussed elsewhere)". Correct.
  - Nit: a quasimode bound locates *an* eigenvalue within 1e-6. Calling it "λ₁" uses other information.
- **Str16 arXiv:1604.02722, p. 22–23:** confirmed verbatim, including the preceding sentence "The reason for it being
  locally maximal is however its degeneracy." **Fix:** quote this sentence too. Str16 explicitly asserts local
  maximality, numerically, and the README novelty text should acknowledge the prior (unproved) assertion.
- **New searches:**
  - web: "Bolza surface local maximum first eigenvalue Laplacian Teichmüller space"; "critical points first
    eigenvalue hyperbolic surfaces moduli space λ1 local maximum genus 2 Bolza 2025";
  - arXiv API: all:Bolza AND all:eigenvalue; abs:"local maximum" AND eigenvalue AND hyperbolic.

  Nothing beyond what LITERATURE.md lists: only the extremal-length systole (2105.03871), the area-normalised problem
  (1704.06384) and the spin bootstrap. I agree that local maximality of λ₁ at B is **not proved** in the literature
  found.

## 5. Labels and honesty
- NUMERICAL and "argument (unreviewed)" labels are appropriate. Novelty UNCHECKED is appropriate.
- The Strohmaier contradiction is handled carefully but slightly too strongly.
  - His remark is correct for every direction outside K, an open dense set.
  - It is contradicted only for directions in K, e.g. F_B, where the analytic branch λ₁(p) has a strict quadratic
    maximum.
  - Suggest: "contradicted as literally stated ('no matter what perturbation'); true for generic directions", rather
    than "treat it as an error".

## Issues (severity)

| # | Severity | Issue | Fix |
|---|---|---|---|
| 1 | Medium | Lemma correct, but the sketch omits: the multi-parameter Kato/Riesz reduction with varying inner product (cite); the rigorous T_B^G = 0 (Bochner + rigidity of the (2,3,8) orbifold); and that tr A ≡ 0 is used in step (ii), not only for c > 0 | Write out the proof as in §1.3 |
| 2 | Medium | The K-Hessian eigenvalues (−14.29, −14.29, −7.29), "q ∈ [7.3, 24]" and the "slowest −7.3ε²" are coordinate-dependent; by irreducibility D²Λ\|_K = c·g_inv | Relabel as coordinate-dependent; add the §1.4 reduction (rank A ∈ {0,3}; sign of one scalar) |
| 3 | Low–Med | refine_k / k_sphere violate the sorted-member Richardson guardrail (λ₁–λ₂ gap below the mesh splitting) | Pair-mean or finer meshes; reduce the stated precision |
| 4 | Low | Observation 2 is wrong: the visible symmetry of F_B (Z/4 ⋊ reflection), mod the hyperelliptic involution (trivial on E), is abelian of order 4 and has only 1-dimensional real irreps, so it cannot force λ₂ ≡ λ₃. Either a hidden isometry is involved (e.g. the extra involution of the y² = x(x⁴+ax²+1) family), or the reason is unexplained | Replace the symmetry explanation. Stationarity along F_B follows from F_B ⊂ K (W-part 1.3e-8, checked) |
| 5 | Low | Str16's "locally maximal" sentence is not quoted; the contradiction wording is too strong | §4/§5 |
| 6 | Low | hessian_sum.py docstring says a negative definite Hessian of Λ implies the maximum, but D²Λ is indefinite on T_B | Fix the docstring |
| 7 | Low | All numerics come from one discretisation family; the off-centre fan is only weakly independent | Already listed as next step 3 |
| 8 | Nit | The SU13 1e-6 bound locates an eigenvalue; identifying it as λ₁ (and the multiplicity) uses FBP | Wording |

## Checks run
- `chk1.py`:
  - rank of the trace Jacobian with the item's words and with 40 independent random words, on the 11-dimensional
    tangent and on E;
  - D₈ characters on T_B, K and W;
  - H_K eigenvalues and invariance under r.
- `chk2.py`:
  - full validation of a random deformed surface;
  - F_B tangent lies in E ⊕ Z with W-component 1.3e-8 (F_B ⊂ K).
- `chk3.py`:
  - 12 mixed K+W directions;
  - D²Λ along K at n = 160 (2 directions) and n = 320 (1 direction);
  - F_B p = 0.02 at n = 80/160/320 with an off-centre mesh;
  - Bolza with an off-centre mesh.
- `chk4.py`: sensitivity to the fan centre (≤ 3e-5).
- Read KMP Thm 1.1, FBP Thm 7.1, the Str16 passage and the SU13 §5 quotes from the arXiv PDFs. Ran 2 web searches
  and 2 arXiv API queries.
- Total compute ≈ 10 min.
