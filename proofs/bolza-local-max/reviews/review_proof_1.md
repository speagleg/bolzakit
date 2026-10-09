# review_proof_1: final adversarial read of PROOF.md (OPEN-1, CERT_G3, whole chain)

Reader: an independent agent that wrote none of item 23. Date 2026-10-08.
Scratch: `scratch_RF/` (re-run output, Felikson–Natanzon text).

**Overall verdict: ACCEPT-WITH-FIXES.** No logical hole was found. OPEN-1's proposed closure is correct and
complete, given inputs that are standard and that I list below. CERT_G3 re-runs bit-identically, and its §3
argument is sound. All the required fixes are wording or status changes. Once they are applied, the headline may
read **PROVED (computer-assisted), conditional on the cited inputs of §Dependencies (weakest: SU13's lower
bound)**.

## Verdict table

| # | Item | Verdict |
|---|---|---|
| 1 | OPEN-1 closure, step 1 (g_S analytic; Hessian invariant at a critical point) | CORRECT |
| 2 | OPEN-1 closure, step 2 (σ∘dg_S = g·σ) | CORRECT |
| 3 | OPEN-1 closure, step 3 (Lemma 11 for all of G, so H\|K₀ = c·g_inv) | CORRECT |
| 4 | OPEN-1 inputs: are they all established? | YES, plus two standard facts that should be listed (F1, F2 below) |
| 5 | cert_dpi_rank.py re-run (`nice -n 19`, work/.venv) | REPRODUCED: stdout identical apart from the runtime line; JSON identical apart from `runtime_s` |
| 6 | The τ_w are functions on T₂ (conjugation-invariant, marking-consistent) | CORRECT |
| 7 | Rank 6 ⇒ slice S, homeomorphism π\|_S, retraction p, σ iso on T S | CORRECT |
| 8 | Felikson–Natanzon quotes (arXiv:1102.4861v1) | VERIFIED verbatim (pdftotext layout only) |
| 9 | Lemma 10 / the final local-max mechanism | VALID |
| 10 | Hypotheses of every lemma discharged | YES |
| 11 | Numbers in PROOF.md vs .out/.json files (spot-check) | MATCH |
| 12 | Labels honest | YES, once the OPEN-1 status lines are updated |
| 13 | Notation | FIX: λ₀ is used for two different things |

## 1. OPEN-1

**The setting.** Λ_S is real-analytic on the 6-dimensional chart S (Lemma 2 on S), and σ: T_{x_B}S → TT is a
linear isomorphism (Lemma 9d2). The closure works entirely in the chart S. It never needs a smooth structure on
T₂ or the analyticity of π into T₂, so the alternative closure (D14) is not needed. This is the route suggested in
my brief ("Λ descends to T₂ … Hessian transforms tensorially"), carried out in the chart that the certificate
actually provides. The fibre question is handled by p: Λ_C = Λ_S∘p, dp|_{TS} = id and dΛ_S(x_B) = 0 give
D²Λ_C(w,w) = D²Λ_S(dp w, dp w). So the Hessian on C is the pull-back of the Hessian on S, and it is degenerate
exactly along ker dp, the fibre directions.

**Step 1.**
- g acts on marked surfaces by precomposing the marking with g⁻¹. Isometric surfaces have equal spectra, so
  Λ∘g = Λ on T₂, and g·B = B.
- For each of the six words w, τ̄_w∘g = τ̄_{g⁻¹w}. I checked the direction: g⁻¹ is an isometry
  ((g⁻¹)*k) → k, so ℓ_{(g⁻¹)*k}(c) = ℓ_k(g⁻¹c).
- g⁻¹w is a free homotopy class. For an orientation-reversing g it is defined only up to an outer automorphism and
  up to inversion. This does no harm: tr_{SO(2,1)} is conjugation-invariant, and tr(A⁻¹) = tr(JAᵀJ) = tr A. Pick
  any word representative of g⁻¹w in the g_k(x_B). Then τ_{g⁻¹w} is an analytic function on C.
- Hence, on a neighbourhood of x_B in S, g_S = (Φ|_S)⁻¹∘(τ_{g⁻¹w}|_S)_w is analytic and fixes x_B.
- This needs g·π(x) ∈ N for x near x_B, so that τ̄ can be inverted there. That holds by continuity of the
  mapping-class action on T₂ (**F1**, standard, not stated).
- Λ_S∘g_S = Λ_S, and dΛ_S(x_B) = 0. Hence D²Λ_S(dg_S·, dg_S·) = D²Λ_S. ✓

**Step 2.**
- Apply Lemma 3 Step 4 to the analytic curve y = g_S∘x ⊂ S ⊂ C. That step needs only a C¹ curve in C.
- Lemma 9d2 holds for *every* word, not just the six, because its proof is word-independent. Hence
  dτ_w(σ(y′)) = dτ_{g⁻¹w}(σ(x′)).
- Next, τ_{g⁻¹w}(k) = τ_w((g⁻¹)*k) on hyperbolic metrics. Differentiating at h₀ = (g⁻¹)*h₀ along k_t gives
  dτ_w((g⁻¹)*σ) + dτ_w((g⁻¹)*L_Y h₀).
- (g⁻¹)*L_Y h₀ = L_{g_*Y}h₀, which is killed by Diff₀-invariance.
- (dτ_w)_{6 words} is injective on TT, so σ(y′) = g·σ(x′). This uses the same convention g·σ = (g⁻¹)*σ as Lemma
  3 Step 5 and Lemma 6. ✓

**Step 3.** H(g·σ(v)) = D²Λ_S(dg_S v, dg_S v) = H(σ(v)) for all g ∈ G, including t. K₀ is absolutely
irreducible (norm 1, EXACT), so H|_{K₀} = c·g_inv. Then c·g_inv(w_B,w_B) = D²Λ_S(dp v_B, dp v_B) ≤ −12.8504, and
w_B ≠ 0, so c < 0. ✓

**Facts needed, and their status.**

| fact | status |
|---|---|
| Λ_S is C² (analytic) on S; dΛ_S(x_B) = 0 | Lemma 2 on S + Lemma 4 + Lemma 9d1: established |
| σ: T S → TT is an isomorphism; (dτ_w) is injective on TT; dim TT = 6 | certificate (rank 6, re-run) + Lemma 9d2 + D6: established |
| d(τ∘x)(v) = dτ(σ(v)) for any word, any C¹ curve in C | Lemma 9d2, with Lemma 3 Step 4 (k_t): established, modulo F2 |
| **F1**: the extended mapping-class group acts on T₂ by homeomorphisms; g fixes B | standard; **add to D12** |
| **F2**: for each loop class, τ_w is a C¹ function on smooth metrics near h₀, with a linear derivative, and is Diff₀-invariant on hyperbolic metrics. CERT_G3 says "holonomy depends smoothly on the metric". A one-line justification suffices: take the holonomy of the so(2,1)-valued Cartan connection ω_LC + solder form, which is defined for every metric, depends smoothly (C¹ in the C¹ topology) on k along a fixed loop, and is flat iff K = −1 | standard; **state it in Lemma 9d2** |
| naturality of length under isometries | trivial |

None of these is new mathematics. So OPEN-1 is **closed**.

## 2. CERT_G3, second read

- **Re-run.** `nice -n 19 .venv/bin/python cert_dpi_rank.py` gave exit code 0. The output is in
  `scratch_RF/cert_dpi_rank.rerun.out`.
  - Every line is identical to `cert_dpi_rank.out` except the runtime (0.1 s against 0.0 s). Key values:
    det = [7872428136359854080.0 ± 3e-6], 392/462 minors certified nonzero, the 11×11 cross-check
    −1.3827971011809151404e+24 ± 2.74e+4, and `CERTIFIED ... True`.
  - The script rewrites `cert_dpi_rank.json`. I saved the original first (`scratch_RF/orig_cert_dpi_rank.json`)
    and restored it byte-for-byte afterwards. The re-run JSON was identical apart from `runtime_s`.
  - `git status` shows only `scratch_RF/` as new.
- **Code read.**
  - `pairings` builds g_k = F(V_{k+1},V_k)·F(V_{k+4},V_{k+5})⁻¹, the convention of octagon.py.
  - The residuals certify g_k V_{k+4} = V_{k+1} and g_k V_{k+5} = V_k, and that g_k ∈ O(2,1). det = +1 holds by
    the frame sign fix, and time orientation holds because tr g_k ≈ 22.3 > 0.
  - The traces are taken of genuine words in the generators. Dual-number derivatives are exact. The kernel basis
    Tb = [−Jc_piv⁻¹Jc_free; I] is a correct parametrisation of ker Jc.
  - The columns of the best minor are Tb-columns for free coordinates {3,6,10,12,13,15}. These are exactly the
    tangent vectors of the coordinate slice S, so rank(JT|_{T S}) = 6 is what is certified.
- **Functions on T₂.**
  - The marking g_k(x_B) ↦ g_k(x) is a group isomorphism Γ_{x_B} → Γ_x, because the single vertex-cycle relation
    holds on C (Poincaré).
  - Points equal in T₂ have holonomies conjugate in PSL₂ℝ ≅ SO⁺(2,1), compatibly with the marking. The
    SO(2,1)-trace needs no sign lift.
  - The marking used by Lemma 3's Ψ_t is the same one, because Ψ̃_t is ρ_t-equivariant. ✓
- **Topology.**
  - π|_S is injective because τ̄∘π|_S = Φ|_S is.
  - Invariance of domain, from a 6-ball into a topological 6-manifold, gives a homeomorphism onto an open N.
  - p = (Φ|_S)⁻¹∘Φ satisfies π∘p = π, using τ̄-injectivity on N and π(x) ∈ N for x near x_B.
  - dp|_{TS} = id, so dΛ_S(x_B) = 0 and D²Λ_C = p*D²Λ_S.
  - σ∘dp = σ holds by injectivity of dτ on TT. ✓
- **Felikson–Natanzon** (PDF fetched; text in `scratch_RF/felikson_natanzon_1102.4861v1.txt`).
  - "So, the Fenchel-Nielsen coordinates on T build from the pants decomposition P consist of 3g − 3 + 2n length
    parameters l(ci ) … and 3g − 3 + n angle parameters" ✓.
  - "The construction establishes the homeomorphism between T and R^{3g−3+2n}_{>0} × R^{3g−3+n}" ✓. The PDF
    layout splits the sub- and superscripts; the content matches.
  - The intro quotes on Wolpert [13]/[14], Schmutz [10] and "9g − 9" are verbatim ✓.
  - For g = 2, n = 0 this gives T₂ ≅ ℝ³_{>0} × ℝ³, which is what §3(c) uses.
- **Minor.** CERT_G3 §3(d)1 uses σ for the section (π|_S)⁻¹, which clashes with Lemma 3's σ. PROOF.md avoids
  the clash. Lemma 5 also uses σ(t) for a radius correction. These are cosmetic.

## 3. Whole-document read

**Chain.**
- Lemma 1 → (H4) of Lemma 2 → Lemma 2 on C and on S.
- Lemma 3 → A intrinsic and equivariant → Lemma 4 (dΛ = 0, tr A ≡ 0).
- Lemma 5 + Lemma 6 → K = K₀, rank 3.
- Lemma 7 → w_B ∈ K₀∖0.
- Lemma 8 → D²Λ_C(v_B,v_B) ≤ −12.85.
- Lemma 9 → chart S, the transfer, σ an isomorphism.
- Lemma 11 (OPEN-1, now closed) → c < 0 → κ > 0.
- Lemma 10 → λ_min M_S < λ₁(B) → π|_S is a homeomorphism → the theorem.

Each hypothesis is discharged by an earlier lemma or a listed dependency. Lemma 1a's run A counts (0, 3.3]
(Lemma A subtracts F(0)), so a bound below 1 correctly excludes eigenvalues there.

**Final mechanism (Lemma 10).** The argument does *not* claim that Λ has a local maximum: off K, Λ can increase at
second order. The mechanism is a case split.
- Case (i), |w| ≥ C|k|²: Weyl's inequality gives λ_min ≤ λ_* − c_A|w|/2 + o, where c_A > 0 because A(w) is
  traceless and nonzero.
- Case (ii), |w| < C|k|²: λ_min ≤ (tr M)/3 = λ_* + tr Q₂(k)/3 + o(|k|²) < λ_*. Here tr A(w) = 0, and the cross
  terms of Q₂ are O(|k|³).
- Then λ₁(π x) = λ_min M_S(x), because M_S's eigenvalues are exactly the cluster λ₁,λ₂,λ₃ (Lemma 2.1) and λ₁ is
  the smallest.
- I checked the constants: with C > 8C′/c_A, the terms 2C′|k|² ≤ 2C′|w|/C < c_A|w|/4. Valid.

**Numbers spot-checked against the output files.** All match the files:
- lp_count.out: c = 0.9816424759937855 and 0.9964283521820437; bounds 0.912612352675290375 and
  3.86966111976973759; sanity windows 3.4803 and 7.5466.
- cert_rank.out: det 46151927.03406868578682884302258; −5.5729711960.
- cert_rr_bolza.out: n = 192, 3.839184804702046, area ±2.23e-13.
- cert_rr_p2_bolza.out: n = 96, 3.8388872829261946, 11.516661838577997, ±1.19e-18.
- cert_first_order.out: −3.0875694043934345, μ, b(u,u)/m, 0.0012, 0.0049, 1.71, β 32.99.
- cert_sign.out: −12.850476214196528, central [−14.2742 ± 2.23e-5], n = 48 value −12.566084…; T₂ 45.6357, coupling
  −59.9099.
- cert_dpi_rank.out: as in §2.
- cert_chars_A.json: dim_T_B_eps_s_line = 1; T_B χ_r = −2, χ_J = 6, χ_t = 0, norm 2, reversing all 0.
- Arithmetic checks: −√(2+2√2) = −2.197368…; D₈ norm (9+9+4+2+4+4)/16 = 2; gap 4.0 − 3.8392 = 0.1608.

## 4. Required fixes (exact text changes to PROOF.md)

1. **Status block (lines 9–15).** Replace the three bullets with:
   "Every component lemma has been second-read; Lemma 9 (CERT_G3) and the OPEN-1 closure were second-read in
   review_proof_1 (ACCEPT-WITH-FIXES). The theorem is PROVED (computer-assisted), conditional on the external
   inputs of §Dependencies; the weakest is the Strohmaier–Uski λ₁ lower bound."
2. **§0, "Conditional on".** Delete item (b) "the closure of OPEN-1".
3. **Lemma 9 header.** Replace "**not independently second-read**, see §Reviews" with "second-read: review_proof_1
   §2, re-run reproduced".
4. **Lemma 9d2.** Append: "(Smoothness of τ_w: τ_w(k) is the trace of the holonomy of the flat so(2,1) Cartan
   connection ω_LC + θ of k, which is defined and C¹ in k for every metric along a fixed loop, and is flat,
   hence free-homotopy-invariant, exactly when K = −1. The identity holds for every word w, not only the six.)"
5. **Lemma 11 header.** Replace "**See OPEN-1.** Proposed proof in §OPEN; not second-read" with "PROVED (§OPEN-1
   closure; second-read: review_proof_1 §1)".
6. **Assembly.** In step 4, delete "**[OPEN-1]**". In step 6, delete "(conditional on OPEN-1)".
7. **§OPEN.** Retitle it "OPEN-1 (CLOSED)", and add after "Inputs of the proposed closure": "• continuity of the
   extended mapping-class action on T₂ (standard), used so that g·π(x) ∈ N for x near x_B." Change "NOT
   second-read" to "second-read in review_proof_1: CORRECT".
8. **Notation clash.** §1 sets λ₀ = 0, and Lemma 8 step 1 uses λ₀ ≡ 0. Lemmas 2, 3 and 10 and Assembly step 5
   use λ₀ for λ₁(B). Rename the latter to λ_* := λ₁(B) throughout Lemmas 2, 3, 10 and the Assembly. This is
   cosmetic, but Assembly step 5 currently reads "λ₀ = λ₁(B)", contradicting §1.
9. **Dependencies.** D12: add "continuity of the extended mapping-class group action on T₂". D14: mark it "not
   needed (OPEN-1 closed in the chart S)".
10. **Reviews table.** Replace the "(none) | CERT_G3" row with "review_proof_1 | CERT_G3, OPEN-1, whole document |
    ACCEPT-WITH-FIXES | fixes 1–9 above". Replace "Integrated-proof read: pending" with "Integrated-proof read:
    review_proof_1."
11. **Optional.** In §0, "Off K, λ₁ decreases linearly" → "Off K, λ₁ − λ₁(B) ≤ −c_A|w|/2 to first order (Lemma
    10 case (i))".
