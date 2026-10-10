# paper/ — drafting notes

Files: `main.tex` (amsart), `refs.bib`, compiled `main.pdf` (27 pages, rebuilt 2026-10-10).
Build: `pdflatex main && bibtex main && pdflatex main && pdflatex main`. It compiles with no errors, no undefined
references and no overfull boxes. BibTeX gives no warnings.

## 1. Places where the exposition needed more than PROOF.md

| Where in main.tex | What was added | Source used / status |
|---|---|---|
| Thm 3.1 (cluster reduction), proof | Steps 1–6 written out. PROOF.md Lemma 2 lists only the ingredients. | F1_PROOF.md §2, condensed. Content is unchanged. |
| Lemma 3.3 (Lemma N) | Lemma N is stated and proved | PROOF.md Lemma 3, nearly verbatim. I renamed the isometry φ_t to ι_t to avoid a clash with quadratic differentials φ. |
| §4.1 algebraic model: P = (det/μ)[[d,b],[c,a]], the y-factor 2+2i for t, χ_T = 2 Re tr Sym²P | Explicit formulas | CERT_P6.md §1 and the cert_chars.py docstring |
| §4.2 characters of K₀, W₀ on r, s_mid, s_vert, t (Table 10) | Values for the two summands | `cert_chars_A.out` (X1/X2 blocks). The identification K₀ = X1 comes from CERT_P6 §2.1 / `v_B_line_lies_in: X1` and LEAN.md. The entries χ(r²) and χ(J) on K₀ and W₀ were left blank because the output does not list them, though they follow trivially. |
| Lemma 4.7 proof, "r_in is the inradius of the regular octagon" | Interpretation of PROOF.md's `sinh r_in` | My reading; PROOF.md does not define r_in. **Check.** |
| §2.2 frame formula: n ∝ 𝖡(P × e) | Written with the Minkowski matrix | CERT_G3 §2 writes "n ∝ J(P × e)", where J is the Lorentz matrix (not the hyperelliptic J) |
| Lemma 6.1 (local-max criterion) constants | A generic C″ with C > 4C″/c_A | PROOF.md uses 2C′(\|w\|²+\|k\|²) and C > 8C′/c_A, which differ only in the constants. Worth a glance. |
| Lemma 6.3(b): the smoothness of τ_w via the Cartan connection | Restated | PROOF.md Lemma 9d2 (review_proof_1 fix 4) |
| Lemma 6.6 (Ky Fan) proof | Reduction to min-max | F1_PROOF.md §5 |
| Lemma 6.7 (formula for g″(0)) proof sketch | One-sentence derivation | CERT_SIGN §1.2. The "3e-39 exact pencil test" is from PROOF.md Lemma 8. |
| §7.1 trace formula statement, FBP Lemma 3.2, Lemma A, test functions, Hermite convention, odd N | Full statements | LITERATURE.md §5.1 and CERT_RESULTS §1.1–1.3. The geodesic weight Λ(γ) is **not** defined in the paper; it is deferred to "the conventions of FBP". |
| §7.2 element integrals (P1: 4 series terms with the remainder formula; P2: 6 terms) | Details | CERT_RESULTS §2.2 and CERT_SIGN §1.4 |
| §7.3 Lemma 7.3 (gap transfer), proof | A 3-line proof from the spectral decomposition, plus the monotonicity remark (needs μ < Λ₄) | The **proof is mine**. PROOF.md and CERT_RESULTS state only the inequalities; review_sign_1 row 3a confirms them ("spectral decomposition … monotone in λ₄"). |
| §7.4 β formulas | β_T ≥ ‖N⁻¹Ṅ‖_F + \|Ḋ/D + Q̇/2Q\| and β_ρ ≥ \|Ḋ/D − 1.5Q̇/Q\| | CERT_RESULTS §4.4 |
| §7.5 corrector construction; the jet formula for Q^{−s} | Details | CERT_SIGN §1.3–1.4 |
| §7.6 "the four generator traces have rank only 3"; the singular values | Context | CERT_G3 §2 remark (the singular values are float) |
| §5 remark "equivariance under r alone suffices for tr A ≡ 0" | Added | LEAN.md finding 3 |
| §9 Lean | Scope and limits | LEAN.md §0–4. The version is written as "toolchain v4.34.1" (from `lean/lean-toolchain`). LEAN.md also says "Mathlib v4.34.1", which I did not repeat. |
| Appendix A tables | Full-precision values | Copied from `lp_count.out`, `cert_rr_bolza.out`, `cert_rr_p2_bolza.out`, `cert_first_order.out`, `cert_sign.out`, `cert_dpi_rank.out`, `cert_rank.out` and `cert_chars_A.out`. The n = 48 row of the A ≠ 0 table is from CERT_RESULTS §4.5. |
| Appendix B commands | Run commands and times | checks/README.md, CERT_SIGN §7 and the cert_chars.py `__main__`. The cert_chars commands (`python cert_chars.py`, `python cert_chars.py B 48 96`) are inferred from the code and were not taken from a doc. |
| MSC codes, keywords | Chosen by me | Not sourced. **Check.** |

## 2. Claims I could not fully source (or that are my own wording or speculation)

1. *(Resolved.)* An earlier draft of Remark 10.1 said "an independent lower bound is in progress". That bound is now
   certified (Lehmann–Goerisch: main.tex `prop:LG` in `subsec:LG`, and PROOF.md Lemma 1(d)), and main.tex no longer contains
   the sentence.
2. Remark 10.3: "A quantitative version … could in principle reduce the global problem to a compact region away
   from B." This is my speculation, marked "we have not attempted this". Delete it if unwanted.
3. Remark 10.4(2): the suggestion to apply the method to the Klein quartic is mine. It is motivated by FBP-LP's
   "conjectured maximizers in these genera" quote.
4. Related work: "for flat tori, local maximality … is classical (see [KLO, §3.3])". NOVELTY.md lists Berger 1973 /
   Nadirashvili 1996 / KLO. I cite only KLO (title unknown). The word "classical" is mine.
5. Remark 10.1 mentions "Lehmann–Goerisch or Liu-type methods" without citation. Liu–Oishi were *not read*
   (LITERATURE §5.3). Only Carstensen–Gedicke is cited.
6. Remark 10.2: the numerical slopes (−5.715|p| along F_A, −21.2p² along F_B) come from README §1 C4. They are
   NUMERICAL, and the paper labels them so. The stabiliser-of-w_B fact is EXACT (CERT_P6 row 8). The explanation of
   λ₂ ≡ λ₃ is unreviewed (CERT_P6 §2.6), and the paper does not claim it.
7. Acknowledgements: "re-ran the key certificates and re-derived the main formulas". This is supported by
   review_cert_1 (re-runs), review_sign_1 (re-derivation, independent P1 code) and review_proof_1 (bit-identical
   re-run of cert_dpi_rank). It is a summary, not a quote.
8. Standard facts are used without a citation: the Poincaré polygon theorem, the Riemannian centre of mass,
   invariance of domain, the finite S₄ ⊂ PSL(2,ℂ) being octahedral, and continuity of the mapping-class action.
   These are the same as PROOF.md D8, D11 and D12.

### Bibliography
- The former `% CHECK` items were resolved on 2026-10-10 against arXiv, Crossref and the publishers' pages: titles,
  first names, published versions (Str16 in Contemp. Math. 700; Hide–Petri in Compositio 162; KMP, Bonifacio
  2107.09674, FBMGVP, NS, FBR, KLO, JLNNP, Mazzoli, Daskalopoulos–Wentworth, El Soufi–Ilias, Duran), volumes, pages
  and the year of the Petri HDR (2024).
- **Not read by the project** (cited second-hand): Kato, Fan 1949/1950 (1949 seen as a page scan only), Jenni 1984,
  Buser, Tromba, and Chavel (via Cook Thm 1.22). Fischer–Marsden [FM75] is attributed in the text via Mazzoli and
  has no bib entry.
- No software citations (python-flint/Arb, sympy, Lean, Mathlib) were added, because none is recorded in
  LITERATURE.md or NOVELTY.md. A journal will probably want them.

## 3. Caveats from the reviews, and where they appear in the paper
- SU13 provenance (one sentence, unpublished code, absolute vs relative): Hypothesis 1.1, Remark 10.1 and §7.5
  sensitivity.
- The script λ-ball sits 4.2e-11 below SU13's ends, with the upper end covered by P2; the re-run with the exact
  lower end gives −12.850504: §7.3.
- β_max is enclosure overestimation, not geometry: §7.4 and §7.5.
- The error-budget rows are not additive: §7.5, Table 8 and §9.
- F1 is not needed for (Λ∘x)″ ≤ g″ itself: §6.4, right after the Ky Fan inequality (Λ∘x)″(0) ≤ g″(0).
- Cai is cited as an analogue only, and Kato as quoted second-hand via Duran: Remark after Thm 3.1.
- The s_mid vs s_vert reflection type is load-bearing: Remark after Lemma 4.7.
- The D₈ group is insufficient, so t is needed (OPEN-1, now closed): §4.4 and Lemma 6.4.
