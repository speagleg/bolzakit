# The Bolza surface is a strict local maximum of λ₁

**Theorem.** Let B be the Bolza surface and T₂ the Teichmüller space of genus 2. There is a neighbourhood N of B
in T₂ such that λ₁(X) < λ₁(B) for every X ∈ N ∖ {B}.

**Status** (copied verbatim from the Status row of math-pursuits item 23's README; there, `work/PROOF.md` is
[PROOF.md](PROOF.md) in this folder):

> **PROVED (computer-assisted), conditional on the cited inputs in work/PROOF.md §Dependencies: FBP Lemmas 7.2–7.3 and §7.4 geodesic data, the Selberg trace formula, the Lehmann–Goerisch theorem (quoted from secondary sources), and standard facts. No longer depends on Strohmaier–Uski (SU13): λ₁ ≥ 3.8388757010 is now certified in-house (review_lb_1: ACCEPT). Every component second-read; integrated proof read: ACCEPT-WITH-FIXES, fixes applied.**

Strohmaier–Uski (SU13) is kept only as an optional cross-check: every SU13-based number below has an SU13-free
twin, and the proof uses the SU13-free one.

**No human expert has read this proof yet.** The "second reads" and "reviews" in [reviews/](reviews/) were done by
AI reviewers, not by referees or specialists. Treat the result as a computer-assisted claim that is awaiting expert
review.

## Contents

| Path | What it is |
|---|---|
| [PROOF.md](PROOF.md) | The integrated proof: theorem, the lemma chain, the assembly, §Dependencies, review verdicts |
| [LB_RESULTS.md](LB_RESULTS.md), [LB_PLAN.md](LB_PLAN.md) | The SU13-free Lehmann–Goerisch lower bound λ₁ ≥ 3.8388757010 (Lemma 1(d)) |
| [NOVELTY.md](NOVELTY.md) | What is and is not new, against the literature |
| [LEAN.md](LEAN.md) | What the Lean 4 formalisation covers (the logical chain, conditional on the computational inputs) |
| [reviews/](reviews/) | The review record (`review_*.md`), one file per review round |
| [paper/](paper/) | The paper: `main.tex`, `refs.bib`, `NOTES.md`, and the built `main.pdf` |
| [certificates/](certificates/) | The golden `bolzakit-cert/1` JSON certificates, regenerated through the bolzakit API |
| [reproduce.py](reproduce.py) | Checks or regenerates the certificates and prints claim / certified / published / provenance / match |

The documents are copied from math-pursuits item 23 (`work/`, `paper/`) as of 2026-10-09 (math-pursuits commit
53bc014). They refer to item-23 script and output names (`lp_count.py`, `cert_sign_48_96.json`, …); the certificates here are their bolzakit
regenerations.

## The chain, in one screen

PROOF.md §2 has the details and §"Theorem: assembly" the assembly. In short:

1. **Lemma 1 — isolation and the λ₁ ball.** A Selberg trace-formula LP shows Spec(B) ∩ (0, 3.3] = ∅ and at most 3
   eigenvalues in [3.3, 4.0] (`lp_A_3.3`, `lp_B_3.3_4.0`). P1 and P2 Rayleigh–Ritz give λ₃ ≤ 3.8391848048 and
   λ₁ ≤ 3.8388872830 (`rr_p1_n192`, `rr_p2_n96`), so λ₁ = λ₂ = λ₃ (with FBP Lemmas 7.2–7.3) and λ₄ > 4.0. A
   Lehmann–Goerisch bound gives λ₁ ≥ 3.8388757010 (LB_RESULTS.md).
2. **Lemmas 2–4.** The cluster λ₁, λ₂, λ₃ is analytic; its first-order matrix A is intrinsic, G-equivariant and
   traceless (dΛ = 0 at B).
3. **Lemma 5 — A ≠ 0.** A certified first-order test along F_A, P1 n = 192: b(e,e)/m ≤ −3.0383475 with the
   Lehmann–Goerisch ball (`first_order_FA_n192_lg`; −3.0875 with the SU13 ball, `first_order_FA_n192_su13`).
4. **Lemmas 6–7 — representation theory (exact).** ker A =: K is 3-dimensional and absolutely irreducible, and the
   equiangular family F_B is tangent to it.
5. **Lemma 8 — the Ky Fan sign bound.** Along F_B, P2 n = 96: (Λ∘x)″(0) ≤ −9.3411 with the Lehmann–Goerisch ball
   (`kyfan_FB_n96_lg`; −12.8504 with the SU13 ball, `kyfan_FB_n96_su13`).
6. **Lemma 9 — the chart.** The trace-coordinate map has rank 6 at B (`rank`), giving a local chart of T₂.
7. **Lemmas 10–11 and assembly.** By G-invariance the Hessian of Λ on K is c·g_inv, and step 5 forces c < 0.
   The local-max criterion (Lemma 10) then gives λ₁ < λ₁(B) near B.

**Dependencies.** The complete list of external inputs (D1–D15, with where each is used and how weak it is) is
PROOF.md [§Dependencies](PROOF.md#dependencies-complete-list-of-external-inputs). The load-bearing ones are FBP
Lemmas 7.2–7.3 and §7.4 geodesic data, the Selberg trace formula in FBP's normalisation, the Lehmann–Goerisch
theorem (quoted from secondary sources), and the trusted software (python-flint/Arb, sympy). SU13 (D1) is
optional.

## Reproduce

From the repository root, with `pip install -e ".[cert]"`:

```bash
python proofs/bolza-local-max/reproduce.py --quick     # a few minutes
python proofs/bolza-local-max/reproduce.py --full      # about 30–45 min on one core; rewrites certificates/
python proofs/bolza-local-max/reproduce.py --full --out /tmp/certs --only rank lp_A_3.3
python -m bolzakit_check proofs/bolza-local-max/certificates/*.json
```

- `--quick` runs the independent checker `bolzakit_check` on every stored certificate, reads the certified value
  out of each, and then re-runs the n ≤ 24 analogues (P1 and P2 Rayleigh–Ritz at n = 24, the n = 24 first-order
  test, which does *not* certify, as published, and the n = 24 Ky Fan test).
- `--full` regenerates the golden set through the bolzakit API.
- Both print a table of claim | certified value | published value | provenance | match, and exit 0 only if every
  row matches. The provenance column names the item-23 file a value is published in, or marks it as a bolzakit
  regression pin (two rows: the Ky Fan LG n = 96 value and the n = 24 first-order value).
- `bolzakit_check` exits 0 if every file is unconditionally VALID, 2 if none is INVALID but some is only
  conditionally VALID, and 1 if any file is INVALID. On the golden set it exits 2: `lp_B_3.3_4.0` is a mode-B
  count, which is valid *given* λ₁ ≥ 3.3, and that hypothesis is discharged by `lp_A_3.3`.

| Certificate | Claim | Certified value |
|---|---|---|
| `lp_A_3.3.json` | trace-formula count in (0, 3.3] < 1 | 0.912612… |
| `lp_B_3.3_4.0.json` | count in [3.3, 4.0] < 4, given λ₁ ≥ 3.3 | 3.869661… |
| `rr_p1_n192.json` | λ₃ ≤ (P1, n = 192) | 3.8391848048 |
| `rr_p2_n96.json` | λ₃ ≤ (P2, n = 96) | 3.8388872830 |
| `first_order_FA_n192_lg.json` | b(e,e)/m < 0 along F_A, λ₁ ∈ [3.8388757010, 3.8388872830] | ≤ −3.0383475 |
| `first_order_FA_n192_su13.json` | the same with the SU13 ball (optional cross-check) | ≤ −3.0875694 |
| `kyfan_FB_n96_lg.json` | D²Λ(v,v) < 0 along F_B, λ₁ ∈ [3.8388757010, 3.8388872830] | ≤ −9.3411481 |
| `kyfan_FB_n96_su13.json` | the same with the SU13 ball (optional cross-check) | ≤ −12.8504761 |
| `rank.json` | a 6×6 minor of the trace-coordinate Jacobian is non-zero | det ≈ 7.8724e18 |

The `_lg` certificates use exactly λ₁ ∈ [3.8388757010, 3.8388872830]. Their bounds are slightly weaker than item
23's (−3.0383476 and −9.3411545), because item 23's scripts (`lb_cert_first_order.py`, `lb_cert_sign.py`) used the
ball arb("3.83888149202636215", "5.791e-6"). That ball starts at 3.8388757010263, slightly above 3.8388757010 but
still below the certified Lehmann–Goerisch bound 3.8388757010527. Both versions certify the sign. The "certified
value" column is the upper end of the ball stored in the certificate (rounded up), i.e. what the certificate attests
and the checker re-derives.

**Known errata in the copied documents.** PROOF.md Lemma 8 (line 495) attributes −9.341154560446736 to the λ ball
[3.8388757010, 3.8388872830]. For exactly that ball, the certified value is −9.3411486 (`kyfan_FB_n96_lg.json`).
The −9.341154560 came from item 23's script ball, which starts at 3.8388757010263 (see above). The proof uses only
the constant −9.3411, so it is unaffected. This will be corrected upstream in item 23. The copies here are
unchanged.

## Trust model

Read this before reading "VALID" as "proved".

**What the checker re-derives.** `bolzakit_check` imports nothing from `bolzakit`. It re-does the mathematics from
the numbers stored in each certificate, in its own Arb code, and a tampered claim fails even if the sha256 is
recomputed.
- **Counts** (`lp_*`) are re-derived fully independently: the Sturm sign checks, the window minimum, the
  integral and the geodesic sum are recomputed from the test-function coefficients and closed-form geodesic
  lengths.
- **Rank** is re-derived from the octagon coordinates and the trace words: the checker builds its own Jacobians
  and recomputes the minors, the stacked determinant and the residual. A stored Jacobian is never trusted.
- **First-order and Ky Fan** are rebuilt from the stored raw sums and Gram matrices (q, m and their
  t-derivatives; the 7×7 Gram blocks, β bounds and the λ ball). The assembly into b(e,e)/m or g₂ is recomputed in
  pure Arb. The checker reports these as "VALID given stored Gram data".
- **Rayleigh–Ritz** (`rr_*`): the checker re-checks positive definiteness of μM − K from the stored Gram matrices
  (also "VALID given stored Gram data").
- `--deep` re-integrates a small sample of elements. That is a consistency check, not a proof.
- The first-order and Ky Fan checker formulas transliterate the producer's derivation. They are independent code,
  not an independent derivation.

**What is trusted.**
- **The element assembly behind the Gram data.** The checker does not re-integrate the finite-element forms over
  the whole mesh, so the stored Gram matrices and raw sums are trusted to be what the bolzakit producer computed.
- **python-flint / Arb** (ball arithmetic, `acb.integral`) and sympy, in both the producer and the checker.
- **The cited inputs** of PROOF.md §Dependencies: FBP's lemmas and geodesic data, the Selberg trace formula, the
  Lehmann–Goerisch theorem, and the standard facts.
- **The λ₁ lower bound itself.** 3.8388757010 enters the `_lg` certificates as an explicit input whose source
  string names item 23's LB_RESULTS. Its Lehmann–Goerisch computation (`lb_bolza.py` in item 23) is not
  regenerated or checked by bolzakit v0.2.
- **The other hand-linked inputs.** The checker takes the λ₁ upper end 3.8388872830 and the gap λ₄ ≥ 4.0 as given
  inputs of the first-order and Ky Fan certificates. They are discharged by `rr_p2_n96` and by `lp_A_3.3` with
  `lp_B_3.3_4.0` respectively, but a person makes that link, not a tool. Likewise, `lp_B`'s hypothesis
  λ₁ ≥ 3.3 is discharged by hand by `lp_A`. The β bounds in the cluster certificates are stored data, trusted
  like the Gram data.
- **The non-certificate steps.** Lemmas 2–4, 6–7, 10–11 and the assembly are written arguments or exact
  computations in item 23. They are not certificates and the checker does not see them. Lean
  (see [LEAN.md](LEAN.md)) machine-checks the logical chain, conditional on these computational inputs.

**What has not happened.** No human expert has read the proof. The reviews are AI second reads.
