# Verification policy

This policy is adapted from the math-pursuits project, where the code was developed. It defines the provenance labels
used in bolzakit and the rules the evidence behind them has to meet.

## Labels

| Label | Meaning | Required evidence |
|---|---|---|
| EXACT | Integer or algebraic arithmetic with no rounding | Computation in Z[zeta_8], integers or exact PARI; no floating point in the chain |
| PROVED | Exact computation plus a written argument, or a Lean theorem | Exact computation plus a written argument, or a Lean theorem; arguments may be AI-written and cited — see the disclosure below. Every external theorem quoted with hypotheses (rule C); adversarial pass (rule A); second-reader review |
| CERTIFIED | Finite computation with rigorous error control | Rigorous enclosure, by interval/ball arithmetic or a rigorous a-priori floating-point error analysis; explicit truncation and tail bounds, no hard-coded flags, reproducible script (rule N) |
| NUMERICAL | Strong numerical evidence | High-precision computation, precision-doubling stability, stated precision or error |
| QUOTED | Value taken from a published source | Exact reference and, for data files, a fixed checksum; the value is not re-derived |

A label is never upgraded relative to the spec §3 table (`docs/superpowers/specs/2026-10-08-bolzakit-design.md`).

## Rule C: citations

1. Quote each external theorem verbatim (or precisely paraphrased), with every hypothesis, plus source and theorem
   number, from a source that has actually been read.
2. Check explicitly that each hypothesis holds in the application.
3. "Classical", "well known" and "standard" are not citations.

## Rule N: numerics

1. "Certified" requires a rigorous enclosure: interval/ball arithmetic (Arb, `mpmath.iv`, exact PARI) or a rigorous
   a-priori floating-point error analysis. A difference between two precisions is not an error bound.
2. Finite checks of a universally quantified statement are NUMERICAL, unless an analytic bound covers the rest and that
   bound is itself proved with explicit constants.
3. Pass/fail values are computed, never typed into an export.
4. Saddle-point and asymptotic arguments state the saddle equation exactly and give explicit error bounds.
5. Normalisations are written out.

## Rule A: adversarial pass

Before PROVED, for each lemma: try to construct a counterexample (small cases, degenerate parameters), test the converse
if it is used, and check every "by symmetry" and "the cross term vanishes". Log the attempts.

## Rule X: circularity

Inputs that presuppose the conclusion (a fit to the formula being tested, an identity that is the definition rewritten)
cannot be evidence. State the data provenance.

## Certified claims (v0.2)

There is one row for each certified claim of the Bolza local-maximum proof. The golden certificates are in
`proofs/bolza-local-max/certificates/`. "Checker status" is the output of
`python -m bolzakit_check proofs/bolza-local-max/certificates/*.json`, which exits 2 on this set, as expected (six
conditional certificates; the cluster rows abbreviate the λ₁ / gap sources the checker prints). For what each status re-derives and what it trusts, see
[docs/certified.md](docs/certified.md#trust-model). "Hand-linked" names the conditions that a person, not a tool,
discharges with another row.

| Claim | Certificate (golden file) | Certified value | Checker status | Hand-linked |
|---|---|---|---|---|
| no eigenvalue in (0, 3.3] (trace-formula count < 1, mode A) | `lp_A_3.3.json` | 0.912612… | VALID | — |
| at most 3 eigenvalues in [3.3, 4.0] (count < 4, mode B) | `lp_B_3.3_4.0.json` | 3.869661… | VALID (conditional on λ₁ ≥ 33/10) | λ₁ ≥ 3.3 ← `lp_A_3.3` |
| λ₃ ≤ 3.8391848048 (P1 Rayleigh–Ritz, n = 192) | `rr_p1_n192.json` | 3.8391848048 | VALID given stored Gram data | — |
| λ₃ ≤ 3.8388872830 and λ₁+λ₂+λ₃ ≤ 11.5166618386 (P2 Rayleigh–Ritz, n = 96) | `rr_p2_n96.json` | 3.8388872830 | VALID given stored Gram data | — |
| λ₁ ≥ 3.83887570105 (Lehmann–Goerisch, n = 48, ρ = 5) | `lg_lower_n48.json` | 3.8388757010527… | VALID given stored Gram data (conditional on λ₄ ≥ 4.0) | λ₄ ≥ 4.0 ← `lp_A_3.3` + `lp_B_3.3_4.0`; the computation cannot detect a wrong ρ |
| A ≠ 0: b(e,e)/m < 0 along F_A (P1, n = 192, LG λ₁ ball) | `first_order_FA_n192_lg.json` | ≤ −3.0383475 | VALID given stored Gram data (conditional on λ₁ ∈ [3.8388757010, 3.8388872830] [D15: …]; λ₄ ≥ 4.0 [gap: …]) | λ₁ ∈ [3.8388757010, 3.8388872830] ← `lg_lower_n48` + `rr_p2_n96`; λ₄ ≥ 4.0 ← counts |
| the same with the SU13 λ₁ ball (optional cross-check) | `first_order_FA_n192_su13.json` | ≤ −3.0875694 | VALID given stored Gram data (conditional on the SU13 λ₁ ball [D1: …]; λ₄ ≥ 4.0 [gap: …]) | λ₁ ball from SU13 (D1, quoted); λ₄ ≥ 4.0 ← counts |
| D²Λ(v,v) < 0 along F_B (Ky Fan, P2, n = 96, LG λ₁ ball) | `kyfan_FB_n96_lg.json` | ≤ −9.3411481 | VALID given stored Gram data (conditional on λ₁ ∈ [3.8388757010, 3.8388872830] [D15: …]; λ₄ ≥ 4.0 [gap: …]) | λ₁ ball ← `lg_lower_n48` + `rr_p2_n96`; λ₄ ≥ 4.0 ← counts |
| the same with the SU13 λ₁ ball (optional cross-check) | `kyfan_FB_n96_su13.json` | ≤ −12.8504761 | VALID given stored Gram data (conditional on the SU13 λ₁ ball [D1: …]; λ₄ ≥ 4.0 [gap: …]) | λ₁ ball from SU13 (D1, quoted); λ₄ ≥ 4.0 ← counts |
| the trace chart has rank 6 at Bolza (non-zero 6×6 minor) | `rank.json` | det ≈ 7.8724e18 | VALID (rank 6 …, surface bolza (Klein vertices checked against the exact octagon)) | — |

The checker's `first_order` and `kyfan` formulas transliterate the producer's derivation (independent code, the same
derivation). The `count` and `lg_lower` re-derivations are derivation-independent. No row covers the element
assembly behind the Gram data: `--deep` samples it as a consistency check only.

## How this was developed, and what has not happened

bolzakit's code was ported from the math-pursuits research project. The derivations there, and the second-reader
reviews of them, were produced by AI agents (Claude) working under Gordon Speagle's direction. No human specialist or
peer review has checked them yet. The "second reader" in these rules was therefore another AI agent, not an independent
human. The Lean proofs are machine-checked, and the test suite compares numbers against independent computations, but
PROVED claims that rest on AI-written arguments should be treated accordingly until a human has reviewed them.

bolzakit's code, tests, documentation and Lean files were also written and reviewed by AI agents (Claude) under Gordon
Speagle's direction; Gordon directed and approved the work.
