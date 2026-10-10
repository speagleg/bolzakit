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

A label is never upgraded relative to the per-function table in [docs/provenance.md](docs/provenance.md).

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

## Sub-project F (0.3.0): environment

Local runs used the bolzakit venv (Python 3.12) with numpy 2.5.3, scipy 1.18.1, scikit-fem 12.0.2; the benchmark also
used pymoto 2.0.1 (a benchmark-only tool, not a dependency).

## Sub-project F (0.3.0): measurements

| Item | Check | Measured | Where |
|---|---|---|---|
| Membrane | rates −8π², −2π² (rtol 2e-3); S'' = 30π², Ky Fan prefix 36π² (rel. 5e-3) | unchanged from 0.2 | `tests/sensitivity/test_skfem.py` |
| Ring rate/λ̄, 4th harmonic (w/R 0.05, Nθ 96) | 3/10 within 5e-4 | 0.29984 | `examples/gyro_trim.py`, `test_trim.py` |
| Ring trim (a) | split reduced ≥ 1e3, dense check | 6.0e-4 → ≤ rel_floor (1.6e-8), dense check also ≤ rel_floor; `stalled` (the floor), 7 eigensolves | same |
| Spoked ring trim (b), 5 seeds | `converged`, ≥ 1e3, dense check | before 3.6e-4 to 8.9e-4; after ≤ rel_floor (2.0e-10) on every seed (below the eigensolver's resolution estimate, so no value is quoted); dense LAPACK check ≤ 1.1e-10 (its own floor); 4 eigensolves each | same |
| Trim mass vs iterated antinode | reported | plain ring: undamped baseline 1.35× our mass (12 eigensolves), 0.5-damped 1.00× (19 eigensolves); ours 7 eigensolves, no mass advantage against the damped baseline. Spoked ring: undamped baseline tied 4/5 seeds at 1.2–5.6× our mass in 13–29 eigensolves; 0.5-damped tied 0/5 within 40 rounds | same |
| Frequency benchmark (budget 200, 5 starts) | reported as measured | within 0.1 % of the best found within budget (not the optimum): SLSQP 3/5 (median 40 eigensolves), ours 1/5 (169), trust-constr 1/5 (194), MMA 1/5 (199; start-3 point mid-trajectory), naive 0/5; decreases: ours 0, SLSQP 138 | `benchmarks/frequency_maximisation.md` |
| Post-supported plate spectrum | λ/π⁴ = 1.58608, 1.83370, 1.83370, 4.37522, 6.76920 | reproduced to 5e-5 | `tests/test_frequency_benchmark.py` |
| `ascent_direction` d = 64, m = 4 (§13) | < 50 ms (a target; the test prints the time and does not assert it) | 12.9 ms in one test run; 21 repeats: median 12.2 ms, max 38.9 ms (local dev machine, under load from another job; timings vary by machine) | `tests/sensitivity/test_optimize.py::test_generic_d64_m4_fast` |
