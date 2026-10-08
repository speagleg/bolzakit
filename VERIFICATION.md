# Verification policy

This policy is adapted from the math-pursuits project, where the code was developed. It defines the provenance labels
used in bolzakit and the rules the evidence behind them has to meet.

## Labels

| Label | Meaning | Required evidence |
|---|---|---|
| EXACT | Integer or algebraic arithmetic with no rounding | Computation in Z[zeta_8], integers or exact PARI; no floating point in the chain |
| PROVED | Complete proof, every step checked | Written or Lean proof; every external theorem quoted with hypotheses (rule C); adversarial pass (rule A); second-reader review |
| CERTIFIED | Finite computation with rigorous error control | Interval or ball arithmetic over every operation, explicit truncation and tail bounds, no hard-coded flags, reproducible script (rule N) |
| NUMERICAL | Strong numerical evidence | High-precision computation, precision-doubling stability, stated precision or error |
| QUOTED | Value taken from a published source | Exact reference and, for data files, a fixed checksum; the value is not re-derived |

A label is never upgraded relative to the provenance table in `docs/provenance.md`.

## Rule C: citations

1. Quote each external theorem verbatim (or precisely paraphrased), with every hypothesis, plus source and theorem
   number, from a source that has actually been read.
2. Check explicitly that each hypothesis holds in the application.
3. "Classical", "well known" and "standard" are not citations.

## Rule N: numerics

1. "Certified" requires enclosures (Arb, `mpmath.iv`, interval arithmetic, exact PARI). A difference between two
   precisions is not an error bound.
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

## How this was developed, and what has not happened

bolzakit's code was ported from the math-pursuits research project. The derivations there, and the second-reader
reviews of them, were produced by AI agents (Claude) working under Gordon Speagle's direction. No human specialist or
peer review has checked them yet. The "second reader" in these rules was therefore another AI agent, not an independent
human. The Lean proofs are machine-checked, and the test suite compares numbers against independent computations, but
PROVED claims that rest on AI-written arguments should be treated accordingly until a human has reviewed them.
