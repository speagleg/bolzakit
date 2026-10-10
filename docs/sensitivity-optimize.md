# Repeated-eigenvalue optimisation (`bolzakit.sensitivity.optimize`)

Label: NUMERICAL (discrete problem). First order only (no Overton SQP); see [Limits](#limits).

## Problem

Maximise `f(p) = λ_k(p)` (`LowestEigenvalue(k)`) or `λ_{k+1}(p) − λ_k(p)` (`Gap(k)`) for a `ParametricPencil(K, M,
dK, dM)` (callables of `p`, dense or sparse, with analytic partial derivatives), subject to `lower ≤ p ≤ upper` and
`A_eq p = b_eq`.

## Steepest ascent at a multiple eigenvalue (Overton; Seyranian–Lund–Olhoff)

At a tie `f'(p; d) = λ_min(A(d))` with `A(d) = Σ d_l A_l` the splitting matrices of the active cluster; this is
concave in `d`, so `λ_k` has directional derivatives and no gradient. The steepest feasible ascent solves
`max_{d∈𝒦, ‖d‖≤1} λ_min(A(d))`, with `𝒦 = {A_eq d = 0, s_l d_l ≥ 0 at bounds}`, whose dual is

`φ* = min_{Y ⪰ 0, tr Y = 1, μ ≥ 0} ‖P_eq(g(Y) + Σ_{l∈B} μ_l s_l e_l)‖`, `g(Y)_l = ⟨Y, A_l⟩`.

`ascent_direction` solves it by FISTA with adaptive restart plus an exact face polish and returns `d`, the primal rate
`phi`, the certified dual bound `phi_upper` (`phi ≤ φ* ≤ phi_upper`), the dual matrices and bound multipliers, and the
predicted branch rates. `φ* = 0` is the SLO/KKT first-order condition.

## The driver `maximize_eigenvalue`

Each iteration: eigensolve with single-linkage cluster growth; status checks; `splitting_matrices`;
`ascent_direction` on the stacks normalised by the eigenvalue scale; a trust-region step clipped at the first bound;
the Poincaré **screen** (a Ritz upper bound on the trial `λ_k` rejects hopeless steps without an eigensolve; it never
accepts a step); and an Armijo test on true eigenvalues that also requires strict increase. An ε-refinement phase
shrinks the clustering tolerance when the tied model stalls on a cluster that is not exactly tied.

Statuses (never exceptions): `"stationary"`, `"subproblem_not_converged"`, `"step_too_small"`, `"max_iter"`,
`"lower_coalescence"`, `"gap_closed"`; see the docstring for what each certifies. Accepted iterates are monotone.

## API

- `ParametricPencil(K, M, dK, dM=None)`: callables of `p`; `dK(p)` (and `dM(p)`) return the list of partial
  derivatives, one per parameter. Objectives: `LowestEigenvalue(k=0)` (maximise `λ_k`) and `Gap(k)` (maximise
  `λ_{k+1} − λ_k`).
- `maximize_eigenvalue(pencil, p0, objective, *, lower=None, upper=None, A_eq=None, b_eq=None, step0=None,
  max_iter=100, mult_rtol=1e-3, stat_tol=1e-8, armijo=1e-4, screen=True, refine=True, mult_rtol_min=1e-11,
  refine_max_iter=100, callback=None) -> EigenOptResult`. `mult_rtol` is the clustering tolerance (relative to the
  eigenvalue scale); `refine`, `mult_rtol_min` and `refine_max_iter` control the ε-refinement (it shrinks the
  clustering tolerance down to `mult_rtol_min`, and the outer budget grows to `max_iter + refine_max_iter` once it
  starts). `callback(IterRecord)` is called for every record. Raises `ValueError` for an infeasible `p0` or an
  out-of-range `objective.k`, `TypeError` for an unknown objective; the core's `PencilError` / `ClusterSplitError`
  propagate.
- `EigenOptResult`: `p`, `value`, `eigvals`, `multiplicity`, `phi`, `phi_upper` (the certified bound on the best
  first-order ascent rate), `Y_upper`, `Y_lower`, `mu` (dual certificate), `status`, `n_iter`, `n_eigensolves`,
  `n_screened`, `history` (`IterRecord`s: `it`, `value`, `multiplicity`, `phi`, `phi_upper`, `tau`, `accepted`,
  `screened`, `rates`, `refined`, the last true for an ε-refinement round), `mult_rtol` (the clustering tolerance in
  force at the end), `eps` (the largest eigenvalue gap inside the active clusters: the reduced model treats
  eigenvalues that far apart as tied), `fallback_from`, `label`, `caveat`. When an ε-refinement does not pay off, the
  result reports the coarser model with its own status (`"stationary"` or `"step_too_small"`), and `fallback_from`
  keeps the stop that actually ended the run (`"max_iter"` means the iteration or eigensolve budget ran out);
  otherwise `fallback_from` is `None`.
- `ascent_direction(A_upper, A_lower=None, *, null_basis=None, bound_signs=None, tol=1e-10, max_iter=10000,
  stat_tol=0.0, f_value=0.0) -> AscentDirection` (`d`, `phi`, `phi_upper`, `duality_gap`, `Y_upper`, `Y_lower`,
  `mu`, `rates`, `rates_lower`, `n_iter`, `stop`). Its stops: `"gap"` (a certified ascent direction of rate `phi`),
  `"stationary"` (`φ* = 0` to accuracy), `"no_progress"` (two consecutive restarts without progress: the best
  certified point seen) and `"max_iter"`; with the last two, `d` may be zero and `phi_upper` then bounds `φ*`.

## Benchmark against standard practice

`benchmarks/frequency_maximisation.py` maximises the fundamental frequency of a Kirchhoff plate (simply supported
edges plus a central post; Morley, 512 dofs; 5×5 thickness zones; volume equality; `0.5 ≤ h ≤ 1.5`). The optima
found have a double or triple λ₀. Every method gets the same problem, eigensolver and budget (200 eigensolves) from
five starts (start 0 uniform, 1–4 random feasible): ours; the bound formulation `max β s.t. λ_i ≥ β (i < 4)` with
simple-eigenvalue gradients, solved by scipy SLSQP, scipy trust-constr and pyMOTO's MMA; and naive fixed-step projected
gradient ascent on λ₀. The table is pasted from `benchmarks/frequency_maximisation.md` as generated (commit e9c0b79; only its prose
wording was edited afterwards, with no rerun):

Plate + central post, 512 dofs, 25 thickness zones (unequal areas, 0.03125 to 0.046875; volume weighted by a_z). Budget 200 eigensolves per run, 5 starts (start 0 uniform, the rest random feasible). Bound formulation: beta / lambda_0(start) in [0.5, 2] for SLSQP, trust-constr and MMA. Versions: scipy 1.18.1, numpy 2.5.3, scikit-fem 12.0.2, pymoto 2.0.1.

Every method reports its best feasible evaluated design within the budget (the 'best λ₀/π⁴ within budget' column). 'Best from a start' is the maximum of those values over methods: the best value found within the 200-eigensolve budget, not the optimum. 'reached' counts eigensolves until the running best is within 0.1 % (99.9 %) or 1 % (99 %) of the best from that start. Feasibility for the best-design record: volume residual ≤ 1e-9 relative, except MMA, whose volume tolerance is 1e-3 relative (pyMOTO MMA takes the volume as an inequality, ≤). 'decreases' counts decreases between successive accepted iterates (ours: accepted steps; SLSQP and trust-constr: callback iterates, where a merit-function method may decrease lambda_0 legitimately; MMA and naive: every iterate, all accepted).

| start | method | best λ₀/π⁴ within budget | mult | eigensolves to 99 % | to 99.9 % | decreases | volume residual | stop |
|---|---|---|---|---|---|---|---|---|
| 0 | ours | 1.9615 | 3 | — | — | 0 | 0e+00 | stationary |
| 0 | SLSQP bound | 1.9615 | 3 | — | — | 35 | 1e-10 | budget |
| 0 | trust-constr bound | 2.0003 | 1 | 186 | 194 | 3 | 2e-16 | budget |
| 0 | MMA bound (pyMOTO) | 1.9036 | 2 | — | — | 86 | 9e-06 | budget |
| 0 | naive fixed-step | 1.8850 | 1 | — | — | 103 | 0e+00 | budget |
| 1 | ours | 2.0515 | 1 | 154 | 169 | 0 | 2e-16 | max_iter |
| 1 | SLSQP bound | 2.0521 | 2 | 39 | 40 | 25 | 0e+00 | converged |
| 1 | trust-constr bound | 2.0074 | 2 | — | — | 2 | 1e-16 | budget |
| 1 | MMA bound (pyMOTO) | 1.8895 | 2 | — | — | 87 | 3e-05 | budget |
| 1 | naive fixed-step | 1.8841 | 1 | — | — | 80 | 0e+00 | budget |
| 2 | ours | 2.0083 | 1 | — | — | 0 | 0e+00 | max_iter |
| 2 | SLSQP bound | 2.0521 | 2 | 36 | 39 | 26 | 2e-16 | converged |
| 2 | trust-constr bound | 2.0076 | 2 | — | — | 2 | 2e-16 | budget |
| 2 | MMA bound (pyMOTO) | 1.8858 | 2 | — | — | 86 | 2e-05 | budget |
| 2 | naive fixed-step | 1.8821 | 1 | — | — | 84 | 0e+00 | budget |
| 3 | ours | 2.0083 | 1 | — | — | 0 | 2e-16 | max_iter |
| 3 | SLSQP bound | 2.0084 | 3 | — | — | 29 | 8e-13 | budget |
| 3 | trust-constr bound | 2.0083 | 2 | — | — | 1 | 7e-16 | budget |
| 3 | MMA bound (pyMOTO) | 2.0415 | 1 | 194 | 199 | 50 | 6e-05 | budget |
| 3 | naive fixed-step | 1.8841 | 1 | — | — | 83 | 0e+00 | budget |
| 4 | ours | 2.0084 | 2 | — | — | 0 | 0e+00 | max_iter |
| 4 | SLSQP bound | 2.0521 | 2 | 43 | 53 | 23 | 0e+00 | converged |
| 4 | trust-constr bound | 1.9699 | 2 | — | — | 1 | 0e+00 | budget |
| 4 | MMA bound (pyMOTO) | 1.8989 | 2 | — | — | 91 | 7e-06 | budget |
| 4 | naive fixed-step | 1.8881 | 1 | — | — | 85 | 1e-16 | budget |

```
summary over starts (reached = within 0.1 % of the best value found within budget from that start):
  ours                 reached 1/5; median eigensolves to reach  169.0; decreases    0; max volume residual 2.2e-16; max bound violation 0.0e+00
  SLSQP bound          reached 3/5; median eigensolves to reach   40.0; decreases  138; max volume residual 1.0e-10; max bound violation 0.0e+00
  trust-constr bound   reached 1/5; median eigensolves to reach  194.0; decreases    9; max volume residual 6.7e-16; max bound violation 0.0e+00
  MMA bound (pyMOTO)   reached 1/5; median eigensolves to reach  199.0; decreases  400; max volume residual 6.4e-05; max bound violation 0.0e+00
  naive fixed-step     reached 0/5; median eigensolves to reach    nan; decreases  435; max volume residual 1.1e-16; max bound violation 0.0e+00
```

Reach counts are relative to the best value found within budget from each start, not to the optimum. From start 3 the best within budget, MMA's 2.0415, is a point mid-way along its trajectory when the budget ran out; 2.052 is reachable from there, since SLSQP started from that design reaches 2.05205 (a review check after this run, not part of the table).

Stop codes: `budget` = eigensolve budget exhausted; `converged` / `failed` = solver's own termination. Full solver messages:
- SLSQP bound: Optimization terminated successfully

What the table shows. "Best from a start" is the best value any method found within the 200-eigensolve budget from
that start, not the optimum, and "reached" means coming within 0.1 % of it. On this problem the SLSQP bound formulation
beat our method on both success rate and cost: it reached the best from 3 of 5 starts (1, 2 and 4, all at
λ₀/π⁴ = 2.0521) in a median of 40 eigensolves. Our method reached it from 1 of 5 starts (start 1, 2.0515, after 169
eigensolves). From starts 2 to 4 it stopped at `max_iter` near 2.008 (where SLSQP from start 3 and trust-constr from
starts 2 and 3 also ended), and from the symmetric start 0 it stopped `stationary` at the symmetric triple point
1.9615, where SLSQP also ended (on a budget stop); only trust-constr left the symmetric subspace there (2.0003, reached at 194 eigensolves).
trust-constr and MMA each reached the best from 1 of 5 starts (medians 194 and 199 eigensolves), and naive fixed-step
ascent from none. MMA's start-3 value, 2.0415, is a point mid-way along its trajectory when the budget ran out, not an
optimum: SLSQP started from that design reaches 2.052.

Where our method wins: its accepted iterates never decreased (0 decreases; the other methods' counts, 9 for
trust-constr, 138 for SLSQP and 400 to 435 for MMA and naive ascent, count different things, see below), and every step comes with a certified bound `phi_upper` on the best first-order
ascent rate, so a `stationary` stop certifies first-order stationarity of the reduced model. It does not win on
success rate or eigensolve count here. Every method met the constraints to rounding (volume residual ≤ 1e-10) except
MMA (volume as an inequality with tolerance 1e-3; residual up to 6.4e-5). The decrease counts are not comparable
across methods, because each method's iterates are counted differently (see the definition above the table): ours
counts accepted steps, SLSQP and trust-constr their callback iterates (a merit-function method may decrease λ₀
legitimately), MMA and naive ascent every iterate. What the counts do show for MMA (pyMOTO defaults) and naive
ascent is that λ₀ went down between successive iterates 50 to 91 and 80 to 103 times per run; the table does not
record whether those decreases are an oscillation.

## Limits

First-order steepest ascent only: it can end at `max_iter` short of the maximum. In the tests this happens on an
ill-conditioned smooth simplex instance, where the iterates zigzag. On the benchmark, four of the five starts ended
with status `max_iter` (there the 200-eigensolve budget, `_max_eigensolves`); why those runs had not converged is not
established. On the benchmark it needed several times more eigensolves than SLSQP. No Overton SQP on the fixed-multiplicity manifold
and no quasi-Newton step (both deferred); no gap opening from a closed gap; `LowestEigenvalue(k>0)` stops at
`lower_coalescence`; analytic partials required; real symmetric pencils only. From a symmetric start the iterates
stay in the symmetric subspace. Scale invariance holds to rounding for converged runs; runs that stop at `max_iter`
or `step_too_small` can differ at rounding level under rescaling.
