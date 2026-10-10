# Frequency trimming of a split mode pair (`bolzakit.sensitivity.trim`)

Label: NUMERICAL (discrete problem). Imports only numpy, scipy and `bolzakit.sensitivity.core`.

## Problem

A nominally degenerate pair `(i, i+1)` of `K x = λ M x` (the n = 2 wineglass pair of a ring gyroscope, a disk
resonator, a symmetric shell) is split by fabrication imperfections. Candidate trims `l = 0..L−1` (a point mass at a
rim site, ablated material, a stiffness patch) change the pencil linearly in their amounts `t_l` (positive adds,
negative removes): `K(t) = K + Σ t_l dK_l`, `M(t) = M + Σ t_l dM_l`. `trim_degenerate_pair` returns the cheapest `t` within bounds that
ties the pair. The cost is the trim magnitude `Σ c_l |t_l|` (`c ≥ 0`, default ones, i.e. the total amount), so removal
(`lower < 0`, e.g. ablation) costs the same as addition.

## Method

In the pair's eigenbasis `U` the reduced model is `diag(−g/2, g/2) + Σ dt_l A_l` with the splitting matrices
`A_l = Uᵀ(dK_l − λ̄ dM_l)U` and the split `g`. The pair is tied to first order iff `G dt = (g, 0)` with
`G = [A_00 − A_11; A_01]`: two equations, whatever the number of candidates. The first step is the linear program
`min Σ c_l |t_l + dt_l|` subject to these equations and the bounds, posed over `t = p − q` with `p, q ≥ 0` (HiGHS); a
vertex uses at most two trims. Later steps are Newton corrections on the selected trims (least squares, exact re-solve
each time); the LP is re-solved only if a correction leaves the bounds or flips the sign of a trim. Convergence is
quadratic until the split reaches the eigensolver's resolution, estimated as `rel_floor = eps·λ_max/|λ̄|`. Splits below
`rel_floor` are reported as "at the resolution floor (≤ rel_floor)", never as a measured value. The estimate is
conservative and can be loose: a run whose split falls below it still reports `"converged"` if it met `tol`.

## API

`trim_degenerate_pair(K, M, dKs, dMs, pair, *, lower=0.0, upper=None, cost=None, tol=1e-10, max_iter=20,
gap_ratio_min=1e-3) -> TrimResult` with fields `t`, `cost`, `split_before`, `split_after`, `rel_split_before`,
`rel_split_after`, `eigvals_before`, `eigvals_after`, `rel_floor`, `status`, `n_iter`, `n_eigensolves`, `history`
(`TrimRecord`s: `it`, `t_cost`, `split`, `rel_split`, `lam_bar`, `step`), `label`, `caveat`. Statuses:
`"converged"` (split ≤ tol·|λ̄|), `"stalled"` (the resolution floor: two iterations without halving the best split,
or out of iterations within a factor 2 of `rel_floor`), `"max_iter"`, `"infeasible"` (no trim pattern within the
bounds ties the pair to first order, e.g. no candidate couples the two modes). The pair must stay isolated
(`ClusterSplitError` otherwise).

## Example (`examples/gyro_trim.py`)

Output of a fresh run (splits are relative, `split/λ̄`; trim masses are relative to the structure's
mass):

```
(a) ring, 4th-harmonic imperfection: rate/lam_bar = 0.29984 (thin-ring theory 0.3)
    split/lam_bar 5.997e-04 -> at floor (<= 1.6e-08) (stalled, 7 eigensolves); dense check at floor (<= 1.6e-08); trim mass 1.9385e-04 (6.17e-04 of the ring) at sites [2, 3]
    history: 6.0e-04 5.2e-07 4.5e-09 3.2e-09 1.4e-09 1.4e-09 4.1e-09
    baseline undamped antinode: tied, split at floor (<= 1.6e-08), 12 eigensolves, mass 2.6219e-04 at sites [3, 10, 11, 12, 13], 1.35x our mass
    baseline damped 0.5 antinode: tied, split at floor (<= 1.6e-08), 19 eigensolves, mass 1.9374e-04 at sites [2, 3, 10, 11], 1.00x our mass

(b) ring on 8 spokes, random E and rho (sigma 0.5 %); nominal pair lambda_6 = lambda_7 = 0.163162 (rel. diff 2.0e-12)
 seed  split before   ours: after  dense check  eigsolves  mass (rel)     sites
    0    8.913e-04     8.07e-13     5.54e-11          4   8.09e-04   [23, 25]
    1    3.816e-04     7.45e-13     3.42e-11          4   3.29e-04     [1, 7]
    2    5.372e-04     1.04e-12     1.85e-11          4   2.35e-04   [22, 23]
    3    3.581e-04     6.19e-13     8.36e-11          4   3.46e-04   [11, 13]
    4    7.140e-04     4.58e-13     1.12e-10          4   7.40e-04     [7, 9]
      single-site antinode baselines (mass ratio = baseline mass / our mass, only if it tied the pair):
    0  undamped: 5.63x mass, 29 eigsolves  |  damped 0.5: did not tie (oscillates / slow; 1.8e-06 after 40 rounds)
    1  undamped: did not tie (oscillates / slow; 2.0e-09 after 40 rounds)  |  damped 0.5: did not tie (oscillates / slow; 7.5e-07 after 40 rounds)
    2  undamped: 3.44x mass, 24 eigsolves  |  damped 0.5: did not tie (oscillates / slow; 1.4e-06 after 40 rounds)
    3  undamped: 2.43x mass, 23 eigsolves  |  damped 0.5: did not tie (oscillates / slow; 3.4e-06 after 40 rounds)
    4  undamped: 1.20x mass, 13 eigsolves  |  damped 0.5: did not tie (oscillates / slow; 2.1e-06 after 40 rounds)
label: NUMERICAL (discrete problem); 325.3 s
```

The `history` line prints the raw relative split of each iterate; entries below the plain ring's resolution
floor (1.6e-8) are not measurements.

What this shows. The baseline is the usual iterated single-site antinode trimming: each round adds mass at the site
where the higher mode moves most, sized (undamped) or half-sized (damped 0.5, the fair standard) to cancel the split
along that one site, up to 40 rounds. Mass ratios count only baseline runs that tied the pair; "did not tie" means
not within the 40-round budget, not never.

- **Spoked ring (no analytic rule), the headline.** Ours tied every seed in 4 eigensolves, from relative splits of
  3.6e-4 to 8.9e-4 down to 4.6e-13 to 1.0e-12 (independent dense check 1.9e-11 to 1.1e-10). Undamped antinode
  trimming tied 4 of 5 seeds, using 1.2 to 5.6 times our mass and 13 to 29 eigensolves; on seed 1 it did not tie.
  The 0.5-damped baseline tied no seed within 40 rounds (residual splits 7.5e-7 to 3.4e-6).
- **Plain ring.** Ours reached the resolution floor (≤ 1.6e-8) in 7 eigensolves (status `stalled`, the floor). Both
  baselines also tied: undamped with 1.35 times our mass in 12 eigensolves, damped with the same mass as ours (1.00×)
  in 19 eigensolves. There is no mass advantage against the damped baseline here; the gain is the eigensolve count.

The masses compare the tool's objective `Σ c|t|` with `c = 1`; every trim in the example is a mass addition, so this
is the added mass. With ablation candidates the same objective charges removal like addition.

## Limits

First-order model of one isolated pair: trims that make a third mode join the pair are not handled. Only the frequency
split is cancelled; the mode orientation (quadrature coupling to the drive axes) is not a target. The trims must act
linearly on (K, M); geometric trims that need remeshing must be linearised by the user. The example's point masses are
lumped nodal masses on an idealised plane-stress model. The `cost` weights are the user's: the tool minimises
`Σ c_l |t_l|`, which equals the added mass only when every trim is a mass addition with `c_l = 1`.
