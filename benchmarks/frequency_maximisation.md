# Fundamental-frequency maximisation benchmark

Label: NUMERICAL (discrete problem). Reported as measured, whoever wins.

Plate + central post, 512 dofs, 25 thickness zones (unequal areas, 0.03125 to 0.046875; volume weighted by a_z). Budget 200 eigensolves per run, 5 starts (start 0 uniform, the rest random feasible). Bound formulation: beta / lambda_0(start) in [0.5, 2] for SLSQP, trust-constr and MMA. Versions: scipy 1.18.1, numpy 2.5.3, scikit-fem 12.0.2, pymoto 2.0.1.

Every method reports its best feasible evaluated design within the budget; 'reached' counts eigensolves until the running best is within 0.1 % (99.9 %) or 1 % (99 %) of the best final value over methods from that start. 'decreases' counts decreases between successive accepted iterates (ours: accepted steps; SLSQP and trust-constr: callback iterates, where a merit-function method may decrease lambda_0 legitimately; MMA and naive: every iterate, all accepted).

| start | method | final λ₀/π⁴ | mult | eigensolves to 99 % | to 99.9 % | decreases | volume residual | stop |
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
summary over starts (reached = within 0.1 % of the best value found from that start):
  ours                 reached 1/5; median eigensolves to reach  169.0; decreases    0; max volume residual 2.2e-16; max bound violation 0.0e+00
  SLSQP bound          reached 3/5; median eigensolves to reach   40.0; decreases  138; max volume residual 1.0e-10; max bound violation 0.0e+00
  trust-constr bound   reached 1/5; median eigensolves to reach  194.0; decreases    9; max volume residual 6.7e-16; max bound violation 0.0e+00
  MMA bound (pyMOTO)   reached 1/5; median eigensolves to reach  199.0; decreases  400; max volume residual 6.4e-05; max bound violation 0.0e+00
  naive fixed-step     reached 0/5; median eigensolves to reach    nan; decreases  435; max volume residual 1.1e-16; max bound violation 0.0e+00
```

Stop codes: `budget` = eigensolve budget exhausted; `converged` / `failed` = solver's own termination. Full solver messages:
- SLSQP bound: Optimization terminated successfully
