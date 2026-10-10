# Sensitivity of repeated eigenvalues: `bolzakit.sensitivity`

`bolzakit.sensitivity` computes derivatives of an eigenvalue **cluster** of a symmetric-definite pencil
K(t) u = λ M(t) u, the kind any finite-element code produces. It works with dense (numpy) or sparse (scipy.sparse)
matrices, it imports only numpy and scipy (nothing else from bolzakit), and all of its results are labelled
NUMERICAL. Adapters assemble the matrices from scikit-fem or FEniCSx forms.

Contents: [Why clusters](#why-individual-repeated-eigenvalues-are-not-differentiable-while-cluster-sums-are) ·
[Formulas](#formulas) · [Near-degenerate clusters](#near-degenerate-clusters) ·
[Discrete vs continuum](#the-discrete-vs-continuum-caveat) · [Membrane example](#the-membrane-example) ·
[Adapters](#adapters) · [0.3.0 additions](#splitting-matrices-optimisation-and-trimming-030) · [API](#api)

## Why individual repeated eigenvalues are not differentiable, while cluster sums are

Suppose λ has multiplicity m at p₀, and the pencil depends smoothly on several parameters p. Near p₀ the m branches
λ₁(p) ≤ … ≤ λ_m(p) split with directional derivatives that depend on the direction v: they are the eigenvalues of
an m×m matrix A(v), linear in v. Sorted eigenvalues of a matrix that depends linearly on v are not linear in v, so
the individual branches have no gradient at p₀. Sorting the pair of branches ±p gives −|p| and |p|, and in two
parameters the branches form the cone λ = ±√(p₁² + p₂²). Along a single analytic one-parameter path the branches can be relabelled to be analytic
(Rellich), but which branch is "λ₁" changes with the path.

The **sum** of the eigenvalues in an isolated cluster is the trace of the pencil restricted to the cluster's spectral
subspace. It is a smooth (analytic) function of p, as long as the cluster stays separated from the rest of the
spectrum. So its gradient tr A and its Hessian are well defined. The same holds for the sum of all eigenvalues up to
the top of a cluster (the prefix sum). Its second derivative has a Ky Fan upper bound, which is what an optimiser
that pushes a repeated eigenvalue up needs.

## Formulas

`cluster_sensitivity(K, M, dK, dM, cluster=(i, j), d2K=None, d2M=None)` uses the M-orthonormal cluster
eigenvectors U (n × m) and the cluster mean λ̄.

- First order: A = Uᵀ(dK − λ̄ dM)U. `rates = eigvalsh(A)` are the splitting rates of an exactly degenerate cluster.
  `tr A` is the first derivative of the cluster sum.
- Cluster-sum Hessian (`trace_hessian`): with B = (dK − λ̄ dM)U and Z solving (K − λ̄M)Z = B − MUUᵀB, UᵀMZ = 0,

      d²/dt² Σ_cluster λ = tr(Uᵀ(d²K − λ̄ d²M)U) − 2 tr(A UᵀdM U) − 2 tr(BᵀZ).

  The constrained system is solved as a bordered system [[K − λ̄M, MU], [UᵀM, 0]].
- **Ky Fan prefix bound** (`kyfan_upper`). Let W be all eigenvectors 0..j (the prefix up to the top of the cluster)
  and Ψ a trial direction. By default Ψ consists of the exact first-order correctors; pass `trial=` to use another.
  By Ky Fan, Σ_{k≤j} λ_k(t) ≤ g(t) = tr[(W_tᵀM(t)W_t)⁻¹ W_tᵀK(t)W_t] with W_t = W + tΨ, and equality holds at
  t = 0. So g''(0) ≥ d²/dt² Σ_{k≤j} λ_k. With G = W_tᵀKW_t and N = W_tᵀMW_t (N(0) = I):

      g''(0) = tr G'' − 2 tr(N'G') + tr((2N'N' − N'')G₀),
      G'  = ΨᵀKW + WᵀdK W + WᵀKΨ,               N'  = ΨᵀMW + WᵀdM W + WᵀMΨ,
      G'' = Wᵀd²K W + 2(ΨᵀdK W + WᵀdK Ψ) + 2ΨᵀKΨ,  N'' = Wᵀd²M W + 2(ΨᵀdM W + WᵀdM Ψ) + 2ΨᵀMΨ.

  This is the exact g''(0). It replaces the formula in the v0.2 plan, which dropped terms that vanish only when the
  prefix equals the cluster, and so was not an upper bound (Task 12 ruling). With the default correctors the bound
  is tight. It stays valid for any trial.
- A gap check guards every call. `cluster_sensitivity` raises `ValueError` if the cluster splits a
  (near-)degenerate group, that is, if g_out < max(gap_ratio_min·|λ̄|, 10·g_in). Here g_in is the spread inside the
  cluster and g_out the distance to the nearest eigenvalue outside it.

## Near-degenerate clusters

If the cluster is only **nearly** degenerate (g_in > 0), `rates` are still, by definition, eig(A) at λ̄. They are
**not** O(g_in)-accurate predictions of the derivatives of the individual branches. Counterexample: take
K = diag(1, 1+δ, 3), M = I, and dK with off-diagonal entries 1 coupling the first two modes. For every δ > 0,
`rates` = ±1, while both branch derivatives at t = 0 are exactly 0 (the diagonal of dK). The rates describe the
splitting once |t| ≫ g_in/|A|, not the derivatives at t = 0. `tests/sensitivity/test_core.py` pins this example.

The cluster **sum** is smooth, so `tr A` and `trace_hessian` are correct: exact for an exactly degenerate cluster,
and accurate to O(g_in) otherwise (exact in the counterexample above). Use the sums, not the individual rates, when
the cluster is not exactly degenerate.

## The discrete-vs-continuum caveat

Every quantity here belongs to the **discrete** pencil (K, M) that you pass in. The rates, the Hessian and the Ky
Fan bound are statements about the eigenvalues of the matrices, not about the eigenvalues of the continuous
operator that the finite-element model discretises. The discretisation error is not controlled, and nothing here is
certified. That is why the label is NUMERICAL. For certified statements about a continuous operator, see
[docs/certified.md](certified.md), which handles the Bolza surface this way, with the transfer from discrete to
exact eigenvectors made rigorous.

## The membrane example

[`examples/membrane_repeated_frequency.py`](../examples/membrane_repeated_frequency.py) stretches a unit-square
Dirichlet membrane to a × 1, with a = 1 + p. At p = 0 the modes (1,2) and (2,1) share λ = 5π². The example
assembles the pencil with the scikit-fem adapter (Q2 elements) and compares the results with the exact values
λ_mn(a) = π²(m²/a² + n²):

| quantity | exact | computed (rel. error) |
|---|---|---|
| rate of (1,2) | −2π² | 1.3e-7 |
| rate of (2,1) | −8π² | 2.1e-6 |
| cluster-sum second derivative | 30π² | 9.8e-5 |
| Ky Fan bound for λ₁₁ + λ₁₂ + λ₂₁ | 36π² | 1.2e-4 |

The residual errors are the discretisation error and the O(h²) finite differences of the adapter.

```
python examples/membrane_repeated_frequency.py
```

## Adapters

Both adapters return `(K, M, dK, dM, d2K, d2M)` as sparse matrices, ready for `cluster_sensitivity`, with the
Dirichlet dofs removed.

- **scikit-fem** (`bolzakit.sensitivity.adapters.skfem.from_forms(basis, stiffness_form, mass_form, p0, h=1e-5,
  dirichlet_dofs=None)`; `pip install "bolzakit[skfem]"`). The forms read the parameter from `w["p"]`. The
  derivatives are central finite differences of the assembled matrices: O(h²) truncation error, plus O(eps/h²)
  round-off in the second derivatives.
- **FEniCSx** (`bolzakit.sensitivity.adapters.fenicsx.from_forms(V, a, m, p, bcs=None)`). The parameter must be
  `p = ufl.variable(dolfinx.fem.Constant(mesh, value))`, and the derivative forms are `ufl.diff`, so the
  derivatives are exact, with no finite differencing. This adapter runs on one MPI rank only. dolfinx is not on
  PyPI: install it with conda or use the `dolfinx/dolfinx` docker image. The weekly `fenicsx.yml` workflow runs its
  test in that image.
- **Your own code.** `from_callables(Kfun, Mfun, p0, v, cluster, h=1e-4)` takes the derivatives by central finite
  differences of any K(p), M(p) along v and emits a `UserWarning`. `directional(...)` takes callables that return
  dK(v), dM(v), and so on.

## Splitting matrices, optimisation and trimming (0.3.0)

- `smallest_eigenpairs(K, M, count)`, `splitting_matrices(K, M, dKs, dMs, cluster, ...)` (one `A_l` per design
  parameter, so `A(v) = Σ v_l A_l` for any direction v) and `ritz_upper_bounds(K, M, W)` (Poincaré separation upper
  bounds on the discrete eigenvalues) extend the core to many design parameters.
- [Optimisation](sensitivity-optimize.md): `ascent_direction` and `maximize_eigenvalue` (steepest ascent of `λ_k` or
  of a gap at a multiple eigenvalue), with a benchmark against the standard bound formulation (which SLSQP won on
  success rate and eigensolve count).
- [Trimming](sensitivity-trim.md): `trim_degenerate_pair`, the inverse problem of tying a split mode pair with the
  cheapest trim, validated in `examples/gyro_trim.py` on a ring and on a ring with 8 spokes.

## API

| Function / field | Meaning |
|---|---|
| `cluster_sensitivity(K, M, dK, dM, cluster, d2K=None, d2M=None, trial=None, gap_ratio_min=1e-3)` | returns a `ClusterSensitivity`. `cluster=(i, j)` is inclusive and 0-based. Without d2K/d2M, `trace_hessian` and `kyfan_upper` are None |
| `.eigvals`, `.basis` | cluster eigenvalues (ascending) and M-orthonormal eigenvectors |
| `.A`, `.rates` | the first-order matrix and eig(A) at λ̄ |
| `.trace_hessian` | d²/dt² of the cluster sum |
| `.kyfan_upper` | an upper bound on d²/dt² of the prefix sum λ₀ + … + λ_j (discrete problem) |
| `.label` | always `"NUMERICAL"` |
| `directional`, `from_callables` | see [Adapters](#adapters) |

The sparse path uses shift-invert `eigsh` with an inertia-checked shift, and `spsolve`. Errors: `ValueError` if K
or M is not symmetric, if M is not positive definite, or if the gap check fails.
