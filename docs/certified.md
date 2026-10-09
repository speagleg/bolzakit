# Certified layer: `bolzakit.certified` and `bolzakit_check`

`bolzakit.certified` computes rigorous enclosures (python-flint / Arb ball arithmetic) for the Laplace spectrum of the
Bolza surface and of small deformations of it. Each call returns a result object and a JSON certificate
(schema `bolzakit-cert/1`). The separate package `bolzakit_check` checks those certificates and imports nothing from
`bolzakit`. Install with `pip install "bolzakit[cert]"`.

These are the computations behind [proofs/bolza-local-max/](../proofs/bolza-local-max/). The proof package keeps the
golden certificates and `reproduce.py`.

Contents: [Method](#method) · [API reference](#api-reference) · [The checker](#the-checker) ·
[Trust model](#trust-model) · [Validity limits](#validity-limits) · [Known limitations](#known-limitations)

## Method

- **Surfaces as glued polygons.** A `Surface` is a hyperbolic octagon, given by its 8 hyperboloid vertices as 96-bit
  Arb balls and glued with the Bolza side pairing. A deformed surface is certified to close up exactly: either in
  closed form (F_B), by a 1-D interval Newton step for the angle-sum correction σ(t) (F_A), or by a Krawczyk
  inclusion in the 5-dimensional normal space of the constraint manifold (vector directions and
  `Surface.from_trace_coords`).
- **Conforming, chart-free elements.** The octagon is triangulated into geodesic triangles (`Mesh(n)`). Each triangle
  is pulled back to the reference simplex through the hyperboloid-barycentric map, so no coordinate chart of the
  surface is needed. The P1 and P2 Lagrange elements are conforming across the glued sides. The element integrals
  of the stiffness q and mass m forms (and of their t-derivatives along a deformation, through order-2 jets of the
  vertices) are enclosed by a binomial series with a Lagrange remainder.
- **Rayleigh–Ritz upper bounds** (`certify_spectrum`). Float eigenvectors are only a proposal. Their exact dyadic
  values span a trial space, and min-max gives λ_k ≤ μ_max(K_p, M_p). μ_max < μ is certified by Sylvester's
  criterion for μM_p − K_p in Arb.
- **Trace-formula counting** (`count_eigenvalues`). This is the Fortier Bourque–Petri linear-programming method
  (arXiv:2111.14699 §3, Lemma 3.2) for the Selberg trace formula. A float LP proposes the test function, and its
  coefficients are rounded to rationals. Then exact Sturm sign checks (sympy), Arb integration with a tail bound,
  and the 48 FBP geodesics give a certified upper bound on the number of eigenvalues in a window.
- **Lehmann–Goerisch lower bound** (`certify_lam1_lower`). P2 trial functions and H(div)-conforming RT₂ fluxes
  give λ₁ ≥ 3.83887570105 at n = 48, with shift ρ = 5. The method needs ρ ≤ λ₄ + 1, which is taken from the
  counts.
- **Cluster perturbation and the η transfer** (`certify_cluster`). At Bolza, λ₁ = λ₂ = λ₃. The exact eigenspace E
  is unknown, so each FE trial vector u stands in for its projection e = P_E u. The gap transfer bounds the error:
  with λ₁ in the input ball and λ₄ ≥ the gap input, m(u − e) ≤ (μ − λ₁)/(λ₄ − λ₁)·m(u), and q(u − e) is
  bounded similarly. Elementwise relative bounds β on the derivative forms carry this through every bilinear term.
  - First order (P1): certifies b(e,e)/m < 0 along F_A, so the first-order matrix A is non-zero.
  - Second order (P2): certifies **Ky Fan**: Λ = λ₁ + λ₂ + λ₃ ≤ g(t) = tr(G_m⁻¹G_q) on span{1, e_i + t ψ_i},
    with equality at t = 0. Hence D²Λ(v,v) ≤ g''(0), and the certificate bounds g''(0) above.
- **Trace chart** (`certify_rank`). Forward-mode dual numbers over Arb give the 6×16 Jacobian of six trace
  functions and the 5×16 constraint Jacobian. A non-zero 6×6 minor on ker of the constraints certifies that the
  trace map has rank 6 at Bolza (a local chart of T₂).

## API reference

All functions set their Arb precision internally and restore it afterwards: 96 bits by default (Rayleigh–Ritz,
cluster tests and `Surface` construction), 128 bits for `certify_rank`, `certify_lam1_lower` and `from_trace_coords`,
and 200 bits for `count_eigenvalues`.
With `strict=True`, a computation that does not separate its enclosures raises `CertificationFailed`. Otherwise
the result has `certified=False`. Every result carries its certificate (`result.certificate`, or
`result.certificates` for clusters). Write it with `cert.save(path)`.

| Function | Returns | Notes |
|---|---|---|
| `Surface.bolza()` | the Bolza octagon | `is_bolza=True`, `r2_invariant=True` |
| `Surface.deform(direction, t)` | a closed deformed octagon | `direction` is `"F_A"`, `"F_B"` or a 16-vector of Klein (x, y) vertex coordinates. Raises `OutOfValidityRange` for \|t\| > 0.05 |
| `Surface.from_trace_coords(tau)` | the octagon with trace coordinates τ (6 floats) | Krawczyk at 128 bits. Raises `OutOfValidityRange` if \|τ − τ_B\|∞ > 0.05·max τ_B |
| `Surface.vertex_jets(direction)`, `.klein_tangent(direction)`, `.klein_coords()` | vertex jets, tangent, coordinates | used by the cluster and chart code |
| `certify_spectrum(surface, k=3, degree=2, n=48)` | `SpectrumResult` | `upper[i] ≥ λ_{i+1}`, attested by the certificate. `upper_inprocess` holds the in-memory bisection values, which the published goldens pin. `cluster_sum_upper` (degree 2) bounds λ₁+…+λ_k |
| `count_eigenvalues("bolza", window=(0, 3.3))` | `CountResult` | `bound` is an upper bound on the count. `window=(0, a)`: mode A, count in (0, a], target 1, assumption-free. `window=(a, b)`: mode B, count in [a, b], target 4, conditional on λ₁ ≥ a |
| `certify_lam1_lower(n=48, rho=5, gap=GapInput(...))` | `LowerBoundResult` | `lam1_lower ≤ λ₁(Bolza)`, conditional on the gap input. About 4 min |
| `certify_cluster(direction, degree, n, lam1=Lam1Input(...), gap=GapInput(...), second_order=False)` | `ClusterResult` | degree 1 (F_A only): `b_ee_upper`, `A_nonzero`. degree 2, `second_order=True`: `kyfan_upper ≥ D²Λ(v,v)`. degree 2, `second_order=False`: float splitting rates `A_float` only, with no certificate. The `*_inprocess` fields keep the in-memory values |
| `certify_rank(surface=None)` | `RankResult` | the rank of the trace map is 6 at the surface (default Bolza) |
| `Lam1Input(lower, upper, source, dep_id="lam1")`, `GapInput(lower, source)` | explicit inputs | recorded as certificate dependencies. `dep_id` names the λ₁ dependency: a PROOF.md §Dependencies number only when the source is that input (`"D1"` for SU13, `"D15"` for the Lehmann–Goerisch route), otherwise the unnumbered `"lam1"` |

`lam1` and `gap` are always explicit. bolzakit never assumes a λ₁ ball or a gap. The values used in the proof are
`Lam1Input("3.8388757010", "3.8388872830", ...)` (the Lehmann–Goerisch lower bound and the P2 Rayleigh–Ritz upper
bound) and `GapInput("4.0", ...)` (the trace-formula counts).

```python
from bolzakit.certified import Surface, certify_spectrum, count_eigenvalues

r = certify_spectrum(Surface.bolza().deform("F_A", 0.02), k=3, degree=2, n=12)
print(r.certified, r.upper)          # True [3.7445..., 3.8416..., 3.9368...]
print(count_eigenvalues("bolza", window=(0, 3.3)).bound)   # 0.9126... < 1: no eigenvalue in (0, 3.3]
```

## The checker

```
python -m bolzakit_check FILE.json... [--deep]
```

The checker prints `VALID`, `VALID given stored Gram data`, `VALID (conditional on …)` or `INVALID: <reason>` for
each file. Exit codes:

- 0: every file is unconditionally VALID;
- 2: none is INVALID, but some is only conditionally VALID;
- 1: some file is INVALID.

On the golden set it exits **2**, as expected, with six conditional files: `lp_B_3.3_4.0` is conditional on
λ₁ ≥ 3.3; `lg_lower_n48` on λ₄ ≥ 4.0; and the four first-order and Ky Fan certificates on their λ₁ ball and on
λ₄ ≥ 4.0. For these the message names the λ₁ ball, its dependency id and source, and the gap with its source.
Never read "VALID" as "proved end to end". See the trust model below.

### Certificate fields by kind

Schema `bolzakit-cert/1` is shared, but the attested quantity and the failure sentinel are named per kind (to be
unified in a later schema version):

| Kind | Attested quantity (in `data`) | On failure |
|---|---|---|
| `rr_upper` | `mu`, a rational string (λ_k ≤ mu) | `mu = "inf"` |
| `count` | `bound` (Arb ball string) and `bound_upper`, a float rounded up from it | `bound_upper` is `"nan"` |
| `cluster_first_order` | `b_ee_upper_over_m` (Arb ball string) | `"nan"` |
| `kyfan_second_order` | `g2_enclosure` (Arb ball string) | `"nan"` |
| `lg_lower` | `lam1_lower_enclosure` (Arb ball string); claim rhs is a rational string | `"nan"`; `A0`/`A1`/`A2` are `null` for non-finite element data, `nubar` is `null` when no ν̄ was found |
| `rank` | `J` (the 6×6 minor), `stacked_det` (Arb ball strings) | `stacked_det` is `null` |

Only `rank` records its working precision (`inputs["prec"]`).

## Trust model

**Trusted** (not checked by anything in this repository):
- python-flint / Arb (ball arithmetic, `acb.integral`) and CPython, in both the producer and the checker; sympy
  (exact Sturm counts, exact Q(ζ₂₄) arithmetic) in the producer only. The checker does its own Sturm counts over
  `Fraction` and never imports sympy;
- the cited inputs D1–D15 of [PROOF.md §Dependencies](../proofs/bolza-local-max/PROOF.md). These include the
  Selberg trace formula in FBP's normalisation, the FBP geodesic data and Lemmas 7.2–7.3, and the Lehmann–Goerisch
  theorem, quoted from secondary sources. D1 (Strohmaier–Uski) is optional.

**Checked by `bolzakit_check`.** The checker re-derives each claim from the certificate's inputs and data, in its
own Arb code. A tampered claim fails even if the sha256 is recomputed.

| Kind | Re-derived | Trusted |
|---|---|---|
| `count` | everything: Sturm sign conditions, the window minimum, the integral and the geodesic sum, from the rational coefficients and closed-form lengths | the trace formula and FBP geodesic data; mode B's hypothesis λ₁ ≥ a |
| `rank` | the Jacobians (own dual numbers), the minors, the stacked determinant and the residual, from the 16 Klein balls and the words; the label "bolza" only if the Klein balls contain the exact Bolza octagon (otherwise reported as unchecked) | nothing beyond Arb |
| `rr_upper` | M ≻ 0 and μM − K ≻ 0; degree 2 also tr(M⁻¹K) ≤ the cluster-sum claim (Ky Fan on the trial space) | the stored Gram balls K, M |
| `cluster_first_order` | μ, B, E and B + E < 0, from the raw element sums | the raw sums; the λ₁ / gap inputs (always conditional) |
| `kyfan_second_order` | g₂ from the raw 7×7 Gram matrices, weighted sums and β bounds | the raw Gram data; the λ₁ / gap inputs (always conditional) |
| `lg_lower` | B ≻ 0, ν* < 0 and ν*B − A ≻ 0, from A0, A1, A2 and ρ (a Lehmann–Goerisch step written from the theorem) | the Gram balls, including the H(div) conformity of the RT₂ fluxes; the gap input |

**Not checked:**
- **Element assembly.** The checker does not re-integrate the finite-element forms over the mesh. The stored Gram
  matrices and raw sums are trusted to be what the producer computed. `--deep` re-integrates a sample of elements
  (1% of the elements, at most 200) with mpmath quadrature and compares them to 1e-8. That is a **consistency check
  only**, not a proof.
- **The links between certificates.** Several conditions are linked by hand, not by a tool (see
  [Known limitations](#known-limitations)).
- **The non-certificate steps of the proof** (written arguments and exact representation theory). The Lean skeleton
  (`lean/BolzaKit/LocalMax`) checks the logical chain, conditional on the computational inputs. Its CI axiom gate
  allows only `propext`, `Classical.choice` and `Quot.sound`.

## Validity limits

- **Deformations: |t| ≤ 0.05, near Bolza only.** `Surface.deform` raises `OutOfValidityRange` beyond that, and
  `from_trace_coords` has the matching limit in trace coordinates. For vector directions the limit applies to t,
  so it depends on the scale of v. Pass unit vectors.
- **Counting is Bolza-only.** `count_eigenvalues` uses the FBP geodesic data of the Bolza surface. Any other
  `surface` raises.
- **Cluster tests are at Bolza.** `certify_cluster` perturbs the isolated cluster of three eigenvalues at Bolza
  (the isolation and the λ₁..λ₃ window are certified; multiplicity 3 is quoted, D4). Its
  degree-1 test is implemented for F_A only.
- **Lower bounds come only from explicit inputs.** The one in-house lower bound is `certify_lam1_lower`, for Bolza,
  and it is itself conditional on a gap input. bolzakit never invents a λ₁ ball or a gap.

## Known limitations

- **Deformation directions near Bolza only** (|t| ≤ 0.05; see above).
- **Counting is Bolza-only.**
- **Arbitrary vector directions certify weakly.** For a 16-vector direction, `Surface.vertex_jets` uses the
  minimum-norm second-order jet w (solving DC w = −D²C(v,v)). This is a valid curve, because dΛ = 0 at Bolza makes
  D²Λ(v,v) independent of the curve. But it bends the vertices about 5× as much as the named F_B family, and the
  β error terms grow with it. F_B as a vector certifies +3.47 (not certified) at n = 24, while the named F_B gives
  −9.71. A better curve choice is deferred to v0.3.
- **Certified runs are slow.** Indicative single-core times are about 4 s for P1 Rayleigh–Ritz at n = 12, 7 min at
  n = 192, and 1 min for P2 at n = 96. A trace-formula count takes 4–20 s, the n = 48 Lehmann–Goerisch bound
  about 4 min, the n = 192 first-order test about 11 min, and the n = 96 Ky Fan bound 2–3 min. Regenerating the
  golden set (`reproduce.py --full`) takes 35–50 min. The default test suite therefore uses n ≤ 24, and the
  golden pins run under `--runslow`.
- **The checker trusts Gram / raw-sum data.** Assembly is outside its trust boundary. `--deep` is a consistency
  check only.
- **The kyfan and b_ee checker formulas share their derivation with the producer.** They are independent code but
  transliterate the same derivation (the g2_upper_v1 and b_ee_upper_v1 assembly formulas documented in
  `bolzakit.certified.cluster`). A derivation error would be reproduced in both. Only the `count` and `lg_lower`
  re-derivations are derivation-independent: the checker works from the FBP lemma and the Lehmann–Goerisch theorem
  statement, not from the producer. `rr_upper` is a direct definiteness check. `rank` uses independent Jacobian code
  on the same trace-word formulation.
- **Hand-linked conditions.** A person must confirm these links. The checker takes each condition as a given input:
  - λ₁ ≥ a, the hypothesis of a mode-B count. For `lp_B_3.3_4.0`, `lp_A_3.3` discharges it.
  - λ₄ ≥ 4.0, the gap input of the Lehmann–Goerisch, Ky Fan and first-order certificates. `lp_A_3.3` with
    `lp_B_3.3_4.0` discharges it. Confirm that the gap in the checker message is 4.0.
  - The λ₁ ball [3.8388757010, 3.8388872830] of the cluster certificates. `lg_lower_n48` and `rr_p2_n96` discharge
    it. Confirm that the checker message of each `_lg` certificate names this ball with dependency id D15.
- **The LG computation cannot detect a wrong ρ.** ρ = 7.5 violates the hypothesis and still returns a plausible
  value below the true λ₁ (review_lb_1 row 4c). The validity of ρ = 5 rests entirely on the gap λ₄ ≥ 4.0. The
  checker checks ρ ≤ gap + 1 against the recorded gap input, and nothing else.
- **No human expert review.** The derivations and reviews behind the certified layer were produced by AI agents.
  See [VERIFICATION.md](../VERIFICATION.md).
