# Changelog

## [0.3.0] — unreleased

### Added

- `bolzakit.sensitivity.errors`: typed, stdlib-only sensitivity errors (`SensitivityError`, with `PencilError`,
  `ClusterRangeError`, `ClusterSplitError`, `NonSymmetricFunctionError` and `MissingExtra`) and
  `FiniteDifferenceWarning`. The `ValueError` subclasses keep existing `except ValueError` handlers working.
- `splitting_matrices`: the per-parameter splitting matrices A_l = U^T (dK_l - lam_bar dM_l) U of a cluster, so that
  A(v) = sum_l v_l A_l for any direction v.
- `ritz_upper_bounds`: Ritz values of the pencil compressed to a trial subspace (Poincaré separation upper bounds on
  the discrete eigenvalues, and the Ky Fan trial bound on their sum).
- `bolzakit.sensitivity.optimize.ascent_direction`: the steepest-ascent unit direction for lambda_min(A_upper(d))
  (optionally minus lambda_max(A_lower(d)), for a gap) from splitting-matrix stacks, with optional null-space and
  bound constraints and a certified duality-gap bound on the optimal rate.

## [0.2.0] - 2026-10-09

### Added

- `bolzakit.certified` (extra `[cert]`: python-flint, sympy). This is a certified (Arb ball arithmetic) layer for
  the Bolza surface and its small deformations, and every result carries a `bolzakit-cert/1` JSON certificate.
  - `Surface`: the Bolza octagon, the deformations `deform("F_A" | "F_B" | 16-vector, t)` for |t| ≤ 0.05, and
    `from_trace_coords`. Closure is certified by interval Newton or Krawczyk.
  - `certify_spectrum`: P1 / P2 Rayleigh–Ritz upper bounds. `upper` is attested by the certificate, and
    `upper_inprocess` keeps the in-memory values.
  - `count_eigenvalues`: Selberg trace-formula counts for Bolza (Fortier Bourque–Petri LP). Mode A is
    assumption-free. Mode B is conditional on λ₁ ≥ a.
  - `certify_lam1_lower`: the Lehmann–Goerisch lower bound λ₁ ≥ 3.83887570105 (n = 48, ρ = 5), conditional on
    λ₄ ≥ 4.0.
  - `certify_cluster`: the first-order test (A ≠ 0 along F_A) and the second-order Ky Fan bound on D²Λ(v,v), with
    explicit `Lam1Input` / `GapInput` and certificate-attested `b_ee_upper` / `kyfan_upper` (plus `*_inprocess`).
  - `certify_rank`: rank 6 of the trace-coordinate chart at Bolza.
- `bolzakit_check`, an independent checker (`python -m bolzakit_check FILE... [--deep]`) that imports nothing from
  `bolzakit`. It re-derives each claim per certificate kind. Exit codes: 0 = all VALID, 2 = some only conditionally
  VALID (mode-B counts, lg_lower, and the first-order / Ky Fan certificates, which rest on hand-linked λ₁ and gap
  inputs), 1 = some INVALID. VALID messages state their trust boundary ("VALID given stored Gram data") and their
  conditions. The rank checker calls the surface "bolza" only after checking the Klein balls against the exact
  Bolza octagon.
- `bolzakit.sensitivity`, a standalone module that depends only on numpy and scipy. `cluster_sensitivity`
  computes the splitting rates, the cluster-sum Hessian and an exact Ky Fan prefix bound for a symmetric pencil.
  It also provides `directional` and `from_callables`, and scikit-fem and FEniCSx adapters.
- `proofs/bolza-local-max/`: the computer-assisted proof that Bolza is a strict local maximum of λ₁, with ten
  golden certificates (an SU13-free route, plus SU13 cross-checks), `reproduce.py --quick | --full`, the paper and
  the review record. The status is conditional on the cited inputs, and no human expert has read the proof.
- Lean 4: the `BolzaKit/LocalMax` skeleton (8 files) of the local-maximum argument, conditional on its computational
  inputs, with a CI axiom gate (`scripts/check_axioms.sh`).
- Docs: `docs/certified.md` (method, trust model, validity limits, known limitations) and `docs/sensitivity.md`.
  The README has new sections for software engineers, physicists and engineers.
- Examples: `bolza_deform.py`, `bolza_strain_rates.py` and `membrane_repeated_frequency.py`.
- CI: the default job installs `.[dev,cert,skfem]` and runs `bolzakit_check` on the golden certificates (exit 0
  or 2 passes). The weekly job runs `--runslow` and `reproduce.py --quick`. A weekly FEniCSx job runs in the dolfinx
  image.

### Changed

- The scikit-fem example is named `membrane_repeated_frequency.py` (it is a membrane, not a plate).
- The `bolzakit.sensitivity.core` docstring no longer claims that near-tie splitting rates are O(spread)-accurate.
  They are eig(A) at the cluster mean. The cluster-sum derivatives are the reliable quantities.

### Fixed

- `scripts/check_axioms.sh` failed open when Lean wrapped a `#print axioms` list across lines. It now joins the log
  before parsing, and it fails on an unterminated list. The Lean CI axiom step runs under `shell: bash` (pipefail),
  so a failing `lean` is no longer masked by `tee`.

### Performance (measured, one core; not asserted in tests)

- Regenerating the golden set (`reproduce.py --full`) takes about 36 min. Per certificate: P1 n = 192 Rayleigh–Ritz
  419 s, P2 n = 96 65 s, first order n = 192 about 680 s each, Ky Fan n = 96 125–150 s, the LP counts 4–21 s,
  rank 1 s, Lehmann–Goerisch n = 48 about 4 min. `reproduce.py --quick` takes about 1–2 min (79 s measured).

### Known limitations

- Deformations are validated near Bolza only (|t| ≤ 0.05), and counting is for Bolza only.
- Arbitrary vector directions certify weakly, because of the minimum-norm second-order jet. This is deferred to
  v0.3.
- The checker trusts the stored Gram / raw-sum data, and `--deep` is a consistency check only. The kyfan and b_ee
  checker formulas share their derivation with the producer. λ₁ ≥ a and λ₄ ≥ 4.0 are linked by hand. The
  Lehmann–Goerisch computation cannot detect a wrong ρ. See `docs/certified.md`.
- No human expert has reviewed the proof or the derivations (see `VERIFICATION.md`).

## 0.1.0

First release.

- `bolzakit.exact`: exact arithmetic in Z[zeta_8] (`zeta8`), the Bolza Fuchsian group with generators, relator,
  enumeration and traces (`group`), the automorphism group of order 48 and the Weierstrass points (`automorphisms`) and
  the length spectrum (`lengths`).
- `bolzakit.geometry`: the regular hyperbolic octagon and its side pairings (`octagon`).
- `bolzakit.spectral`: Laplace eigenvalues (bundled FEM cross-check clusters; Strohmaier-Uski data downloaded on demand,
  sha256-checked), Selberg trace check, ζ_Δ(0) and the tau data.
- `bolzakit.green`: Green's function, regular part and Robin mass in Weierstrass charts, Bergman-kernel data, heat
  kernel and the Green's-function constants, computed in closed form (`data/green_constants.json` is a reference copy,
  checked against the closed form by a regression test).
- `bolzakit.vortex`: linearised ring-mode stability of regular vortex rings (flat, sinh, tanh and Bolza models) and the
  λ₃ series for the Boatto-Koiller model.
- `bolzakit.arith`: the arithmetic facts: curve data (`curve`, optional cypari2), the arithmetic core of the
  golden-geodesic criterion (`golden`) and the Salem-norm identity (`salem`).
- `bolzakit.provenance` and `bolzakit.errors`: provenance labels (EXACT, PROVED, CERTIFIED, NUMERICAL, QUOTED), available
  with `with_provenance=True`, and the package exceptions.
- Bundled data: `data/*.json` with a `MANIFEST.json`.
- Lean 4: three proofs under `lean/` (the group relation, the arithmetic core of the golden-geodesic criterion and the
  Salem-norm identity), with no `sorry`, `admit`, `native_decide` or `axiom` declarations.
- Docs: README, per-module pages, provenance table, `CITATIONS_CHECK.md`, `limitations.md`, `VERIFICATION.md`,
  `CITATION.cff`.

### Known limitations

- It covers one surface, the Bolza surface. There are no covers or other genera, hence no codes or larger lattices.
- The FEM eigenvalues are a cross-check only: 9 clusters, agreement with Strohmaier-Uski to about 3e-5, labelled with a
  2e-4 relative error. The Strohmaier-Uski data is not bundled.
- There is no tight-binding layer (adjacency or graph-Laplacian Hamiltonians, lattice patches, Bloch phases).
- The Green's function G(z, w) is valid only within Poincare radius 0.45 of the same Weierstrass point; other pairs
  raise `OutOfValidityRange`. `robin_mass(z)` works everywhere. Only (T, D_ball) = (2.6, 15) or (2.6, 13) is accepted.
- `selberg_check(t)` is validated for about 0.03 ≤ t ≤ 1 (bundled data, T_max = 80) and raises outside it.
- The lambda_3 series holds only for N = 7 in the SO(2)-invariant degree-7 truncation of the Boatto-Koiller model.
- No human or peer review yet (see `VERIFICATION.md`).
