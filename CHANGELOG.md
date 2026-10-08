# Changelog

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
  raise `OutOfValidityRange`. `robin_mass(z)` works everywhere.
- The lambda_3 series holds only for N = 7 in the SO(2)-invariant degree-7 truncation of the Boatto-Koiller model.
- No human or peer review yet (see `VERIFICATION.md`).
