# Limitations

These limits are deliberate. They were first recorded for v0.1 in an earlier assessment in the original research code (private archive) and still apply,
except that the certified layer (v0.2) adds deformations near Bolza and the sensitivity module works with any
finite-element model.

## It covers one surface

Everything in bolzakit is about the Bolza surface. Codes, larger lattice experiments and genus-3 work need families of
surfaces (covers, other genera), and bolzakit has none. The Γ(p) golden-geodesic function concerns congruence covers of
the modular surface, not covers of Bolza.

## The FEM is a cross-check only

The bundled finite-element eigenvalues are 9 clusters. They agree with Strohmaier–Uski to about 3e-5 and are labelled
with a conservative relative error of 2e-4. Their job is to check the quoted data independently and to give an offline
answer for the lowest clusters. The real spectral work is carried by the Strohmaier–Uski eigenvalues, which are QUOTED
and downloaded on demand (sha256-checked, never bundled). Without network access or a cached copy, `su2013`,
`eigenvalues("su2013")` and `selberg_check` raise `DataUnavailable`. See [spectral.md](spectral.md).

## There is no tight-binding layer

There are no adjacency or graph-Laplacian Hamiltonians, no finite lattice patches and no Bloch-phase machinery, which
most physics users build. bolzakit gives the continuum data such a layer would be compared against.

## The Green's function is local to the Weierstrass charts

`green().G(z, w)` and `regular(z, w)` are validated only when both points lie within Poincaré radius 0.45 of the same
Weierstrass point, in a chart centred there. Other pairs raise `OutOfValidityRange`; a general-position G_Σ is not
provided. The charts cover the surface (covering radius 0.4056), so every single point is valid and `robin_mass` works
everywhere. The omitted heat-kernel tail limits accuracy to about 4e-6 in value at T = 2.6. Only the validated
parameter pairs (T, D_ball) = (2.6, 15) and (2.6, 13) are accepted; other pairs raise `OutOfValidityRange` (un-validated
pairs can be off by up to 25%).

## Vortex results are model-conditional

- The λ₃ series holds only for N = 7 in the SO(2)-invariant degree-7 truncation of the Boatto–Koiller model.
- In the full Bolza Hamiltonian the ring is only an approximate relative equilibrium, with an O(ξ⁴) residual.
- The Boatto–Koiller model freezes the harmonic degrees of freedom (recorded in the original research code). See [vortex.md](vortex.md).

## Other

- `lengths` is bundled up to half-trace 120. Beyond that, `compute=True` enumerates afresh, at a cost that grows like
  e^ℓ.
- `tau_bounds` data is bundled. Regenerating it takes about 30 minutes and 8 GB of RAM.

## Roadmap

1. **Congruence covers of Bolza.** Reduce the exact ℤ[ζ₈] generators modulo prime ideals 𝔭 of ℤ[ζ₈] and enumerate the
   finite quotients Γ/Γ(𝔭). From each quotient, build the genus-g cover:
   - its group;
   - its tiling as a cell complex;
   - its homological code parameters [[n, 2g, d]], with d from the combinatorial systole.

   This is new work, not a port of the Γ(p) material, which is about SL(2,ℤ).
2. **A tight-binding layer:**
   - adjacency and graph-Laplacian Hamiltonians on the covers and on finite patches;
   - abelian Bloch phases on the genus-2 cell [MR21];
   - comparison of graph spectra with the continuum spectrum.

[MR21] J. Maciejko, S. Rayan, *Hyperbolic band theory*, Sci. Adv. 7, eabe9170 (2021); arXiv:2008.05489.
