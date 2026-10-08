# bolzakit

Exact and provenance-labelled computations on the Bolza surface, the most symmetric closed hyperbolic surface of genus 2.

## Introduction

The Bolza surface is the closed hyperbolic surface (curvature −1, area 4π) of genus 2 with the most symmetry. You get it by
gluing opposite sides of a regular hyperbolic octagon with interior angles π/4, so eight octagons meet at every vertex.
Its orientation-preserving automorphism group is GL(2,3), of order 48, or 96 if you include reflections. Its systole
(shortest closed geodesic) is 2·arccosh(1+√2) ≈ 3.0571, the largest in genus 2 [SU13]. In the physics literature it is
also called the regular octagon or the Hadamard–Gutzwiller model [SU13].

Hyperbolic lattices are now built in the lab. Networks of coplanar-waveguide resonators in circuit QED form lattices in an
effective hyperbolic space [KFH19], and electric-circuit networks have been used to measure the eigenstates of a
"hyperbolic drum" [L+22]. Hyperbolic band theory describes such lattices with Bloch-like phases under the
non-commutative Fuchsian group of a higher-genus surface [MR21]. Hyperbolic surface codes trade encoding rate against
protection more favourably than the surface code [BT16]. For all of these, the Bolza surface is the smallest
non-trivial example.

bolzakit gives you the exact Fuchsian group (matrix entries in ℤ[ζ₈]), the automorphism group, the length spectrum,
Laplace eigenvalues, Green's-function and Bergman-kernel data at the Weierstrass points, vortex-ring stability data and a
few arithmetic facts. It ships with tests and three Lean 4 proofs (the group relation, the golden-geodesic criterion and
the Salem-norm identity). Every value bolzakit returns carries a provenance label (EXACT, PROVED, CERTIFIED, NUMERICAL or QUOTED) — ask for it with `with_provenance=True`.

## Contents

- [Install](#install)
- [Quick-start](#quick-start)
- [Use cases](#use-cases)
  - [(a) Hyperbolic band theory](#a-hyperbolic-band-theory)
  - [(b) Circuit-QED and electric-circuit lattices](#b-circuit-qed-and-electric-circuit-lattices)
  - [(c) Hyperbolic codes](#c-hyperbolic-codes)
  - [(d) Spectral geometry and quantum chaos](#d-spectral-geometry-and-quantum-chaos)
  - [(e) Vortex dynamics](#e-vortex-dynamics)
  - [(f) Arithmetic](#f-arithmetic)
  - [(g) Teaching](#g-teaching)
- [Provenance labels](#provenance-labels)
- [Results and corrections](#results-and-corrections)
- [Limitations](#limitations)
- [Roadmap (v0.2)](#roadmap-v02)
- [Documentation](#documentation)
- [References](#references)
- [Citing](#citing)
- [Licence](#licence)
- [AI-assisted development](#ai-assisted-development)

## Install

```
pip install bolzakit                 # core: numpy, scipy, mpmath
pip install "bolzakit[pari]"         # curve arithmetic via cypari2
pip install "bolzakit[sym,flint]"    # sympy and python-flint extras
```

You need Python 3.10 or later. The Strohmaier–Uski eigenvalue file is not bundled. It is downloaded from arXiv on first
use, checked against a fixed sha256 and cached in `~/.cache/bolzakit` (or `$BOLZAKIT_CACHE`).

## Quick-start

```python
import bolzakit as bk

G = bk.group()
print(len(G.GENS))                          # 8 side-pairing generators
print(bk.automorphisms().order)             # 48 = |GL(2,3)|
print(bk.lengths(max_half_trace=5)[0])      # the systole class: half-trace 1+√2, 24 oriented geodesics
print(bk.systole())                         # 3.057141838961996
print(bk.systole(with_provenance=True).label)   # EXACT
```

## Use cases

Every snippet below runs as written: `tests/test_readme_snippets.py` executes each one. Blocks marked
`# requires: network` need the Strohmaier–Uski file (downloaded once). Blocks marked `# requires: slow` take from
seconds to a few minutes. Blocks marked `# requires: pari` need `bolzakit[pari]`.

### (a) Hyperbolic band theory

Bloch theory on a hyperbolic lattice starts from the Fuchsian group. Its generators play the role of lattice
translations, and the automorphism group plays the role of the point group. bolzakit gives the generators exactly and
the automorphism group with an explicit isomorphism to GL(2,3), so you can build symmetry-adapted bases and read off
multiplicities.

```python
import bolzakit as bk

G = bk.group()
print(G.to_matrix(G.GENS[0]))               # SU(1,1) matrix of g0 (translation length = systole)
print(G.half_trace(G.GENS[0]))              # (1, 1), i.e. 1 + 1*sqrt(2), exact
assert G.word(G.RELATOR) == G.word([])      # g0 g1^-1 g2 g3^-1 g0^-1 g1 g2^-1 g3 = 1, exactly
A = bk.automorphisms()
phi = A.isomorphism_to_GL23()               # element index -> 2x2 matrix over F_3
print(A.order, A.full_order, A.element_order_counts())   # 48 96 {1: 1, 2: 13, 3: 8, 4: 6, 6: 8, 8: 12}
```

Outputs: the `group()` arithmetic is exact in ℤ[ζ₈] (its helper functions carry no label of their own; `elements` and
`systole` are EXACT). `automorphisms()` is PROVED.

### (b) Circuit-QED and electric-circuit lattices

A finite lattice samples the Laplacian of the continuum surface, so its low-lying graph spectrum should approach the
continuum spectrum. bolzakit ships 9 finite-element eigenvalue clusters as a quick offline reference. It also fetches
the 1000 Strohmaier–Uski eigenvalues on demand.

```python
import bolzakit as bk

for lam, mult in bk.spectral.eigenvalues("fem"):          # 9 clusters, bundled
    print(f"{lam:9.4f}  x{mult}")
print(bk.provenance.info(bk.spectral.eigenvalues).error)  # FEM: relative error <= 2e-4 vs Strohmaier–Uski
```

```python
# requires: network
import bolzakit as bk

su = bk.spectral.eigenvalues("su2013", n=9)               # QUOTED values (Strohmaier–Uski 2013)
fem = bk.spectral.eigenvalues("fem")
print(max(abs(a - b) / b for (a, _), (b, _) in zip(fem, su)))   # about 3e-5
print(bk.spectral.selberg_check(0.5))                     # |geometric - spectral side| of the heat trace, ~1e-16
```

Outputs: `eigenvalues` is NUMERICAL, but the `su2013` values themselves are QUOTED (see
[docs/spectral.md](docs/spectral.md)). `selberg_check` is NUMERICAL.

### (c) Hyperbolic codes

For a homological code on a surface, the shortest non-trivial cycles set the code distance, and the automorphism group
supplies candidate symmetries. bolzakit gives the systole, the short geodesics with their multiplicities, and the
automorphism group. As a separate arithmetic fact, it decides when the golden geodesic is a systole of the congruence
cover Γ(p)\H of the modular surface. That result is about SL(2,ℤ), not about Bolza.

```python
import bolzakit as bk

print(bk.systole())                                       # 2*arccosh(1+sqrt2)
for e in bk.lengths(10):                                  # half-trace p + q*sqrt2 <= 10
    print(e["halftrace"], round(e["length"], 6), e["n_unoriented"])
print(bk.automorphisms().order)                           # 48
print([p for p in range(2, 100) if bk.arith.golden_geodesic_is_systole(p)])   # [3, 7, 47]
```

Outputs: `systole` is EXACT. `lengths` is NUMERICAL: the half-traces are exact and the multiplicities are rounded
from sums with integrality residual ≤ 1.1e-8. `automorphisms` is PROVED. `golden_geodesic_is_systole` is PROVED, with a
Lean proof.

### (d) Spectral geometry and quantum chaos

The Bolza surface is a classic test case for the Selberg trace formula and for quantum chaos. bolzakit gives the
primitive length spectrum up to half-trace 120, the trace-formula balance, and the spectral zeta function at 0. It also
gives certified bounds on a von Neumann trace in the group algebra.

```python
import bolzakit as bk

L = bk.lengths()                                          # all lengths with half-trace <= 120
print(len(L), L[-1]["length"])                          # 42 distinct lengths, longest 10.934
print(bk.spectral.zeta_zero())                            # -4/3 = chi/6 - 1, as a Fraction
tb = bk.spectral.tau_bounds("systole_shell")
print(tb["lower"], tb["upper"])                           # certified interval for tau(P_-(h))
```

```python
# requires: network
import bolzakit as bk

for t in (0.25, 0.5, 1.0):
    print(t, bk.spectral.selberg_check(t))                # residual: ~1e-16 at t = 0.25, 0.5; ~5e-11 at t = 1
```

Outputs: `lengths` is NUMERICAL, `zeta_zero` is PROVED, `tau_bounds` is CERTIFIED and `selberg_check` is NUMERICAL.

### (e) Vortex dynamics

A ring of N identical point vortices around a Weierstrass point is a natural test of how curvature and topology change
vortex stability. bolzakit gives the surface Green's function and Robin mass near the Weierstrass points. It gives the
Bergman-kernel constants that control small rings in closed form. It also gives the energy–momentum mode eigenvalues of
a ring in the plane, in H² and on the Bolza surface (Boatto–Koiller model [BK08]).

```python
import bolzakit as bk

c = bk.green.constants()                                  # closed forms at a Weierstrass point
print(c["B_Sigma_00"]["value"], c["beta3"]["value"])      # 0.69660196..., -0.16248381...
print(bk.vortex.lambda3_series()["coefficients"])         # lambda_3 = 12 xi - 4 xi^2 - 94.8177 xi^3 + O(xi^4)
for model in ("flat", "sinh", "tanh"):
    print(model, bk.vortex.ring(7, 0.01, model).modes()[3])   # flat: 0 (marginal); sinh 0.1224; tanh 0.1210
```

```python
# requires: slow
import bolzakit as bk

g = bk.green()                         # builds the heat-kernel image sum once (seconds to a minute)
print(g.robin_mass(0j))                # Robin mass at the Weierstrass point 0, about -0.02799
print(g.G(0.1 + 0j, 0.2j))             # both points within Poincaré radius 0.45 of one Weierstrass point
ring = bk.vortex.ring(7, 0.01, "bolza", D_ball=13.0)     # about 20 s to a few minutes
print(ring.modes()[3], ring.split(3))  # ~0.1195; split ~ 1000 xi^4 from the C8 symmetry of the point
```

Outputs: `green.constants` is PROVED. `lambda3_series` is PROVED, valid for N = 7 in the truncated Boatto–Koiller
model. `ring` is NUMERICAL. `green()` and its `G`, `regular` and `robin_mass` are NUMERICAL.

### (f) Arithmetic

The systole length is arithmetic. The number e^{ℓ_sys/2} = 1+√2+√(2+2√2) is a Salem number: it is the largest root of
x⁴−4x³−2x²−4x+1. The curve itself is y² = x⁵−x. With PARI you can compute its Euler factors and check that its Jacobian
is simple over ℚ.

```python
import math
import bolzakit as bk

print(bk.arith.bolza_quartic())                           # [1, -4, -2, -4, 1]
s = bk.arith.salem_number(dps=30)
print(s, abs(float(s) - math.exp(bk.systole() / 2)))      # largest root = e^{l_sys/2}
```

```python
# requires: pari
import bolzakit as bk

print(bk.arith.curve.euler_factor(3))                     # [1, 0, -2, 0, 9]: x^4 - 2x^2 + 9
print(bk.arith.curve.jacobian_is_Q_simple())              # True
```

Outputs: `bolza_quartic` is EXACT, `salem_number` is PROVED, `euler_factor` is EXACT and `jacobian_is_Q_simple` is
PROVED.

### (g) Teaching

The Bolza surface is a hyperbolic surface you can compute with. You can reduce any point of the disk into the octagon,
check the defining relation exactly, and measure distances.

```python
import bolzakit as bk

G = bk.group()
z, word, M = bk.geometry.reduce_point(0.9 + 0.1j)         # z = M(0.9+0.1i) lies in the octagon
print(z, word, bk.geometry.in_F(z))
print(G.word(G.RELATOR))                                  # ((1,0,0,0),(0,0,0,0)): the identity, exactly
print(bk.geometry.dist(0j, 0.5 + 0j))                     # 2 artanh(0.5)
print(bk.geometry.disk_to_upper(0j))                      # 1j
```

Outputs: the geometry functions are NUMERICAL (closed-form float64); `group().word` is exact.

## Provenance labels

Every labelled function accepts `with_provenance=True` and then returns a `Result(value, label, source, error)`.
`bolzakit.provenance.info(fn)` returns the label without calling the function. The full table of labelled functions is
in [docs/provenance.md](docs/provenance.md).

| label | meaning |
|---|---|
| EXACT | integer or algebraic arithmetic, with no rounding |
| PROVED | an exact computation plus a written argument, or a Lean theorem |
| CERTIFIED | a rigorous enclosure (interval or ball arithmetic, or exact tail bounds) |
| NUMERICAL | floating point with a stated error estimate |
| QUOTED | taken from the literature, with a reference |

```python
import bolzakit as bk
from bolzakit.provenance import info

r = bk.lengths(5, with_provenance=True)
print(r.label, r.error)                     # NUMERICAL  integrality residual <= 1.1e-8 ...
print(info(bk.green.constants).label)       # PROVED
print(info(bk.spectral.su2013).label)       # QUOTED
```

Some low-level helpers carry no label of their own: the arithmetic functions of the `group()` module, and the methods of
the `Automorphisms` and `Ring` objects.

## Results and corrections

Details and sources are in [docs/results-and-corrections.md](docs/results-and-corrections.md).

- **Aut⁺(Bolza) ≅ GL(2,3), not S₄×ℤ₂.** Strohmaier–Uski 2013 §7.1 say the Bolza surface "has symmetry group S4 ×Z2" of
  order 48 [SU13]. The two groups have the same order but are not isomorphic: GL(2,3) has 12 elements of order 8, and
  S₄×ℤ₂ has none. bolzakit computes the group exactly, as the quotient of the (2,3,8) triangle group by Γ, and gives an
  explicit isomorphism to GL(2,3).
- **Strict fractionality (math-pursuits item 04).** Let Γ be a torsion-free cocompact Fuchsian group, such as the
  Bolza group, and let h = h* ∈ C*_r(Γ) be nonzero with τ(h) = 0 (in particular, any such h in ℂΓ). Then
  0 < τ(P₋(h)) < 1, P₋(h) ∉ C*_r(Γ), and spec(h) = [a, b] with a < 0 < b. The first statement holds for every discrete
  group and uses only faithfulness of τ. The others use the Baum–Connes conjecture for hyperbolic groups [MY02, Laf12]
  through the Kadison–Kaplansky conjecture [MY02, Thm 21]. The ingredients are known, and the source note calls the
  corollary folklore-level.
- **Smaller results.** The length-spectrum multiplicities round to integers with residual ≤ 1.1e-8. The Weierstrass
  charts of radius 0.45 cover the surface (covering radius 0.4056). The m = 3 mode of a 7-ring at a Weierstrass point
  splits by about 1000·ξ⁴.

## Limitations

More detail is in [docs/limitations.md](docs/limitations.md).

- **It covers one surface.** Everything is for the Bolza surface. Codes, larger lattice experiments and genus-3 work
  need families of surfaces, which v0.1 does not have.
- **The FEM is a cross-check only.** The bundled finite-element eigenvalues are 9 clusters. They agree with
  Strohmaier–Uski to about 3e-5 and carry a conservative 2e-4 error label. They exist to check the quoted data. The real
  spectral work is carried by the Strohmaier–Uski data, which is downloaded on demand and not bundled.
- **There is no tight-binding layer.** There are no adjacency or graph-Laplacian Hamiltonians and no lattice patches,
  which most physics users build.
- **The Green's function is local to the Weierstrass charts.** `G(z, w)` is validated only when both points lie within
  Poincaré radius 0.45 of the same Weierstrass point; other pairs raise `OutOfValidityRange`. Every single point is
  covered, so `robin_mass(z)` works everywhere.
- **The λ₃ series** holds only for N = 7 in the SO(2)-invariant degree-7 truncation of the Boatto–Koiller model. In
  the full Bolza Hamiltonian the ring is only an approximate relative equilibrium.

## Roadmap (v0.2)

- **Congruence covers of Bolza.** Reduce the exact ℤ[ζ₈] generators modulo prime ideals of ℤ[ζ₈] and enumerate the
  finite quotients. For each one, build the cover: its group, its tiling as a cell complex, and its homological code
  parameters [[n, 2g, d]].
- **A tight-binding layer.** This means adjacency and graph-Laplacian Hamiltonians on the covers and on finite patches,
  abelian Bloch phases on the genus-2 cell, and a comparison of graph spectra with the continuum spectrum.

## Documentation

[quick-start](docs/quickstart.md) · [exact group](docs/exact.md) · [geometry](docs/geometry.md) ·
[spectral](docs/spectral.md) · [Green's function](docs/green.md) · [vortex](docs/vortex.md) ·
[arithmetic](docs/arith.md) · [results and corrections](docs/results-and-corrections.md) ·
[limitations](docs/limitations.md) · [provenance table](docs/provenance.md) ·
[citation checks](docs/CITATIONS_CHECK.md)

## References

Each reference was checked against its source; see [docs/CITATIONS_CHECK.md](docs/CITATIONS_CHECK.md).

- [BK08] S. Boatto, J. Koiller, *Vortices on closed surfaces*, arXiv:0802.4313 (2008).
- [BT16] N. P. Breuckmann, B. M. Terhal, *Constructions and noise threshold of hyperbolic surface codes*, IEEE Trans.
  Inf. Theory 62(6), 3731–3744 (2016); arXiv:1506.04029.
- [KFH19] A. J. Kollár, M. Fitzpatrick, A. A. Houck, *Hyperbolic lattices in circuit quantum electrodynamics*, Nature 571,
  45–50 (2019); arXiv:1802.09549.
- [L+22] P. M. Lenggenhager et al., *Simulating hyperbolic space on a circuit board*, Nat. Commun. 13, 4373 (2022);
  arXiv:2109.01148.
- [Laf12] V. Lafforgue, *La conjecture de Baum–Connes à coefficients pour les groupes hyperboliques*, J. Noncommut.
  Geom. 6 (2012) 1–197; arXiv:1201.4653.
- [MR21] J. Maciejko, S. Rayan, *Hyperbolic band theory*, Sci. Adv. 7, eabe9170 (2021); arXiv:2008.05489.
- [MY02] I. Mineyev, G. Yu, *The Baum–Connes conjecture for hyperbolic groups*, Invent. Math. 149, 97–122 (2002);
  arXiv:math/0105086.
- [SU13] A. Strohmaier, V. Uski, *An algorithm for the computation of eigenvalues, spectral zeta functions and
  zeta-determinants on hyperbolic surfaces*, Commun. Math. Phys. 317, 827–869 (2013); arXiv:1110.2150.

## Citing

If you use bolzakit, please cite it with the metadata in `CITATION.cff`. If you use the quoted eigenvalues, please
also cite [SU13].

## Licence

MIT. See [LICENSE](LICENSE).

## AI-assisted development

bolzakit was written with substantial help from AI coding agents, working from Gordon Speagle's mathematical notes (the
math-pursuits project) under his direction. To keep that honest:
- every returned value carries a provenance label;
- every literature reference in these docs was checked against its source;
- the Lean proofs build with no `sorry`, `admit` or `native_decide`;
- the tests check the numbers against independent computations.

The verification policy is in [VERIFICATION.md](VERIFICATION.md).
