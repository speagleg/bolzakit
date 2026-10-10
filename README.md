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
protection more favourably than the surface code [BT16]. The lattices in these experiments are not the Bolza surface
(the circuit-QED lattices of [KFH19], for instance, are {7,3}-type), but for this kind of physics the Bolza surface is a
natural minimal, maximally symmetric test case.

bolzakit gives you the exact Fuchsian group (matrix entries in ℤ[ζ₈]), the automorphism group, the length spectrum,
Laplace eigenvalues, Green's-function and Bergman-kernel data at the Weierstrass points, vortex-ring stability data and a
few arithmetic facts. It ships with tests and three Lean 4 proofs (the group relation, the arithmetic core of the
golden-geodesic criterion and the Salem-norm identity). Every public function listed in
[docs/provenance.md](docs/provenance.md) carries a provenance label (EXACT, PROVED, CERTIFIED, NUMERICAL or QUOTED) —
ask for it with `with_provenance=True`.

Version 0.2 adds three things:
- a **certified layer** (`bolzakit.certified`): Arb-ball Rayleigh–Ritz upper bounds, trace-formula eigenvalue counts,
  a Lehmann–Goerisch lower bound, certified cluster perturbation tests and a certified trace chart, for Bolza and
  for small deformations of it. Every result comes with a JSON certificate, and the independent checker
  `bolzakit_check` re-derives the claims ([docs/certified.md](docs/certified.md));
- a standalone **sensitivity** module (`bolzakit.sensitivity`) for repeated eigenvalues of any finite-element pencil,
  with scikit-fem and FEniCSx adapters ([docs/sensitivity.md](docs/sensitivity.md));
- the proof package [proofs/bolza-local-max/](proofs/bolza-local-max/): a computer-assisted proof that the Bolza
  surface is a strict local maximum of λ₁, with its golden certificates and `reproduce.py`, and a Lean 4 skeleton
  of its logical chain ([below](#the-bolza-surface-is-a-strict-local-maximum-of-λ₁-computer-assisted-conditional-on-cited-inputs)).

## Contents

- [Install](#install)
- [Quick-start](#quick-start)
- [The Bolza surface is a strict local maximum of λ₁ (computer-assisted, conditional on cited inputs)](#the-bolza-surface-is-a-strict-local-maximum-of-λ₁-computer-assisted-conditional-on-cited-inputs)
- [For software engineers](#for-software-engineers)
- [For physicists](#for-physicists)
- [For engineers](#for-engineers)
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
- [Roadmap](#roadmap)
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
pip install "bolzakit[cert]"         # certified layer: python-flint and sympy
pip install "bolzakit[skfem]"        # scikit-fem adapter for bolzakit.sensitivity
```

You need Python 3.10 or later. The Strohmaier–Uski eigenvalue file is not bundled. It is downloaded from arXiv on first
use, checked against a fixed sha256 and cached in `~/.cache/bolzakit` (or `$BOLZAKIT_CACHE`).

## Quick-start

```python
import bolzakit as bk

G = bk.group()
print(len(G.GENS))                          # 8 side-pairing generators
assert bk.automorphisms().order == 48       # = |GL(2,3)|
print(bk.lengths(max_half_trace=5)[0])      # the systole class: half-trace 1+√2, 24 oriented geodesics
print(bk.systole())                         # 3.057141838961996
print(bk.systole(with_provenance=True).label)   # NUMERICAL (float64 of the exact closed form 2·arccosh(1+√2))
```

## The Bolza surface is a strict local maximum of λ₁ (computer-assisted, conditional on cited inputs)

[proofs/bolza-local-max/](proofs/bolza-local-max/) holds a computer-assisted proof, conditional on cited inputs, that
λ₁(X) < λ₁(B) for every surface X ≠ B in a neighbourhood of the Bolza surface B in the Teichmüller space of genus 2. The package contains
the integrated proof (PROOF.md), the paper, the review record, ten golden certificates and `reproduce.py`.

**Status.** The proof is computer-assisted and conditional on the cited inputs listed in PROOF.md §Dependencies,
including Fortier Bourque–Petri's Lemmas 7.2–7.3 and geodesic data, the Selberg trace formula, the
Lehmann–Goerisch theorem (quoted from secondary sources), and the trusted software (python-flint/Arb, sympy in the
producer). The Strohmaier–Uski
eigenvalue (SU13) is optional: λ₁ ≥ 3.8388757010 is certified in-house, and each SU13-based number has an
SU13-free twin. No human expert has read the proof yet. Its "second reads" were done by AI reviewers. The checker
re-derives the certified claims from stored Gram data, but it does not re-check the element assembly or the links
between certificates, which are made by hand ([trust model](docs/certified.md#trust-model)). The Lean 4 skeleton
under `lean/BolzaKit/LocalMax` checks the logical chain, conditional on the computational inputs.

```
python -m bolzakit_check proofs/bolza-local-max/certificates/*.json   # exit 2: all VALID, six conditionally
python proofs/bolza-local-max/reproduce.py --quick                     # about 1-2 min
```

## For software engineers

**Use case: checking a certificate you did not compute.** Each certified computation writes a JSON certificate
(schema `bolzakit-cert/1`, with a sha256 of its canonical core). `bolzakit_check` imports nothing from `bolzakit`. It
re-derives every claim from the stored inputs and data in its own Arb code, so a tampered claim fails even if the
sha256 is recomputed. On the command line it exits 0 (all VALID), 2 (some only conditionally VALID) or 1 (some
INVALID). CI runs it on the golden set and accepts 0 or 2. Run this from the repository root:

```python
# requires: cert
import glob
from bolzakit_check.check import check_file
for f in sorted(glob.glob("proofs/bolza-local-max/certificates/*.json")):
    ok, msg = check_file(f)                   # e.g. "VALID given stored Gram data"
    assert ok, msg
```

"VALID given stored Gram data" means that the checker trusts the element integrals. `--deep` re-integrates a sample
of elements, as a consistency check only. See [docs/certified.md](docs/certified.md#trust-model).

## For physicists

**Use case: rigorous eigenvalue bounds for a hyperbolic surface.** `certify_spectrum` gives certified upper bounds
on λ₁..λ_k for Bolza or a deformation of it (|t| ≤ 0.05), and `count_eigenvalues` gives certified Selberg
trace-formula counts for Bolza. The snippet below certifies that a small F_A deformation lowers λ₁ below the
certified lower bound λ₁(Bolza) ≥ 3.8388757010 (about 15 s):

```python
# requires: cert
from bolzakit.certified import Surface, certify_spectrum, count_eigenvalues
assert count_eigenvalues("bolza", window=(0, 3.3)).bound < 1         # certified: no eigenvalue in (0, 3.3]
r = certify_spectrum(Surface.bolza().deform("F_A", 0.02), k=3, degree=2, n=12)
print(r.upper)                    # certified upper bounds on lambda_1..3: [3.74456, 3.84164, 3.93680]
assert r.certified and r.upper[0] < 3.8388757010
```

[examples/bolza_deform.py](examples/bolza_deform.py) does this for three directions, and
[examples/bolza_strain_rates.py](examples/bolza_strain_rates.py) prints the splitting rates of the triple
eigenvalue and a certified Ky Fan bound. The validity limits (deformations near Bolza only, counting for Bolza
only, weak certificates for arbitrary vector directions) are in
[docs/certified.md](docs/certified.md#known-limitations).

## For engineers

**Use case: design sensitivity of a repeated natural frequency.** Symmetric structures (square plates and
membranes, rings, rotors) have repeated eigenvalues. Their individual branches are not differentiable, but the sum
over the cluster is. `cluster_sensitivity` takes the stiffness K, the mass M and their design derivatives from any
FE code, and returns the splitting rates, the Hessian of the cluster sum and a Ky Fan bound. It needs only numpy
and scipy:

```python
import numpy as np
from bolzakit.sensitivity import cluster_sensitivity
K, M = np.diag([2.0, 2.0, 5.0]), np.eye(3)            # a repeated eigenvalue 2
dK = np.array([[1, 1, 0], [1, -1, 1], [0, 1, 0.0]])    # derivative of K along a design change
r = cluster_sensitivity(K, M, dK, 0 * M, (0, 1), d2K=0 * M, d2M=0 * M)
print(r.rates, r.trace_hessian)                       # [-1.414 1.414] (+-sqrt 2) -0.667 (= 2/(2-5))
```

The results are NUMERICAL and belong to the discrete pencil, not to the continuum. For a nearly (not exactly)
repeated eigenvalue, the rates are eig(A) and not the branch derivatives, but the cluster-sum derivatives remain
correct. [examples/membrane_repeated_frequency.py](examples/membrane_repeated_frequency.py) runs the scikit-fem
adapter on a stretched membrane and checks the results against the exact values. There is also a FEniCSx adapter
with exact UFL derivatives ([docs/sensitivity.md](docs/sensitivity.md)).

**Use case: trimming a split mode pair.** A ring gyroscope or any other degenerate-mode resonator needs its mode
pair at one frequency, and fabrication splits it. `trim_degenerate_pair` takes the FE model and the per-unit effect of
candidate trims (point masses, ablation, stiffness patches) and returns the cheapest trim (least Σ c|t|, so ablation
costs like addition) that ties the pair: a linear program on the splitting matrices, then Newton re-solves.
[examples/gyro_trim.py](examples/gyro_trim.py) checks it with an independent dense eigensolve. On a ring on 8 spokes
(no analytic rule) it tied every one of 5 random seeds in 4 eigensolves (status `converged`), from relative splits of
3.6e-4 to 8.9e-4 to below the eigensolver's resolution estimate, ≤ rel_floor (2.0e-10); the independent dense LAPACK
check gives ≤ 1.1e-10, its own rounding floor, so it confirms the tie to that level and no further. Iterated antinode
trimming tied 4 of the 5 seeds, using 1.2–5.6× our mass and 13–29 eigensolves,
and its 0.5-damped variant tied none within 40 rounds. On a plain ring both reach the resolution floor (≤ 1.6e-8);
there the damped baseline matched our mass (1.00×) but needed 19 eigensolves to our 7
([docs/sensitivity-trim.md](docs/sensitivity-trim.md)):

```python
import numpy as np
from bolzakit.sensitivity.trim import trim_degenerate_pair
K, E = np.diag([1.0, 2.0, 2.01, 5.0]), np.eye(4)      # lambda_1, lambda_2 split by 0.5 %
r = trim_degenerate_pair(K, np.eye(4), [np.outer(E[1], E[1]), np.outer(E[2], E[2])], None, (1, 2))
print(r.status, r.t)                                  # converged [0.01 0.  ]: stiffen dof 1 by 0.01
assert r.status == "converged" and np.allclose(r.t, [0.01, 0.0])
```

**Use case: maximising a repeated frequency.** `bolzakit.sensitivity.optimize.maximize_eigenvalue` raises the lowest
eigenvalue (or a gap) of a parametrised pencil with steepest ascent that is correct at multiple eigenvalues.
[benchmarks/frequency_maximisation.py](benchmarks/frequency_maximisation.py) compares it with the standard bound
formulation (SLSQP, trust-constr, MMA) on a plate whose optimum has a repeated fundamental frequency. On that
benchmark SLSQP beat it on success rate and cost: SLSQP came within 0.1 % of the best value found within the
200-eigensolve budget (not the optimum) from 3 of 5 starts in a median of 40 eigensolves, ours from 1 of 5 in 169.
Ours wins only on monotone progress (its accepted iterates never decreased; SLSQP's callback iterates decreased 138
times, and the benchmark notes define the count per method) and on its certified bound `phi_upper` on the ascent
rate. It is a first-order method and can end at `max_iter` short of the maximum: in the tests on an ill-conditioned
instance, and on four of the five benchmark starts (there the eigensolve budget; the cause is not established).
The table is in [docs/sensitivity-optimize.md](docs/sensitivity-optimize.md).

## Use cases

Every Python snippet in this README, including those in the audience sections above, runs as written:
`tests/test_readme_snippets.py` executes each one (from the repository root). Blocks marked
`# requires: network` need the Strohmaier–Uski file (downloaded once); the test suite runs them only when that file is
already cached or with `pytest --runnetwork`. Blocks marked `# requires: slow` take from seconds to a few minutes
(`pytest --runslow`). Blocks marked `# requires: pari` need `bolzakit[pari]`, and blocks marked `# requires: cert`
need `bolzakit[cert]`.

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
assert A.element_order_counts()[8] == 12    # 12 elements of order 8: GL(2,3), not S4 x Z2
```

Outputs: the `group()` arithmetic is exact in ℤ[ζ₈] (its helper functions carry no label of their own; `elements` is EXACT;
`systole` is NUMERICAL, the float64 evaluation of the exact closed form 2·arccosh(1+√2)). `automorphisms()` is PROVED.

### (b) Circuit-QED and electric-circuit lattices

Under refinement of the tiling, low-lying graph spectra are expected to approximate the continuum spectrum; bolzakit
provides the continuum side for such comparisons. It ships 9 finite-element eigenvalue clusters as a quick offline reference. It also fetches
the 1000 Strohmaier–Uski eigenvalues on demand.

```python
import bolzakit as bk

fem = bk.spectral.eigenvalues("fem")
assert len(fem) == 9                                      # 9 clusters, bundled
for lam, mult in fem:
    print(f"{lam:9.4f}  x{mult}")
print(bk.provenance.info(bk.spectral.eigenvalues).error)  # FEM: relative error <= 2e-4 vs Strohmaier–Uski
```

```python
# requires: network
import bolzakit as bk

su = bk.spectral.eigenvalues("su2013", n=9)               # QUOTED values (Strohmaier–Uski 2013)
fem = bk.spectral.eigenvalues("fem")
rel = max(abs(a - b) / b for (a, _), (b, _) in zip(fem, su))
assert rel < 3e-5                                         # measured 2.98e-5
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
odd_primes = [p for p in range(3, 100) if all(p % d for d in range(2, int(p ** 0.5) + 1))]
hits = [p for p in odd_primes if bk.arith.golden_geodesic_is_systole(p)]
assert hits == [3, 7, 47]                                 # prime Lucas numbers L_{2^j}
```

Outputs: `systole` is NUMERICAL (float64 of the exact closed form). `lengths` is NUMERICAL: the half-traces are exact
and the multiplicities are rounded from sums with integrality residual ≤ 1.1e-8. `automorphisms` is PROVED.
`golden_geodesic_is_systole` is PROVED: the arithmetic step (L_{n(p)} = p²−2 ⇔ p = L_{2^j}) is a Lean theorem; the
geometric part is item 05's written proof. It accepts odd primes only and raises `ValueError` otherwise (other p are
outside the theorem's scope).

### (d) Spectral geometry and quantum chaos

The Bolza surface is a classic test case for the Selberg trace formula and for quantum chaos. bolzakit gives the
primitive length spectrum up to half-trace 120, the trace-formula balance, and the spectral zeta function at 0. It also
gives certified bounds on a von Neumann trace in the group algebra.

```python
from fractions import Fraction
import bolzakit as bk

L = bk.lengths()                                          # all lengths with half-trace <= 120
print(len(L), L[-1]["length"])                          # 42 half-traces, longest 10.934
assert sum(e["n_oriented"] == 0 for e in L) == 2          # 5+4√2 and 33+24√2: no primitive geodesic
assert bk.spectral.zeta_zero() == Fraction(-4, 3)         # chi/6 - 1
tb = bk.spectral.tau_bounds("systole_shell")
print(tb["lower"], tb["upper"])                           # certified interval for tau(P_-(h))
```

```python
# requires: network
import bolzakit as bk

for t in (0.25, 0.5, 1.0):
    print(t, bk.spectral.selberg_check(t))                # residual: ~1e-16 at t = 0.25, 0.5; ~5e-11 at t = 1
```

Outputs: `lengths` is NUMERICAL, `zeta_zero` is PROVED, `tau_bounds` is CERTIFIED and `selberg_check` is NUMERICAL;
`selberg_check` raises `OutOfValidityRange` outside its tested window, about 0.03 ≤ t ≤ 1.

### (e) Vortex dynamics

A ring of N identical point vortices around a Weierstrass point is a natural test of how curvature and topology change
vortex stability. bolzakit gives the surface Green's function and Robin mass near the Weierstrass points. It gives the
Bergman-kernel constants that control small rings in closed form. It also gives the energy–momentum mode eigenvalues of
a ring in the plane, in H² and on the Bolza surface (Boatto–Koiller model [BK08]).

```python
import bolzakit as bk

c = bk.green.constants()                                  # closed forms at a Weierstrass point
print(c["B_Sigma_00"]["value"], c["beta3"]["value"])      # 0.69660196..., -0.16248381...
coeffs = bk.vortex.lambda3_series()["coefficients"]       # lambda_3 = 12 xi - 4 xi^2 - 94.8177 xi^3 + O(xi^4)
assert coeffs[1] == 12 and coeffs[2] == -4 and abs(coeffs[3] + 94.8177) < 1e-4
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
assert 900 < ring.split(3) / 0.01**4 < 1100
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

assert bk.arith.bolza_quartic() == [1, -4, -2, -4, 1]    # x^4 - 4x^3 - 2x^2 - 4x + 1
s = bk.arith.salem_number(dps=30)
print(s, abs(float(s) - math.exp(bk.systole() / 2)))      # largest root = e^{l_sys/2}
```

```python
# requires: pari
import bolzakit as bk

assert bk.arith.curve.euler_factor(3) == [1, 0, -2, 0, 9]  # x^4 - 2x^2 + 9
print(bk.arith.curve.jacobian_is_Q_simple())              # True (proved by L_3 irreducible)
assert bk.arith.curve.jacobian_is_Q_simple(primes=(5,)) is None   # inconclusive: L_5 is reducible
```

Outputs: `bolza_quartic` is EXACT, `salem_number` is PROVED, `euler_factor` is EXACT and `jacobian_is_Q_simple` is
PROVED when it returns True; it returns None (inconclusive) when L_p is reducible at every given prime.

### (g) Teaching

The Bolza surface is a hyperbolic surface you can compute with. You can reduce any point of the disk into the octagon,
check the defining relation exactly, and measure distances.

```python
import bolzakit as bk

G = bk.group()
z, word, M = bk.geometry.reduce_point(0.9 + 0.1j)         # z = M(0.9+0.1i) lies in the octagon
assert bk.geometry.in_F(z)
print(z, word)
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
| PROVED | an exact computation plus a written argument, or a Lean theorem; arguments may be AI-written and cited (see [AI-assisted development](#ai-assisted-development)) |
| CERTIFIED | a rigorous enclosure (interval/ball arithmetic, or a rigorous a-priori floating-point error analysis) |
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
  S₄×ℤ₂ has none. bolzakit computes the group as the quotient of the (2,3,8) triangle group by Γ and gives an explicit
  isomorphism to GL(2,3): the normaliser check is exact; the 48 cosets and their multiplication table are enumerated
  with floating-point Möbius maps (keyed by rounding), and the index 48 agrees with the area formula 4π/(π/12).
- **Strict fractionality (math-pursuits item 04).** Let Γ be a torsion-free cocompact Fuchsian group, such as the
  Bolza group, and let h = h* ∈ C*_r(Γ) be nonzero with τ(h) = 0 (in particular, any such h in ℂΓ). Then
  0 < τ(P₋(h)) < 1, P₋(h) ∉ C*_r(Γ), and spec(h) = [a, b] with a < 0 < b. The first statement holds for every discrete
  group and uses only faithfulness of τ. The others use the Baum–Connes conjecture for hyperbolic groups [MY02, Laf12]
  through the Kadison–Kaplansky conjecture [MY02, Thm 21]. The ingredients are known, and the source note calls the
  corollary folklore-level.
- **Smaller results.** The length-spectrum multiplicities round to integers with residual ≤ 1.1e-8. The Weierstrass
  charts of radius 0.45 cover the surface (covering radius 0.4056). The m = 3 mode of a 7-ring at a Weierstrass point
  splits by about 1000·ξ⁴.
- **Ring mode-pair split coefficient 3/10, not 1/4.** Research note 02 derived split rates ±ω²ε/4 for a 4th-harmonic
  density imperfection of a thin ring; the correct rates are ±(3/10)λ̄ε (measured 0.29984 by `examples/gyro_trim.py`).

## Limitations

More detail is in [docs/limitations.md](docs/limitations.md).

- **It covers one surface.** Everything is for the Bolza surface. The certified layer adds only deformations near
  it (|t| ≤ 0.05), and its eigenvalue counts are for Bolza only. Codes, larger lattice experiments and genus-3 work
  need families of surfaces, which bolzakit does not have.
- **Certified runs are slow, and the checker has a trust boundary.** See
  [docs/certified.md](docs/certified.md#known-limitations).
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

## Roadmap

- **Vector directions in the certified layer.** Replace the minimum-norm second-order jet with a lower-curvature
  curve, so that arbitrary deformation directions certify as well as the named families.
- **Congruence covers of Bolza.** Reduce the exact ℤ[ζ₈] generators modulo prime ideals of ℤ[ζ₈] and enumerate the
  finite quotients. For each one, build the cover: its group, its tiling as a cell complex, and its homological code
  parameters [[n, 2g, d]].
- **A tight-binding layer.** This means adjacency and graph-Laplacian Hamiltonians on the covers and on finite patches,
  abelian Bloch phases on the genus-2 cell, and a comparison of graph spectra with the continuum spectrum.

## Documentation

[quick-start](docs/quickstart.md) · [exact group](docs/exact.md) · [geometry](docs/geometry.md) ·
[spectral](docs/spectral.md) · [certified layer and checker](docs/certified.md) ·
[repeated-eigenvalue sensitivity](docs/sensitivity.md) ([optimisation](docs/sensitivity-optimize.md),
[trimming](docs/sensitivity-trim.md)) · [Green's function](docs/green.md) · [vortex](docs/vortex.md) ·
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

bolzakit's code was ported from the math-pursuits research project. The derivations and notes there, and their
second-reader reviews, were produced by AI agents (Claude) working under Gordon Speagle's direction. No human specialist
or peer review has checked them yet. The provenance labels, including PROVED, reflect machine-checked computation, Lean
proofs where stated, and AI-written arguments that are cited. Treat PROVED claims that rest on such arguments
accordingly. bolzakit's code, tests, documentation and Lean files were also written and reviewed by AI agents (Claude)
under Gordon Speagle's direction; Gordon directed and approved the work.

What is in place:
- every labelled function (the public factories and the functions listed in [docs/provenance.md](docs/provenance.md))
  can return its provenance label; low-level helpers and object methods carry none (see
  [Provenance labels](#provenance-labels));
- every literature reference in these docs was checked against its source;
- the Lean proofs build with no `sorry`, `admit` or `native_decide`;
- the tests check the numbers against independent computations.

The verification policy is in [VERIFICATION.md](VERIFICATION.md).
