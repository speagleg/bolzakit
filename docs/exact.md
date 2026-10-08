# Exact group data

Modules: `bolzakit.exact.group` (returned by `bk.group()`), `bolzakit.exact.automorphisms`, `bolzakit.exact.lengths`.

## The Fuchsian group

The Bolza group Γ ⊂ SU(1,1) acts on the Poincaré disk. Its fundamental domain is the regular octagon F centred at 0,
with interior angles π/4. The side midpoints are at angles kπ/4 and the vertices at (2k+1)π/8.

The generator g_k (k = 0..7) is the hyperbolic translation of length ℓ_sys = 2·arccosh(1+√2) along the axis at angle
kπ/4. It pairs side k+4 with side k, and g_{k+4} = g_k⁻¹. With μ = √(2+2√2), every group element has the form

  γ = [[A, μB], [μ·conj B, conj A]],   A, B ∈ ℤ[ζ₈],   |A|² − μ²|B|² = 1,

and is stored exactly as the pair (A, B) of integer 4-tuples (coefficients of 1, ζ, ζ², ζ³). Since μ² ∈ ℤ[ζ₈],
products stay in this form. Half-traces Re A lie in ℤ[√2] and are returned as pairs (p, q), meaning p + q√2.

| name | what it is | label |
|---|---|---|
| `GENS` | the 8 generators, `GENS[k] = g_k` | exact data |
| `RELATOR` | `(0, 5, 2, 7, 4, 1, 6, 3)`, i.e. g0 g1⁻¹ g2 g3⁻¹ g0⁻¹ g1 g2⁻¹ g3 = 1 | exact data |
| `word(ws)` | the product of `GENS[k]` for k in `ws` | exact (unlabelled helper) |
| `det(g)`, `half_trace(g)` | exact values in ℤ[√2] as (p, q) | exact (unlabelled helpers) |
| `to_matrix(g)`, `act(g, z)` | complex float matrix; Möbius action on the disk | float (unlabelled helpers) |
| `elements(max_word_length)` | all distinct elements of word length ≤ L | EXACT |
| `ball_by_displacement(Rmax)` | all elements with d(0, γ0) ≤ Rmax, as integer rows | exact enumeration (unlabelled) |
| `L_SYS`, `systole()` | 3.0571418389619960 = 2·arccosh(1+√2) | NUMERICAL (float64 evaluation of the exact closed form) |

The relation is also proved in Lean (`lean/BolzaKit/Group.lean`, theorem `bolza_relation`) by kernel evaluation.

## Automorphisms

`bk.automorphisms()` (PROVED) returns an `Automorphisms` object. Aut⁺ = N/Γ, where N = ⟨Γ, ρ, σ⟩ is the (2,3,8)
triangle group: ρ is the rotation by π/4 about 0, and σ is the half-turn about the point at distance r_in/2 from 0
towards the midpoint of side 0. bolzakit checks exactly that the
generators of N conjugate Γ into Γ, then enumerates the quotient.

| member | value |
|---|---|
| `.order` | 48 |
| `.full_order` | 96 (with the reflection z ↦ z̄) |
| `.element_order_counts()` | `{1: 1, 2: 13, 3: 8, 4: 6, 6: 8, 8: 12}` |
| `.isomorphism_to_GL23()` | dict: element index → 2×2 matrix over F₃, `((a, b), (c, d))` |
| `.weierstrass_points()` | the 6 Weierstrass points (fixed points of the hyperelliptic involution), as points of the disk |
| `.stabiliser(z)` | indices of the elements fixing z on the surface |
| `.chart_to(p)` | a 2×2 matrix M of N with M(p) = 0 mod Γ (a chart centred at p) |

The element orders alone show Aut⁺ ≇ S₄×ℤ₂; see [results-and-corrections.md](results-and-corrections.md).

## Length spectrum

`bk.lengths(max_half_trace=120, compute=False)` (NUMERICAL) returns a list of dicts sorted by length:

- `halftrace`: (p, q), exact, meaning p + q√2 = cosh(ℓ/2);
- `length`: ℓ = 2·arccosh(p + q√2), as a float;
- `n_oriented`: the number of oriented primitive closed geodesics of that length (the count the Selberg trace formula
  uses);
- `n_unoriented`: `n_oriented / 2`.

The data up to half-trace 120 is bundled (42 distinct lengths). Asking for more raises `OutOfValidityRange` unless you
pass `compute=True`, which enumerates afresh; the cost grows like e^ℓ.

**Method.** For each exact half-trace t, bolzakit sums the length of axis(γ) ∩ F over all γ with |half-trace| = t. That
sum equals Σ ℓ(primitive root) over the oriented classes with that trace. Subtracting the contributions of powers and
dividing by ℓ(t) gives a number that must be an integer. The largest distance from an integer is 1.09e-8 (at
59 + 42√2), which is why the label is NUMERICAL rather than EXACT. The multiplicities are confirmed by the Selberg trace
formula ([spectral.md](spectral.md)).

The first entries are: ℓ = 3.0571 (24 oriented), 4.8969 (24) and 5.8281 (48).
