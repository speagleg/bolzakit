# Arithmetic companions

Module: `bolzakit.arith`. Sources: math-pursuits items 05 and 15.

## The Bolza quartic and the Salem number

- `bolza_quartic()` (EXACT) returns `[1, -4, -2, -4, 1]`, i.e. x⁴ − 4x³ − 2x² − 4x + 1. It is the norm from ℚ(√2) of
  x² − (2+2√2)x + 1, the characteristic polynomial of a systole translation g_k (trace 2+2√2, eigenvalues e^{±ℓ_sys/2}). The identity
  (x² − (2+2√2)x + 1)(x² − (2−2√2)x + 1) = x⁴ − 4x³ − 2x² − 4x + 1 is proved in Lean (`lean/BolzaKit/Salem.lean`,
  `item15_bolza_norm`).
- `salem_number(dps=50)` (PROVED) returns 1 + √2 + √(2+2√2) = e^{ℓ_sys/2} = 4.61158178930871…. The quartic is
  reciprocal and irreducible, with this root, its inverse, and two roots on |x| = 1. So it is a Salem number (item 15b,
  checked in PARI).

## Golden geodesic on Γ(p)\H

This concerns congruence covers of the **modular surface** (subgroups Γ(p) of SL(2,ℤ)), not the Bolza group.

- `lucas(n)` (EXACT) is the Lucas number Lₙ (L₀ = 2, L₁ = 1).
- `golden_geodesic_is_systole(p)` (PROVED) decides whether the golden geodesic (the axis of [[1,1],[1,0]]) is a systole
  of Γ(p)\H.
  - Let k be the rank of apparition of p (the least k ≥ 1 with p | F_k), and set n(p) = k if k is even, else 2k. Then
    the geodesic is a systole iff L_{n(p)} = p² − 2, iff p is a prime Lucas number L_{2^j}: 3, 7, 47, 2207, ….
  - The theorem covers odd primes only. Non-primes, p = 2 and non-int input (including bool) raise `ValueError`: they
    are outside the theorem's scope.
  - The equivalence "L_{n(p)} = p² − 2 ⇔ p = L_{2^j}" is the Lean theorem `golden_systolic_iff`
    (`lean/BolzaKit/GoldenGeodesic.lean`). The geometric reading is item 05's proof.

## The curve y² = x⁵ − x (needs `bolzakit[pari]`)

Without cypari2 these functions raise `OptionalDependencyMissing`.

- `curve.euler_factor(p)` (EXACT) returns the coefficients of PARI's `hyperellcharpoly` for a good prime p, for example
  `[1, 0, -2, 0, 9]` for p = 3, i.e. x⁴ − 2x² + 9. p = 2 has bad reduction and raises `ValueError`, as do non-primes.
- `curve.jacobian_is_Q_simple(primes=(3,))` (PROVED) returns True. An abelian surface whose Frobenius characteristic
  polynomial at a good prime is irreducible over ℚ is ℚ-simple, and x⁴ − 2x² + 9 is irreducible (item 15c).
  - The criterion is one-sided. It returns True if L_p is irreducible at some given good prime, and None
    (inconclusive) otherwise; it never returns False. True is proved; None means inconclusive at the given primes.
    For example `primes=(5,)` returns None (L_5 is reducible), while `(3,)` and `(5, 3)` return True.
