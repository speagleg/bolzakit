# Arithmetic companions

Module: `bolzakit.arith`. Sources: the original research code (private archive).

## The Bolza quartic and the Salem number

- `bolza_quartic()` (EXACT) returns `[1, -4, -2, -4, 1]`, i.e. x⁴ − 4x³ − 2x² − 4x + 1. It is the norm from ℚ(√2) of
  x² − (2+2√2)x + 1, the characteristic polynomial of a systole translation g_k (trace 2+2√2, eigenvalues e^{±ℓ_sys/2}). The identity
  (x² − (2+2√2)x + 1)(x² − (2−2√2)x + 1) = x⁴ − 4x³ − 2x² − 4x + 1 is proved in Lean (`lean/BolzaKit/Salem.lean`,
  `item15_bolza_norm`).
- `salem_number(dps=50)` (PROVED) returns 1 + √2 + √(2+2√2) = e^{ℓ_sys/2} = 4.61158178930871…. The quartic is
  reciprocal and irreducible, with this root, its inverse, and two roots on |x| = 1. So it is a Salem number
  (proof below; PARI agrees).

### Proof (Salem number)

The quartic P(x) = x⁴ − 4x³ − 2x² − 4x + 1 is reciprocal, so with y = x + 1/x, P(x) = x²(y² − 4y − 4), whose roots
are y = 2 ± 2√2.
- y = 2 + 2√2 > 2: x² − (2+2√2)x + 1 has two real roots x₊ > 1 and 1/x₊, and
  x₊ = (1+√2) + √((1+√2)² − 1) = 1 + √2 + √(2+2√2) = 4.6115…. Since the systole translations g_k have trace 2+2√2
  (`exact/group.py`), their eigenvalues are x₊^{±1} = e^{±ℓ_sys/2}.
- y = 2 − 2√2 ∈ (−2, 2): the roots are a conjugate pair of non-real numbers with |x| = 1.

Irreducibility over ℚ: P has no rational root (the only candidates ±1 give P(1) = −8, P(−1) = 8). A factorisation
into two rational quadratics would put x₊ in a rational quadratic factor with another root of P. That root cannot be a
non-real root (the factor would have non-real coefficients) and cannot be 1/x₊ (then the sum of the pair, 2+2√2, would
be rational). So P is irreducible, x₊ is an algebraic integer > 1 whose conjugates are 1/x₊ and two numbers of modulus
1, i.e. a Salem number. □

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
    (`lean/BolzaKit/GoldenGeodesic.lean`). The geometric reading (Step 4) is proved below.

### Proof (golden geodesic)

Notation: F_n, L_n are the Fibonacci and Lucas numbers (L₀ = 2, L₁ = 1), Q = [[1,1],[1,0]] (det Q = −1),
Qⁿ = [[F_{n+1}, F_n],[F_n, F_{n−1}]], tr Qⁿ = L_n, det Qⁿ = (−1)ⁿ. For an odd prime p, Γ(p) = ker(SL₂(ℤ) → SL₂(𝔽_p)),
±Γ(p) = {A ∈ SL₂(ℤ) : A ≡ ±I mod p}, and X(p) = Γ̄(p)\H where Γ̄(p) is the common image of Γ(p) and ±Γ(p) in PSL₂(ℤ).
α(p) is the least k ≥ 1 with p | F_k, and n(p) is the least even n ≥ 2 with Qⁿ ∈ ±Γ(p). Only even powers of Q lie in
SL₂(ℤ); Q² is hyperbolic with axis from −1/φ to φ (the golden geodesic).

**Step 0 (n(p) = lcm(2, α(p))).** Qᵏ is scalar mod p iff p | F_k (off-diagonal entries); then F_{k+1} ≡ F_{k−1} =: c and
c² ≡ det Qᵏ = (−1)ᵏ. The k with Qᵏ scalar form the subgroup α(p)ℤ (Q is invertible over ℤ), and for even k, c ≡ ±1.
So the even n with Qⁿ ≡ ±I are the even multiples of α(p).

**Step 1 (systole trace).** If A ∈ SL₂(ℤ), A = ε(I + pX), ε = ±1, X ∈ M₂(ℤ), then 1 = det A = 1 + p tr X + p² det X, so
tr X = −p d with d = det X, and tr A = ε(2 − p²d). If A is hyperbolic then d ≠ 0 and |tr A| ≥ p² − 2 (d = 1),
with p² + 2 for d = −1 and ≥ 2p² − 2 for |d| ≥ 2. Equality is attained by [[1−p², p],[−p, 1]]. Also Γ̄(p) is
torsion-free (elliptic elements of SL₂(ℤ) have |tr| ≤ 1, impossible in ±Γ(p) apart from ±I). Since ℓ = 2 arccosh(|tr|/2),
the systole of X(p) is 2 arccosh((p²−2)/2). (Known: Murillo, Commun. Math. 34 (2026), Lemma 5.5.)

**Step 2 (Lucas lemmas).** (a) L_{2k} = L_k² − 2(−1)ᵏ: apply tr M² = (tr M)² − 2 det M (Cayley–Hamilton) to M = Qᵏ.
(b) If m is odd then L_r | L_{rm}: with M = Q^r, t = L_r, δ = det M, u_n = tr Mⁿ = L_{rn} satisfies
u_{n+2} = t u_{n+1} − δ u_n, u₁ = t, so t | u_n for all odd n by induction. (c) If k is even and L_k is prime, then
k = 2ʲ with j ≥ 1: write k = 2ʲm, m odd; by (b) L_{2ʲ} | L_k with 3 ≤ L_{2ʲ}, and L is strictly increasing on n ≥ 1, so
m > 1 would make L_k composite.

**Step 3 (theorem).** For an odd prime p: L_{n(p)} = p² − 2 iff p = L_{2ʲ} for some j ≥ 1, and then
n(p) = α(p) = 2^{j+1}. Qⁿ⁽ᵖ⁾ is hyperbolic (L_n ≥ 3 for even n ≥ 2) and lies in ±Γ(p), so L_{n(p)} ≥ p² − 2 by Step 1.
(⇒) Let n(p) = 2k and L_{2k} = p² − 2. If k is odd, (a) gives L_k² = p² − 4, so 0 < L_k < p and
L_k² ≤ (p−1)² < p² − 4 (as p ≥ 3), a contradiction. So k is even, L_k² = p², L_k = p is prime and k = 2ʲ by (c).
(⇐) Let p = L_k prime, k = 2ʲ, j ≥ 1, M = Qᵏ (det M = 1, tr M = p). Then Q²ᵏ = M² = pM − I ≡ −I (mod p), so
n(p) ≤ 2k. By Step 1 and (a), L_{n(p)} ≥ p² − 2 = L_k² − 2 = L_{2k}, so n(p) ≥ 2k by monotonicity. Hence
n(p) = 2k = 2^{j+1} and L_{n(p)} = p² − 2. Finally n(p) = lcm(2, α(p)): if α(p) were odd it would equal 2ʲ,
which is even, so α(p) = 2^{j+1}. The equivalence "L_{n(p)} = p² − 2 ⇔ p = L_{2ʲ}" is the Lean theorem
`golden_systolic_iff`.

**Step 4 (the geodesic, not only the lift).**
*Stabiliser of φ.* If A = [[a,b],[c,d]] ∈ SL₂(ℤ) fixes φ, then cφ² + (d−a)φ − b = 0; with φ² = φ + 1 and φ irrational,
b = c and a = c + d, so A = dI + cQ. The map ℤ[φ] → ℤ[Q], φ ↦ Q, is an injective ring map (Q² = Q + I) sending
u = d + cφ to A, and N(u) = det A = 1, so u is a unit. Replace A by −A to get u > 0, pick m with φᵐ ≤ u < φ^{m+1}
and put v = uφ^{−m} = x + yφ, so 1 ≤ v < φ and |v̄| = 1/v ≤ 1. Then y√5 = v − v̄ ∈ [0, φ+1), so y ∈ {0, 1}; y = 1 would
need an integer x = v − φ ∈ [1−φ, 0), which does not exist; hence v = 1, u = φᵐ, and N(u) = (−1)ᵐ = 1 forces m even.
So A = ±Q^{m} with m even; this A also fixes −1/φ.
*Stabiliser of the axis in Γ̄(p).* An element swapping the endpoints would fix a point of the axis, hence be elliptic,
which is impossible in the torsion-free Γ̄(p). So it fixes φ and equals ±Q^{2m} ∈ ±Γ(p); by Step 0 (applied to |2m|,
since ±Γ(p) is a group) n(p) | 2m. Hence the stabiliser is generated by Qⁿ⁽ᵖ⁾ and the golden geodesic projects to a
closed geodesic of length 2 arccosh(L_{n(p)}/2). Combined with Step 3, it is a systole of X(p) iff p = L_{2ʲ}. □

`golden_geodesic_is_systole` evaluates the criterion L_{n(p)} = p² − 2 exactly (integers); the Lean theorem checks
Steps 2–3 and the algebraic part of Step 4. The geometric facts used (length dictionary 2 arccosh(|tr|/2), closed
geodesics ↔ hyperbolic conjugacy classes, elliptic elements of SL₂(ℤ) have |tr| ≤ 1) are standard.

## The curve y² = x⁵ − x (needs `bolzakit[pari]`)

Without cypari2 these functions raise `OptionalDependencyMissing`.

- `curve.euler_factor(p)` (EXACT) returns the coefficients of PARI's `hyperellcharpoly` for a good prime p, for example
  `[1, 0, -2, 0, 9]` for p = 3, i.e. x⁴ − 2x² + 9. p = 2 has bad reduction and raises `ValueError`, as do non-primes.
- `curve.jacobian_is_Q_simple(primes=(3,))` (PROVED) returns True. An abelian surface whose Frobenius characteristic
  polynomial at a good prime is irreducible over ℚ is ℚ-simple, and x⁴ − 2x² + 9 is irreducible (proof below).
  - The criterion is one-sided. It returns True if L_p is irreducible at some given good prime, and None
    (inconclusive) otherwise; it never returns False. True is proved; None means inconclusive at the given primes.
    For example `primes=(5,)` returns None (L_5 is reducible), while `(3,)` and `(5, 3)` return True.

### Proof (Q-simplicity)

Let J be the Jacobian of C: y² = x⁵ − x (genus 2, good reduction at odd p, since disc(x⁵ − x) = −2⁸). For a good prime p,
the characteristic polynomial L_p of Frobenius on the ℓ-adic Tate module of J (ℓ ≠ p) has integer coefficients and is
an isogeny invariant that is multiplicative on products (Tate). If J were ℚ-isogenous to a product A × B of positive
dimension, L_p would factor over ℤ as L_p(A)·L_p(B) with factors of degree 2dim A and 2dim B. So irreducibility of L_p
over ℚ at one good prime proves that J is ℚ-simple.

For p = 3, PARI's `hyperellcharpoly` gives L₃ = x⁴ − 2x² + 9 (the data computed by `euler_factor`; this is the only
computed input). It is irreducible over ℚ: as a quadratic in x² its discriminant 4 − 36 = −32 < 0, so it has no
rational root x², hence no linear factor; and a factorisation (x² + ax + b)(x² − ax + b) needs b² = 9 and
2b − a² = −2, i.e. a² = 8 (b = 3) or a² = −4 (b = −3), both impossible over ℚ. □
