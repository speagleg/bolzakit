# Green's function, Robin mass and Bergman constants

Module: `bolzakit.green`. `bk.green` is both a module and a callable: `bk.green(...)` is `bk.green.green(...)`, and
`info(bk.green)` and `inspect.signature(bk.green)` describe that function. Static type checkers see `bk.green` as a
module, so `bk.green.green(T, D_ball)` is the type-checker-friendly spelling.

## API

| call | returns | label |
|---|---|---|
| `bk.green(T=2.6, D_ball=15.0)` | a `SurfaceGreen` (cached per `(T, D_ball)`) | NUMERICAL |
| `g.G(z, w)` | G_Σ(z, w), with −Δ G_Σ = δ − 1/A and ∫ G_Σ dA = 0 | NUMERICAL |
| `g.regular(z, w)` | R(z, w) = G_Σ(z, w) − G₁(d(z, w)) | NUMERICAL |
| `g.robin_mass(z)` | m(z) = R(z, z) + (ln 2)/2π = lim [G_Σ + (1/2π) ln d] | NUMERICAL |
| `g.chart(z)` | `(index, z_local)` of the Weierstrass chart nearest z | (unlabelled helper) |
| `bk.green.constants(dps=30)` | `B_Sigma_00`, `beta1`, `beta2`, `beta3`, each `{value, closed_form, label}` | PROVED |
| `bk.green.bergman_at_weierstrass(dps=50, M=24, npts=50)` | B_Σ(0,0) from a Hejhal-type 1-form computation | NUMERICAL |
| `bk.green.VALIDITY_RADIUS` | `0.45` | |

G₁(d) = −(1/2π) ln tanh(d/2) is the free Green's function on H². The area is A = 4π.

**Method.** The Green's function uses a heat-kernel splitting at time T. The short-time part is an image sum over the
exact group elements with d(0, γ0) ≤ `D_ball`, and the spectral tail Σ e^{−λₙT} φₙφₙ/λₙ is omitted. At T = 2.6 the
omitted tail is ≤ 4e-6 in value, and ΔG = 1/A holds to 2e-5 (the label's error field). The first call for a given
`(T, D_ball)` builds the image ball. On the development machine this took about 8 s at the defaults; budget up to a
minute on slower machines.

**Validity.** Both points must lie within Poincaré radius 0.45 of the same Weierstrass point, measured in a chart
centred there. `G`, `regular` and `robin_mass` pick the chart that minimises max(|z_local|, |w_local|).
- Only the validated parameter pairs `(T, D_ball)` = (2.6, 15.0) (default) and (2.6, 13.0) are accepted
  (`bk.green.VALIDATED_T_D`); any other pair raises `OutOfValidityRange` before anything is built. At (2.6, 13),
  G(0.1, 0.2i) and robin_mass(0) differ from the defaults by 1.1e-7 and 8.2e-8. Smaller T leaves a larger omitted
  tail (T = 2.0: −3.1e-5 in G, −3.9e-5 in the Robin mass, so it is not accepted), and larger T needs a larger image
  ball; un-validated pairs can be off by up to 25% (T = 10 at D_ball = 15).
- A pair of points with no common chart raises `OutOfValidityRange`.
- `G(z, w)` raises `ValueError` when z and w coincide on the surface (local Poincaré distance < 1e-9, e.g. w = z or
  w = γz for γ ∈ Γ); use `robin_mass(z)` or `regular(z, z)` instead. Distances are computed as
  d = 2 asinh(|z − w| / √((1 − |z|²)(1 − |w|²))), which is exact for nearby points (no cancellation in 1 + u).
- The six charts cover the surface: the covering radius is ρ = 0.4056164 < 0.45 (Poincaré distance 0.8607). The worst
  point is at argument π/8 mod π/4 in the octagon, equidistant from the centre chart and two edge-midpoint charts. So
  every *single* point is valid, and `robin_mass(z)` works everywhere.

```python
import bolzakit as bk
g = bk.green()                    # slow the first time
g.robin_mass(0j)                  # -0.027989475...
g.G(0.1 + 0j, 0.2j)               # 0.104134729...
```

## Closed form of the Bergman constants (the PROVED label)

`constants()` is labelled PROVED, and its source field cites this page. The derivation below is self-contained: it
restates, in condensed form, the argument of the original research code (private archive, not public), where it was
written by an AI agent and then re-derived by an independent AI second reader (Claude). It has not been checked by a
human specialist. The standard facts it quotes are named where they are used.

### 1. What the constants are

Lift the Weierstrass point to 0 in the disk (curvature −1, dA = 4 dx dy/(1−|z|²)², Δ = (1−|z|²)² ∂_z∂_z̄, A = 4π).
G_Σ is the Green's function with −ΔG_Σ(·, w) = δ_w − 1/A, ∫G_Σ dA = 0. Its free part is
G₁(d) = −(1/2π) ln tanh(d/2) = −(1/4π)[ln|z−w|² − ln|1−z w̄|²], and R := G_Σ − G₁ is smooth near the diagonal.

*Lemma (two identities).* (i) ∂_z∂_z̄ R(z, w) = 1/(A(1−|z|²)²): off the diagonal ΔG_Σ = 1/A and ΔG₁ = 0, and the
identity extends across the diagonal because R is smooth. (ii) For z ≠ w, ∂_z∂_w̄ G_Σ(z, w) = −¼ B_Σ(z, w): the function
K = ∂_z∂_w̄ G_Σ is holomorphic in z and antiholomorphic in w for z ≠ w (apply ∂_z̄ or ∂_w and use (i) and its symmetric
version); near the diagonal ∂_z∂_w̄[−(1/4π) ln|z−w|²] = ¼ δ is a delta, not a pole, so K(·, w) dz is a holomorphic
1-form; for holomorphic f, Stokes gives ∫(∂_zG) f̄ dz∧dz̄ = 0, and applying ∂_w̄ distributionally gives
∫K(z, w) conj f(z) dx dy = −¼ conj f(w), i.e. K̄ = −¼ B by the reproducing property. The same identity on the disk gives
∂_z∂_w̄ G₁ = −1/(4π(1−z w̄)²). (This identity is also in Gustafsson, arXiv:2405.19215, eqs. (5.23)–(5.24), and
Grotta-Ragazzo, J. Nonlinear Sci. 2024, Thm 2.2; we have not re-read those papers.) Hence
∂_z∂_w̄ R = −¼(B_Σ − B_𝔻). Here B_Σ(z, w) = Σᵢ fᵢ(z) conj fᵢ(w) for an orthonormal basis fᵢ dz of holomorphic
1-forms, in the inner product ⟨f, g⟩ = ∫_F f ḡ dx dy. Expanding the regular part to Taylor degree 7 gives

  R(z, w) = R₀ + (u(z) + u(w))/A + Σ_{p=1}^{3} β_p (z^p w̄^p + z̄^p w^p) + O(deg 8),   u = −ln(1 − |z|²),

  p² β_p = −¼ [(z w̄)^{p−1}] (B_Σ − B_𝔻),   B_𝔻 = 1/(π(1 − z w̄)²).

The notation [(z w̄)^k] means the coefficient of (z w̄)^k.

*Proof of the expansion.* Since ∂∂̄u = (1−|z|²)^{−2}, P := R − (u(z) + u(w))/A is harmonic in z and in w by (i) and the
symmetry R(z, w) = R(w, z); such a function is a sum of monomials z^p w^q, z^p w̄^q and their conjugates. The rotation
ρ: z ↦ e^{iπ/4} z and the reflection z ↦ z̄ normalise Γ (exact, `exact/group.py`), so G_Σ, G₁ and R are invariant under both:
only monomials z^a z̄^b w^c w̄^d with a − b + c − d ≡ 0 (mod 8) occur and the coefficients are real. In total degree ≤ 7
this forces charge 0, leaving constants and z^p w̄^p + c.c. The coefficients follow from the second identity because
∂_z∂_w̄(z^p w̄^p) = p²(z w̄)^{p−1}. □

At the Bolza Weierstrass point the two eigen-1-forms of the rotation ρ: z ↦ e^{iπ/4}z are f₁ = Σ cₙ z^{8n} and
f₂ = Σ dₙ z^{8n+2}. On y² = x⁵ − x these are dx/y and x dx/y, with gap sequence {1, 3}. Hence, by the expansion above:

- β₁ = ¼ (1/π − B_Σ(0,0)), the p = 1 case;
- β₂ = 1/(8π). B_Σ has no z w̄ term because there is no holomorphic 1-form of character 1 mod 8. So
  4β₂ = −¼(0 − 2/π), where 2/π is the z w̄ coefficient of B_𝔻 (shown above). This holds at any C₈ point where no
  holomorphic 1-form has character k₀ ≡ 1 (mod 8). It is not universal in higher genus;
- β₃ = −(1/36) (|d₀|²/‖f₂‖² − 3/π).

(Here "B_Σ has no z w̄ term" holds because f₁, f₂ contain only the powers 8n and 8n+2; for β₂ the p = 2 coefficient is
−¼ times the z w̄ coefficient of B_Σ − B_𝔻, and [z w̄](B_𝔻) = 2/π.)

So everything reduces to two numbers: B_Σ(0,0) = |c₀|²/‖f₁‖² and |d₀|²/‖f₂‖².

### 2. Local coordinate from the (2,3,8) Schwarz map

Work on the curve y² = x⁵ − x. The Weierstrass point x = 0 has stabiliser C₈ (x ↦ ix, y ↦ ζ₈ y), and y is a local
parameter there, with x = −y² + O(y¹⁰).

- The octahedral quotient is τ = 108 f⁴/W³, with f = x(x⁴ − 1) and W = x⁸ + 14x⁴ + 1. Klein's identity
  W³ − T² = 108 f⁴, with T = x¹² − 33x⁸ − 33x⁴ + 1, was checked in sympy. Near x = 0, τ = 108x⁴ + … = 108y⁸ + ….
- The (2,3,8) Schwarz triangle map is

  s(τ) = τ^{1/8} F(a + ⅛, b + ⅛; 9/8; τ) / F(a, b; 7/8; τ),   a = 17/48,  b = 1/48,

  where F = ₂F₁. The disk coordinate is z = K·s(τ), with the Weierstrass point at 0 and τ = (z/K)⁸(1 + O(z⁸)).
- K is fixed by putting the edge vertex τ = 1 at hyperbolic distance c from 0, where

  cosh c = cos(π/3)/sin(π/8) = 1/(2 sin(π/8)).

  So K = tanh(c/2)/s(1). Check: the image of the third side is orthogonal to |z| = 1, to 1e−31.
- s(1) follows from Gauss's formula F(A, B; C; 1) = Γ(C)Γ(C−A−B)/(Γ(C−A)Γ(C−B)), applied to numerator and
  denominator. So s(1) is a ratio of Γ-products. (Our check: C − A − B = ½ in both cases, so the formula applies.)
- Comparing τ = 108y⁸ + … with τ = (z/K)⁸ + … gives

  |a| := |dy/dz|(0) = 1/(K·108^{1/8}) = s(1) / (tanh(c/2)·108^{1/8}) = 1.767870456679306753677…

### 3. The forms and their norms

- In the local parameter, dx/y = −2 dy (1 + O(y⁸)) and x dx/y = 2y² dy (1 + O(y⁸)) (from x = −y² + O(y¹⁰): dx = −2y dy(1 + O(y⁸))). In the disk
  coordinate this gives |c₀| = 2|a| and |d₀| = 2|a|³.
- The norms are ‖ω‖² = (i/2)∫ ω ∧ ω̄ = 2 ∫_ℂ |g(x)|² d²x for ω = g(x) dx; the factor 2 comes from the two sheets.
  This gives ‖dx/y‖² = ‖x dx/y‖² = 2I₀ with

  I₀ = ∫_ℂ d²x / |x⁵ − x|.

  That the second norm equals the first: with w = x⁴ the integrals are ¼∫|w|^{−7/4}|1−w|^{−1} d²w (for dx/y) and
  ¼∫|w|^{−5/4}|1−w|^{−1} d²w (for x dx/y), and the substitution w ↦ 1/w turns the second into the first
  (d²w ↦ |w|^{−4}d²w, |w|^{−5/4} ↦ |w|^{5/4}, |1−w|^{−1} ↦ |w||1−w|^{−1}; total exponent 5/4 + 1 − 4 = −7/4).
- With the substitution w = x⁴, I₀ = ¼ ∫_ℂ |w|^{−7/4} |1 − w|^{−1} d²w. The complex beta integral then gives

  I₀ = (π/4) γ(⅛) γ(½) γ(⅜),   γ(x) = Γ(x)/Γ(1 − x),

  with I₀ = 8.973175814411088600661542424577714345855…. The evaluation uses the standard complex beta integral
  ∫_ℂ |w|^{2a−2}|1−w|^{2b−2} d²w = π γ(a)γ(b)γ(c), a + b + c = 1 (a = ⅛, b = ½, c = ⅜ here), which is quoted from
  the literature (Dotsenko–Fateev / Shapiro–Virasoro; we did not re-read a source) and not re-proved here. The resulting
  value was checked against direct quadrature to 40 digits.

### 4. The constants

f₁ and f₂ are ρ-eigenforms with different characters, so they are orthogonal, and B_Σ(0, 0) is computed from the
orthonormal ρ-eigenforms. Therefore

  B_Σ(0,0) = |c₀|²/‖f₁‖² = 4|a|²/(2I₀) = **2|a|²/I₀** = 2 s(1)² / (tanh²(c/2)·108^{1/4}·I₀),

  |d₀|²/‖f₂‖² = 4|a|⁶/(2I₀),

and so

- B_Σ(0,0) = 0.6966019648428384295921231301626841213…
- β₁ = (1/π − B)/4 = −0.09457301966476193951…
- β₂ = 1/(8π) = 0.03978873577297383394…
- β₃ = −(4|a|⁶/(2I₀) − 3/π)/36 = −0.16248381288685044074…
- the cubic coefficient of the λ₃ series, −61/2 + 126πβ₃ = −94.8177020644660778…

`bolzakit.green.constants._core` evaluates exactly these formulas in mpmath.

### 5. Independent agreement

Two independent numerical methods in the original research code agree with the closed form:
- the Hejhal-type computation of the holomorphic 1-forms, available here as `bergman_at_weierstrass()`, agrees to
  about 18 digits;
- the heat-kernel image-sum fit for β₁ agrees to 8e-12.

That is why the label is PROVED rather than NUMERICAL. The proof is the derivation above, and the numerics confirm it.
