# Green's function, Robin mass and Bergman constants

Module: `bolzakit.green`. `bk.green` is both a module and a callable.

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
- A pair of points with no common chart raises `OutOfValidityRange`.
- `G(z, z)` raises `ValueError`; use `robin_mass(z)` or `regular(z, z)` instead.
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

`constants()` is labelled PROVED, and its source field cites this page. The derivation below reproduces the source.
It is due to an independent AI second reader (Claude) in math-pursuits item 03. No human specialist has checked it yet.
Sources:
- `items/03-compact-surface-n7-stability/work/review_1.md` §3.A;
- `note.tex`: Proposition "structure through degree 7" in §2 and the paragraph "Closed form" in §5;
- the second reader's scripts `review_1_scripts/schwarz_bergman.py` and `closed_form.py`.

Steps that the source does not spell out are marked *see source*. bolzakit adds nothing to the argument.

### 1. What the constants are

Lift the Weierstrass point to 0 in the disk. Item 03 (note.tex, Lemma "two PDE identities") shows that off the
diagonal ∂_z∂_w̄ G_Σ = −¼ B_Σ. Here B_Σ(z, w) = Σᵢ fᵢ(z) conj fᵢ(w) for an orthonormal basis fᵢ dz of holomorphic
1-forms, in the inner product ⟨f, g⟩ = ∫_F f ḡ dx dy. Expanding the regular part to Taylor degree 7 gives

  R(z, w) = R₀ + (u(z) + u(w))/A + Σ_{p=1}^{3} β_p (z^p w̄^p + z̄^p w^p) + O(deg 8),   u = −ln(1 − |z|²),

  p² β_p = −¼ [(z w̄)^{p−1}] (B_Σ − B_𝔻),   B_𝔻 = 1/(π(1 − z w̄)²).

The notation [(z w̄)^k] means the coefficient of (z w̄)^k.

At the Bolza Weierstrass point the two eigen-1-forms of the rotation ρ: z ↦ e^{iπ/4}z are f₁ = Σ cₙ z^{8n} and
f₂ = Σ dₙ z^{8n+2}. On y² = x⁵ − x these are dx/y and x dx/y, with gap sequence {1, 3}. Hence (note.tex, Prop. 3):

- β₁ = ¼ (1/π − B_Σ(0,0)), the p = 1 case;
- β₂ = 1/(8π). B_Σ has no z w̄ term because there is no holomorphic 1-form of character 1 mod 8. So
  4β₂ = −¼(0 − 2/π), where 2/π is the z w̄ coefficient of B_𝔻 (review_1 §1.2). This holds at any C₈ point where no
  holomorphic 1-form has character k₀ ≡ 1 (mod 8). It is not universal in higher genus;
- β₃ = −(1/36) (|d₀|²/‖f₂‖² − 3/π).

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

- In the local parameter, dx/y = −2 dy (1 + O(y⁸)) and x dx/y = 2y² dy (1 + O(y⁸)) (`schwarz_bergman.py`). In the disk
  coordinate this gives |c₀| = 2|a| and |d₀| = 2|a|³.
- The norms are ‖ω‖² = (i/2)∫ ω ∧ ω̄ = 2 ∫_ℂ |g(x)|² d²x for ω = g(x) dx; the factor 2 comes from the two sheets.
  This gives ‖dx/y‖² = ‖x dx/y‖² = 2I₀ with

  I₀ = ∫_ℂ d²x / |x⁵ − x|.

  That the second norm equals the first is shown in the source by the substitution r ↦ 1/r (*see source*).
- With the substitution w = x⁴, I₀ = ¼ ∫_ℂ |w|^{−7/4} |1 − w|^{−1} d²w. The complex beta integral then gives

  I₀ = (π/4) γ(⅛) γ(½) γ(⅜),   γ(x) = Γ(x)/Γ(1 − x),

  with I₀ = 8.973175814411088600661542424577714345855…. The evaluation of the complex beta integral itself is not
  written out in the source (*see source*). It was checked against direct quadrature to 40 digits.

### 4. The constants

f₁ and f₂ are ρ-eigenforms with different characters, so they are orthogonal; the source works with the orthonormal
ρ-eigenforms (review_1 §1.2). Therefore

  B_Σ(0,0) = |c₀|²/‖f₁‖² = 4|a|²/(2I₀) = **2|a|²/I₀** = 2 s(1)² / (tanh²(c/2)·108^{1/4}·I₀),

  |d₀|²/‖f₂‖² = 4|a|⁶/(2I₀),

and so

- B_Σ(0,0) = 0.6966019648428384295921231301626841213…
- β₁ = (1/π − B)/4 = −0.09457301966476193951…
- β₂ = 1/(8π) = 0.03978873577297383394…
- β₃ = −(4|a|⁶/(2I₀) − 3/π)/36 = −0.16248381288685044074…
- the cubic coefficient of the λ₃ series, −61/2 + 126πβ₃ = −94.8177020644660778…

`bolzakit.green.constants._core` evaluates exactly these formulas in mpmath. They match the AI second reader's
`closed_form.py` line for line.

### 5. Independent agreement

Two independent numerical methods in item 03 agree with the closed form:
- the Hejhal-type computation of the holomorphic 1-forms, available here as `bergman_at_weierstrass()`, agrees to
  about 18 digits;
- the heat-kernel image-sum fit for β₁ agrees to 8e-12.

That is why the label is PROVED rather than NUMERICAL. The proof is the derivation above, and the numerics confirm it.
