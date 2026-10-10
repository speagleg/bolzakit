# Vortex rings

Module: `bolzakit.vortex`. Source: ported from the original research code (private archive, not public); the arguments are restated here.

## Model

N identical unit point vortices sit on a regular N-gon of Euclidean radius r = √ξ in the disk, centred at 0. For the
`bolza` model, 0 is a Weierstrass point of the surface. The Hamiltonian is the Boatto–Koiller one [BK08]: the Green's
function with the −1/A background for each pair, plus ½ the Robin function for each vortex. In their words, "a single
vortex so drifts according to the Hamiltonian system (Ωg , Rg ), where Ωg is the area form of g and Rg is Robins's
function", and the collective Hamiltonian is H = Σ_{i<j} κᵢκⱼ G(sᵢ, sⱼ) + Σ ½ κₗ² R(sₗ). bolzakit uses G_Σ = −G_BK,
which flips the sign of H; that is time reversal and does not affect spectral stability (the AI second reader's review in the original research code, §1.1).

Units are those of the original research code, "C1" (c_norm = 1/2π).

| `model` | kernel |
|---|---|
| `"flat"` | the plane |
| `"sinh"` | −(1/4π) ln sinh²(d/2): the Montaldi–Tokieda kernel [MT13, eq. (1.12)], log(\|z−w\|²/((1+λ\|z\|²)(1+λ\|w\|²))), at λ = −1, up to the factor −1/4π |
| `"tanh"` | G₁ = −(1/2π) ln tanh(d/2), the free Green's function of H², i.e. the local model of G_Σ |
| `"bolza"` | the Bolza surface, H_Σ = H_tanh + the smooth part from `bk.green(T, D_ball)` |

Montaldi and Tokieda construct "a smooth family of Hamiltonian systems, together with a family of group symmetries and
momentum maps, for the dynamics of point vortices on surfaces parametrized by the curvature of the surface" [MT13].

## API

- `bk.vortex.ring(N, xi, model="bolza", dps=30, *, T=2.6, D_ball=15.0)` (NUMERICAL) returns a `Ring`. Arguments are
  validated before any expensive work: 0 < ξ < 1 and N ≥ 2. For `bolza` it also requires √ξ ≤ 0.45 and a validated
  `(T, D_ball)` pair, (2.6, 15.0) or (2.6, 13.0) (see [green.md](green.md)), and otherwise raises `OutOfValidityRange`.
- `Ring.modes()` returns `{m: λ_m}` for m = 1..N−1. λ_m is the radial-mode Rayleigh quotient of Hess(H − ΩJ), with
  λ_{N−m} = λ_m.
  - flat, sinh and tanh: computed in mpmath by central differences with step 10^(−dps/3).
  - bolza: the original code's construction. The tanh part is done in mpmath, the float64 smooth part by central differences
    with step 2e-3, plus the Ω correction. It returns the **mean of the cos(mθ_k) and sin(mθ_k) copies**, i.e. the
    C₇-symmetrised eigenvalue that the series below describes.
- `Ring.split(m)` returns cos copy − sin copy. It is 0.0 for flat, sinh and tanh, which are rotation invariant. For
  bolza it is the O(ξ⁴) breaking between the surface's C₈ symmetry at the Weierstrass point and the ring's C₇: about
  1000·ξ⁴ for m = 3.
- `bk.vortex.lambda3_series()` (PROVED) returns `{"coefficients": [0, 12, −4, −61/2 + 126πβ₃], ...}`, i.e.

  λ₃ = 12ξ − 4ξ² + (−61/2 + 126πβ₃)ξ³ + O(ξ⁴) = 12ξ − 4ξ² − 94.8177ξ³ + O(ξ⁴),

  for N = 7 at a Bolza Weierstrass point. β₃ is in closed form; see [green.md](green.md).

## Numbers (N = 7, ξ = 0.01, mode m = 3)

| model | λ₃ |
|---|---|
| flat | 0 (marginal: (N−1) − m(N−m)/2 = 0) |
| sinh | 0.1224365 |
| tanh | 0.1209974 |
| bolza (D_ball = 13) | 0.1195036; split(3) = 1.006e-5 ≈ 1006 ξ⁴ |
| bolza (D_ball = 15) | 0.1195036; split(3) = 1.015e-5 |
| series | 0.1195052 |

The bolza value differs from the series by −1.6e-6. That is the O(ξ) heat-kernel-tail offset at T = 2.6 (about −1.5e-6
at ξ = 0.01, as recorded in the label's error field). The O(h²) differencing error is negligible. On the development
machine a bolza ring took about 20 s at D_ball = 13 and about 3 min at the default D_ball = 15, including the build of
the Green's function.

## Derivation of the λ₃ series

This is the argument behind the PROVED label of `lambda3_series`, written out here. It is an AI-written argument
that no human specialist has checked; the pieces it quotes are [BK08] and the Montaldi–Tokieda formula [MT13].

**Setup.** Write z_k = e^{s_k + iθ_k} for the vortices, a ring s_k = s₀, θ_k = 2πk/N, ξ = e^{2s₀}. The symplectic form is
Σ κ(s_k) ds_k ∧ dθ_k with κ = 4ξ/(1−ξ)², and the rotation momentum is J = Σ j(s_k) with j = 2e^{2s}/(1−e^{2s}), j' = κ.
For a rotation-invariant H the ring is a relative equilibrium of F = H − ΩJ with Ω = ∂_sH/j'. The radial and tangential
eigenvalues of Hess F in mode m are α_m and β_m; λ_m = α_m/c_norm with c_norm = 1/2π. In a mode block the only coupling
allowed by the reflection θ ↦ −θ is Coriolis-type, so mode m is spectrally stable iff α_mβ_m > 0, and β_m =
c_norm m(N−m)/2 · (1 + O(ξ)) > 0 for small ξ.

**Hamiltonian.** With G_Σ = −G_BK, H = Σ_{i<j} G_Σ(z_i, z_j) + ½ Σ m(z_i) = H_tanh + S + const, where G₁ = −(1/2π) ln tanh(d/2)
is the free Green's function, H_tanh = Σ_{i<j} G₁(d_ij), and S = Σ_{i<j} R(z_i, z_j) + ½ Σ_i R(z_i, z_i) with R = G_Σ − G₁
(Robin mass m(z) = R(z, z) + ln 2/2π). At a point with stabiliser C₈ (the Bolza Weierstrass point), [green.md](green.md)
§1 gives R through Taylor degree 7, and hence
S = const + (N/A) Σ_k u(z_k) + Σ_{p=1}^{3} β_p |Σ_k z_k^p|² + O(deg 8), with u = −ln(1 − |z|²) and A = 4π.

**Lemma (momentum map).** Let H be rotation invariant with the ring a relative equilibrium, and Φ = Σ_k φ(s_k). Then the ring
(same radius) is a relative equilibrium of H + εΦ with Ω' = Ω + εφ'/j', each α_m changes by ε(φ'' − φ'j''/j'), and the
β_m and the mode couplings are unchanged (Hess F' = Hess F + ε diag(φ'') − ε(φ'/j') diag(j'') on the s-block, a multiple
of the identity on every mode). For φ = ln(1 + j/2) = u the change is exactly −4εξ²/(1−ξ)².

**Lemma (power sums).** Φ_p = |Σ_k z_k^p|² has vanishing gradient at the ring when N ∤ p, and its Hessian is nonzero only
on modes m ≡ ±p (mod N), of size O(ξᵖ) in C1 units (δΣ z_k^p = p r^p Σ_k e^{ipθ_k}(δs_k + iδθ_k) picks only the Fourier
components ∓p).

**Exact identity for the tanh model.** Expanding ln|1 − z w̄|² = −Σ_{p≥1} ((z w̄)^p + c.c.)/p gives, for |z_k| < 1,
H_tanh = H_sinh + (N/4π) Σ_k u(z_k) − Σ_{p≥1} (1/4πp) |Σ_k z_k^p|².
For H_sinh [MT13] gives λ_m = C₁ − m(N−m)/2 with C₁ = (N−1)(1 + ξ²)/(1 − ξ)². Apply the two lemmas to the u term and to
the power-sum terms (only p ≡ ±3 mod 7 matter for m = 3, N = 7; the p ≡ 0 mod 7 terms are O(ξ⁷)). For N = 7, m = 3,
λ₃^tanh = 6(1+ξ²)/(1−ξ)² − 6 − 14ξ²/(1−ξ)² − (7/2) Σ_{p≡±3 (7)} p ξᵖ + O(ξ⁷) = 12ξ + 10ξ² − (5/2)ξ³ + O(ξ⁴).

**The Bolza correction.** The truncation H^{≤7} = H_tanh + const + (N/A) Σ u(z_k) + Σ_{p≤3} β_p |Σ z_k^p|² is the
SO(2)-average of H up to Taylor degree 7 of S, and the regular N-ring is an exact relative equilibrium of it. The two
lemmas give, for m ∉ {1, 2, 3, N−3, N−2, N−1}, λ_m = λ_m^tanh − (8πN/A) ξ²/(1−ξ)², and for m ∈ {p, N−p}, p ≤ 3, the
extra term 2πN p² β_p ξᵖ. For N = 7, m = 3, A = 4π:

  λ₃ = 12ξ + (10 − 14)ξ² + (−5/2 − 28 + 126πβ₃)ξ³ + O(ξ⁴) = 12ξ − 4ξ² + (−61/2 + 126πβ₃)ξ³ + O(ξ⁴).

The O(ξ) coefficient 12 is the same as for the sinh model. Since α₁, α₂ → 3, 1 and β_m > 0, the 7-ring of the truncated
model is spectrally stable modulo rotations for small ξ (α_mβ_m > 0, m = 1, 2, 3). The theorem concerns only the
truncation and the Boatto–Koiller model; see the next section. □

## What is and is not proved

The λ₃ series is a theorem (proved in the derivation above) for the **SO(2)-invariant degree-7 truncation of the
Boatto–Koiller model**, in which the 7-ring is an exact relative equilibrium. In the full Bolza Hamiltonian the
charge-8 Taylor terms make the ring only an approximate relative equilibrium, with relative residual O(ξ⁴); `split` is
the visible trace of this. Persistence and stability there are open. The
Boatto–Koiller model also freezes the 2g harmonic degrees of freedom (recorded in the original research code; we have not checked this against the literature, see [CITATIONS_CHECK.md](CITATIONS_CHECK.md)), so every stability
statement here is conditional on that model.

## References

- [BK08] S. Boatto, J. Koiller, *Vortices on closed surfaces*, arXiv:0802.4313.
- [MT13] J. Montaldi, T. Tokieda, *Deformation of geometry and bifurcations of vortex rings*, in Recent Trends in
  Dynamical Systems, Springer Proc. Math. Stat., 335–370 (2013); arXiv:1210.5662.
