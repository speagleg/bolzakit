# Vortex rings

Module: `bolzakit.vortex`. Source: math-pursuits item 03 (`note.tex`, `ringmodes.py`, `dynamics.py`).

## Model

N identical unit point vortices sit on a regular N-gon of Euclidean radius r = √ξ in the disk, centred at 0. For the
`bolza` model, 0 is a Weierstrass point of the surface. The Hamiltonian is the Boatto–Koiller one [BK08]: the Green's
function with the −1/A background for each pair, plus ½ the Robin function for each vortex. In their words, "a single
vortex so drifts according to the Hamiltonian system (Ωg , Rg ), where Ωg is the area form of g and Rg is Robins's
function", and the collective Hamiltonian is H = Σ_{i<j} κᵢκⱼ G(sᵢ, sⱼ) + Σ ½ κₗ² R(sₗ). bolzakit uses G_Σ = −G_BK,
which flips the sign of H; that is time reversal and does not affect spectral stability (item 03 review_1 §1.1, by an AI second reader).

Units are those of item 03, "C1" (c_norm = 1/2π).

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
  - bolza: item 03's construction. The tanh part is done in mpmath, the float64 smooth part by central differences
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

## What is and is not proved

The λ₃ series is a theorem (the main theorem of item 03's note) for the **SO(2)-invariant degree-7 truncation of the
Boatto–Koiller model**, in which the 7-ring is an exact relative equilibrium. In the full Bolza Hamiltonian the
charge-8 Taylor terms make the ring only an approximate relative equilibrium, with relative residual O(ξ⁴); `split` is
the visible trace of this. Persistence and stability there are open (note, Remark "full Bolza Hamiltonian: open"). The
Boatto–Koiller model also freezes the 2g harmonic degrees of freedom (as item 03's note records), so every stability
statement here is conditional on that model.

## References

- [BK08] S. Boatto, J. Koiller, *Vortices on closed surfaces*, arXiv:0802.4313.
- [MT13] J. Montaldi, T. Tokieda, *Deformation of geometry and bifurcations of vortex rings*, in Recent Trends in
  Dynamical Systems, Springer Proc. Math. Stat., 335–370 (2013); arXiv:1210.5662.
