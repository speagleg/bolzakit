# Spectral data

Module: `bolzakit.spectral`. The eigenvalue functions live in `bolzakit.spectral.eigen`.

The Bolza surface here has curvature −1 and area 4π, and Δ is the positive Laplace–Beltrami operator. All eigenvalues
are the non-zero ones; λ₀ = 0 is simple.

## API

| call | returns | label |
|---|---|---|
| `eigenvalues(source="fem", n=None)` | list of `(value, multiplicity)` clusters, ascending | NUMERICAL (see below) |
| `su2013(path=None, dps=50)` | the 1000 Strohmaier–Uski eigenvalues as mpmath numbers | QUOTED |
| `selberg_check(t, T_max=80.0)` | \|geometric − spectral side\| of the heat-trace Selberg formula at time t | NUMERICAL |
| `zeta_zero()` | ζ_Δ(0) = −4/3, as a `Fraction` | PROVED |
| `tau_bounds(case="systole_shell")` | `{"lower", "upper", "h", "details"}` for τ(P₋(h)) | CERTIFIED |
| `weyl_counting(lam)` | Weyl-law leading term A·λ/4π = λ (asymptotic only) | NUMERICAL |

`eigenvalues("fem", n)` with n > 9 raises `DataUnavailable`. So does `eigenvalues("su2013", n)` with n beyond the number
of clusters in the file. An unknown `source` raises `ValueError`.

## Two sources, one label: why `su2013` values are QUOTED

`eigenvalues` is a single function, so its decorator gives it a single label, NUMERICAL. That label describes the
**FEM** branch, the bundled data. The `source="su2013"` branch returns numbers we did not compute. They are the
eigenvalues published by Strohmaier and Uski [SU13], and they are **QUOTED**:

- `su2013()` downloads the authors' ancillary file
  `https://arxiv.org/src/1110.2150v4/anc/eig-bolza-refined0-1000.txt`. It checks the sha256
  (`869f48bc…4520345`) and caches the file in `~/.cache/bolzakit`, or in `$BOLZAKIT_CACHE` if that is set. The file is
  never bundled. If the download or the checksum fails, it raises `DataUnavailable` with instructions; you can then
  pass `path=` to a file you downloaded yourself (the checksum is still enforced).
- `eigenvalues("su2013")` groups those values into clusters (relative gap ≤ 1e-9) and returns them as floats.
  Grouping does not change the values, so the values remain QUOTED even though the function's label says NUMERICAL.
  If you need the label to travel with the data, call `su2013(with_provenance=True)`.
- **Accuracy as stated by the authors.** [SU13] §7.1 says that "To compute the first 500 eigenvalues of the Bolza surface
  to a precision of 12 digits about 10000 λ-evaluations ... were needed". For λ₁ they give
  3.8388872588421995185866224504354645970819150157 "where we believe that all digits are correct". They checked the
  completeness of their eigenvalue lists with the Selberg trace formula applied to the heat trace. bolzakit makes no
  stronger claim than theirs.

## The FEM data: a 9-cluster cross-check

The bundled `eigenvalues_fem.json` holds 9 clusters (multiplicities 3, 4, 2, 4, 3, 3, 4, 1, 3). They come from P1
finite elements on meshes n = 40 and 80, with h² Richardson extrapolation of the cluster means (item 02).

- **Role.** The FEM data exists to cross-check the quoted data independently, and to give an offline answer for the
  lowest clusters. It is not meant to replace [SU13].
- **Agreement.** All 9 cluster multiplicities match [SU13], and the values agree to a relative 2.98e-5 (≤ 3e-5).
- **Label.** The NUMERICAL label carries a more conservative error bound, "relative error ≤ 2e-4 vs Strohmaier–Uski".
  The tests check that bound on the first 5 clusters.

## Selberg trace formula

`selberg_check(t)` evaluates both sides of the heat-trace Selberg formula, [SU13] eq. (18):
- the **spectral side** is 1 + Σ e^{−tλ} over the 1000 quoted eigenvalues (so it needs the downloaded file);
- the **geometric side** is the identity term plus the hyperbolic term over bolzakit's length spectrum, with
  half-traces up to `T_max` (default 80, enumerated on the fly).

The sum runs over *oriented* primitive geodesics. With unoriented counts the residual is O(1e-2), which confirms the
convention. Typical residuals are ~1e-16 at t = 0.25 and 0.5, and ~5e-11 at t = 1. A call takes a few seconds.

**Validity window.** `selberg_check` raises `OutOfValidityRange` unless both truncation errors are ≤ 1e-10: the Weyl
bound e^{−λ₁₀₀₀t}/t on the spectral tail beyond the 1000 eigenvalues, and the prime-geodesic estimate of the
geometric tail beyond half-trace `T_max`. With the bundled data at T_max = 80 the window is about **0.03 ≤ t ≤ 1**:

| t | 0.02 | 0.025 | 0.03 | 0.04 | 0.5 | 1.0 | 1.2 | 2.0 |
|---|---|---|---|---|---|---|---|---|
| spectral tail bound | 1.1e-7 | 5.8e-10 | 3.3e-12 | 1.2e-16 | ~0 | 0 | 0 | 0 |
| geometric tail estimate | 0 | 0 | 0 | ~0 | 4.9e-22 | 9.8e-11 | 7.6e-9 | 4.6e-5 |
| residual | 9.9e-8 | 5.3e-10 | 3.6e-12 | 2.5e-14 | 6.6e-17 | 5.3e-11 | 4.1e-9 | 2.4e-5 |

## ζ_Δ(0)

The heat trace has constant term χ/6 = −1/3. This is the identity term of the Selberg formula. Removing the zero mode
(dim ker Δ = 1) gives ζ_Δ(0) = −1/3 − 1 = −4/3. The value is exact.

## τ bounds

`tau_bounds` returns certified two-sided bounds on τ(P₋(h)), the von Neumann trace of the negative spectral projection.
Here h is a self-adjoint element of the group algebra of the Bolza group.
- `case="systole_shell"` uses h = c·Σ_{k=0}^{7} u_{g_k} with c = G₁(ℓ_sys). The bounds are [0.40993…, 0.5].
- The integer cases 0..3 use G₁-weighted ball operators; see `details`.

The method has two parts:
- exact group moments τ(h^k). These are integers for the systole shell; for float coefficients there is a rigorous
  round-off budget;
- self-certifying polynomial bounds, checked on all of [−‖h‖, ‖h‖].

The bundled numbers can be regenerated with `python scripts/regen_tau.py` (about 30 min and 8 GB of RAM). The
strict-fractionality theorem in [results-and-corrections.md](results-and-corrections.md) shows that the true value lies
strictly between 0 and 1.

## Reference

[SU13] A. Strohmaier, V. Uski, *An algorithm for the computation of eigenvalues, spectral zeta functions and
zeta-determinants on hyperbolic surfaces*, Commun. Math. Phys. 317, 827–869 (2013); arXiv:1110.2150 (v4 carries the
ancillary data file).
