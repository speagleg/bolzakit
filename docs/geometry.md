# Geometry of the octagon

Module: `bolzakit.geometry` (implementation in `bolzakit.geometry.octagon`). Everything uses the Poincaré disk with
curvature −1. All functions are NUMERICAL: closed-form hyperbolic geometry in float64, with error at the level of
float64 rounding.

| function | meaning |
|---|---|
| `mobius(M, z)` | (a z + b)/(c z + d) for a 2×2 matrix M |
| `dist(z, w)`, `cosh_dist(z, w)` | hyperbolic distance and its cosh |
| `reduce_point(z)` | `(z_reduced, word, M)`: z_reduced = M(z) lies in the octagon F; `word` lists the generator indices used |
| `in_F(z, tol=1e-12)` | whether z lies in the closed octagon |
| `disk_to_upper(z)`, `upper_to_disk(w)` | Cayley maps between the disk and the upper half-plane (0 ↦ i) |
| `poincare_to_klein(z)`, `klein_to_poincare(k)` | maps between the Poincaré and Klein disk models |
| `VERTICES`, `MIDPOINTS` | the 8 vertices (Poincaré radius 2^{−1/4}) and the 8 side midpoints (radius √(√2−1)) |

`reduce_point` works greedily: while some g_k(0) is closer to z than 0 is, it replaces z by g_k⁻¹(z). The result is the
Dirichlet-domain representative.

```python
import bolzakit as bk
z, word, M = bk.geometry.reduce_point(0.9 + 0.1j)
bk.geometry.in_F(z)                 # True
abs(bk.geometry.mobius(M, 0.9 + 0.1j) - z) < 1e-12
```

Useful constants in `bk.group()`: `R_C` (circumradius, cosh R_c = 3 + 2√2), `R_IN` (inradius, cosh r_in = 1 + √2),
`AREA` = 4π, `GENUS` = 2.
