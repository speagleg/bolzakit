# Quick-start

```
pip install bolzakit
```

```python
import bolzakit as bk

G = bk.group()
len(G.GENS)                                   # 8 generators, exact in Z[zeta8]
bk.automorphisms().order                      # 48 (GL(2,3)); .full_order is 96
bk.lengths(max_half_trace=5)[0]               # {'halftrace': (1, 1), 'length': 3.0571..., 'n_oriented': 24, 'n_unoriented': 12}
bk.systole()                                  # 3.057141838961996
bk.spectral.eigenvalues("fem")[0]             # (3.83889..., 3)
bk.green.constants()["B_Sigma_00"]["value"]   # 0.6966019648428384...
```

## Provenance

```python
r = bk.systole(with_provenance=True)          # Result(value, label, source, error)
r.label                                       # 'EXACT'
from bolzakit.provenance import info
info(bk.lengths)                              # Info(label='NUMERICAL', source=..., error='integrality residual <= 1.1e-8 ...')
```

The labels are EXACT, PROVED, CERTIFIED, NUMERICAL and QUOTED; see the README. Every labelled function is listed in
[provenance.md](provenance.md).

## Errors

All errors derive from `bolzakit.errors.BolzakitError`.

| error | when |
|---|---|
| `OutOfValidityRange` (also a `ValueError`) | `lengths` beyond half-trace 120 without `compute=True`; `green().G` for a pair with no common Weierstrass chart; a bolza ring with √ξ > 0.45 |
| `DataUnavailable` | the Strohmaier–Uski download or checksum failed; more FEM clusters requested than bundled |
| `OptionalDependencyMissing` (also an `ImportError`) | `arith.curve` without cypari2 |

## Slow and networked calls

- `bk.green()` builds a heat-kernel image ball on first use (seconds to a minute) and caches it per `(T, D_ball)`.
- `bk.vortex.ring(..., "bolza")` builds the Green's function and differentiates it. It takes about 20 s at
  `D_ball=13` and a few minutes at the default 15.
- `bk.spectral.su2013()`, `eigenvalues("su2013")` and `selberg_check` need the Strohmaier–Uski file. It is downloaded
  once into `~/.cache/bolzakit` (override with `$BOLZAKIT_CACHE`).
