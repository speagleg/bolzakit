import time

import mpmath as mp
import pytest
from bolzakit import vortex
from bolzakit.errors import OutOfValidityRange


def test_flat_modes_formula():
    lam = vortex.ring(7, 0.04, model="flat").modes()
    for m in range(1, 7):
        assert abs(lam[m] - ((7 - 1) - m * (7 - m) / 2)) < 1e-6      # (N-1) - m(N-m)/2


def test_sinh_lambda3_O_xi():
    xi = mp.mpf("1e-3")
    lam = vortex.ring(7, xi, model="sinh").modes()
    assert abs(lam[3] / xi - 12) < 0.1                                  # Montaldi–Tokieda


def test_lambda3_series():
    s = vortex.lambda3_series()
    c = s["coefficients"]
    assert c[:3] == [0, 12, -4] and abs(c[3] - mp.mpf("-94.8177")) < 1e-3


def test_bolza_ring_out_of_range():
    t0 = time.time()
    with pytest.raises(OutOfValidityRange):
        vortex.ring(7, 0.3, model="bolza").modes()       # ring radius sqrt(0.3) = 0.55 > 0.45
    assert time.time() - t0 < 2.0                         # raised before building the Green function


def test_bad_model_and_N():
    with pytest.raises(ValueError):
        vortex.ring(7, 0.01, model="cosh")
    with pytest.raises(ValueError):
        vortex.ring(1, 0.01, model="flat")


@pytest.mark.parametrize("xi", [-0.1, 0, 1.5])
def test_xi_out_of_unit_interval(xi):
    with pytest.raises(ValueError):
        vortex.ring(7, xi, model="tanh")


def test_xi_validation_complex_and_bolza_before_green():
    with pytest.raises(ValueError):
        vortex.ring(7, mp.mpc(0.01, 0.1), model="tanh")
    with pytest.raises(ValueError):
        vortex.ring(7, -0.1, model="bolza")
    r = vortex.ring(7, mp.mpc(0.01, 0), model="flat")
    assert r.xi == mp.mpf(0.01) and r.split(3) == 0.0


@pytest.mark.slow
def test_bolza_ring_matches_series():
    xi = 0.01
    r = vortex.ring(7, xi, model="bolza", D_ball=13.0)           # item 03's setting; lambda3 moves 3.7e-9 vs 15
    lam = r.modes()
    s = vortex.lambda3_series()["coefficients"]
    pred = sum(float(c) * xi ** k for k, c in enumerate(s))
    print(f"lam3 (symmetrised) = {lam[3]!r}  series = {pred!r}  diff = {lam[3] - pred:+.4e}  "
          f"split(3) = {r.split(3):+.4e}  split/xi^4 = {r.split(3) / xi ** 4:.2f}")
    assert abs(lam[3] - pred) < 5e-6                                # (a)
    assert 900 < r.split(3) / xi ** 4 < 1100                        # (b) C8-vs-C7 symmetry breaking


class _FakeGreen:
    """Cheap smooth, non-rotation-invariant stand-in for BolzaGreen.S_ring (exercises every term)."""

    def S_ring(self, zs):
        import numpy as np
        zs = np.asarray(zs)
        return float(np.sum(np.abs(zs) ** 4) + np.real(np.sum(zs ** 8)) + 0.3 * np.abs(np.sum(zs)) ** 2)


def test_bolza_lam_sigma_equals_item03_analyse():
    from bolzakit.vortex import dynamics as D
    g, xi = _FakeGreen(), 0.02
    fast = D.lam_sigma(g, xi, (1, 3), dps=15)
    for m in (1, 3):
        ref = D.analyse(g, xi, m=m, dps=15)
        assert abs(fast[m][0] - ref["lam_sigma"]) < 1e-12
        assert abs(fast[m][1] - ref["lam_sigma_sin"]) < 1e-12


@pytest.mark.parametrize("kw", [dict(T=5.0), dict(D_ball=9.0), dict(T=10.0)])
def test_bolza_ring_unvalidated_T_D_raises_fast(kw):
    t0 = time.time()
    with pytest.raises(OutOfValidityRange, match="allowed pairs"):
        vortex.ring(7, 0.01, "bolza", **kw)
    assert time.time() - t0 < 1.0
    vortex.ring(7, 0.01, "tanh", **kw)        # T, D_ball are irrelevant for the closed-form models
