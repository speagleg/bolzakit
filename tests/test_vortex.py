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


@pytest.mark.slow
def test_bolza_ring_matches_series():
    xi = 0.01
    lam = vortex.ring(7, xi, model="bolza").modes()
    s = vortex.lambda3_series()["coefficients"]
    pred = sum(float(c) * xi ** k for k, c in enumerate(s))
    assert abs(lam[3] - pred) < 5e-6


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
        assert abs(fast[m] - D.analyse(g, xi, m=m, dps=15)["lam_sigma"]) < 1e-12
