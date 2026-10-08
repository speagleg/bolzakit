import mpmath as mp
import pytest

from bolzakit.exact import lengths as LS

pytestmark = pytest.mark.network


@pytest.fixture(scope="module")
def lspec():
    return LS.length_spectrum(80.0)


def test_selberg_balance(lspec):
    from bolzakit.spectral import selberg as SE
    rows = SE.balance([0.1, 0.25, 0.5, 1.0], spectrum_rows=lspec, dps=25)
    for r in rows:
        tol = max(1e-14, 10 * float(r["geometric_tail_estimate"]))
        assert abs(float(r["residual"])) < tol, (r["t"], r["residual"])


def test_selberg_detects_wrong_convention(lspec):
    """Unoriented counts (half) must FAIL the balance -> confirms the oriented convention."""
    from bolzakit.spectral import selberg as SE
    half = [dict(r, n_oriented=r["n_oriented"] // 2) for r in lspec]
    r = SE.balance([0.5], spectrum_rows=half, dps=20)[0]
    assert abs(float(r["residual"])) > 1e-3


def test_heat_constant_term():
    from bolzakit.spectral import selberg as SE
    # identity term = 1/t - 1/3 + O(t): chi/6 = -1/3, hence zeta(0) = -1/3 - 1 = -4/3
    for t in (1e-3, 2e-3):
        it = SE.identity_term(t)
        assert abs(float(it - 1 / mp.mpf(t)) + 1 / 3) < 2 * t
    rows = SE.heat_constant_term(ts=(0.04,))
    assert abs(float(rows[0]["theta_minus_inv_t"]) + 1 / 3) < 0.04 * 0.1 + 1e-6


def test_selberg_check_small_residual():
    from bolzakit.spectral.selberg import selberg_check
    assert selberg_check(0.5) < 1e-12
