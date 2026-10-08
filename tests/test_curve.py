import pytest
from bolzakit.arith import curve
from bolzakit.errors import OptionalDependencyMissing


def test_missing_pari_message(monkeypatch):
    import bolzakit.provenance as P
    real = P.require
    monkeypatch.setattr(curve, "require", lambda m, e: real("no_such_mod_zz", e))
    with pytest.raises(OptionalDependencyMissing, match=r"bolzakit\[pari\]"):
        curve.euler_factor(3)


def test_euler_factor_p3():
    pytest.importorskip("cypari2")
    assert curve.euler_factor(3) == [1, 0, -2, 0, 9]          # x^4 - 2x^2 + 9 (item 15c, PARI)


def test_jacobian_Q_simple():
    pytest.importorskip("cypari2")
    assert curve.jacobian_is_Q_simple() is True


@pytest.mark.parametrize("bad", [2, 0, -3, 1, 4, 9, 15, 3.0, "3", True])
def test_jacobian_invalid_prime_raises(bad):
    with pytest.raises(ValueError):
        curve.jacobian_is_Q_simple(primes=(bad,))


@pytest.mark.parametrize("bad", [2, 1, 4, 9, 15, True])
def test_euler_factor_invalid_prime_raises(bad):
    with pytest.raises(ValueError):
        curve.euler_factor(bad)
