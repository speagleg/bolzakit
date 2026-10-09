import pytest
flint = pytest.importorskip("flint"); pytest.importorskip("sympy")
from bolzakit.certified import count_eigenvalues


def test_no_eigenvalues_below_3_3():
    r = count_eigenvalues("bolza", window=(0, 3.3))
    assert r.certified and r.bound < 1 and abs(r.bound - 0.912612) < 5e-3


def test_at_most_three_in_3_3_4_0():
    r = count_eigenvalues("bolza", window=(3.3, 4.0))
    assert r.certified and r.bound < 4 and abs(r.bound - 3.869661) < 5e-3


def test_window_containing_lambda4_fails():
    r = count_eigenvalues("bolza", window=(3.3, 5.4), target=4)
    assert not r.certified and r.bound >= 4


def test_other_surfaces_rejected():
    with pytest.raises(ValueError, match="Bolza only"):
        count_eigenvalues("klein_quartic", window=(0, 1))


def test_precision_not_changed_globally():
    before = flint.ctx.prec
    count_eigenvalues("bolza", window=(0, 3.3))
    assert flint.ctx.prec == before


def test_mode_b_carries_lambda1_hypothesis_mode_a_does_not():
    rb = count_eigenvalues("bolza", window=(3.3, 4.0))
    ids = [d.id for d in rb.certificate.dependencies]
    assert "H_a" in ids and rb.certificate.inputs["assumes"] == "lambda_1 >= 33/10"
    ra = count_eigenvalues("bolza", window=(0, 3.3))
    assert "H_a" not in [d.id for d in ra.certificate.dependencies]
    assert "assumes" not in ra.certificate.inputs


def test_certified_window_contains_requested_window():
    from fractions import Fraction
    r = count_eigenvalues("bolza", window=(0, 3.8388757))
    lo, hi = r.certificate.inputs["window_exact"]
    assert Fraction(lo) == 0 and Fraction(hi) >= Fraction("3.8388757")
    assert Fraction(r.certificate.data["a"]) >= Fraction("3.8388757")


def test_nan_piece_forces_failure(monkeypatch):
    import flint
    from bolzakit.certified import counting
    real = counting.arb_poly
    calls = [0]

    class Wrap:
        def __init__(self, c):
            self.p = real(c)

        def __call__(self, x):
            calls[0] += 1
            if calls[0] == 10:
                return flint.arb("nan")
            return self.p(x)

    monkeypatch.setattr(counting, "arb_poly", Wrap)
    r = count_eigenvalues("bolza", window=(0, 3.3))
    assert not r.certified
