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
