import pytest
from bolzakit import provenance as P
from bolzakit.errors import BolzakitError, OutOfValidityRange, DataUnavailable, OptionalDependencyMissing


@P.provenance("EXACT", "test source", error=None)
def _seven():
    return 7


def test_plain_value_by_default():
    assert _seven() == 7


def test_with_provenance_returns_result():
    r = _seven(with_provenance=True)
    assert r == P.Result(7, "EXACT", "test source", None)


def test_info_and_registry():
    assert P.info(_seven) == P.Info("EXACT", "test source", None)
    assert any(k.endswith("_seven") for k in P.registry())


def test_bad_label_rejected():
    with pytest.raises(ValueError):
        P.provenance("PROBABLY", "x")


def test_require_missing_module():
    with pytest.raises(OptionalDependencyMissing, match=r"pip install bolzakit\[pari\]"):
        P.require("definitely_not_a_module_xyz", "pari")


def test_error_hierarchy():
    assert issubclass(OutOfValidityRange, BolzakitError) and issubclass(OutOfValidityRange, ValueError)
    assert issubclass(DataUnavailable, BolzakitError)
    assert issubclass(OptionalDependencyMissing, ImportError)
