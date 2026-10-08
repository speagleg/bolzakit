import math
import pytest
from bolzakit.exact import lengths as LS
from bolzakit.errors import OutOfValidityRange
from bolzakit.provenance import info


def test_first_entries_match_item02():
    L = LS.lengths()
    assert [(e["halftrace"], e["n_oriented"]) for e in L[:5]] == [
        ((1, 1), 24), ((3, 2), 24), ((5, 3), 48), ((5, 4), 0), ((7, 5), 96)]
    assert abs(L[0]["length"] - 2 * math.acosh(1 + math.sqrt(2))) < 1e-12
    assert len(L) == 42


def test_lengths_truncates_to_request():
    L = LS.lengths(max_half_trace=10)
    assert all(e["halftrace"][0] + e["halftrace"][1] * math.sqrt(2) <= 10 + 1e-9 for e in L)


def test_lengths_beyond_bundle():
    with pytest.raises(OutOfValidityRange, match="bundled up to half-trace 120"):
        LS.lengths(max_half_trace=500)


def test_bundle_max_residual():
    from bolzakit import data
    assert data.load("lengths")["max_residual"] < 1.1e-8


def test_label():
    assert info(LS.lengths).label == "NUMERICAL"


@pytest.mark.slow
def test_bundle_regenerates_small():
    rows = LS.length_spectrum(40.0)
    bundle = [e for e in LS.lengths(max_half_trace=40)]
    assert [(tuple(r["halftrace"]), r["n_oriented"]) for r in rows] == [(e["halftrace"], e["n_oriented"]) for e in bundle]


@pytest.fixture(scope="module")
def lspec():
    return LS.length_spectrum(80.0)


def test_length_spectrum(lspec):
    for r in lspec:
        assert r["residual"] < 1e-3
        assert r["n_oriented"] % 2 == 0
    head = [(r["halftrace"], r["n_oriented"] // 2) for r in lspec[:8]]
    assert head == [((1, 1), 12), ((3, 2), 12), ((5, 3), 24), ((5, 4), 0), ((7, 5), 48), ((9, 6), 24),
                    ((9, 7), 24), ((11, 8), 4)]
    # 5+4sqrt2 = T_2(1+sqrt2): only squares of systoles
    assert lspec[3]["power_of"] == [((1, 1), 2)]
    assert LS.cheb_T(2, (1, 1)) == (5, 4)
