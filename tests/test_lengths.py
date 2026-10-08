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


def test_compute_rejects_bad_integrality(monkeypatch):
    # final review M1: compute=True raises if any multiplicity residual > 1e-3 (stubbed; no enumeration)
    from bolzakit.errors import BolzakitError
    rows = [dict(halftrace=(121, 0), length=1.0, n_oriented=24, residual=1e-9),
            dict(halftrace=(100, 20), length=2.0, n_oriented=48, residual=0.02)]
    monkeypatch.setattr(LS, "length_spectrum", lambda T: rows)
    with pytest.raises(BolzakitError, match=r"residual 0.02 > 1e-3 at half-trace 100\+20\*sqrt2"):
        LS.lengths(max_half_trace=130, compute=True)
    rows[1]["residual"] = 1e-6
    assert [e["n_unoriented"] for e in LS.lengths(max_half_trace=130, compute=True)] == [12, 24]


def test_two_halftraces_without_primitive_geodesic():
    # final review M6: of the 42 bundled half-traces, 5+4√2 and 33+24√2 have n_oriented = 0
    L = LS.lengths()
    assert len(L) == 42
    assert [e["halftrace"] for e in L if e["n_oriented"] == 0] == [(5, 4), (33, 24)]
