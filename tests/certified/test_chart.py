# tests/certified/test_chart.py
import pytest
flint = pytest.importorskip("flint")
from flint import arb
from bolzakit.certified import Surface
from bolzakit.certified.chart import certify_rank, traces
from bolzakit.errors import OutOfValidityRange

def test_rank6_golden():
    r = certify_rank()
    assert r.certified and r.n_nonzero_minors == 392 and r.cols == [1, 3, 5, 7, 8, 10]
    assert arb(r.det).overlaps(arb("7872428136359854080", "1e-6"))

def test_trace_coords_roundtrip_bolza():
    from bolzakit.certified.chart import bolza_traces
    tau = bolza_traces()
    assert abs(tau[0] - 340.70562748477140586) < 1e-9 and abs(tau[1] - 134.88225099390856234) < 1e-9
    s = Surface.from_trace_coords(tau)
    assert s.is_bolza or max(abs(a - b) for a, b in zip(s.klein_coords(), Surface.bolza().klein_coords())) < 1e-20

def test_trace_coords_out_of_range():
    from bolzakit.certified.chart import bolza_traces
    tau = bolza_traces(); tau[0] *= 1.5
    with pytest.raises(OutOfValidityRange):
        Surface.from_trace_coords(tau)


# ------------------------------------------------------------------ additional checks (not in the brief)
def test_rank_certificate_shape():
    c = certify_rank().certificate
    assert c.kind == "rank" and c.certified and c.data["cols"] == [1, 3, 5, 7, 8, 10]
    assert len(c.data["J"]) == 6 and all(len(row) == 6 for row in c.data["J"])
    assert [d.id for d in c.dependencies] == ["D12", "D13"]


def test_trace_coords_perturbed_contains_tau():
    from bolzakit.certified.chart import bolza_traces
    tau = bolza_traces(); tau[0] += 0.5; tau[3] -= 0.25
    s = Surface.from_trace_coords(tau)
    assert not s.is_bolza and not s.r2_invariant and repr(tau[0]) in s.label
    from bolzakit.certified._arb import prec
    with prec(128):
        tv = traces(s._klein_arb())
    assert all(t.contains(arb(x)) for t, x in zip(tv, tau))
    # certify_rank on the deformed surface (balls) still certifies rank 6
    assert certify_rank(s).certified
