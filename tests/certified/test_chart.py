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
    # Replaces the plan's "is_bolza or within 1e-20" check: 1e-20 is unreachable, because the float tau_B are not the
    # exact (irrational) Bolza traces and their certified preimage lies ~6e-17 from Bolza.  So we assert closeness
    # 1e-15 and that the certified traces of the result contain tau.
    from bolzakit.certified.chart import bolza_traces, bolza_klein
    from bolzakit.certified._arb import prec
    tau = bolza_traces()
    assert abs(tau[0] - 340.70562748477140586) < 1e-9 and abs(tau[1] - 134.88225099390856234) < 1e-9
    s = Surface.from_trace_coords(tau)
    assert not s.is_bolza and not s.r2_invariant and s.trace_coords == tau
    with prec(128):
        xs = s._klein_arb()
        assert all(abs(float((a - b).mid())) < 1e-15 and a.rad() < 1e-20 for a, b in zip(xs, bolza_klein()))
        assert all(t.contains(arb(x)) for t, x in zip(traces(xs), tau))


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
    assert len(c.claims) == 4 and len(c.inputs["klein"]) == 16
    deps = {d.id: d.statement for d in c.dependencies}
    assert list(deps) == ["D12", "D13"]
    assert "6-manifold" in deps["D12"] and "tau_w are C^1" in deps["D12"]
    assert "mapping-class group action on T2 is continuous" in deps["D12"]
    assert deps["D13"] == "software: python-flint (Arb) and sympy are trusted"


def test_off_manifold_surface_not_certified():
    from bolzakit.certified.chart import bolza_klein
    from bolzakit.certified._arb import prec
    from bolzakit.certified.surface import _tohyp
    with prec(128):
        x = bolza_klein(); x[0] = x[0] * arb("1.001")   # breaks the side-length / angle-sum constraints
        s = Surface([_tohyp(x[2 * k], x[2 * k + 1]) for k in range(8)], label="off-manifold")
    r = certify_rank(s)
    assert r.certified is False and r.certificate.certified is False
    assert not next(c for c in r.certificate.claims if c.lhs == "constraints(x)").holds


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
    r = certify_rank(s)
    assert r.certified and r.certificate.certified and r.certificate.inputs["tau"] == tau
