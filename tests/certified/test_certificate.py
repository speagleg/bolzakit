import json, math, pytest
flint = pytest.importorskip("flint")
from flint import arb
from bolzakit.certified import _arb as A
from bolzakit.certified.certificate import Certificate, Claim, Dependency

def test_certainly_lt_rejects_nan_and_straddle():
    assert A.certainly_lt(arb(1), 2)
    assert not A.certainly_lt(arb(1, 2), 2)          # [−1, 3] straddles 2
    nan = arb("nan")
    assert not A.certainly_lt(nan, 2) and not A.certainly_gt(nan, 0)
    assert not A.certainly_gt(arb(0, 1).sqrt(), -1)  # sqrt of a ball straddling 0 is indeterminate

def test_ball_roundtrip_and_outward_rounding():
    x = arb("3.8388872830", "1e-12")
    assert A.parse_ball(A.ball_str(x)).overlaps(x)
    assert float(A.rational_upper(x)) >= 3.8388872830 + 1e-12 - 1e-18
    assert float(A.rational_lower(x)) <= 3.8388872830 - 1e-12 + 1e-18

def test_prec_context_restores():
    old = flint.ctx.prec
    with A.prec(200):
        assert flint.ctx.prec == 200
    assert flint.ctx.prec == old

def _cert():
    return Certificate(kind="rank", inputs={"surface": "bolza"}, claims=[Claim("det", "!=0", "0", True)],
                       data={"J": [["1"]], "cols": [0]}, dependencies=[Dependency("D12", "T2 is a 6-manifold", "FN")])

def test_roundtrip_and_sha(tmp_path):
    c = _cert(); p = tmp_path / "c.json"; c.save(p)
    d = Certificate.load(p)
    assert d.sha256() == c.sha256() and d.certified
    raw = json.loads(p.read_text()); assert raw["schema"] == "bolzakit-cert/1" and raw["sha256"] == c.sha256()

def test_not_certified_if_any_claim_fails_or_empty():
    c = _cert(); c.claims.append(Claim("x", "<", "1", False)); assert not c.certified
    c.claims = []; assert not c.certified
