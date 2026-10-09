import math

import pytest
flint = pytest.importorskip("flint")
from bolzakit.certified.cluster import certify_cluster, Lam1Input, GapInput

# Published numbers (SRC cert_sign_{24_48,48_96,48_nosym}.json) were produced with LAM = arb("3.8388872588",
# "1.0000001e-6"); LEGACY is that interval, so goldens reproduce them to 1e-6 (controller ruling, ledger).
LEGACY = Lam1Input(lower="3.8388862587999999", upper="3.8388882588000001",
                   source="SU13 as used in item-23 scripts: arb('3.8388872588','1.0000001e-6')")
SU13 = Lam1Input(lower="3.83888625884", upper="3.83888825884", source="SU13 (D1)")
LB = Lam1Input("3.8388757010", "3.8388872830", "Lehmann–Goerisch lower bound + P2 upper bound")
GAP = GapInput(lower="4.0", source="trace-formula LP count (certified)")

FB24 = -9.713754495605823          # cert_sign_24_48.json n=24
FB48 = -12.566084055230016         # cert_sign_48_96.json n=48
FB48_NOSYM = -12.566084126010532   # cert_sign_48_nosym.json n=48
FB96 = -12.850476214196528         # cert_sign_48_96.json n=96
# lb_cert_sign_96.json publishes -9.341154560446736, but its ball arb("3.83888149202636215", "5.791e-6") has lower end
# 3.83887570102634 > the certified Lehmann-Goerisch bound 3.8388757010 (by 2.6e-11; d g2 / d lam_low ~ 2e5), so it
# does not contain the certified interval's low end.  LB below is exactly the certified interval; its bound is
# 5.9e-6 weaker (measured; the same systematic shift appears at n=24: -6.7971273 vs SRC -6.7971333, SRC's own
# run-to-run spread there 8e-8).  Pin = this implementation (fixed ARPACK v0).
FB96_LB = -9.341148662380872


@pytest.fixture(scope="module")
def fb24():
    return certify_cluster("F_B", degree=2, n=24, lam1=LEGACY, gap=GAP, second_order=True)


def test_fb_n24_golden(fb24):
    r = fb24
    assert r.certified and abs(r.kyfan_upper - FB24) < 1e-6
    assert r.b_ee_upper is None and len(r.A_float) == 3


def test_fb_n24_certificate(fb24):
    c = fb24.certificates[0]
    assert c.kind == "kyfan_second_order" and c.certified
    assert c.inputs["direction"] == "F_B" and c.inputs["degree"] == 2 and c.inputs["n"] == 24
    assert c.inputs["sym"] is True and c.inputs["lam1"]["lower"] == LEGACY.lower
    assert c.claims[0].lhs == "g2_upper" and c.claims[0].op == "<" and c.claims[0].rhs == "0"
    d = c.data
    assert d["assembly"] == "g2_upper_v1"
    assert set(d["terms"]) == {"trNX", "A_term", "c1_term", "g2_central", "transfer_err"}
    t = {k: flint.arb(v) for k, v in d["terms"].items()}
    g2 = flint.arb(d["g2_enclosure"])
    # g2_upper_v1: g2 = trNX + A_term + c1_term, and g2_central + transfer_err encloses it
    assert (t["trNX"] + t["A_term"] + t["c1_term"]).overlaps(g2)
    # (ball_str strings are short until T10 widens them, so compare with the certified float, not string ends)
    s2 = t["g2_central"] + t["transfer_err"]
    assert s2.lower() <= fb24.kyfan_upper <= s2.upper() + 1e-3
    assert fb24.kyfan_upper <= float(g2.upper())
    # one-at-a-time error rows are diagnostics only (non-additive)
    diag = d["diagnostics"]
    assert diag["non_additive"] is True
    assert set(diag["error_breakdown_upper_minus_central"]) == {"A", "M0", "M1", "X1", "X2", "c1"}
    assert "error_breakdown_upper_minus_central" not in d["terms"]
    ids = {x.id for x in c.dependencies}
    srcs = {x.source for x in c.dependencies}
    assert {"D10", "F1", "Lemma4"} <= ids and LEGACY.source in srcs and GAP.source in srcs


def test_fb_n24_su13_exact_certifies():
    r = certify_cluster("F_B", degree=2, n=24, lam1=SU13, gap=GAP, second_order=True)
    assert r.certified and abs(r.kyfan_upper - FB24) < 1e-4


@pytest.mark.slow
def test_fb_n48_golden_and_nosym():
    r = certify_cluster("F_B", degree=2, n=48, lam1=LEGACY, gap=GAP, second_order=True)
    r0 = certify_cluster("F_B", degree=2, n=48, lam1=LEGACY, gap=GAP, second_order=True, sym=False)
    assert r.certified and r0.certified
    assert abs(r.kyfan_upper - FB48) < 1e-6 and abs(r0.kyfan_upper - FB48_NOSYM) < 1e-6
    assert abs(r0.kyfan_upper - r.kyfan_upper) < 1e-6


@pytest.mark.slow
def test_fb_n96_golden():
    r = certify_cluster("F_B", degree=2, n=96, lam1=LEGACY, gap=GAP, second_order=True)
    assert r.certified and abs(r.kyfan_upper - FB96) < 1e-6


@pytest.mark.slow
def test_fb_n96_lehmann_goerisch_lam1_certifies():
    r = certify_cluster("F_B", degree=2, n=96, lam1=LB, gap=GAP, second_order=True)
    assert r.certified and r.kyfan_upper <= -9.3411 and abs(r.kyfan_upper - FB96_LB) < 1e-6


def test_bad_inputs_never_certify():
    wide = Lam1Input(lower="3.8387872588", upper="3.8389872588", source="radius 1e-4")
    r = certify_cluster("F_B", degree=2, n=24, lam1=wide, gap=GAP, second_order=True)
    assert not r.certified and not r.certificates[0].certified
    with pytest.raises(ValueError):
        certify_cluster("F_B", degree=2, n=24, lam1=SU13, gap=GapInput("3.8", "bad"), second_order=True)
    excl = Lam1Input(lower="3.9", upper="3.91", source="excludes truth")
    r = certify_cluster("F_B", degree=2, n=12, lam1=excl, gap=GAP, second_order=True)
    assert not r.certified and not r.certificates[0].certified


def test_fa_degree2_rates_match_p1_limit():
    r = certify_cluster("F_A", degree=2, n=24, lam1=SU13, gap=GAP, second_order=False)
    rates = sorted(r.A_float)
    assert abs(rates[0] + 4.806) < 0.02 and abs(rates[1]) < 1e-6 and abs(rates[2] - 4.806) < 0.02
    assert r.kyfan_upper is None and not r.certified and r.certificates == []


def test_second_order_needs_degree2():
    with pytest.raises(NotImplementedError, match="degree=2"):
        certify_cluster("F_A", degree=1, n=12, lam1=SU13, gap=GAP, second_order=True)


def test_nan_in_element_data_never_passes(monkeypatch):
    import bolzakit.certified.forms as F
    real = F.element_p2
    calls = [0]

    def poisoned(*a, **k):
        Ks, Ms, bet = real(*a, **k)
        calls[0] += 1
        if calls[0] == 5:
            Ks = [[row[:] for row in Mx] for Mx in Ks]
            Ks[1][0][1] = Ks[1][1][0] = flint.arb("nan")
        return Ks, Ms, bet
    monkeypatch.setattr(F, "element_p2", poisoned)
    r = certify_cluster("F_B", degree=2, n=12, lam1=SU13, gap=GAP, second_order=True)
    assert not r.certified and r.kyfan_upper is None and not r.certificates[0].certified


def test_element_p2_nonfinite_beta_fails_closed():
    # NaN derivative seeds reach the relative-bound boxes; the former bare asserts must not raise or vanish
    # (python -O) but mark the element non-finite
    from bolzakit.certified import forms
    from bolzakit.certified._arb import prec
    from bolzakit.certified.jets import J
    from bolzakit.certified.mesh import Mesh, positions_jet
    from bolzakit.certified.surface import Surface
    with prec(96):
        mesh = Mesh(12)
        X = positions_jet(mesh, Surface.bolza().vertex_jets("F_B"))
        i0, i1, i2 = mesh.T[0]
        P0 = [J(c.v, flint.arb("nan"), c.e) for c in X[i0]]
        _, _, bet = forms.element_p2(P0, X[i1], X[i2], order=2)
    assert not all(b.is_finite() for b in bet)


def test_path_independence_named_vs_vector_fb(fb24):
    # dLambda = 0 at Bolza, so the bound on D2 Lambda(v, v) depends only on the tangent v
    from bolzakit.certified.surface import Surface
    v = Surface.bolza().klein_tangent("F_B")
    r = certify_cluster(v, degree=2, n=24, lam1=LEGACY, gap=GAP, second_order=True, sym=True)
    c = r.certificates[0]
    assert c.inputs["sym"] is False                         # forced for vector directions
    assert c.inputs["direction"] == [repr(float(x)) for x in v]
    g_named = flint.arb(fb24.certificates[0].data["g2_enclosure"])
    g_vec = flint.arb(c.data["g2_enclosure"])
    assert g_named.overlaps(g_vec)
    assert abs(r.kyfan_upper - fb24.kyfan_upper) < float(g_named.rad() + g_vec.rad())
    # the central values (transfer errors set to 0) agree far more tightly
    c_named = flint.arb(fb24.certificates[0].data["terms"]["g2_central"])
    c_vec = flint.arb(c.data["terms"]["g2_central"])
    assert c_named.overlaps(c_vec)
    # NOT certified at n=24: the min-norm second-order Klein jet (surface.vertex_jets for vectors) has ~5x the
    # vertex curvature of the F_B family, so beta2 and the X2 transfer error are much larger (measured: +3.47)
    assert not r.certified
