import math

import pytest
flint = pytest.importorskip("flint")
from bolzakit.certified.cluster import certify_cluster, Lam1Input, GapInput

# Published numbers (SRC cert_first_order_{48,96_192}.json) were produced with LAM = arb("3.8388872588",
# "1.0000001e-6"); LEGACY is that interval, so the goldens reproduce them to 1e-6 (controller ruling, ledger).
LEGACY = Lam1Input(lower="3.8388862587999999", upper="3.8388882588000001",
                   source="SU13 as used in item-23 scripts: arb('3.8388872588','1.0000001e-6')")
SU13 = Lam1Input(lower="3.83888625884", upper="3.83888825884", source="SU13 (D1)")
LB = Lam1Input("3.8388757010", "3.8388872830", "Lehmann–Goerisch lower bound + P2 upper bound")
GAP = GapInput(lower="4.0", source="trace-formula LP count (certified)")


def test_inputs_required():
    with pytest.raises(ValueError, match="required"):
        certify_cluster("F_A", degree=1, n=24)
    with pytest.raises(ValueError, match="required"):
        certify_cluster("F_A", degree=1, n=24, lam1=SU13)


def test_gap_must_exceed_lam1():
    with pytest.raises(ValueError, match="gap"):
        certify_cluster("F_A", degree=1, n=12, lam1=SU13, gap=GapInput(lower="3.8", source="bad"))
    with pytest.raises(ValueError, match="gap"):
        certify_cluster("F_A", degree=1, n=12, lam1=SU13, gap=GapInput(lower="3.83888825884", source="bad"))


def test_other_directions_not_published_at_degree1():
    with pytest.raises(NotImplementedError, match="published for F_A"):
        certify_cluster("F_B", degree=1, n=12, lam1=SU13, gap=GAP)


def test_n24_not_certified_and_certificate():
    r = certify_cluster("F_A", degree=1, n=24, lam1=LEGACY, gap=GAP)
    assert not r.A_nonzero and not r.certified
    # regression pin (this implementation, fixed ARPACK v0). No published value exists for the current script at
    # n=24: cert_first_order_24_48.json (93.418...) predates the energy-weighted transfer bound (no beta_rms keys).
    assert abs(r.b_ee_upper_inprocess - 21.700387418284663) < 1e-6
    assert abs(abs(r.A_float[0]) - 4.815547668200662) < 1e-8
    c = r.certificates[0]
    assert c.kind == "cluster_first_order" and not c.certified
    assert c.data["assembly"] == "b_ee_upper_v1"
    assert set(c.data["terms"]) == {"b_uu_over_m", "transfer_err_over_m"}
    assert "beta_a" in c.data["diagnostics"] and "beta_rho" in c.data["diagnostics"]
    assert c.claims[0].lhs == "b_ee_upper_over_m" and c.claims[0].op == "<" and c.claims[0].rhs == "0"
    srcs = [d.source for d in c.dependencies]
    assert LEGACY.source in srcs and GAP.source in srcs and any(d.id == "F1" for d in c.dependencies)
    # b_ee_upper_v1: the claim ball is the sum of the two terms
    s = flint.arb(c.data["terms"]["b_uu_over_m"]) + flint.arb(c.data["terms"]["transfer_err_over_m"])
    assert s.overlaps(flint.arb(c.data["b_ee_upper_over_m"]))
    # the stored (short) ball string is a wider enclosure of the certified ball
    assert float(flint.arb(c.data["b_ee_upper_over_m"]).upper()) >= r.b_ee_upper_inprocess


@pytest.mark.slow  # ~90 s on one core
def test_n48_not_certified_golden():
    # SRC cert_first_order_48.json (the current script; cert_first_order_24_48.json's 35.339 is from an older
    # version without the energy-weighted Sa/Sr transfer bound)
    r = certify_cluster("F_A", degree=1, n=48, lam1=LEGACY, gap=GAP)
    assert not r.A_nonzero and abs(r.b_ee_upper_inprocess - 4.268784508109093) < 1e-6
    assert abs(abs(r.A_float[0]) - 4.808400696505938) < 1e-8


def test_nan_in_element_data_never_passes(monkeypatch):
    import bolzakit.certified.cluster as C
    real = C.element_dual
    calls = [0]

    def poisoned(*a):
        Ke, Me, dKe, dMe, ba, br = real(*a)
        calls[0] += 1
        if calls[0] == 5:
            dKe = [row[:] for row in dKe]
            dKe[0][1] = dKe[1][0] = flint.arb("nan")
        return Ke, Me, dKe, dMe, ba, br
    monkeypatch.setattr(C, "element_dual", poisoned)
    r = certify_cluster("F_A", degree=1, n=12, lam1=SU13, gap=GAP)
    assert not r.A_nonzero and not r.certified and not r.certificates[0].certified


def test_nan_from_straddling_sqrt_never_passes():
    # lam1 ball reaching above mu makes mu - lam1 straddle 0, so sqrt(q_w) is NaN: never a pass
    wide = Lam1Input(lower="3.5", upper="3.99", source="test")
    r = certify_cluster("F_A", degree=1, n=12, lam1=wide, gap=GapInput("4.5", "test"))
    assert not r.A_nonzero and (r.b_ee_upper is None or math.isnan(r.b_ee_upper) or r.b_ee_upper >= 0)


@pytest.mark.slow
def test_n96_certified_golden():
    r = certify_cluster("F_A", degree=1, n=96, lam1=LEGACY, gap=GAP)
    assert r.A_nonzero and r.certified and abs(r.b_ee_upper_inprocess - (-1.0282409843057394)) < 1e-6


@pytest.mark.slow
def test_n96_su13_exact_certifies():
    r = certify_cluster("F_A", degree=1, n=96, lam1=SU13, gap=GAP)
    assert r.A_nonzero and abs(r.b_ee_upper_inprocess - (-1.0282409843057394)) < 1e-4


@pytest.mark.slow
def test_n192_certified_golden():
    r = certify_cluster("F_A", degree=1, n=192, lam1=LEGACY, gap=GAP)
    assert r.A_nonzero and abs(r.b_ee_upper_inprocess - (-3.0875694043934345)) < 1e-6


@pytest.mark.slow
def test_n192_lehmann_goerisch_lam1_certifies():
    # published lb_cert_first_order (LB_RESULTS.md): b(e,e)/m <= -3.0383476242423058
    r = certify_cluster("F_A", degree=1, n=192, lam1=LB, gap=GAP)
    assert r.A_nonzero and r.b_ee_upper_inprocess <= -3.0383 and abs(r.b_ee_upper_inprocess - (-3.0383476242423058)) < 1e-6


@pytest.mark.slow
def test_weak_lam1_lower_does_not_certify():
    # PROOF.md Lemma 5's robustness remark ("survives any lambda_1 lower bound >= 3.83788") does NOT hold for this
    # P1 n=192 test: measured b_ee_upper/m = +0.01027 (mu - lam1 grows to 1.2e-3, so q(w)/m ~ 0.03). The plan's
    # test asserted it certifies; reported to the controller instead of loosening anything.
    weak = Lam1Input(lower="3.83788", upper="3.83888825884", source="weak")
    r = certify_cluster("F_A", degree=1, n=192, lam1=weak, gap=GAP)
    assert not r.A_nonzero and abs(r.b_ee_upper_inprocess - 0.0102731604129076) < 1e-4


def test_b_ee_upper_is_certificate_attested():
    from fractions import Fraction
    r = certify_cluster("F_A", degree=1, n=24, lam1=LEGACY, gap=GAP)    # finite (not certified) at n=24
    mid, rad = r.certificates[0].data["b_ee_upper_over_m"].strip()[1:-1].split("+/-")
    stored = Fraction(mid.strip()) + Fraction(rad.strip())
    assert Fraction(r.b_ee_upper) >= stored
    assert Fraction(r.b_ee_upper) >= Fraction(r.b_ee_upper_inprocess)
    assert r.b_ee_upper - float(stored) < 1e-12 * max(1.0, abs(float(stored)))
