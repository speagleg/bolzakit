import math

import pytest
flint = pytest.importorskip("flint")
from bolzakit.certified.cluster import certify_cluster, Lam1Input, GapInput

# Published numbers (SRC cert_sign_{24_48,48_96,48_nosym}.json) were produced with LAM = arb("3.8388872588",
# "1.0000001e-6"); LEGACY is that interval, so goldens reproduce them to 1e-6 (controller ruling, ledger).
LEGACY = Lam1Input(lower="3.8388862587999999", upper="3.8388882588000001",
                   source="SU13 as used in item-23 scripts: arb('3.8388872588','1.0000001e-6')", dep_id="D1")
SU13 = Lam1Input(lower="3.83888625884", upper="3.83888825884", source="SU13 (D1)")
LB = Lam1Input("3.8388757010", "3.8388872830",
               "the decimal-rounded (outward) certified LG interval: Lehmann–Goerisch lower bound + P2 upper bound")
GAP = GapInput(lower="4.0", source="trace-formula LP count (certified)")

FB24 = -9.713754495605823          # cert_sign_24_48.json n=24
FB48 = -12.566084055230016         # cert_sign_48_96.json n=48
FB48_NOSYM = -12.566084126010532   # cert_sign_48_nosym.json n=48
FB96 = -12.850476214196528         # cert_sign_48_96.json n=96
# lb_cert_sign_96.json publishes -9.341154560446736 with arb("3.83888149202636215", "5.791e-6"), whose lower end
# 3.838875701026... lies below the certified Lehmann-Goerisch bound 3.8388757010527243... (item-23 LB_RESULTS.md), so
# that published value is sound.  LB here is the decimal-rounded (outward) certified LG interval, whose lower end
# 3.8388757010 is lower still; d g2 / d lam_low ~ 2e5 makes its bound 5.9e-6 weaker (same shift at n=24: -6.7971273
# vs SRC -6.7971333; SRC run-to-run spread there 8e-8).  Pin = this implementation (fixed ARPACK v0).
FB96_LB = -9.341148662380872


@pytest.fixture(scope="module")
def fb24():
    return certify_cluster("F_B", degree=2, n=24, lam1=LEGACY, gap=GAP, second_order=True)


def test_fb_n24_golden(fb24):
    r = fb24
    assert r.certified and abs(r.kyfan_upper_inprocess - FB24) < 1e-6
    assert r.b_ee_upper is None and len(r.A_float) == 3


def test_fb_n24_certificate(fb24):
    c = fb24.certificates[0]
    assert c.kind == "kyfan_second_order" and c.certified
    assert c.inputs["direction"] == "F_B" and c.inputs["degree"] == 2 and c.inputs["n"] == 24
    assert c.inputs["sym"] is True and c.inputs["lam1"]["lower"] == LEGACY.lower
    assert c.claims[0].lhs == "g2_upper" and c.claims[0].op == "<" and c.claims[0].rhs == "0"
    d = c.data
    assert d["assembly"] == "g2_upper_v1"
    assert set(d["terms"]) == {"trNX", "A_term", "c1_term"}
    t = {k: flint.arb(v) for k, v in d["terms"].items()}
    g2 = flint.arb(d["g2_enclosure"])
    # g2_upper_v1: the terms are additive, g2 = trNX + A_term + c1_term
    assert (t["trNX"] + t["A_term"] + t["c1_term"]).overlaps(g2)
    assert fb24.kyfan_upper_inprocess <= float(g2.upper())
    assert "g2_central" in d["diagnostics"]
    # one-at-a-time error rows are diagnostics only (non-additive)
    diag = d["diagnostics"]
    assert diag["non_additive"] is True
    assert set(diag["error_breakdown_upper_minus_central"]) == {"A", "M0", "M1", "X1", "X2", "c1"}
    assert "error_breakdown_upper_minus_central" not in d["terms"]
    ids = {x.id for x in c.dependencies}
    srcs = {x.source for x in c.dependencies}
    assert {"D1", "gap", "D4", "D9", "F1", "Lemma4"} == ids and LEGACY.source in srcs and GAP.source in srcs
    assert "F1_PROOF" not in "".join(srcs) and "D10" not in ids     # PROOF.md §Dependencies: Ky Fan is D9


def test_fb_n24_su13_exact_certifies():
    r = certify_cluster("F_B", degree=2, n=24, lam1=SU13, gap=GAP, second_order=True)
    assert r.certified and abs(r.kyfan_upper_inprocess - FB24) < 1e-4


@pytest.mark.slow
def test_fb_n48_golden_and_nosym():
    r = certify_cluster("F_B", degree=2, n=48, lam1=LEGACY, gap=GAP, second_order=True)
    r0 = certify_cluster("F_B", degree=2, n=48, lam1=LEGACY, gap=GAP, second_order=True, sym=False)
    assert r.certified and r0.certified
    assert abs(r.kyfan_upper_inprocess - FB48) < 1e-6 and abs(r0.kyfan_upper_inprocess - FB48_NOSYM) < 1e-6
    assert abs(r0.kyfan_upper_inprocess - r.kyfan_upper_inprocess) < 1e-6


@pytest.mark.slow
def test_fb_n96_golden():
    r = certify_cluster("F_B", degree=2, n=96, lam1=LEGACY, gap=GAP, second_order=True)
    assert r.certified and abs(r.kyfan_upper_inprocess - FB96) < 1e-6


@pytest.mark.slow
def test_fb_n96_lehmann_goerisch_lam1_certifies():
    r = certify_cluster("F_B", degree=2, n=96, lam1=LB, gap=GAP, second_order=True)
    assert r.certified and r.kyfan_upper_inprocess <= -9.3411 and abs(r.kyfan_upper_inprocess - FB96_LB) < 1e-6


def _rebuild_g2(cert):
    """Independent g2_upper_v1 checker, written only from the _kyfan docstring: g2 from cert.data["raw"] alone."""
    from flint import arb, arb_mat
    from bolzakit.certified._arb import prec
    with prec(96):
        R = cert.data["raw"]
        PM = arb(0, 1)
        lam = arb(R["lam"][0]).union(arb(R["lam"][1])); lam4 = arb(R["lam4"]); lamU = arb(float(R["lamU"]))
        assert lamU >= lam.upper()
        G = {k: [[arb(x) for x in row] for row in M] for k, M in R["G"].items()}
        sqrt_up = lambda x: arb(x.upper()).sqrt()
        W = {k: [sqrt_up(arb(x)) for x in v] for k, v in R["Wdiag"].items()}
        W["mps"] = [sqrt_up(G["M0"][f][f]) for f in range(7)]
        b1, b1r, b2, b2r = (arb(x) for x in R["bmax"])
        exc = [G["K0"][i][i] - lam * G["M0"][i][i] for i in range(3)]
        sm = [sqrt_up(e / (lam4 - lam)) for e in exc]
        sq = [sqrt_up(e * lam4 / (lam4 - lam)) for e in exc]
        B1 = [[G["K1"][a][b] - lam * G["M1"][a][b] for b in range(7)] for a in range(7)]
        B2 = [[2 * (G["K2"][a][b] - lam * G["M2"][a][b]) for b in range(7)] for a in range(7)]
        S = [[G["K0"][a][b] - lam * G["M0"][a][b] for b in range(7)] for a in range(7)]
        W1q, W1m, W2q, W2m, mps = (W[k] for k in ("W1q", "W1m", "W2q", "W2m", "mps"))
        M0e = arb_mat(3, 3); Ae = arb_mat(3, 3); Xe = arb_mat(3, 3); M1e = arb_mat(3, 3); c1e = arb_mat(3, 1)
        ONE = 6
        for i in range(3):
            p = 3 + i
            for j in range(3):
                P = 3 + j
                M0e[i, j] = G["M0"][i][j] + sm[i] * sm[j] * PM
                Ae[i, j] = B1[i][j] + (W1q[i] * sq[j] + lamU * W1m[i] * sm[j] + W1q[j] * sq[i] + lamU * W1m[j] * sm[i]
                                       + b1 * sq[i] * sq[j] + lamU * b1r * sm[i] * sm[j]) * PM
                Xe[i, j] = (B2[i][j] + 2 * B1[p][j] + 2 * B1[i][P] + 2 * S[p][P]
                            + (W2q[i] * sq[j] + lamU * W2m[i] * sm[j] + W2q[j] * sq[i] + lamU * W2m[j] * sm[i]
                               + b2 * sq[i] * sq[j] + lamU * b2r * sm[i] * sm[j]) * PM
                            + (2 * (W1q[p] * sq[j] + lamU * W1m[p] * sm[j])
                               + 2 * (W1q[P] * sq[i] + lamU * W1m[P] * sm[i])) * PM)
                M1e[i, j] = (G["M1"][i][j] + G["M0"][p][j] + G["M0"][i][P]
                             + (W1m[i] * sm[j] + W1m[j] * sm[i] + b1r * sm[i] * sm[j] + mps[p] * sm[j]
                                + mps[P] * sm[i]) * PM)
            c1e[i, 0] = G["M1"][i][ONE] + G["M0"][p][ONE] + W1m[ONE] * sm[i] * PM
        N = M0e.inv()
        m11 = G["M0"][ONE][ONE]
        return ((N * Xe).trace() - 2 * (N * M1e * N * Ae).trace()
                + 2 * lam * (c1e.transpose() * N * c1e)[0, 0] / m11)


def test_g2_rebuilds_from_raw_data(fb24):
    c = fb24.certificates[0]
    g2 = _rebuild_g2(c)
    assert g2.overlaps(flint.arb(c.data["g2_enclosure"]))
    assert g2.upper() < 0                      # the claim also holds from the stored raw inputs alone
    assert abs(float(g2.upper()) - fb24.kyfan_upper) < 1e-3


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
    assert len(c.inputs["tangent_projected"]) == 16 and float(c.inputs["projection_residual"]) < 1e-12
    g_named = flint.arb(fb24.certificates[0].data["g2_enclosure"])
    g_vec = flint.arb(c.data["g2_enclosure"])
    assert g_named.overlaps(g_vec)
    assert abs(r.kyfan_upper - fb24.kyfan_upper) < float(g_named.rad() + g_vec.rad())
    # the central values (transfer errors set to 0) agree far more tightly
    c_named = flint.arb(fb24.certificates[0].data["diagnostics"]["g2_central"])
    c_vec = flint.arb(c.data["diagnostics"]["g2_central"])
    assert c_named.overlaps(c_vec)
    # NOT certified at n=24: the min-norm second-order Klein jet (surface.vertex_jets for vectors) has ~5x the
    # vertex curvature of the F_B family, so beta2 and the X2 transfer error are much larger (measured: +3.47)
    assert not r.certified


def _stored_upper(ball):
    """Exact upper end of a stored '[mid +/- rad]' ball string."""
    from fractions import Fraction
    mid, rad = ball.strip()[1:-1].split("+/-")
    return Fraction(mid.strip()) + Fraction(rad.strip())


def test_kyfan_upper_is_certificate_attested(fb24):
    from fractions import Fraction
    stored = _stored_upper(fb24.certificates[0].data["g2_enclosure"])
    assert Fraction(fb24.kyfan_upper) >= stored
    assert Fraction(fb24.kyfan_upper) >= Fraction(fb24.kyfan_upper_inprocess)
    assert fb24.kyfan_upper - float(stored) < 1e-12 * max(1.0, abs(float(stored)))    # rounded up, not loosened
