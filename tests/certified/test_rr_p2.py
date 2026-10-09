import pytest
flint = pytest.importorskip("flint")
from bolzakit.certified import Surface, certify_spectrum


# lambda_1 is degenerate (x1,x2 are an arbitrary rotation in the eigenspace), so Arb radii depend on the float
# eigenvector. n=24 over 6 start-vector seeds: 3.838893607026..3.838893607287 (spread 2.6e-10; the SRC value
# 3.838893606900 lies 1.3e-10 below the lowest); n=48 is identical to 1e-15 on all seeds. Hence 1e-9 at n=24.
def test_p2_n24_golden():
    r = certify_spectrum(Surface.bolza(), k=3, degree=2, n=24)
    assert r.certified and abs(r.upper_inprocess[-1] - 3.838893606899935) < 1e-9
    assert abs(r.cluster_sum_upper - 11.516677957238421) < 1e-8


def test_p2_n48_golden():
    r = certify_spectrum(Surface.bolza(), k=3, degree=2, n=48)
    assert r.certified
    assert abs(r.upper_inprocess[-1] - 3.8388876433743824) < 1e-11


@pytest.mark.slow
def test_p2_n96_golden():
    r = certify_spectrum(Surface.bolza(), k=3, degree=2, n=96)
    assert r.certified
    assert abs(r.upper_inprocess[-1] - 3.8388872829261946) < 1e-12
    assert r.upper[-1] >= 3.83888725884                         # SU13 midpoint: an upper bound must not undercut it by more than SU13's radius


def test_cluster_sum_upper_is_rigorous_and_claimed(monkeypatch):
    import flint
    from flint import arb
    import bolzakit.certified.rayleigh_ritz as rr
    orig, seen = rr._project_p2, {}

    def spy(mesh, surface, k):
        out = orig(mesh, surface, k)
        seen["GK"], seen["GM"] = out[2], out[3]
        return out
    monkeypatch.setattr(rr, "_project_p2", spy)
    r = certify_spectrum(Surface.bolza(), k=3, degree=2, n=12)
    flint.ctx.prec = 96
    try:
        kyfan = seen["GM"].solve(seen["GK"]).trace()
        rhs_ok = arb(r.certificate.claims[-1].rhs) >= kyfan.upper()
    finally:
        flint.ctx.prec = 53
    assert arb(r.cluster_sum_upper) >= kyfan.upper()      # certain comparison: float is an outward-rounded bound
    c = r.certificate.claims[-1]
    assert c.lhs == "lambda_1+lambda_2+lambda_3" and c.op == "<=" and c.holds
    assert rhs_ok


def test_cluster_sum_failed_solve_is_none(monkeypatch):
    from flint import arb_mat
    import bolzakit.certified.rayleigh_ritz as rr
    orig = rr._project_p2

    def singular_gm(mesh, surface, k):
        Kp, Mp, GK, GM, *rest = orig(mesh, surface, k)
        return (Kp, Mp, GK, arb_mat(k + 1, k + 1), *rest)   # zero matrix: solve cannot succeed
    monkeypatch.setattr(rr, "_project_p2", singular_gm)
    r = certify_spectrum(Surface.bolza(), k=3, degree=2, n=12)
    assert r.cluster_sum_upper is None
    assert r.certificate.claims[-1].holds is False


@pytest.mark.parametrize("degree", [1, 2])
def test_upper_is_certificate_attested(degree):
    from fractions import Fraction
    r = certify_spectrum(Surface.bolza(), k=3, degree=degree, n=12)
    assert r.certified and len(r.upper) == len(r.upper_inprocess) == 3
    assert all(Fraction(u) >= Fraction(t) for u, t in zip(r.upper, r.upper_inprocess))
    assert Fraction(r.upper[-1]) >= Fraction(r.certificate.data["mu"])
    assert r.certificate.claims[0].rhs == r.certificate.data["mu"]
