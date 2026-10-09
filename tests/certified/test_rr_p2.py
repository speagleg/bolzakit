import pytest
flint = pytest.importorskip("flint")
from bolzakit.certified import Surface, certify_spectrum


# lambda_1 is degenerate (x1,x2 are an arbitrary rotation in the eigenspace), so Arb radii depend on the float
# eigenvector. n=24 over 6 start-vector seeds: 3.838893607026..3.838893607287 (spread 2.6e-10; the SRC value
# 3.838893606900 lies 1.3e-10 below the lowest); n=48 is identical to 1e-15 on all seeds. Hence 1e-9 at n=24.
def test_p2_n24_golden():
    r = certify_spectrum(Surface.bolza(), k=3, degree=2, n=24)
    assert r.certified and abs(r.upper[-1] - 3.838893606899935) < 1e-9
    assert abs(r.cluster_sum_upper - 11.516677957238421) < 1e-8


def test_p2_n48_golden():
    r = certify_spectrum(Surface.bolza(), k=3, degree=2, n=48)
    assert abs(r.upper[-1] - 3.8388876433743824) < 1e-11


@pytest.mark.slow
def test_p2_n96_golden():
    r = certify_spectrum(Surface.bolza(), k=3, degree=2, n=96)
    assert abs(r.upper[-1] - 3.8388872829261946) < 1e-12
    assert r.upper[-1] >= 3.83888725884                         # SU13 midpoint: an upper bound must not undercut it by more than SU13's radius
