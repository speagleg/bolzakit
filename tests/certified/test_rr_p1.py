import pytest
flint = pytest.importorskip("flint")
from bolzakit.certified import Surface, certify_spectrum

# NOTE: lambda_1 has multiplicity 2 on the FE level, so the float eigenvectors x1,x2 are an arbitrary rotation inside the
# eigenspace; the certified bound depends on that rotation through the Arb radii (n=12 spread: 3.9797..3.9830 over
# seeds; the SRC json itself is not reproducible to 1e-9). certify_spectrum fixes the ARPACK start vector, so results
# are deterministic, and the goldens are checked to the spread, not to 1e-9.
def test_p1_n24_golden():
    r = certify_spectrum(Surface.bolza(), k=3, degree=1, n=24)
    assert r.certified
    assert abs(r.upper[-1] - 3.857985420933427) < 1e-5          # SRC cert_rr_bolza_12_24.json; see note below
    assert r.upper[-1] >= 3.8388872588                          # contains true lambda_1 from above
    assert r.certificate.kind == "rr_upper" and len(r.certificate.data["K"]) == 4

def test_p1_n12_golden():
    r = certify_spectrum(Surface.bolza(), k=3, degree=1, n=12)
    assert abs(r.upper[-1] - 3.9827723626133213) < 5e-3

@pytest.mark.slow
def test_p1_n192_golden():
    r = certify_spectrum(Surface.bolza(), k=3, degree=1, n=192)
    assert abs(r.upper[-1] - 3.839184804702046) < 1e-9

def test_bisection_never_certifies_below_truth():
    from bolzakit.certified.rayleigh_ritz import is_pd
    from flint import arb
    # mu*M - K with mu just below the RR max eigenvalue must NOT be PD
    r = certify_spectrum(Surface.bolza(), k=3, degree=1, n=12)
    K = [[arb(s) for s in row] for row in r.certificate.data["K"]]
    M = [[arb(s) for s in row] for row in r.certificate.data["M"]]
    mu = arb(r.float_eigs[3]) - arb("1e-6")
    assert not is_pd([[mu * M[i][j] - K[i][j] for j in range(4)] for i in range(4)])
