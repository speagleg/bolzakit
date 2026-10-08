import cmath, math
import numpy as np
import pytest
from bolzakit import green
from bolzakit.errors import OutOfValidityRange
from bolzakit.exact import group as G
from bolzakit.geometry import octagon as Geo


@pytest.fixture(scope="module")
def g():
    return green.green()          # T=2.6, D=15 (setup may take ~1 min; cached)


def test_laplacian_off_diagonal(g):          # -Δ_z G = -1/A away from w, i.e. Δ_z G = 1/A
    w = 0.1 + 0.05j
    for z in (0.3 + 0.1j, -0.2 + 0.25j, 0.05 - 0.3j):
        h = 1e-3
        lap = (g.G(z + h, w) + g.G(z - h, w) + g.G(z + 1j * h, w) + g.G(z - 1j * h, w) - 4 * g.G(z, w)) / h ** 2
        lap *= (1 - abs(z) ** 2) ** 2 / 4
        assert abs(lap - 1 / (4 * math.pi)) < 5e-5


def test_symmetry_and_c8(g):
    z, w = 0.2 + 0.1j, 0.1 + 0.05j
    assert abs(g.G(z, w) - g.G(w, z)) < 1e-14
    r = cmath.exp(1j * math.pi / 4)
    assert abs(g.G(r * z, r * w) - g.G(z, w)) < 1e-13


def test_robin_mass_expansion(g):
    # The residual (~3.85e-5) is dominated by the omitted spectral tail (~3.1e-5): it is deterministic, so the 5e-5
    # tolerance must not be loosened.
    # m(r) - m(0) = c2 r^2 + O(r^4) along the real axis; Richardson on f(r) = (m(r) - m(0))/r^2 removes the r^2 term of f
    beta1 = -0.0945730196647619
    m0 = g.robin_mass(0.0)
    f = lambda r: (g.robin_mass(r) - m0) / r ** 2
    c2 = (4 * f(0.05) - f(0.1)) / 3
    assert abs(c2 - 2 * (1 / (4 * math.pi) + beta1)) < 5e-5


def test_green_gamma_invariance(g):
    z, w = 0.2 + 0.1j, 0.1 + 0.05j
    gam = G.GENS[3]
    zt = G.act(gam, z)               # a Γ-translate far outside the octagon
    assert abs(g.G(zt, w) - g.G(z, w)) < 1e-10
    assert abs(g.robin_mass(zt) - g.robin_mass(z)) < 1e-10


def test_green_other_weierstrass_chart(g):
    # G is Aut-invariant, so the same configuration near another Weierstrass point gives the same value
    from bolzakit.exact.automorphisms import automorphisms
    A = automorphisms()
    p = [q for q in A.weierstrass_points() if abs(q) > 1e-6][0]
    M = np.linalg.inv(A.chart_to(p))
    z, w = 0.2 + 0.1j, 0.1 + 0.05j
    assert abs(g.G(Geo.mobius(M, z), Geo.mobius(M, w)) - g.G(z, w)) < 1e-10


def test_green_pair_without_common_chart(g):
    """Every single point is valid (covering radius 0.4056 < 0.45), so the only failure is a pair with no common
    chart: zfar is valid in chart 1 (|z_local| = 0.315, also 0.425 in chart 5) but shares no chart with w = 0."""
    zfar = 0.5 * Geo.VERTICES[0] + 0.5 * Geo.MIDPOINTS[0]
    with pytest.raises(OutOfValidityRange, match="0.45"):
        g.G(zfar, 0.0)


def test_green_coincident_points(g):
    with pytest.raises(ValueError, match="robin_mass"):
        g.G(0.1 + 0.1j, 0.1 + 0.1j)


def test_green_overlapping_charts(g):
    # The 0.45-discs about adjacent Weierstrass points overlap (distance 1.5286 < 2*2 atanh 0.45).  z is nearer the
    # centre, but both z and w lie within 0.45 of the edge-midpoint Weierstrass point at -MIDPOINT_RADIUS: valid.
    z = -math.tanh(0.70 / 2)            # d(0, z) = 0.70, d(z, p) = 0.83
    w = -math.tanh(1.20 / 2)            # d(0, w) = 1.20 (> 2 atanh 0.45), d(w, p) = 0.33
    assert g.chart(z)[0] == 0 and abs(g.chart(w)[1]) < 0.45 and abs(w) > green.VALIDITY_RADIUS
    assert abs(g.G(z, w) - g.G(w, z)) < 1e-14

    # Same-point-in-two-charts consistency: this pair is within 0.45 of the same point in >= 2 charts (0 and 3); the
    # core G must agree between charts (the image sum runs over all of Gamma), here to 4e-11.
    z2, w2 = 0.07191337862080915 + 0.3738994124692343j, -0.1551718461079472 - 0.36792601406121417j
    zs, ws = g._locals(z2), g._locals(w2)
    ok = [i for i in range(len(zs)) if abs(zs[i]) < 0.45 and abs(ws[i]) < 0.45]
    assert len(ok) >= 2
    vals = [g._core.G(zs[i], ws[i]) for i in ok]
    assert max(vals) - min(vals) < 1e-9


# ---- final review C1: (T, D_ball) outside the validated set raise before anything is built

@pytest.mark.parametrize("kw", [dict(T=10.0), dict(D_ball=9.0), dict(T=1.0), dict(T=5.0), dict(T=2.0, D_ball=12.0),
                                dict(T=2.6, D_ball=13.5)])
def test_green_unvalidated_T_D_raises_fast(kw):
    import time
    t0 = time.time()
    with pytest.raises(OutOfValidityRange, match=r"\(2\.6, 15\.0\), \(2\.6, 13\.0\).*25%"):
        green.green(**kw)
    with pytest.raises(OutOfValidityRange):
        green.SurfaceGreen(**kw)
    assert time.time() - t0 < 1.0            # validated before the image ball is built


def test_green_validated_pairs_accepted(monkeypatch):
    import bolzakit.green.surface as S
    monkeypatch.setattr(S, "_green", lambda T, D: (T, D))
    assert green.green() == (2.6, 15.0)
    assert green.green(T=2.6 + 1e-13, D_ball=13) == (2.6, 13.0)
    assert set(green.VALIDATED_T_D) == {(2.6, 15.0), (2.6, 13.0)}
