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


@pytest.mark.slow                            # 15 G evaluations, ~27 s
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
    beta1 = -0.0945730196647619
    pred = 2 * (1 / (4 * math.pi) + beta1) * 0.01
    assert abs((g.robin_mass(0.1) - g.robin_mass(0.0)) - pred) < 1e-4 * abs(pred) + 1e-7


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


def test_green_out_of_range(g):
    zfar = 0.5 * Geo.VERTICES[0] + 0.5 * Geo.MIDPOINTS[0]     # between charts
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
