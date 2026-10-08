import cmath
import math

import numpy as np
import pytest

from bolzakit import geometry
from bolzakit.exact import group as G
from bolzakit.geometry import octagon as Geo
from bolzakit.exact.automorphisms import same_point_mod_gamma


def test_octagon_geometry():
    assert abs(G.VERTEX_RADIUS - math.tanh(G.R_C / 2)) < 1e-14
    assert abs(G.MIDPOINT_RADIUS - math.tanh(G.R_IN / 2)) < 1e-14
    assert abs(math.cosh(G.R_C) - 1 / math.tan(math.pi / 8) ** 2) < 1e-12
    # side k (midpoint angle k pi/4) joins V_{k-1}, V_k; g_k maps side k+4 = [V_{k+3}, V_{k+4}] onto it
    V = Geo.VERTICES
    for k in range(8):
        a, b = G.act(G.GENS[k], V[(k + 3) % 8]), G.act(G.GENS[k], V[(k + 4) % 8])
        assert abs(a - V[k]) < 1e-9 and abs(b - V[(k - 1) % 8]) < 1e-9
        assert abs(G.act(G.GENS[k], Geo.MIDPOINTS[(k + 4) % 8]) - Geo.MIDPOINTS[k]) < 1e-9
    # all 8 vertices are identified (single vertex cycle, angle sum 8 * pi/4 = 2 pi)
    v0 = V[0]
    for v in V:
        assert same_point_mod_gamma(v0, v)


def test_cayley_roundtrip():
    for z in (0j, 0.3 + 0.2j, -0.7j, 0.99):
        assert abs(geometry.upper_to_disk(geometry.disk_to_upper(z)) - z) < 1e-12
    assert abs(geometry.disk_to_upper(0) - 1j) < 1e-15


def test_reduce_point_lands_in_F_and_preserves_orbit():
    rng = np.random.default_rng(1)
    for _ in range(50):
        r = 0.97 * math.sqrt(rng.random()); z = r * cmath.exp(2j * math.pi * rng.random())
        zr, word, M = geometry.reduce_point(z)
        assert geometry.in_F(zr)
        assert abs(geometry.mobius(M, z) - zr) < 1e-9


@pytest.mark.parametrize("z", [1.0, 1j, 0.6 + 0.8j, 2.0, -1.5j])
def test_reduce_point_rejects_outside_disk(z):
    # final review M4
    with pytest.raises(ValueError, match="open unit disk"):
        geometry.reduce_point(z)
