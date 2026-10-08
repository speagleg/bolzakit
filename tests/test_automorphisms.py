from fractions import Fraction

import numpy as np
import pytest

import bolzakit.exact.automorphisms as AU
from bolzakit.exact import automorphisms as Au, group as G
from bolzakit.geometry import octagon as Geo


@pytest.fixture(scope="module")
def autq():
    return Au.quotient_group()


def test_normaliser_exact():
    r = Au.normalises_exactly()
    assert all(v["in_Gamma"] and v["dets_one"] for v in r.values())


def test_aut_group_gl23(autq):
    E, Tb = autq
    assert len(E) == 48
    o = Au.element_orders(Tb)
    assert {k: o.count(k) for k in set(o)} == {1: 1, 2: 13, 3: 8, 4: 6, 6: 8, 8: 12}
    assert len(Au.center(Tb)) == 2
    idx = {Au._key(M): i for i, M in enumerate(E)}
    gi = [idx[Au._key(Au._normalise(Au.M_RHO))], idx[Au._key(Au._normalise(Au.M_SIGMA))]]
    assert Au.find_isomorphism_to_gl23(Tb, gi) is not None


def test_weierstrass(autq):
    E, Tb = autq
    W = Au.weierstrass_points_from_involution()
    assert len(W) == 6
    radii = sorted(round(abs(z), 6) for z in W)
    assert radii == [0.0] + [round(G.MIDPOINT_RADIUS, 6)] * 4 + [round(G.VERTEX_RADIUS, 6)]
    for z in W:
        st = Au.stabiliser(E, Tb, z)
        assert len(st) == 8 and Au.is_cyclic(Tb, st)
    assert len(Au.orbit_on_surface(E, 0j)) == 6
    # signature (0; 2, 3, 8): orbit sizes 24, 16, 6 ; Riemann-Hurwitz 2g-2 = 48(-2 + 1/2 + 2/3 + 7/8)
    assert len(Au.orbit_on_surface(E, Au.Q_ORDER2 + 0j)) == 24
    assert len(Au.orbit_on_surface(E, Au.P_ORDER3)) == 16
    assert 48 * (-2 + Fraction(1, 2) + Fraction(2, 3) + Fraction(7, 8)) == 2


@pytest.fixture(scope="module")
def A():
    return AU.automorphisms()


def test_orders(A):
    assert A.order == 48 and A.full_order == 96
    assert A.element_order_counts() == {1: 1, 2: 13, 3: 8, 4: 6, 6: 8, 8: 12}   # GL(2,3)


def test_isomorphism_is_homomorphism(A):
    phi = A.isomorphism_to_GL23()
    def m3(x, y):
        return tuple(tuple(sum(x[i][k] * y[k][j] for k in range(2)) % 3 for j in range(2)) for i in range(2))
    for i in range(48):
        for j in range(48):
            assert phi[int(A.table[i, j])] == m3(phi[i], phi[j])
    assert len(set(phi.values())) == 48


def test_six_weierstrass_points_one_orbit(A):
    W = A.weierstrass_points()
    assert len(W) == 6
    orbit = AU.orbit_on_surface(A.elements, 0j)
    assert len(orbit) == 6          # Aut+ acts transitively on the Weierstrass points
    # the orbit of 0 equals the Weierstrass set modulo Gamma, in both directions
    assert all(any(AU.same_point_mod_gamma(o, w) for w in W) for o in orbit)
    assert all(any(AU.same_point_mod_gamma(w, o) for o in orbit) for w in W)


def test_chart_to_sends_point_to_origin(A):
    for p in A.weierstrass_points():
        M = A.chart_to(p)
        assert abs(Geo.mobius(M, p)) < 1e-9
