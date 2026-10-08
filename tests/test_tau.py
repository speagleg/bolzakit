import math
from fractions import Fraction

import mpmath as mp
import numpy as np
import pytest

from bolzakit.spectral import tau as TA
from bolzakit.exact import group as G

T = TA


def test_closed_walk_counts():
    X = np.array([G.row(g) for g in G.GENS])
    m, _, _ = T.moments(X, np.ones(8, dtype=np.int64), 5, integer=True)
    assert m == [1, 0, 8, 0, 120, 0, 2192, 0, 44264, 0, 950608]
    # 8-regular tree: 1, 8, 120, 2192, 44248 ; difference 16 = 8-cycles through e (8 faces x 2 orientations)
    assert m[8] - 44248 == 16


def test_float_moments_match_integer():
    X, c, d = T.systole_shell()
    mf, _, _ = T.moments(X, c, 4)
    mi, _, _ = T.moments(X, np.ones(8, dtype=np.int64), 4, integer=True)
    for k in range(len(mi)):
        assert abs(mf[k] - c[0] ** k * mi[k]) <= 1e-13 * c[0] ** k * mi[k] + 1e-300


def test_moments_self_adjoint_general():
    X, c, d = T.build_h(4.3)
    assert len(c) == 24
    keys = {tuple(r) for r in X.tolist()}
    for r in X:
        assert tuple(G.row(G.inv(G.elt(r))).tolist()) in keys
    m, sizes, _ = T.moments(X, c, 3)
    assert m[1] == 0.0 and m[3] > 0       # tau(h) = 0; tau(h^3) > 0 (even elements present)
    # Hankel positivity
    H = np.array([[m[i + j] for j in range(4)] for i in range(4)])
    assert np.all(np.linalg.eigvalsh(H) > 0)


@pytest.mark.slow
def test_christoffel_and_certificates():
    m = [1, 0, 8, 0, 120, 0, 2192, 0, 44264]
    lam = T.christoffel_bound(m)
    assert lam == Fraction(1, 1) / (1 / T.christoffel_bound(m))   # exact Fraction
    assert 0 < lam < 1
    # mu = (delta_{-1} + delta_{1})/2 has no atom at 0: lambda_2(0) = 1/(1 + 0^2/1) ... exact value 0 for n=2
    # (P(x) = 1 - x^2 vanishes on the support) -- needs n = 3 (moments to 4)
    mm = [1 if k % 2 == 0 else 0 for k in range(5)]
    assert T.christoffel_bound(mm[:3]) == 1          # degree-1 polys: best is P = 1
    # certificate sanity on a known measure: uniform on [-1,1] (mu((-1,0)) = 1/2)
    K = 12
    mu = [float(mp.quad(lambda x: mp.chebyt(k, x) / 2, [-1, 1])) for k in range(K + 1)]
    b = T.rigorous_bounds(mu, [1e-15] * (K + 1))
    assert b["lower_lt0"] <= 0.5 <= b["upper_le0"]
    assert b["upper_le0"] - b["lower_lt0"] < 0.4
    # asymmetric known measure: 0.3 delta_{-0.5} + 0.7 delta_{0.6}
    mu2 = [0.3 * math.cos(k * math.acos(-0.5)) + 0.7 * math.cos(k * math.acos(0.6)) for k in range(K + 1)]
    b2 = T.rigorous_bounds(mu2, [1e-15] * (K + 1))
    assert b2["lower_lt0"] <= 0.3 + 1e-9 and 0.3 - 1e-9 <= b2["upper_le0"]


def test_tau_bounds_bundle():
    from bolzakit.spectral.tau import tau_bounds
    from bolzakit.provenance import info
    import pytest
    b = tau_bounds()                         # systole shell
    assert abs(b["lower"] - 0.40993338889993414) < 1e-12 and b["upper"] == 0.5
    assert info(tau_bounds).label == "CERTIFIED"
    for i in range(4):
        g = tau_bounds(i)
        assert 0 < g["lower"] <= g["upper"] < 1
    with pytest.raises(ValueError, match="systole_shell"):
        tau_bounds("nope")
