import pytest
flint = pytest.importorskip("flint")
from flint import arb
from bolzakit.certified import Surface, certify_spectrum
from bolzakit.errors import OutOfValidityRange

def test_fb_closed_form_at_zero_is_bolza():
    J = Surface.bolza().vertex_jets("F_B")
    C = J[0][0]
    assert (C.c[0] - (3 + 2 * arb(2).sqrt())).contains(0)

def test_deform_too_large_raises():
    with pytest.raises(OutOfValidityRange, match="0.05"):
        Surface.bolza().deform("F_B", 0.2)

def test_deformed_surface_spectrum_moves_down_along_FB():
    s = Surface.bolza().deform("F_B", 0.01)
    r = certify_spectrum(s, k=3, degree=2, n=24)
    assert r.certified
    # Lambda'' ~ -14.27 along F_B, so sum drops by ~ 7.1e-4 at t=0.01 (float check, generous)
    assert r.cluster_sum_upper < 11.516677957238421

def test_vector_direction_matches_named_FB_tangent():
    b = Surface.bolza()
    v = b.klein_tangent("F_B")                      # 16 floats: d/dt of Klein (x, y) of V0..V7
    Jv = b.vertex_jets(v)
    Jn = b.vertex_jets("F_B")
    for a, c in zip(Jv, Jn):
        for x, y in zip(a, c):
            assert abs(float(x.c[1].mid()) - float(y.c[1].mid())) < 1e-9


# ------------------------------------------------------------------ additional tests (beyond the brief)
from bolzakit.certified import surface as surf  # noqa: E402
from bolzakit.certified._arb import prec  # noqa: E402


def _residuals(s):
    with prec(96):
        return surf._constraints_hyp(s.vertices_hyp)


def test_deformed_surfaces_are_certified_closed():
    """deform() output contains an exact zero of the 5 gluing constraints (L_{k+4}=L_k, angle sum 2 pi)."""
    import numpy as np
    b = Surface.bolza()
    v = 0.05 * np.random.RandomState(1).randn(16)   # Klein steps ~1e-3: Bolza vertices sit at Klein radius 0.985
    for s in (b.deform("F_A", 0.03), b.deform("F_B", -0.04), b.deform(v, 0.02)):
        for c in _residuals(s):
            assert c.contains(0)
        assert max(x.rad() for V in s.vertices_hyp for x in V) < 1e-20


def test_FA_linear_family_alone_is_not_closed():
    """Guard for the F_A sigma(t) correction: R_k = R + t(-1)^k alone misses the angle sum by ~2.83 t^2."""
    from flint import arb as A
    with prec(96):
        R = (3 + 2 * A(2).sqrt()).acosh()
        H = []
        for k in range(8):
            Rk, th = R + A("0.03") * (-1) ** k, A.pi() * k / 4
            H.append([Rk.cosh(), Rk.sinh() * th.cos(), Rk.sinh() * th.sin()])
        assert abs(float(surf._constraints_hyp(H)[4].mid())) > 1e-3
    assert Surface.bolza().deform("F_A", 0.03).r2_invariant


def test_FA_jets_first_order_match_vertices_dual():
    from bolzakit.certified.mesh import vertices_dual
    with prec(96):
        Vd = vertices_dual()
    for a, d in zip(Surface.bolza().vertex_jets("F_A"), Vd):
        for x, y in zip(a, d):
            assert (x.d - y.d).contains(0) and (x.v - y.v).contains(0)


@pytest.mark.parametrize("direction", ["F_A", "F_B", "vector"])
def test_vertex_jets_satisfy_constraints_to_second_order(direction):
    import numpy as np
    b = Surface.bolza()
    d = np.random.RandomState(2).randn(16) if direction == "vector" else direction
    with prec(96):
        for c in surf._constraints_hyp(b.vertex_jets(d)):
            # radii are Arb dependency growth of the (cancellation-heavy) check itself, from 1e-28 vertex balls
            for part, tol in zip(c.c, (1e-15, 1e-12, 1e-7)):
                assert part.contains(0) and part.rad() < tol


def test_vector_tangent_is_projected():
    """A direction normal to the constraint manifold has zero tangent after projection."""
    import numpy as np
    b = Surface.bolza()
    with prec(96):
        Jm = surf._jacobian(b._klein_arb())
    normal = np.array([float(Jm[0, i].mid()) for i in range(16)])
    t = b.klein_tangent(normal)
    assert np.max(np.abs(t)) < 1e-12


def test_r2_invariance_flags():
    import numpy as np
    b = Surface.bolza()
    assert b.r2_invariant and b.deform("F_B", 0.01).r2_invariant and b.deform("F_A", 0.01).r2_invariant
    assert not b.deform(np.ones(16), 0.01).r2_invariant


def test_vector_leaving_the_disk_raises():
    import numpy as np
    with pytest.raises(OutOfValidityRange):
        Surface.bolza().deform(np.random.RandomState(3).randn(16), 0.01)   # pushes a vertex outside the Klein disk


def test_vector_deformed_surface_runs_all_sectors(monkeypatch):
    import numpy as np
    from bolzakit.certified import forms
    seen = []
    orig = forms.compute_elements

    def spy(*a, **kw):
        seen.append(kw.get("sym"))
        return orig(*a, **kw)
    monkeypatch.setattr(forms, "compute_elements", spy)
    s = Surface.bolza().deform(0.05 * np.random.RandomState(3).randn(16), 0.01)
    r = certify_spectrum(s, k=3, degree=2, n=12)
    assert seen == [False] and r.certified
    seen.clear()
    certify_spectrum(Surface.bolza(), k=3, degree=2, n=12)
    assert seen == [True]


def test_klein_coords_bolza():
    import math
    kc = Surface.bolza().klein_coords()
    rho = math.tanh(math.acosh(3 + 2 * math.sqrt(2)))
    assert len(kc) == 16 and abs(kc[0] - rho) < 1e-15 and abs(kc[1]) < 1e-15
