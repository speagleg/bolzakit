import numpy as np, pytest
skfem = pytest.importorskip("skfem")
from skfem import MeshQuad, ElementQuad2, Basis, BilinearForm
from bolzakit.sensitivity import cluster_sensitivity
from bolzakit.sensitivity.adapters.skfem import from_forms

@BilinearForm
def kform(u, v, w):
    p = w["p"]; return (1 / (1 + p)) * u.grad[0] * v.grad[0] + (1 + p) * u.grad[1] * v.grad[1]

@BilinearForm
def mform(u, v, w):
    return (1 + w["p"]) * u * v

def test_rectangle_cluster_rates_and_curvature():
    mesh = MeshQuad().refined(5); basis = Basis(mesh, ElementQuad2())
    K, M, dK, dM, d2K, d2M = from_forms(basis, kform, mform, 0.0, dirichlet_dofs=basis.get_dofs())
    r = cluster_sensitivity(K, M, dK, dM, (1, 2), d2K=d2K, d2M=d2M)
    pi2 = np.pi ** 2
    assert np.allclose(np.sort(r.rates), [-8 * pi2, -2 * pi2], rtol=2e-3)
    assert abs(r.trace_hessian - 30 * pi2) / (30 * pi2) < 5e-3
    assert abs(r.kyfan_upper - 36 * pi2) / (36 * pi2) < 5e-3


def _small():
    return Basis(MeshQuad().refined(3), ElementQuad2())


@pytest.mark.parametrize("make", [
    lambda b: b.get_dofs(),
    lambda b: b.get_dofs().all(),
    lambda b: [int(i) for i in b.get_dofs().all()],
    lambda b: {"all": b.get_dofs()},
    lambda b: {k: v for k, v in b.get_dofs().nodal.items()} | {"f": b.get_dofs().facet},
])
def test_dirichlet_forms_keep_interior(make):
    b = _small()
    K, M, *_ = from_forms(b, kform, mform, 0.0, dirichlet_dofs=make(b))
    assert K.shape == M.shape == (225, 225) and b.N == 289


def test_dirichlet_none_keeps_all_and_validates():
    b = _small()
    assert from_forms(b, kform, mform, 0.0)[0].shape == (289, 289)
    with pytest.raises(ValueError):
        from_forms(b, kform, mform, 0.0, dirichlet_dofs=[289])
    with pytest.raises(TypeError):
        from_forms(b, kform, mform, 0.0, dirichlet_dofs=np.array([True, False]))
    with pytest.raises(ValueError):
        from_forms(b, kform, mform, 0.0, h=0.0)
