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
