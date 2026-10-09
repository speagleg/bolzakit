import numpy as np, pytest
pytest.importorskip("dolfinx")
import ufl
from mpi4py import MPI
from dolfinx import fem, mesh
from bolzakit.sensitivity import cluster_sensitivity
from bolzakit.sensitivity.adapters.fenicsx import from_forms


def _setup(n):
    msh = mesh.create_unit_square(MPI.COMM_WORLD, n, n, mesh.CellType.quadrilateral)
    V = fem.functionspace(msh, ("Lagrange", 2))
    p = ufl.variable(fem.Constant(msh, 0.0))
    u, v = ufl.TrialFunction(V), ufl.TestFunction(V)
    a = ((1 / (1 + p)) * u.dx(0) * v.dx(0) + (1 + p) * u.dx(1) * v.dx(1)) * ufl.dx
    m = (1 + p) * u * v * ufl.dx
    tdim = msh.topology.dim
    msh.topology.create_connectivity(tdim - 1, tdim)
    facets = mesh.exterior_facet_indices(msh.topology)
    dofs = fem.locate_dofs_topological(V, tdim - 1, facets)
    bc = fem.dirichletbc(fem.Constant(msh, 0.0), dofs, V)
    return V, a, m, p, bc


def test_rectangle_cluster_rates_and_curvature():
    V, a, m, p, bc = _setup(32)
    K, M, dK, dM, d2K, d2M = from_forms(V, a, m, p, bcs=[bc])
    r = cluster_sensitivity(K, M, dK, dM, (1, 2), d2K=d2K, d2M=d2M)
    pi2 = np.pi ** 2
    assert np.allclose(np.sort(r.rates), [-8 * pi2, -2 * pi2], rtol=2e-3)
    assert abs(r.trace_hessian - 30 * pi2) / (30 * pi2) < 5e-3
    assert abs(r.kyfan_upper - 36 * pi2) / (36 * pi2) < 5e-3


def test_requires_variable_and_removes_bc_dofs():
    V, a, m, p, bc = _setup(4)
    with pytest.raises(TypeError):
        from_forms(V, a, m, p.ufl_operands[0], bcs=[bc])
    n = V.dofmap.index_map.size_global
    K = from_forms(V, a, m, p, bcs=[bc])[0]
    assert K.shape == (n - len(bc.dof_indices()[0]),) * 2
