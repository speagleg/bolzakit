import pytest
flint = pytest.importorskip("flint")
from flint import arb
from bolzakit.certified._arb import prec
from bolzakit.certified.mesh import Mesh, bolza_vertices_hyp, positions_arb, pairing_residual
from bolzakit.certified._hyp import det3


@pytest.fixture(autouse=True)
def _prec():
    with prec(128):
        yield


@pytest.mark.parametrize("n", [6, 12])
def test_euler_and_orientation(n):
    m = Mesh(n)
    E = len({tuple(sorted((m.lab[a], m.lab[b]))) for t in m.T for a, b in ((t[0], t[1]), (t[1], t[2]), (t[0], t[2]))})
    assert m.Nf - E + len(m.T) == -2
    X = positions_arb(m, bolza_vertices_hyp())
    assert all(det3(X[i0], X[i1], X[i2]) > 0 for (i0, i1, i2) in m.T)


def test_paired_nodes_coincide():
    m = Mesh(12); X = positions_arb(m, bolza_vertices_hyp())
    assert pairing_residual(m, X) < arb("1e-25")


def _mesh_pos():
    m = Mesh(6)
    return m, positions_arb(m, bolza_vertices_hyp())


def test_nan_node_is_not_a_pass():
    m, X = _mesh_pos()
    X[m.node(4, 6, 2)] = [arb("nan")] * 3
    assert not pairing_residual(m, X) < arb("1e-25")


def test_misplaced_node_is_not_a_pass():
    m, X = _mesh_pos()
    X[m.node(4, 6, 2)] = X[m.node(4, 6, 3)]
    assert not pairing_residual(m, X) < arb("1e-25")
