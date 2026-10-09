import builtins, sys
import numpy as np, pytest
from bolzakit.sensitivity import cluster_sensitivity, from_callables, ClusterSensitivity, DISCRETE_CAVEAT
from bolzakit.sensitivity.errors import (SensitivityError, PencilError, ClusterRangeError, ClusterSplitError,
                                         NonSymmetricFunctionError, MissingExtra, FiniteDifferenceWarning)


def test_hierarchy():
    for cls in (PencilError, ClusterRangeError, ClusterSplitError, NonSymmetricFunctionError):
        assert issubclass(cls, SensitivityError) and issubclass(cls, ValueError)
    assert issubclass(MissingExtra, SensitivityError) and issubclass(MissingExtra, ImportError)
    assert issubclass(FiniteDifferenceWarning, UserWarning)
    assert SensitivityError.__mro__[1] is Exception


def test_pencil_error_nonsymmetric_and_not_pd():
    with pytest.raises(PencilError, match="K and M must be symmetric with M positive definite") as e:
        cluster_sensitivity(np.array([[1.0, 2], [0, 1]]), np.eye(2), np.eye(2), None, (0, 0))
    assert type(e.value) is PencilError
    with pytest.raises(PencilError) as e:
        cluster_sensitivity(np.eye(2), np.diag([1.0, 0.0]), np.eye(2), None, (0, 0))
    assert type(e.value) is PencilError
    with pytest.raises(PencilError) as e:
        cluster_sensitivity(np.eye(2), np.eye(3), np.eye(2), None, (0, 0))
    assert type(e.value) is PencilError


@pytest.mark.parametrize("cl", [(-1, 0), (1, 0), (0, 3)])
def test_cluster_range_error(cl):
    with pytest.raises(ClusterRangeError) as e:
        cluster_sensitivity(np.diag([1.0, 2, 3]), np.eye(3), np.eye(3), None, cl)
    assert type(e.value) is ClusterRangeError


def test_cluster_split_error_keeps_word_splits():
    with pytest.raises(ClusterSplitError, match="splits") as e:
        cluster_sensitivity(np.diag([1.0, 1.0, 3.0]), np.eye(3), np.eye(3), None, (0, 0))
    assert type(e.value) is ClusterSplitError


def test_from_callables_emits_fd_warning():
    with pytest.warns(FiniteDifferenceWarning, match="finite-difference"):
        from_callables(lambda p: np.diag([1.0 + p[0], 2.0, 3.0]), lambda p: np.eye(3), np.zeros(1), np.ones(1), (0, 0))


def test_skfem_adapter_emits_fd_warning():
    skfem = pytest.importorskip("skfem")
    from skfem import MeshQuad, ElementQuad1, Basis, BilinearForm
    from bolzakit.sensitivity.adapters.skfem import from_forms

    @BilinearForm
    def a(u, v, w):
        return (1 + w["p"]) * (u.grad[0] * v.grad[0] + u.grad[1] * v.grad[1])

    @BilinearForm
    def m(u, v, w):
        return u * v

    with pytest.warns(FiniteDifferenceWarning):
        from_forms(Basis(MeshQuad().refined(2), ElementQuad1()), a, m, 0.0)


def test_skfem_missing_raises_missing_extra(monkeypatch):
    from bolzakit.sensitivity.adapters import skfem as ad
    real_import = builtins.__import__

    def fake(name, *a, **k):
        if name == "skfem" or name.startswith("skfem."):
            raise ImportError("no skfem")
        return real_import(name, *a, **k)

    monkeypatch.setattr(builtins, "__import__", fake)
    with pytest.raises(MissingExtra, match=r"bolzakit\[skfem\]") as e:
        ad._require_skfem()
    assert isinstance(e.value, ImportError)


def test_fenicsx_missing_raises_missing_extra(monkeypatch):
    from bolzakit.sensitivity.adapters import fenicsx as ad
    monkeypatch.setitem(sys.modules, "dolfinx", None)
    with pytest.raises(MissingExtra, match="dolfinx"):
        ad._require_dolfinx()


def test_new_result_fields_and_positional_construction():
    r = cluster_sensitivity(np.diag([1.0, 1.0, 3.0]), np.eye(3), np.eye(3), None, (0, 1))
    assert r.cluster == (0, 1) and r.lam_bar == 1.0 and r.gap_ratio >= 1 and r.caveat == DISCRETE_CAVEAT
    old = ClusterSensitivity(np.ones(1), np.ones((1, 1)), np.ones((1, 1)), np.ones(1), None, None)
    assert old.label == "NUMERICAL" and old.cluster is None
