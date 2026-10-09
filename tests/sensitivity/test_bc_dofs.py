import numpy as np, pytest
from bolzakit.sensitivity.adapters._util import _dof_indices


class _BC:
    def __init__(self, d): self._d = np.array(d)
    def dof_indices(self): return self._d, 1


def _bc_dofs():
    # fenicsx.py imports dolfinx lazily, so the helper is importable without it
    from bolzakit.sensitivity.adapters.fenicsx import _bc_dofs
    return _bc_dofs


def test_bc_dofs_normalisation_without_dolfinx():
    f = _bc_dofs()
    assert f(None, 10).size == 0
    assert f(np.array([], dtype=int), 10).size == 0
    assert list(f(np.array([3, 1, 3, 7]), 10)) == [1, 3, 7]
    assert list(f([5, 2], 10)) == [2, 5]
    assert list(f([_BC([4, 2]), _BC([2, 9])], 10)) == [2, 4, 9]
    assert f([], 10).size == 0
    with pytest.raises(ValueError):
        f(np.array([10]), 10)
