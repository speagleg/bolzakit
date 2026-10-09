"""FEniCSx (dolfinx) adapter: assemble K, M and their exact parameter derivatives for ``cluster_sensitivity``.

Imports nothing from bolzakit outside ``bolzakit.sensitivity``. The parameter must enter the UFL forms through
``p = ufl.variable(dolfinx.fem.Constant(mesh, value))``; the derivative forms are ``ufl.diff(a, p)`` and
``ufl.diff(ufl.diff(a, p), p)``, so no finite differencing is involved.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from bolzakit.sensitivity.adapters._util import _dof_indices


def _require_dolfinx():
    try:
        import dolfinx  # noqa: F401
        import dolfinx.fem.petsc  # noqa: F401
        import ufl
    except ImportError as e:
        raise ImportError(
            "dolfinx not available: install via conda or the dolfinx/dolfinx docker image") from e
    return ufl


def _bc_dofs(bcs, n: int):
    """Normalise ``bcs`` (None, an index array/list, or DirichletBC-like objects with ``dof_indices()``) to dofs."""
    if bcs is None:
        return np.array([], dtype=int)
    if isinstance(bcs, (list, tuple)) and len(bcs) > 0 and all(hasattr(b, "dof_indices") for b in bcs):
        return _dof_indices(np.concatenate([np.asarray(b.dof_indices()[0], dtype=int) for b in bcs]), n)
    return _dof_indices(bcs, n)


def from_forms(V, a, m, p, bcs=None):
    """Return ``(K, M, dK, dM, d2K, d2M)`` as scipy CSR matrices at the current value of ``p``.

    ``a`` and ``m`` are UFL bilinear forms on ``V`` built from ``p``, which must be a ``ufl.variable`` wrapping a
    ``dolfinx.fem.Constant``. Rows and columns of the Dirichlet dofs in ``bcs`` (a list of
    ``dolfinx.fem.DirichletBC``, or an index array) are removed; None keeps all dofs. Serial (single MPI rank) only."""
    ufl = _require_dolfinx()
    from dolfinx import fem
    from dolfinx.fem.petsc import assemble_matrix
    if not isinstance(p, ufl.classes.Variable):
        raise TypeError("p must be a ufl.variable wrapping a dolfinx.fem.Constant")
    if V.mesh.comm.size != 1:
        raise NotImplementedError("fenicsx adapter supports a single MPI rank only")

    def build(form):
        if not form.integrals():  # identically zero after differentiation
            n = V.dofmap.index_map.size_global * V.dofmap.index_map_bs
            return sp.csr_matrix((n, n))
        A = assemble_matrix(fem.form(form))
        A.assemble()
        indptr, indices, data = A.getValuesCSR()
        out = sp.csr_matrix((data, indices, indptr), shape=A.getSize())
        A.destroy()
        return out

    forms = []
    for f in (a, m):
        d1 = ufl.diff(f, p)
        forms += [f, d1, ufl.diff(d1, p)]
    K, dK, d2K, M, dM, d2M = (build(f) for f in forms)
    n = K.shape[0]
    drop = _bc_dofs(bcs, n)
    keep = np.setdiff1d(np.arange(n), drop)
    return tuple(X[keep][:, keep].tocsr() for X in (K, M, dK, dM, d2K, d2M))
