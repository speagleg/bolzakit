"""scikit-fem adapter: assemble K, M and their parameter derivatives for ``cluster_sensitivity``.

Imports nothing from the rest of bolzakit. The derivatives are central finite differences of the assembled
matrices (first derivative O(h^2), second derivative O(h^2) truncation plus O(eps/h^2) round-off).
"""
from __future__ import annotations

import numpy as np

from bolzakit.sensitivity.adapters._util import _dof_indices


def _require_skfem():
    try:
        import skfem
    except ImportError as e:  # local stand-in for bolzakit.provenance.require
        raise ImportError("pip install bolzakit[skfem]") from e
    return skfem


def from_forms(basis, stiffness_form, mass_form, p0: float, h: float = 1e-5, dirichlet_dofs=None):
    """Return ``(K, M, dK, dM, d2K, d2M)`` as sparse matrices at ``p0``.

    The forms are ``skfem.BilinearForm`` objects reading the parameter from ``w["p"]``. Dirichlet dofs are removed
    (only the interior dofs are kept). ``dirichlet_dofs`` may be an index array or the dict returned by
    ``basis.get_dofs()``; None keeps all dofs."""
    _require_skfem()
    from skfem import asm
    if not h > 0:
        raise ValueError("h must be positive")

    def at(p):
        return asm(stiffness_form, basis, p=p), asm(mass_form, basis, p=p)

    K0, M0 = at(p0)
    Kp, Mp = at(p0 + h)
    Km, Mm = at(p0 - h)
    n = K0.shape[0]
    drop = _dof_indices(dirichlet_dofs, n)
    keep = np.setdiff1d(np.arange(K0.shape[0]), drop)

    def cond(X):
        return X.tocsr()[keep][:, keep].tocsc()

    K0, M0, Kp, Mp, Km, Mm = (cond(X) for X in (K0, M0, Kp, Mp, Km, Mm))
    dK, dM = (Kp - Km) / (2 * h), (Mp - Mm) / (2 * h)
    d2K, d2M = (Kp - 2 * K0 + Km) / h**2, (Mp - 2 * M0 + Mm) / h**2
    return K0, M0, dK, dM, d2K, d2M
