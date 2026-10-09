"""Helpers shared by the FEM adapters (numpy only)."""
from __future__ import annotations

import numpy as np


def _dof_indices(spec, n: int):
    """Normalise a Dirichlet-dof specification to a sorted unique int array, validated against range(n)."""
    if spec is None:
        return np.array([], dtype=int)
    if isinstance(spec, dict):
        parts = [_dof_indices(v, n) for v in spec.values()]
        return np.unique(np.concatenate(parts)) if parts else np.array([], dtype=int)
    if not isinstance(spec, (np.ndarray, list, tuple)) and hasattr(spec, "all"):  # skfem DofsView
        spec = spec.all()
    arr = np.asarray(spec)
    if arr.size == 0:
        return np.array([], dtype=int)
    if arr.dtype.kind not in "iu":
        raise TypeError("dirichlet_dofs must hold integer dof indices")
    arr = arr.astype(int).ravel()
    if arr.min() < 0 or arr.max() >= n:
        raise ValueError(f"dirichlet_dofs indices must lie in [0, {n})")
    return np.unique(arr)
