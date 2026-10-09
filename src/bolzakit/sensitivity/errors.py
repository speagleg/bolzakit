"""Exception and warning types of ``bolzakit.sensitivity`` (stdlib bases only).

This hierarchy is deliberately separate from ``bolzakit.errors``: ``bolzakit.sensitivity`` imports nothing from the
rest of bolzakit (layering rule), so it cannot subclass ``BolzakitError``. Every class also subclasses the builtin
that v0.2 raised (``ValueError`` or ``ImportError``), so ``except ValueError`` keeps working.
"""
from __future__ import annotations

__all__ = ["SensitivityError", "PencilError", "ClusterRangeError", "ClusterSplitError",
           "NonSymmetricFunctionError", "MissingExtra", "FiniteDifferenceWarning"]


class SensitivityError(Exception):
    """Base class of all bolzakit.sensitivity errors."""


class PencilError(SensitivityError, ValueError):
    """K or M is not square/symmetric, M is not positive definite, or a trial basis W is rank-deficient."""


class ClusterRangeError(SensitivityError, ValueError):
    """``cluster = (i, j)`` is out of range or has i > j."""


class ClusterSplitError(SensitivityError, ValueError):
    """The cluster is not isolated under the applicable gap rule (it splits a (near-)degenerate group)."""


class NonSymmetricFunctionError(SensitivityError, ValueError):
    """``cluster_function``'s symmetry probe found f not symmetric in the cluster eigenvalues."""


class MissingExtra(SensitivityError, ImportError):
    """An optional dependency (scikit-fem, dolfinx) is missing; the message names the install route."""


class FiniteDifferenceWarning(UserWarning):
    """Derivatives were obtained by finite differences (truncation and round-off error apply)."""
