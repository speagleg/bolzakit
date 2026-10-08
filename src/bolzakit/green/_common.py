"""Shared helpers for item 03.  Imports item 02's exact Bolza group read-only (sys.path), never edits it.

Conventions (fixed for all of item 03):
  * Poincare disk, curvature -1, area form dA = 4 dx dy/(1-|z|^2)^2, Delta = (1-|z|^2)^2/4 * 4 d_z d_zbar.
  * A = area(Bolza) = 4 pi.
  * G1(d) = -(1/2pi) ln tanh(d/2)  (free Green function, -Delta G1 = delta).
  * u-variable: u = cosh d - 1 = 2|z-w|^2/((1-|z|^2)(1-|w|^2)).
"""
from __future__ import annotations
import math
import numpy as np

from ..exact import group as BG          # exact group (ported from item 02)
from ..provenance import provenance

A_AREA = 4 * math.pi
SQ2 = math.sqrt(2.0)
MU = math.sqrt(2 + 2 * SQ2)
ZETA = np.exp(1j * np.pi / 4)


def rows_to_ab(X: np.ndarray):
    """Rows (a0..a3, b0..b3) of exact group elements -> complex arrays (a, b) with gamma = [[a, b],[conj b, conj a]]."""
    X = np.asarray(X)
    zp = ZETA ** np.arange(4)
    A = X[:, :4].astype(float) @ zp
    B = (X[:, 4:].astype(float) @ zp) * MU
    return A, B


def ball(R: float):
    """All gamma with d(0, gamma 0) <= R (exact enumeration from item 02), as complex (a, b) arrays, identity first removed flag."""
    X = BG.ball_by_displacement(R)
    a, b = rows_to_ab(X)
    ident = (np.abs(b) < 1e-12) & (np.abs(a - 1) < 1e-12)
    return a, b, ident


def mobius(a, b, z):
    return (a * z + b) / (np.conj(b) * z + np.conj(a))


def dmobius(a, b, z):
    return 1.0 / (np.conj(b) * z + np.conj(a)) ** 2


def u_of(z, w):
    return 2 * np.abs(z - w) ** 2 / ((1 - np.abs(z) ** 2) * (1 - np.abs(w) ** 2))


def dist(z, w):
    # stability fix vs item 03 original (final review I1): 2 asinh(sqrt(u/2)) == arccosh(1 + u), no cancellation
    return 2 * np.arcsinh(np.sqrt(u_of(z, w) / 2))


@provenance("NUMERICAL", "closed form -(1/2pi) ln tanh(d/2) in float64", "float64 rounding")
def G1(d):
    return -np.log(np.tanh(d / 2)) / (2 * np.pi)
