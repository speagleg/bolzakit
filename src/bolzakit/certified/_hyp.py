"""Hyperboloid-model primitives (Minkowski form diag(1,-1,-1)); ported from the Bolza lambda1 research code.

mink/tohyp/tok/hdist are float (numpy) helpers (SRC octagon.py); Bm/det3/hyp_interp are generic and work on arb
triples as well as on jet (jets.J) and dual (jets.D) triples (SRC cert_rr_p2.py version).
"""
import numpy as np


def mink(a, b):
    return a[..., 0] * b[..., 0] - a[..., 1] * b[..., 1] - a[..., 2] * b[..., 2]


def tohyp(p):
    p = np.asarray(p, float)
    t = 1 / np.sqrt(1 - np.sum(p * p, -1))
    return np.stack([t, t * p[..., 0], t * p[..., 1]], -1)


def tok(X):
    return X[..., 1:] / X[..., :1]


def hdist(a, b):
    return np.arccosh(np.maximum(mink(a, b), 1.0))


def Bm(x, y):
    return x[0] * y[0] - x[1] * y[1] - x[2] * y[2]


def det3(c0, c1, c2):
    return (c0[0] * (c1[1] * c2[2] - c1[2] * c2[1]) - c1[0] * (c0[1] * c2[2] - c0[2] * c2[1])
            + c2[0] * (c0[1] * c1[2] - c0[2] * c1[1]))


def hyp_interp(A, Bp, s):
    """point at fraction s (arb) of the hyperbolic arclength from A to Bp (arb, jet or dual triples)."""
    c = Bm(A, Bp)
    L = c.acosh(); sh = L.sinh()
    wa = (L * (1 - s)).sinh() / sh; wb = (L * s).sinh() / sh
    return [wa * A[i] + wb * Bp[i] for i in range(3)]
