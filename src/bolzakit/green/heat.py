"""Heat-kernel splitting for the Green function of the Bolza surface (method 2, independent of the Bergman route).

    G_Sigma(z, w) = int_0^inf (K_Sigma(t, z, w) - 1/A) dt
                  = sum_gamma F_T(d(z, gamma w)) - T/A + Tail_T(z, w),
    F_T(d) = int_0^T k_t(d) dt,   Tail_T(z, w) = int_T^inf (K_Sigma(t, z, w) - 1/A) dt
                  = sum_{n>=1} e^{-lambda_n T} lambda_n^{-1} phi_n(z) phi_n(w),
with k_t the heat kernel of H^2 (curvature -1).  The image sum converges like a Gaussian (k_t ~ e^{-d^2/4t}), so
it is absolutely convergent and NOT ball-truncation-sensitive once D_max is large enough (checked).

Free heat kernel (McKean):  k_t(d) = sqrt2 e^{-t/4}/(4 pi t)^{3/2} int_d^inf s e^{-s^2/4t} / sqrt(cosh s - cosh d) ds.
Time integral in closed form:  int_0^T t^{-3/2} e^{-t/4 - s^2/(4t)} dt
                              = (sqrt(pi)/s) e^{-x^2-y^2} [erfcx(x - y) + erfcx(x + y)],  x = s/(2 sqrt T), y = sqrt(T)/2.
(T -> inf gives F_inf = G1 = -(1/2pi) ln tanh(d/2) via Mehler's formula for Q_0 -- checked numerically below.)

Regular part:  R(z, w) = G_Sigma - G1 = sum_{gamma != 1} F_T(d(z, gamma w)) - E_T(d(z, w)) - T/A + Tail_T,
E_T = G1 - F_T = int_T^inf k_t dt (smooth at d = 0).
"""
from __future__ import annotations
import math
import numpy as np
from scipy.special import erfcx, erfc
from numpy.polynomial import chebyshev as C
from ._common import A_AREA

CNST = math.sqrt(2) / (4 * math.pi) ** 1.5
_GX, _GW = np.polynomial.legendre.leggauss(240)


def inner_F(s, T):
    """int_0^T t^{-3/2} e^{-t/4 - s^2/4t} dt."""
    x = s / (2 * math.sqrt(T)); y = math.sqrt(T) / 2
    return math.sqrt(math.pi) / s * np.exp(-x * x - y * y) * (erfcx(x - y) + erfcx(x + y))


def inner_k(s, t):
    """t^{-3/2} e^{-t/4 - s^2/4t}  (so that k_t = CNST * Abel[s * inner_k])."""
    return t ** -1.5 * np.exp(-t / 4 - s * s / (4 * t))


def abel(d, inner, par, Tscale):
    """CNST * int_d^inf s inner(s) / sqrt(cosh s - cosh d) ds, via s = d + v^2 and
    cosh s - cosh d = 2 sinh((s+d)/2) sinh(v^2/2).  Vectorised over d (array)."""
    d = np.atleast_1d(np.asarray(d, float))
    Vmax = math.sqrt(math.sqrt(200.0 * Tscale) + 1.0)
    v = (_GX + 1) / 2 * Vmax
    wv = _GW * Vmax / 2
    s = d[:, None] + v[None, :] ** 2
    # 2v / sqrt(2 sinh(v^2/2)) -> smooth; at v -> 0 equals 2/1 * ... (limit 2)
    vv = v[None, :]
    fac = np.where(vv > 1e-8, 2 * vv / np.sqrt(2 * np.sinh(vv * vv / 2)), 2.0)
    integrand = fac * s * inner(s, par) / np.sqrt(np.sinh((s + d[:, None]) / 2))
    return CNST * integrand @ wv


def F_T(d, T):
    return abel(d, inner_F, T, T)


def k_t(d, t):
    return abel(d, inner_k, t, t)


def k_t0(t):
    """k_t(0) (Abel integral at d = 0 has integrand s/sqrt(cosh s - 1) -> sqrt2; handled by same formula)."""
    return float(abel(np.array([0.0]), inner_k, t, t)[0])


class ChebTable:
    """Piecewise Chebyshev interpolant of f on [lo, hi] (unit pieces) in variable d; returns f, f', f''."""

    def __init__(self, f, lo, hi, pieces=None, deg=28):
        self.lo, self.hi = lo, hi
        n = pieces or int(math.ceil(hi - lo))
        self.edges = np.linspace(lo, hi, n + 1)
        xc = np.cos(np.pi * (np.arange(deg + 1) + 0.5) / (deg + 1))
        self.coef = []
        for a, b in zip(self.edges[:-1], self.edges[1:]):
            dd = (a + b) / 2 + (b - a) / 2 * xc
            self.coef.append(C.chebfit(xc, f(dd), deg))
        self.coef = np.array(self.coef)
        self.c1 = np.array([C.chebder(c) for c in self.coef])
        self.c2 = np.array([C.chebder(c, 2) for c in self.coef])

    def __call__(self, d, der=0):
        d = np.asarray(d, float)
        if d.size and (d.min() < self.lo - 1e-12 or d.max() > self.hi + 1e-12):
            raise ValueError("ChebTable out of range [%g,%g]: %g..%g" % (self.lo, self.hi, d.min(), d.max()))
        idx = np.clip(np.searchsorted(self.edges, d) - 1, 0, len(self.edges) - 2)
        a = self.edges[idx]; b = self.edges[idx + 1]
        x = (2 * d - a - b) / (b - a)
        cs = [self.coef, self.c1, self.c2][der][idx]
        # Clenshaw per point (vectorised)
        out = np.zeros_like(x)
        b1 = np.zeros_like(x); b2 = np.zeros_like(x)
        for j in range(cs.shape[1] - 1, 0, -1):
            b1, b2 = 2 * x * b1 - b2 + cs[:, j], b1
        out = x * b1 - b2 + cs[:, 0]
        return out * (2 / (b - a)) ** der


def dd_mixed_at_origin(ftab, a, b):
    """d_z d_wbar f(d(z, gamma w)) at z = w = 0 for each gamma = [[a,b],[conj b, conj a]] (vectorised), using
    D = -conj(gamma'(0)) ((cosh d + 1)^2/2) [u phi''(u) + phi'(u)],  phi(u) = f(arccosh(1+u))."""
    q = b / np.conj(a)
    r2 = np.abs(q) ** 2
    u = 2 * r2 / (1 - r2)
    d = np.arccosh(1 + u)
    f1 = ftab(d, 1); f2 = ftab(d, 2)
    sh = np.sinh(d); ch = np.cosh(d)
    p1 = f1 / sh
    p2 = (f2 - ch / sh * f1) / sh ** 2
    gp = 1.0 / np.conj(a) ** 2
    return -np.conj(gp) * ((ch + 1) ** 2 / 2) * (u * p2 + p1)


def self_check():
    """F_T -> G1 as T -> large; k_t integrates to 1 over H^2; Delta F_T = k_T (radial ODE) at a few d."""
    from ._common import G1
    d = np.array([0.5, 1.0, 3.0, 6.0])
    print("F_T(d) for T=200 vs G1:", np.max(np.abs(F_T(d, 200.0) - G1(d))))
    for t in (0.5, 2.0):
        dd = np.linspace(0, 40 * math.sqrt(t) + 8 * t, 4001)[1:]
        mass = np.trapezoid(k_t(dd, t) * 2 * np.pi * np.sinh(dd), dd)
        print(f"int k_t dA (t={t}) = {mass:.12f}")
    T = 2.0
    tab = ChebTable(lambda x: F_T(x, T), 1.0, 8.0)
    dd = np.array([1.5, 3.2, 5.7])
    lap = tab(dd, 2) + np.cosh(dd) / np.sinh(dd) * tab(dd, 1)
    print("Delta F_T - k_T (should be 0):", lap - k_t(dd, T))


if __name__ == '__main__':
    self_check()
