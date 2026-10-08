"""G_Sigma, regular part R and Robin mass on the Bolza surface (heat-kernel splitting, ported from math-pursuits
item 03).

Method (heat-kernel splitting, heat.py):  for points z, w in a neighbourhood of the Weierstrass point 0,
    R(z, w)      = G_Sigma(z, w) - G1(d(z, w))
                 = sum_{gamma != 1} F_T(d(z, gamma w)) - E_T(d(z, w)) - T/A + Tail_T(z, w),
    G_Sigma(z,w) = G1(d(z,w)) + R(z,w),     m(z) = R(z, z) + (1/2pi) ln 2   (Robin mass, m = lim[G + (1/2pi) ln d]).
Tail_T = sum_{n>=1} e^{-lambda_n T} phi_n(z) phi_n(w)/lambda_n is NOT included (it is C8-invariant and smooth; at
T = 2.6 its size is <= (K(T,0,0) - 1/A)/lambda_1 ~ 4e-6 in value; its effect on the ring is quantified by varying T).
Normalisation:  -Delta G_Sigma = delta - 1/A, int G_Sigma dA = 0 (the image sum has int = T exactly; tail has mean 0).
Validity: |z|, |w| <= 0.45 (table ranges); images with d(0, gamma 0) <= D_ball.
"""
from __future__ import annotations
import math
import numpy as np
from ._common import ball, A_AREA, G1, mobius
from . import heat as H


COINCIDE_TOL = 1e-9


def _pdist(z, w):
    """Poincaré distance, array version: d = 2 asinh(|z-w| / sqrt((1-|z|^2)(1-|w|^2))).

    Algebraically identical to acosh(1 + 2|z-w|^2/((1-|z|^2)(1-|w|^2))) (cosh d = 1 + 2 sinh^2(d/2)), but without the
    cancellation in 1 + u for nearby points (final review I1)."""
    return 2 * np.arcsinh(np.abs(z - w) / np.sqrt((1 - np.abs(z) ** 2) * (1 - np.abs(w) ** 2)))


def _pdist1(z, w) -> float:
    """Scalar version of _pdist."""
    return 2 * math.asinh(abs(z - w) / math.sqrt((1 - abs(z) ** 2) * (1 - abs(w) ** 2)))


class BolzaGreen:
    def __init__(self, T=2.6, D_ball=15.0, zmax=0.45):
        self.T = T
        a, b, idn = ball(D_ball)
        self.a = a[~idn]; self.b = b[~idn]
        rh = 2 * math.atanh(zmax)
        lo = max(0.05, 2 * math.acosh(1 + math.sqrt(2)) - 2 * rh - 0.05)
        self.Ftab = H.ChebTable(lambda x: H.F_T(x, T), lo, D_ball + 2 * rh + 0.1, pieces=int(math.ceil(D_ball + 2 * rh + 0.1 - lo)) * 2)
        E = lambda x: G1(x) - H.F_T(x, T)
        self.Etab = H.ChebTable(lambda x: np.array([self._E_direct(xx) for xx in x]), 0.0, 2 * rh + 0.05, pieces=4)

    def _E_direct(self, d):
        """E_T(d) = int_T^inf k_t(d) dt; for d > 0.3 via G1 - F_T, for small d via the Abel integral of J_T."""
        T = self.T
        if d > 0.3:
            return float(G1(d) - H.F_T(np.array([d]), T)[0])
        # J_T(s) = int_T^inf t^{-3/2} e^{-t/4 - s^2/4t} dt = total - I_T
        from scipy.special import erfc
        def inner_E(s, TT):
            s = np.maximum(s, 1e-300)
            x = s / (2 * math.sqrt(TT)); y = math.sqrt(TT) / 2
            tot = 2 * math.sqrt(math.pi) * np.exp(-s / 2) / s
            return tot - H.inner_F(s, TT)
        # for tiny s the difference loses precision; use direct quadrature in t instead (smooth for t >= T)
        from scipy.integrate import quad
        def innerE2(s, TT):
            return np.array([quad(lambda t: t ** -1.5 * math.exp(-t / 4 - ss * ss / (4 * t)), TT, np.inf, epsabs=1e-18, epsrel=1e-13)[0] for ss in np.atleast_1d(s).ravel()]).reshape(np.shape(s))
        return float(H.abel(np.array([d]), innerE2, T, 40.0)[0])

    def R(self, z, w):
        """Regular part R(z, w) (Tail omitted) for scalar z, w."""
        gw = mobius(self.a, self.b, w)
        # stability fix vs item 03 original (final review I1): 2 asinh(|z-w|/sqrt(...)) instead of acosh(1 + 2|z-w|^2/...)
        d = _pdist(z, gw)
        d0 = _pdist1(z, w)
        return float(self.Ftab(d).sum() - self.Etab(np.array([d0]))[0] - self.T / A_AREA)

    def G(self, z, w):
        # stability fix vs item 03 original (final review I1)
        d0 = _pdist1(z, w)
        if d0 < COINCIDE_TOL:
            raise ValueError(f"points coincide on the surface (distance < {COINCIDE_TOL:g}); use robin_mass(z) or "
                             "regular(z, z)")
        return float(G1(d0)) + self.R(z, w)

    def robin_mass(self, z):
        return self.R(z, z) + math.log(2) / (2 * math.pi)

    def S_ring(self, zs):
        """Smooth part of the Hamiltonian: sum_{i<j} R(z_i, z_j) + (1/2) sum_i R(z_i, z_i)."""
        zs = np.asarray(zs)
        n = len(zs)
        # images of all z_j at once
        gw = mobius(self.a[:, None], self.b[:, None], zs[None, :])        # (images, n)
        tot = 0.0
        for i in range(n):
            zi = zs[i]
            d = _pdist(zi, gw[:, i:])        # stability fix vs item 03 original (final review I1)
            vals = self.Ftab(d.ravel()).reshape(d.shape).sum(axis=0)   # per j >= i
            for jj, j in enumerate(range(i, n)):
                d0 = _pdist1(zi, zs[j]) if j != i else 0.0      # stability fix vs item 03 original (final review I1)
                r = vals[jj] - self.Etab(np.array([d0]))[0] - self.T / A_AREA
                tot += r if j != i else r / 2
        return tot



# ---------------------------------------------------------------- bolzakit: Weierstrass charts and validity

import functools  # noqa: E402
from ..errors import OutOfValidityRange  # noqa: E402
from ..geometry import octagon as Geo  # noqa: E402
from ..provenance import provenance  # noqa: E402

VALIDITY_RADIUS = 0.45
# Covering radius of the 6 Weierstrass charts: rho = max_z min_i |z_local,i| = 0.4056164 (Poincare distance 0.8607),
# attained at z = -0.155223+0.374741j (|z| = rho, arg 5pi/8): equidistant (0.40562) from the centre chart and the two
# adjacent edge-midpoint charts; it is not equidistant from the vertex class.
# rho < VALIDITY_RADIUS, so every single point is valid; only a pair with no common chart raises.
_ERR = "omitted tail <= 4e-6 in value at T=2.6; ΔG = 1/A to 2e-5"

# (T, D_ball) pairs whose accuracy is validated: G(0.1, 0.2j) and robin_mass(0) agree with the (2.6, 15) values to
# 1.2e-7 and 8.2e-8 at (2.6, 13).  Rejected after measurement (final review C1): T = 2.0 with D_ball = 12, 13 or 15
# (-3.1e-5 in G, -3.9e-5 in the Robin mass: the omitted spectral tail grows like e^{-lambda_1 T}); and the review's
# T=1 (-1.6e-3), T=5 (-5.6e-5), D_ball=9 (-1.5e-4), T=10 (-2.6e-2, 25%).
VALIDATED_T_D = ((2.6, 15.0), (2.6, 13.0))


def check_T_D(T, D_ball) -> tuple[float, float]:
    """(float(T), float(D_ball)) if the pair is in VALIDATED_T_D (to 1e-12), else OutOfValidityRange. Cheap."""
    T, D_ball = float(T), float(D_ball)
    for t, d in VALIDATED_T_D:
        if abs(T - t) <= 1e-12 and abs(D_ball - d) <= 1e-12:
            return t, d
    raise OutOfValidityRange(
        f"(T, D_ball) = ({T!r}, {D_ball!r}) is not validated; allowed pairs: "
        + ", ".join(f"({t}, {d})" for t, d in VALIDATED_T_D)
        + ". Smaller T leaves a larger omitted spectral tail and larger T needs a larger image ball; un-validated "
        "pairs can be off by up to 25% while the error label (<= 4e-6) would not apply.")


class SurfaceGreen:
    """Green's function of the Bolza surface: -Δ G = δ - 1/A, ∫ G dA = 0 (heat-kernel splitting, item 03).
    Validated only when both points lie within Poincaré radius 0.45 of the SAME Weierstrass point (chart coordinates)."""

    def __init__(self, T: float = 2.6, D_ball: float = 15.0):
        from ..exact.automorphisms import automorphisms
        T, D_ball = check_T_D(T, D_ball)
        self._core = BolzaGreen(T=T, D_ball=D_ball, zmax=VALIDITY_RADIUS)
        A = automorphisms()
        self._charts = [A.chart_to(p) for p in A.weierstrass_points()]

    def _locals(self, z: complex) -> list[complex]:
        """Local coordinate of z in every Weierstrass chart (Γ-reduced, so |z_local| = tanh(d(z, p_i)/2))."""
        return [Geo.reduce_point(Geo.mobius(M, z))[0] for M in self._charts]

    def chart(self, z: complex) -> tuple[int, complex]:
        """(index, z_local) of the Weierstrass chart in which z has the smallest |z_local|."""
        best = None
        for i, zl in enumerate(self._locals(z)):
            if best is None or abs(zl) < abs(best[1]):
                best = (i, zl)
        return best

    def _local_pair(self, z, w):
        # The 0.45-discs of adjacent Weierstrass points overlap (their distance 1.5286 < 2 * 2 atanh 0.45), so the
        # chart nearest z alone may reject a valid pair.  Use the chart minimising max(|z_local|, |w_local|).
        # Within a chart, reduce_point returns the Γ-image of w nearest the chart centre, the unique one within
        # 0.45 of it; the core image sum runs over all of Γ, so this is the correct representative.
        zs, ws = self._locals(z), self._locals(w)
        i = min(range(len(zs)), key=lambda k: (max(abs(zs[k]), abs(ws[k])), k))
        zl, wl0 = zs[i], ws[i]
        if abs(zl) > VALIDITY_RADIUS or abs(wl0) > VALIDITY_RADIUS:
            raise OutOfValidityRange(
                f"points must lie within Poincaré radius {VALIDITY_RADIUS} of the same Weierstrass point "
                f"(best chart {i}: |z|={abs(zl):.3f}, |w|={abs(wl0):.3f}); general-position G_Σ is not in v0.1")
        return zl, wl0

    @provenance("NUMERICAL", "heat-kernel splitting (item 03) in a Weierstrass chart", _ERR)
    def G(self, z: complex, w: complex) -> float:
        zl, wl = self._local_pair(z, w)
        return self._core.G(zl, wl)       # raises ValueError if d(zl, wl) < COINCIDE_TOL (coincident on the surface)

    @provenance("NUMERICAL", "heat-kernel splitting (item 03) in a Weierstrass chart", _ERR)
    def regular(self, z: complex, w: complex) -> float:
        zl, wl = self._local_pair(z, w)
        return self._core.R(zl, wl)

    @provenance("NUMERICAL", "heat-kernel splitting (item 03) in a Weierstrass chart", _ERR)
    def robin_mass(self, z: complex) -> float:
        zl, _ = self._local_pair(z, z)
        return self._core.robin_mass(zl)


@functools.lru_cache(maxsize=4)
def _green(T: float, D_ball: float) -> SurfaceGreen:
    return SurfaceGreen(T, D_ball)


@provenance("NUMERICAL", "heat-kernel image sum over exact group elements + omitted spectral tail (item 03)",
            _ERR)
def green(T: float = 2.6, D_ball: float = 15.0) -> SurfaceGreen:
    """Green's function of the Bolza surface (cached per (T, D_ball)).

    Only the validated pairs in VALIDATED_T_D, (2.6, 15.0) and (2.6, 13.0), are accepted; any other pair raises
    OutOfValidityRange before anything is built."""
    return _green(*check_T_D(T, D_ball))
