"""Holomorphic 1-forms on the Bolza surface in the disk coordinate centred at the Weierstrass point z0 = 0
(octagon centre), by collocation of the automorphy condition in high precision (mpmath), and the Bergman kernel
B(z, w) near 0.

A holomorphic 1-form is f(z) dz with f(g z) g'(z) = f(z) for all g in Gamma.  rho: z -> zeta8 z is an automorphism
(item 02), so the 2-dim space splits into rho-eigenforms, f(z) = z^k sum_n c_n z^{8n} (residue class k).  Because
rho g0 rho^{-1} = g1 etc., for a rho-eigenform the automorphy under g0 implies it under all g_j, and the g_j generate
Gamma; so it suffices to collocate f(g0 z) g0'(z) = f(z) on points z of side 4 (g0 maps side 4 onto side 0).
The collocation matrix is very ill-conditioned in floating point (high powers live near the vertices), so it is
solved in mpmath at 50-60 digits; the null vector is separated from the rest of the spectrum by ~15 orders.

Bergman kernel, orthonormal in <f, g> = int_F f conj(g) dx dy:  B(z, w) = sum_i f_i(z) conj(f_i(w)) / ||f_i||^2.
Item 03 uses:  d_z d_wbar G_Sigma(z, w) = -(1/4) B(z, w) for z != w (note.tex, Lemma 2).
"""
from __future__ import annotations
import sys, math, time
import mpmath as mp


def setup(dps):
    mp.mp.dps = dps
    s2 = mp.sqrt(2)
    RM = mp.sqrt(s2 - 1)
    C = (RM + 1 / RM) / 2
    return dict(s2=s2, RM=RM, C=C, ch=1 + s2, sh=mp.sqrt(2 + 2 * s2), RV=mp.mpf(2) ** mp.mpf(-0.25))


def r_F(th, P):
    """Boundary radius of the octagon in direction th, for th in [-pi/8, pi/8] (side 0)."""
    return P['C'] * mp.cos(th) - mp.sqrt(P['C'] ** 2 * mp.cos(th) ** 2 - 1)


def solve_class(k, M, npts, P):
    pw = [8 * n + k for n in range(M)]
    rows = []
    for i in range(1, npts + 1):
        th = -mp.pi / 8 + mp.pi / 4 * i / (npts + 1)
        z = -r_F(th, P) * mp.expj(th)                      # point of side 4 (direction pi + th)
        gz = (P['ch'] * z + P['sh']) / (P['sh'] * z + P['ch'])
        dg = 1 / (P['sh'] * z + P['ch']) ** 2
        rows.append([(gz ** p * dg - z ** p) / P['RV'] ** p for p in pw])
    A = mp.matrix(rows)
    U, S, V = mp.svd_c(A)
    svals = [abs(S[i]) for i in range(len(S))]
    j = min(range(len(svals)), key=lambda i: svals[i])
    v = [mp.conj(V[j, n]) / P['RV'] ** pw[n] for n in range(M)]
    v0 = v[0]
    c = [x / v0 for x in v]
    return c, sorted(svals)


def norm2(k, c, P, nth=120):
    """int_F |f|^2 dx dy: exact radial integration of monomial products, Gauss-Legendre in theta on one sector
    (|f|^2 is C8-invariant since f is a rho-eigenform), times 8."""
    pw = [8 * n + k for n in range(len(c))]
    xs, ws = mp.gauss_quadrature(nth, 'legendre') if hasattr(mp, 'gauss_quadrature') else (None, None)
    if xs is None:
        # fall back: mp.quad
        def integrand(th):
            rf = r_F(th, P)
            s = mp.mpf(0)
            for a, p in zip(c, pw):
                for b, q in zip(c, pw):
                    s += mp.re(a * mp.conj(b) * mp.expj((p - q) * th)) * rf ** (p + q + 2) / (p + q + 2)
            return s
        return 8 * mp.quad(integrand, [-mp.pi / 8, 0, mp.pi / 8])


def norm2_fast(k, c, P):
    """Same integral, organised as int dth sum_{p,q} ... with mp.quad (degree-adaptive)."""
    pw = [8 * n + k for n in range(len(c))]
    cc = [(a, p) for a, p in zip(c, pw)]

    def integrand(th):
        rf = r_F(th, P)
        # g(r) = sum_p a_p r^p e^{i p th};  |f|^2 r dr integrated: sum_{p,q} a_p conj(a_q) e^{i(p-q)th} rf^{p+q+2}/(p+q+2)
        s = mp.mpf(0)
        terms = [(a * mp.expj(p * th) * rf ** p, p) for a, p in cc]
        for (x, p) in terms:
            for (y, q) in terms:
                s += mp.re(x * mp.conj(y)) / (p + q + 2)
        return s * rf ** 2
    return 8 * mp.quad(integrand, [-mp.pi / 8, -mp.pi / 16, 0, mp.pi / 16, mp.pi / 8])


def automorphy_check(k, c, P, ntest=30):
    """max relative |f(g z) g'(z) - f(z)| over z in F and a few group elements g with g z in F-ish (|g z|<=RV):
    uses g = g_j and g_j g_l (exact matrices via mp)."""
    import random
    random.seed(3)
    pw = [8 * n + k for n in range(len(c))]
    f = lambda z: sum(a * z ** p for a, p in zip(c, pw))
    zeta = mp.expjpi(mp.mpf(1) / 4)
    gens = []
    for j in range(8):
        e = zeta ** j
        gens.append(mp.matrix([[P['ch'], P['sh'] * e], [P['sh'] / e, P['ch']]]))
    worst = mp.mpf(0)
    count = 0
    while count < ntest:
        th = random.uniform(0, 2 * math.pi)
        thr = (th + math.pi / 8) % (math.pi / 4) - math.pi / 8
        r = random.uniform(0.3, 1.0) * r_F(mp.mpf(thr), P)
        z = r * mp.expj(th)
        for g in gens + [gens[a] * gens[b] for a in range(8) for b in range(8) if (a - b) % 8 != 4]:
            den = g[1, 0] * z + g[1, 1]
            gz = (g[0, 0] * z + g[0, 1]) / den
            if abs(gz) <= P['RV']:
                worst = max(worst, abs(f(gz) / den ** 2 - f(z)) / abs(f(z)))
                count += 1
    return worst


from ..provenance import provenance  # noqa: E402


@provenance("NUMERICAL", "Hejhal-type computation of holomorphic 1-forms (item 03 bergman.py)",
            "automorphy residual < 1e-23 at dps=70, M=40")
def bergman_at_weierstrass(dps: int = 50, M: int = 24, npts: int = 50):
    """B_Σ(0, 0) at the Weierstrass point 0, i.e. 1/||f_0||² with f_0 the normalised 1-form of residue class 0.
    Same call sequence as item 03 bergman.py's __main__ driver. An independent cross-check of the closed form."""
    saved = mp.mp.dps                       # setup() sets the global precision; restore it for the caller
    try:
        P = setup(dps)
        c, _sv = solve_class(0, M, npts, P)
        return 1 / norm2_fast(0, c, P)
    finally:
        mp.mp.dps = saved
