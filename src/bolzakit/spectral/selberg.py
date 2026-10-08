"""Selberg trace formula consistency test (README C5) with the heat-kernel test function.

Formula used (Strohmaier-Uski 2013, eq. (18), quoted from arXiv:1110.2150v4 p.29, attributed there to
[M04] eq. (193)); for a compact hyperbolic surface M:

  tr e^{-t Delta} = Vol(M) e^{-t/4}/(4 pi t) * int_0^inf pi e^{-r^2 t} / cosh^2(pi r) dr
                    + sum_{n>=1} sum_{gamma} e^{-t/4} l(gamma) e^{-n^2 l(gamma)^2/(4t)} / ( sqrt(4 pi t) * 2 sinh(n l(gamma)/2) ),

where gamma runs over PRIMITIVE closed geodesics.  The text does not say "oriented" explicitly; the
standard form (test function h(r)=e^{-t(1/4+r^2)}, sum over all hyperbolic conjugacy classes of Gamma)
counts gamma and gamma^-1 separately, i.e. oriented geodesics.  We use oriented counts and let the
numerical balance confirm the convention (with unoriented counts the residual is O(1e-2), see tests).

Spectral side: 1 (lambda_0 = 0) + sum over the 1000 S-U eigenvalues; tail beyond lambda_1000 ~ 998 is
bounded by Weyl (N(lam) ~ lam): int_{998}^inf e^{-t lam} d lam = e^{-998 t}/t (negligible for t >= 0.05).
Geometric side truncated at the enumerated half-trace bound; the omitted tail is estimated (not bounded)
using the prime geodesic theorem count ~ e^l / l per unit length (labelled ESTIMATE).
"""
from __future__ import annotations

import math

import mpmath as mp

from ..exact import lengths as LS
from . import eigen as SP


def identity_term(t):
    t = mp.mpf(t)
    f = lambda r: mp.pi * mp.exp(-r * r * t) / mp.cosh(mp.pi * r) ** 2
    I = mp.quad(f, [0, 1, 3, 8, mp.inf])
    return SP.AREA * mp.exp(-t / 4) / (4 * mp.pi * t) * I


def hyperbolic_term(t, spectrum_rows, n_max=40):
    """spectrum_rows: output of lengths.length_spectrum (oriented primitive counts)."""
    t = mp.mpf(t)
    s = mp.mpf(0)
    for r in spectrum_rows:
        n_or = r["n_oriented"]
        if n_or == 0:
            continue
        l = mp.mpf(2) * mp.acosh(mp.mpf(r["halftrace"][0]) + mp.mpf(r["halftrace"][1]) * mp.sqrt(2))
        for n in range(1, n_max + 1):
            term = n_or * mp.exp(-t / 4) * l * mp.exp(-(n * l) ** 2 / (4 * t)) / (mp.sqrt(4 * mp.pi * t) * 2 * mp.sinh(n * l / 2))
            s += term
            if term < mp.mpf(10) ** (-40):
                break
    return s


def hyperbolic_tail_estimate(t, l_max):
    """ESTIMATE of omitted geodesics with l > l_max: count density ~ e^l/l (oriented: 2 e^l / l? -- we use
    the prime geodesic theorem pi(x) ~ e^x/x for unoriented-or-oriented ambiguity, multiplied by 2 to be safe)."""
    t = mp.mpf(t)
    f = lambda l: 2 * mp.exp(l) / l * l * mp.exp(-t / 4) * mp.exp(-l * l / (4 * t)) / (mp.sqrt(4 * mp.pi * t) * 2 * mp.sinh(l / 2))
    return mp.quad(f, [l_max, l_max + 5, l_max + 20, mp.inf])


def spectral_side(t, eigs=None):
    t = mp.mpf(t)
    if eigs is None:
        eigs = SP.su_eigenvalues()
    s = 1 + mp.fsum(mp.exp(-t * lam) for lam in eigs)
    tail = mp.exp(-eigs[-1] * t) / t
    return s, tail


def balance(ts, T_max=80.0, dps=30, spectrum_rows=None):
    """Return rows (t, spectral, identity, hyperbolic, residual, spectral_tail_bound, geometric_tail_est)."""
    with mp.workdps(dps):
        if spectrum_rows is None:
            spectrum_rows = LS.length_spectrum(T_max)
        l_max = 2 * math.acosh(T_max)
        eigs = SP.su_eigenvalues(dps + 10)
        out = []
        for t in ts:
            spec, stail = spectral_side(t, eigs)
            idt = identity_term(t)
            hyp = hyperbolic_term(t, spectrum_rows)
            gtail = hyperbolic_tail_estimate(t, l_max)
            out.append(dict(t=t, spectral=spec, identity=idt, hyperbolic=hyp, residual=spec - idt - hyp,
                            spectral_tail_bound=stail, geometric_tail_estimate=gtail))
        return out


def heat_constant_term(ts=(0.01, 0.02, 0.04), dps=30):
    """Theta(t) - 1/t from the S-U eigenvalues (should tend to chi/6 = -1/3 with O(t) correction; exact
    identity-term expansion: 1/t - 1/3 + t/15 + O(t^2)?) -- we report Theta(t) - identity_term(t), which is
    exponentially small (hyperbolic part), and identity_term(t) - 1/t + 1/3, which is O(t)."""
    with mp.workdps(dps):
        eigs = SP.su_eigenvalues(dps + 10)
        rows = []
        for t in ts:
            spec, stail = spectral_side(t, eigs)
            idt = identity_term(t)
            rows.append(dict(t=t, theta_minus_inv_t=spec - 1 / mp.mpf(t), identity_minus_inv_t=idt - 1 / mp.mpf(t),
                             spectral_tail_bound=stail))
        return rows


from ..provenance import provenance  # noqa: E402


@provenance("NUMERICAL", "Selberg trace formula, heat-kernel test function (Strohmaier–Uski eq. (18)); geometric side "
            "from bolzakit lengths, spectral side from Strohmaier–Uski eigenvalues")
def selberg_check(t: float, T_max: float = 80.0) -> float:
    """|geometric side − spectral side| of the heat-trace Selberg formula at time t (should be ≲ 1e-12)."""
    rows = balance([t], T_max=T_max)
    return float(abs(rows[0]["residual"]))
