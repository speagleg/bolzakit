"""Length spectrum of the Bolza surface with multiplicities of primitive conjugacy classes.

Method (exact traces, float lengths, integer counts obtained by rounding a sum that must be an integer).
For a hyperbolic conjugacy class [gamma] with primitive root gamma_0, the conjugates eta gamma eta^-1 have
axes eta(axis gamma); summing over all conjugates the length of (axis cap F) gives exactly l(gamma_0)
(F = closed Dirichlet octagon; a segment lying ON a side of F is shared by two tiles and gets weight 1/2).
So for each exact half-trace t:

    S(t) := sum_{gamma in Gamma, |halftr gamma| = t} len(axis(gamma) cap F)
          = sum_{oriented classes [gamma] with |halftr| = t} l(primitive root of gamma).

Elements whose axis meets F satisfy d(0, axis) <= R_c, hence sinh(d(0,gamma 0)/2) <= cosh(R_c) sinh(l/2);
we enumerate the displacement ball with that radius (complete, see group.ball_by_displacement).
Then the number of ORIENTED primitive classes with half-trace t is

    n(t) = ( S(t) - sum_{t0, m>=2 : T_m(t0) = t} n(t0) l(t0) ) / l(t),

which must come out an integer (checked).  gamma and gamma^-1 are distinct oriented classes with the same
trace and are never conjugate in a torsion-free Fuchsian group, so UNORIENTED primitive closed geodesics
number n(t)/2.  Conventions: 'oriented' = what the Selberg trace formula sums over (S-U 2013 eq. (18)).
"""
from __future__ import annotations

import math
from collections import defaultdict

import numpy as np

from . import group as G
from . import zeta8 as Z
from ..geometry import octagon as Geo


def cheb_T(m, t):
    """T_m(t) for t = (p, q) in Z[sqrt2], exact (half-trace of gamma^m)."""
    def zmul(x, y):
        return (x[0] * y[0] + 2 * x[1] * y[1], x[0] * y[1] + x[1] * y[0])
    T0, T1 = (1, 0), t
    if m == 0:
        return T0
    for _ in range(m - 1):
        T0, T1 = T1, tuple(2 * a - b for a, b in zip(zmul(t, T1), T0))
    return T1


def length_of(t):
    return 2 * math.acosh(Z.zs_float(t))


def length_spectrum(T_max: float, return_raw=False):
    """Primitive oriented multiplicities for all half-traces t <= T_max (exact t).

    Returns list of dicts sorted by t: {'halftrace': (p,q), 'halftrace_float', 'length', 'S', 'n_oriented',
    'n_unoriented', 'primitive_power_of': [...], 'residual'} where residual = |n - round(n)|."""
    l_max = 2 * math.acosh(T_max)
    D = 2 * math.asinh(G.COSH_RC * math.sinh(l_max / 2))
    X = G.ball_by_displacement(D)
    P, Q = G.rows_trace(X)
    tr = P.astype(float) + Q.astype(float) * math.sqrt(2)
    ht = np.abs(tr) / 2
    ell = 2 * np.arccosh(np.maximum(ht, 1.0))
    ch = G.rows_cosh_disp0(X)
    sh_half = np.sqrt(np.maximum((ch - 1) / 2, 0))   # sinh(d/2)
    keep = (ht > 1 + 1e-9) & (ht <= T_max * (1 + 1e-12)) & (sh_half <= G.COSH_RC * np.sinh(ell / 2) * (1 + 1e-9))
    idx = np.nonzero(keep)[0]
    S = defaultdict(float)
    count = defaultdict(int)
    for i in idx:
        g = G.elt(X[i])
        t = G.abs_half_trace(g)
        L = Geo.axis_length_in_F(G.to_matrix(g))
        if L > 0:
            S[t] += L
            count[t] += 1
    ts = sorted(S, key=Z.zs_float)
    out = []
    n_or = {}
    for t in ts:
        lt = length_of(t)
        corr = 0.0
        powers = []
        for t0 in n_or:
            for m in range(2, 64):
                tm = cheb_T(m, t0)
                if Z.zs_cmp(tm, t) > 0:
                    break
                if tm == t:
                    corr += n_or[t0] * length_of(t0)
                    powers.append((t0, m))
        nval = (S[t] - corr) / lt
        n = int(round(nval))
        n_or[t] = n
        out.append(dict(halftrace=t, halftrace_float=Z.zs_float(t), length=lt, S=S[t],
                        n_oriented=n, n_unoriented=n / 2, power_of=powers, residual=abs(nval - n),
                        elements_meeting_F=count[t]))
    if return_raw:
        return out, dict(ball_radius=D, ball_size=len(X))
    return out


def classes_by_counting_check(T_max: float):
    """Independent sanity check of the multiplicities via the prime geodesic theorem is too weak; instead
    the Selberg trace formula (selberg.py) is the acceptance test."""
    raise NotImplementedError


from .. import data as _data  # noqa: E402
from ..errors import OutOfValidityRange  # noqa: E402
from ..provenance import provenance  # noqa: E402


@provenance("NUMERICAL", "exact half-traces; multiplicities rounded from axis sums and confirmed by the "
            "Selberg trace formula (item 02)",
            "integrality residual <= 1.1e-8 (max 1.09e-8 at half-trace 59+42*sqrt2)")
def lengths(max_half_trace: float = 120.0, compute: bool = False) -> list[dict]:
    """Primitive length spectrum: entries with half-trace p+q√2 <= max_half_trace, sorted by length."""
    bundle = _data.load("lengths")
    bound = bundle["max_half_trace"]
    if max_half_trace > bound:
        if not compute:
            raise OutOfValidityRange(
                f"lengths are bundled up to half-trace {bound:g}; pass compute=True to enumerate further "
                f"(cost grows like exp(length); ~minutes at 120)")
        rows = length_spectrum(max_half_trace)
        return [{"halftrace": tuple(r["halftrace"]), "length": r["length"], "n_oriented": r["n_oriented"],
                 "n_unoriented": r["n_oriented"] // 2} for r in rows]
    s2 = math.sqrt(2)
    return [{"halftrace": tuple(e["halftrace"]), "length": e["length"], "n_oriented": e["n_oriented"],
             "n_unoriented": e["n_unoriented"]}
            for e in bundle["entries"] if e["halftrace"][0] + e["halftrace"][1] * s2 <= max_half_trace * (1 + 1e-12)]
