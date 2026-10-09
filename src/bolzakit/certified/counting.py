"""Certified eigenvalue counts for the Bolza surface via the Selberg trace formula (linear programming).

Port of math-pursuits item 23 ``lp_count.py``. Method: Fortier Bourque--Petri, arXiv:2111.14699, Section 3
(Lemma 3.2, Cohn--Elkies strategy), with the variant "Lemma A" (window (0, a] instead of [a, b]).

Conventions (FBP (3.1)): fhat(xi) = (2pi)^(-1/2) int f(x) e^{-i xi x} dx;  r(lam)^2 = lam - 1/4;
  sum_{lam in spec} fhat(r(lam)) = 2(g-1) int_0^inf r fhat(r) tanh(pi r) dr
                                   + (2pi)^(-1/2) sum_gamma Lam(gamma) f(l(gamma)) / (2 sinh(l/2)).
Test functions: f(x) = u(x^2/s^2) exp(-x^2/(2 s^2)),  fhat(xi) = s v(s^2 xi^2) exp(-s^2 xi^2/2), with
  u(y) = sum_n c_n H_{2n}(sqrt y),  v(y) = sum_n (-1)^n c_n H_{2n}(sqrt y).

Step 1 (float, untrusted): a linear program proposes c_n.  Step 2: c_n are rounded to rationals.  Step 3
(rigorous): exact sign conditions on u, v (sympy root counting); Arb enclosure of min F on the window; exact
Gaussian moments plus rigorous Arb integration and a tail bound for the integral term; the 48 FBP geodesics
of lengths l1, l2, l3 (all with f <= 0) for the geometric term.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

import numpy as np
from scipy.integrate import quad
from scipy.optimize import linprog
import flint
from flint import acb, arb, arb_poly

from ..errors import CertificationFailed
from ._arb import ball_str, certainly_gt, certainly_lt, json_float, prec
from .certificate import Certificate, Claim, Dependency, make_producer

PRECISION_BITS = 200



def _lengths():
    """(length ball, multiplicity) of the 48 FBP-verified Bolza geodesics (FBP Sec. 7.4); needs prec set."""
    sq2 = arb(2).sqrt()
    return [(2 * (1 + sq2).acosh(), 12), (2 * (3 + 2 * sq2).acosh(), 12), (2 * (5 + 3 * sq2).acosh(), 24)]


# ---------------- float LP (proposal only) ----------------
def psi(nmax, x):
    """orthonormal Hermite functions psi_0..psi_nmax at x (float array)."""
    x = np.asarray(x, float)
    P = np.zeros((nmax + 1,) + x.shape)
    P[0] = np.pi ** -0.25 * np.exp(-x * x / 2)
    if nmax > 0:
        P[1] = np.sqrt(2) * x * P[0]
    for k in range(2, nmax + 1):
        P[k] = np.sqrt(2 / k) * x * P[k - 1] - np.sqrt((k - 1) / k) * P[k - 2]
    return P


def psi_even_at_y(N, y):
    """psi_{2n}(t), n=0..N, at t = sqrt(y) (y may be negative: even functions of t are functions of y)."""
    t = np.sqrt(complex(y))
    P = [np.pi ** -0.25 * np.exp(-t * t / 2), 0]
    P[1] = np.sqrt(2) * t * P[0]
    for k in range(2, 2 * N + 1):
        P.append(np.sqrt(2 / k) * t * P[k - 1] - np.sqrt((k - 1) / k) * P[k - 2])
    return np.array([P[2 * n].real for n in range(N + 1)])


def propose(N, s, mode, a, X, lengths_f, delta=1e-5, eta=0.02):
    idx = np.arange(0, 2 * N + 1, 2)
    sgn = (-1.0) ** (idx // 2)
    l1f = lengths_f[0][0]
    fcol = lambda x: psi(2 * N, np.asarray(x, float) / s)[idx].T
    Fcol = lambda lam: np.array([s * sgn * psi_even_at_y(N, s * s * (l - 0.25)) for l in np.atleast_1d(lam)])
    env_f = lambda x: np.abs(fcol(x)).sum(1)
    env_F = lambda lam: np.abs(Fcol(lam)).sum(1)
    Icol = np.array([2 * quad(lambda r: r * s * sgn[j] * psi(2 * N, np.array(s * r))[idx[j]] * np.tanh(np.pi * r),
                              0, 60 / s, limit=400)[0] for j in range(len(idx))])
    Gcol = sum(m * l * fcol([l])[0] / np.sinh(l / 2) for l, m in lengths_f) / np.sqrt(2 * np.pi)
    obj = Icol + Gcol - Fcol(0.0)[0]
    lo = 0.0 if mode == "A" else a
    hi = a if mode == "A" else X
    lamW = np.linspace(lo, hi, 400)
    lamP = np.linspace(hi, hi + 160 / s ** 2, 8000)
    xs = np.linspace(l1f - eta, l1f + 16 * s, 8000)
    rows = [-Fcol(lamW), -Fcol(lamP), fcol(xs)]
    rhs = [-np.ones(len(lamW)), -delta * env_F(lamP), -delta * env_f(xs)]
    res = linprog(obj, A_ub=np.vstack(rows), b_ub=np.concatenate(rhs), bounds=[(-1e6, 1e6)] * len(idx),
                  method="highs")
    if res.status != 0:
        raise RuntimeError("LP failed: " + res.message)
    return res.x, res.fun


# ---------------- exact / rigorous verification ----------------
def hermite_even_in_y(n):
    """H_{2n}(sqrt y) as an exact polynomial in y."""
    import sympy as sp
    _y, _x = sp.symbols("y"), sp.symbols("x")
    p = sp.Poly(sp.hermite(2 * n, _x), _x)
    coeffs = {}
    for (k,), c in p.terms():
        assert k % 2 == 0
        coeffs[k // 2] = sp.Rational(c)
    return sp.Poly(sum(c * _y ** k for k, c in coeffs.items()), _y)


def rationalize(dvec, N):
    """c_n = d_n / sqrt(2^{2n} (2n)! sqrt(pi)) rounded to 40 significant digits (exact rationals from here on)."""
    import mpmath as mp
    cs = []
    with mp.workdps(60):
        for n, d in enumerate(dvec):
            val = mp.mpf(float(d)) / mp.sqrt(mp.mpf(2) ** (2 * n) * mp.factorial(2 * n) * mp.sqrt(mp.pi))
            cs.append(Fraction(mp.nstr(val, 40)))
    return cs


def arb_of(q):
    q = Fraction(q)
    return arb(q.numerator) / arb(q.denominator)


def no_root_on_ray(poly, y0, want_sign):
    """Exact: poly has no real root in [y0, inf) and has sign want_sign (+1/-1) at y0.  y0 rational."""
    import sympy as sp
    _y = sp.symbols("y")
    nroots = poly.count_roots(y0, None)  # exact count of distinct real roots in [y0, inf)
    val = poly.eval(y0)
    ok = (nroots == 0) and (sp.sign(val) == want_sign)
    locs = [float(r) for r in sp.Poly(poly, _y).real_roots() if r >= y0][:8] if nroots else []
    return ok, f"real roots in [y0,inf): {nroots} {locs}; sign at y0 = {sp.sign(val)}"


def acb_poly_eval(coeffs_low_to_high, z):
    acc = acb(0)
    for c in reversed(coeffs_low_to_high):
        acc = acc * z + acb(c)
    return acc


def _coef_arb(c):
    import sympy as sp
    c = sp.Rational(c)
    return arb_of(Fraction(int(c.p), int(c.q)))


def certify(mode, a, X, N, s, dvec, target, lengths):
    """Rigorous verification (caller sets the working precision). Returns a dict of results."""
    import sympy as sp
    _y = sp.symbols("y")
    s_q = Fraction(s).limit_denominator(1000)
    # exact endpoints: the certified window is exactly the requested decimal window
    a_q = Fraction(str(a))
    X_q = Fraction(str(X)) if X is not None else None
    cs = rationalize(dvec, N)
    Hs = [hermite_even_in_y(n) for n in range(N + 1)]
    u = sum((sp.Rational(c.numerator, c.denominator) * H for c, H in zip(cs, Hs)), sp.Poly(0, _y))
    v = sum((sp.Rational((-1) ** n * c.numerator, c.denominator) * H for n, (c, H) in enumerate(zip(cs, Hs))),
            sp.Poly(0, _y))
    S = sp.Rational(s_q.numerator, s_q.denominator)
    out = {"mode": mode, "N": N, "s": str(s_q), "a": str(a_q), "X": None if X_q is None else str(X_q)}
    # (1) f <= 0 on [sys, inf): u(y) <= 0 for y >= l1^2/s^2.  Use a rational y0 <= l1^2/s^2 (safe: larger set).
    y0_arb = (lengths[0][0] ** 2) / (arb_of(s_q) ** 2)
    y0 = sp.Rational(Fraction(float(y0_arb.lower())).limit_denominator(10 ** 8)) - sp.Rational(1, 10 ** 6)
    y0_ok = certainly_lt(arb_of(Fraction(int(y0.p), int(y0.q))), y0_arb)
    ok_u, msg_u = no_root_on_ray(u, y0, -1)
    ok_u = bool(ok_u and y0_ok)
    out["u_nonpos_on_ray"] = bool(ok_u)
    out["u_msg"] = msg_u
    # (2) F >= 0 beyond the window: v(y) >= 0 for y >= s^2 (hi - 1/4)
    hi = a_q if mode == "A" else X_q
    y1 = S ** 2 * (sp.Rational(hi.numerator, hi.denominator) - sp.Rational(1, 4))
    ok_v, msg_v = no_root_on_ray(v, y1, +1)
    out["v_nonneg_beyond_window"] = bool(ok_v)
    out["v_msg"] = msg_v
    vc = [_coef_arb(c) for c in reversed(v.all_coeffs())]
    uc = [_coef_arb(c) for c in reversed(u.all_coeffs())]
    vpol = arb_poly(vc)
    upol = arb_poly(uc)
    sA = arb_of(s_q)
    F = lambda lam: sA * vpol(sA ** 2 * (lam - arb(0.25))) * (-(sA ** 2) * (lam - arb(0.25)) / 2).exp()
    # (3) c = rigorous lower bound of F on the window [lo, hi]
    lo = Fraction(0) if mode == "A" else a_q
    M = 4000
    cmin = None
    pieces_ok = True
    for i in range(M):
        l0 = lo + (hi - lo) * Fraction(i, M)
        l1_ = lo + (hi - lo) * Fraction(i + 1, M)
        ball = arb((arb_of(l0) + arb_of(l1_)) / 2, float(Fraction(l1_ - l0) / 2) * (1 + 1e-12))
        # the ball must contain the closed subinterval [l0, l1_] (checked in ball arithmetic, not assumed)
        if not (ball.lower() <= arb_of(l0).lower() and ball.upper() >= arb_of(l1_).upper()):
            pieces_ok = False  # subinterval not covered by its ball: cannot trust this piece
        lb = F(ball).lower()
        if not lb.is_finite():
            pieces_ok = False
            continue
        cmin = lb if cmin is None or lb < cmin else cmin
    if cmin is None:
        cmin = arb(0)
        pieces_ok = False
    out["c_lower"] = float(cmin)
    # (4) integral term I = 2(g-1) int_0^inf r F~(r) tanh(pi r) dr, g = 2;  F~(r) = s v(s^2 r^2) e^{-s^2 r^2/2}
    #   part 1 (tanh -> 1): exact Gaussian moments
    part1 = sum(sp.Rational(vk) * 2 ** k * sp.factorial(k) for k, vk in enumerate(reversed(v.all_coeffs()))) * 2 / S
    part1_arb = arb_of(Fraction(int(part1.p), int(part1.q)))
    #   part 2: -4 s int_0^inf r v(s^2 r^2) e^{-s^2 r^2/2} / (e^{2 pi r} + 1) dr ; rigorous on [0,R] + tail
    R = 12

    def integrand(z, analytic):
        return (4 * acb(sA) * z * acb_poly_eval(vc, acb(sA) ** 2 * z * z) * (-(acb(sA) ** 2) * z * z / 2).exp()
                / ((2 * acb.pi() * z).exp() + 1))

    part2_main = -acb.integral(integrand, 0, R, rel_tol=arb(2) ** -60).real
    absmom = sum(abs(sp.Rational(vk)) * 2 ** k * sp.factorial(k)
                 for k, vk in enumerate(reversed(v.all_coeffs()))) * 4 / S
    tail = arb_of(Fraction(int(absmom.p), int(absmom.q))) * (-2 * arb.pi() * R).exp()
    I = part1_arb + part2_main + arb(0, float(tail.upper()) * 1.0000001)
    out["I"] = ball_str(I)
    # (5) geometric term over the 48 FBP geodesics (f(l) <= 0 there, by (1))
    G = arb(0)
    for L, m in lengths:
        fL = upol(L ** 2 / sA ** 2) * (-(L ** 2) / (2 * sA ** 2)).exp()
        G += m * L * fL / (L / 2).sinh()
    G = G / (2 * arb.pi()).sqrt()
    out["G_P"] = ball_str(G)
    # (6) F(0) = fhat(i/2)
    F0 = F(arb(0))
    out["F0"] = ball_str(F0)
    bound = (I + G - F0) / cmin
    out["bound"] = ball_str(bound)
    out["bound_upper"] = float(bound.upper())
    out["target"] = target
    out["c_positive"] = bool(pieces_ok and certainly_gt(cmin, 0))
    out["bound_below_target"] = bool(out["c_positive"] and certainly_lt(bound, target))
    out["c_coeffs"] = [str(c) for c in cs]
    return out


@dataclass
class CountResult:
    bound: float
    target: int
    certified: bool
    certificate: Certificate


def count_eigenvalues(surface="bolza", window=(0, 3.3), N=21, strict=False, target=None, s=None) -> CountResult:
    """Certified upper bound on the number of Laplace eigenvalues in a window, by the Selberg trace formula.

    window=(0, a): mode A, count in (0, a]; default target 1 (bound < 1 means there are none).
    Otherwise mode B, count in [a, b]; default target 4.

    Mode A (Lemma A variant) is assumption-free. Mode B is FBP Lemma 3.2 verbatim and is conditional on the
    hypothesis lambda_1 >= a, which this function does not prove: the certificate records it as dependency
    "H_a" and in inputs["assumes"] (discharge it with a mode-A certificate for (0, a)).
    The certified window is exactly the requested decimal window (endpoints taken as exact rationals,
    recorded in inputs["window_exact"]).
    """
    if surface != "bolza":
        raise ValueError("counting needs exact length data; Bolza only in v0.2")
    a, b = float(window[0]), float(window[1])
    mode = "A" if a == 0 else "B"
    if mode == "A":
        a_win, X = b, None
    else:
        a_win, X = a, b
    if target is None:
        target = 1 if mode == "A" else 4
    if s is None:
        s = 1.0 if mode == "A" else 1.3
    with prec(PRECISION_BITS):
        lengths = _lengths()
        lengths_f = [(float(l.mid()), m) for l, m in lengths]
        d, fun = propose(N, s, mode, a_win, X, lengths_f)
        res = certify(mode, a_win, X, N, s, d, target, lengths)
        lengths_data = [[ball_str(l), m] for l, m in lengths]
    ok_u, ok_v = res["u_nonpos_on_ray"], res["v_nonneg_beyond_window"]
    bound = res["bound_upper"]
    claims = [
        Claim("u(y) on [l1^2/s^2, inf)", "<=", "0", ok_u),
        Claim("v(y) on [s^2(hi-1/4), inf)", ">=", "0", ok_v),
        Claim("min F on window", ">", "0", res["c_positive"]),
        Claim("count_upper_bound", "<", str(target), res["bound_below_target"]),
    ]
    certified = all(c.holds for c in claims)
    data = {"N": N, "s": res["s"], "a": res["a"], "X": res["X"], "c_coeffs": res["c_coeffs"],
            "lengths": lengths_data, "target": target, "mode": mode,
            "I": res["I"], "G_P": res["G_P"], "F0": res["F0"], "c_lower": json_float(res["c_lower"]),
            "bound": res["bound"], "bound_upper": json_float(bound), "u_msg": res["u_msg"], "v_msg": res["v_msg"]}
    lo_exact = "0" if mode == "A" else res["a"]
    hi_exact = res["a"] if mode == "A" else res["X"]
    inputs = {"surface": "bolza", "window": [window[0], window[1]], "N": N,
              "window_exact": [lo_exact, hi_exact]}
    deps = [Dependency("D2", "FBP trace formula (3.1), " + ("Lemma 3.2 (mode B)" if mode == "B"
                                                              else "Lemma A variant (mode A)"),
                       "Fortier Bourque-Petri, arXiv:2111.14699"),
            Dependency("D3", "Bolza geodesic length data (48 geodesics of lengths l1, l2, l3)",
                       "Fortier Bourque-Petri, arXiv:2111.14699, Sec. 7.4")]
    if mode == "B":
        inputs["assumes"] = f"lambda_1 >= {res['a']}"
        deps.append(Dependency("H_a", "lambda_1 >= a (e.g. from a mode-A count certificate for (0, a))",
                               "FBP Lemma 3.2 hypothesis"))
    cert = Certificate(
        kind="count",
        inputs=inputs,
        claims=claims, data=data,
        dependencies=deps,
        producer=make_producer())
    if strict and not certified:
        raise CertificationFailed("could not certify the eigenvalue count bound")
    return CountResult(bound=bound, target=target, certified=certified, certificate=cert)
