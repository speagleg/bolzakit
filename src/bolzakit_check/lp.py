"""Independent re-evaluation of a trace-formula eigenvalue count (kind "count").

Written from the formulas only (FBP arXiv:2111.14699 (3.1), Lemma 3.2 and the "Lemma A" variant; PROOF.md Lemma 1):
  F(lam) = s v(s^2 (lam - 1/4)) exp(-s^2 (lam - 1/4)/2),   f(x) = u(x^2/s^2) exp(-x^2/(2 s^2)),
  u(y) = sum_n c_n H_2n(sqrt y),  v(y) = sum_n (-1)^n c_n H_2n(sqrt y)   (physicists' Hermite, c_n and s rational).
  Trace formula (g = 2): sum_{lam in spec} F(lam) = I + (2 pi)^(-1/2) sum_gamma Lam(gamma) f(l(gamma)) / (2 sinh(l/2)),
  I = 2 int_0^inf r F~(r) tanh(pi r) dr,  F~(r) = s v(s^2 r^2) exp(-s^2 r^2/2).
If u <= 0 on [l1^2/s^2, inf) (f <= 0 on every closed geodesic, l1 = systole), F >= 0 for lam > hi and F >= c > 0 on
the window [lo, hi], then  #(spec in window) <= (I + G_P - F(0)) / c,  G_P = the geodesic sum restricted to the 48
FBP geodesics (12 of length l1, 12 of l2, 24 of l3; both orientations of each).  Mode A: window (0, a] (lo = 0).
Mode B: window [a, X], valid under the hypothesis lambda_1 >= a (so (0, a) holds no eigenvalue).

Sign conditions are decided exactly (Sturm sequences over Fraction); everything else is Arb.
"""
from __future__ import annotations

from fractions import Fraction
from math import factorial

from flint import acb, arb

from ._cert import CheckError, arb_q, ball, frac, gt, lt, precision

PREC = 200
MULTS = (12, 12, 24)
LHS = ("u(y) on [l1^2/s^2, inf)", "v(y) on [s^2(hi-1/4), inf)", "min F on window", "count_upper_bound")


# ---------------------------------------------------------------- exact polynomials (coefficient lists, low -> high)
def hermite_even(n):
    """H_2n(sqrt y) in y:  H_2n(x) = (2n)! sum_k (-1)^(n-k) (2x)^(2k) / ((2k)! (n-k)!)."""
    return [Fraction(factorial(2 * n) * (-1) ** (n - k) * 4 ** k, factorial(2 * k) * factorial(n - k))
            for k in range(n + 1)]


def _trim(p):
    p = list(p)
    while p and p[-1] == 0:
        p.pop()
    return p


def _eval(p, y):
    acc = Fraction(0)
    for c in reversed(p):
        acc = acc * y + c
    return acc


def _rem(a, b):
    a = list(a)
    while len(a) >= len(b) and a:
        q = a[-1] / b[-1]
        sh = len(a) - len(b)
        for i, c in enumerate(b):
            a[sh + i] -= q * c
        a = _trim(a)
    return a


def _sign(x):
    return (x > 0) - (x < 0)


def sturm_roots_above(p, y0):
    """Number of distinct real roots of p in (y0, inf), for p(y0) != 0 (Sturm's theorem, exact)."""
    p = _trim(p)
    seq = [p, _trim([i * c for i, c in enumerate(p)][1:])]
    while seq[-1]:
        r = _rem(seq[-2], seq[-1])
        if not r:
            break
        lc = abs(r[-1])
        seq.append([-c / lc for c in r])   # positive rescaling keeps the signs of the Sturm sequence
    seq = [s for s in seq if s]

    def changes(signs):
        signs = [s for s in signs if s != 0]
        return sum(1 for a, b in zip(signs, signs[1:]) if a != b)
    return changes([_sign(_eval(s, y0)) for s in seq]) - changes([_sign(s[-1]) for s in seq])


def no_root_on_ray(p, y0, want):
    """p has sign `want` on the closed ray [y0, inf) (exact)."""
    p = _trim(p)
    if not p:
        return False
    return _sign(_eval(p, y0)) == want and sturm_roots_above(p, y0) == 0


# ---------------------------------------------------------------- Arb evaluation
def _horner(p_arb, z):
    acc = z * 0
    for c in reversed(p_arb):
        acc = acc * z + c
    return acc


def lengths():
    """Closed-form Bolza lengths l1, l2, l3 (FBP Sec. 7.4)."""
    r2 = arb(2).sqrt()
    return [2 * (1 + r2).acosh(), 2 * (3 + 2 * r2).acosh(), 2 * (5 + 3 * r2).acosh()]


def _rational_below(x: arb) -> Fraction:
    """A rational certainly < x (finite x)."""
    q = Fraction(float(x.lower())) - Fraction(1, 10 ** 6)
    if not lt(arb_q(q), x):
        raise CheckError("could not place a rational below l1^2/s^2")
    return q


def verify(data: dict, inputs: dict, dependencies: list):
    """Re-derive the four count claims.  Returns ([(lhs, op, rhs, truth)], conditional note or None)."""
    mode = data["mode"]
    if mode not in ("A", "B"):
        raise CheckError("unknown count mode")
    N = data["N"]
    cs = [frac(c) for c in data["c_coeffs"]]
    if not isinstance(N, int) or len(cs) != N + 1:
        raise CheckError("c_coeffs length does not match N")
    s = frac(data["s"])
    a = frac(data["a"])
    if s <= 0 or a <= 0:
        raise CheckError("s and a must be positive")
    target = data["target"]
    if isinstance(target, bool) or not isinstance(target, int):
        raise CheckError("target must be an integer")
    win = inputs["window"]
    win_q = [Fraction(str(float(w))) for w in win]
    if mode == "A":
        if data["X"] is not None or win_q != [0, a]:
            raise CheckError("mode A window does not match (0, a]")
        lo, hi = Fraction(0), a
    else:
        X = frac(data["X"])
        if win_q != [a, X] or not 0 < a < X:
            raise CheckError("mode B window does not match [a, X]")
        lo, hi = a, X
        if inputs.get("assumes") != f"lambda_1 >= {data['a']}":
            raise CheckError("mode B certificate does not record its hypothesis lambda_1 >= a in inputs['assumes']")
    if inputs.get("window_exact") != [str(lo), str(hi)]:
        raise CheckError("window_exact does not match the certified window")
    Hs = [hermite_even(n) for n in range(N + 1)]
    u = [sum((c * H[k] for c, H in zip(cs, Hs) if k < len(H)), Fraction(0)) for k in range(N + 1)]
    v = [sum(((-1) ** n * c * H[k] for n, (c, H) in enumerate(zip(cs, Hs)) if k < len(H)), Fraction(0))
         for k in range(N + 1)]
    with precision(PREC):
        Ls = lengths()
        st = data["lengths"]
        if len(st) != 3 or any(int(m) != mm or not ball(l_).overlaps(L) for (l_, m), mm, L
                                in zip(st, MULTS, Ls)):
            raise CheckError("stored geodesic lengths / multiplicities differ from the FBP data")
        sA = arb_q(s)
        # (1) u < 0 on [y0, inf), y0 rational <= l1^2/s^2
        y0 = _rational_below(Ls[0] ** 2 / sA ** 2)
        ok_u = no_root_on_ray(u, y0, -1)
        # (2) v > 0 on [s^2 (hi - 1/4), inf)
        ok_v = no_root_on_ray(v, s * s * (hi - Fraction(1, 4)), +1)
        ua = [arb_q(c) for c in u]
        va = [arb_q(c) for c in v]

        def F(lam):
            y = sA ** 2 * (lam - arb(1) / 4)
            return sA * _horner(va, y) * (-y / 2).exp()
        # (3) lower bound of F on the window: union balls of 4000 exact rational pieces
        M = 4000
        c = None
        for i in range(M):
            l0 = lo + (hi - lo) * Fraction(i, M)
            l1 = lo + (hi - lo) * Fraction(i + 1, M)
            fb = F(arb_q(l0).union(arb_q(l1))).lower()
            if not fb.is_finite():
                c = arb("nan")
                break
            c = fb if c is None or fb < c else c
        c_pos = gt(c, 0)
        # (4) I = 2 int_0^inf r F~(r) dr - 4 int_0^inf r F~(r) / (e^{2 pi r} + 1) dr   (tanh = 1 - 2/(e^{2 pi r}+1))
        #     first part exactly: int_0^inf r s v(s^2 r^2) e^{-s^2 r^2/2} dr = (1/s) sum_k v_k 2^k k!
        mom = sum((vk * 2 ** k * factorial(k) for k, vk in enumerate(v)), Fraction(0))
        absmom = sum((abs(vk) * 2 ** k * factorial(k) for k, vk in enumerate(v)), Fraction(0))
        R = 12
        vac = [acb(x) for x in va]
        sC = acb(sA)

        def g(z, analytic):
            y = sC * sC * z * z
            return 4 * sC * z * _horner(vac, y) * (-y / 2).exp() / ((2 * acb.pi() * z).exp() + 1)
        main = acb.integral(g, 0, R, rel_tol=arb(2) ** -60).real
        # tail on [R, inf): 1/(e^{2 pi r}+1) <= e^{-2 pi R};  int_0^inf 4 r s |v|(s^2 r^2) e^{-s^2 r^2/2} dr
        tail = arb_q(4 * absmom / s) * (-2 * arb.pi() * R).exp()
        I = arb_q(2 * mom / s) - main + arb(0, 1) * tail
        # (5) G_P = (2 pi)^(-1/2) sum_{48 geodesics, 2 orientations} l f(l) / (2 sinh(l/2))
        G = arb(0)
        for L, m in zip(Ls, MULTS):
            fL = _horner(ua, L ** 2 / sA ** 2) * (-(L ** 2) / (2 * sA ** 2)).exp()
            G += 2 * m * L * fL / (2 * (L / 2).sinh())
        G = G / (2 * arb.pi()).sqrt()
        F0 = F(arb(0))
        for key, val in (("I", I), ("G_P", G), ("F0", F0)):
            if not (val.is_finite() and ball(data[key]).overlaps(val)):
                raise CheckError(f"stored {key} does not overlap the re-derived enclosure")
        bound = (I + G - F0) / c if c_pos else arb("nan")
        ok_b = c_pos and lt(bound, target)
    claims = [(LHS[0], "<=", "0", ok_u), (LHS[1], ">=", "0", ok_v), (LHS[2], ">", "0", c_pos),
              (LHS[3], "<", str(target), ok_b)]
    note = None
    if mode == "B":
        if not any(d.get("id") == "H_a" for d in dependencies):
            raise CheckError("mode B certificate lacks the H_a dependency")
        note = inputs["assumes"]
    return claims, note
