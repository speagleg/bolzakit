"""NaN-safe Arb ball helpers shared by the certified modules."""
from __future__ import annotations

import math
from contextlib import contextmanager
from decimal import ROUND_CEILING, Decimal, localcontext
from fractions import Fraction

import flint
from flint import arb

PM = arb(0, 1)  # the ball [0 +/- 1]


@contextmanager
def prec(bits):
    """Temporarily set the working precision (bits) of the flint context."""
    old = flint.ctx.prec
    flint.ctx.prec = bits
    try:
        yield
    finally:
        flint.ctx.prec = old


def _cmp(x, q, op):
    # Arb comparisons are certain-only: False when the ball cannot decide. Non-finite x is rejected.
    return bool(x.is_finite() and op(x, arb(q)))


def certainly_lt(x, q) -> bool:
    return _cmp(x, q, lambda a, b: a < b)


def certainly_le(x, q) -> bool:
    return _cmp(x, q, lambda a, b: a <= b)


def certainly_gt(x, q) -> bool:
    return _cmp(x, q, lambda a, b: a > b)


def ball_str(x) -> str:
    """Outward-enclosing ball string '[mid +/- rad]': 30 mid digits, and a 7-digit radius covering the true radius
    plus the decimal rounding of the mid (arb's own str rounds the radius to 3 digits, which can widen a ball by
    0.1% of its radius).  arb(ball_str(x)) contains x."""
    if not x.is_finite():
        return x.str(radius=True)
    m = x.mid()
    ms = m.str(30, radius=False, more=True)
    pm = arb(ms)
    exact = pm.rad() == 0 and pm == m                   # the 30-digit mid is exactly the arb midpoint
    r = x.rad() if exact else (x.rad() + abs(m) * arb("1e-29")).upper()   # |rounding| <= 0.5e-29 |m|
    if r == 0:
        return ms
    with localcontext() as ctx:                         # radius: 8 significant digits, rounded UP past r
        ctx.prec = 8
        ctx.rounding = ROUND_CEILING
        rd = +(Decimal(r.str(20, radius=False, more=True)) * (1 + Decimal("1e-15")))
    return f"[{ms} +/- {rd:E}]"


def json_float(x):
    """x as a float if finite, else the string 'nan' / 'inf' / '-inf' (certificates are strict JSON)."""
    x = float(x)
    if x != x:
        return "nan"
    if x in (float("inf"), float("-inf")):
        return "inf" if x > 0 else "-inf"
    return x


def attested_upper(ball: str):
    """Smallest-ish float >= the upper end of a stored ball string '[mid +/- rad]' or exact decimal string (exact
    Fraction arithmetic, as the checker reads it); None for a non-finite or unparsable string."""
    s = ball.strip()
    try:
        if s.startswith("["):
            mid, rad = s[1:-1].split("+/-")
            q = Fraction(mid.strip()) + Fraction(rad.strip())
        else:
            q = Fraction(s)
    except (ValueError, ZeroDivisionError, OverflowError):
        return None
    f = float(q)
    while Fraction(f) < q:
        f = math.nextafter(f, math.inf)
    return f


def parse_ball(s) -> arb:
    return arb(s)


def _endpoint_str(end, sign) -> str:
    # end: arb with (almost) exact value; round outward in `sign` direction using Decimal arithmetic.
    if not end.is_finite():
        raise ValueError("cannot take a rational bound of a non-finite ball")
    with localcontext() as ctx:
        ctx.prec = 80
        mid = Decimal(end.mid().str(30, radius=False))
        nudge = abs(mid) * Decimal("1e-25") + Decimal("1e-40") + Decimal(end.rad().str(5, radius=False) or 0)
        return format(mid + sign * nudge, "f")


def rational_upper(x) -> str:
    """Decimal string >= the upper end of x."""
    return _endpoint_str(x.upper(), 1)


def rational_lower(x) -> str:
    """Decimal string <= the lower end of x."""
    return _endpoint_str(x.lower(), -1)
