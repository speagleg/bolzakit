"""NaN-safe Arb ball helpers shared by the certified modules."""
from __future__ import annotations

from contextlib import contextmanager
from decimal import Decimal, localcontext

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
    """Round-trippable ball string, e.g. '[3.83 +/- 1.01e-12]'."""
    return x.str(radius=True)


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
