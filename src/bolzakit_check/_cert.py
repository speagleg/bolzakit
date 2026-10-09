"""Minimal certificate loading for the checker: strict JSON, the sha256 rule, and Arb parsing helpers.

Deliberately self-contained (stdlib + python-flint): it re-states the bolzakit-cert/1 hashing rule instead of importing
it, so the checker never trusts producer code.
"""
from __future__ import annotations

import hashlib
import json
from contextlib import contextmanager
from fractions import Fraction
from pathlib import Path

import flint
from flint import arb

SCHEMA = "bolzakit-cert/1"


class CheckError(Exception):
    """The certificate is malformed or a re-derivation contradicts it (message = reason)."""


def _no_constant(name):
    raise CheckError(f"non-strict JSON constant {name}")


def load(path) -> dict:
    """Parse a certificate file as strict JSON (NaN / Infinity are rejected)."""
    raw = json.loads(Path(path).read_text(encoding="utf-8"), parse_constant=_no_constant)
    if not isinstance(raw, dict):
        raise CheckError("top level is not a JSON object")
    return raw


def core_hash(raw: dict) -> str:
    """sha256 of the canonical core {schema, kind, inputs, data} (sorted keys, compact separators, UTF-8)."""
    core = {k: raw[k] for k in ("schema", "kind", "inputs", "data")}
    s = json.dumps(core, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(s.encode("utf-8")).hexdigest()


@contextmanager
def precision(bits):
    old = flint.ctx.prec
    flint.ctx.prec = bits
    try:
        yield
    finally:
        flint.ctx.prec = old


def ball(s) -> arb:
    """An Arb ball from its stored string (non-strings are malformed)."""
    if not isinstance(s, str):
        raise CheckError(f"expected a ball string, got {type(s).__name__}")
    try:
        return arb(s)
    except (ValueError, TypeError) as e:
        raise CheckError(f"unparsable ball string {s!r}") from e


def matrix(rows, n=None) -> list:
    """Square matrix of balls from a list of lists of ball strings (size n if given)."""
    if not isinstance(rows, list) or (n is not None and len(rows) != n):
        raise CheckError("matrix has the wrong size")
    out = [[ball(x) for x in r] for r in rows]
    if any(len(r) != len(out) for r in out):
        raise CheckError("matrix is not square")
    return out


def frac(s) -> Fraction:
    """Exact rational from a decimal or p/q string (or an int)."""
    if isinstance(s, bool) or not isinstance(s, (str, int)):
        raise CheckError(f"expected an exact rational, got {s!r}")
    try:
        return Fraction(s)
    except (ValueError, ZeroDivisionError) as e:
        raise CheckError(f"not an exact rational: {s!r}") from e


def arb_q(q) -> arb:
    """Ball enclosing the rational q (exact when the denominator is a power of two)."""
    q = Fraction(q)
    return arb(q.numerator) / arb(q.denominator)


def lt(x: arb, y) -> bool:
    """x < y certainly (False for non-finite x)."""
    return bool(x.is_finite() and x < y)


def gt(x: arb, y) -> bool:
    return bool(x.is_finite() and x > y)
