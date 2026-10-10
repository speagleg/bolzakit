"""Provenance labels. Every public function is decorated with @provenance(label, source, error).

Labels (see VERIFICATION.md):
  EXACT      integer / algebraic arithmetic
  PROVED     exact computation plus an argument, or a Lean theorem
  CERTIFIED  rigorous enclosure (interval / ball arithmetic, exact tail bounds)
  NUMERICAL  floating point with an error estimate
  QUOTED     taken from the literature, with a reference
"""
from __future__ import annotations

import functools
import importlib
from dataclasses import dataclass
from typing import Any, NamedTuple

from .errors import OptionalDependencyMissing

LABELS = ("EXACT", "PROVED", "CERTIFIED", "NUMERICAL", "QUOTED")


@dataclass(frozen=True)
class Info:
    label: str
    source: str
    error: str | None = None


class Result(NamedTuple):
    value: Any
    label: str
    source: str
    error: str | None


_REGISTRY: dict[str, Info] = {}


def provenance(label: str, source: str, error: str | None = None):
    if label not in LABELS:
        raise ValueError(f"unknown provenance label {label!r}; expected one of {LABELS}")
    meta = Info(label, source, error)

    def deco(fn):
        @functools.wraps(fn)
        def wrapper(*args, with_provenance: bool = False, **kwargs):
            value = fn(*args, **kwargs)
            return Result(value, label, source, error) if with_provenance else value

        wrapper.__provenance__ = meta
        _REGISTRY[f"{fn.__module__}.{fn.__qualname__}"] = meta
        return wrapper

    return deco


def info(fn) -> Info:
    try:
        return fn.__provenance__
    except AttributeError:
        raise ValueError(f"{getattr(fn, '__name__', fn)!r} carries no provenance label") from None


def registry() -> dict[str, Info]:
    return dict(_REGISTRY)


def require(module_name: str, extra: str):
    try:
        return importlib.import_module(module_name)
    except ImportError as e:
        raise OptionalDependencyMissing(
            f"{module_name} is required for this function: pip install bolzakit[{extra}]") from e
