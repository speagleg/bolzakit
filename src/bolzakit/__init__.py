"""bolzakit: exact and provenance-labelled computations on the Bolza surface."""
__version__ = "0.1.0.dev0"
from . import errors, provenance, geometry, spectral, green, vortex, arith  # noqa: F401
from .exact import group as _group_module
from .exact.group import systole  # noqa: F401
from .exact.lengths import lengths  # noqa: F401
from .exact.automorphisms import automorphisms  # noqa: F401


def group():
    """The exact Bolza Fuchsian group module (generators, relator, enumeration, traces)."""
    return _group_module
