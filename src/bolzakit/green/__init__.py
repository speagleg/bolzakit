"""Green's function, Robin mass and Bergman kernel of the Bolza surface."""
from .surface import green, SurfaceGreen, VALIDITY_RADIUS, VALIDATED_T_D  # noqa: F401
from ._common import G1, A_AREA as AREA  # noqa: F401
from .constants import constants  # noqa: F401
from .bergman import bergman_at_weierstrass  # noqa: F401

import sys as _sys, types as _types  # noqa: E402


class _CallableModule(_types.ModuleType):
    """bk.green(T, D_ball) builds the Green's function; bk.green.constants() etc. still work."""
    def __call__(self, *args, **kwargs):
        return green(*args, **kwargs)


_sys.modules[__name__].__class__ = _CallableModule

# so that info(bk.green) and inspect.signature(bk.green) describe the callable (final review M2);
# bk.green.green(...) is the same function, spelled for type checkers
import inspect as _inspect  # noqa: E402
_sys.modules[__name__].__provenance__ = green.__provenance__
_sys.modules[__name__].__signature__ = _inspect.signature(green)
