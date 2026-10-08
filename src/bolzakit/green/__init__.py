"""Green's function, Robin mass and Bergman kernel of the Bolza surface."""
from .surface import green, SurfaceGreen, VALIDITY_RADIUS  # noqa: F401
from ._common import G1, A_AREA as AREA  # noqa: F401
from .constants import constants  # noqa: F401
from .bergman import bergman_at_weierstrass  # noqa: F401

import sys as _sys, types as _types  # noqa: E402


class _CallableModule(_types.ModuleType):
    """bk.green(T, D_ball) builds the Green's function; bk.green.constants() etc. still work."""
    def __call__(self, *args, **kwargs):
        return green(*args, **kwargs)


_sys.modules[__name__].__class__ = _CallableModule
