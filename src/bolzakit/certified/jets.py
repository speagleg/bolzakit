"""Order-2 jets (J) and order-1 duals (D) over arb; ported verbatim from the Bolza lambda1 research code
(cert_rr_p2.py and cert_first_order.py)."""
from flint import arb


# ------------------------------------------------------------------ order-2 jets over arb
class J:
    """Taylor jet  v + d p + e p^2  (e = half the second derivative)."""
    __slots__ = ("v", "d", "e")

    def __init__(self, v, d=None, e=None):
        self.v = v if isinstance(v, arb) else arb(v)
        self.d = arb(0) if d is None else d
        self.e = arb(0) if e is None else e

    @staticmethod
    def c(x):
        return x if isinstance(x, J) else J(x)

    def __add__(s, o):
        o = J.c(o); return J(s.v + o.v, s.d + o.d, s.e + o.e)
    __radd__ = __add__

    def __sub__(s, o):
        o = J.c(o); return J(s.v - o.v, s.d - o.d, s.e - o.e)

    def __rsub__(s, o):
        return J.c(o) - s

    def __neg__(s):
        return J(-s.v, -s.d, -s.e)

    def __mul__(s, o):
        if not isinstance(o, J):
            o = arb(o) if not isinstance(o, arb) else o
            return J(s.v * o, s.d * o, s.e * o)
        return J(s.v * o.v, s.v * o.d + s.d * o.v, s.v * o.e + s.d * o.d + s.e * o.v)
    __rmul__ = __mul__

    def __truediv__(s, o):
        if not isinstance(o, J):
            o = arb(o) if not isinstance(o, arb) else o
            return J(s.v / o, s.d / o, s.e / o)
        q0 = s.v / o.v
        q1 = (s.d - q0 * o.d) / o.v
        q2 = (s.e - q0 * o.e - q1 * o.d) / o.v
        return J(q0, q1, q2)

    def __rtruediv__(s, o):
        return J.c(o) / s

    def _f(s, f0, f1, f2h):
        """compose with f: f0 = f(v), f1 = f'(v), f2h = f''(v)/2."""
        return J(f0, f1 * s.d, f1 * s.e + f2h * s.d * s.d)

    def sqrt(s):
        r = s.v.sqrt(); return s._f(r, 1 / (2 * r), -1 / (8 * r * s.v))

    def sinh(s):
        a, b = s.v.sinh(), s.v.cosh(); return s._f(a, b, a / 2)

    def cosh(s):
        a, b = s.v.cosh(), s.v.sinh(); return s._f(a, b, a / 2)

    def acosh(s):
        w = (s.v * s.v - 1); r = w.sqrt()
        return s._f(s.v.acosh(), 1 / r, -s.v / (2 * w * r))

    def tan(s):
        t = s.v.tan(); u = 1 + t * t; return s._f(t, u, t * u)

    def cos(s):
        c, sn = s.v.cos(), s.v.sin(); return s._f(c, -sn, -c / 2)

    def sin(s):
        c, sn = s.v.cos(), s.v.sin(); return s._f(sn, c, -sn / 2)




class D:
    """dual number over arb: value + eps * derivative."""
    __slots__ = ("v", "d")

    def __init__(self, v, d=None):
        self.v = v if isinstance(v, arb) else arb(v)
        self.d = d if d is not None else arb(0)

    def __add__(s, o):
        o = o if isinstance(o, D) else D(o); return D(s.v + o.v, s.d + o.d)
    __radd__ = __add__

    def __sub__(s, o):
        o = o if isinstance(o, D) else D(o); return D(s.v - o.v, s.d - o.d)

    def __rsub__(s, o):
        return D(o) - s

    def __neg__(s):
        return D(-s.v, -s.d)

    def __mul__(s, o):
        o = o if isinstance(o, D) else D(o); return D(s.v * o.v, s.v * o.d + s.d * o.v)
    __rmul__ = __mul__

    def __truediv__(s, o):
        o = o if isinstance(o, D) else D(o); return D(s.v / o.v, (s.d * o.v - s.v * o.d) / (o.v * o.v))

    def __rtruediv__(s, o):
        return D(o) / s

    def sqrt(s):
        r = s.v.sqrt(); return D(r, s.d / (2 * r))

    def acosh(s):
        return D(s.v.acosh(), s.d / (s.v * s.v - 1).sqrt())

    def sinh(s):
        return D(s.v.sinh(), s.d * s.v.cosh())

    def cosh(s):
        return D(s.v.cosh(), s.d * s.v.sinh())


