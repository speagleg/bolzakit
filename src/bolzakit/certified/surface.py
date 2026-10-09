"""Deformed Bolza surfaces for the certified layer: named families F_A, F_B and general (projected) directions.

A surface is a glued hyperbolic octagon given by its 8 vertices on the hyperboloid (Arb balls, built at 96 bits).
Directions of deformation are
  "F_A"  equilateral, alternating angles:  R_k(t) = R + t (-1)^k + sigma(t), angles k pi/4 fixed (SRC
         cert_first_order).  sigma(t) = O(t^2) restores the angle sum; it is computed rigorously here (the linear
         family R + t(-1)^k alone misses the angle sum by ~2.83 t^2, so it is not a closed surface for t != 0).
  "F_B"  equiangular, alternating sides (SRC cert_rr_p2.fb_vertices, closed form, exact in Arb).
  a 16-vector v: a direction in Klein coordinates (x, y) of V_0..V_7, projected onto the 11-dim manifold
         {constraints = 0} (SRC octagon.constraints: L_{k+4} = L_k, k = 0..3, and angle sum = 2 pi).

All Arb work is done at 96 bits (53-bit vertex balls measurably degrade downstream bounds).
"""
from __future__ import annotations

import numpy as np
from flint import arb, arb_mat

from ..errors import CertificationFailed, OutOfValidityRange
from ._arb import certainly_lt, prec
from ._hyp import Bm
from .jets import J

T_MAX = 0.05          # chart neighbourhood of Bolza in which deform() is validated
_PREC = 96
_RADII = (1e-26, 1e-24, 1e-22, 1e-20, 1e-18)   # candidate inclusion radii, tightest first


# ------------------------------------------------------------------ generic (arb / jet) octagon geometry
def _acos(x):
    v = x.v
    w = 1 - v * v
    r = w.sqrt()
    return x._f(v.acos(), -1 / r, -v / (2 * w * r))


def _tohyp(x, y):
    t = 1 / (1 - x * x - y * y).sqrt()
    return [t, t * x, t * y]


def _constraints_hyp(H):
    """SRC octagon.constraints on hyperboloid vertices (arb or J triples): [L_{k+4} - L_k (k=0..3), sum angles - 2pi].
    Returns arb if the input is arb, else J."""
    jet = isinstance(H[0][0], J)
    H = [[J.lift(c) for c in V] for V in H]
    L = [Bm(H[k], H[(k + 1) % 8]).acosh() for k in range(8)]
    ang = J(0)
    for k in range(8):
        a, p, q = H[k], H[(k - 1) % 8], H[(k + 1) % 8]
        bp, bq = Bm(a, p), Bm(a, q)
        w1 = [p[i] - bp * a[i] for i in range(3)]
        w2 = [q[i] - bq * a[i] for i in range(3)]
        ang = ang + _acos(-Bm(w1, w2) / (Bm(w1, w1) * Bm(w2, w2)).sqrt())
    out = [L[k + 4] - L[k] for k in range(4)] + [ang - 2 * arb.pi()]
    return out if jet else [c.v for c in out]


def _constraints_klein(xs):
    """constraints of 16 Klein coordinates (arb or J)."""
    return _constraints_hyp([_tohyp(xs[2 * k], xs[2 * k + 1]) for k in range(8)])


def _jacobian(xk, dirs=None):
    """5 x m Arb matrix of directional derivatives of the constraints at Klein point xk (arb balls) along the
    columns of dirs (16 x m arb_mat; default the identity, i.e. the full 5 x 16 Jacobian).  Ball inputs give an
    enclosure of the derivative over the ball."""
    m = 16 if dirs is None else dirs.ncols()
    cols = []
    for j in range(m):
        d = [arb(1 if i == j else 0) for i in range(16)] if dirs is None else [dirs[i, j] for i in range(16)]
        cols.append([c.d for c in _constraints_klein([J(xk[i], d[i]) for i in range(16)])])
    return arb_mat([[cols[j][r] for j in range(m)] for r in range(5)])


def _min_norm(Jm, b):
    """min-norm solution Jm^T (Jm Jm^T)^-1 b of Jm w = b (Jm: 5 x 16 arb_mat, b: list of 5)."""
    Jt = Jm.transpose()
    lam = (Jm * Jt).solve(arb_mat([[x] for x in b]))
    w = Jt * lam
    return [w[i, 0] for i in range(16)]


def _mid(x):
    return arb(x.mid())


# ------------------------------------------------------------------ named families
def _fb_vertices(p0=0):
    """SRC cert_rr_p2.fb_vertices: family F_B as order-2 jets in p = phi - pi/4 at p = p0 (exact, Arb).
    Vertices at common radius R(p); V_{2j} at angle j pi/2, V_{2j+1} at angle phi + j pi/2; cosh R = C solves
    C^2 x y - C (x + y) - 1 = 0, x = tan(phi/2), y = tan(pi/4 - phi/2)  (angle at V_0 = pi/4)."""
    pi8 = arb.pi() / 8
    p0 = arb(p0)
    x = J(pi8 + p0 / 2, arb(1) / 2).tan()
    y = J(pi8 - p0 / 2, arb(-1) / 2).tan()
    s = x + y; pr = x * y
    C = (s + (s * s + 4 * pr).sqrt()) / (2 * pr)
    S = (C * C - 1).sqrt()
    out = []
    for k in range(8):
        j = k // 2
        th = J(arb.pi() * j / 2) if k % 2 == 0 else J(arb.pi() / 4 + p0 + arb.pi() * j / 2, arb(1))
        out.append([C, S * th.cos(), S * th.sin()])
    return out


def _bolza_R():
    return (3 + 2 * arb(2).sqrt()).acosh()


def _fa_hyp(Rk):
    """F_A vertices from radii R_k (arb or J), angles k pi/4."""
    out = []
    for k in range(8):
        th = arb.pi() * k / 4
        out.append([Rk[k].cosh(), Rk[k].sinh() * th.cos(), Rk[k].sinh() * th.sin()])
    return out


def _fa_anglesum(t, sig, dsig):
    """angle-sum residual (J in sigma, derivative seed dsig) of R_k = R + t(-1)^k + sigma."""
    R = _bolza_R()
    return _constraints_hyp(_fa_hyp([J(R + t * (-1) ** k + sig, dsig) for k in range(8)]))[4]


def _fa_sigma(t):
    """Certified enclosure of the unique sigma near 0 with angle sum 2 pi (1-D interval Newton)."""
    m = arb(0)
    for _ in range(30):
        f = _fa_anglesum(t, m, arb(1))
        step = f.v / f.d
        m = _mid(m - step)
        if abs(float(step.mid())) < 1e-26:
            break
    f0 = _fa_anglesum(t, m, arb(0)).v
    for r in _RADII:   # smallest inclusion box that verifies: N(box) inside box => unique zero in box
        box = arb(m.mid(), r)
        N = m - f0 / _fa_anglesum(t, box, arb(1)).d
        if certainly_lt(abs(N - m), r):
            return N
    raise OutOfValidityRange(f"F_A closure could not be certified at t={t} (validated for |t| <= {T_MAX})")


def _fa_jets():
    """F_A order-2 vertex jets at t = 0: R_k = R + t(-1)^k + sigma2 t^2, sigma2 = -(S_tt/2) / S_sigma."""
    R = _bolza_R()
    s_tt_half = _constraints_hyp(_fa_hyp([J(R, arb((-1) ** k)) for k in range(8)]))[4].e
    s_sig = _constraints_hyp(_fa_hyp([J(R, arb(1)) for k in range(8)]))[4].d
    sig2 = -s_tt_half / s_sig
    return _fa_hyp([J(R, arb((-1) ** k), sig2) for k in range(8)])


# ------------------------------------------------------------------ Surface
class Surface:
    """A closed hyperbolic genus-2 surface glued from an octagon (Bolza side pairing pattern).

    vertices_hyp: 8 hyperboloid vertices [arb, arb, arb].  r2_invariant: the octagon is invariant under rotation
    by pi/2 (V_k -> V_{k+2}); certify_spectrum then computes only 2 of the 8 sectors."""

    def __init__(self, vertices_hyp, is_bolza=False, label="custom", r2_invariant=False):
        self.trace_coords = None   # the tau of Surface.from_trace_coords, else None
        self.vertices_hyp = vertices_hyp
        self.is_bolza = is_bolza
        self.label = label
        self.r2_invariant = r2_invariant

    @classmethod
    def bolza(cls) -> "Surface":
        from .mesh import bolza_vertices_hyp
        # Vertex balls must be built at the certification precision (96 bits); default-precision balls
        # (53 bits) would cap the whole certified computation at ~1e-8 absolute radius.
        with prec(_PREC):
            return cls(bolza_vertices_hyp(), is_bolza=True, label="bolza", r2_invariant=True)

    @classmethod
    def from_trace_coords(cls, tau, check=True) -> "Surface":
        """The surface with Teichmueller trace coordinates tau (6 floats, chart.WORDS order) near Bolza.

        The octagon is the unique zero, certified by a Krawczyk inclusion at 128 bits, of [constraints; traces - tau]
        in 11 Klein coordinates (the other 5 fixed at their Bolza values; chart.solve_trace_coords).  Raises
        OutOfValidityRange if |tau - tau_B|_inf > 0.05 max(tau_B) or the inclusion fails.  The result is never
        flagged is_bolza, even for tau = chart.bolza_traces(): the floats are not the exact (irrational) Bolza traces,
        and their preimage is ~6e-17 from Bolza.  tau is kept in .trace_coords and in the label (certificates record
        it).  check=True also certifies rank d(pi) = 6 at the solution (CertificationFailed otherwise)."""
        from . import chart
        tau = [float(t) for t in tau]
        xs = chart.solve_trace_coords(tau)
        with prec(chart._PREC):
            V = [_tohyp(xs[2 * k], xs[2 * k + 1]) for k in range(8)]
        s = cls(V, label=f"trace_coords(tau={tau!r})", r2_invariant=False)
        s.trace_coords = tau
        if check and not chart.certify_rank(s).certified:
            raise CertificationFailed(f"rank d(pi) = 6 could not be certified at tau={tau!r}")
        return s

    # -------------------------------------------------------------- coordinates
    def _klein_arb(self):
        return [c for V in self.vertices_hyp for c in (V[1] / V[0], V[2] / V[0])]

    def klein_coords(self) -> list:
        """16 floats: Klein (x, y) of V_0..V_7."""
        with prec(_PREC):
            return [float(c.mid()) for c in self._klein_arb()]

    @staticmethod
    def _direction(direction):
        if isinstance(direction, str):
            if direction not in ("F_A", "F_B"):
                raise ValueError(f"unknown family {direction!r}; expected 'F_A', 'F_B' or a length-16 vector")
            return direction
        v = np.asarray(direction, dtype=float)
        if v.shape != (16,) or not np.all(np.isfinite(v)):
            raise ValueError("a vector direction must be 16 finite floats (Klein x, y of V_0..V_7)")
        return v

    def _require_bolza(self, name):
        if not self.is_bolza:
            raise ValueError(f"the named family {name} is parametrised from the Bolza surface only")

    # -------------------------------------------------------------- deformation
    def deform(self, direction, t) -> "Surface":
        """A new Surface: the family evaluated at (the exact binary value of) t, vertices as 96-bit Arb balls.

        Named families start at Bolza.  A vector v is applied in Klein coordinates, x = x_self + t v, then
        projected onto {constraints = 0} by min-norm Newton (SRC octagon.project) and certified by a Krawczyk
        inclusion in the 5-dim normal space: the returned balls contain an exactly closed octagon."""
        d = self._direction(direction)
        if not abs(float(t)) <= T_MAX:
            raise OutOfValidityRange(f"|t| = {abs(float(t))} exceeds the validated chart neighbourhood |t| <= {T_MAX}")
        with prec(_PREC):
            ta = arb(float(t))
            if isinstance(d, str):
                self._require_bolza(d)
                if d == "F_B":
                    V = [[c.v for c in W] for W in _fb_vertices(ta)]
                else:
                    sig = _fa_sigma(ta)
                    R = _bolza_R()
                    V = _fa_hyp([R + ta * (-1) ** k + sig for k in range(8)])
                return Surface(V, label=f"{self.label}+{d}(t={float(t)!r})", r2_invariant=True)
            V = self._project([x + ta * arb(vi) for x, vi in zip(self._klein_arb(), d)], float(t))
            return Surface(V, label=f"{self.label}+vector(t={float(t)!r})", r2_invariant=False)

    @staticmethod
    def _project(x, t):
        fail = OutOfValidityRange(f"projection onto the constraint manifold could not be certified at t={t} "
                                  f"(validated for |t| <= {T_MAX} with all vertices inside the Klein disk)")
        try:
            return Surface._project_certified([_mid(c) for c in x], fail)
        except ZeroDivisionError:   # singular / non-finite Jacobian, e.g. a vertex outside the Klein disk
            raise fail from None

    @staticmethod
    def _project_certified(x, fail):
        if any(not certainly_lt(x[2 * k] ** 2 + x[2 * k + 1] ** 2, 1) for k in range(8)):
            raise fail
        for _ in range(40):
            c = _constraints_klein(x)
            if max(abs(float(ci.mid())) for ci in c) < 1e-25:
                break
            dx = _min_norm(_jacobian(x), c)
            x = [_mid(xi - di) for xi, di in zip(x, dx)]
        # Krawczyk on g(lam) = C(x + A lam), A = mid(Jacobian)^T (16 x 5, exact): K(Lam) in int(Lam) => unique zero
        Jm = _jacobian(x)
        A = arb_mat([[_mid(Jm[r, i]) for r in range(5)] for i in range(16)])
        g0 = _constraints_klein(x)
        Gm = _jacobian(x, A)
        Ginv = arb_mat([[_mid(Gm[i, j]) for j in range(5)] for i in range(5)]).inv()
        Y = arb_mat([[_mid(Ginv[i, j]) for j in range(5)] for i in range(5)])
        yg = Y * arb_mat([[c] for c in g0])
        I5 = arb_mat(5, 5, [1 if i == j else 0 for i in range(5) for j in range(5)])
        for r in _RADII:   # smallest box Lam that verifies K(Lam) in int(Lam)
            Lam = arb_mat([[arb(0, r)] for _ in range(5)])
            AL = A * Lam
            X = [x[i] + AL[i, 0] for i in range(16)]
            K = -yg + (I5 - Y * _jacobian(X, A)) * Lam
            if all(certainly_lt(abs(K[i, 0]), r) for i in range(5)):
                return [_tohyp(X[2 * k], X[2 * k + 1]) for k in range(8)]
        raise fail

    # -------------------------------------------------------------- jets
    def vertex_jets(self, direction) -> list:
        """Order-2 jets J (value, d/dt, half d2/dt2) of the 8 hyperboloid vertices at t = 0 along direction.

        Vector v (Klein coordinates): the tangent is the orthogonal projection v_T of v onto ker(DC) (DC = constraint
        Jacobian).  The second-order term is the min-norm w with DC w = -D2C(v_T, v_T), so the Klein curve
        x(t) = x0 + t v_T + t^2 w / 2 satisfies C(x(t)) = O(t^3).  Any C^2 curve in the manifold with tangent v_T is
        valid downstream (dLambda = 0 at Bolza, so second-derivative bounds depend on the tangent only); this is the
        simplest rigorous choice.  D2C(v_T, v_T) = 2 * (t^2 coefficient of C(x0 + t v_T)), computed with J jets."""
        d = self._direction(direction)
        with prec(_PREC):
            if isinstance(d, str):
                self._require_bolza(d)
                return _fb_vertices(0) if d == "F_B" else _fa_jets()
            x0 = self._klein_arb()
            Jm = _jacobian(x0)
            v = [arb(vi) for vi in d]
            Jv = [sum((Jm[r, i] * v[i] for i in range(16)), arb(0)) for r in range(5)]
            n = _min_norm(Jm, Jv)
            vt = [v[i] - n[i] for i in range(16)]
            c2 = [2 * c.e for c in _constraints_klein([J(x0[i], vt[i]) for i in range(16)])]
            w = _min_norm(Jm, [-c for c in c2])
            xs = [J(x0[i], vt[i], w[i] / 2) for i in range(16)]
            return [_tohyp(xs[2 * k], xs[2 * k + 1]) for k in range(8)]

    def klein_tangent(self, direction) -> np.ndarray:
        """16 floats: first-order Taylor coefficient of Klein (x, y) = (X1/X0, X2/X0) of the vertex jets."""
        jets = self.vertex_jets(direction)
        with prec(_PREC):
            return np.array([float((V[i] / V[0]).d.mid()) for V in jets for i in (1, 2)])
