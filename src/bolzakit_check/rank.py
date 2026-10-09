"""Re-derivation of a `rank` certificate from its inputs alone (the 16 Klein balls and the trace words).

Mathematics (written for the checker, independent of the producer's code):
  * Klein point (x, y) -> hyperboloid point H = (1, x, y) / sqrt(1 - x^2 - y^2), B(u, v) = u0 v0 - u1 v1 - u2 v2.
  * Octagon constraints c(x) (5 values): side_{k+4} - side_k (k = 0..3), side_k = acosh B(H_k, H_{k+1}), and
    (sum of the 8 interior angles) - 2 pi, the angle at H_k between the geodesics to H_{k-1} and H_{k+1}.
  * Side pairing g_k (k = 0..3): the orientation-preserving isometry of the hyperboloid (SO+(2,1)) with
    H_{k+4} -> H_{k+1} and the direction toward H_{k+5} -> the direction toward H_k.  With the positively oriented
    Lorentz-orthonormal frame F(P, Q) = [P | e | n] (e the unit tangent at P toward Q, n completing it, det F = +1),
    g_k = F(H_{k+1}, H_k) F(H_{k+4}, H_{k+5})^-1, where F^-1 = J F^T J, J = diag(1, -1, -1).
  * Trace coordinates: tau_w = tr(w(g)) for the words w ("g2" = g_2, "G2" = g_2^-1, products left to right).
Forward-mode dual numbers over Arb give the Jacobians Jc (5 x 16) and Jf (6 x 16).  With the pivot columns p of Jc
(det Jc[:, p] != 0), the kernel basis Tb (16 x 11) has x_free = e_j, x_p = -Jc[:, p]^-1 Jc[:, free] e_j.  Claims:
  det (Jf Tb)[:, cols] != 0;  det [Jc; Jf][:, stacked_cols] != 0;  every entry of Jc Tb encloses 0 (finite);
  c(x) encloses 0 (finite).
"""
from __future__ import annotations

import re

from flint import arb, arb_mat

from ._cert import CheckError
from .hermitian import bareiss_det


class D:
    """First-order dual number v + d eps over Arb."""
    __slots__ = ("v", "d")

    def __init__(self, v, d=0):
        self.v, self.d = arb(v) if not isinstance(v, arb) else v, arb(d) if not isinstance(d, arb) else d

    @staticmethod
    def lift(x):
        return x if isinstance(x, D) else D(x, 0)

    def __add__(self, o):
        o = D.lift(o); return D(self.v + o.v, self.d + o.d)
    __radd__ = __add__

    def __sub__(self, o):
        o = D.lift(o); return D(self.v - o.v, self.d - o.d)

    def __rsub__(self, o):
        return D.lift(o) - self

    def __neg__(self):
        return D(-self.v, -self.d)

    def __mul__(self, o):
        o = D.lift(o); return D(self.v * o.v, self.v * o.d + self.d * o.v)
    __rmul__ = __mul__

    def __truediv__(self, o):
        o = D.lift(o); return D(self.v / o.v, (self.d * o.v - self.v * o.d) / (o.v * o.v))

    def __rtruediv__(self, o):
        return D.lift(o) / self


def _sqrt(x):
    if isinstance(x, D):
        r = x.v.sqrt(); return D(r, x.d / (2 * r))
    return x.sqrt()


def _acosh(x):
    if isinstance(x, D):
        return D(x.v.acosh(), x.d / (x.v * x.v - 1).sqrt())
    return x.acosh()


def _acos(x):
    if isinstance(x, D):
        return D(x.v.acos(), -x.d / (1 - x.v * x.v).sqrt())
    return x.acos()


def _val(x):
    return x.v if isinstance(x, D) else x


def B(u, v):
    return u[0] * v[0] - u[1] * v[1] - u[2] * v[2]


def hyperboloid(xs):
    out = []
    for k in range(8):
        x, y = xs[2 * k], xs[2 * k + 1]
        t = 1 / _sqrt(1 - x * x - y * y)
        out.append([t, t * x, t * y])
    return out


def constraints(xs):
    H = hyperboloid(xs)
    side = [_acosh(B(H[k], H[(k + 1) % 8])) for k in range(8)]
    total = 0
    for k in range(8):
        a, p, q = H[k], H[(k - 1) % 8], H[(k + 1) % 8]
        t1 = [p[i] - B(a, p) * a[i] for i in range(3)]      # tangents at a (spacelike: B < 0)
        t2 = [q[i] - B(a, q) * a[i] for i in range(3)]
        total = total + _acos(-B(t1, t2) / _sqrt(B(t1, t1) * B(t2, t2)))
    return [side[k + 4] - side[k] for k in range(4)] + [total - 2 * arb.pi()]


def _det3(F):
    return (F[0][0] * (F[1][1] * F[2][2] - F[1][2] * F[2][1]) - F[0][1] * (F[1][0] * F[2][2] - F[1][2] * F[2][0])
            + F[0][2] * (F[1][0] * F[2][1] - F[1][1] * F[2][0]))


def frame(P, Q):
    """Columns [P, e, n]: Lorentz-orthonormal, det = +1 (CheckError if the orientation is not certain)."""
    e = [Q[i] - B(P, Q) * P[i] for i in range(3)]
    s = _sqrt(-B(e, e))
    e = [c / s for c in e]
    # n: Lorentz-orthogonal to P and e.  J(P x e) is B-orthogonal to both (B(u, J w) = u . w).
    cr = [P[1] * e[2] - P[2] * e[1], P[2] * e[0] - P[0] * e[2], P[0] * e[1] - P[1] * e[0]]
    n = [cr[0], -cr[1], -cr[2]]
    s = _sqrt(-B(n, n))
    n = [c / s for c in n]
    F = [[P[i], e[i], n[i]] for i in range(3)]
    dv = _val(_det3(F))
    if not dv.is_finite() or dv.contains(0):
        raise CheckError("pairing frame orientation is not certain")
    if dv < 0:
        for i in range(3):
            F[i][2] = -F[i][2]
    return F


def _mm(A, Bm):
    return [[A[i][0] * Bm[0][j] + A[i][1] * Bm[1][j] + A[i][2] * Bm[2][j] for j in range(3)] for i in range(3)]


SIG = (1, -1, -1)


def _linv(F):
    """Inverse of a Lorentz-orthonormal (or SO(2,1)) matrix: J F^T J."""
    return [[SIG[i] * F[j][i] * SIG[j] for j in range(3)] for i in range(3)]


def pairings(xs):
    H = hyperboloid(xs)
    return [_mm(frame(H[k + 1], H[k]), _linv(frame(H[k + 4], H[(k + 5) % 8]))) for k in range(4)]


_TOKEN = re.compile(r"([gG])([0-3])")


def parse_word(w):
    if not isinstance(w, str) or not w or _TOKEN.sub("", w):
        raise CheckError(f"bad trace word {w!r}")
    return [(int(i), c == "g") for c, i in _TOKEN.findall(w)]


def traces(xs, words):
    g = pairings(xs)
    out = []
    for w in words:
        M = None
        for i, fwd in w:
            X = g[i] if fwd else _linv(g[i])
            M = X if M is None else _mm(M, X)
        out.append(M[0][0] + M[1][1] + M[2][2])
    return out


def jacobians(x0, words):
    Jc = [[None] * 16 for _ in range(5)]
    Jf = [[None] * 16 for _ in range(len(words))]
    cval = None
    for j in range(16):
        xs = [D(x0[i], 1 if i == j else 0) for i in range(16)]
        c = constraints(xs)
        f = traces(xs, words)
        cval = [ci.v for ci in c]
        for i in range(5):
            Jc[i][j] = c[i].d
        for i in range(len(words)):
            Jf[i][j] = f[i].d
    return Jc, Jf, cval


def _finite0(v):
    return bool(v.is_finite() and v.contains(0))


def _nonzero(d):
    return bool(d.is_finite() and not d.contains(0))


def pairing_selfcheck(xs):
    """Independent validation of the pairing construction at the Klein balls xs: each g_k is in SO(2,1)
    (g^T J g contains J, det g contains 1) and maps H_{k+4} -> H_{k+1}, H_{k+5} -> H_k (enclosures contain)."""
    H = hyperboloid(xs)
    for k, g in enumerate(pairings(xs)):
        for i in range(3):
            for j in range(3):
                gtjg = sum((g[a][i] * SIG[a] * g[a][j] for a in range(3)), arb(0))
                if not (gtjg.is_finite() and gtjg.contains(SIG[i] if i == j else 0)):
                    raise CheckError(f"side pairing g{k} is not Lorentzian (g^T J g != J)")
        d = _det3(g)
        if not (d.is_finite() and d.contains(1)):
            raise CheckError(f"side pairing g{k} does not have determinant 1")
        for src, dst in ((H[k + 4], H[k + 1]), (H[(k + 5) % 8], H[k])):
            img = [sum((g[i][a] * src[a] for a in range(3)), arb(0)) for i in range(3)]
            if not all(v.is_finite() and v.overlaps(w) for v, w in zip(img, dst)):
                raise CheckError(f"side pairing g{k} does not map the octagon side as required")


def _int_cols(cols, n, size):
    return (isinstance(cols, list) and len(cols) == n and len(set(cols)) == n
            and all(isinstance(c, int) and not isinstance(c, bool) and 0 <= c < size for c in cols))


CANONICAL_WORDS = ["g0g1", "g0G2", "g1g2", "g1g3", "g1G3", "g2g3"]


def bolza_vertices():
    """The exact Bolza octagon in the Klein model (computed here, independently of the producer): the regular octagon
    with vertex k at R (cos k pi/4, sin k pi/4), R = tanh(d), cosh d = 3 + 2 sqrt 2, i.e. R = sqrt(c^2 - 1) / c.
    The angular factors are exact (0, +-1, +-sqrt(2)/2), so axis coordinates are the exact 0."""
    c = 3 + 2 * arb(2).sqrt()
    R = (c * c - 1).sqrt() / c
    h = arb(2).sqrt() / 2
    trig = [(1, 0), (h, h), (0, 1), (-h, h), (-1, 0), (-h, -h), (0, -1), (h, -h)]
    return [R * t for cs in trig for t in cs]


def is_bolza(x0) -> bool:
    """Every stored Klein ball contains the corresponding exact Bolza vertex coordinate."""
    return len(x0) == 16 and all(x.is_finite() and x.contains(v) for x, v in zip(x0, bolza_vertices()))


def describe(inputs, ball):
    """Message tail.  The label "bolza" is reported only when the stored Klein balls contain the exact Bolza octagon
    (checked here); any other label is echoed as unchecked."""
    w = inputs["words"]
    tag = "" if w == CANONICAL_WORDS else "non-canonical words "
    label = inputs.get("surface")
    if label == "bolza":
        if is_bolza([ball(x) for x in inputs["klein"]]):
            surf = "surface bolza (Klein vertices checked against the exact octagon)"
        else:
            surf = "surface custom (label 'bolza' rejected: the Klein balls do not contain the exact Bolza octagon)"
    else:
        surf = f"surface {label} (label unchecked)"
    return f"rank 6 for {tag}{','.join(w)}, {surf}"


def derive(inputs, data, ball):
    """[(lhs, op, rhs, truth)] for the four rank claims; stored J / dets / residual only cross-checked."""
    x0 = [ball(x) for x in inputs["klein"]]
    if len(x0) != 16:
        raise CheckError("inputs.klein must hold 16 balls")
    words = [parse_word(w) for w in inputs["words"]]
    if len(words) != 6:
        raise CheckError("expected 6 trace words")
    piv = inputs["pivot_columns_Jc"]
    if not _int_cols(piv, 5, 16):
        raise CheckError("bad pivot_columns_Jc")
    pairing_selfcheck(x0)
    Jc, Jf, cval = jacobians(x0, words)
    ok_con = all(_finite0(v) for v in cval)
    free = [j for j in range(16) if j not in piv]
    Ap = [[Jc[i][c] for c in piv] for i in range(5)]
    if not _nonzero(bareiss_det(Ap)):
        raise CheckError("pivot minor of Jc does not certainly exclude 0")
    X = arb_mat(Ap).solve(arb_mat([[Jc[i][c] for c in free] for i in range(5)]))
    Tb = arb_mat(16, 11)
    for jj, c in enumerate(free):
        Tb[c, jj] = 1
        for ii, p in enumerate(piv):
            Tb[p, jj] = -X[ii, jj]
    R = arb_mat(Jc) * Tb
    res = [R[i, j] for i in range(5) for j in range(11)]
    ok_ker = all(_finite0(v) for v in res)
    JT = arb_mat(Jf) * Tb
    cols = data["cols"]
    ok_det = False
    if _int_cols(cols, 6, 11):
        Jmin = [[JT[i, c] for c in cols] for i in range(6)]
        ok_det = _nonzero(bareiss_det(Jmin))
        stored = data["J"]
        if (not isinstance(stored, list) or len(stored) != 6
                or any(len(r) != 6 or any(not ball(s).overlaps(v) for s, v in zip(r, rv))
                       for r, rv in zip(stored, Jmin))):
            raise CheckError("stored J differs from the re-derived minor (Jf Tb)[:, cols]")
    sc = data["stacked_cols"]
    ok_st = False
    if _int_cols(sc, 11, 16):
        d = bareiss_det([[(Jc[i] if i < 5 else Jf[i - 5])[c] for c in sc] for i in range(11)])
        ok_st = _nonzero(d)
        if not (isinstance(data["stacked_det"], str) and ball(data["stacked_det"]).overlaps(d)):
            raise CheckError("stored stacked_det differs from the re-derived determinant")
    if ok_ker:
        up = max((abs(v).upper() for v in res), key=lambda a: float(a))
        st = ball(data["kernel_residual_upper"])
        # the residual is rounding noise around 0.  The producer works from exact-ish Klein points, the checker from
        # the stored 30-digit balls, so its maximum is larger; the stored value must be a nonnegative number no
        # larger than the re-derived one (up to a factor 1e3)
        if not (st.is_finite() and st >= 0 and st <= 1000 * up + arb(2) ** -90):
            raise CheckError("stored kernel_residual_upper differs from the re-derived residual")
    if data["constraints"] is not None:
        sv = data["constraints"]
        if len(sv) != 5 or any(not ball(s).overlaps(v) for s, v in zip(sv, cval)):
            raise CheckError("stored constraints differ from the re-derived values")
    return [("det J[:, cols]", "!=0", "0", ok_det), ("det [Jc; Jf][:, stacked_cols]", "!=0", "0", ok_st),
            ("Jc . Tb (all entries)", "encloses", "0", ok_ker), ("constraints(x)", "encloses", "0", ok_con)]
