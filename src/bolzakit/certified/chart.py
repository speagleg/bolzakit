"""Teichmueller trace chart near Bolza: certified rank of d(pi) and the inverse chart tau -> octagon.

Ported from the Bolza lambda1 research code (cert_rank.py: constraints, bolza_klein; cert_dpi_rank.py: pairings,
traces, the rank-6 certificate).  The side pairings g_k(x) in SO+(2,1) are real-analytic in the 16 Klein vertex
coordinates x; the traces tau_w(x) = tr w(x) of six words w factor through T_2.  Forward-mode dual numbers (jets.D)
over Arb give exact-formula enclosures of the 6 x 16 trace Jacobian Jf and the 5 x 16 constraint Jacobian Jc.
A 6 x 6 minor of Jf restricted to ker Jc whose determinant excludes 0 certifies rank d(pi)(x) = 6 = dim T_2.

All work runs at 128 bits (as in the SRC), set locally with _arb.prec.
"""
from __future__ import annotations

import itertools
from dataclasses import dataclass

from flint import arb, arb_mat

from ..errors import OutOfValidityRange
from ._arb import ball_str, certainly_lt, prec
from .certificate import Certificate, Claim, Dependency, make_producer
from .jets import D

_PREC = 128
PIV = (0, 2, 5, 7, 8)  # certified nonzero 5x5 minor of Jc at Bolza (SRC cert_rank.out)
# the six trace functions: words in the side pairings (G = inverse); chosen from a float scan
# (SRC scratch_G3/explore.py) as the best-conditioned 6-subset of the length-2 words
WORDS = [("g0g1", [(0, 1), (1, 1)]),
         ("g0G2", [(0, 1), (2, -1)]),
         ("g1g2", [(1, 1), (2, 1)]),
         ("g1g3", [(1, 1), (3, 1)]),
         ("g1G3", [(1, 1), (3, -1)]),
         ("g2g3", [(2, 1), (3, 1)])]
SIG = [1, -1, -1]
TAU_REL_RANGE = 0.05      # from_trace_coords is validated for |tau - tau_B|_inf <= 0.05 max(tau_B)
# Krawczyk box radii, tightest first, x10 up to 1e-3.  Naive interval evaluation of J(X) amplifies box width by
# ~1e7 (hyperboloid lift + frames), so at 128 bits the inclusion verifies for radii ~1e-30 .. 1e-11 only; boxes of
# radius >= 1e-9 cannot even certify the frame orientation (measured at tau_B + (0.5, 0, 0, -0.25, 0, 0)).
_KRAWCZYK_RADII = tuple(10.0 ** -k for k in range(28, 2, -1))


# ------------------------------------------------------------------ SRC cert_rank.py
def acos_d(x):
    return D(x.v.acos(), -x.d / (1 - x.v * x.v).sqrt())


def Bm(a, b):
    return a[0] * b[0] - a[1] * b[1] - a[2] * b[2]


def lift(p):
    t = 1 / (1 - p[0] * p[0] - p[1] * p[1]).sqrt()
    return [t, t * p[0], t * p[1]]


def constraints(xs):
    """xs: 16 dual numbers (Klein coords of V_0..V_7).  Returns the 5 constraint values (dual)."""
    H = [lift((xs[2 * k], xs[2 * k + 1])) for k in range(8)]
    L = [Bm(H[k], H[(k + 1) % 8]).acosh() for k in range(8)]
    ang = []
    for k in range(8):
        a, p, q = H[k], H[(k - 1) % 8], H[(k + 1) % 8]
        w1 = [p[i] - Bm(a, p) * a[i] for i in range(3)]
        w2 = [q[i] - Bm(a, q) * a[i] for i in range(3)]
        Bp = lambda u, v: -Bm(u, v)  # noqa: E731
        ang.append(acos_d(Bp(w1, w2) / (Bp(w1, w1) * Bp(w2, w2)).sqrt()))
    S = ang[0]
    for k in range(1, 8):
        S = S + ang[k]
    return [L[k + 4] - L[k] for k in range(4)] + [S - 2 * arb.pi()]


def bolza_klein():
    R = (3 + 2 * arb(2).sqrt()).acosh(); rho = R.tanh()
    out = []
    for k in range(8):
        th = arb.pi() * k / 4
        out += [rho * th.cos(), rho * th.sin()]
    return out


# ------------------------------------------------------------------ SRC cert_dpi_rank.py
def val(z):
    return z.v if isinstance(z, D) else z


def frame(P, Q):
    """Lorentz-orthonormal frame [P, e, n] (columns) at P, e = unit tangent toward Q, det > 0.
    Same frame as octagon.frame, with the normal written explicitly: n = J (P x e), normalised."""
    e = [Q[i] - Bm(P, Q) * P[i] for i in range(3)]
    s = (-Bm(e, e)).sqrt()
    e = [ei / s for ei in e]
    c = [P[1] * e[2] - P[2] * e[1], P[2] * e[0] - P[0] * e[2], P[0] * e[1] - P[1] * e[0]]
    n = [SIG[i] * c[i] for i in range(3)]
    s = (-Bm(n, n)).sqrt()
    n = [ni / s for ni in n]
    F = [[P[i], e[i], n[i]] for i in range(3)]
    det = (F[0][0] * (F[1][1] * F[2][2] - F[1][2] * F[2][1]) - F[0][1] * (F[1][0] * F[2][2] - F[1][2] * F[2][0])
           + F[0][2] * (F[1][0] * F[2][1] - F[1][1] * F[2][0]))
    dv = val(det)
    if dv.contains(0):
        raise ZeroDivisionError("frame determinant not certified nonzero")
    if dv < 0:
        for i in range(3):
            F[i][2] = -F[i][2]
    return F


def matmul(A, B):
    return [[A[i][0] * B[0][j] + A[i][1] * B[1][j] + A[i][2] * B[2][j] for j in range(3)] for i in range(3)]


def linv(F):
    """Inverse of a Lorentz-orthonormal frame: J F^T J (exact identity since F^T J F = J)."""
    return [[SIG[i] * F[j][i] * SIG[j] for j in range(3)] for i in range(3)]


def linv_iso(g):
    """Inverse of g in SO(2,1): J g^T J."""
    return [[SIG[i] * g[j][i] * SIG[j] for j in range(3)] for i in range(3)]


def pairings(xs):
    """g_k (k = 0..3): orientation-preserving isometry with V_{k+4} -> V_{k+1}, V_{k+5} -> V_k."""
    H = [lift((xs[2 * k], xs[2 * k + 1])) for k in range(8)]
    g = []
    for k in range(4):
        src = frame(H[k + 4], H[(k + 5) % 8])
        dst = frame(H[k + 1], H[k])
        g.append(matmul(dst, linv(src)))
    return g, H


def traces(xs):
    g, _ = pairings(xs)
    gi = [linv_iso(x) for x in g]
    out = []
    for _, w in WORDS:
        M = None
        for (a, s) in w:
            X = g[a] if s == 1 else gi[a]
            M = X if M is None else matmul(M, X)
        out.append(M[0][0] + M[1][1] + M[2][2])
    return out


def jacobians(x0):
    Jc = [[None] * 16 for _ in range(5)]
    Jf = [[None] * 16 for _ in range(len(WORDS))]
    cval = fval = None
    for j in range(16):
        xs = [D(x0[i], arb(1) if i == j else arb(0)) for i in range(16)]
        c = constraints(xs)
        f = traces(xs)
        cval = [ci.v for ci in c]
        fval = [fi.v for fi in f]
        for i in range(5):
            Jc[i][j] = c[i].d
        for i in range(len(WORDS)):
            Jf[i][j] = f[i].d
    return Jc, Jf, cval, fval


def bolza_traces() -> list:
    """The six trace coordinates (floats, midpoints at 128 bits) of the Bolza octagon, in WORDS order."""
    with prec(_PREC):
        return [float(t.mid()) for t in traces(bolza_klein())]


# ------------------------------------------------------------------ certified rank 6
@dataclass
class RankResult:
    det: str                 # ball string of the chosen 6 x 6 minor determinant
    cols: list               # tangent-basis columns (0..10) of that minor
    n_nonzero_minors: int    # 6 x 6 minors (of 462) whose determinant certainly excludes 0
    certified: bool
    certificate: Certificate


def _encloses_0(v) -> bool:
    """v is a finite ball containing 0 (arb nan / infinite balls also 'contain' 0 in python-flint)."""
    return bool(v.is_finite() and v.contains(0))


def _pivot(Jc):
    """PIV if its 5x5 minor of Jc certainly excludes 0, else the first such minor (None if there is none)."""
    for cols in [PIV] + list(itertools.combinations(range(16), 5)):
        d = arb_mat([[Jc[i][c] for c in cols] for i in range(5)]).det()
        if d.is_finite() and not d.contains(0):
            return tuple(cols)
    return None


def certify_rank(surface=None) -> RankResult:
    """Certify rank d(pi)(x) = 6 at the octagon x of `surface` (default: Bolza), SRC cert_dpi_rank.main.

    certified = (best 6x6 minor of Jf.Tb excludes 0) and (an 11x11 minor of [Jc; Jf] excludes 0) and
    (|Jc.Tb| encloses 0) and (the constraints at x enclose 0)."""
    with prec(_PREC):
        x0 = bolza_klein() if surface is None else surface._klein_arb()
        label = "bolza" if surface is None else surface.label
        try:
            Jc, Jf, cval, _ = jacobians(x0)
        except ZeroDivisionError:
            Jc, cval = None, None
        piv = None if Jc is None else _pivot(Jc)
        best, ncert, resid_upper, stacked = None, 0, None, None
        kernel_ok = False
        if piv is not None:
            free = [j for j in range(16) if j not in piv]
            Ap = arb_mat([[Jc[i][c] for c in piv] for i in range(5)])
            Af = arb_mat([[Jc[i][c] for c in free] for i in range(5)])
            X = Ap.solve(Af)  # 5 x 11 ; kernel vectors: x_free = e_j, x_piv = -X e_j
            Tb = [[arb(0)] * 11 for _ in range(16)]
            for jj, c in enumerate(free):
                Tb[c][jj] = arb(1)
                for ii, p in enumerate(piv):
                    Tb[p][jj] = -X[ii, jj]
            Tb = arb_mat(Tb)
            resid = arb_mat(Jc) * Tb
            kernel_ok = all(_encloses_0(resid[i, j]) for i in range(5) for j in range(11))
            resid_upper = max(abs(resid[i, j]).upper() for i in range(5) for j in range(11))
            JT = arb_mat(Jf) * Tb  # 6 x 11
            for cols in itertools.combinations(range(11), 6):
                d = arb_mat([[JT[i, c] for c in cols] for i in range(6)]).det()
                if d.is_finite() and not d.contains(0):
                    ncert += 1
                    if best is None or abs(d).lower() > abs(best[1]).lower():
                        best = (cols, d)
            # cross-check: an 11 x 11 minor of the stacked [Jc; Jf]
            S = arb_mat([Jc[i] for i in range(5)] + [Jf[i] for i in range(6)])
            for cols in itertools.combinations(range(16), 11):
                d = arb_mat([[S[i, c] for c in cols] for i in range(11)]).det()
                if d.is_finite() and not d.contains(0):
                    stacked = (list(cols), d)
                    break
        if best is None:
            cols, det, Jmin = [], arb("nan"), []
        else:
            cols, det = list(best[0]), best[1]
            Jmin = [[ball_str(JT[i, c]) for c in cols] for i in range(6)]
        constraints_ok = cval is not None and all(_encloses_0(v) for v in cval)
        # the certificate's claims are exactly the conditions of RankResult.certified
        claims = [
            Claim(lhs="det J[:, cols]", op="!=0", rhs="0", holds=bool(det.is_finite() and not det.contains(0))),
            Claim(lhs="det [Jc; Jf][:, stacked_cols]", op="!=0", rhs="0", holds=stacked is not None),
            Claim(lhs="Jc . Tb (all entries)", op="encloses", rhs="0", holds=bool(kernel_ok)),
            Claim(lhs="constraints(x)", op="encloses", rhs="0", holds=bool(constraints_ok))]
        ok = all(c.holds for c in claims)
        data = {"J": Jmin, "cols": cols,
                "stacked_cols": stacked[0] if stacked else None,
                "stacked_det": ball_str(stacked[1]) if stacked else None,
                "kernel_residual_upper": ball_str(resid_upper) if resid_upper is not None else None,
                "constraints": [ball_str(v) for v in cval] if cval is not None else None}
        inputs = {"surface": label, "klein": [ball_str(c) for c in x0], "prec": _PREC,
                  "words": [w for w, _ in WORDS], "pivot_columns_Jc": list(piv) if piv else None}
        tau = getattr(surface, "trace_coords", None)
        if tau is not None:
            inputs["tau"] = list(tau)
        cert = Certificate(
            kind="rank", inputs=inputs, claims=claims, data=data,
            dependencies=[
                Dependency("D12", "T2 is a topological 6-manifold; invariance of domain; continuity of the "
                                  "quotient to T2; the trace functions tau_w are C^1; the extended mapping-class "
                                  "group action on T2 is continuous",
                           "Fenchel-Nielsen; Felikson-Natanzon arXiv:1102.4861 sec. 2; PROOF.md Lemma 9d2"),
                Dependency("D13", "software: python-flint (Arb) and sympy are trusted", "python-flint; sympy")],
            producer=make_producer())
        return RankResult(det=ball_str(det), cols=cols, n_nonzero_minors=ncert, certified=bool(ok),
                          certificate=cert)


# ------------------------------------------------------------------ inverse chart tau -> octagon (Krawczyk)
def _unknowns():
    """The 11 unknown Klein coordinates: PIV plus the free coordinates of the Bolza rank-6 minor columns
    (tangent columns (1, 3, 5, 7, 8, 10) -> free coords (3, 6, 10, 12, 13, 15)).  The stacked Jacobian on these
    columns has det = det(Jc[:, PIV]) * det(6x6 rank minor) != 0 at Bolza (Schur complement)."""
    free = [j for j in range(16) if j not in PIV]
    return sorted(list(PIV) + [free[c] for c in (1, 3, 5, 7, 8, 10)])


def _system(x, tau, cols):
    """F(x) = [constraints(x); traces(x) - tau] (11 arb values) and its 11 x 11 Jacobian in the unknowns `cols`."""
    rows = [[None] * 11 for _ in range(11)]
    F = None
    for jj, j in enumerate(cols):
        xs = [D(x[i], arb(1) if i == j else arb(0)) for i in range(16)]
        out = constraints(xs) + [t - ti for t, ti in zip(traces(xs), tau)]
        F = [o.v for o in out]
        for i in range(11):
            rows[i][jj] = out[i].d
    return F, arb_mat(rows)


def _mid(x):
    return arb(x.mid())


def solve_trace_coords(tau) -> list:
    """Klein balls (16 arb) of the unique octagon near Bolza with traces tau (and the 5 unknowns
    outside _unknowns() fixed at their Bolza values), certified by a Krawczyk inclusion.

    Raises OutOfValidityRange if |tau - tau_B|_inf > 0.05 max(tau_B) or the inclusion cannot be verified."""
    tau = [float(t) for t in tau]
    if len(tau) != len(WORDS):
        raise ValueError(f"expected {len(WORDS)} trace coordinates ({', '.join(w for w, _ in WORDS)})")
    tB = bolza_traces()
    if not max(abs(a - b) for a, b in zip(tau, tB)) <= TAU_REL_RANGE * max(tB):
        raise OutOfValidityRange(f"|tau - tau_B|_inf exceeds the validated range {TAU_REL_RANGE} * max(tau_B)")
    fail = OutOfValidityRange(f"the trace chart could not be certified at tau={tau!r}")
    cols = _unknowns()
    with prec(_PREC):
        ta = [arb(t) for t in tau]
        x = bolza_klein()
        x = [_mid(c) if i in cols else c for i, c in enumerate(x)]
        try:
            for _ in range(40):   # Newton (untrusted proposal)
                F, Jm = _system(x, ta, cols)
                step = Jm.solve(arb_mat([[f] for f in F]))
                for jj, j in enumerate(cols):
                    x[j] = _mid(x[j] - step[jj, 0])
                if max(abs(float(step[i, 0].mid())) for i in range(11)) < 1e-30:
                    break
            # Krawczyk: K(X) = x0 - Y F(x0) + (I - Y J(X)) (X - x0), Y = mid(J(x0))^-1; K(X) in int X => unique zero
            F0, J0 = _system(x, ta, cols)
            Yi = arb_mat([[_mid(J0[i, j]) for j in range(11)] for i in range(11)]).inv()
            Y = arb_mat([[_mid(Yi[i, j]) for j in range(11)] for i in range(11)])
            yF = Y * arb_mat([[f] for f in F0])
            I11 = arb_mat(11, 11, [1 if i == j else 0 for i in range(11) for j in range(11)])
        except (ZeroDivisionError, ValueError):   # singular / non-finite Jacobian along the Newton path
            raise fail from None
        for r in _KRAWCZYK_RADII:   # smallest radius that verifies
            dX = arb_mat([[arb(0, r)] for _ in range(11)])
            X = list(x)
            for jj, j in enumerate(cols):
                X[j] = x[j] + dX[jj, 0]
            try:
                _, JX = _system(X, ta, cols)
            except (ZeroDivisionError, ValueError):
                continue
            K = -yF + (I11 - Y * JX) * dX   # K(X) - x0
            if all(certainly_lt(abs(K[i, 0]), r) for i in range(11)):
                out = list(x)
                for jj, j in enumerate(cols):
                    out[j] = x[j] + K[jj, 0]
                return out
        raise fail
