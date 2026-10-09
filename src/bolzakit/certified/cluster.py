"""Certified perturbation tests of the multiplicity-3 lambda_1 cluster at Bolza.

First order (degree=1): a port of the Bolza lambda1 research code (cert_first_order.py: element_dual,
series_integrator, ball_simplex_eval and the eta transfer). Its hard-coded lambda_1 ball and lambda_4 >= 4 are replaced
by the explicit inputs ``lam1`` and ``gap``.

Deformation F_A: hyperbolic vertex radii R_k(t) = R + t (-1)^k + sigma(t), angles k pi/4 fixed, sigma'(0) = 0, so the
first-order vertex derivative is mesh.vertices_dual().  On each geodesic triangle (hyperboloid-barycentric pullback):
   q_t(f) = sum_T int grad f^T a_T(t) grad f,  a_T = adj(L L^T - Q G) / (D sqrt Q)
   m_t(f) = sum_T int f^2 rho_T(t),            rho_T = D Q^{-3/2}
First variation b(f,g) = q'(f,g) - lambda1 m'(f,g).  The first-order matrix A = [b(e_i,e_j)] on the exact
lambda_1-eigenspace E is traceless, so A != 0 iff b(e,e) < 0 for some e in E.  We certify b(e,e) < 0 for e = P_E u,
u an (untrusted) FE trial function made exactly m-orthogonal to the constants, via the gap transfer
   m(u - e) <= (mu - lambda1)/(lambda4 - lambda1) m(u),  q(u - e) <= (mu - lambda1) lambda4/(lambda4 - lambda1) m(u)
and |b(f,g)| <= beta_a q(f)^1/2 q(g)^1/2 + lambda1 beta_rho m(f)^1/2 m(g)^1/2 (elementwise relative bounds in Arb).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla
from flint import arb

from ..errors import CertificationFailed, OutOfValidityRange
from . import forms
from ._arb import ball_str, certainly_gt, certainly_lt, prec
from ._hyp import Bm
from .certificate import Certificate, Claim, Dependency, make_producer
from .jets import D
from .mesh import Mesh, positions_dual, vertices_dual


@dataclass
class Lam1Input:
    lower: str
    upper: str
    source: str          # decimal strings; e.g. SU13


@dataclass
class GapInput:
    lower: str
    source: str          # lower bound for lambda_4; e.g. "4.0", "LP count (certified)"


@dataclass
class ClusterResult:
    A_float: list                      # eigenvalues of the float first-order matrix (splitting rates per unit t)
    b_ee_upper: float | None           # certified upper bound of min b(e,e)/m (first order)
    A_nonzero: bool
    kyfan_upper: float | None          # Task 9
    certified: bool
    certificates: list = field(default_factory=list)


# ------------------------------------------------------------------ element forms with t-derivative (dual numbers)
def _pv(p):
    return {k: c.v for k, c in p.items()}


def _pd(p):
    return {k: c.d for k, c in p.items()}


def series_integrator(Qv):
    """returns integ(p, s) = enclosure of int_simplex p Q^{-s} for arb polynomials p (Q fixed, arb)."""
    third = arb(1) / 3
    Qc = sum((c * third ** (i + j) for (i, j), c in Qv.items()), arb(0))
    eps = forms.padd(forms.pscale(Qv, 1 / Qc), {(0, 0): arb(1)}, -1)
    emax = forms.psup_centered(eps)
    if not certainly_lt(emax, 0.5):
        raise OutOfValidityRange("element too large for the series bound (emax >= 0.5)")
    K = forms.K_SERIES + 1
    epows = [{(0, 0): arb(1)}]
    for k in range(K):
        epows.append(forms.pmul(epows[-1], eps))

    def integ(p, s):
        tot = arb(0)
        for k in range(K + 1):
            tot += forms.binom_neg(s, k) * forms.pint(forms.pmul(p, epows[k]))
        rem = abs(forms.binom_neg(s, K + 1)) * emax ** (K + 1) * (1 - emax) ** (-s - K - 1)
        err = forms.psup_simplex(p) * rem / 2
        return (tot + arb(0, float(err.upper()) * (1 + 1e-9) + 1e-300)) * Qc ** (-s)
    return integ


def ball_simplex_eval(p, a, b):
    return sum((c * a ** i * b ** j for (i, j), c in p.items()), arb(0))


def element_dual(P0, P1, P2):
    """returns Ke, Me, dKe, dMe (3x3 arb) and local relative-derivative bounds beta_a, beta_rho (arb upper ends).

    P0..P2: dual (jets.D) hyperboloid triples.  -beta_a a_T <= a_T' <= beta_a a_T and |rho_T'| <= beta_rho rho_T hold
    pointwise on the simplex.  A non-certifiable bound gives beta = +inf (never a pass)."""
    W1 = [P0[i] - P2[i] for i in range(3)]; W2 = [P1[i] - P2[i] for i in range(3)]
    G11, G12, G22 = Bm(W1, W1), Bm(W1, W2), Bm(W2, W2)
    b1, b2 = Bm(W1, P2), Bm(W2, P2)
    Q = {(0, 0): Bm(P2, P2), (1, 0): 2 * b1, (0, 1): 2 * b2, (2, 0): G11, (1, 1): 2 * G12, (0, 2): G22}
    L1 = {(0, 0): b1, (1, 0): G11, (0, 1): G12}
    L2 = {(0, 0): b2, (1, 0): G12, (0, 1): G22}
    det = (P0[0] * (P1[1] * P2[2] - P1[2] * P2[1]) - P1[0] * (P0[1] * P2[2] - P0[2] * P2[1])
           + P2[0] * (P0[1] * P1[2] - P0[2] * P1[1]))
    if certainly_gt(det.v, 0):
        sgn = 1
    elif certainly_lt(det.v, 0):
        sgn = -1
    else:
        raise OutOfValidityRange("degenerate element: orientation not certain")
    Dd = D(det.v * sgn, det.d * sgn)
    N11 = forms.padd(forms.pmul(L2, L2), {k: v * G22 for k, v in Q.items()}, -1)
    N12 = forms.padd({k: v * G12 for k, v in Q.items()}, forms.pmul(L1, L2), -1)
    N22 = forms.padd(forms.pmul(L1, L1), {k: v * G11 for k, v in Q.items()}, -1)
    Qv, Qd = _pv(Q), _pd(Q)
    integ = series_integrator(Qv)
    half, th = arb(1) / 2, arb(3) / 2

    def I_and_dI(Np, s):
        """int Np Q^{-s} and its derivative int Np' Q^{-s} - s int Np Q' Q^{-s-1}."""
        v = integ(_pv(Np), s)
        d = integ(_pd(Np), s) - s * integ(forms.pmul(_pv(Np), Qd), s + 1)
        return v, d
    I11, dI11 = I_and_dI(N11, half); I12, dI12 = I_and_dI(N12, half); I22, dI22 = I_and_dI(N22, half)
    Nm = [[I11, I12], [I12, I22]]; dNm = [[dI11, dI12], [dI12, dI22]]
    grads = [(1, 0), (0, 1), (-1, -1)]
    Ke = [[None] * 3 for _ in range(3)]; dKe = [[None] * 3 for _ in range(3)]
    for p in range(3):
        for q in range(p, 3):
            gp, gq = grads[p], grads[q]
            s0 = sum((gp[i] * gq[j] * Nm[i][j] for i in range(2) for j in range(2)), arb(0))
            s1 = sum((gp[i] * gq[j] * dNm[i][j] for i in range(2) for j in range(2)), arb(0))
            Ke[p][q] = Ke[q][p] = s0 / Dd.v
            dKe[p][q] = dKe[q][p] = s1 / Dd.v - s0 * Dd.d / (Dd.v * Dd.v)
    phis = [{(1, 0): arb(1)}, {(0, 1): arb(1)}, {(0, 0): arb(1), (1, 0): arb(-1), (0, 1): arb(-1)}]
    Me = [[None] * 3 for _ in range(3)]; dMe = [[None] * 3 for _ in range(3)]
    for p in range(3):
        for q in range(p, 3):
            pp = forms.pmul(phis[p], phis[q])
            v = integ(pp, th); d = -th * integ(forms.pmul(pp, Qd), th + 1)
            Me[p][q] = Me[q][p] = Dd.v * v
            dMe[p][q] = dMe[q][p] = Dd.d * v + Dd.v * d
    # pointwise relative bounds over the box a, b in [0,1] (contains the simplex):
    #  a^{-1} a' = N^{-1} N' - (D'/D + Q'/(2Q)) I;  beta_a >= ||N^{-1} N'||_F + |D'/D + Q'/2Q|
    A_ = arb(0.5, 0.5); B_ = arb(0.5, 0.5)
    Qb = ball_simplex_eval(Qv, A_, B_); Qdb = ball_simplex_eval(Qd, A_, B_)
    Nb = [[ball_simplex_eval(_pv(N11), A_, B_), ball_simplex_eval(_pv(N12), A_, B_)],
          [ball_simplex_eval(_pv(N12), A_, B_), ball_simplex_eval(_pv(N22), A_, B_)]]
    Ndb = [[ball_simplex_eval(_pd(N11), A_, B_), ball_simplex_eval(_pd(N12), A_, B_)],
           [ball_simplex_eval(_pd(N12), A_, B_), ball_simplex_eval(_pd(N22), A_, B_)]]
    detN = Nb[0][0] * Nb[1][1] - Nb[0][1] * Nb[1][0]
    if not (certainly_gt(detN, 0) and certainly_gt(Nb[0][0], 0)):
        # N not certainly positive definite on the box: no relative bound (SRC asserted here)
        inf = arb("inf")
        return Ke, Me, dKe, dMe, inf, inf
    Ninv = [[Nb[1][1] / detN, -Nb[0][1] / detN], [-Nb[1][0] / detN, Nb[0][0] / detN]]
    Mprod = [[sum((Ninv[i][k] * Ndb[k][j] for k in range(2)), arb(0)) for j in range(2)] for i in range(2)]
    fro = sum((Mprod[i][j] ** 2 for i in range(2) for j in range(2)), arb(0)).sqrt()
    scal = Dd.d / Dd.v + Qdb / (2 * Qb)
    beta_a = (fro + abs(scal)).upper()
    beta_rho = (abs(Dd.d / Dd.v - th * Qdb / Qb)).upper()
    return Ke, Me, dKe, dMe, beta_a, beta_rho


def _all_finite(elem):
    Ke, Me, dKe, dMe, ba, br = elem
    return (ba.is_finite() and br.is_finite()
            and all(x.is_finite() for Mx in (Ke, Me, dKe, dMe) for row in Mx for x in row))


def _float_up(x):
    """smallest-ish float certainly >= the upper end of the finite ball x."""
    up = x.upper()
    f = float(up)
    while not arb(f) >= up:
        f = math.nextafter(f, math.inf)
    return f


# ------------------------------------------------------------------ first-order test
def _first_order(n, lam1_in, gap_in, sym):
    """Certified upper bound of b(e,e)/m(u) for the F_A deformation, P1 trial space on Mesh(n).

    Certificate data (assembly "b_ee_upper_v1"):
      data["terms"] = {"b_uu_over_m": B, "transfer_err_over_m": E}   (Arb ball strings)
      data["b_ee_upper_over_m"] = ball of  B + E   (Arb ball addition)
      claim: b_ee_upper_over_m < 0, i.e. the upper end of B + E is certainly negative.
    Here B = b(u,u)/m(u) and E = err/m(u), where b(e,e) = b(u,u) - 2 b(u,w) + b(w,w) (w = u - e) and
      E = 2 (sqrt(Sa/m) sqrt(qw) + lam1 sqrt(Sr/m) sqrt(mw)) + beta_a qw + lam1 beta_rho mw,
      mw = (mu - lam1)/(gap - lam1),  qw = (mu - lam1) gap/(gap - lam1),
    with lam1 the input ball [lam1.lower, lam1.upper], gap = gap.lower, and, all under data["transfer_terms"]:
    mu = q(u)/m(u), Sa/m = sum_T beta_a,T^2 q_T(u)/m(u), Sr/m = sum_T beta_rho,T^2 m_T(u)/m(u), beta_a = max_T beta_a,T,
    beta_rho = max_T beta_rho,T (mw, qw as computed are also stored there).  Recomputing E from transfer_terms gives
    another enclosure of the same quantity.  NaN or infinite balls anywhere never pass.
    """
    mesh = Mesh(n)
    X = positions_dual(mesh, vertices_dual())
    T = mesh.T; lab = mesh.lab
    elems = []; rows = []; cols = []; kv = []; mv = []; dkv = []; dmv = []
    # pi/2 rotation shortcut (SRC SYM=1): F_A and the mesh construction are invariant under V_k -> V_{k+2}, and
    # mesh.T lists the 8 sectors as consecutive congruent blocks, so sector k+2 has the exact element data of sector k.
    per = len(T) // 8
    cache = {}
    finite = True
    for ti, (i0, i1, i2) in enumerate(T):
        key = (divmod(ti, per)[0] % 2, ti % per)
        if sym and key in cache:
            el = cache[key]
        else:
            el = element_dual(X[i0], X[i1], X[i2])
            if sym:
                cache[key] = el
        finite = finite and _all_finite(el)
        elems.append(el)
        Ke, Me, dKe, dMe, _, _ = el
        L = [lab[i0], lab[i1], lab[i2]]
        for p in range(3):
            for q in range(3):
                rows.append(L[p]); cols.append(L[q])
                kv.append(float(Ke[p][q].mid())); mv.append(float(Me[p][q].mid()))
                dkv.append(float(dKe[p][q].mid())); dmv.append(float(dMe[p][q].mid()))
    lam1 = arb(lam1_in.lower).union(arb(lam1_in.upper))
    lam4_low = arb(gap_in.lower)
    inputs = {"direction": "F_A", "degree": 1, "n": n, "sym": sym, "surface": "bolza",
              "lam1": {"lower": lam1_in.lower, "upper": lam1_in.upper, "source": lam1_in.source},
              "gap": {"lower": gap_in.lower, "source": gap_in.source}}
    deps = [Dependency("D1", f"lambda_1 in [{lam1_in.lower}, {lam1_in.upper}]", lam1_in.source),
            Dependency("gap", f"lambda_4 >= {gap_in.lower} (no eigenvalue in (lambda_1, {gap_in.lower}))",
                       gap_in.source),
            Dependency("F1", "C^1 regularity of the lambda_1 cluster under the deformation (A is the first-order "
                             "matrix of the cluster, traceless)", "PROOF.md Lemma 2 (F1_PROOF.md)")]
    if not finite:
        # NaN / infinite element data: no float model and no certified bound (fail closed)
        cert = Certificate(kind="cluster_first_order", inputs=inputs,
                           claims=[Claim(lhs="b_ee_upper_over_m", op="<", rhs="0", holds=False)],
                           data={"assembly": "b_ee_upper_v1", "terms": {}, "b_ee_upper_over_m": "nan",
                                 "transfer_terms": {}, "diagnostics": {"error": "non-finite element data"}},
                           dependencies=deps, producer=make_producer())
        return [], None, False, cert
    sh = (mesh.Nf, mesh.Nf)
    Kf = sp.csr_matrix((kv, (rows, cols)), shape=sh); Mf = sp.csr_matrix((mv, (rows, cols)), shape=sh)
    dKf = sp.csr_matrix((dkv, (rows, cols)), shape=sh); dMf = sp.csr_matrix((dmv, (rows, cols)), shape=sh)
    w, U = sla.eigsh(Kf, k=6, M=Mf, sigma=-0.5,
                     v0=np.random.RandomState(0).rand(mesh.Nf))  # fixed start: reproducible eigenvectors
    o = np.argsort(w); w, U = w[o], U[:, o]
    Phi = U[:, 1:4]; lamf = w[1:4].mean()
    Af = Phi.T @ (dKf - lamf * dMf) @ Phi; Af = (Af + Af.T) / 2
    ev, evec = np.linalg.eigh(Af)
    xvec = Phi @ evec[:, 0]                          # most negative direction (float, untrusted)
    # exact quantities for vectors [1, x] (arb)
    xa = [arb(float(v)) for v in xvec]
    S = {k: arb(0) for k in ["q_xx", "m_xx", "m_x1", "m_11", "dq_xx", "dm_xx", "dm_x1", "dm_11"]}
    for t, (i0, i1, i2) in enumerate(T):
        Ke, Me, dKe, dMe, _, _ = elems[t]
        lx = [xa[lab[i0]], xa[lab[i1]], xa[lab[i2]]]
        for p in range(3):
            for q in range(3):
                S["q_xx"] += lx[p] * Ke[p][q] * lx[q]; S["m_xx"] += lx[p] * Me[p][q] * lx[q]
                S["m_x1"] += lx[p] * Me[p][q]; S["m_11"] += Me[p][q]
                S["dq_xx"] += lx[p] * dKe[p][q] * lx[q]; S["dm_xx"] += lx[p] * dMe[p][q] * lx[q]
                S["dm_x1"] += lx[p] * dMe[p][q]; S["dm_11"] += dMe[p][q]
    # u = x - c 1 with c = m(x,1)/m(1,1): u is m-orthogonal to constants exactly
    c = S["m_x1"] / S["m_11"]
    q_u = S["q_xx"]
    m_u = S["m_xx"] - c * S["m_x1"]
    b_xx = S["dq_xx"] - lam1 * S["dm_xx"]
    b_x1 = -lam1 * S["dm_x1"]                         # q'(x,1) = 0
    b_11 = -lam1 * S["dm_11"]
    b_uu = b_xx - 2 * c * b_x1 + c * c * b_11
    # Sa = sum_T beta_a,T^2 q_T(u),  Sr = sum_T beta_rho,T^2 m_T(u);  beta = max_T (all elements finite here)
    Sa = arb(0); Sr = arb(0); beta_a = arb(0); beta_rho = arb(0)
    for t, (i0, i1, i2) in enumerate(T):
        Ke, Me, dKe, dMe, bta, btr = elems[t]
        lu = [xa[lab[i0]] - c, xa[lab[i1]] - c, xa[lab[i2]] - c]
        qT = sum((lu[p] * Ke[p][q] * lu[q] for p in range(3) for q in range(3)), arb(0))
        mT = sum((lu[p] * Me[p][q] * lu[q] for p in range(3) for q in range(3)), arb(0))
        Sa += bta * bta * qT; Sr += btr * btr * mT
        beta_a = bta if bta > beta_a else beta_a; beta_rho = btr if btr > beta_rho else beta_rho
    mu = q_u / m_u
    gap = lam4_low - lam1
    exc = mu - lam1                                     # >= 0 (lambda1 = min on 1-perp)
    m_w = exc / gap * m_u
    q_w = exc * lam4_low / gap * m_u
    # |b(u,w)| <= Sa^{1/2} q(w)^{1/2} + lam1 Sr^{1/2} m(w)^{1/2};  |b(w,w)| <= beta_a q(w) + lam1 beta_rho m(w)
    err = 2 * (Sa.sqrt() * q_w.sqrt() + lam1 * Sr.sqrt() * m_w.sqrt()) + beta_a * q_w + lam1 * beta_rho * m_w
    B = b_uu / m_u
    E = err / m_u
    bound = B + E
    ok = certainly_lt(bound, 0)                         # NaN / infinite never pass
    b_ee_upper = _float_up(bound) if bound.is_finite() else None
    data = {
        "assembly": "b_ee_upper_v1",
        "terms": {"b_uu_over_m": ball_str(B), "transfer_err_over_m": ball_str(E)},
        "b_ee_upper_over_m": ball_str(bound),
        "transfer_terms": {"mu": ball_str(mu), "m_u": ball_str(m_u), "Sa_over_m": ball_str(Sa / m_u),
                           "Sr_over_m": ball_str(Sr / m_u), "beta_a": ball_str(beta_a),
                           "beta_rho": ball_str(beta_rho), "m_w_over_m": ball_str(m_w / m_u),
                           "q_w_over_m": ball_str(q_w / m_u)},
        "diagnostics": {"non_additive": True, "beta_a": float(beta_a.upper()), "beta_rho": float(beta_rho.upper()),
                        "beta_rms_a": ball_str((Sa / q_u).sqrt()), "beta_rms_rho": ball_str((Sr / m_u).sqrt()),
                        "float_cluster": w[1:4].tolist(), "float_A_eigs_per_unit_t": ev.tolist()},
    }
    cert = Certificate(kind="cluster_first_order", inputs=inputs,
                       claims=[Claim(lhs="b_ee_upper_over_m", op="<", rhs="0", holds=ok)],
                       data=data, dependencies=deps, producer=make_producer())
    return ev.tolist(), b_ee_upper, ok, cert


def certify_cluster(direction="F_A", degree=1, n=96, lam1: Lam1Input = None, gap: GapInput = None,
                    second_order=False, sym=True, strict=False) -> ClusterResult:
    """Certified first-order test of the lambda_1 cluster at Bolza (A != 0).

    lam1, gap: required explicit inputs (lambda_1 ball and lambda_4 lower bound), recorded with their sources.
    degree=1: the published P1 test for direction "F_A"; see _first_order for the certified quantity and the
    certificate's b_ee_upper_v1 assembly formula.  Precision is set internally (96 bits)."""
    if lam1 is None or gap is None:
        raise ValueError("lam1 and gap inputs are required (spec §3.4)")
    if float(lam1.lower) > float(lam1.upper):
        raise ValueError(f"lam1 interval is empty: lower {lam1.lower} > upper {lam1.upper}")
    if float(gap.lower) <= float(lam1.upper):
        raise ValueError(f"gap lower bound {gap.lower} must exceed lam1 upper {lam1.upper}")
    if degree != 1 or second_order:
        raise NotImplementedError("degree=2 / second-order cluster test is Task 9")
    if direction != "F_A":
        raise NotImplementedError("first-order P1 test is published for F_A; use degree=2 via Task 9 for other "
                                  "directions")
    with prec(96):
        A_float, b_ee_upper, ok, cert = _first_order(n, lam1, gap, sym)
    if strict and not ok:
        raise CertificationFailed("could not certify b(e,e) < 0 (A != 0) at first order")
    return ClusterResult(A_float=A_float, b_ee_upper=b_ee_upper, A_nonzero=ok, kyfan_upper=None, certified=ok,
                         certificates=[cert])
