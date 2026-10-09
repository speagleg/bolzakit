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
from flint import arb, arb_mat

from ..errors import CertificationFailed, OutOfValidityRange
from . import forms
from ._arb import PM, attested_upper, ball_str, certainly_gt, certainly_lt, json_float, prec
from ._hyp import Bm
from .certificate import Certificate, Claim, Dependency, make_producer, trial_provenance
from .jets import D
from .mesh import Mesh, positions_dual, positions_jet, vertices_dual
from .surface import Surface


@dataclass
class Lam1Input:
    lower: str
    upper: str
    source: str          # decimal strings; e.g. SU13
    # dep_id: the dependency id recorded in the certificate.  A PROOF.md §Dependencies number only when the source IS
    # that input: "D1" for SU13, "D15" for the Lehmann-Goerisch route (lower end from the lg_lower certificate);
    # otherwise "lam1" (an unnumbered hand-linked input)
    dep_id: str = "lam1"


@dataclass
class GapInput:
    lower: str
    source: str          # lower bound for lambda_4; e.g. "4.0", "LP count (certified)"


@dataclass
class ClusterResult:
    A_float: list                      # eigenvalues of the float first-order matrix A (splitting rates per unit t;
                                       # degree 2: A_ij = b(x_i, x_j) on the P2 float eigenvectors)
    b_ee_upper: float | None           # certified upper bound of min b(e,e)/m (first order), certificate-attested:
                                       # >= the upper end of the stored ball data["b_ee_upper_over_m"]
    A_nonzero: bool
    kyfan_upper: float | None          # certified upper bound of D2 Lambda(v,v) (degree 2, second order),
                                       # certificate-attested: >= the upper end of the stored data["g2_enclosure"]
    certified: bool
    certificates: list = field(default_factory=list)
    b_ee_upper_inprocess: float | None = None    # in-memory ball's upper end (<= b_ee_upper; fidelity pins)
    kyfan_upper_inprocess: float | None = None   # in-memory ball's upper end (<= kyfan_upper; fidelity pins)


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


def _cluster_deps(lam1_in, gap_in):
    """Dependencies shared by both cluster kinds: the lambda_1 ball (id from its source, see Lam1Input.dep_id), the gap,
    and D4 (FBP Lemmas 7.2-7.3: the lambda_1 eigenspace E has dimension 3, used for dim E = 3 in the gap transfer)."""
    return [Dependency(lam1_in.dep_id, f"lambda_1 in [{lam1_in.lower}, {lam1_in.upper}]", lam1_in.source),
            Dependency("gap", f"lambda_4 >= {gap_in.lower} (no eigenvalue in (lambda_1, {gap_in.lower}))",
                       gap_in.source),
            Dependency("D4", "lambda_1 has multiplicity exactly 3 (no 1- or 2-dimensional irreducible representation "
                             "in the lambda_1-eigenspace), so dim E = 3 in the gap transfer",
                       "FBP Lemmas 7.2 and 7.3 (quoted; PROOF.md Lemma 1c, §Dependencies D4)")]


def _float_up(x):
    """smallest-ish float certainly >= the upper end of the finite ball x."""
    up = x.upper()
    f = float(up)
    while not arb(f) >= up:
        f = math.nextafter(f, math.inf)
    return f


_attested_upper = attested_upper     # float >= the upper end of a stored ball string (see _arb.attested_upper)


# ------------------------------------------------------------------ first-order test
def _first_order(n, lam1_in, gap_in, sym):
    """Certified upper bound of b(e,e)/m(u) for the F_A deformation, P1 trial space on Mesh(n).

    Certificate data (assembly "b_ee_upper_v1"):
      data["raw"] (the checker's input; ball strings): the exact sums over elements for the float vector x and the
        constants 1:  q_xx, m_xx, m_x1, m_11 (forms q, m) and dq_xx, dm_xx, dm_x1, dm_11 (their t-derivatives q', m'),
        c = m_x1/m_11, Sa = sum_T beta_a,T^2 q_T(u), Sr = sum_T beta_rho,T^2 m_T(u) (u = x - c 1), beta_a, beta_rho.
        Rebuild:  m_u = m_xx - c m_x1,  q_u = q_xx (q annihilates constants),  mu = q_u/m_u,
          b_uu = (dq_xx - lam1 dm_xx) + 2 c lam1 dm_x1 - c^2 lam1 dm_11   (b = q' - lam1 m', q'(., 1) = 0),
          B = b_uu/m_u,  Sa/m = Sa/m_u,  Sr/m = Sr/m_u,  then E as below.
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
    deps = _cluster_deps(lam1_in, gap_in) + [
            Dependency("F1", "C^1 regularity of the lambda_1 cluster under the deformation (A is the first-order "
                             "matrix of the cluster, traceless)", "PROOF.md Lemma 2")]
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
    inputs["trial"] = trial_provenance([xvec], "RandomState(0).rand(Nf)")
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
        "raw": {**{k: ball_str(v) for k, v in S.items()}, "c": ball_str(c), "Sa": ball_str(Sa), "Sr": ball_str(Sr),
                "beta_a": ball_str(beta_a), "beta_rho": ball_str(beta_rho)},
        "terms": {"b_uu_over_m": ball_str(B), "transfer_err_over_m": ball_str(E)},
        "b_ee_upper_over_m": ball_str(bound),
        "transfer_terms": {"mu": ball_str(mu), "m_u": ball_str(m_u), "Sa_over_m": ball_str(Sa / m_u),
                           "Sr_over_m": ball_str(Sr / m_u), "beta_a": ball_str(beta_a),
                           "beta_rho": ball_str(beta_rho), "m_w_over_m": ball_str(m_w / m_u),
                           "q_w_over_m": ball_str(q_w / m_u)},
        "diagnostics": {"non_additive": True, "beta_a": json_float(beta_a.upper()),
                        "beta_rho": json_float(beta_rho.upper()),
                        "beta_rms_a": ball_str((Sa / q_u).sqrt()), "beta_rms_rho": ball_str((Sr / m_u).sqrt()),
                        "float_cluster": [json_float(x) for x in w[1:4]],
                        "float_A_eigs_per_unit_t": [json_float(x) for x in ev]},
    }
    cert = Certificate(kind="cluster_first_order", inputs=inputs,
                       claims=[Claim(lhs="b_ee_upper_over_m", op="<", rhs="0", holds=ok)],
                       data=data, dependencies=deps, producer=make_producer())
    return ev.tolist(), b_ee_upper, ok, cert


# ------------------------------------------------------------------ second-order Ky Fan bound (P2)
def _ub(x):
    """float >= the upper end of x (rounded outward); NaN for a non-finite ball (never a pass downstream)."""
    return _float_up(x) if x.is_finite() else math.nan


def _sqrt_up(x):
    """exact float point (as arb) >= sqrt(upper(x)); NaN if x is non-finite or certainly-not-nonnegative at the top."""
    return arb(_ub(arb(_ub(x)).sqrt()))


def _kyfan(direction, n, lam1_in, gap_in, sym, second_order):
    """Certified upper bound g2_upper >= g''(0) >= D2 Lambda(v, v), Lambda = lambda1 + lambda2 + lambda3 along the
    deformation `direction` ("F_A", "F_B" or a 16-vector in Klein coordinates) at Bolza; P2 trial space on Mesh(n).

    Port of the Bolza lambda1 research code (cert_sign.main), the F_B vertex source generalised to
    Surface.bolza().vertex_jets(direction).  Ky Fan: lambda_0 = 0, so for V(p) = span{1, e_i + p psi_i} (e_i a basis
    of the exact lambda1-eigenspace E, psi_i fixed P2 correctors) Lambda(p) <= g(p) = tr(G_m(p)^-1 G_q(p)) with
    equality at p = 0; with Lambda and g C^2 (F1), g - Lambda >= 0 has a minimum at 0, so Lambda''(0) <= g''(0);
    dLambda = 0 at Bolza (Lemma 4) makes Lambda''(0) = D2 Lambda(v, v) for every curve with tangent v.  Here
       g''(0) = tr(N X) - 2 tr(N M1 N A) + 2 lam c1^T N c1 / m(1,1),   N = M0^-1,
       X_ij  = (q''-lam m'')(e_i,e_j) + 2 b(psi_i,e_j) + 2 b(e_i,psi_j) + 2 (q-lam m)(psi_i,psi_j),
       A_ij  = b(e_i,e_j),  b = q' - lam m',  M0_ij = m(e_i,e_j),
       M1_ij = m'(e_i,e_j) + m(psi_i,e_j) + m(e_i,psi_j),   c1_i = m'(e_i,1) + m(psi_i,1).
    Each matrix entry at e_i = P_E u_i (u_i = float FE eigenvector made exactly m-orthogonal to 1) is the same entry
    at u plus a +-ball bounding the cross terms in w_i = u_i - e_i (gap transfer with lam in the lam1 input ball and
    lambda4 >= gap.lower:  m(w_i) <= (q(u_i) - lam m(u_i))/(lambda4 - lam),  q(w_i) <= that * lambda4, and the
    elementwise relative derivative bounds beta of forms.element_p2).

    Certificate data (assembly "g2_upper_v1").  All values are Arb ball strings unless noted; the checker works at
    96 bits.  PM = [0 +/- 1].
      data["raw"]: the complete input of the final assembly
        lam  = [lo, hi]  (decimal strings; lam := arb(lo).union(arb(hi)), the lam1 input ball)
        lam4 = decimal string (lam4 := arb(lam4), the gap input)
        lamU = float string with lamU >= upper(lam)
        G[k], k in K0 K1 K2 M0 M1 M2: 7x7 Gram matrices of the forms q, q', q''/2, m, m', m''/2 on the vectors
             f = 0..6 = [u_1, u_2, u_3, psi_1, psi_2, psi_3, 1]   (u_i = x_i - c_i 1 exactly m-orthogonal to 1)
        Wdiag[k], k in W1q W1m W2q W2m: 7-vectors, sum_T beta_T^2 q_T(f) resp. m_T(f) (beta1, beta1r, beta2, beta2r)
        bmax = [b1, b1r, b2, b2r]: max over elements of the four relative bounds (float points)
        W, sq, sm: the producer's float-point values of the derived radii below (informational; each is >= the
             checker's value up to the rounding of the stored strings)
      Rebuild, pure Arb from the above.  Derived radii (sqrt(upper(.)) of a ball = arb(upper).sqrt()):
        W[k]_f = sqrt(upper(Wdiag[k]_f))  (k in W1q W1m W2q W2m),   mps_f = sqrt(upper(G[M0]_ff)),
        exc_i = G[K0]_ii - lam G[M0]_ii,  sm_i = sqrt(upper(exc_i / (lam4 - lam))),
        sq_i = sqrt(upper(exc_i lam4 / (lam4 - lam)))   (exc_i < 0 gives NaN: never a pass)
      Then (i, j = 0..2;  u = i, p = 3 + i, P = 3 + j, ONE = 6;  B1 = G[K1] - lam G[M1],
      B2 = 2 (G[K2] - lam G[M2]),  S = G[K0] - lam G[M0]):
        M0e_ij = G[M0]_ij + sm_i sm_j PM
        Ae_ij  = B1_ij + (W1q_i sq_j + lamU W1m_i sm_j + W1q_j sq_i + lamU W1m_j sm_i + b1 sq_i sq_j
                          + lamU b1r sm_i sm_j) PM
        Xe_ij  = B2_ij + 2 B1_pj + 2 B1_iP + 2 S_pP
                 + (W2q_i sq_j + lamU W2m_i sm_j + W2q_j sq_i + lamU W2m_j sm_i + b2 sq_i sq_j + lamU b2r sm_i sm_j) PM
                 + (2 (W1q_p sq_j + lamU W1m_p sm_j) + 2 (W1q_P sq_i + lamU W1m_P sm_i)) PM
        M1e_ij = G[M1]_ij + G[M0]_pj + G[M0]_iP + (W1m_i sm_j + W1m_j sm_i + b1r sm_i sm_j + mps_p sm_j + mps_P sm_i) PM
        c1e_i  = G[M1]_i,ONE + G[M0]_p,ONE + W1m_ONE sm_i PM
        N = M0e^-1,  m11 = G[M0]_ONE,ONE
        trNX = tr(N Xe),  A_term = -2 tr(N M1e N Ae),  c1_term = 2 lam (c1e^T N c1e) / m11
        g2 = trNX + A_term + c1_term
      data["terms"] = {trNX, A_term, c1_term}: an additive set, sum = data["g2_enclosure"] = g2 (no double counting;
        the transfer errors are already inside each term).
      claim: g2_upper < 0, i.e. the upper end of g2 is certainly negative.
    data["diagnostics"] holds untrusted / non-additive information: g2_central (all PM terms dropped), the
    one-at-a-time error rows error_breakdown_upper_minus_central (each row includes ONE error group; they do NOT sum to
    g2_upper - g2_central), eta (gap-transfer) data, float estimates.  NaN or infinite data never pass.

    Vector directions: only the 2-jet of the vertex curve enters g''(0); (v_T, w) is the 2-jet of a curve in C
    (implicit function theorem, rank DC = 5).  inputs record v, the projected tangent v_T and |v - v_T|.

    second_order=False: float pass only (A_float = eigenvalues of the float A = [b(x_i, x_j)]), no certificate.
    Returns (A_float, g2_upper or None, ok, certificate or None)."""
    surface = Surface.bolza()
    d = Surface._direction(direction)
    named = isinstance(d, str)
    sym = bool(sym) and named          # the pi/2 element cache is valid only for r2-invariant (named) families
    mesh = Mesh(n); space = forms.P2Space(mesh)
    Vj = surface.vertex_jets(d)
    X = positions_jet(mesh, Vj)
    elems = forms.compute_elements(mesh, X, order=2, sym=sym)
    F = {}
    for wh, nm in ((0, "K"), (1, "M")):
        for o in range(3):
            F[nm + str(o)] = forms.assemble_float(space, elems, o, wh)
    inputs = {"direction": d if named else [repr(float(x)) for x in d], "degree": 2, "n": n, "sym": sym,
              "surface": "bolza",
              "lam1": {"lower": lam1_in.lower, "upper": lam1_in.upper, "source": lam1_in.source},
              "gap": {"lower": gap_in.lower, "source": gap_in.source}}
    if not named:
        vT = np.array([float((V[i] / V[0]).d.mid()) for V in Vj for i in (1, 2)])   # Klein tangent of the jets
        inputs["tangent_projected"] = [repr(float(t)) for t in vT]
        inputs["projection_residual"] = repr(float(np.linalg.norm(d - vT)))
    deps = _cluster_deps(lam1_in, gap_in) + [
            Dependency("D9", "Ky Fan: lambda_1+lambda_2+lambda_3 <= tr(G_m^-1 G_q) on any 4-dim space containing "
                             "the constants (lambda_0 = 0)",
                       "Ky Fan (PROOF.md §Dependencies D9: Moslehian arXiv:1108.1467, Fan 1950)"),
            Dependency("F1", "C^2 regularity of Lambda = lambda_1+lambda_2+lambda_3 along the deformation",
                       "PROOF.md Lemma 2"),
            Dependency("Lemma4", "dLambda = 0 at Bolza, so Lambda''(0) (<= g''(0) by Ky Fan) is the second derivative "
                                 "D^2 Lambda(v, v): independent of the curve's second-order term",
                       "PROOF.md Lemma 4 (dLambda = 0)")]
    claim_false = [Claim(lhs="g2_upper", op="<", rhs="0", holds=False)]
    if not all(np.isfinite(F[k].data).all() for k in F):
        # NaN / infinite element data: no float model and no certified bound (fail closed)
        if not second_order:
            return [], None, False, None
        cert = Certificate(kind="kyfan_second_order", inputs=inputs, claims=claim_false,
                           data={"assembly": "g2_upper_v1", "raw": {}, "terms": {}, "g2_enclosure": "nan",
                                 "diagnostics": {"non_additive": True, "error": "non-finite element data"}},
                           dependencies=deps, producer=make_producer())
        return [], None, False, cert
    K0, K1, K2, M0, M1, M2 = (F[k] for k in ("K0", "K1", "K2", "M0", "M1", "M2"))
    one = np.ones(space.ndof)
    w, U = sla.eigsh(K0, k=6, M=M0, sigma=-0.5,
                     v0=np.random.RandomState(0).rand(space.ndof))  # fixed start: reproducible eigenvectors
    o = np.argsort(w); w, U = w[o], U[:, o]
    x = U[:, 1:4].copy()
    x -= np.outer(one, (one @ (M0 @ x)) / (one @ (M0 @ one)))
    G = x.T @ (M0 @ x); L = np.linalg.cholesky(G); x = x @ np.linalg.inv(L).T    # M-orthonormalise (float)
    lamf = float(np.trace(x.T @ (K0 @ x)) / 3)
    B1 = (K1 - lamf * M1).tocsr()
    Af = x.T @ (B1 @ x)
    A_float = sorted(np.linalg.eigvalsh((Af + Af.T) / 2).tolist())
    if not second_order:
        return A_float, None, False, None
    # ---- correctors (float, untrusted): (K - lam M) psi_i = -(K' - lam M') x_i  on the M-complement of {x, 1}
    Cb = M0 @ np.column_stack([x, one])
    Abig = sp.bmat([[K0 - lamf * M0, sp.csr_matrix(Cb)], [sp.csr_matrix(Cb.T), None]]).tocsc()
    lu = sla.splu(Abig)
    psi = np.zeros((space.ndof, 3))
    for i in range(3):
        rhs = np.concatenate([-(B1 @ x[:, i]), np.zeros(4)])
        psi[:, i] = lu.solve(rhs)[:space.ndof]
    alpha = -0.5 * (x.T @ (M1 @ x))
    psi += x @ alpha.T
    m11f = one @ (M0 @ one)
    gam = -(x.T @ (M1 @ one)) / m11f
    psi += np.outer(one, gam)
    # float estimate of g''(0) (u = x), diagnostics only
    Xf = (2 * x.T @ ((K2 - lamf * M2) @ x) + 2 * (psi.T @ (B1 @ x) + x.T @ (B1 @ psi))
          + 2 * psi.T @ ((K0 - lamf * M0) @ psi))
    M1f = x.T @ (M1 @ x) + psi.T @ (M0 @ x) + x.T @ (M0 @ psi)
    c1f = x.T @ (M1 @ one) + psi.T @ (M0 @ one)
    T2f = float(np.trace(2 * x.T @ ((K2 - lamf * M2) @ x)))
    g2f = float(np.trace(Xf) - 2 * np.trace(M1f @ Af) + 2 * lamf * c1f @ c1f / m11f)
    # ---- certified pass: exact Gram matrices of the 7 float vectors [x_1..x_3, psi_1..psi_3, 1]
    V = np.column_stack([x, psi, one])
    inputs["trial"] = trial_provenance([V[:, j] for j in range(V.shape[1])], "RandomState(0).rand(ndof)")
    nv = 7
    Gs = {k: arb_mat(nv, nv) for k in ("K0", "K1", "K2", "M0", "M1", "M2", "W1q", "W1m", "W2q", "W2m")}
    bmax = [arb(0)] * 4
    finite = True
    cache = {}
    for t, e in enumerate(elems):
        key = id(e)
        if key not in cache:
            Ks, Ms, bet = e
            cache[key] = ([arb_mat(Ks[o_]) for o_ in range(3)], [arb_mat(Ms[o_]) for o_ in range(3)],
                          [arb(b) for b in bet])
        Km, Mm, bet = cache[key]
        Vl = arb_mat(V[space.loc[t]].tolist()); Vt = Vl.transpose()
        gk = [Vt * (Km[o_] * Vl) for o_ in range(3)]
        gm = [Vt * (Mm[o_] * Vl) for o_ in range(3)]
        for o_ in range(3):
            Gs["K" + str(o_)] += gk[o_]; Gs["M" + str(o_)] += gm[o_]
        Gs["W1q"] += (bet[0] * bet[0]) * gk[0]; Gs["W1m"] += (bet[1] * bet[1]) * gm[0]
        Gs["W2q"] += (bet[2] * bet[2]) * gk[0]; Gs["W2m"] += (bet[3] * bet[3]) * gm[0]
        for j in range(4):
            if not bet[j].is_finite():
                finite = False                      # SRC asserted; a NaN would be skipped by the max below
            elif bet[j] > bmax[j]:
                bmax[j] = bet[j]
    del cache
    if not finite:
        bmax = [arb("inf")] * 4                     # fail closed: every bound built on bmax is infinite
    # transform to u_i = x_i - c_i 1, c_i = m(x_i,1)/m(1,1) exactly
    c = [Gs["M0"][i, 6] / Gs["M0"][6, 6] for i in range(3)]
    Tm = arb_mat(nv, nv)
    for i in range(nv):
        Tm[i, i] = arb(1)
    for i in range(3):
        Tm[i, 6] = -c[i]
    Gu = {k: Tm * (g * Tm.transpose()) for k, g in Gs.items()}
    lam = arb(lam1_in.lower).union(arb(lam1_in.upper))
    lam4 = arb(gap_in.lower)
    B1u = Gu["K1"] - lam * Gu["M1"]                    # b = q' - lam m'
    B2u = 2 * (Gu["K2"] - lam * Gu["M2"])              # q'' - lam m''
    Su = Gu["K0"] - lam * Gu["M0"]                     # q - lam m
    U_, PS, ONE = [0, 1, 2], [3, 4, 5], 6
    # gap argument (a negative or NaN excess gives NaN square roots: never a pass)
    sq, sm, eta = [], [], []
    for i in range(3):
        exc = Gu["K0"][i, i] - lam * Gu["M0"][i, i]
        mw = exc / (lam4 - lam); qw = exc * lam4 / (lam4 - lam)
        sm.append(_sqrt_up(mw)); sq.append(_sqrt_up(qw))
        eta.append(dict(mu=ball_str(Gu["K0"][i, i] / Gu["M0"][i, i]), m_u=ball_str(Gu["M0"][i, i]),
                        sqrt_q_w_le=json_float(_ub(sq[-1])), sqrt_m_w_le=json_float(_ub(sm[-1]))))
    nrm = lambda G_, f: _sqrt_up(G_[f, f])           # point >= sqrt of the upper end of a (nonneg) weighted sum
    W1q = [nrm(Gu["W1q"], f) for f in range(7)]; W1m = [nrm(Gu["W1m"], f) for f in range(7)]
    W2q = [nrm(Gu["W2q"], f) for f in range(7)]; W2m = [nrm(Gu["W2m"], f) for f in range(7)]
    mps = [nrm(Gu["M0"], f) for f in range(7)]
    b1, b1r, b2, b2r = bmax
    lamU = arb(_ub(lam))

    def build(groups):
        """g2 and its three top-level parts; groups = the transfer-error groups included (diagnostic rows)."""
        M0e = arb_mat(3, 3); Ae = arb_mat(3, 3); Xe = arb_mat(3, 3); M1e = arb_mat(3, 3)
        c1e = arb_mat(3, 1)
        for i in range(3):
            for j in range(3):
                ui, uj, pi_, pj = U_[i], U_[j], PS[i], PS[j]
                r = sm[i] * sm[j]
                M0e[i, j] = Gu["M0"][ui, uj] + (r * PM if "M0" in groups else 0)
                r = (W1q[ui] * sq[j] + lamU * W1m[ui] * sm[j] + W1q[uj] * sq[i] + lamU * W1m[uj] * sm[i]
                     + b1 * sq[i] * sq[j] + lamU * b1r * sm[i] * sm[j])
                Ae[i, j] = B1u[ui, uj] + (r * PM if "A" in groups else 0)
                r2 = (W2q[ui] * sq[j] + lamU * W2m[ui] * sm[j] + W2q[uj] * sq[i] + lamU * W2m[uj] * sm[i]
                      + b2 * sq[i] * sq[j] + lamU * b2r * sm[i] * sm[j])
                r1 = (2 * (W1q[pi_] * sq[j] + lamU * W1m[pi_] * sm[j])
                      + 2 * (W1q[pj] * sq[i] + lamU * W1m[pj] * sm[i]))
                Xe[i, j] = (B2u[ui, uj] + 2 * B1u[pi_, uj] + 2 * B1u[ui, pj] + 2 * Su[pi_, pj]
                            + (r2 * PM if "X2" in groups else 0) + (r1 * PM if "X1" in groups else 0))
                r = (W1m[ui] * sm[j] + W1m[uj] * sm[i] + b1r * sm[i] * sm[j] + mps[pi_] * sm[j] + mps[pj] * sm[i])
                M1e[i, j] = (Gu["M1"][ui, uj] + Gu["M0"][pi_, uj] + Gu["M0"][ui, pj]
                             + (r * PM if "M1" in groups else 0))
            c1e[i, 0] = (Gu["M1"][U_[i], ONE] + Gu["M0"][PS[i], ONE]
                         + (W1m[ONE] * sm[i] * PM if "c1" in groups else 0))
        m11 = Gu["M0"][ONE, ONE]
        try:
            N = M0e.inv()
        except ZeroDivisionError:                    # M0e not provably invertible (e.g. NaN): never a pass
            nan = arb("nan")
            return nan, dict(trNX=nan, A_term=nan, c1_term=nan)
        parts = dict(trNX=(N * Xe).trace(), A_term=-2 * (N * M1e * N * Ae).trace(),
                     c1_term=2 * lam * (c1e.transpose() * N * c1e)[0, 0] / m11)
        return parts["trNX"] + parts["A_term"] + parts["c1_term"], parts
    ALL = {"M0", "A", "X2", "X1", "M1", "c1"}
    g2, parts = build(ALL)
    g2_0, _ = build(set())
    ok = certainly_lt(g2, 0)                           # NaN / infinite never pass
    g2_upper = _float_up(g2) if g2.is_finite() else None
    breakdown = {grp: json_float(_ub(build({grp})[0]) - float(g2_0.mid())) for grp in sorted(ALL)}
    # raw inputs keep all 30 mid digits (short ball_str folds the mid rounding into a much wider radius)
    bs = ball_str
    vec = lambda xs: [bs(x) for x in xs]
    raw = {"lam": [lam1_in.lower, lam1_in.upper], "lam4": gap_in.lower, "lamU": repr(_ub(lam)),
           "G": {k: [[bs(Gu[k][a, b]) for b in range(nv)] for a in range(nv)]
                 for k in ("K0", "K1", "K2", "M0", "M1", "M2")},
           "Wdiag": {k: vec(Gu[k][f, f] for f in range(nv)) for k in ("W1q", "W1m", "W2q", "W2m")},
           "W": {"W1q": vec(W1q), "W1m": vec(W1m), "W2q": vec(W2q), "W2m": vec(W2m), "mps": vec(mps)},
           "bmax": vec(bmax), "sq": vec(sq), "sm": vec(sm)}
    data = {
        "assembly": "g2_upper_v1",
        "raw": raw,
        "terms": {k: ball_str(v) for k, v in parts.items()},
        "g2_enclosure": ball_str(g2),
        "diagnostics": {
            "non_additive": True,
            "g2_central": ball_str(g2_0),
            "error_breakdown_upper_minus_central": breakdown,
            "beta_max": {k: json_float(_ub(b)) for k, b in zip(("beta1", "beta1r", "beta2", "beta2r"), bmax)},
            "eta": eta, "area_over_4pi": ball_str(Gu["M0"][ONE, ONE] / (4 * arb.pi())),
            "float_cluster": [json_float(v) for v in w[1:4]], "float_g2": json_float(g2f),
            "float_T2": json_float(T2f), "float_A_eigs": [json_float(v) for v in A_float],
            "A_e_enclosure_diag": [ball_str(B1u[i, i]) for i in range(3)],
            "B2_uu_diag": [ball_str(B2u[i, i]) for i in range(3)],
        },
    }
    cert = Certificate(kind="kyfan_second_order", inputs=inputs,
                       claims=[Claim(lhs="g2_upper", op="<", rhs="0", holds=ok)],
                       data=data, dependencies=deps, producer=make_producer())
    return A_float, g2_upper, ok, cert


def certify_cluster(direction="F_A", degree=1, n=96, lam1: Lam1Input = None, gap: GapInput = None,
                    second_order=False, sym=True, strict=False) -> ClusterResult:
    """Certified perturbation tests of the lambda_1 cluster at Bolza.

    lam1, gap: required explicit inputs (lambda_1 ball and lambda_4 lower bound), recorded with their sources.
    degree=1: the published P1 first-order test (A != 0) for direction "F_A"; see _first_order for the certified
    quantity and the certificate's b_ee_upper_v1 assembly formula.
    degree=2, second_order=True: the certified Ky Fan bound kyfan_upper >= D2 Lambda(v, v) along direction ("F_A",
    "F_B" or a 16-vector of Klein coordinates); see _kyfan for the g2_upper_v1 assembly formula.  sym=True (pi/2
    element cache) applies to the named, r2-invariant families only; it is forced False for vector directions.
    degree=2, second_order=False: float splitting rates A_float only (no certificate, certified=False).
    b_ee_upper / kyfan_upper are certificate-attested: floats rounded up from the stored enclosure the checker
    re-derives; the *_inprocess fields keep the (marginally tighter) in-memory values that the published goldens pin.
    Precision is set internally (96 bits)."""
    if lam1 is None or gap is None:
        raise ValueError("lam1 and gap inputs are required (spec §3.4)")
    if float(lam1.lower) > float(lam1.upper):
        raise ValueError(f"lam1 interval is empty: lower {lam1.lower} > upper {lam1.upper}")
    if float(gap.lower) <= float(lam1.upper):
        raise ValueError(f"gap lower bound {gap.lower} must exceed lam1 upper {lam1.upper}")
    if degree == 2:
        with prec(96):
            A_float, g2_upper, ok, cert = _kyfan(direction, n, lam1, gap, sym, second_order)
        if strict and second_order and not ok:
            raise CertificationFailed("could not certify D2 Lambda(v, v) < 0 (Ky Fan, second order)")
        attested = None if g2_upper is None else _attested_upper(cert.data["g2_enclosure"])
        return ClusterResult(A_float=A_float, b_ee_upper=None, A_nonzero=False, kyfan_upper=attested,
                             certified=ok, certificates=[cert] if cert is not None else [],
                             kyfan_upper_inprocess=g2_upper)
    if degree != 1:
        raise NotImplementedError("degree must be 1 or 2")
    if second_order:
        raise NotImplementedError("the second-order (Ky Fan) test needs degree=2")
    if direction != "F_A":
        raise NotImplementedError("first-order P1 test is published for F_A; use degree=2 via Task 9 for other "
                                  "directions")
    with prec(96):
        A_float, b_ee_upper, ok, cert = _first_order(n, lam1, gap, sym)
    if strict and not ok:
        raise CertificationFailed("could not certify b(e,e) < 0 (A != 0) at first order")
    attested = None if b_ee_upper is None else _attested_upper(cert.data["b_ee_upper_over_m"])
    return ClusterResult(A_float=A_float, b_ee_upper=attested, A_nonzero=ok, kyfan_upper=None, certified=ok,
                         certificates=[cert], b_ee_upper_inprocess=b_ee_upper)
