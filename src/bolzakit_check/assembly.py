"""Re-assembly of the cluster perturbation bounds from their stored inputs (PROOF.md Lemma 5, §4.4 gap transfer;
Lemma 8, step 5).  The Gram / transfer inputs are trusted data of the certificate; the final bound is recomputed here
and never read from the stored sum.

b_ee_upper_v1 (kind cluster_first_order), with lam1 = [lam1.lower, lam1.upper] and lam4 = gap.lower from the inputs:
    mw = (mu - lam1)/(lam4 - lam1),   qw = (mu - lam1) lam4/(lam4 - lam1)       (gap transfer, w = u - P_E u)
    E  = 2 (sqrt(Sa/m) sqrt(qw) + lam1 sqrt(Sr/m) sqrt(mw)) + beta_a qw + lam1 beta_rho mw
    claim: B + E < 0,  B = b(u,u)/m(u).
g2_upper_v1 (kind kyfan_second_order): Ky Fan second-order bound g''(0), rebuilt from data["raw"] (see g2()).
"""
from __future__ import annotations

from flint import arb, arb_mat

from ._cert import CheckError, ball, lt, precision

PM = arb(0, 1)
PREC = 96


def _lam_inputs(inputs):
    lam1, gap = inputs["lam1"], inputs["gap"]
    lam = ball(lam1["lower"]).union(ball(lam1["upper"]))
    lam4 = ball(gap["lower"])
    if not lt(lam, lam4):
        raise CheckError("gap lower bound does not exceed lam1")
    return lam, lam4


def b_ee(data, inputs):
    """Recomputed B + E (ball) for the first-order certificate."""
    if data.get("assembly") != "b_ee_upper_v1":
        raise CheckError("unknown assembly id for cluster_first_order")
    with precision(PREC):
        lam, lam4 = _lam_inputs(inputs)
        t = data["transfer_terms"]
        B = ball(data["terms"]["b_uu_over_m"])
        mu, Sa, Sr = ball(t["mu"]), ball(t["Sa_over_m"]), ball(t["Sr_over_m"])
        ba, br = ball(t["beta_a"]), ball(t["beta_rho"])
        exc = mu - lam
        mw = exc / (lam4 - lam)
        qw = exc * lam4 / (lam4 - lam)
        E = 2 * (Sa.sqrt() * qw.sqrt() + lam * Sr.sqrt() * mw.sqrt()) + ba * qw + lam * br * mw
        return B + E


def g2(data, inputs):
    """g2 rebuilt from data["raw"], pure Arb at 96 bits (formula: bolzakit cluster._kyfan docstring, Lemma 8)."""
    if data.get("assembly") != "g2_upper_v1":
        raise CheckError("unknown assembly id for kyfan_second_order")
    raw = data["raw"]
    with precision(PREC):
        if raw["lam"] != [inputs["lam1"]["lower"], inputs["lam1"]["upper"]] or raw["lam4"] != inputs["gap"]["lower"]:
            raise CheckError("raw lam / lam4 differ from the recorded lam1 / gap inputs")
        lam, lam4 = _lam_inputs(inputs)
        lamU = arb(float(raw["lamU"]))          # exact float point
        if not (lamU.is_finite() and lamU >= lam.upper()):
            raise CheckError("lamU is not >= the upper end of lam")
        G = {k: [[ball(x) for x in r] for r in raw["G"][k]] for k in ("K0", "K1", "K2", "M0", "M1", "M2")}
        if any(len(m) != 7 or any(len(r) != 7 for r in m) for m in G.values()):
            raise CheckError("raw Gram matrices must be 7 x 7")
        Wd = {k: [ball(x) for x in raw["Wdiag"][k]] for k in ("W1q", "W1m", "W2q", "W2m")}
        if any(len(w) != 7 for w in Wd.values()) or len(raw["bmax"]) != 4:
            raise CheckError("raw Wdiag / bmax have the wrong length")
        b1, b1r, b2, b2r = (ball(x) for x in raw["bmax"])
        U, P, ONE = (0, 1, 2), (3, 4, 5), 6
        for i in U:   # u_i is exactly m-orthogonal to the constants
            if not G["M0"][i][ONE].contains(0):
                raise CheckError("raw G[M0] does not enclose m(u_i, 1) = 0")
        sqrt_up = lambda x: arb(x.upper()).sqrt()  # noqa: E731  (nan for a negative upper end: never a pass)
        W = {k: [sqrt_up(x) for x in w] for k, w in Wd.items()}
        mps = [sqrt_up(G["M0"][f][f]) for f in range(7)]
        sm, sq = [], []
        for i in U:
            exc = G["K0"][i][i] - lam * G["M0"][i][i]
            sm.append(sqrt_up(exc / (lam4 - lam)))
            sq.append(sqrt_up(exc * lam4 / (lam4 - lam)))
        B1 = lambda a, b: G["K1"][a][b] - lam * G["M1"][a][b]  # noqa: E731
        B2 = lambda a, b: 2 * (G["K2"][a][b] - lam * G["M2"][a][b])  # noqa: E731
        S = lambda a, b: G["K0"][a][b] - lam * G["M0"][a][b]  # noqa: E731
        W1q, W1m, W2q, W2m = W["W1q"], W["W1m"], W["W2q"], W["W2m"]
        M0e, Ae, Xe, M1e = (arb_mat(3, 3) for _ in range(4))
        c1e = arb_mat(3, 1)
        for i in range(3):
            for j in range(3):
                u, v, p, q = U[i], U[j], P[i], P[j]
                M0e[i, j] = G["M0"][u][v] + sm[i] * sm[j] * PM
                Ae[i, j] = B1(u, v) + (W1q[u] * sq[j] + lamU * W1m[u] * sm[j] + W1q[v] * sq[i]
                                       + lamU * W1m[v] * sm[i] + b1 * sq[i] * sq[j] + lamU * b1r * sm[i] * sm[j]) * PM
                Xe[i, j] = (B2(u, v) + 2 * B1(p, v) + 2 * B1(u, q) + 2 * S(p, q)
                            + (W2q[u] * sq[j] + lamU * W2m[u] * sm[j] + W2q[v] * sq[i] + lamU * W2m[v] * sm[i]
                               + b2 * sq[i] * sq[j] + lamU * b2r * sm[i] * sm[j]) * PM
                            + (2 * (W1q[p] * sq[j] + lamU * W1m[p] * sm[j])
                               + 2 * (W1q[q] * sq[i] + lamU * W1m[q] * sm[i])) * PM)
                M1e[i, j] = (G["M1"][u][v] + G["M0"][p][v] + G["M0"][u][q]
                             + (W1m[u] * sm[j] + W1m[v] * sm[i] + b1r * sm[i] * sm[j]
                                + mps[p] * sm[j] + mps[q] * sm[i]) * PM)
            c1e[i, 0] = G["M1"][U[i]][ONE] + G["M0"][P[i]][ONE] + W1m[ONE] * sm[i] * PM
        m11 = G["M0"][ONE][ONE]
        try:
            N = M0e.inv()
        except ZeroDivisionError:
            return arb("nan")
        trNX = (N * Xe).trace()
        A_term = -2 * (N * M1e * N * Ae).trace()
        c1_term = 2 * lam * (c1e.transpose() * N * c1e)[0, 0] / m11
        return trNX + A_term + c1_term
