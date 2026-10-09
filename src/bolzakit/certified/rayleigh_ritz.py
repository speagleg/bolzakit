"""Certified Rayleigh-Ritz upper bounds for Laplace eigenvalues (P1 / P2 trial space, any certified Surface).

Ported from the Bolza lambda1 research code (cert_rr.py).  The float eigenvectors are untrusted; only their exact
dyadic values enter the certified part: lambda_k <= mu_max(Kp, Mp) by min-max, and mu_max < t is certified by
Sylvester's criterion for t*Mp - Kp in Arb.
"""
from __future__ import annotations

import math
from fractions import Fraction
import time
from dataclasses import dataclass, field

import numpy as np
import scipy.linalg as sl
import scipy.sparse.linalg as sla
from flint import arb, arb_mat

from ..errors import CertificationFailed
from . import forms
from ._arb import ball_str, certainly_gt, prec, rational_upper
from .certificate import Certificate, Claim, Dependency, make_producer, trial_provenance
from .jets import J
from .mesh import Mesh, positions_arb, positions_jet


def is_pd(A):
    """Sylvester: all leading principal minors > 0 (certified in Arb)."""
    m = len(A)
    for k in range(1, m + 1):
        d = arb_mat([[A[i][j] for j in range(k)] for i in range(k)]).det()
        if not certainly_gt(d, 0):
            return False
    return True


def certified_upper(Kp, Mp, guess, step0=1e-6, iters=40):
    """smallest t (by bisection from guess) with t*Mp - Kp certified positive definite => mu_max < t."""
    m = len(Kp)
    def ok(t):
        ta = arb(t)
        return is_pd([[ta * Mp[i][j] - Kp[i][j] for j in range(m)] for i in range(m)])
    step = step0
    hi = guess * (1 + step)
    while not ok(hi):
        step *= 2
        hi = guess * (1 + step)
        if step > 1:
            raise RuntimeError("could not certify an upper bound")
    lo = guess * (1 - step0)
    for _ in range(iters):
        mid = (lo + hi) / 2
        if ok(mid): hi = mid
        else: lo = mid
    return hi


def _stored_mu(data, t0, m=None):
    """Decimal mu >= t0 for which mu M - K (leading m x m block, default all) is certified positive definite from
    the balls as STORED (parsed at 128 bits), so that an independent checker reading the certificate can re-derive
    the claim.  The bisection point t0 is certifiable from the in-memory balls only; the stored balls are marginally
    wider, so mu starts 1e-12 (relative) above t0 and grows geometrically if needed."""
    with prec(128):
        K = [[arb(x) for x in r] for r in data["K"]]
        M = [[arb(x) for x in r] for r in data["M"]]
        m = len(K) if m is None else m
        rel = arb(10) ** -12
        for _ in range(30):
            mu = rational_upper(arb(t0) * (1 + rel))
            q = Fraction(mu)
            t = arb(q.numerator) / arb(q.denominator)
            if is_pd([[t * M[i][j] - K[i][j] for j in range(m)] for i in range(m)]):
                return mu
            rel *= 10
    raise RuntimeError("could not certify an upper bound from the stored balls")


def _float_at_least(q: str) -> float:
    """Smallest-ish float >= the exact decimal q."""
    f = float(Fraction(q))
    while Fraction(f) < Fraction(q):
        f = math.nextafter(f, math.inf)
    return f


@dataclass
class SpectrumResult:
    upper: list               # certificate-attested: upper[j-1] >= the stored-ball bound for lambda_j (upper[-1] >= mu)

    float_eigs: list
    area_over_4pi: str
    certified: bool
    certificate: Certificate
    cluster_sum_upper: float | None = None   # degree 2: upper bound of lambda_1+..+lambda_k (Ky Fan)
    upper_inprocess: list = field(default_factory=list)   # in-memory bisection values (<= upper; fidelity pins)


def _project_p2(mesh, surface, k):
    """P2 Rayleigh-Ritz data: exact Gram matrices on span{1, x_1..x_k} of float FE eigenvectors (untrusted)."""
    space = forms.P2Space(mesh)
    X = positions_jet(mesh, [[J(c) for c in V] for V in surface.vertices_hyp])
    # pi/2-rotation element cache: valid only for octagons invariant under V_k -> V_{k+2}
    elems = forms.compute_elements(mesh, X, order=0, sym=surface.r2_invariant)
    Kf = forms.assemble_float(space, elems, 0, 0)
    Mf = forms.assemble_float(space, elems, 0, 1)
    w, U = sla.eigsh(Kf, k=k + 3, M=Mf, sigma=-0.5,
                     v0=np.random.RandomState(0).rand(space.ndof))  # fixed start: reproducible eigenvectors
    o = np.argsort(w); U = U[:, o]
    vecs = [np.ones(space.ndof)] + [U[:, i] for i in range(1, k + 1)]
    GK, GM = forms.gram(space, forms.to_arbmats(elems, 1, 1), vecs, 0, 0)
    m = k + 1
    Kp = [[GK[i, j] for j in range(m)] for i in range(m)]
    Mp = [[GM[i, j] for j in range(m)] for i in range(m)]
    sample_src = [([c.v for c in X[i0]], [c.v for c in X[i1]], [c.v for c in X[i2]], e[0][0], e[1][0])
                  for (i0, i1, i2), e in zip(mesh.T, elems)]
    return Kp, Mp, GK, GM, vecs, sample_src


def element_sample(cert, elements):
    """Deterministic 1% sample (at least 1, at most 200) of the elements, for the checker's --deep consistency mode.

    elements[t] = (P0, P1, P2, Ke, Me): arb hyperboloid vertices and the element matrices.  The indices are
    random.Random(seed).sample(range(len(elements)), size) with seed = int(sha256[:8], 16), sha256 the hash of the
    certificate before element_sample is added (so a checker can re-derive the selection)."""
    import random
    size = max(1, min(200, len(elements) // 100))
    seed = int(cert.sha256()[:8], 16)
    idx = sorted(random.Random(seed).sample(range(len(elements)), size))
    mat = lambda A: [[ball_str(x) for x in row] for row in A]
    return {"n_elements": len(elements), "seed_rule": "int(sha256_without_element_sample[:8], 16)",
            "items": [{"index": t, "vertices": [[ball_str(c) for c in P] for P in elements[t][:3]],
                       "K": mat(elements[t][3]), "M": mat(elements[t][4])} for t in idx]}


def certify_spectrum(surface, k=3, degree=2, n=48, strict=False) -> SpectrumResult:
    """Certified upper bounds for lambda_1..lambda_k of a Surface (degree 1 or 2 trial space)."""
    if degree not in (1, 2):
        raise NotImplementedError("only degree=1 or 2 is supported in this version")
    cluster = None
    with prec(96):
        mesh = Mesh(n)
        if degree == 1:
            X = positions_arb(mesh, surface.vertices_hyp)
            elems, Kf, Mf, _ = forms.assemble_p1(mesh, X)
            w, U = sla.eigsh(Kf, k=k + 3, M=Mf, sigma=-0.5,
                             v0=np.random.RandomState(0).rand(mesh.Nf))  # fixed start: reproducible eigenvectors
            o = np.argsort(w); w, U = w[o], U[:, o]
            vecs = [np.ones(mesh.Nf)] + [U[:, i] for i in range(1, k + 1)]
            Kp, Mp = forms.project_p1(mesh, elems, vecs)
            sample_src = [(X[i0], X[i1], X[i2], e[0], e[1]) for (i0, i1, i2), e in zip(mesh.T, elems)]
            step0, iters = 1e-6, 40
            v0_rule = "RandomState(0).rand(Nf)"
        else:
            Kp, Mp, GK, GM, vecs, sample_src = _project_p2(mesh, surface, k)
            step0, iters = 1e-9, 50
            v0_rule = "RandomState(0).rand(ndof)"
        Kpf = np.array([[float(x.mid()) for x in r] for r in Kp])
        Mpf = np.array([[float(x.mid()) for x in r] for r in Mp])
        mus = sl.eigh(Kpf, Mpf, eigvals_only=True)
        area = Mp[0][0] / (4 * arb.pi())
        data = {"K": [[ball_str(x) for x in r] for r in Kp], "M": [[ball_str(x) for x in r] for r in Mp]}
        upper, upper_inprocess, ok = [], [], True
        try:
            # a bound for each lambda_j, j <= k, from the leading (j+1)-dim sub-space {1, x1..xj}; the bisection
            # value (in-memory balls) is re-certified from the stored balls, and upper reports that attested bound
            for j in range(1, k + 1):
                sub_mus = sl.eigh(Kpf[:j + 1, :j + 1], Mpf[:j + 1, :j + 1], eigvals_only=True)
                t0 = certified_upper([r[:j + 1] for r in Kp[:j + 1]], [r[:j + 1] for r in Mp[:j + 1]],
                                     float(sub_mus[-1]), step0, iters)
                upper_inprocess.append(t0)
                mu = _stored_mu(data, t0, j + 1)
                upper.append(_float_at_least(mu))
        except RuntimeError:
            ok = False
            mu = "inf"
        data["mu"] = mu
        claims_extra = []
        if degree == 2:
            # Ky Fan: lambda_1 + .. + lambda_k <= tr(M^-1 K) on the (k+1)-dim space (lambda_0 = 0 is exact)
            cs_rhs, cs_ok = "inf", False
            try:
                kyfan = GM.solve(GK).trace()
                # the claim is about the balls actually stored: a re-solve from the parsed K / M strings at 128 bits
                # (what an independent checker does) is joined to the in-memory enclosure
                with prec(128):
                    stored = arb_mat([[arb(x) for x in r] for r in data["M"]]).solve(
                        arb_mat([[arb(x) for x in r] for r in data["K"]])).trace()
                kyfan = kyfan.union(stored)
                if ok and kyfan.is_finite():
                    # outward margin 1e-6 of the radius: a checker on another flint version may enclose the trace
                    # of the stored balls slightly differently
                    cs_rhs = rational_upper(kyfan + kyfan.rad() * arb("1e-6"))
                    cluster, cs_ok = _float_at_least(cs_rhs), True   # never tighter than the claim rhs
                    data["cluster_sum"] = ball_str(kyfan)
            except ZeroDivisionError:  # GM not provably invertible
                pass
            claims_extra.append(Claim(lhs="+".join(f"lambda_{i}" for i in range(1, k + 1)), op="<=",
                                      rhs=cs_rhs, holds=cs_ok))
        claim = Claim(lhs=f"lambda_{k}", op="<=", rhs=mu, holds=ok)
        cert = Certificate(
            kind="rr_upper",
            inputs={"surface": surface.label, "degree": degree, "n": n, "k": k,
                    "trial": trial_provenance(vecs, v0_rule)},
            claims=[claim] + claims_extra, data=data,
            dependencies=[Dependency("D9", "min-max principle for the Laplacian on a closed surface",
                                     "Courant–Hilbert")],
            producer=make_producer())
        cert.data["element_sample"] = element_sample(cert, sample_src)
        area_s = ball_str(area)
    if not ok and strict:
        raise CertificationFailed("could not certify a Rayleigh-Ritz upper bound")
    return SpectrumResult(upper=upper, float_eigs=[float(x) for x in mus],
                          area_over_4pi=area_s, certified=ok, certificate=cert, cluster_sum_upper=cluster,
                          upper_inprocess=upper_inprocess)
