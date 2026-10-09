"""Certified Rayleigh-Ritz upper bounds for Laplace eigenvalues (P1 trial space, Bolza surface).

Ported from the Bolza lambda1 research code (cert_rr.py).  The float eigenvectors are untrusted; only their exact
dyadic values enter the certified part: lambda_k <= mu_max(Kp, Mp) by min-max, and mu_max < t is certified by
Sylvester's criterion for t*Mp - Kp in Arb.
"""
from __future__ import annotations

import time
from dataclasses import dataclass

import numpy as np
import scipy.linalg as sl
import scipy.sparse.linalg as sla
from flint import arb, arb_mat

from ..errors import CertificationFailed
from . import forms
from ._arb import ball_str, certainly_gt, prec, rational_upper
from .certificate import Certificate, Claim, Dependency, make_producer
from .mesh import Mesh, positions_arb


def is_pd(A):
    """Sylvester: all leading principal minors > 0 (certified in Arb)."""
    m = len(A)
    for k in range(1, m + 1):
        d = arb_mat([[A[i][j] for j in range(k)] for i in range(k)]).det()
        if not certainly_gt(d, 0):
            return False
    return True


def certified_upper(Kp, Mp, guess):
    """smallest t (by bisection from guess) with t*Mp - Kp certified positive definite => mu_max < t."""
    m = len(Kp)
    def ok(t):
        ta = arb(t)
        return is_pd([[ta * Mp[i][j] - Kp[i][j] for j in range(m)] for i in range(m)])
    step = 1e-6
    hi = guess * (1 + step)
    while not ok(hi):
        step *= 2
        hi = guess * (1 + step)
        if step > 1:
            raise RuntimeError("could not certify an upper bound")
    lo = guess * (1 - 1e-6)
    for _ in range(40):
        mid = (lo + hi) / 2
        if ok(mid): hi = mid
        else: lo = mid
    return hi


@dataclass
class SpectrumResult:
    upper: list
    float_eigs: list
    area_over_4pi: str
    certified: bool
    certificate: Certificate


def certify_spectrum(surface, k=3, degree=2, n=48, strict=False) -> SpectrumResult:
    """Certified upper bounds for lambda_1..lambda_k.  Only degree=1 on the Bolza surface is supported so far."""
    if degree != 1 or not surface.is_bolza:
        raise NotImplementedError("only degree=1 on the Bolza surface is supported in this version")
    with prec(96):
        mesh = Mesh(n)
        X = positions_arb(mesh, surface.vertices_hyp)
        elems, Kf, Mf, _ = forms.assemble_p1(mesh, X)
        w, U = sla.eigsh(Kf, k=k + 3, M=Mf, sigma=-0.5,
                         v0=np.random.RandomState(0).rand(mesh.Nf))  # fixed start: reproducible eigenvectors
        o = np.argsort(w); w, U = w[o], U[:, o]
        vecs = [np.ones(mesh.Nf)] + [U[:, i] for i in range(1, k + 1)]
        Kp, Mp = forms.project_p1(mesh, elems, vecs)
        Kpf = np.array([[float(x.mid()) for x in r] for r in Kp])
        Mpf = np.array([[float(x.mid()) for x in r] for r in Mp])
        mus = sl.eigh(Kpf, Mpf, eigvals_only=True)
        area = Mp[0][0] / (4 * arb.pi())
        m = k + 1
        data = {"K": [[ball_str(x) for x in r] for r in Kp], "M": [[ball_str(x) for x in r] for r in Mp]}
        upper, ok = [], True
        try:
            # a bound for each lambda_j, j <= k, from the leading (j+1)-dim sub-space {1, x1..xj}
            for j in range(1, k + 1):
                sub_mus = sl.eigh(Kpf[:j + 1, :j + 1], Mpf[:j + 1, :j + 1], eigvals_only=True)
                upper.append(certified_upper([r[:j + 1] for r in Kp[:j + 1]], [r[:j + 1] for r in Mp[:j + 1]],
                                             float(sub_mus[-1])))
            mu = rational_upper(arb(upper[-1]))
        except RuntimeError:
            ok = False
            mu = "inf"
        data["mu"] = mu
        claim = Claim(lhs=f"lambda_{k}", op="<=", rhs=mu, holds=ok)
        cert = Certificate(
            kind="rr_upper",
            inputs={"surface": surface.label, "degree": degree, "n": n, "k": k},
            claims=[claim], data=data,
            dependencies=[Dependency("D9", "min-max principle for the Laplacian on a closed surface",
                                     "Courant–Hilbert")],
            producer=make_producer())
        area_s = ball_str(area)
    if not ok and strict:
        raise CertificationFailed("could not certify a Rayleigh-Ritz upper bound")
    return SpectrumResult(upper=upper, float_eigs=[float(x) for x in mus],
                          area_over_4pi=area_s, certified=ok, certificate=cert)
