"""Bolza-surface ring dynamics: port of the original `dynamics.py` (mode-m energy-momentum Hessian of an N-ring centred at
the Weierstrass point z0 = 0, FULL surface Green function).

H_Sigma = sum_{i<j} G1(d_ij) + S,  S = sum_{i<j} R(z_i,z_j) + (1/2) sum_i R(z_i,z_i)   (Boatto-Koiller (3.8), our sign).
H_tanh part: ringmodes.py (mpmath, exact model).  S part: double precision, central differences in (s, theta) with
step h = 2e-3 (S is float64, so the mpmath step 10^(-dps/3) of ringmodes.analyse would be rounding-dominated).
Restricted to the 4-dim mode-m subspace {s cos, s sin, theta cos, theta sin}(m theta_k).

bolzakit edits relative to the original (porting rule): imports; N is a parameter (the original hard-codes the module constant
N = 7); `mp.mp.dps = 30` (global side effect) became `mp.workdps(dps)`; main() dropped; `lam_sigma` added -- it
evaluates only the (0,0) and (1,1) entries of the same restricted Hessian (the cos and sin copies `lam_sigma`,
`lam_sigma_sin` of `analyse`), which needs 4 S-evaluations per mode instead of 33. Formulas are unchanged.
"""
from __future__ import annotations
import math
import numpy as np
import mpmath as mp
from . import ringmodes as RM

CN = 1 / (2 * math.pi)


def jmom_derivs(s0):
    r2 = math.exp(2 * s0)
    j1 = 4 * r2 / (1 - r2) ** 2            # kappa
    j2 = 8 * r2 * (1 + r2) / (1 - r2) ** 3
    return j1, j2


def mode_basis(m, N=7):
    th = 2 * np.pi * np.arange(N) / N
    c, s = np.cos(m * th), np.sin(m * th)
    Z = np.zeros(N)
    B = [np.r_[c, Z], np.r_[s, Z], np.r_[Z, c], np.r_[Z, s]]
    return [v / np.linalg.norm(v) for v in B]


def zs_of(X, N=7):
    return np.exp(X[:N] + 1j * X[N:])


def restricted_hessian(f, X0, basis, h):
    k = len(basis)
    Hm = np.zeros((k, k))
    f0 = f(X0)
    for i in range(k):
        Hm[i, i] = (f(X0 + h * basis[i]) - 2 * f0 + f(X0 - h * basis[i])) / h ** 2
        for j in range(i + 1, k):
            v = (f(X0 + h * (basis[i] + basis[j])) - f(X0 + h * (basis[i] - basis[j]))
                 - f(X0 - h * (basis[i] - basis[j])) + f(X0 - h * (basis[i] + basis[j]))) / (4 * h ** 2)
            Hm[i, j] = Hm[j, i] = v
    return Hm


def analyse(g, xi, m=3, h=2e-3, N=7, dps=30):
    """The original `analyse`, verbatim up to the N / dps parameters. Requires 2m != N (else the sin basis vanishes)."""
    s0 = 0.5 * math.log(xi)
    X0 = np.r_[np.full(N, s0), 2 * np.pi * np.arange(N) / N]
    fS = lambda X: g.S_ring(zs_of(X, N))
    # dS/ds (uniform dilation) -> Omega correction
    dil = np.r_[np.ones(N), np.zeros(N)]
    dS = (fS(X0 + h * dil) - fS(X0 - h * dil)) / (2 * h) / N
    # tanh model (mpmath)
    with mp.workdps(dps):
        rt = RM.analyse(RM.H_tanh, N, mp.mpf(xi), modes=False)
    Ht = np.array(rt['hess'].tolist(), dtype=float)            # already includes -Omega_t j''
    j1, j2 = jmom_derivs(s0)
    Omega_t = float(rt['Omega'])
    dOmega = dS / j1
    B = mode_basis(m, N)
    Ht_r = np.array([[bi @ Ht @ bj for bj in B] for bi in B])
    HS_r = restricted_hessian(fS, X0, B, h)
    corr = np.zeros((4, 4)); corr[0, 0] = corr[1, 1] = -dOmega * j2
    Hsig = Ht_r + HS_r + corr
    Jm = np.array([[0, 0, 1, 0], [0, 0, 0, 1], [-1, 0, 0, 0], [0, -1, 0, 0]]) / j1
    ev_t = np.linalg.eigvals(Jm @ Ht_r)
    ev_s = np.linalg.eigvals(Jm @ Hsig)
    return dict(xi=xi, lam_tanh=Ht_r[0, 0] / CN, lam_sigma=Hsig[0, 0] / CN, lam_sigma_sin=Hsig[1, 1] / CN,
                offdiag=float(np.max(np.abs(Hsig - np.diag(np.diag(Hsig))))) / CN,
                omega2_tanh=float(np.max(-(ev_t ** 2).real)), omega2_sigma=float(np.max(-(ev_s ** 2).real)),
                maxRe_sigma=float(np.max(np.abs(ev_s.real))), Omega_t=Omega_t, dOmega=dOmega)


def lam_sigma(g, xi, ms, h=2e-3, N=7, dps=30):
    """{m: (lam_sigma, lam_sigma_sin)} for m in ms: exactly `analyse(g, xi, m, h, N, dps)['lam_sigma']` and
    `['lam_sigma_sin']` (Hsig[0,0] / CN cos copy, Hsig[1,1] / CN sin copy), computing only those diagonal entries.
    When 2m = N the sin vector vanishes and the cos value is returned for both.
    Cost: 3 + 4 len(ms) S-evaluations (f0 shared, dilation pair shared)."""
    s0 = 0.5 * math.log(xi)
    X0 = np.r_[np.full(N, s0), 2 * np.pi * np.arange(N) / N]
    fS = lambda X: g.S_ring(zs_of(X, N))
    dil = np.r_[np.ones(N), np.zeros(N)]
    dS = (fS(X0 + h * dil) - fS(X0 - h * dil)) / (2 * h) / N
    with mp.workdps(dps):
        rt = RM.analyse(RM.H_tanh, N, mp.mpf(xi), modes=False)
    Ht = np.array(rt['hess'].tolist(), dtype=float)
    j1, j2 = jmom_derivs(s0)
    dOmega = dS / j1
    f0 = fS(X0)
    th = 2 * np.pi * np.arange(N) / N

    def diag(v):                                   # = Hsig[i,i] / CN for the normalised radial vector v
        b = np.r_[v, np.zeros(N)]
        b = b / np.linalg.norm(b)
        HS = (fS(X0 + h * b) - 2 * f0 + fS(X0 - h * b)) / h ** 2
        return (b @ Ht @ b + HS - dOmega * j2) / CN

    out = {}
    for m in ms:
        c = diag(np.cos(m * th))                   # = mode_basis(m, N)[0]
        sn = np.sin(m * th)
        out[m] = (c, diag(sn) if np.linalg.norm(sn) > 1e-12 else c)    # = mode_basis(m, N)[1]
    return out
