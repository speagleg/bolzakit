"""Energy-momentum (Lagrangian) Hessian and linearised dynamics of a regular N-ring of identical unit vortices
centred at 0 in the Poincare disk (or the plane), for an arbitrary rotation-invariant Hamiltonian H(z_1..z_N).

Coordinates: z_k = exp(s_k + i theta_k).  Symplectic form omega = sum kappa(s_k) ds_k ^ dtheta_k with
kappa(s) = rho^2 r^2,  rho^2 = 4/(1-r^2)^2 (hyperbolic area form, curvature -1) or 1 (flat).
Rotation momentum j(s) with j' = kappa:  j = 2 e^{2s}/(1-e^{2s}) (hyp), e^{2s}/2 (flat).
Ring: s_k = s0, theta_k = 2 pi k/N; Omega fixed by dH/ds_k = Omega j'(s0) (relative-equilibrium condition).
F = H - Omega J.  Radial mode eigenvalue  alpha_m = v^T Hess(F) v / v^T v,  v_k = cos(m theta_k) in s-directions.
"C1 units": lambda_m := alpha_m / c_norm with c_norm fixed by the FLAT model so that lambda_m^flat = (N-1) - m(N-m)/2.
Linearised dynamics: d/dt (ds, dth) = kappa(s0)^{-1} [[0, 1], [-1, 0]] Hess(F) (ds, dth); eigenvalues returned.
All in mpmath (precision set by caller); Hessian by central differences with step h (error O(h^2)).
"""
from __future__ import annotations
import mpmath as mp


def kappa(s, flat=False):
    r2 = mp.e ** (2 * s)
    return r2 if flat else 4 * r2 / (1 - r2) ** 2


def jmom(s, flat=False):
    r2 = mp.e ** (2 * s)
    return r2 / 2 if flat else 2 * r2 / (1 - r2)


def ring_point(N, xi):
    s0 = mp.log(xi) / 2       # xi = r^2  (Euclidean disk radius squared = tanh^2(r_hyp/2))
    return [s0] * N, [2 * mp.pi * k / N for k in range(N)]


def _H_of(Hz, S, TH):
    return Hz([mp.e ** (s + 1j * t) for s, t in zip(S, TH)])


def analyse(Hz, N, xi, flat=False, h=None, modes=True):
    S, TH = ring_point(N, xi)
    X0 = S + TH
    n = 2 * N
    h = h or mp.mpf(10) ** (-(mp.mp.dps // 3))
    f = lambda X: _H_of(Hz, X[:N], X[N:])
    # gradient (for Omega)
    def grad_i(i):
        Xp = list(X0); Xm = list(X0); Xp[i] += h; Xm[i] -= h
        return (f(Xp) - f(Xm)) / (2 * h)
    dHds = [grad_i(i) for i in range(N)]
    dHdt = [grad_i(N + i) for i in range(N)]
    s0 = S[0]
    jp = mp.diff(lambda s: jmom(s, flat), s0)
    Omega = sum(dHds) / N / jp
    jpp = mp.diff(lambda s: jmom(s, flat), s0, 2)
    # Hessian of H
    Hs = mp.matrix(n, n)
    f0 = f(X0)
    for i in range(n):
        for j in range(i, n):
            if i == j:
                Xp = list(X0); Xm = list(X0); Xp[i] += h; Xm[i] -= h
                Hs[i, i] = (f(Xp) - 2 * f0 + f(Xm)) / h ** 2
            else:
                vals = []
                for si, sj in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
                    X = list(X0); X[i] += si * h; X[j] += sj * h
                    vals.append(si * sj * f(X))
                Hs[i, j] = Hs[j, i] = sum(vals) / (4 * h ** 2)
    for i in range(N):
        Hs[i, i] -= Omega * jpp
    res = dict(Omega=Omega, resid_s=max(abs(x - sum(dHds) / N) for x in dHds), resid_t=max(abs(x) for x in dHdt), hess=Hs)
    if modes:
        alpha = {}
        for m in range(1, N):
            v = [mp.cos(m * TH[k]) for k in range(N)] + [0] * N
            num = sum(v[i] * Hs[i, j] * v[j] for i in range(n) for j in range(n))
            alpha[m] = num / sum(x * x for x in v)
        res['alpha'] = alpha
        # tangential
        beta = {}
        for m in range(1, N):
            v = [0] * N + [mp.sin(m * TH[k]) for k in range(N)]
            num = sum(v[i] * Hs[i, j] * v[j] for i in range(n) for j in range(n))
            beta[m] = num / sum(x * x for x in v)
        res['beta'] = beta
    # linearisation
    k0 = kappa(s0, flat)
    Jm = mp.matrix(n, n)
    for i in range(N):
        Jm[i, N + i] = 1 / k0
        Jm[N + i, i] = -1 / k0
    L = Jm * Hs
    ev = mp.eig(L, left=False, right=False)
    res['eig'] = sorted(ev, key=lambda z: (abs(z), mp.im(z)))
    return res


# ---------------- model Hamiltonians -----------------
def H_flat(zs):
    N = len(zs)
    return -sum(mp.log(abs(zs[i] - zs[j]) ** 2) for i in range(N) for j in range(i + 1, N)) / (4 * mp.pi)


def H_tanh(zs):
    """sum_{i<j} G1(d_ij), G1 = -(1/2pi) ln tanh(d/2) = -(1/4pi)[ln|z-w|^2 - ln|1 - z wbar|^2]."""
    N = len(zs)
    return -sum(mp.log(abs(zs[i] - zs[j]) ** 2) - mp.log(abs(1 - zs[i] * mp.conj(zs[j])) ** 2)
                for i in range(N) for j in range(i + 1, N)) / (4 * mp.pi)


def H_sinh(zs):
    """sum_{i<j} -(1/2pi) ln sinh(d_ij/2) = -(1/4pi)[ln|z-w|^2 - ln(1-|z|^2) - ln(1-|w|^2)] (Montaldi-Tokieda (1.12))."""
    N = len(zs)
    return -sum(mp.log(abs(zs[i] - zs[j]) ** 2) - mp.log(1 - abs(zs[i]) ** 2) - mp.log(1 - abs(zs[j]) ** 2)
                for i in range(N) for j in range(i + 1, N)) / (4 * mp.pi)


def c_norm(N):
    """alpha_1^flat / ((N-1)/2): fixes C1 units."""
    r = analyse(H_flat, N, mp.mpf('0.04'), flat=True)
    return r['alpha'][1] / (mp.mpf(N - 1) - mp.mpf(N - 1) / 2)
