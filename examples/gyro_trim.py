"""Frequency trimming of ring-gyroscope resonators with bolzakit.sensitivity.trim (NUMERICAL, discrete problem).

A vibratory ring gyroscope needs its two n = 2 wineglass modes at the same frequency; fabrication imperfections split
them, and the split is trimmed by adding (or ablating) small masses at sites on the rim. This example computes a
minimal-mass trim from the FE model and checks it with an independent dense eigensolve.

(a) Free-free thin ring (mean radius 1, width 0.05, plane stress, nu = 0.3, E = rho0 = 1, Quad2 elements, 96 cells
    around, 2 across) with a 4th-harmonic density imperfection rho0 * eps * cos(4(theta - phi)), eps = 1e-3,
    phi = 0.3. Thin-ring theory: the pair (3, 4) (indices 0-2 are rigid modes) splits at rates +-(3/10) lam_bar eps
    (spec 5.10, correcting research note 02). Measured here: 0.29984.
(b) The same ring on 8 tapered spokes clamped at a hub (r = 0.3): a geometry with no closed-form split rule. Random
    element-wise Young's modulus and density (1 + 0.005 xi, xi ~ N(0, 1)), five seeds. The n = 2 pair is (6, 7).

Trim sites: point masses (lumped, both displacement components) at outer-rim nodes, 16 sites for (a), 32 for (b).
Baseline: iterated antinode trimming. Each round adds, at the site where the higher mode moves most (largest
first-order split reduction per unit mass), the mass that cancels the split to first order along that one site, then
re-solves (10 rounds). It ignores the coupling term, which is why it can oscillate.
"""
from __future__ import annotations

import time

import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp
from skfem import Basis, BilinearForm, ElementQuad2, ElementVector, MeshQuad, asm
from skfem.helpers import ddot, dot, sym_grad, trace

from bolzakit.sensitivity import smallest_eigenpairs, splitting_matrices
from bolzakit.sensitivity.trim import trim_degenerate_pair

NU = 0.3
N_THETA = 96


@BilinearForm
def _stiffness(u, v, w):
    E = w["E"]
    mu, lam = E / (2 * (1 + NU)), E * NU / (1 - NU ** 2)        # plane stress
    eu, ev = sym_grad(u), sym_grad(v)
    return lam * trace(eu) * trace(ev) + 2 * mu * ddot(eu, ev)


@BilinearForm
def _mass(u, v, w):
    return w["rho"] * dot(u, v)


def ring_mesh(width=0.05, n_theta=N_THETA, n_r=2, n_spokes=0, spoke_cells=2, r_hub=0.3, n_s=12):
    """Structured Quad mesh of a ring (node (a, b) = a * n_theta + b, a = 0 inner edge, a = n_r outer edge), with
    optional tapered spokes (spoke_cells ring cells wide, n_s cells long) from the inner edge to radius r_hub.
    Returns (mesh, hub_nodes)."""
    r = 1 + width * (np.arange(n_r + 1) / n_r - 0.5)
    th = 2 * np.pi * np.arange(n_theta) / n_theta
    R, TH = np.meshgrid(r, th, indexing="ij")
    P = [np.vstack([(R * np.cos(TH)).ravel(), (R * np.sin(TH)).ravel()])]
    a, b = (x.ravel() for x in np.meshgrid(np.arange(n_r), np.arange(n_theta), indexing="ij"))
    nid = lambda a_, b_: a_ * n_theta + (b_ % n_theta)
    T = [np.vstack([nid(a, b), nid(a + 1, b), nid(a + 1, b + 1), nid(a, b + 1)])]
    nnode, hub = (n_r + 1) * n_theta, []
    rs = r_hub + (r[0] - r_hub) * np.arange(n_s) / n_s
    for k in range(n_spokes):
        bs = k * (n_theta // n_spokes) - spoke_cells // 2 + np.arange(spoke_cells + 1)
        ids = np.empty((n_s + 1, spoke_cells + 1), dtype=int)
        for cc in range(n_s):
            for q, bb in enumerate(bs):
                ids[cc, q] = nnode
                P.append([[rs[cc] * np.cos(th[bb % n_theta])], [rs[cc] * np.sin(th[bb % n_theta])]])
                nnode += 1
        ids[n_s] = [nid(0, bb) for bb in bs]
        hub.extend(ids[0])
        for cc in range(n_s):
            for q in range(spoke_cells):
                T.append([[ids[cc, q]], [ids[cc + 1, q]], [ids[cc + 1, q + 1]], [ids[cc, q + 1]]])
    return MeshQuad(np.hstack(P), np.hstack(T)), np.array(hub, dtype=int)


class Model:
    """Assembled resonator: K(E), M(rho) with element-wise fields, Dirichlet hub dofs removed, rim trim sites."""

    def __init__(self, n_spokes=0, n_sites=16):
        self.mesh, hub = ring_mesh(n_spokes=n_spokes)
        self.basis = Basis(self.mesh, ElementVector(ElementQuad2()))
        self.ne, self.nqp = self.mesh.t.shape[1], self.basis.X.shape[-1]
        fixed = self.basis.nodal_dofs[:, hub].ravel()
        self.keep = np.setdiff1d(np.arange(self.basis.N), fixed)
        pos = np.full(self.basis.N, -1)
        pos[self.keep] = np.arange(len(self.keep))
        n = len(self.keep)
        outer = np.arange(2 * N_THETA, 3 * N_THETA)                  # a = n_r = 2: outer-rim nodes
        self.sites = outer[:: N_THETA // n_sites]
        self.dMs = [sp.csc_matrix((np.ones(2), (pos[self.basis.nodal_dofs[:, s]],) * 2), shape=(n, n))
                    for s in self.sites]
        X = self.basis.mapping.F(self.basis.X)                       # (2, elements, quadrature points)
        self.theta = np.arctan2(X[1], X[0])                          # at the quadrature points

    def _cond(self, X):
        return X.tocsr()[self.keep][:, self.keep].tocsc()

    def element_field(self, values):
        return np.repeat(np.asarray(values, float)[:, None], self.nqp, axis=1)

    def K(self, E=None):
        return self._cond(asm(_stiffness, self.basis, E=np.ones((self.ne, self.nqp)) if E is None else E))

    def M(self, rho=None):
        return self._cond(asm(_mass, self.basis, rho=np.ones((self.ne, self.nqp)) if rho is None else rho))


def naive_trim(K, M, dMs, pair, rounds=10):
    """Iterated antinode trimming (the baseline; see the module docstring). Returns (t, rel_splits, eigensolves)."""
    i, j = pair
    t = np.zeros(len(dMs))
    zero = [sp.csc_matrix(K.shape)] * len(dMs)
    rel = []
    for r in range(rounds + 1):
        Mt = M + sum(tl * d for tl, d in zip(t, dMs) if tl)
        w, V = smallest_eigenpairs(K, Mt, j + 2)
        g = w[j] - w[i]
        rel.append(g / w[i:j + 1].mean())
        if r == rounds:
            break
        A = splitting_matrices(K, Mt, zero, dMs, pair, eig=(w, V), validate=False).mats
        rate = A[:, 0, 0] - A[:, 1, 1]                               # split reduction per unit mass at each site
        s = int(np.argmax(rate))
        t[s] += g / rate[s]
    return t, np.array(rel), rounds + 1


def dense_split(K, M, pair):
    """Independent check: dense LAPACK generalized eigensolve (full spectrum, divide and conquer), relative split of
    the pair. (The subset driver gvx can return a resolved pair as an exact tie, so it is not used.)"""
    w = sla.eigh(K.toarray(), M.toarray(), eigvals_only=True, driver="gvd")[list(pair)]
    return (w[1] - w[0]) / w.mean()


def part_a():
    m = Model(n_spokes=0, n_sites=16)
    K, M0 = m.K(), m.M()
    dM4 = asm(_mass, m.basis, rho=np.cos(4 * (m.theta - 0.3)))
    dM4 = m._cond(dM4)
    pair = (3, 4)
    w, V = smallest_eigenpairs(K, M0, 6)
    sm = splitting_matrices(K, M0, [sp.csc_matrix(K.shape)], [dM4], pair, eig=(w, V))
    rate = float(np.abs(np.linalg.eigvalsh(sm.mats[0])).max() / sm.lam_bar)
    M = M0 + 1e-3 * dM4
    res = trim_degenerate_pair(K, M, None, m.dMs, pair)
    Mt = M + sum(tl * d for tl, d in zip(res.t, m.dMs) if tl)
    tn, reln, nn = naive_trim(K, M, m.dMs, pair)
    return dict(rate=rate, res=res, dense=dense_split(K, Mt, pair), naive_t=tn, naive_rel=reln, naive_eigs=nn,
                mass=float(M0.sum() / 2))


def part_b(seeds=range(5)):
    m = Model(n_spokes=8, n_sites=32)
    pair = (6, 7)
    w0, _ = smallest_eigenpairs(m.K(), m.M(), 9)
    out = []
    for seed in seeds:
        rng = np.random.default_rng(seed)
        K = m.K(m.element_field(1 + 5e-3 * rng.standard_normal(m.ne)))
        M = m.M(m.element_field(1 + 5e-3 * rng.standard_normal(m.ne)))
        t0 = time.perf_counter()
        res = trim_degenerate_pair(K, M, None, m.dMs, pair)
        sec = time.perf_counter() - t0
        Mt = M + sum(tl * d for tl, d in zip(res.t, m.dMs) if tl)
        tn, reln, nn = naive_trim(K, M, m.dMs, pair)
        out.append(dict(seed=seed, res=res, sec=sec, dense=dense_split(K, Mt, pair), naive_t=tn, naive_rel=reln,
                        naive_eigs=nn))
    return dict(nominal=w0, runs=out, mass=float(m.M().sum() / 2))


def main(quick: bool = False):
    t_start = time.perf_counter()
    a = part_a()
    r = a["res"]
    print(f"(a) ring, 4th-harmonic imperfection: rate/lam_bar = {a['rate']:.5f} (thin-ring theory 0.3)")
    print(f"    split/lam_bar {r.rel_split_before:.3e} -> {r.rel_split_after:.3e} ({r.status}, "
          f"{r.n_eigensolves} eigensolves); dense check {a['dense']:.2e}; trim mass {r.cost:.4e} "
          f"({r.cost / a['mass']:.2e} of the ring) at sites {np.flatnonzero(r.t).tolist()}")
    print(f"    history: {' '.join(f'{h.rel_split:.1e}' for h in r.history)}")
    print(f"    baseline (antinode, {a['naive_eigs']} eigensolves): split {a['naive_rel'][-1]:.3e}, "
          f"mass {a['naive_t'].sum():.4e} at sites {np.flatnonzero(a['naive_t']).tolist()}")
    b = part_b(range(2) if quick else range(5))
    print(f"\n(b) ring on 8 spokes, random E and rho (sigma 0.5 %); nominal pair lambda_6 = lambda_7 = "
          f"{b['nominal'][6]:.6f} (rel. diff {(b['nominal'][7] - b['nominal'][6]) / b['nominal'][6]:.1e})")
    print(" seed  split before   ours: after  dense check  eigsolves  mass (rel)     sites   |  "
          "antinode: after  mass (rel)  eigsolves")
    for x in b["runs"]:
        r = x["res"]
        print(f" {x['seed']:4d}  {r.rel_split_before:11.3e}  {r.rel_split_after:11.2e}  {x['dense']:11.2e}  "
              f"{r.n_eigensolves:9d}  {r.cost / b['mass']:9.2e}  {str(np.flatnonzero(r.t).tolist()):>9s}  |  "
              f"{x['naive_rel'][-1]:15.2e}  {x['naive_t'].sum() / b['mass']:10.2e}  {x['naive_eigs']:9d}")
    print(f"label: NUMERICAL (discrete problem); {time.perf_counter() - t_start:.1f} s")
    # self-checks (spec 10.1)
    assert abs(a["rate"] - 0.3) < 5e-4
    assert a["res"].rel_split_after <= 1e-3 * a["res"].rel_split_before and a["dense"] <= 1e-3 * r.rel_split_before
    for x in b["runs"]:
        r = x["res"]
        assert r.status == "converged" and r.rel_split_after <= 1e-3 * r.rel_split_before
        assert x["dense"] <= 1e-3 * r.rel_split_before
    return dict(a=a, b=b)


if __name__ == "__main__":
    main()
