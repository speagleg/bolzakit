"""Fundamental-frequency maximisation: bolzakit's repeated-eigenvalue steepest ascent against standard practice
(NUMERICAL, discrete problem). Reported as measured, whoever wins.

Problem. Kirchhoff plate on the unit square, simply supported on its edges and on a central point support (a post),
Morley elements on MeshTri.init_symmetric().refined(3) (512 dofs), nu = 0.3. Thickness h_z on a 5 x 5 grid of zones
(25 design variables): K = sum h_z^3 K_z, M = sum h_z M_z. Constraints: volume sum a_z h_z = sum a_z (equality),
0.5 <= h_z <= 1.5. Objective: maximise lambda_0. At the uniform design lambda_0 = 1.586 pi^4 is simple and
lambda_1 = lambda_2 = 1.834 pi^4; the optima found have lambda_0 double or triple, the case where methods that
treat eigenvalues as differentiable chatter.

Methods, identical problem and eigensolve budget (an eigensolve = one call of the core eigensolver; every method uses
it through bolzakit.sensitivity, counted by a wrapper installed only while this benchmark runs):
  ours           maximize_eigenvalue(LowestEigenvalue(0)) with its defaults
  SLSQP bound    max beta s.t. lambda_i(h) >= beta (i < 4), volume, bounds; scipy SLSQP with the simple-eigenvalue
                 (Hellmann-Feynman) gradients u_i^T (dK - lambda_i dM) u_i
  trust-constr   the same bound formulation with scipy trust-constr
  MMA bound      the same bound formulation with pyMOTO's MMA (pymoto.minimize_mma, defaults, volume as <=); run only
                 if pymoto is importable (a benchmark-only tool, never a bolzakit dependency)
  naive          projected gradient ascent on lambda_0 (gradient of the eigenvector eigh returned), fixed step 0.02
Every method reports the best feasible design it evaluated within the eigensolve budget (ours: its monotone best
iterate; the others: the best feasible evaluation, line-search trials included) and its value is re-solved. "Best
from a start" is the maximum of those finals over methods. Evaluations are cached by design (last 8), so no method
is charged twice for a point it already evaluated; every real eigensolve is counted.
Metrics per run: final lambda_0 / pi^4 at the reported design, its multiplicity (eigenvalues within 1e-4 relative),
eigensolves and wall time until the running best reaches 99 % and 99.9 % of the best value from that start, the
number of decreases between successive accepted iterates (ours: accepted steps; SLSQP / trust-constr: callback
iterates, where a merit-function method may legitimately decrease lambda_0; MMA and naive: every iterate, all
accepted), constraint residuals, and the stop reason.
Zones: the 5 x 5 grid of x- and y-intervals does not align with the mesh's 1/16 grid, so the zones are unequal
(areas 0.03125 to 0.046875, summing to 1); volume constraints and gradients use the zone areas a_z.
The bound formulation uses beta normalised by lambda_0(start), bounded in [0.5, 2], for all three solvers.
"""
from __future__ import annotations

import collections
import contextlib
import time
import warnings

import numpy as np
import scipy.optimize as so
from skfem import Basis, BilinearForm, ElementTriMorley, MeshTri, asm
from skfem.helpers import dd, ddot, trace

from bolzakit.sensitivity import core as _core
from bolzakit.sensitivity import smallest_eigenpairs
from bolzakit.sensitivity.optimize import LowestEigenvalue, ParametricPencil, maximize_eigenvalue

NU = 0.3
LO, HI = 0.5, 1.5
NZ = 5
BETA_LO, BETA_HI = 0.5, 2.0   # beta / lambda_0(start), the same box for SLSQP, trust-constr and MMA
MB = 4                      # bound formulation: lambda_i >= beta for i < MB
P4 = np.pi ** 4
METHODS = ("ours", "SLSQP bound", "trust-constr bound", "MMA bound (pyMOTO)", "naive fixed-step")


@BilinearForm
def _bending(u, v, w):
    return ((1 - NU) * ddot(dd(u), dd(v)) + NU * trace(dd(u)) * trace(dd(v))) / 12 * w["chi"]


@BilinearForm
def _mass(u, v, w):
    return u * v * w["chi"]


def build(refine: int = 3, nz: int = NZ):
    """Zone matrices with the edge supports and the central post removed, and the zone areas."""
    mesh = MeshTri.init_symmetric().refined(refine)
    basis = Basis(mesh, ElementTriMorley())
    c = mesh.p[:, mesh.t].mean(axis=1)
    zone = np.minimum((c[1] * nz).astype(int), nz - 1) * nz + np.minimum((c[0] * nz).astype(int), nz - 1)
    centre = int(np.argmin(np.linalg.norm(mesh.p - 0.5, axis=0)))
    drop = np.union1d(basis.get_dofs().nodal["u"], basis.nodal_dofs[0, [centre]])
    keep = np.setdiff1d(np.arange(basis.N), drop)
    nqp = basis.X.shape[-1]
    Ks, Ms = [], []
    for z in range(nz * nz):
        chi = np.repeat((zone == z).astype(float)[:, None], nqp, axis=1)
        Ks.append(asm(_bending, basis, chi=chi).tocsr()[keep][:, keep].tocsc())
        Ms.append(asm(_mass, basis, chi=chi).tocsr()[keep][:, keep].tocsc())
    p, t = mesh.p, mesh.t
    a, b, cc = p[:, t[0]], p[:, t[1]], p[:, t[2]]
    area = 0.5 * np.abs((b[0] - a[0]) * (cc[1] - a[1]) - (b[1] - a[1]) * (cc[0] - a[0]))
    return Ks, Ms, np.bincount(zone, weights=area, minlength=nz * nz)


class Benchmark:
    def __init__(self, budget: int = 200):
        self.Ks, self.Ms, self.areas = build()
        self.d, self.V, self.budget = len(self.Ks), float(self.areas.sum()), budget
        Ks, Ms = self.Ks, self.Ms
        self.pen = ParametricPencil(K=lambda h: sum(hz ** 3 * Kz for hz, Kz in zip(h, Ks)),
                                    M=lambda h: sum(hz * Mz for hz, Mz in zip(h, Ms)),
                                    dK=lambda h: [3 * hz ** 2 * Kz for hz, Kz in zip(h, Ks)],
                                    dM=lambda h: list(Ms))
        self.count, self._mma = 0, False
        self._cache = collections.OrderedDict()
        self._reset_best()

    def _reset_best(self):
        self.best_val, self.best_h, self.bestlog, self._t0 = -np.inf, None, [], time.perf_counter()

    def _feasible(self, h):
        tol_v = 1e-3 if self._mma else 1e-9
        return (float(np.max(LO - h)) <= 1e-12 and float(np.max(h - HI)) <= 1e-12
                and abs(self.areas @ h - self.V) / self.V <= tol_v)

    @contextlib.contextmanager
    def counting(self):
        orig = _core._eigs

        def wrapped(*a, **k):
            self.count += 1
            return orig(*a, **k)
        _core._eigs = wrapped
        try:
            yield
        finally:
            _core._eigs = orig

    def project(self, y):
        """Euclidean projection onto {areas . h = V, LO <= h <= HI}."""
        f = lambda nu: self.areas @ np.clip(y - nu * self.areas, LO, HI) - self.V
        return np.clip(y - so.brentq(f, -1e3, 1e3, xtol=1e-15) * self.areas, LO, HI)

    def eig(self, h):
        """(lambda_0..MB-1, simple-eigenvalue gradients) at h, one counted eigensolve per new h."""
        key = h.tobytes()
        if key in self._cache:
            self._cache.move_to_end(key)
        else:
            if self.count >= self.budget:
                raise _Budget
            w, U = smallest_eigenpairs(self.pen.K(h), self.pen.M(h), MB)
            G = np.array([[U[:, i] @ (3 * h[l] ** 2 * (self.Ks[l] @ U[:, i]) - w[i] * (self.Ms[l] @ U[:, i]))
                           for l in range(self.d)] for i in range(MB)])
            self._cache[key] = (w, G)
            if len(self._cache) > 8:
                self._cache.popitem(last=False)
            if w[0] > self.best_val and self._feasible(h):
                self.best_val, self.best_h = w[0], h.copy()
            self.bestlog.append((self.count, time.perf_counter() - self._t0, self.best_val))
        return self._cache[key]

    # ------------------------------------------------------------------ methods: (h0, trace) -> (h, stop reason)
    def ours(self, h0, trace):
        t0 = time.perf_counter()
        res = maximize_eigenvalue(self.pen, h0, LowestEigenvalue(0), lower=LO, upper=HI, A_eq=self.areas[None, :],
                                  b_eq=[self.V], max_iter=10_000, _max_eigensolves=self.budget,
                                  callback=lambda r: trace.append((self.count, time.perf_counter() - t0, r.value)))
        trace.append((self.count, time.perf_counter() - t0, res.value))
        return res.p, res.status

    def bound(self, h0, trace, method):
        t0, d = time.perf_counter(), self.d
        lam0 = self.eig(h0)[0][0]
        x0 = np.r_[h0, 1.0]
        fun = lambda x: self.eig(x[:d])[0] / lam0 - x[d]
        jac = lambda x: np.hstack([self.eig(x[:d])[1] / lam0, -np.ones((MB, 1))])

        def cb(xk):
            trace.append((self.count, time.perf_counter() - t0, self.eig(xk[:d])[0][0]))
        obj, gobj = (lambda x: -x[d]), (lambda x: np.r_[np.zeros(d), -1.0])
        try:
            if method == "SLSQP":
                r = so.minimize(obj, x0, jac=gobj, method="SLSQP", bounds=[(LO, HI)] * d + [(BETA_LO, BETA_HI)],
                                constraints=[{"type": "ineq", "fun": fun, "jac": jac},
                                             {"type": "eq", "fun": lambda x: np.array([self.areas @ x[:d] - self.V]),
                                              "jac": lambda x: np.r_[self.areas, 0.0][None, :]}],
                                callback=cb, options={"maxiter": 10_000, "ftol": 1e-12})
            else:
                r = so.minimize(obj, x0, jac=gobj, method="trust-constr",
                                bounds=so.Bounds([LO] * d + [BETA_LO], [HI] * d + [BETA_HI]),
                                constraints=[so.NonlinearConstraint(fun, 0, np.inf, jac=jac),
                                             so.LinearConstraint(np.r_[self.areas, 0.0][None, :], self.V, self.V)],
                                callback=lambda xk, st: cb(xk),
                                options={"maxiter": 10_000, "gtol": 1e-10, "xtol": 1e-12})
            return None, str(r.message) if r.success else "failed: " + str(r.message)
        except _Budget:
            return None, "budget"

    def mma(self, h0, trace):
        import pymoto as pym
        t0, d, bm = time.perf_counter(), self.d, self
        lam0 = self.eig(h0)[0][0]

        class Responses(pym.Module):
            def __call__(self, x):
                self.x = x.copy()
                w, self.G = bm.eig(x[:d])
                trace.append((bm.count, time.perf_counter() - t0, w[0]))
                return (np.array([-x[d]]), *[np.array([x[d] - w[i] / lam0]) for i in range(MB)],
                        np.array([(bm.areas @ x[:d] - bm.V) / bm.V]))

            def _sensitivity(self, *dys):
                g = np.zeros_like(self.x)
                if dys[0] is not None:
                    g[d] -= dys[0][0]
                for i in range(MB):
                    if dys[1 + i] is not None:
                        g[d] += dys[1 + i][0]
                        g[:d] -= dys[1 + i][0] * self.G[i] / lam0
                if dys[-1] is not None:
                    g[:d] += dys[-1][0] * bm.areas / bm.V
                return g

        sx = pym.Signal("x", state=np.r_[h0, 1.0])
        with pym.Network() as fn:
            outs = Responses()(sx)
        try:
            pym.minimize_mma(sx, list(outs), function=fn, maxit=10_000, xmin=np.r_[[LO] * d, BETA_LO],
                             xmax=np.r_[[HI] * d, BETA_HI], verbosity=0, tolx=1e-9, tolf=1e-12)
            return None, "converged"
        except _Budget:
            return None, "budget"

    def naive(self, h0, trace, alpha=0.02):
        t0, h = time.perf_counter(), h0.copy()
        try:
            while True:
                w, G = self.eig(h)
                trace.append((self.count, time.perf_counter() - t0, w[0]))
                g = G[0] - self.areas * (self.areas @ G[0]) / (self.areas @ self.areas)
                h = self.project(h + alpha * g / np.linalg.norm(g))
        except _Budget:
            return None, "budget"

    def run(self, name, h0):
        fn = {"ours": self.ours, "SLSQP bound": lambda h, t: self.bound(h, t, "SLSQP"),
              "trust-constr bound": lambda h, t: self.bound(h, t, "trust-constr"),
              "MMA bound (pyMOTO)": self.mma, "naive fixed-step": self.naive}[name]
        self.count, self._cache, trace = 0, collections.OrderedDict(), []
        self._mma = "MMA" in name
        self._reset_best()
        t0 = time.perf_counter()
        with self.counting(), warnings.catch_warnings():
            warnings.simplefilter("ignore")
            h, stop = fn(h0.copy(), trace)
        sec, n_eig = time.perf_counter() - t0, self.count
        if h is None:                                  # best feasible evaluated design, as ours reports
            h = self.best_h if self.best_h is not None else h0
            reach = self.bestlog
        else:
            reach = trace
        short = stop if name == "ours" else ("budget" if stop == "budget" else
                                             "failed" if stop.startswith("failed") else "converged")
        w, _ = smallest_eigenpairs(self.pen.K(h), self.pen.M(h), 5)     # final re-solve, not counted
        vals = np.array([x[2] for x in trace])
        return dict(method=name, final=w[0] / P4, mult=int(np.sum(w - w[0] <= 1e-4 * w[0])),
                    decreases=int(np.sum(np.diff(vals) < -1e-12 * np.abs(vals).max())), eig=n_eig, sec=sec,
                    vol=abs(self.areas @ h - self.V) / self.V,
                    bnd=max(0.0, float(np.max(LO - h)), float(np.max(h - HI))), stop=short, stop_full=stop, trace=reach)


class _Budget(Exception):
    pass


def starts(n_seeds: int, bm: Benchmark):
    """Seed 0: the uniform (D4-symmetric) design. Seeds 1..: projections of 1 + U(-0.25, 0.25)."""
    rng = np.random.default_rng(0)
    return [np.ones(bm.d)] + [bm.project(1 + 0.25 * rng.uniform(-1, 1, bm.d)) for _ in range(n_seeds - 1)]


def _hit(row, target):
    k = next((x for x in row["trace"] if x[2] / P4 >= target), None)
    return None if k is None else (k[0], k[1])


def main(quick: bool = False, budget: int | None = None, n_seeds: int | None = None, markdown: str | None = None):
    t_start = time.perf_counter()
    budget = budget or (100 if quick else 200)
    n_seeds = n_seeds or (2 if quick else 5)
    bm = Benchmark(budget)
    try:
        import pymoto  # noqa: F401
        methods = METHODS
    except ImportError:
        methods = tuple(m for m in METHODS if "pyMOTO" not in m)
        print("MMA bound (pyMOTO): not run (pymoto not installed)")
    head = (f"plate + post, {bm.pen.K(np.ones(bm.d)).shape[0]} dofs, {bm.d} zones (unequal areas, volume weighted by a_z); "
            f"budget {budget} eigensolves; {n_seeds} starts")
    print(head)
    lines = ["| start | method | final λ₀/π⁴ | mult | eigensolves to 99 % | to 99.9 % | decreases | "
             "volume residual | stop |", "|---|---|---|---|---|---|---|---|---|"]
    rows = []
    for s, h0 in enumerate(starts(n_seeds, bm)):
        res = [bm.run(m, h0) for m in methods]
        best = max(r["final"] for r in res)
        lam_start = smallest_eigenpairs(bm.pen.K(h0), bm.pen.M(h0), 1)[0][0] / P4
        print(f"\nstart {s}: lambda_0/pi^4 = {lam_start:.5f}; best found {best:.6f}")
        print("  method               final     mult eigs  time   99%:eigs/s   99.9%:eigs/s  decr  vol      stop")
        for r in res:
            r["seed"], r["best"] = s, best
            h99, h999 = _hit(r, 0.99 * best), _hit(r, 0.999 * best)
            r["to99"], r["to999"] = h99, h999
            f = lambda h: "-" if h is None else f"{h[0]}/{h[1]:.1f}s"
            print(f"  {r['method']:20s} {r['final']:.6f} {r['mult']:4d} {r['eig']:4d} {r['sec']:5.1f}s "
                  f"{f(h99):>12s} {f(h999):>13s} {r['decreases']:5d}  {r['vol']:.1e}  {r['stop']}")
            lines.append(f"| {s} | {r['method']} | {r['final']:.4f} | {r['mult']} | "
                         f"{'—' if h99 is None else h99[0]} | {'—' if h999 is None else h999[0]} | "
                         f"{r['decreases']} | {r['vol']:.0e} | {r['stop']} |")
            rows.append(r)
    summ = ["summary over starts (reached = within 0.1 % of the best value found from that start):"]
    for m in methods:
        rs = [r for r in rows if r["method"] == m]
        reach = [r for r in rs if r["to999"] is not None]
        med = np.median([r["to999"][0] for r in reach]) if reach else float("nan")
        summ.append(f"  {m:20s} reached {len(reach)}/{len(rs)}; median eigensolves to reach {med:6.1f}; "
                    f"decreases {sum(r['decreases'] for r in rs):4d}; max volume residual "
                    f"{max(r['vol'] for r in rs):.1e}; max bound violation {max(r['bnd'] for r in rs):.1e}")
    print("\n" + "\n".join(summ))
    if markdown:
        import importlib.metadata as md
        ver = {k: md.version(k) for k in ("scipy", "numpy", "scikit-fem")}
        try:
            ver["pymoto"] = md.version("pymoto")
        except md.PackageNotFoundError:
            ver["pymoto"] = "not installed"
        msgs = sorted({f"{r['method']}: {r['stop_full']}" for r in rows if r["method"] != "ours"
                       and r["stop"] != "budget"})
        out = ["# Fundamental-frequency maximisation benchmark", "",
               "Label: NUMERICAL (discrete problem). Reported as measured, whoever wins.", "",
               f"Plate + central post, 512 dofs, 25 thickness zones (unequal areas, 0.03125 to 0.046875; volume "
               f"weighted by a_z). Budget {budget} eigensolves per run, {n_seeds} starts (start 0 uniform, the rest "
               f"random feasible). Bound formulation: beta / lambda_0(start) in [0.5, 2] for SLSQP, trust-constr and "
               f"MMA. Versions: scipy {ver['scipy']}, numpy {ver['numpy']}, scikit-fem {ver['scikit-fem']}, "
               f"pymoto {ver['pymoto']}.", "",
               "Every method reports its best feasible evaluated design within the budget; 'reached' counts "
               "eigensolves until the running best is within 0.1 % (99.9 %) or 1 % (99 %) of the best final value "
               "over methods from that start. 'decreases' counts decreases between successive accepted iterates "
               "(ours: accepted steps; SLSQP and trust-constr: callback iterates, where a merit-function method "
               "may decrease lambda_0 legitimately; MMA and naive: every iterate, all accepted).", ""]
        out += lines + ["", "```"] + summ + ["```", ""]
        out += ["Stop codes: `budget` = eigensolve budget exhausted; `converged` / `failed` = solver's own "
                "termination. Full solver messages:"] + [f"- {m}" for m in msgs] + [""]
        with open(markdown, "w", encoding="utf-8") as fh:
            fh.write("\n".join(out))
    print(f"label: NUMERICAL (discrete problem); {time.perf_counter() - t_start:.1f} s")
    # self-checks (spec 10.1): constraints hold and our accepted iterates never decrease; nothing about who wins
    for r in rows:
        assert r["bnd"] <= 1e-12
        assert r["vol"] <= (1e-3 if "MMA" in r["method"] else 1e-9)
        if r["method"] == "ours":
            assert r["decreases"] == 0
    return rows


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--markdown", default=None)
    a = ap.parse_args()
    main(quick=a.quick, markdown=a.markdown)
