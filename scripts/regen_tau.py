"""B12 driver: moments and bounds for tau(P_-(h_{R,z0})).  Writes src/bolzakit/data/tau.json (--quick: scripts/outputs/tau_results_quick.json).

Usage:  python3 -B run_tau.py [--quick] [--rebound]
  --rebound: recompute bounds from the moments already stored in outputs/tau_results.json
  --quick: smaller J (used by the test-suite smoke test); default reproduces the LOG.md numbers.
"""
import json
import math
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bolzakit.exact import group as G  # noqa: E402
from bolzakit.spectral import tau as T    # noqa: E402

QUICK = "--quick" in sys.argv
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = (os.path.join(_ROOT, "scripts", "outputs", "tau_results_quick.json") if QUICK
       else os.path.join(_ROOT, "src", "bolzakit", "data", "tau.json"))


def bounds_from(m, M, rel_err, label):
    mu, err = T.cheb_moments(m, M, rel_err)
    res = {}
    for J in range(2, (len(m) - 1) // 2 + 1):
        b = T.rigorous_bounds(mu[:2 * J + 1], err[:2 * J + 1])
        res[f"moments_to_{2*J}"] = dict(lower_mu_lt0=float(b["lower_lt0"]), upper_mu_le0=float(b["upper_le0"]),
                                       width=float(b["upper_le0"] - b["lower_lt0"]),
                                       cert_shift=float(max(b["eps_upper"], b["eps_lower"])),
                                       kpm_jackson_estimate=float(T.kpm_estimate(mu[:2 * J + 1])))
        print(label, 2 * J, res[f"moments_to_{2*J}"], flush=True)
    # every row is a valid rigorous bound, so the intersection is too
    res["best_rigorous"] = dict(lower_mu_lt0=max(r["lower_mu_lt0"] for r in res.values()),
                                upper_mu_le0=min(r["upper_mu_le0"] for r in res.values()))
    return res, mu


def case_systole(J):
    X, c, d = T.systole_shell()
    t = time.time()
    m, sizes, times = T.moments(X, np.ones(8, dtype=np.int64), J, integer=True, verbose=True)
    cval = float(c[0])
    # spectral measure of h = c*A has the same mass on (-inf,0) as that of A; work with A, M = 8 = ||A||_l1.
    res, mu = bounds_from(m, 8, [0.0] * len(m), "systole")
    chris = {}
    for n in range(2, J + 2):
        lam = T.christoffel_bound(m[:2 * n - 1])
        chris[f"n={n} (moments to {2*n-2})"] = dict(exact=str(lam), float=float(lam))
    lam = T.christoffel_bound(m)
    sym_lower = max((1 - float(lam)) / 2, res["best_rigorous"]["lower_mu_lt0"])
    return dict(
        h="c * sum_{k=0}^{7} u_{g_k},  c = G1(l_sys) = -(1/2pi) ln tanh(arccosh(1+sqrt2)) = %.15g; z0 = 0, R = l_sys" % cval,
        support_size=8, l1_norm_bound=8 * cval,
        exact_closed_walk_counts_W_k=[str(x) for x in m],
        moments_note="m_k(h) = c^k W_k; W_k exact integers (closed walks of length k at e in the Cayley graph)",
        support_sizes=sizes, step_times_s=times, total_time_s=time.time() - t,
        christoffel_mu_atom_at_0_upper=chris,
        lp_certificate_bounds=res,
        symmetry="supp h is in the odd part of the parity character (relator has even length) => mu_h symmetric, "
                 "so tau(P_-) = (1 - mu{0})/2 exactly",
        tau_P_minus_rigorous_interval=[sym_lower, 0.5],
        tau_P_minus_if_no_atom=0.5,
    )


def case_general(R, z0, J):
    X, c, d = T.build_h(R, z0)
    t = time.time()
    m, sizes, times = T.moments(X, c, J, verbose=True)
    rel = T.moment_rel_errors(m, len(c), sizes)
    M = float(c.sum()) * (1 + 1e-12)
    res, mu = bounds_from(m, M, rel, f"R={R} z0={z0}")
    best = res["best_rigorous"]
    shells = {}
    for dd in np.round(d, 6):
        shells[str(dd)] = shells.get(str(dd), 0) + 1
    return dict(
        h=f"sum over gamma != 1 with d(z0,gamma z0) <= {R} of G1(d(z0,gamma z0)) u_gamma, z0 = {z0}",
        support_size=len(c), shells_distance_count=shells, l1_norm_bound=M,
        moments=[float(x) for x in m], moment_rel_err_bounds=rel, chebyshev_moments=[float(x) for x in mu],
        support_sizes=sizes, step_times_s=times, total_time_s=time.time() - t,
        lp_certificate_bounds=res,
        tau_P_minus_rigorous_interval=[best["lower_mu_lt0"], best["upper_mu_le0"]],
        note="interval is [lower bound on mu((-inf,0)), upper bound on mu((-inf,0])]; tau(P_-) = mu((-inf,0)) lies in it",
    )


def rebound(path):
    """Recompute all bounds from the moments stored in `path` (no group computation)."""
    with open(path) as f:
        out = json.load(f)
    sy = out["systole_shell"]
    m = [int(x) for x in sy["exact_closed_walk_counts_W_k"]]
    res, _ = bounds_from(m, 8, [0.0] * len(m), "systole")
    sy["lp_certificate_bounds"] = res
    lam = T.christoffel_bound(m)
    sy["tau_P_minus_rigorous_interval"] = [max((1 - float(lam)) / 2, res["best_rigorous"]["lower_mu_lt0"]), 0.5]
    for r in out["general"]:
        res, _ = bounds_from(r["moments"], r["l1_norm_bound"], r["moment_rel_err_bounds"], f"R={r['R']} z0={r['z0']}")
        r["lp_certificate_bounds"] = res
        r["tau_P_minus_rigorous_interval"] = [res["best_rigorous"]["lower_mu_lt0"], res["best_rigorous"]["upper_mu_le0"]]
    out["rebounded"] = time.strftime("%Y-%m-%d %H:%M:%S")
    with open(path, "w") as f:
        json.dump(out, f, indent=1, default=str)
    print("rebounded", path)


if __name__ == "__main__":
    if "--rebound" in sys.argv:
        rebound(OUT)
        sys.exit(0)
    out = dict(generated=time.strftime("%Y-%m-%d %H:%M:%S"), quick=QUICK,
               G1="G1(d) = -(1/(2 pi)) ln tanh(d/2)")
    out["systole_shell"] = case_systole(5 if QUICK else 8)
    cases = [(4.3, 0j, 3 if QUICK else 6), (5.0, 0j, 3 if QUICK else 5)]
    if not QUICK:
        cases += [(4.3, 0.15 + 0.07j, 5), (5.0, 0.15 + 0.07j, 4)]
    out["general"] = []
    for R, z0, J in cases:
        r = case_general(R, z0, J)
        r["R"], r["z0"], r["J"] = R, [z0.real, z0.imag], J
        out["general"].append(r)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(out, f, indent=1, default=str)
    print("wrote", OUT)
