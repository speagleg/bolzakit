"""trim_degenerate_pair: first-order trim LP plus Newton re-solves."""
import numpy as np
import pytest
import scipy.linalg as sla

from bolzakit.sensitivity import ClusterSplitError
from bolzakit.sensitivity.trim import TrimResult, trim_degenerate_pair

E = np.eye(4)


def _outer(v):
    v = np.asarray(v, float)
    return np.outer(v, v)


def test_diagonal_split_is_cancelled_exactly_by_the_cheapest_site():
    K = np.diag([1.0, 2.0, 2.01, 5.0])
    dKs = [_outer(E[1]), _outer(E[2])]           # stiffening dof 2 would widen the split; dof 1 closes it
    r = trim_degenerate_pair(K, np.eye(4), dKs, None, (1, 2))
    assert isinstance(r, TrimResult)
    assert r.status == "converged" and r.n_iter == 1 and r.n_eigensolves == 2
    assert np.allclose(r.t, [0.01, 0.0], atol=1e-14)
    assert r.rel_split_after <= 1e-12 and abs(r.rel_split_before - 0.01 / 2.005) < 1e-12
    assert r.label == "NUMERICAL" and "discrete" in r.caveat


def test_coupled_pair_needs_a_coupling_site_else_infeasible():
    eps = 1e-3
    K = np.diag([1.0, 2.0, 2.0, 5.0])
    K[1, 2] = K[2, 1] = eps                       # pair 2 -+ eps, eigenvectors (e1 -+ e2)/sqrt 2
    diag_only = [_outer(E[1]), _outer(E[2])]      # cannot touch the off-diagonal: split >= 2 eps for any t
    r = trim_degenerate_pair(K, np.eye(4), diag_only, None, (1, 2))
    assert r.status == "infeasible" and np.all(r.t == 0) and r.n_eigensolves == 1
    sites = diag_only + [_outer(E[1] - E[2]) / 2]  # stiffens the lower mode (e1 - e2)/sqrt 2 only
    r = trim_degenerate_pair(K, np.eye(4), sites, None, (1, 2))
    assert r.status == "converged" and np.allclose(r.t, [0, 0, 2 * eps], atol=1e-13)


def test_upper_bound_makes_it_infeasible():
    K = np.diag([1.0, 2.0, 2.01, 5.0])
    r = trim_degenerate_pair(K, np.eye(4), [_outer(E[1])], None, (1, 2), upper=0.005)
    assert r.status == "infeasible" and r.t[0] == 0


def test_cost_weights_choose_the_site():
    K = np.diag([1.0, 2.0, 2.01, 5.0])
    dKs = [_outer(E[1]), 2 * _outer(E[1])]        # the second does the same job with half the amount
    r1 = trim_degenerate_pair(K, np.eye(4), dKs, None, (1, 2))
    r2 = trim_degenerate_pair(K, np.eye(4), dKs, None, (1, 2), cost=[1.0, 3.0])
    assert np.allclose(r1.t, [0, 0.005]) and np.allclose(r2.t, [0.01, 0])


def _random_mass_problem(seed):
    rng = np.random.default_rng(seed)
    n = 10
    K = np.diag([1.0, 2.0, 2.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0])
    X = rng.standard_normal((n, n))
    M = np.eye(n) + 2e-3 * (X + X.T)              # imperfect mass: the pair (1, 2) splits by O(1e-3)
    dMs = [1e-2 * _outer(rng.standard_normal(n)) for _ in range(8)]   # point-mass-like rank-one trims
    return K, M, dMs


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_newton_converges_quadratically_and_matches_independent_eigh(seed):
    K, M, dMs = _random_mass_problem(seed)
    r = trim_degenerate_pair(K, M, None, dMs, (1, 2))
    assert r.status == "converged" and r.n_iter <= 6
    assert np.all(r.t >= 0) and np.count_nonzero(r.t) <= 2  # an LP vertex: at most two trims
    Mt = M + sum(tl * d for tl, d in zip(r.t, dMs))
    w = sla.eigh(K, Mt, eigvals_only=True)
    assert (w[2] - w[1]) / w[1:3].mean() <= 1e-9
    rel = [h.rel_split for h in r.history]
    assert rel[0] > 1e-4 and rel[-1] <= 1e-10
    steps = [h.step for h in r.history]
    assert steps[0] == "lp" and steps[-1] == "" and set(steps[1:-1]) <= {"newton"}


def test_negative_lower_bound_costs_removal_like_addition():
    K = np.diag([1.0, 2.0, 2.01, 5.0])
    ablate = trim_degenerate_pair(K, np.eye(4), [_outer(E[2])], None, (1, 2), lower=-1.0)
    add = trim_degenerate_pair(K, np.eye(4), [_outer(E[1])], None, (1, 2), lower=-1.0)
    assert ablate.status == "converged" and np.allclose(ablate.t, [-0.01], atol=1e-13)   # not maximal ablation
    assert add.status == "converged" and np.allclose(add.t, [0.01], atol=1e-13)         # the mirrored-sign trim
    assert abs(ablate.cost - add.cost) < 1e-13 and abs(ablate.cost - 0.01) < 1e-13
    both = trim_degenerate_pair(K, np.eye(4), [_outer(E[1]), _outer(E[2])], None, (1, 2), lower=-1.0)
    assert both.status == "converged" and abs(np.abs(both.t).sum() - 0.01) < 1e-13 and both.cost == pytest.approx(0.01)
    assert np.count_nonzero(both.t) == 1


def test_resolve_lp_when_newton_leaves_the_bounds_dk_and_dm_together():
    rng = np.random.default_rng(29)               # found by search: lp, lp (Newton left the bounds), newton
    n = 10
    K = np.diag([1.0, 2.0, 2.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0])
    X = rng.standard_normal((n, n))
    M = np.eye(n) + 2e-3 * (X + X.T)
    dMs = [1e-2 * _outer(rng.standard_normal(n)) for _ in range(8)]
    dKs = [0.05 * _outer(rng.standard_normal(n)) for _ in range(8)]
    r = trim_degenerate_pair(K, M, dKs, dMs, (1, 2), lower=-0.5, upper=0.3)
    assert [h.step for h in r.history] == ["lp", "lp", "newton", ""] and r.status == "converged"
    assert np.all(r.t >= -0.5) and np.all(r.t <= 0.3) and np.any(r.t < 0)
    Kt = K + sum(tl * d for tl, d in zip(r.t, dKs))
    Mt = M + sum(tl * d for tl, d in zip(r.t, dMs))
    w = sla.eigh(Kt, Mt, eigvals_only=True)
    assert (w[2] - w[1]) / w[1:3].mean() <= 1e-9


def test_status_max_iter_and_stalled_at_the_resolution_floor():
    K, M, dMs = _random_mass_problem(0)
    r = trim_degenerate_pair(K, M, None, dMs, (1, 2), max_iter=1)
    assert r.status == "max_iter" and r.n_iter == 1 and r.rel_split_after > 10 * r.rel_floor
    rng = np.random.default_rng(0)                # lambda_max / lambda_pair = 1e6: floor ~ 2e-10 > tol below
    Q, _ = np.linalg.qr(rng.standard_normal((6, 6)))
    K = Q @ np.diag([1e-3, 1.0, 1.0 + 2e-3, 3, 50, 1e6]) @ Q.T
    K = (K + K.T) / 2
    sites = [_outer(Q[:, 1]), _outer(Q[:, 2]), _outer(Q[:, 1] + Q[:, 2]) / 2]
    r = trim_degenerate_pair(K, np.eye(6), sites, None, (1, 2), tol=1e-14)
    assert r.status == "stalled" and r.rel_split_after <= 2 * r.rel_floor and 1e-11 < r.rel_floor < 1e-9


def test_max_iter_just_above_tol_is_not_labelled_stalled():
    # N3: floor < tol < split < 2 tol at the end of the budget is budget exhaustion, not the resolution floor
    K, M, dMs = _random_mass_problem(0)
    r1 = trim_degenerate_pair(K, M, None, dMs, (1, 2), max_iter=1)
    tol = 0.7 * r1.rel_split_after                # same trajectory (tol only decides convergence)
    r = trim_degenerate_pair(K, M, None, dMs, (1, 2), max_iter=1, tol=tol)
    assert r.rel_floor < tol < r.rel_split_after < 2 * tol
    assert r.status == "max_iter" and r.rel_split_after == r1.rel_split_after


def test_infeasible_after_progress_keeps_the_best_iterate():
    # Mass at dof 2 lowers lambda_2 = 2.01 / (1 + t); the first-order LP step t1 = 0.01 / 2.005 undershoots the exact
    # 0.005. The Newton correction then needs t > upper, and the re-solved LP is infeasible within the bounds.
    K = np.diag([1.0, 2.0, 2.01, 5.0])
    t1 = 0.01 / 2.005
    r = trim_degenerate_pair(K, np.eye(4), None, [_outer(E[2])], (1, 2), upper=0.004995)
    assert t1 < 0.004995 < 0.005
    assert r.status == "infeasible" and [h.step for h in r.history] == ["lp", ""] and r.n_eigensolves == 2
    assert abs(r.t[0] - t1) < 1e-12 and r.rel_split_after < 1e-2 * r.rel_split_before


def test_highs_failure_is_lp_failed_not_infeasible(monkeypatch):
    import types
    import bolzakit.sensitivity.trim as trim_mod
    K = np.diag([1.0, 2.0, 2.01, 5.0])
    for code in (1, 4):                           # iteration limit, numerical difficulties
        monkeypatch.setattr(trim_mod, "linprog", lambda *a, code=code, **k: types.SimpleNamespace(status=code, x=None))
        r = trim_degenerate_pair(K, np.eye(4), [_outer(E[1])], None, (1, 2))
        assert r.status == "lp_failed" and r.t[0] == 0 and r.n_eigensolves == 1
    monkeypatch.setattr(trim_mod, "linprog", lambda *a, **k: types.SimpleNamespace(status=2, x=None))
    assert trim_degenerate_pair(K, np.eye(4), [_outer(E[1])], None, (1, 2)).status == "infeasible"


def test_lambda_max_fallback_never_underestimates(monkeypatch):
    # M-2: when ARPACK fails the floor must not come from a lower bound on lambda_max (that would shrink rel_floor)
    import scipy.sparse as sp
    import bolzakit.sensitivity.trim as trim_mod

    def boom(*a, **k):
        raise RuntimeError("ARPACK did not converge")
    monkeypatch.setattr(trim_mod.spla, "eigsh", boom)
    for n in (500, 5000):                         # dense fallback; Gershgorin upper bound
        d = np.linspace(1.0, 1e4, n)
        K = sp.diags(d).tocsc()
        M = sp.diags(np.full(n, 2.0)).tocsc() + sp.eye(n, k=1, format="csc") * 0.1 + sp.eye(n, k=-1, format="csc") * 0.1
        lam = trim_mod._lambda_max(K, M)
        assert lam >= 1e4 / 2.2 * (1 - 1e-12)     # lambda_max(K, M) >= max(d) / lambda_max(M) >= 1e4 / 2.2
        if n == 500:
            ref = sla.eigh(K.toarray(), M.toarray(), eigvals_only=True)[-1]
            assert abs(lam - ref) <= 1e-10 * ref
        else:
            assert lam >= 1e4 / 2.0               # an upper bound: ||K||_inf / gamma(M) = 1e4 / 1.8


def test_argument_validation():
    K = np.diag([1.0, 2.0, 2.01, 5.0])
    with pytest.raises(ValueError, match="pair"):
        trim_degenerate_pair(K, np.eye(4), [E], None, (1, 3))
    with pytest.raises(ValueError, match="both be None"):
        trim_degenerate_pair(K, np.eye(4), None, None, (1, 2))
    with pytest.raises(ValueError, match="same length"):
        trim_degenerate_pair(K, np.eye(4), [E, E], [E], (1, 2))
    with pytest.raises(ValueError, match="lower"):
        trim_degenerate_pair(K, np.eye(4), [E], None, (1, 2), lower=1.0, upper=0.0)
    with pytest.raises(ValueError, match="at least one"):
        trim_degenerate_pair(K, np.eye(4), [], None, (1, 2))
    for kw in ({"lower": [0.0, 0.0]}, {"upper": [1.0, 1.0, 1.0]}, {"cost": [1.0, 1.0]}):
        with pytest.raises(ValueError, match="length 1"):
            trim_degenerate_pair(K, np.eye(4), [E], None, (1, 2), **kw)
    with pytest.raises(ValueError, match="non-negative"):
        trim_degenerate_pair(K, np.eye(4), [E], None, (1, 2), cost=[-1.0])


def test_pair_must_be_isolated():
    K = np.diag([1.0, 2.0, 2.01, 2.02])           # lambda_3 is as close as the split itself
    with pytest.raises(ClusterSplitError):
        trim_degenerate_pair(K, np.eye(4), [_outer(E[1])], None, (1, 2))
    K = np.diag([1.0, 2.0, 2.0005, 2.0006])       # a third eigenvalue within gap_ratio_min |lam_bar|: the core message
    with pytest.raises(ClusterSplitError, match="splits"):
        trim_degenerate_pair(K, np.eye(4), [_outer(E[1])], None, (1, 2))
    K = np.diag([1.0, 2.0, 2.01, 5.0])            # M-8: (2, 3) is not a near-degenerate pair at all
    with pytest.raises(ClusterSplitError, match="not nearly degenerate"):
        trim_degenerate_pair(K, np.eye(4), [_outer(E[1])], None, (2, 3))


def test_peak_objective_spreads_the_trim_over_equivalent_sites():
    K = np.diag([1.0, 2.0, 2.03, 5.0])
    dKs = [_outer(E[1])] * 3                      # three identical candidates (e.g. electrodes of one group)
    s = trim_degenerate_pair(K, np.eye(4), dKs, None, (1, 2))
    p = trim_degenerate_pair(K, np.eye(4), dKs, None, (1, 2), objective="peak")
    assert s.status == p.status == "converged"
    assert np.count_nonzero(s.t) == 1 and s.t.sum() == pytest.approx(0.03)
    assert np.allclose(p.t, 0.01, atol=1e-9)
    with pytest.raises(ValueError, match="objective"):
        trim_degenerate_pair(K, np.eye(4), dKs, None, (1, 2), objective="max")


def test_backtracking_rescues_a_strongly_nonlinear_trim():
    rng = np.random.default_rng(120)                            # found by search: without halving the run stalls
    K = np.diag([1.0, 2.0, 2.02, 2.6, 5.0, 9.0])
    dKs, dMs = [], []
    for _ in range(6):
        v = rng.standard_normal(6)
        v[3] *= 3                                               # strong coupling to the mode at 2.6
        dKs.append(0.3 * _outer(v))
        dMs.append(0.3 * 0.47 * _outer(v))                      # first-order pair effect nearly cancels
    kw = dict(lower=-1.0, upper=1.0, max_iter=12)
    r0 = trim_degenerate_pair(K, np.eye(6), dKs, dMs, (1, 2), backtrack=0, **kw)
    r4 = trim_degenerate_pair(K, np.eye(6), dKs, dMs, (1, 2), **kw)
    assert r0.status == "stalled" and r0.rel_split_after > 1e-4
    assert r4.status == "converged" and "halve" in [h.step for h in r4.history]


# ------------------------------------------------------------------ FE resonators (examples/gyro_trim.py; skfem)
def _gyro():
    pytest.importorskip("skfem")
    import pathlib, sys
    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2] / "examples"))
    import gyro_trim
    return gyro_trim


def test_ring_rate_is_three_tenths_and_trim_gains_three_orders():
    a = _gyro().part_a()
    assert abs(a["rate"] - 0.3) < 5e-4                       # 3/10 (w/R = 0.05, N_theta = 96: 0.29984)
    r = a["res"]
    assert r.rel_split_after <= 1e-3 * r.rel_split_before
    assert a["dense"] <= 1e-3 * r.rel_split_before            # independent dense eigensolve
    assert np.count_nonzero(r.t) == 2 and np.all(r.t >= 0)


def test_spoked_ring_trim_converges_and_beats_antinode_mass():
    b = _gyro().part_b(seeds=[1])
    x = b["runs"][0]
    r = x["res"]
    assert r.status == "converged" and r.rel_split_after <= 1e-10
    assert x["dense"] <= 1e-9                                 # independent dense eigensolve
    for k in ("naive", "damped"):                             # single-site antinode baselines, 40 rounds
        tied = x[f"{k}_rel"][-1] <= x["stop"]
        assert (not tied) or r.cost < x[f"{k}_t"].sum()       # a baseline either fails to tie or costs more mass
    assert x["damped_rel"][-1] > x["stop"]                    # seed 1: the damped baseline has not tied the pair
