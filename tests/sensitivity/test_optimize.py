import time
import numpy as np, pytest
from bolzakit.sensitivity.optimize import ascent_direction


def _lmin(Au, d):
    return np.linalg.eigvalsh(np.einsum("l,lab->ab", d, Au))[0]


# ------------------------------------------------------------------ ascent_direction, analytic
def test_m1_is_normalised_projected_gradient():
    g = np.array([3.0, -4.0, 0.0])
    r = ascent_direction(g.reshape(3, 1, 1))
    assert np.allclose(r.d, g / 5) and abs(r.phi - 5) < 1e-12 and r.duality_gap <= 1e-10
    N = np.array([[1.0], [0.0], [0.0]])
    r = ascent_direction(g.reshape(3, 1, 1), null_basis=N)
    assert np.allclose(r.d, [1, 0, 0]) and abs(r.phi - 3) < 1e-12


def test_split_case_phi_one():
    A = np.array([np.diag([1.0, -1.0]), np.diag([1.0, 1.0])])
    r = ascent_direction(A)
    assert np.allclose(r.d, [0, 1], atol=1e-6) and abs(r.phi - 1) < 1e-9 and r.duality_gap <= 1e-10
    assert np.allclose(r.rates, [1, 1], atol=1e-6)


def test_stationary_case_Y_not_unique():
    A = np.array([np.diag([1.0, -1.0]), np.diag([-1.0, 1.0])])
    r = ascent_direction(A, stat_tol=1e-12, f_value=1.0)
    Y = r.Y_upper
    assert np.all(np.linalg.eigvalsh(Y) >= -1e-12) and abs(np.trace(Y) - 1) < 1e-12
    g = np.einsum("lab,ab->l", A, Y)
    assert np.linalg.norm(g) <= 1e-10 and abs(Y[0, 0] - 0.5) < 1e-10 and abs(Y[1, 1] - 0.5) < 1e-10
    assert r.phi_upper <= 1e-10 and np.all(r.d == 0)


def test_gap_closed_form_stationary():
    r = ascent_direction(np.array([np.diag([1.0, 0.0])]), np.array([np.diag([-1.0, 0.0])]), stat_tol=1e-12)
    assert r.phi_upper <= 1e-10 and r.Y_lower is not None and r.rates_lower is not None


def test_bounds_cone_R1():
    A = np.array([[[1.0]]])
    r = ascent_direction(A, bound_signs=np.array([1.0]))
    assert abs(r.phi - 1) < 1e-12 and np.allclose(r.d, [1.0]) and r.duality_gap <= 1e-10
    r = ascent_direction(A, bound_signs=np.array([-1.0]), stat_tol=1e-12)
    assert r.phi_upper <= 1e-10 and abs(r.mu[0] - 1) < 1e-8 and np.all(r.d == 0)


def test_rank_deficient_stationary_stops_by_rule_b_R4():
    rng = np.random.default_rng(0)
    d, m = 64, 4
    A = rng.standard_normal((d, m, m)); A = (A + A.transpose(0, 2, 1)) / 2
    A[:, 0, 0] = 0.0                                   # Y* = e0 e0^T (rank 1) gives g(Y*) = 0
    t0 = time.perf_counter()
    r = ascent_direction(A, stat_tol=1e-8, f_value=1.0)
    elapsed = time.perf_counter() - t0
    assert r.stop == "stationary" and r.n_iter < 2000 and r.phi_upper <= 1e-8
    assert abs(r.Y_upper[0, 0] - 1) < 1e-6
    print(f"rank-1 stationary: {r.n_iter} iterations, {elapsed * 1e3:.1f} ms")


def test_generic_d64_m4_fast():
    rng = np.random.default_rng(1)
    A = rng.standard_normal((64, 4, 4)); A = (A + A.transpose(0, 2, 1)) / 2
    t0 = time.perf_counter()
    r = ascent_direction(A)
    elapsed = time.perf_counter() - t0
    assert r.stop == "gap" and r.phi > 0
    print(f"generic d=64 m=4: {r.n_iter} iterations, {elapsed * 1e3:.1f} ms")


@pytest.mark.parametrize("m", [2, 3])
def test_optimality_of_dstar(m):
    rng = np.random.default_rng(10 + m)
    A = rng.standard_normal((6, m, m)); A = (A + A.transpose(0, 2, 1)) / 2
    N = np.linalg.qr(rng.standard_normal((6, 4)))[0]       # D = a random 4-dim subspace
    r = ascent_direction(A, null_basis=N)
    G = np.linalg.norm(A.reshape(6, -1), 2)
    best = r.phi
    X = rng.standard_normal((10_000, 4)); X /= np.linalg.norm(X, axis=1, keepdims=True)
    D = X @ N.T
    vals = np.linalg.eigvalsh(np.einsum("nl,lab->nab", D, A))[:, 0]
    assert np.all(vals <= best + r.duality_gap + 1e-12 * max(1, G))
    assert abs(_lmin(A, r.d) - r.phi) < 1e-12 and abs(np.linalg.norm(r.d) - 1) < 1e-12


# ------------------------------------------------------------------ fix round 1: stationarity, restarts, census
def _sym(rng, *shape):
    A = rng.standard_normal(shape)
    return (A + A.transpose(0, 2, 1)) / 2


def _gscale(*stacks):
    return max(1.0, np.linalg.norm(np.hstack([S.reshape(S.shape[0], -1) for S in stacks]), 2))


@pytest.mark.parametrize("case", ["free_d3_m4", "subspace"])
def test_default_args_stationary_stops_quickly(case):
    """phi* = 0 (d < m(m+1)/2): with the default stat_tol = 0 the subproblem must still stop, as stationary, d = 0."""
    rng = np.random.default_rng(42)
    n_stat = 0
    for _ in range(20):
        if case == "free_d3_m4":
            A, kw = _sym(rng, 3, 4, 4), {}
        else:
            A, kw = _sym(rng, 8, 4, 4), {"null_basis": np.linalg.qr(rng.standard_normal((8, 3)))[0]}
        t0 = time.perf_counter()
        r = ascent_direction(A, **kw)
        assert time.perf_counter() - t0 < 2.0 and r.stop in ("gap", "stationary")
        assert r.phi >= 0 and r.phi <= r.phi_upper + 1e-12 * _gscale(A) and r.duality_gap <= 1e-10 * _gscale(A)
        if r.stop == "stationary":
            n_stat += 1
            assert np.all(r.d == 0) and r.phi == 0 and r.phi_upper <= 1e-10 * _gscale(A)
    assert n_stat >= 10                                   # most such instances are stationary


def test_census_no_max_iter():
    rng = np.random.default_rng(2026)
    times, stops = [], []
    for i in range(100):
        kind = i % 4
        if kind == 0:
            args, kw = (_sym(rng, 6, 2, 2),), {}
        elif kind == 1:
            args, kw = (_sym(rng, 4, 2, 2), _sym(rng, 4, 2, 2)), {}
        elif kind == 2:
            args, kw = (_sym(rng, 6, 2, 2),), {"bound_signs": rng.choice([-1.0, 0.0, 1.0], 6)}
        else:
            args, kw = (_sym(rng, 6, 3, 3),), {"null_basis": np.linalg.qr(rng.standard_normal((6, 4)))[0]}
        t0 = time.perf_counter()
        r = ascent_direction(*args, **kw)
        times.append(time.perf_counter() - t0)
        stops.append(r.stop)
        assert r.stop != "max_iter"
        assert r.phi <= r.phi_upper + 1e-12 * _gscale(*args)
        if r.stop in ("gap", "stationary"):
            assert r.duality_gap <= 1e-10 * _gscale(*args)
        if r.stop == "stationary":
            assert np.all(r.d == 0)
    med = float(np.median(times))
    print(f"census: median {med * 1e3:.1f} ms, max {max(times) * 1e3:.1f} ms, stops "
          f"{ {s: stops.count(s) for s in set(stops)} }")
    assert med < 0.05


def test_null_basis_must_be_orthonormal():
    A = np.array([np.diag([1.0, -1.0]), np.diag([1.0, 1.0])])
    with pytest.raises(ValueError, match="orthonormal"):
        ascent_direction(A, null_basis=np.array([[2.0], [0.0]]))


def test_unreachable_tol_stops_no_progress_with_best_certificate():
    """tol = 0 cannot be met in floating point: the run must stop by the stall rule, not spin to max_iter, and return
    the best certified point seen."""
    rng = np.random.default_rng(5)
    stops = []
    for _ in range(6):
        Au, Al = _sym(rng, 4, 2, 2), _sym(rng, 4, 2, 2)
        r = ascent_direction(Au, Al, tol=0.0)
        stops.append(r.stop)
        assert r.stop in ("gap", "no_progress") and r.n_iter < 1000
        assert r.duality_gap <= 1e-12 * _gscale(Au, Al)
    assert "no_progress" in stops
