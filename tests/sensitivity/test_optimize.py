import time
import numpy as np, pytest
import scipy.linalg as sla
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


def _plate_like(rng, m=4):
    """Area equality (null_basis of ones, non-square) plus ~90% of the variables at bounds: the Task 7/9 plate set-up."""
    return (_sym(rng, 64, m, m),), {"null_basis": sla.null_space(np.ones((1, 64))),
                                    "bound_signs": rng.choice([-1.0, 0.0, 1.0], 64, p=[0.45, 0.1, 0.45])}


def test_census_no_max_iter():
    """Correctness census. Wall-clock limits are loose guards (~10x the unloaded times, so shared CI runners pass);
    the strict spec 13 targets (< 50 ms) are measured by the benchmark / perf docs, not asserted. The load-independent
    guard is the iteration count."""
    times, stops, iters = [], [], []
    instances = []
    rng = np.random.default_rng(2026)
    for i in range(112):
        instances.append((rng, i % 4 if i < 100 else 4 + i % 3))
    for seed in (47, 101, 7):                     # round-4 kind: small bounds + non-square subspace + Gap
        srng = np.random.default_rng(seed)
        instances += [(srng, 7)] * 30
    for rng, kind in instances:
        if kind == 0:
            args, kw = (_sym(rng, 6, 2, 2),), {}
        elif kind == 1:
            args, kw = (_sym(rng, 4, 2, 2), _sym(rng, 4, 2, 2)), {}
        elif kind == 2:
            args, kw = (_sym(rng, 6, 2, 2),), {"bound_signs": rng.choice([-1.0, 0.0, 1.0], 6)}
        elif kind == 3:
            args, kw = (_sym(rng, 6, 3, 3),), {"null_basis": np.linalg.qr(rng.standard_normal((6, 4)))[0]}
        elif kind == 4:
            args, kw = (_sym(rng, 64, 4, 4),), {"bound_signs": rng.choice([-1.0, 0.0, 1.0], 64, p=[0.45, 0.1, 0.45])}
        elif kind == 5:
            args, kw = (_sym(rng, 64, 2, 2), _sym(rng, 64, 2, 2)), {"bound_signs": rng.choice([-1.0, 0.0, 1.0], 64)}
        elif kind == 6:
            args, kw = _plate_like(rng)
        else:
            args, kw = _small_bounds_subspace_gap(rng)
        t0 = time.perf_counter()
        r = ascent_direction(*args, **kw)
        times.append(time.perf_counter() - t0)
        stops.append(r.stop)
        iters.append(r.n_iter)
        assert r.stop != "max_iter"
        _assert_certificate(r, *args, **kw)
    med = float(np.median(times))
    print(f"census: median {med * 1e3:.1f} ms, max {max(times) * 1e3:.1f} ms, iterations median "
          f"{np.median(iters):.0f} max {max(iters)}, stops { {s: stops.count(s) for s in set(stops)} }")
    assert np.median(iters) <= 100 and max(iters) <= 2000
    assert med < 0.25


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


@pytest.mark.parametrize("kind", ["m4_bounds", "gap_bounds", "plate_m4"])
def test_d64_with_bounds_fast(kind):
    """d = 64, m <= 4 with many variables at bounds (and, for plate_m4, the area equality via a non-square
    null_basis): the configuration of Task 9. The strict spec 13 target (< 50 ms) is measured by the benchmark /
    perf docs, not asserted; here wall clock is only a loose guard (~10x the unloaded times, for shared CI
    runners) and the load-independent guard is the iteration count."""
    rng = np.random.default_rng(64)
    times, iters = [], []
    for _ in range(15):
        if kind == "m4_bounds":
            args, kw = (_sym(rng, 64, 4, 4),), {"bound_signs": rng.choice([-1.0, 0.0, 1.0], 64, p=[0.45, 0.1, 0.45])}
        elif kind == "gap_bounds":
            args, kw = (_sym(rng, 64, 2, 2), _sym(rng, 64, 2, 2)), {"bound_signs": rng.choice([-1.0, 0.0, 1.0], 64)}
        else:
            args, kw = _plate_like(rng)
        t0 = time.perf_counter()
        r = ascent_direction(*args, **kw)
        times.append(time.perf_counter() - t0)
        iters.append(r.n_iter)
        assert r.stop in ("gap", "stationary") and r.duality_gap <= 1e-10 * _gscale(*args)
        assert np.all(kw["bound_signs"] * r.d >= 0)
        N = kw.get("null_basis")
        if N is not None:
            assert np.linalg.norm(r.d - N @ (N.T @ r.d)) <= 1e-12
    med = float(np.median(times))
    print(f"d=64 {kind}: median {med * 1e3:.1f} ms, max {max(times) * 1e3:.1f} ms, iterations median "
          f"{np.median(iters):.0f} max {max(iters)}")
    assert np.median(iters) <= 100 and max(iters) <= 2000
    assert med < 0.25 and max(times) < 2.5


# ------------------------------------------------------------------ fix round 4: full-rank faces, certificate checks
def _assert_certificate(r, Au, Al=None, null_basis=None, bound_signs=None, stat_tol=0.0, f_value=0.0, tol=1e-10):
    """Recompute every certificate property from the returned fields: d unit-or-zero, in span N, sign feasible;
    phi = lambda_min(A_upper(d)) - lambda_max(A_lower(d)); Y on the spectraplex, mu >= 0 on the bounds only;
    phi_upper = |P(g(Y) + mu s)|; weak duality phi <= phi_upper; the stop contract."""
    dim = Au.shape[0]
    g = _gscale(Au, *([] if Al is None else [Al]))
    N = np.eye(dim) if null_basis is None else null_basis
    s = np.zeros(dim) if bound_signs is None else bound_signs
    nd = np.linalg.norm(r.d)
    assert nd == 0 or abs(nd - 1) <= 1e-12
    assert np.linalg.norm(r.d - N @ (N.T @ r.d)) <= 1e-10
    assert np.all(s * r.d >= 0)
    lm = np.linalg.eigvalsh(np.einsum("l,lab->ab", r.d, Au))[0]
    if Al is not None:
        lm -= np.linalg.eigvalsh(np.einsum("l,lab->ab", r.d, Al))[-1]
    assert abs(lm - r.phi) <= 1e-12 * g
    for Y in (r.Y_upper, r.Y_lower):
        if Y is not None:
            assert np.allclose(Y, Y.T) and np.linalg.eigvalsh(Y)[0] >= -1e-13 and abs(np.trace(Y) - 1) <= 1e-12
    assert np.all(r.mu >= 0) and np.all(r.mu[s == 0] == 0)
    gv = np.einsum("lab,ab->l", Au, r.Y_upper) + r.mu * s
    if Al is not None:
        gv -= np.einsum("lab,ab->l", Al, r.Y_lower)
    assert abs(np.linalg.norm(N @ (N.T @ gv)) - r.phi_upper) <= 1e-12 * g
    assert r.phi <= r.phi_upper + 1e-12 * g
    assert abs(r.duality_gap - max(r.phi_upper - r.phi, 0.0)) <= 1e-15 * g
    if r.stop == "gap":
        assert r.phi > 0 and nd > 0 and r.duality_gap <= tol * g
    if r.stop == "stationary":
        assert np.all(r.d == 0) and r.phi == 0
        assert r.phi_upper <= max(tol * g, stat_tol * max(1.0, abs(f_value)))


def _small_bounds_subspace_gap(rng, all_bounded=True):
    """d = 8, a random 5-dim null_basis, a Gap stack, bounds on every variable (or on a random subset)."""
    signs = rng.choice([-1.0, 1.0], 8) if all_bounded else rng.choice([-1.0, 0.0, 1.0], 8)
    return (_sym(rng, 8, 2, 2), _sym(rng, 8, 2, 2)), {"null_basis": np.linalg.qr(rng.standard_normal((8, 5)))[0],
                                                     "bound_signs": signs}


def test_project_face_full_rank_returns_exact_zeros():
    """The fixed coordinates pin all of span N (rank N[fixed] = q): the projection is exactly 0, not rounding noise."""
    from bolzakit.sensitivity.optimize import _feasible_direction, _project_face
    rng = np.random.default_rng(3)
    for d, q, nfix in ((8, 5, 5), (8, 5, 7), (64, 10, 40), (12, 3, 12)):
        N = np.linalg.qr(rng.standard_normal((d, q)))[0]
        fixed = np.zeros(d, bool)
        fixed[rng.choice(d, nfix, replace=False)] = True
        for scale in (1.0, 1e-15, 1e8):
            r = scale * rng.standard_normal(d)
            out = _project_face(N, r, fixed)
            assert np.all(out == 0)
            signs = np.where(fixed, 1.0, 0.0)
            assert np.all(_feasible_direction(N, r, signs, signs) == 0)
    # rank < q: a genuine projection, in span N with the fixed coordinates at rounding level
    N = np.linalg.qr(rng.standard_normal((8, 5)))[0]
    fixed = np.zeros(8, bool)
    fixed[:3] = True
    out = _project_face(N, rng.standard_normal(8), fixed)
    assert np.linalg.norm(out) > 0.1 and np.linalg.norm(out - N @ (N.T @ out)) <= 1e-14
    assert np.linalg.norm(out[fixed]) <= 1e-14


def _sub_bounds_gap_seed47_instance0():
    """The review's reproducer: the reviewer census (seed 47) draws 50 'sub_bounds' instances, then this one."""
    rng = np.random.default_rng(47)
    for _ in range(50):
        d = int(rng.choice([6, 8, 12]))
        _sym(rng, d, 4, 4)
        np.linalg.qr(rng.standard_normal((d, d - 3)))
        rng.choice([-1.0, 0.0, 1.0], d, p=[0.3, 0.4, 0.3])
    args = (_sym(rng, 8, 2, 2), _sym(rng, 8, 2, 2))
    kw = {"null_basis": np.linalg.qr(rng.standard_normal((8, 5)))[0], "bound_signs": rng.choice([-1.0, 0.0, 1.0], 8)}
    return args, kw


def test_full_rank_face_reproducer_is_certified():
    """Round-3 regression: the QR face projection returned ~1e-31 noise on a face where the fixed coordinates pin all
    of span N; it was normalised into a 'unit' d outside span N and returned as stop='gap' with phi = 0.31 >
    phi_upper = 7e-16."""
    args, kw = _sub_bounds_gap_seed47_instance0()
    assert np.all(kw["bound_signs"] != 0)                   # the instance as reported: every variable at a bound
    r = ascent_direction(*args, **kw)
    _assert_certificate(r, *args, **kw)
    assert r.stop in ("gap", "stationary")


def test_uncertified_direction_is_never_returned(monkeypatch):
    """Soundness guard: even if the face projection produced a unit d outside span N (the round-3 failure mode),
    ascent_direction must return a certified result (here d = 0, phi = 0), never an uncertified 'gap'."""
    import bolzakit.sensitivity.optimize as O
    args, kw = _sub_bounds_gap_seed47_instance0()
    off_span = np.zeros(8)
    off_span[[2, 5]] = [-0.35498954, -0.93487027]           # the direction round 3 returned: rate 0.31 > phi_upper
    off_span /= np.linalg.norm(off_span)
    assert np.linalg.norm(off_span - kw["null_basis"] @ (kw["null_basis"].T @ off_span)) > 0.1
    monkeypatch.setattr(O, "_feasible_direction", lambda N, r, signs, mu: off_span.copy())
    r = ascent_direction(*args, **kw)
    _assert_certificate(r, *args, **kw)
    assert np.all(r.d == 0) and r.phi == 0 and r.stop in ("stationary", "no_progress")
