import time
import numpy as np, pytest
import scipy.linalg as sla
from bolzakit.sensitivity.optimize import (ParametricPencil, LowestEigenvalue, Gap, EigenOptResult, IterRecord,
                                           ascent_direction, maximize_eigenvalue)


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


# ------------------------------------------------------------------ maximize_eigenvalue
def _linear_pencil(Ss, M=None):
    Ss = [np.asarray(S, float) for S in Ss]
    n = Ss[0].shape[0]
    Mm = np.eye(n) if M is None else M
    return ParametricPencil(K=lambda p: sum(pl * S for pl, S in zip(p, Ss)), M=lambda p: Mm, dK=lambda p: Ss)


def test_bounds_driver_R1():
    pen = ParametricPencil(K=lambda p: np.diag([p[0], 2 - p[0]]), M=lambda p: np.eye(2),
                           dK=lambda p: [np.diag([1.0, -1.0])])
    seen = []
    res = maximize_eigenvalue(pen, [0.0], LowestEigenvalue(0), lower=[0.0], upper=[2.0], callback=seen.append)
    assert isinstance(res, EigenOptResult) and res.status == "stationary"
    assert len(seen) == len(res.history) > 0 and all(a is b for a, b in zip(seen, res.history))   # callback
    assert abs(res.p[0] - 1) < 1e-3 and abs(res.value - 1) < 1e-3
    assert all(isinstance(h, IterRecord) for h in res.history)


def _simplex_instance(seed):
    rng = np.random.default_rng(seed)
    Ss = []
    for _ in range(3):
        X = rng.standard_normal((3, 3)); Ss.append((X + X.T) / 2 + 2 * np.eye(3))
    return Ss


def _grid_simplex(fun, h=1e-3):
    a = np.arange(0, 1 + h / 2, h)
    A, B = np.meshgrid(a, a, indexing="ij")
    mask = A + B <= 1 + 1e-12
    P = np.stack([A[mask], B[mask], 1 - A[mask] - B[mask]], axis=1)
    vals = fun(P)
    i = int(np.argmax(vals))
    return P[i], vals[i]


def _nm_polish(fun1, p_start):
    from scipy.optimize import minimize
    def neg(x):
        p = np.array([x[0], x[1], 1 - x[0] - x[1]])
        if np.any(p < 0):
            return 1e3 + float(np.sum(np.minimum(p, 0) ** 2))
        return -fun1(p)
    r = minimize(neg, p_start[:2], method="Nelder-Mead",
                 options={"xatol": 1e-12, "fatol": 1e-14, "maxiter": 20000})
    return np.array([r.x[0], r.x[1], 1 - r.x[0] - r.x[1]]), -r.fun


def test_simplex_lowest_eigenvalue_matches_reference():
    Ss = _simplex_instance(3)
    lmin = lambda P: np.linalg.eigvalsh(np.einsum("nl,lab->nab", P, np.array(Ss)))[:, 0]
    p_ref, _ = _grid_simplex(lmin)
    _, v_ref = _nm_polish(lambda p: lmin(p[None])[0], p_ref)
    res = maximize_eigenvalue(_linear_pencil(Ss), np.ones(3) / 3, LowestEigenvalue(0), lower=np.zeros(3),
                              upper=np.ones(3), A_eq=np.ones((1, 3)), b_eq=[1.0], max_iter=200)
    assert abs(res.value - v_ref) < 1e-5, (res.status, res.value, v_ref)
    assert abs(res.p.sum() - 1) < 1e-10 and np.all(res.p >= 0)
    vals = [h.value for h in res.history if h.accepted]
    assert all(b >= a - 1e-14 for a, b in zip(vals, vals[1:]))       # monotone accepted iterates


def test_simplex_gap_polished_from_endpoint():
    Ss = _simplex_instance(4)
    gap = lambda p: (lambda w: w[1] - w[0])(np.linalg.eigvalsh(sum(pl * S for pl, S in zip(p, Ss))))
    res = maximize_eigenvalue(_linear_pencil(Ss), np.ones(3) / 3, Gap(0), lower=np.zeros(3), upper=np.ones(3),
                              A_eq=np.ones((1, 3)), b_eq=[1.0], max_iter=200)
    assert res.status in ("stationary", "step_too_small") and res.Y_lower is not None
    _, v_pol = _nm_polish(gap, res.p)
    assert v_pol - res.value <= 1e-6 and abs(res.p.sum() - 1) < 1e-10 and np.all(res.p >= 0)


def _concave_pencil():
    # lambda_0(p) = p - p^2 (eigenvector e0 for all p near [0, 1]); maximum 1/4 at p = 1/2
    return ParametricPencil(K=lambda p: np.diag([p[0] - p[0] ** 2, 5.0]), M=lambda p: np.eye(2),
                            dK=lambda p: [np.diag([1 - 2 * p[0], 0.0])])


def test_screen_soundness_forced():
    # step0 = 10 overshoots: theta_0 = tau - tau^2 << target, so the screen must fire. _debug_screen eigensolves
    # every screened trial and asserts that the Armijo test would have failed too.
    res = maximize_eigenvalue(_concave_pencil(), [0.0], LowestEigenvalue(0), step0=10.0, _debug_screen=True)
    assert res.n_screened >= 1 and res.status == "stationary"
    assert abs(res.p[0] - 0.5) < 1e-4 and abs(res.value - 0.25) < 1e-8


def test_status_lower_coalescence():
    pen = ParametricPencil(K=lambda p: np.diag([1.0, 1.0 + p[0], 3.0]), M=lambda p: np.eye(3),
                           dK=lambda p: [np.diag([0.0, 1.0, 0.0])])
    assert maximize_eigenvalue(pen, [0.0], LowestEigenvalue(1)).status == "lower_coalescence"


def test_status_gap_closed():
    pen = ParametricPencil(K=lambda p: np.diag([1.0, 1.0 + p[0], 3.0]), M=lambda p: np.eye(3),
                           dK=lambda p: [np.diag([0.0, 1.0, 0.0])])
    assert maximize_eigenvalue(pen, [0.0], Gap(0)).status == "gap_closed"


def test_status_step_too_small():
    # armijo = 1 demands f(tau) >= f + tau * phi, impossible for the strictly concave p - p^2
    res = maximize_eigenvalue(_concave_pencil(), [0.0], LowestEigenvalue(0), armijo=1.0)
    assert res.status == "step_too_small" and res.n_screened > 0


def test_status_max_iter():
    Ss = _simplex_instance(3)
    res = maximize_eigenvalue(_linear_pencil(Ss), np.ones(3) / 3, LowestEigenvalue(0), lower=np.zeros(3),
                              upper=np.ones(3), A_eq=np.ones((1, 3)), b_eq=[1.0], max_iter=1)
    assert res.status == "max_iter" and res.n_iter == 1


def test_status_stationary_at_start():
    pen = ParametricPencil(K=lambda p: np.diag([p[0], 2 - p[0]]), M=lambda p: np.eye(2),
                           dK=lambda p: [np.diag([1.0, -1.0])])
    res = maximize_eigenvalue(pen, [1.0], LowestEigenvalue(0))
    assert res.status == "stationary" and res.multiplicity == 2 and res.n_iter == 0


def test_status_subproblem_not_converged():
    # Adapted from the brief: its 2x2 instance is solved exactly by the subproblem's first polish (stop "stationary"
    # even at max_iter = 1). Here a triple eigenvalue at p = 0 with all four parameters at their lower bound is a
    # stationary point (phi* = 0) that one subproblem iteration cannot certify: ascent_direction stops "max_iter"
    # with d = 0, phi = 0, which the driver must report as subproblem_not_converged, never as stationary.
    X = np.random.default_rng(0).standard_normal((3, 4, 3, 3))[2]
    big = np.zeros((4, 4, 4))
    big[:, :3, :3] = (X + X.transpose(0, 2, 1)) / 2
    pen = ParametricPencil(K=lambda p: np.diag([1.0, 1.0, 1.0, 5.0]) + np.einsum("l,lab->ab", p, big),
                           M=lambda p: np.eye(4), dK=lambda p: list(big))
    res = maximize_eigenvalue(pen, np.zeros(4), LowestEigenvalue(0), lower=0.0, stat_tol=0.0, _subproblem_max_iter=1)
    assert res.status == "subproblem_not_converged" and res.multiplicity == 3
    assert res.phi == 0.0 and res.phi_upper > 0.1                     # uncertified: phi_upper is far from 0
    res = maximize_eigenvalue(pen, np.zeros(4), LowestEigenvalue(0), lower=0.0)      # default budget
    assert res.status == "stationary" and res.phi_upper <= 1e-8


def test_infeasible_p0_raises():
    pen = _linear_pencil(_simplex_instance(3))
    with pytest.raises(ValueError, match="bounds"):
        maximize_eigenvalue(pen, [-0.1, 0.6, 0.5], LowestEigenvalue(0), lower=np.zeros(3))
    with pytest.raises(ValueError, match="A_eq"):
        maximize_eigenvalue(pen, [0.2, 0.2, 0.2], LowestEigenvalue(0), A_eq=np.ones((1, 3)), b_eq=[1.0])


def test_cluster_growth_without_split_error_R3():
    # lambda_j - lambda_k = 0.99e-3 s and lambda_{j+1} - lambda_k = 1.01e-3 s: single linkage chains them.
    # s = max(|lambda_0|, |lambda_1|) ~ 1, so mult_rtol * s ~ 1e-3; the core gap rule would reject cluster (0, 1)
    lam = [1.0, 1.0 + 0.99e-3, 1.0 + 1.01e-3, 1.5, 2.0]
    rng = np.random.default_rng(5)
    S1 = rng.standard_normal((5, 5)); S1 = (S1 + S1.T) / 2
    pen = ParametricPencil(K=lambda p: np.diag(lam) + p[0] * S1, M=lambda p: np.eye(5), dK=lambda p: [S1])
    # refine=False: the fixed-eps clustering itself (eps-refinement would split the cluster at once)
    res = maximize_eigenvalue(pen, [0.0], LowestEigenvalue(0), max_iter=3, refine=False)
    assert res.status in ("stationary", "max_iter", "step_too_small")   # a status, never ClusterSplitError
    assert res.multiplicity == 3
    # with eps-refinement the eps-stationary start is not final: the cluster is split and lambda_0 improves
    ref = maximize_eigenvalue(pen, [0.0], LowestEigenvalue(0), max_iter=3)
    assert ref.status in ("stationary", "max_iter", "step_too_small") and ref.value > res.value
    assert ref.mult_rtol < 1e-3


# ------------------------------------------------------------------ review focus (driver)
def test_sparse_pencil_with_screen_review_focus():
    # Review Focus 2: FE users pass scipy.sparse K, M, dK; the driver (eigensolve, splitting, screen) must accept them.
    import scipy.sparse as sp
    n = 60
    e = np.ones(n)
    L = sp.diags([-e[:-1], 2 * e, -e[:-1]], [-1, 0, 1]).tocsc()
    B = sp.block_diag([L, L], format="csc")                       # every eigenvalue exactly double
    rng = np.random.default_rng(0)
    D = [sp.diags(rng.uniform(0.5, 1.5, 2 * n)).tocsc() for _ in range(3)]
    pen = ParametricPencil(K=lambda p: B + sum(pl * Dl for pl, Dl in zip(p, D)), M=lambda p: sp.identity(2 * n,
                           format="csc"), dK=lambda p: D)
    res = maximize_eigenvalue(pen, np.zeros(3), LowestEigenvalue(0), lower=-1.0, upper=1.0, step0=0.5, max_iter=5)
    assert res.status in ("max_iter", "stationary", "step_too_small") and res.value > res.history[0].value
    assert res.history[0].multiplicity == 2


def test_scalar_bounds_and_integer_p0_review_focus():
    # Review Focus 4: p0 as a Python list of ints and scalar bounds are broadcast, not rejected.
    pen = ParametricPencil(K=lambda p: np.diag([p[0], 2 - p[0]]), M=lambda p: np.eye(2),
                           dK=lambda p: [np.diag([1.0, -1.0])])
    res = maximize_eigenvalue(pen, [0], LowestEigenvalue(0), lower=0, upper=2)
    assert res.status == "stationary" and abs(res.p[0] - 1) < 1e-3 and res.p.dtype == float


def test_gap_with_upper_cluster_at_top_of_spectrum_review_focus():
    # Review Focus 5: Gap(k) whose upper cluster ends at lambda_{n-1}: no out-of-range access, a status is returned.
    pen = ParametricPencil(K=lambda p: np.diag([1.0, 2.0 + p[0], 2.0 + p[1]]), M=lambda p: np.eye(3),
                           dK=lambda p: [np.diag([0.0, 1.0, 0.0]), np.diag([0.0, 0.0, 1.0])])
    res = maximize_eigenvalue(pen, [0.0, 0.0], Gap(0), lower=-0.5, upper=0.5, max_iter=20)
    assert res.multiplicity >= 1 and res.status in ("max_iter", "stationary", "step_too_small")
    assert res.value >= 1.0 - 1e-12


# ------------------------------------------------------------------ scaling (controller ruling, Task 6 review)
def _concave_family(seed):
    """lambda_min of K(p) = A0 + sum p_l S_l - sum p_l^2 D_l (D_l >= 0 diagonal) is concave in p (fix round 2 family);
    M = I for even seeds, a fixed SPD M for odd ones. n = 6, d = 4, box [-1, 1]^4."""
    rng = np.random.default_rng(1000 + seed)
    n, d = 6, 4
    A0 = np.diag(np.sort(rng.uniform(1, 3, n)))
    Ss = [(lambda X: (X + X.T) / 2)(rng.standard_normal((n, n))) for _ in range(d)]
    Ds = [np.diag(rng.uniform(0.0, 1.0, n)) for _ in range(d)]
    if seed % 2:
        X = rng.standard_normal((n, n))
        M0 = X @ X.T / n + np.eye(n)
    else:
        M0 = np.eye(n)
    K = lambda p: A0 + sum(pl * S for pl, S in zip(p, Ss)) - sum(pl * pl * D for pl, D in zip(p, Ds))
    dK = lambda p: [S - 2 * pl * D for pl, S, D in zip(p, Ss, Ds)]
    return K, dK, M0


_SCALES = [(1e6, 1.0), (1e-6, 1.0), (1e6, 1e-6), (1e-8, 1e4), (3.7, 0.21)]


@pytest.mark.parametrize("seed", [3, 0])        # 3: simple optimum on an edge; 0: double eigenvalue (eps-refinement)
def test_driver_scale_invariant(seed):
    # The driver divides the splitting stacks by the eigenvalue scale and every refinement threshold is relative to
    # s, so scaling K by cK and M by cM (eigenvalues by cK / cM) gives the same iterates: the same final value
    # relative to scale and the same p. (Status and screen counts are not compared: they can flip on rounding.)
    Ss = _simplex_instance(seed)
    kw = dict(lower=np.zeros(3), upper=np.ones(3), A_eq=np.ones((1, 3)), b_eq=[1.0], max_iter=200)
    ref = maximize_eigenvalue(_linear_pencil(Ss), np.ones(3) / 3, LowestEigenvalue(0), **kw)
    for cK, cM in _SCALES:
        res = maximize_eigenvalue(_linear_pencil([cK * S for S in Ss], M=cM * np.eye(3)), np.ones(3) / 3,
                                  LowestEigenvalue(0), **kw)
        assert abs(res.value / (cK / cM) - ref.value) <= 1e-10 * abs(ref.value), (cK, cM)
        assert np.allclose(res.p, ref.p, atol=1e-6)


@pytest.mark.parametrize("seed", [2, 7])        # 2: M = I; 7: SPD M. Both optima are double eigenvalues.
def test_driver_scale_invariant_concave(seed):
    K, dK, M0 = _concave_family(seed)
    kw = dict(lower=-1.0, upper=1.0, max_iter=200)
    ref = maximize_eigenvalue(ParametricPencil(K=K, M=lambda p: M0, dK=dK), np.zeros(4), LowestEigenvalue(0), **kw)
    assert ref.eigvals[1] - ref.eigvals[0] <= 1e-6 * abs(ref.value)           # refined onto the crossing
    for cK, cM in _SCALES:
        pen = ParametricPencil(K=lambda p: cK * K(p), M=lambda p: cM * M0, dK=lambda p: [cK * G for G in dK(p)])
        res = maximize_eigenvalue(pen, np.zeros(4), LowestEigenvalue(0), **kw)
        assert abs(res.value / (cK / cM) - ref.value) <= 1e-10 * abs(ref.value), (cK, cM)


def test_power_of_two_scaling_is_bit_identical():
    # Scaling by powers of two is exact in floating point, so no threshold may see any difference at all.
    K, dK, M0 = _concave_family(4)
    kw = dict(lower=-1.0, upper=1.0, max_iter=200)
    ref = maximize_eigenvalue(ParametricPencil(K=K, M=lambda p: M0, dK=dK), np.zeros(4), LowestEigenvalue(0), **kw)
    pen = ParametricPencil(K=lambda p: 1024.0 * K(p), M=lambda p: M0 / 64.0,
                           dK=lambda p: [1024.0 * G for G in dK(p)])
    res = maximize_eigenvalue(pen, np.zeros(4), LowestEigenvalue(0), **kw)
    assert res.value == ref.value * 65536.0 and np.array_equal(res.p, ref.p) and res.n_iter == ref.n_iter


# ------------------------------------------------------------------ fix round 1
def test_rounding_level_armijo_does_not_cycle_I1():
    # lambda_0 = 1 - |p| with p1 >= 0.1: maximum 0.9 at (0.1, 0), smooth there. Once armijo * tau * phi < ulp(f) the
    # Armijo target rounds to f; equal-value steps must not be accepted (they 2-cycled p2 = +-4e-9 to max_iter).
    X = np.array([[0.0, 1.0], [1.0, 0.0]])
    pen = ParametricPencil(K=lambda p: np.eye(2) + p[0] * np.diag([1.0, -1.0]) + p[1] * X, M=lambda p: np.eye(2),
                           dK=lambda p: [np.diag([1.0, -1.0]), X])
    for refine in (True, False):
        res = maximize_eigenvalue(pen, [0.3, 0.4], LowestEigenvalue(0), lower=[0.1, -1.0], upper=[1.0, 1.0],
                                  refine=refine)
        assert res.status in ("stationary", "step_too_small") and res.n_iter <= 60, (res.status, res.n_iter)
        assert abs(res.value - 0.9) < 1e-12 and abs(res.p[0] - 0.1) < 1e-12 and abs(res.p[1]) < 1e-6
        acc = [h.value for h in res.history if h.accepted]
        vals = acc + [res.value]
        assert all(b > a for a, b in zip(vals, vals[1:]))                   # strictly increasing accepted values


def _kelley_simplex(Ss, tol=1e-7, max_cuts=500):
    """Certified bracket LB <= max_{p in simplex} lambda_min(sum p_l S_l) <= UB with UB - LB <= tol: lambda_min is
    concave, so v^T S_l v (v a unit eigenvector of lambda_min) is a supergradient and every cut bounds f from above
    (Kelley's cutting-plane LP gives UB); LB is the best evaluated feasible point."""
    from scipy.optimize import linprog
    A = np.array(Ss)

    def fg(p):
        w, V = np.linalg.eigh(np.einsum("l,lab->ab", p, A))
        return w[0], np.einsum("a,lab,b->l", V[:, 0], A, V[:, 0])

    cuts, LB = [], -np.inf
    p = np.ones(3) / 3
    for _ in range(max_cuts):
        f, g = fg(p)
        cuts.append((f, g, p))
        LB = max(LB, f)
        r = linprog(np.r_[0.0, 0.0, 0.0, -1.0], A_ub=np.array([np.r_[-g, 1.0] for f, g, q in cuts]),
                    b_ub=np.array([f - g @ q for f, g, q in cuts]), A_eq=np.array([[1.0, 1.0, 1.0, 0.0]]),
                    b_eq=[1.0], bounds=[(0, 1)] * 3 + [(None, None)], method="highs")
        UB = -r.fun
        if UB - LB <= tol:
            return LB, UB
        p = np.clip(r.x[:3], 0.0, None)
        p /= p.sum()
    raise AssertionError("Kelley did not converge")


@pytest.mark.parametrize("seed", [0, 11])
def test_double_eigenvalue_optimum_eps_refinement_I2(seed):
    # Simplex instances whose maximiser has a double lambda_min. At mult_rtol = 1e-3 alone the driver is stationary
    # for the eps-reduced model only, 1e-4 short; the eps-refinement must reach the certified bracket to 1e-5.
    Ss = _simplex_instance(seed)
    LB, UB = _kelley_simplex(Ss)
    kw = dict(lower=np.zeros(3), upper=np.ones(3), A_eq=np.ones((1, 3)), b_eq=[1.0], max_iter=200)
    res = maximize_eigenvalue(_linear_pencil(Ss), np.ones(3) / 3, LowestEigenvalue(0), **kw)
    assert res.value <= UB + 1e-12 and LB - res.value <= 1e-5, (res.status, LB - res.value)
    # the refined run ends "stationary" or, once the tied model stalls at the finest eps, "step_too_small"
    assert res.status in ("stationary", "step_too_small") and res.mult_rtol < 1e-3
    assert res.eigvals[1] - res.eigvals[0] <= 1e-6 * abs(res.value)          # on the crossing
    unrefined = maximize_eigenvalue(_linear_pencil(Ss), np.ones(3) / 3, LowestEigenvalue(0), refine=False, **kw)
    assert unrefined.status == "stationary" and LB - unrefined.value > 1e-5     # the defect the refinement fixes


def test_exact_crossing_1d_I2():
    # lambda_0(p) = min(1.3 + p/2, 1.55 - p/2): an exact crossing (kink) at p = 0.25 with value 1.425
    pen = ParametricPencil(K=lambda p: np.diag([1.3 + 0.5 * p[0], 1.55 - 0.5 * p[0]]), M=lambda p: np.eye(2),
                           dK=lambda p: [np.diag([0.5, -0.5])])
    res = maximize_eigenvalue(pen, [0.0], LowestEigenvalue(0), lower=0.0, upper=1.0)
    assert res.status == "stationary" and res.multiplicity == 2
    assert abs(res.value - 1.425) < 1e-8 and abs(res.p[0] - 0.25) < 1e-7


@pytest.mark.parametrize("obj", [Gap(2), Gap(-1), LowestEigenvalue(3), LowestEigenvalue(-1)])
def test_objective_k_out_of_range_M4(obj):
    pen = ParametricPencil(K=lambda p: np.diag([1.0, 2.0, 3.0 + p[0]]), M=lambda p: np.eye(3),
                           dK=lambda p: [np.diag([0.0, 0.0, 1.0])])
    with pytest.raises(ValueError, match="out of range"):
        maximize_eigenvalue(pen, [0.0], obj)


def test_direction_hook_nonpositive_M5():
    pen = ParametricPencil(K=lambda p: np.diag([p[0], 2 - p[0]]), M=lambda p: np.eye(2),
                           dK=lambda p: [np.diag([1.0, -1.0])])
    res = maximize_eigenvalue(pen, [0.0], LowestEigenvalue(0), _direction=lambda sm, N, sg: (np.zeros(1), 0.0))
    assert res.status == "direction_nonpositive"


def test_screen_trial_assembles_once_M3():
    calls = {"K": 0}
    base = _concave_pencil()

    def K(p):
        calls["K"] += 1
        return base.K(p)

    pen = ParametricPencil(K=K, M=base.M, dK=base.dK)
    res = maximize_eigenvalue(pen, [0.0], LowestEigenvalue(0), step0=10.0)
    assert res.n_screened >= 1
    # one assembly per trial (screened or eigensolved) plus the start; none twice
    trials = sum(1 for h in res.history)
    assert calls["K"] == trials + 1, (calls["K"], trials)


# ------------------------------------------------------------------ fix round 2
def _concave_kelley(K, dK, M0, tol=1e-7, max_cuts=600):
    """Certified bracket for max over [-1, 1]^4 of the concave lambda_min (cuts from v^T dK_l v, v M-normalised)."""
    from scipy.optimize import linprog
    cuts, LB, p = [], -np.inf, np.zeros(4)
    for _ in range(max_cuts):
        w, V = sla.eigh(K(p), M0)
        f, g = w[0], np.array([V[:, 0] @ G @ V[:, 0] for G in dK(p)])
        cuts.append((f, g, p.copy()))
        LB = max(LB, f)
        r = linprog(np.r_[np.zeros(4), -1.0], A_ub=np.array([np.r_[-g, 1.0] for f, g, q in cuts]),
                    b_ub=np.array([f - g @ q for f, g, q in cuts]), bounds=[(-1, 1)] * 4 + [(None, None)],
                    method="highs")
        if -r.fun - LB <= tol:
            return LB, -r.fun
        p = r.x[:4]
    raise AssertionError("Kelley did not converge")


@pytest.mark.parametrize("seed", [0, 10])
def test_concave_double_optimum_refines_on_stall(seed):
    # The eps-model stalls at phi ~ 7e-8 > stat_tol * s, so the trust radius collapses; these runs used to end
    # step_too_small ~1e-4 short without ever refining. Refinement on the stall must reach the certified bracket.
    K, dK, M0 = _concave_family(seed)
    LB, UB = _concave_kelley(K, dK, M0)
    res = maximize_eigenvalue(ParametricPencil(K=K, M=lambda p: M0, dK=dK), np.zeros(4), LowestEigenvalue(0),
                              lower=-1.0, upper=1.0, max_iter=200)
    assert res.value <= UB + 1e-12 and LB - res.value <= 1e-5, (res.status, LB - res.value)
    assert any(h.refined for h in res.history)


def test_refinement_respects_eigensolve_cap_and_is_recorded():
    lam = [1.0, 1.0 + 0.99e-3, 1.0 + 1.01e-3, 1.5, 2.0]
    S1 = (lambda X: (X + X.T) / 2)(np.random.default_rng(5).standard_normal((5, 5)))
    pen = ParametricPencil(K=lambda p: np.diag(lam) + p[0] * S1, M=lambda p: np.eye(5), dK=lambda p: [S1])
    res = maximize_eigenvalue(pen, [0.0], LowestEigenvalue(0), _max_eigensolves=1)    # eps-stationary start
    assert res.status == "max_iter" and res.n_eigensolves == 1
    seen = []
    res = maximize_eigenvalue(pen, [0.0], LowestEigenvalue(0), max_iter=5, callback=seen.append)
    refined = [h for h in res.history if h.refined]
    assert refined and refined[0].it == 0 and refined[0].tau == 0.0 and not refined[0].accepted
    assert len(seen) == len(res.history) and all(a is b for a, b in zip(seen, res.history))
