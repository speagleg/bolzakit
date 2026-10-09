# tests/sensitivity/test_core.py
import warnings, numpy as np, scipy.sparse as sp, pytest
from bolzakit.sensitivity import cluster_sensitivity, from_callables

def _family(t):
    return np.array([[1 + t, t * t, t], [t * t, 1 - t, 0], [t, 0, 3.0]])

def test_analytic_rates_and_hessian():
    K = _family(0); M = np.eye(3)
    dK = np.array([[1, 0, 1], [0, -1, 0], [1, 0, 0.0]]); d2K = np.array([[0, 2, 0], [2, 0, 0], [0, 0, 0.0]])
    r = cluster_sensitivity(K, M, dK, np.zeros((3, 3)), (0, 1), d2K=d2K, d2M=np.zeros((3, 3)))
    assert np.allclose(sorted(r.rates), [-1, 1])
    assert abs(r.trace_hessian - (-1.0)) < 1e-12        # 2*(1)^2/(1-3) from coupling to lambda=3; tr(U^T K'' U)=0
    assert r.label == "NUMERICAL"

def test_random_pencil_vs_finite_differences():
    rng = np.random.default_rng(0); n = 9
    Q = np.linalg.qr(rng.standard_normal((n, n)))[0]
    D = np.diag([1.0, 2.0, 2.0, 2.0, 3.5, 4, 5, 6, 7]); K0 = Q @ D @ Q.T
    S1 = rng.standard_normal((n, n)); S1 = S1 + S1.T; S2 = rng.standard_normal((n, n)); S2 = S2 + S2.T
    R = rng.standard_normal((n, n)); M0 = R @ R.T / n + np.eye(n); M0 = Q @ np.eye(n) @ Q.T
    K = lambda t: K0 + t * S1 + t * t * S2; Mf = lambda t: M0 + 0.1 * t * np.eye(n)
    r = cluster_sensitivity(K(0), Mf(0), S1, 0.1 * np.eye(n), (1, 3), d2K=2 * S2, d2M=np.zeros((n, n)))
    import scipy.linalg as sl
    h = 1e-4; s = lambda t: sl.eigh(K(t), Mf(t), eigvals_only=True)[1:4].sum()
    fd = (s(h) - 2 * s(0) + s(-h)) / h**2
    assert abs(r.trace_hessian - fd) < 1e-4 * max(1, abs(fd))
    ev = lambda t: np.sort(sl.eigh(K(t), Mf(t), eigvals_only=True)[1:4])
    # Brief had (ev(1e-6) - ev(-1e-6)) / 2e-6. That is wrong: the branches 2 + t*r_k sort in REVERSED order at
    # t < 0, so a sorted central difference yields (r_k + r_{m-1-k})/2, not r_k. For t > 0 the sorted order is
    # the order of the rates, so a forward difference is correct (error O(h)).
    assert np.allclose(np.sort(r.rates), (ev(1e-6) - ev(0)) / 1e-6, atol=1e-4)

def test_kyfan_upper_bounds_prefix_hessian_and_is_tight_with_correctors():
    rng = np.random.default_rng(1); n = 8
    K0 = np.diag([1.0, 2, 2, 4, 5, 6, 7, 8]); S1 = rng.standard_normal((n, n)); S1 += S1.T
    S2 = rng.standard_normal((n, n)); S2 += S2.T; M = np.eye(n)
    r = cluster_sensitivity(K0, M, S1, 0 * M, (1, 2), d2K=2 * S2, d2M=0 * M)
    import scipy.linalg as sl
    h = 1e-4; s = lambda t: sl.eigvalsh(K0 + t * S1 + t * t * S2)[:3].sum()
    fd = (s(h) - 2 * s(0) + s(-h)) / h**2
    assert r.kyfan_upper >= fd - 1e-4 and abs(r.kyfan_upper - fd) < 1e-4
    r2 = cluster_sensitivity(K0, M, S1, 0 * M, (1, 2), d2K=2 * S2, d2M=0 * M,
                             trial=rng.standard_normal((n, 3)) * 0.01)
    assert r2.kyfan_upper >= fd - 1e-6

def test_sparse_matches_dense():
    n = 60; main = 2 + np.arange(n) * 0.0; K = sp.diags([-np.ones(n - 1), 2 * np.ones(n), -np.ones(n - 1)], [-1, 0, 1]).tocsc()
    M = sp.identity(n, format="csc"); dK = sp.diags([np.linspace(0, 1, n)], [0]).tocsc()
    rs = cluster_sensitivity(K, M, dK, 0 * M, (0, 0))
    rd = cluster_sensitivity(K.toarray(), M.toarray(), dK.toarray(), 0 * M.toarray(), (0, 0))
    assert np.allclose(rs.rates, rd.rates)

def test_cluster_splits_degenerate_group_raises():
    with pytest.raises(ValueError, match="splits"):
        cluster_sensitivity(np.diag([1.0, 1.0, 3.0]), np.eye(3), np.eye(3), 0 * np.eye(3), (0, 0))

def test_rejects_nonsymmetric_and_singular():
    K = np.array([[1.0, 2], [0, 1]])
    with pytest.raises(ValueError, match="symmetric"):
        cluster_sensitivity(K, np.eye(2), np.eye(2), 0 * np.eye(2), (0, 0))
    with pytest.raises(ValueError, match="positive definite"):
        cluster_sensitivity(np.eye(2), np.diag([1.0, 0.0]), np.eye(2), 0 * np.eye(2), (0, 0))

def test_from_callables_warns():
    with pytest.warns(UserWarning, match="finite-difference"):
        r = from_callables(lambda p: _family(p[0]), lambda p: np.eye(3), np.zeros(1), np.ones(1), (0, 1))
    assert np.allclose(sorted(r.rates), [-1, 1], atol=1e-6)
