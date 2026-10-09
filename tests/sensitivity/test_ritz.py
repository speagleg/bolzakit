import numpy as np, scipy.linalg as sl, scipy.sparse as sp, pytest
from bolzakit.sensitivity import ritz_upper_bounds
from bolzakit.sensitivity.errors import PencilError


def test_poincare_property_200_random_pencils():
    rng = np.random.default_rng(0)
    for _ in range(200):
        n = int(rng.integers(3, 31)); p = int(rng.integers(1, n + 1))
        S = rng.standard_normal((n, n)); K = S + S.T
        R = rng.standard_normal((n, n)); M = R @ R.T + 0.5 * np.eye(n)
        W = rng.standard_normal((n, p))
        lam = sl.eigh(K, M, eigvals_only=True)
        th = ritz_upper_bounds(K, M, W)
        scale = max(1.0, np.abs(lam).max())
        assert th.shape == (p,) and np.all(np.diff(th) >= -1e-12 * scale)
        assert np.all(lam[:p] <= th + 1e-12 * scale)
        assert lam[:p].sum() <= np.trace(np.linalg.solve(W.T @ M @ W, W.T @ K @ W)) + 1e-10 * scale * p


def test_equality_on_eigenvectors():
    rng = np.random.default_rng(1); n = 12
    S = rng.standard_normal((n, n)); K = S + S.T
    R = rng.standard_normal((n, n)); M = R @ R.T + np.eye(n)
    lam, V = sl.eigh(K, M)
    W = V[:, :4] @ np.linalg.qr(rng.standard_normal((4, 4)))[0]   # same span, rotated basis
    assert np.allclose(ritz_upper_bounds(K, M, W), lam[:4], atol=1e-10 * np.abs(lam).max())


def test_sparse_inputs():
    n = 40
    K = sp.diags([-np.ones(n - 1), 2 * np.ones(n), -np.ones(n - 1)], [-1, 0, 1]).tocsc()
    M = sp.identity(n, format="csc")
    W = np.random.default_rng(2).standard_normal((n, 3))
    assert np.allclose(ritz_upper_bounds(K, M, W), ritz_upper_bounds(K.toarray(), M.toarray(), W))


def test_rank_deficient_raises():
    W = np.zeros((5, 2)); W[0, 0] = W[0, 1] = 1.0
    with pytest.raises(PencilError, match="full column rank") as e:
        ritz_upper_bounds(np.eye(5), np.eye(5), W)
    assert type(e.value) is PencilError
