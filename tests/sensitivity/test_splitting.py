import numpy as np, scipy.sparse as sp, pytest
from bolzakit.sensitivity import cluster_sensitivity, smallest_eigenpairs, splitting_matrices, SplittingMatrices
from bolzakit.sensitivity.errors import ClusterSplitError


def _pencil(seed=0, n=10):
    rng = np.random.default_rng(seed)
    Q = np.linalg.qr(rng.standard_normal((n, n)))[0]
    lam = np.array([0.5, 2.0, 2.0, 2.0] + list(np.linspace(4, 9, n - 4)))
    R = rng.standard_normal((n, n)); M = R @ R.T / n + np.eye(n)
    L = np.linalg.cholesky(M)
    K = L @ (Q @ np.diag(lam) @ Q.T) @ L.T          # K x = lam M x has spectrum lam (exact triple at 2)
    K = (K + K.T) / 2
    d = 4
    dKs = [(lambda S: S + S.T)(rng.standard_normal((n, n))) for _ in range(d)]
    dMs = [(lambda S: (S + S.T) / 10)(rng.standard_normal((n, n))) for _ in range(d)]
    return K, M, dKs, dMs, rng


@pytest.mark.parametrize("sparse", [False, True])
def test_stack_reproduces_cluster_sensitivity(sparse):
    K, M, dKs, dMs, rng = _pencil()
    conv = (lambda X: sp.csc_matrix(X)) if sparse else (lambda X: X)
    sm = splitting_matrices(conv(K), conv(M), [conv(X) for X in dKs], [conv(X) for X in dMs], (1, 3))
    assert isinstance(sm, SplittingMatrices) and sm.mats.shape == (4, 3, 3) and sm.label == "NUMERICAL"
    assert np.allclose(sm.mats, np.transpose(sm.mats, (0, 2, 1)), atol=0)
    for _ in range(5):
        v = rng.standard_normal(4)
        dK = sum(vl * X for vl, X in zip(v, dKs)); dM = sum(vl * X for vl, X in zip(v, dMs))
        r = cluster_sensitivity(conv(K), conv(M), conv(dK), conv(dM), (1, 3))
        Av = np.einsum("l,lab->ab", v, sm.mats)
        # same eigenvectors up to a rotation inside the tie: compare the invariant spectrum and the trace
        assert np.allclose(np.linalg.eigvalsh(Av), r.rates, atol=1e-12 * max(1, np.abs(r.rates).max()))
    assert abs(sm.lam_bar - 2.0) < 1e-12 and sm.gap_ratio >= 1


def test_same_basis_gives_identical_A():
    K, M, dKs, dMs, rng = _pencil(1)
    w, V = smallest_eigenpairs(K, M, 5)
    sm = splitting_matrices(K, M, dKs, dMs, (1, 3), eig=(w, V))
    v = rng.standard_normal(4)
    U = V[:, 1:4]
    dK = sum(vl * X for vl, X in zip(v, dKs)); dM = sum(vl * X for vl, X in zip(v, dMs))
    A = U.T @ (dK - sm.lam_bar * dM) @ U
    assert np.allclose(np.einsum("l,lab->ab", v, sm.mats), (A + A.T) / 2, atol=1e-12)


def test_eig_reuse_is_identical():
    K, M, dKs, dMs, _ = _pencil(2)
    a = splitting_matrices(K, M, dKs, dMs, (1, 3))
    w, V = smallest_eigenpairs(K, M, 5)
    b = splitting_matrices(K, M, dKs, dMs, (1, 3), eig=(w, V))
    assert np.array_equal(a.mats, b.mats) and np.array_equal(a.eigvals, b.eigvals)
    with pytest.raises(ValueError, match="at least 5"):
        splitting_matrices(K, M, dKs, dMs, (1, 3), eig=(w[:4], V[:, :4]))


def test_eig_reuse_still_validates():
    K, M, dKs, dMs, _ = _pencil(3)
    w, V = smallest_eigenpairs(K, M, 5)
    Kbad = K.copy(); Kbad[0, 1] += 1.0
    with pytest.raises(ValueError, match="symmetric"):
        splitting_matrices(Kbad, M, dKs, dMs, (1, 3), eig=(w, V))
    splitting_matrices(Kbad, M, dKs, dMs, (1, 3), eig=(w, V), validate=False)  # no raise


def test_check_gap_false_reports_ratio():
    K = np.diag([1.0, 1.0005, 1.001, 3.0])
    with pytest.raises(ValueError, match="splits"):
        splitting_matrices(K, np.eye(4), [np.eye(4)], None, (0, 1))
    with pytest.raises(ClusterSplitError) as e:                    # exact type
        splitting_matrices(np.diag([1.0, 1.0, 3.0]), np.eye(3), [np.eye(3)], None, (0, 0))
    assert type(e.value) is ClusterSplitError
    sm = splitting_matrices(K, np.eye(4), [np.eye(4)], None, (0, 1), check_gap=False)
    assert sm.gap_ratio < 1 and sm.mats.shape == (1, 2, 2)


def test_dMs_none_means_constant_M():
    K, M, dKs, _, _ = _pencil(4)
    a = splitting_matrices(K, M, dKs, None, (1, 3))
    b = splitting_matrices(K, M, dKs, [np.zeros_like(K)] * 4, (1, 3))
    assert np.allclose(a.mats, b.mats, atol=1e-14)


@pytest.mark.parametrize("sparse", [False, True])
def test_smallest_eigenpairs_count_and_orthonormality(sparse):
    n = 80
    K = sp.diags([-np.ones(n - 1), 2 * np.ones(n), -np.ones(n - 1)], [-1, 0, 1]).tocsc()
    M = sp.diags([1 + 0.2 * np.cos(np.arange(n))], [0]).tocsc()
    if not sparse:
        K, M = K.toarray(), M.toarray()
    w, V = smallest_eigenpairs(K, M, 6)
    assert w.shape == (6,) and V.shape == (n, 6) and np.all(np.diff(w) > 0)
    Md = M.toarray() if sparse else M
    assert np.allclose(V.T @ Md @ V, np.eye(6), atol=1e-10)
    import scipy.linalg as sl
    assert np.allclose(w, sl.eigh(K.toarray() if sparse else K, Md, eigvals_only=True)[:6], rtol=1e-10)
    with pytest.raises(ValueError):
        smallest_eigenpairs(K, M, 0)


def test_mixed_sparse_K_dense_M_review_focus():
    # Review Focus 1: sparse K with a dense M (or the reverse) is accepted and agrees with the all-dense result.
    K, M, dKs, dMs, _ = _pencil(5, n=12)
    a = splitting_matrices(sp.csc_matrix(K), M, dKs, dMs, (1, 3))
    b = splitting_matrices(K, sp.csc_matrix(M), [sp.csc_matrix(X) for X in dKs], dMs, (1, 3))
    c = splitting_matrices(K, M, dKs, dMs, (1, 3))
    for r in (a, b):
        for l in range(4):
            assert np.allclose(np.linalg.eigvalsh(r.mats[l]), np.linalg.eigvalsh(c.mats[l]), atol=1e-10)
    w, V = smallest_eigenpairs(sp.csc_matrix(K), M, 5)
    assert np.allclose(w[1:4], c.eigvals, atol=1e-10)
