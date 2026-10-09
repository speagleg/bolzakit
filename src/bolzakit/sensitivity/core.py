"""Repeated-eigenvalue (cluster) sensitivity for a symmetric pencil K(t) u = lambda M(t) u.

The pencil comes from any finite-element model. It may be dense (numpy) or sparse (scipy.sparse). Given the
first (and optionally second) parameter derivatives of K and M, this module returns, for a cluster of eigenvalues
lambda_i..lambda_j:

* ``A = U^T (dK - lam_bar dM) U``. Its eigenvalues are the first-order splitting rates of an exactly degenerate
  cluster (for a nearly degenerate one, see "Near-degenerate clusters" below).
* ``trace_hessian``: the second derivative of the sum of the cluster eigenvalues,
  ``tr(U^T (d2K - lam_bar d2M) U) - 2 tr(A U^T dM U) - 2 tr(B^T Z)``. Here ``B = (dK - lam_bar dM) U`` and Z is
  the solution of the deflated system ``(K - lam_bar M) Z = B - M U U^T B`` with ``U^T M Z = 0``.
* ``kyfan_upper``: an upper bound on the second derivative of the PREFIX sum lambda_0 + ... + lambda_j.

Caveat (discrete vs. continuum): every quantity here belongs to the DISCRETE pencil (K, M) that you pass in. The
Ky Fan bound and the derivatives are statements about the eigenvalues of the matrices. They are not statements
about the eigenvalues of the continuous operator that the finite-element model discretises, and the
discretisation error is not controlled. All results carry the label "NUMERICAL".

Near-degenerate clusters. If the cluster is only nearly degenerate (spread g_in > 0), the formulas use the cluster
mean lam_bar, and ``rates`` are by definition eigvalsh(A) at lam_bar: the slopes of the reduced model
    lambda_{i+r}(p + t v) = lam_bar + eig_r(Lambda_C - lam_bar I + t A) + O(t^2 + t g_in),
i.e. the branch slopes once |t| >> g_in / ||A||. They are NOT O(g_in)-accurate predictions of the t = 0 derivatives
of the individual branches. Counterexample: K = diag(1, 1 + delta, 3), M = I, and dK = e1 e2^T + e2 e1^T. The rates
are +-1 for every delta > 0, but both branch derivatives at t = 0 are 0. The cluster SUM is smooth, so ``tr A`` (its
first derivative) and ``trace_hessian`` (its second derivative) are correct: exact for an exactly degenerate
cluster, and accurate to O(g_in) otherwise (exact in the counterexample).

This module imports only numpy and scipy, and nothing from the rest of bolzakit.
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass

import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from bolzakit.sensitivity.errors import (ClusterRangeError, ClusterSplitError, FiniteDifferenceWarning,
                                         PencilError)

__all__ = ["ClusterSensitivity", "cluster_sensitivity", "directional", "from_callables", "smallest_eigenpairs",
           "SplittingMatrices", "splitting_matrices", "ritz_upper_bounds", "DISCRETE_CAVEAT"]

_PENCIL_MSG = "K and M must be symmetric with M positive definite"

DISCRETE_CAVEAT = ("NUMERICAL: computed for the discrete pencil (K, M) as given; discretisation error relative to "
                   "any continuum operator is not controlled.")


@dataclass
class ClusterSensitivity:
    eigvals: np.ndarray          # cluster eigenvalues (ascending)
    basis: np.ndarray            # n x m, M-orthonormal cluster eigenvectors
    A: np.ndarray                # m x m symmetric first-order matrix U^T (dK - lam_bar dM) U
    rates: np.ndarray            # eigvalsh(A) at lam_bar (not branch derivatives of a near-tie; module docstring)
    trace_hessian: float | None  # d^2/dt^2 of the sum of the cluster eigenvalues (needs d2K/d2M)
    kyfan_upper: float | None    # upper bound for d^2/dt^2 (lambda_0 + ... + lambda_j), discrete problem only
    label: str = "NUMERICAL"
    cluster: tuple[int, int] | None = None   # (i, j) as passed
    lam_bar: float | None = None             # mean of the cluster eigenvalues
    gap_ratio: float | None = None           # g_out / max(gap_ratio_min |lam_bar|, 10 g_in); >= 1 when accepted
    caveat: str = DISCRETE_CAVEAT


# ----------------------------------------------------------------------------------------------- helpers
def _maxabs(X) -> float:
    if sp.issparse(X):
        return float(abs(X).max()) if X.nnz else 0.0
    return float(np.max(np.abs(X))) if X.size else 0.0


def _is_symmetric(X) -> bool:
    return _maxabs(X - X.T) <= 1e-12 * _maxabs(X)


def _sparse_pd_by_lu(X) -> bool:
    """Certify positive definiteness by an LU with diagonal pivoting only, under a symmetric ordering.

    If perm_r == perm_c, the factorisation is P X P^T = L U with L unit lower triangular. For symmetric X this is
    L D L^T with D = diag(U), and Sylvester's law makes sign(D) the inertia. SuperLU still swaps rows off the
    diagonal when a diagonal pivot is exactly zero, even with diag_pivot_thresh=0 (e.g. [[0,1],[1,0]] blocks).
    Then perm_r != perm_c, diag(U) says nothing about the inertia, and the result is "not certified" (False)."""
    try:
        lu = spla.splu(sp.csc_matrix(X), permc_spec="MMD_AT_PLUS_A", diag_pivot_thresh=0.0,
                       options={"SymmetricMode": True})
    except RuntimeError:  # exactly singular
        return False
    if not np.array_equal(lu.perm_r, lu.perm_c):
        return False
    d = lu.U.diagonal()
    return bool(np.all(np.isfinite(d)) and np.all(d > 0))


def _validate(K, M, sparse: bool) -> None:
    if K.shape != M.shape or K.shape[0] != K.shape[1]:
        raise PencilError(_PENCIL_MSG)
    if not (_is_symmetric(K) and _is_symmetric(M)):
        raise PencilError(_PENCIL_MSG)
    if sparse:
        ok = _sparse_pd_by_lu(M)
        if ok:
            try:
                mu = spla.eigsh(M, k=1, sigma=0, which="LM", return_eigenvectors=False)[0]
            except Exception:  # noqa: BLE001 - any failure here means we cannot certify PD
                mu = -1.0
            ok = mu > 0
    else:
        try:
            np.linalg.cholesky(M)
            ok = True
        except np.linalg.LinAlgError:
            ok = False
    if not ok:
        raise PencilError(_PENCIL_MSG)


def _eigs(K, M, count: int, sparse: bool):
    """The smallest ``count`` eigenpairs (or all of them, if n is smaller), ascending, M-orthonormal."""
    n = K.shape[0]
    count = min(int(count), n)
    if not sparse or count >= n - 1:
        Kd = K.toarray() if sp.issparse(K) else K
        Md = M.toarray() if sp.issparse(M) else M
        w, V = sla.eigh(Kd, Md)
        return w[:count], V[:, :count]
    v0 = np.random.default_rng(0).standard_normal(n)  # fixed start vector: reproducible runs (ARPACK's is random)
    # Find a shift below lambda_0: K - s M must be positive definite (checked by inertia).
    scale = max(_maxabs(K), 1e-300)
    s = -1e-8 * scale
    while not _sparse_pd_by_lu(K - s * M):
        s = 4.0 * s
        if s < -1e12 * scale:
            raise PencilError(_PENCIL_MSG)
    w, V = spla.eigsh(K, k=count, M=M, sigma=s, which="LM", v0=v0)
    o = np.argsort(w)
    w, V = w[o], V[:, o]
    # Move the shift to just below lambda_0, then re-solve.
    s2 = w[0] - 1e-6 * max(abs(w[0]), 1e-8 * scale)
    if s2 > s and _sparse_pd_by_lu(K - s2 * M):
        w, V = spla.eigsh(K, k=count, M=M, sigma=s2, which="LM", v0=v0)
        o = np.argsort(w)
        w, V = w[o], V[:, o]
    # M-orthonormalise. eigsh does this only up to its tolerance, and it does not do it inside degenerate groups.
    G = V.T @ (M @ V)
    L = np.linalg.cholesky((G + G.T) / 2)
    V = sla.solve_triangular(L, V.T, lower=True).T
    return w, V


def _prepare(K, M):
    """Normalise (K, M): sparse K selects the sparse path (both CSC); otherwise both dense float arrays."""
    sparse = sp.issparse(K)
    if sparse:
        return sp.csc_matrix(K), sp.csc_matrix(M), True
    K = np.asarray(K, dtype=float)
    M = M.toarray() if sp.issparse(M) else np.asarray(M, dtype=float)
    return K, M, False


def _check_range(cluster, n: int) -> tuple[int, int]:
    i, j = (int(c) for c in cluster)
    if not (0 <= i <= j < n):
        raise ClusterRangeError(f"cluster {tuple(cluster)} out of range for n = {n}")
    return i, j


def _gap(w, i: int, j: int, gap_ratio_min: float):
    """Core gap rule (spec 5.1). Returns (lam_bar, g_in, g_out, thresh, ratio)."""
    lam = w[i:j + 1]
    lam_bar = float(lam.mean())
    g_in = float(lam[-1] - lam[0])
    outside = []
    if i > 0:
        outside.append(lam[0] - w[i - 1])
    if j + 1 < len(w):
        outside.append(w[j + 1] - lam[-1])
    g_out = float(min(outside)) if outside else np.inf
    thresh = max(gap_ratio_min * abs(lam_bar), 10 * g_in)
    ratio = g_out / thresh if thresh > 0 else np.inf
    return lam_bar, g_in, g_out, thresh, float(ratio)


def _split_error(cluster, g_in, g_out, ratio):
    return ClusterSplitError(f"cluster {tuple(cluster)} splits a (near-)degenerate group: gap ratio "
                             f"g_out / max(gap_ratio_min*|lam_bar|, 10*g_in) = {ratio:.3g} < 1 "
                             f"(g_in = {g_in:.3g}, g_out = {g_out:.3g})")


def smallest_eigenpairs(K, M, count: int) -> tuple[np.ndarray, np.ndarray]:
    """The ``count`` smallest eigenpairs of K x = lambda M x, ascending, with M-orthonormal eigenvectors.

    Dense inputs use scipy.linalg.eigh; sparse K uses shift-invert eigsh with an inertia-certified shift below
    lambda_0. Raises PencilError if K, M are not symmetric with M positive definite."""
    K, M, sparse = _prepare(K, M)
    n = K.shape[0]
    if not 1 <= int(count) <= n:
        raise ValueError(f"count must lie in [1, {n}]")
    _validate(K, M, sparse)
    return _eigs(K, M, int(count), sparse)


def _as_op(X, sparse: bool, n: int):
    if X is None:
        return sp.csc_matrix((n, n)) if sparse else np.zeros((n, n))
    if sparse:
        return sp.csc_matrix(X)
    return X.toarray() if sp.issparse(X) else np.asarray(X, dtype=float)


def _deflated_solve(K, M, shift: float, W, RHS, sparse: bool):
    """Solve (K - shift M) Z = RHS with W^T M Z = 0, via the bordered system [[K - shift M, M W], [W^T M, 0]].
    This equals the brief's deflated solve (K - shift M + M W W^T M) Z = RHS when RHS is M-projected off W, and it
    is nonsingular whenever shift lies in the prefix spanned by W. RHS must satisfy W^T RHS = 0."""
    n, p = W.shape
    MW = np.asarray(M @ W)
    rhs = np.vstack([RHS, np.zeros((p, RHS.shape[1]))])
    if sparse:
        S = sp.bmat([[K - shift * M, sp.csc_matrix(MW)], [sp.csc_matrix(MW.T), None]], format="csc")
        sol = spla.spsolve(S, rhs)
        sol = sol.reshape(n + p, -1) if sol.ndim == 1 else np.asarray(sol)
    else:
        S = np.block([[K - shift * M, MW], [MW.T, np.zeros((p, p))]])
        sol = np.linalg.solve(S, rhs)
    return sol[:n]


def _proj(M, W, B):
    """B - M W (W^T B): removes the M W components, so that W^T result = 0 (for W M-orthonormal)."""
    return B - np.asarray(M @ (W @ (W.T @ B)))


# ----------------------------------------------------------------------------------------------- public API
def cluster_sensitivity(K, M, dK, dM, cluster, d2K=None, d2M=None, trial=None,
                        gap_ratio_min=1e-3) -> ClusterSensitivity:
    """Sensitivity of the eigenvalue cluster ``cluster = (i, j)`` (inclusive, 0-based, ascending order) of the
    symmetric-definite pencil K u = lambda M u, along one parameter direction.

    ``dK, dM`` are the first derivatives and ``d2K, d2M`` the second derivatives along that direction. If both
    second derivatives are None, ``trace_hessian`` and ``kyfan_upper`` are None. If only one is given, the other is
    taken as zero. Inputs may be dense arrays or scipy.sparse matrices. A sparse K selects the sparse path
    (shift-invert ``eigsh`` and ``spsolve``).

    ``kyfan_upper`` bounds d^2/dt^2 of the PREFIX sum lambda_0 + ... + lambda_j. It is valid for the DISCRETE
    pencil only. Let W be all eigenvectors 0..j and Psi a trial direction (``trial``, an n x (j+1) array, or by
    default the exact first-order correctors -Z_a of each prefix eigenvector). By the Ky Fan principle,
    sum_{k<=j} lambda_k(t) <= g(t) = tr[(W_t^T M(t) W_t)^{-1} W_t^T K(t) W_t] with W_t = W + t Psi, and equality
    holds at t = 0. Hence g''(0) >= d^2/dt^2 sum_{k<=j} lambda_k at 0. The bound is exact (tight) with the default
    correctors, and it remains valid for any trial. The prefix is a union of complete groups, because the gap
    check below guarantees a gap above lambda_j.

    Raises PencilError if K or M is not symmetric or M is not positive definite, ClusterRangeError for a bad
    ``cluster``, and ClusterSplitError if the cluster splits a (near-)degenerate group (all subclass ValueError).
    The split check fails when g_out < max(gap_ratio_min |lam_bar|, 10 g_in), where g_in is the spread inside the
    cluster and g_out the distance to the nearest eigenvalue outside it.
    """
    K, M, sparse = _prepare(K, M)
    n = K.shape[0]
    i, j = _check_range(cluster, n)
    _validate(K, M, sparse)
    dK, dM = _as_op(dK, sparse, n), _as_op(dM, sparse, n)
    have_d2 = d2K is not None or d2M is not None
    d2K, d2M = _as_op(d2K, sparse, n), _as_op(d2M, sparse, n)

    w, V = _eigs(K, M, j + 2, sparse)  # indices 0..j+1: lambda_{j+1} is needed for the gap check
    lam = w[i:j + 1]
    lam_bar, g_in, g_out, thresh, ratio = _gap(w, i, j, gap_ratio_min)
    if g_out < thresh:
        raise _split_error(cluster, g_in, g_out, ratio)

    U = V[:, i:j + 1]
    B = np.asarray(dK @ U) - lam_bar * np.asarray(dM @ U)
    A = U.T @ B
    A = (A + A.T) / 2
    rates = np.linalg.eigvalsh(A)

    trace_hessian = kyfan = None
    if have_d2:
        Z = _deflated_solve(K, M, lam_bar, U, _proj(M, U, B), sparse)
        UdMU = U.T @ np.asarray(dM @ U)
        H2 = U.T @ (np.asarray(d2K @ U) - lam_bar * np.asarray(d2M @ U))
        trace_hessian = float(np.trace(H2) - 2 * np.trace(A @ UdMU) - 2 * np.sum(B * Z))
        kyfan = _kyfan_upper(K, M, dK, dM, d2K, d2M, w[:j + 1], V[:, :j + 1], trial, sparse)

    return ClusterSensitivity(eigvals=lam.copy(), basis=U.copy(), A=A, rates=rates,
                              trace_hessian=trace_hessian, kyfan_upper=kyfan,
                              cluster=(i, j), lam_bar=lam_bar, gap_ratio=ratio)


def _kyfan_upper(K, M, dK, dM, d2K, d2M, lamW, W, trial, sparse) -> float:
    """Second derivative at t = 0 of g(t) = tr[N(t)^{-1} G(t)], with G = W_t^T K W_t, N = W_t^T M W_t,
    W_t = W + t Psi. With G0 = W^T K W and N0 = W^T M W = I (W is M-orthonormal),
        g'' = tr(G'') - 2 tr(N' G') + tr((2 N'N' - N'') G0),
    where (Psi enters linearly, so W'' = 0)
        G'  = Psi^T K W + W^T dK W + W^T K Psi,       N'  = Psi^T M W + W^T dM W + W^T M Psi,
        G'' = W^T d2K W + 2(Psi^T dK W + W^T dK Psi) + 2 Psi^T K Psi,
        N'' = W^T d2M W + 2(Psi^T dM W + W^T dM Psi) + 2 Psi^T M Psi.
    This is the exact g''(0). It is invariant under K -> K - sigma M (that only shifts g by a constant)."""
    n, p = W.shape
    if trial is None:
        Psi = np.empty((n, p))
        for a in range(p):  # exact first-order corrector of eigenvector a, off the prefix span
            b = np.asarray(dK @ W[:, a]) - lamW[a] * np.asarray(dM @ W[:, a])
            Psi[:, a] = -_deflated_solve(K, M, float(lamW[a]), W, _proj(M, W, b[:, None]), sparse)[:, 0]
    else:
        Psi = np.asarray(trial, dtype=float)
        if Psi.shape != (n, p):
            raise ValueError(f"trial must have shape {(n, p)} (all eigenvectors 0..j of the prefix)")
    KW, MW = np.asarray(K @ W), np.asarray(M @ W)
    KP, MP = np.asarray(K @ Psi), np.asarray(M @ Psi)
    dKW, dMW = np.asarray(dK @ W), np.asarray(dM @ W)
    G0 = W.T @ KW
    t1 = Psi.T @ KW
    G1 = t1 + t1.T + W.T @ dKW
    t2 = Psi.T @ MW
    N1 = t2 + t2.T + W.T @ dMW
    c = Psi.T @ dKW
    G2 = W.T @ np.asarray(d2K @ W) + 2 * (c + c.T) + 2 * Psi.T @ KP
    c = Psi.T @ dMW
    N2 = W.T @ np.asarray(d2M @ W) + 2 * (c + c.T) + 2 * Psi.T @ MP
    return float(np.trace(G2) - 2 * np.trace(N1 @ G1) + np.trace((2 * N1 @ N1 - N2) @ G0))


def directional(K, M, dK, dM, cluster, v, d2K=None, d2M=None, **kw) -> ClusterSensitivity:
    """Like cluster_sensitivity, but dK, dM, d2K and d2M are callables of a direction vector v. For example,
    d2K(v) is the second derivative of K along v."""
    return cluster_sensitivity(K, M, dK(v), dM(v), cluster,
                               d2K=None if d2K is None else d2K(v),
                               d2M=None if d2M is None else d2M(v), **kw)


def from_callables(Kfun, Mfun, p0, v, cluster, h=1e-4, **kw) -> ClusterSensitivity:
    """Sensitivity along p0 + t v, with dK, dM, d2K and d2M taken from central finite differences of K(p) and M(p)
    (step h). This emits a FiniteDifferenceWarning (a UserWarning): finite-difference derivatives carry O(h^2)
    truncation error and O(eps/h^2) round-off in the second derivatives."""
    warnings.warn("finite-difference derivatives used in from_callables (O(h^2) truncation error)",
                  FiniteDifferenceWarning, stacklevel=2)
    p0 = np.asarray(p0, dtype=float)
    v = np.asarray(v, dtype=float)
    K0, Kp, Km = Kfun(p0), Kfun(p0 + h * v), Kfun(p0 - h * v)
    M0, Mp, Mm = Mfun(p0), Mfun(p0 + h * v), Mfun(p0 - h * v)
    dK = (Kp - Km) / (2 * h)
    dM = (Mp - Mm) / (2 * h)
    d2K = (Kp - 2 * K0 + Km) / h**2
    d2M = (Mp - 2 * M0 + Mm) / h**2
    return cluster_sensitivity(K0, M0, dK, dM, cluster, d2K=d2K, d2M=d2M, **kw)


# ----------------------------------------------------------------------------------------------- multi-parameter
@dataclass(frozen=True)
class SplittingMatrices:
    eigvals: np.ndarray      # cluster eigenvalues, ascending
    basis: np.ndarray        # n x m, M-orthonormal
    lam_bar: float
    mats: np.ndarray         # (d, m, m): A_l = U^T (dK_l - lam_bar dM_l) U, symmetrised
    gap_ratio: float
    label: str = "NUMERICAL"
    caveat: str = DISCRETE_CAVEAT


def splitting_matrices(K, M, dKs, dMs, cluster, *, gap_ratio_min=1e-3, eig=None, check_gap: bool = True,
                       validate: bool = True) -> SplittingMatrices:
    """Splitting matrices A_l = U^T (dK_l - lam_bar dM_l) U of the cluster (i, j), one per parameter l.

    ``dKs`` (and ``dMs``; None means M does not depend on p) are length-d sequences of partial derivatives, dense or
    sparse. A(v) = sum_l v_l mats[l] equals cluster_sensitivity(..., dK=sum v_l dK_l, dM=sum v_l dM_l).A.
    ``eig=(w, V)`` reuses eigenpairs from smallest_eigenpairs (at least min(j + 2, n) of them); it skips the
    eigensolve but not the validation. ``validate=False`` skips the (K, M) validation; ``check_gap=False`` skips
    the core gap rule (gap_ratio is still reported). The rates eig(A(v)) are reduced-model slopes at a near-tie
    (see the module docstring)."""
    K, M, sparse = _prepare(K, M)
    n = K.shape[0]
    i, j = _check_range(cluster, n)
    if validate:
        _validate(K, M, sparse)
    need = min(j + 2, n)
    if eig is None:
        w, V = _eigs(K, M, need, sparse)
    else:
        w, V = np.asarray(eig[0], dtype=float), np.asarray(eig[1], dtype=float)
        if len(w) < need or V.shape != (n, len(w)):
            raise ValueError(f"eig must hold at least {need} eigenpairs of the n = {n} pencil")
    lam_bar, g_in, g_out, thresh, ratio = _gap(w, i, j, gap_ratio_min)
    if check_gap and g_out < thresh:
        raise _split_error(cluster, g_in, g_out, ratio)
    U = V[:, i:j + 1]
    dKs = list(dKs)
    dMs = [None] * len(dKs) if dMs is None else list(dMs)
    if len(dMs) != len(dKs):
        raise ValueError("dKs and dMs must have the same length")
    mats = np.empty((len(dKs), j - i + 1, j - i + 1))
    for l, (dK, dM) in enumerate(zip(dKs, dMs)):
        B = np.asarray(_as_op(dK, sparse, n) @ U)
        if dM is not None:
            B = B - lam_bar * np.asarray(_as_op(dM, sparse, n) @ U)
        A = U.T @ B
        mats[l] = (A + A.T) / 2
    return SplittingMatrices(eigvals=w[i:j + 1].copy(), basis=U.copy(), lam_bar=lam_bar, mats=mats,
                             gap_ratio=ratio)


def ritz_upper_bounds(K, M, W) -> np.ndarray:
    """Ritz values theta_0 <= ... <= theta_{p-1} of the pencil compressed to span W (W: n x p, full column rank).

    Poincare separation: lambda_i(K, M) <= theta_i for i < p, hence sum_{i<p} lambda_i <= tr[(W^T M W)^{-1} W^T K W]
    (the Ky Fan trial bound). Valid for the DISCRETE pencil. Raises PencilError if cond(W^T M W) > 1e12."""
    W = np.asarray(W, dtype=float)
    if W.ndim != 2 or W.shape[0] != K.shape[0]:
        raise ValueError(f"W must have shape (n, p) with n = {K.shape[0]}")
    KW, MW = np.asarray(K @ W), np.asarray(M @ W)
    G = W.T @ MW
    G = (G + G.T) / 2
    H = W.T @ KW
    H = (H + H.T) / 2
    c = np.linalg.cond(G)
    if not np.isfinite(c) or c > 1e12:
        raise PencilError(f"W must have full column rank: cond(W^T M W) = {c:.3g} > 1e12")
    return sla.eigh(H, G, eigvals_only=True)
