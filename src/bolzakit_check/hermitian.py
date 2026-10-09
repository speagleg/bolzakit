"""Ball linear algebra for the checker: symmetric LDL^T positive-definiteness, Bareiss determinant, Ky Fan trace."""
from __future__ import annotations

from flint import arb, arb_mat

from ._cert import CheckError


def check_symmetric(A):
    """A itself, after checking that A_ij and A_ji overlap (else it cannot enclose a symmetric Gram matrix).

    Each stored entry encloses the true symmetric value, so A as stored (lower triangle, or the full matrix in a
    solve) is a valid enclosure; no symmetrising union is needed."""
    n = len(A)
    for i in range(n):
        for j in range(i):
            if not A[i][j].overlaps(A[j][i]):
                raise CheckError(f"stored matrix is not symmetric at ({i},{j})")
    return A


def ldl_pd(A) -> bool:
    """True iff every symmetric matrix in the ball matrix A is positive definite, certified by an LDL^T
    factorisation (no pivoting) whose pivots d_j are all certainly > 0.  Uses the lower triangle of A."""
    n = len(A)
    L = [[arb(0)] * n for _ in range(n)]
    d = [arb(0)] * n
    for j in range(n):
        dj = A[j][j] - sum((L[j][k] * L[j][k] * d[k] for k in range(j)), arb(0))
        if not (dj.is_finite() and dj > 0):
            return False
        d[j] = dj
        for i in range(j + 1, n):
            L[i][j] = (A[i][j] - sum((L[i][k] * L[j][k] * d[k] for k in range(j)), arb(0))) / dj
    return True


def is_pd(A) -> bool:
    """Certified positive definiteness: own LDL^T pivots, or else Sylvester's criterion (leading principal minors
    certainly > 0, flint determinants).  Each criterion alone is a proof; the second only guards against the
    different overestimation of ball arithmetic in the two algorithms."""
    if ldl_pd(A):
        return True
    n = len(A)
    for k in range(1, n + 1):
        d = arb_mat([[A[i][j] for j in range(k)] for i in range(k)]).det()
        if not (d.is_finite() and d > 0):
            return False
    return True


def bareiss_det(A) -> arb:
    """Determinant of a square ball matrix by fraction-free (Bareiss) elimination with row pivoting on the entry
    of largest certain magnitude.  Returns nan if no pivot certainly excludes 0 (determinant undecided)."""
    M = [list(r) for r in A]
    n = len(M)
    sign, prev = 1, arb(1)
    for k in range(n - 1):
        best, bmag = None, None
        for r in range(k, n):
            x = M[r][k]
            if x.is_finite() and not x.contains(0):
                mag = abs(x).lower()
                if bmag is None or mag > bmag:
                    best, bmag = r, mag
        if best is None:
            return arb("nan")
        if best != k:
            M[k], M[best] = M[best], M[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                M[i][j] = (M[i][j] * M[k][k] - M[i][k] * M[k][j]) / prev
        prev = M[k][k]
    return sign * M[n - 1][n - 1]


def trace_minv_k(K, M) -> arb:
    """tr(M^-1 K) for ball matrices (nan if M is not provably invertible)."""
    try:
        return arb_mat(M).solve(arb_mat(K)).trace()
    except ZeroDivisionError:
        return arb("nan")
