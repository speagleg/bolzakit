"""Independent re-derivation of a Lehmann-Goerisch lower bound (kind "lg_lower").

Written from the theorem statement (Liu-Plum arXiv:2512.23182 Thm 4.1/4.3 with A1-A5, as quoted in LB_PLAN.md §1(c)
and review_lb_1 rows 1a-1f; Behnke-Mertins-Plum-Wieners 2000 Thm 2.1 in shifted form), not from the producer:
  A0 = (M(v_i, v_j)), A1 = (N(v_i, v_j)), A2 = (b_G(w_i, w_j)),  A = A0 - rho A1,  B = A0 - 2 rho A1 + rho^2 A2.
  If rho > 0, rho <= lambda~_{m+1} with m >= n, B is positive definite and nu_1 <= .. <= nu_n are the eigenvalues of
  A z = nu B z, then lambda~_{m+1-k} >= rho - rho/(1 - nu_k) for every k with nu_k < 0.
Instantiation (shift s = 1): lambda~ = lambda + s, lambda~_1 = s (constants), m = n + 1, so k = n bounds
lambda~_2 = lambda_1 + s, and rho <= lambda~_{n+2} = lambda_{n+1} + s is the gap hypothesis lambda_{n+1} >= rho - s
(an input of the certificate: the checker reports it as a condition).

The eigenvalue step, done directly at the claimed value: for a rational L put nu* = 1 - rho/(rho - s - L), so that
rho - rho/(1 - nu*) - s = L.  If nu* < 0 and nu* B - A is positive definite, then z^T A z < nu* z^T B z for all
z != 0, so every nu_k < nu* < 0 and, the map nu -> rho - rho/(1 - nu) being decreasing,
lambda_1 >= rho - rho/(1 - nu_n) - s > L.  No eigenvalue solver and no stored nubar is used.

L is the lower end of the stored enclosure data["lam1_lower_enclosure"] (so the float the producer attests from it
is covered), and the claim's rhs must be <= that L.  Trusted: the stored Gram balls A0, A1, A2 (element integration
and the H(div) conformity of the fluxes, which the checker cannot see), and the gap input.
"""
from __future__ import annotations

from fractions import Fraction

from . import hermitian
from ._cert import CheckError, arb_q, frac, matrix, precision

PREC = 256


def _ball_lower(s) -> Fraction:
    """Exact lower end of a stored ball string '[mid +/- rad]' (or an exact decimal)."""
    if not isinstance(s, str):
        raise CheckError("lam1_lower_enclosure is not a ball string")
    t = s.strip()
    try:
        if t.startswith("[") and t.endswith("]"):
            mid, rad = t[1:-1].split("+/-")
            return Fraction(mid.strip()) - Fraction(rad.strip())
        return Fraction(t)
    except (ValueError, ZeroDivisionError) as e:
        raise CheckError(f"unparsable enclosure {s!r}") from e


def verify(raw: dict):
    """([(lhs, op, rhs, truth)], conditional note) for an lg_lower certificate."""
    data, inputs = raw["data"], raw["inputs"]
    k = data["k"]
    if isinstance(k, bool) or not isinstance(k, int) or k < 1:
        raise CheckError("bad lg_lower k")
    if data["shift"] != 1 or isinstance(data["shift"], bool):
        raise CheckError("lg_lower shift must be 1 (M = q + m, N = m)")
    rho = frac(data["rho"])
    if rho <= 0:
        raise CheckError("rho must be > 0 (Liu-Plum A5)")
    gap_s = inputs["gap"]["lower"]
    gap = frac(gap_s)
    if gap + 1 < rho:
        raise CheckError(f"Lehmann-Goerisch hypothesis rho <= lambda_(k+1) + 1 fails: rho = {rho}, gap = {gap_s}")
    deps = {d.get("id") for d in raw["dependencies"]}
    if not {"D15", "gap"} <= deps:
        raise CheckError("lg_lower certificate lacks the D15 / gap dependencies")
    rhs = raw["claims"][0]["rhs"]
    note = f"lambda_{k + 1} >= {gap_s}"
    if rhs == "nan":                                   # the producer certified nothing: nothing to re-derive
        return [("lambda_1", ">=", rhs, False)], note
    L = _ball_lower(data["lam1_lower_enclosure"])
    ok = False
    with precision(PREC):
        A0 = hermitian.check_symmetric(matrix(data["A0"], k))
        A1 = hermitian.check_symmetric(matrix(data["A1"], k))
        A2 = hermitian.check_symmetric(matrix(data["A2"], k))
        r = arb_q(rho)
        A = [[A0[i][j] - r * A1[i][j] for j in range(k)] for i in range(k)]
        B = [[A0[i][j] - 2 * r * A1[i][j] + r * r * A2[i][j] for j in range(k)] for i in range(k)]
        den = rho - 1 - L
        if frac(rhs) <= L and den > 0:
            nu = 1 - rho / den                         # exact rational
            if nu < 0 and hermitian.is_pd(B):
                nuA = arb_q(nu)
                ok = hermitian.is_pd([[nuA * B[i][j] - A[i][j] for j in range(k)] for i in range(k)])
    return [("lambda_1", ">=", rhs, ok)], note
