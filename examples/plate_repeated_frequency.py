"""Repeated frequency of a stretched Dirichlet membrane (NUMERICAL, discrete problem).

The unit square is stretched to a x 1 with a = 1 + p. Pulled back to the unit square (xi-scaling):
    K(p) = int (1+p)^-1 u_xi v_xi + (1+p) u_eta v_eta,    M(p) = int (1+p) u v.
Exact eigenvalues: lambda_mn(a) = pi^2 (m^2/a^2 + n^2). At p = 0 the modes (1,2) and (2,1) share lambda = 5 pi^2.

Why the cluster is treated as a unit: with several parameters, the repeated eigenvalue splits into branches whose
directional derivatives depend on the direction, so the individual eigenvalues are not differentiable at the
degenerate point (only directionally, with rates that are the eigenvalues of A). The SUM of the cluster
eigenvalues is smooth (analytic), so its gradient and Hessian are well defined. Here, with one parameter, the
branches are still differentiable: rates -2 pi^2 and -8 pi^2.
"""
import numpy as np
from skfem import MeshQuad, ElementQuad2, Basis, BilinearForm

from bolzakit.sensitivity import cluster_sensitivity
from bolzakit.sensitivity.adapters.skfem import from_forms


@BilinearForm
def kform(u, v, w):
    p = w["p"]
    return (1 / (1 + p)) * u.grad[0] * v.grad[0] + (1 + p) * u.grad[1] * v.grad[1]


@BilinearForm
def mform(u, v, w):
    return (1 + w["p"]) * u * v


def main():
    basis = Basis(MeshQuad().refined(5), ElementQuad2())
    K, M, dK, dM, d2K, d2M = from_forms(basis, kform, mform, 0.0, dirichlet_dofs=basis.get_dofs())
    r = cluster_sensitivity(K, M, dK, dM, (1, 2), d2K=d2K, d2M=d2M)
    pi2 = np.pi ** 2
    rows = [("rate (1,2)", r.rates[1], -2 * pi2), ("rate (2,1)", r.rates[0], -8 * pi2),
            ("cluster-sum second derivative", r.trace_hessian, 30 * pi2),
            ("Ky Fan bound, prefix lam11+lam12+lam21", r.kyfan_upper, 36 * pi2)]
    worst = 0.0
    for name, val, ex in rows:
        err = abs(val - ex) / abs(ex)
        worst = max(worst, err)
        print(f"{name:42s} {val:12.5f}  exact {ex:12.5f}  rel.err {err:.2e}")
    print(f"label: {r.label}; worst relative error {worst:.2e}")
    assert worst < 5e-3


if __name__ == "__main__":
    main()
