"""Fundamental-frequency maximisation benchmark (Task 9'): the model, and a quick run of every method."""
import pathlib, sys
import numpy as np
import pytest

pytest.importorskip("skfem")
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "benchmarks"))
import frequency_maximisation as fm  # noqa: E402
from bolzakit.sensitivity import smallest_eigenpairs  # noqa: E402


def test_post_supported_plate_spectrum():
    Ks, Ms, areas = fm.build()
    w, _ = smallest_eigenpairs(sum(Ks), sum(Ms), 5)                  # uniform design h = 1
    assert np.allclose(w / np.pi ** 4, [1.58608, 1.83370, 1.83370, 4.37522, 6.76920], atol=5e-5)
    assert abs(w[2] - w[1]) <= 1e-9 * w[1]                            # the (1,2)/(2,1)-type pair, tied by symmetry
    assert abs(areas.sum() - 1.0) < 1e-12 and len(Ks) == 25


@pytest.mark.slow
def test_benchmark_quick_runs_and_reports_honestly():
    rows = fm.main(quick=True)
    by = {(r["seed"], r["method"]): r for r in rows}
    assert {r["seed"] for r in rows} == {0, 1}
    for s in (0, 1):
        assert by[(s, "ours")]["decreases"] == 0                       # monotone (also asserted inside main)
        assert by[(s, "naive fixed-step")]["decreases"] >= 10          # the naive baseline chatters
        assert by[(s, "ours")]["final"] > 1.9                          # lambda_0 rises from 1.586 / 1.677 pi^4
    assert by[(0, "ours")]["mult"] == 3                                # symmetric start: a triple eigenvalue
