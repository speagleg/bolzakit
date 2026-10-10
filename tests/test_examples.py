"""The 0.3.0 examples run and pass their own assertions. The v0.2 examples (membrane, bolza_*) run
in tests/test_docs_v02.py::test_examples_run; the frequency benchmark in tests/test_frequency_benchmark.py."""
import importlib, pathlib, sys
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "examples"))


@pytest.mark.slow
def test_gyro_trim_quick():
    pytest.importorskip("skfem")
    out = importlib.import_module("gyro_trim").main(quick=True)
    assert abs(out["a"]["rate"] - 0.3) < 5e-4
    assert [x["res"].status for x in out["b"]["runs"]] == ["converged", "converged"]
