"""Every ```python block in README.md must run.

Blocks may start with requirement markers:
  # requires: pari     skipped unless cypari2 is importable
  # requires: network  skipped unless the Strohmaier–Uski file is already in the cache with the right sha256, or
                       pytest is run with --runnetwork (so plain `pytest` never downloads anything)
  # requires: slow     skipped unless pytest is run with --runslow
"""
import importlib.util
import pathlib
import re

import pytest

README = pathlib.Path(__file__).resolve().parents[1] / "README.md"
BLOCKS = re.findall(r"```python\n(.*?)```", README.read_text(encoding="utf-8"), re.S)


def _su_cached() -> bool:
    from conftest import su_cached
    return su_cached()


def test_readme_has_blocks():
    assert len(BLOCKS) >= 8


def test_readme_blocks_assert_values():
    # the README doubles as a doctest: key printed values are asserted inside the blocks
    assert sum(b.count("assert ") for b in BLOCKS) >= 4


@pytest.mark.parametrize("code", BLOCKS, ids=[f"block{i}" for i in range(len(BLOCKS))])
def test_readme_block(code, request):
    if "# requires: pari" in code and importlib.util.find_spec("cypari2") is None:
        pytest.skip("needs cypari2 (pip install 'bolzakit[pari]')")
    if "# requires: network" in code and not (_su_cached() or request.config.getoption("--runnetwork")):
        pytest.skip("needs --runnetwork (Strohmaier–Uski file not cached)")
    if "# requires: slow" in code and not request.config.getoption("--runslow"):
        pytest.skip("slow (run with --runslow)")
    exec(compile(code, "README", "exec"), {})
