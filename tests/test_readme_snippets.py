"""Every ```python block in README.md must run.

Blocks may start with requirement markers:
  # requires: pari     skipped unless cypari2 is importable
  # requires: network  skipped unless the Strohmaier–Uski file is already in the cache with the right sha256, or
                       pytest is run with --runnetwork (so plain `pytest` never downloads anything)
  # requires: slow     skipped unless pytest is run with --runslow
  # requires: cert     skipped unless python-flint and sympy are importable (bolzakit[cert])
Blocks run with the repository root as the working directory (the README's snippets assume it).

The three audience sections ("For software engineers", "For physicists", "For engineers") must each have at least one
python block of at most 6 code lines (requirement markers excluded); those blocks run here like every other block.
"""
import importlib.util
import pathlib
import re

import pytest

README = pathlib.Path(__file__).resolve().parents[1] / "README.md"
TEXT = README.read_text(encoding="utf-8")
BLOCKS = re.findall(r"```python\n(.*?)```", TEXT, re.S)
AUDIENCES = ["For software engineers", "For physicists", "For engineers"]


def _section(title):
    m = re.search(rf"^## {re.escape(title)}\n(.*?)(?=^## )", TEXT, re.S | re.M)
    assert m, title
    return m.group(1)


def _su_cached() -> bool:
    from conftest import su_cached
    return su_cached()


def test_readme_has_blocks():
    assert len(BLOCKS) >= 8


def test_readme_blocks_assert_values():
    # the README doubles as a doctest: key printed values are asserted inside the blocks
    assert sum(b.count("assert ") for b in BLOCKS) >= 4


@pytest.mark.parametrize("title", AUDIENCES)
def test_audience_section_has_short_snippet(title):
    blocks = re.findall(r"```python\n(.*?)```", _section(title), re.S)
    assert blocks, title
    for b in blocks:
        assert b in BLOCKS                                   # so test_readme_block executes it
        code = [ln for ln in b.splitlines() if ln.strip() and not ln.startswith("# requires:")]
        assert len(code) <= 6, (title, len(code))


@pytest.mark.parametrize("code", BLOCKS, ids=[f"block{i}" for i in range(len(BLOCKS))])
def test_readme_block(code, request, monkeypatch):
    if "# requires: pari" in code and importlib.util.find_spec("cypari2") is None:
        pytest.skip("needs cypari2 (pip install 'bolzakit[pari]')")
    if "# requires: network" in code and not (_su_cached() or request.config.getoption("--runnetwork")):
        pytest.skip("needs --runnetwork (Strohmaier–Uski file not cached)")
    if "# requires: slow" in code and not request.config.getoption("--runslow"):
        pytest.skip("slow (run with --runslow)")
    if "# requires: cert" in code and not all(importlib.util.find_spec(m) for m in ("flint", "sympy")):
        pytest.skip("needs python-flint and sympy (pip install 'bolzakit[cert]')")
    monkeypatch.chdir(README.parent)
    exec(compile(code, "README", "exec"), {})
