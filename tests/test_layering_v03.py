import ast, importlib.util, os, pathlib, subprocess, sys
import pytest
SRC = pathlib.Path(__file__).resolve().parents[1] / "src"

def _imports(pkgdir):
    out = set()
    for f in pkgdir.rglob("*.py"):
        tree = ast.parse(f.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                out |= {(f, a.name) for a in node.names}
            elif isinstance(node, ast.ImportFrom):
                mod = node.module or ""
                if node.level:  # relative import -> resolve against package
                    rel = f.relative_to(SRC).with_suffix("").parts[:-node.level]
                    mod = ".".join(rel + ((mod,) if mod else ()))
                out.add((f, mod))
    return out

def test_sensitivity_is_standalone():
    bad = [(f, m) for f, m in _imports(SRC / "bolzakit" / "sensitivity")
           if m.startswith("bolzakit") and not m.startswith("bolzakit.sensitivity")]
    assert bad == []

def test_checker_never_imports_bolzakit():
    bad = [(f, m) for f, m in _imports(SRC / "bolzakit_check")
           if m == "bolzakit" or m.startswith("bolzakit.")]
    assert bad == []

def test_certified_layering():
    allowed = ("bolzakit.certified", "bolzakit.exact", "bolzakit.geometry", "bolzakit.provenance", "bolzakit.errors")
    bad = [(f, m) for f, m in _imports(SRC / "bolzakit" / "certified")
           if m.startswith("bolzakit") and not m.startswith(allowed)]
    assert bad == []

def test_v01_modules_do_not_import_v02():
    v01 = [p for p in (SRC / "bolzakit").iterdir() if p.is_dir() and p.name not in ("certified", "sensitivity", "__pycache__", "data")]
    bad = []
    for d in v01:
        bad += [(f, m) for f, m in _imports(d) if m.startswith(("bolzakit.certified", "bolzakit.sensitivity"))]
    assert bad == []

def test_version_and_error():
    import bolzakit
    from bolzakit.errors import CertificationFailed, BolzakitError
    assert bolzakit.__version__ == "0.3.0"
    assert issubclass(CertificationFailed, BolzakitError)


# ------------------------------------------------------------------ sub-project F (spec 3.1)
SENS = SRC / "bolzakit" / "sensitivity"
NO_FRAMEWORK = [SENS / "core.py", SENS / "errors.py", SENS / "optimize.py", SENS / "__init__.py",
                SENS / "adapters" / "_util.py"]


def _mods(path):
    return {m for f, m in _imports_file(path)}


def _imports_file(f):
    tree = ast.parse(f.read_text(encoding="utf-8"))
    out = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            out |= {(f, a.name) for a in node.names}
        elif isinstance(node, ast.ImportFrom):
            mod = node.module or ""
            if node.level:
                rel = f.relative_to(SRC).with_suffix("").parts[:-node.level]
                mod = ".".join(rel + ((mod,) if mod else ()))
            out.add((f, mod))
    return out


def _is(m, root):
    return m == root or m.startswith(root + ".")


def test_core_errors_optimize_import_no_framework():
    for f in NO_FRAMEWORK:
        if not f.exists():            # optimize.py arrives in Task 6
            continue
        bad = [m for m in _mods(f) if _is(m, "torch") or _is(m, "jax") or _is(m, "jaxlib")]
        assert bad == [], (f.name, bad)


def test_errors_module_is_stdlib_only():
    if not (SENS / "errors.py").exists():   # errors.py arrives in Task 1
        pytest.skip("errors.py not written yet")
    assert _mods(SENS / "errors.py") <= {"__future__"}


def test_nothing_outside_sensitivity_imports_the_new_layers():
    new = ("bolzakit.sensitivity.optimize",)
    bad = []
    for f in (SRC / "bolzakit").rglob("*.py"):
        if SENS in f.parents:
            continue
        bad += [(f, m) for _, m in _imports_file(f) if m.startswith(new)]
    for f in (SRC / "bolzakit_check").rglob("*.py"):
        bad += [(f, m) for _, m in _imports_file(f) if m.startswith(new)]
    assert bad == []


def _core_imports():
    """The core modules that exist so far (optimize.py arrives in Task 6)."""
    mods = ["bolzakit.sensitivity"] + (["bolzakit.sensitivity.optimize"] if (SENS / "optimize.py").exists() else [])
    return ", ".join(mods)


def _run(code):
    return subprocess.run([sys.executable, "-c", code], capture_output=True, text=True,
                          env={**os.environ, "PYTHONPATH": str(SRC)})


def test_importable_without_torch_and_jax():
    r = _run('import sys; sys.modules["torch"] = None; sys.modules["jax"] = None\n'
             f'import {_core_imports()}\nprint("ok")')
    assert r.returncode == 0 and "ok" in r.stdout, r.stderr


@pytest.mark.skipif(importlib.util.find_spec("torch") is None and importlib.util.find_spec("jax") is None,
                    reason="neither extra installed (the CI torch/jax jobs run this)")
def test_core_import_leaves_frameworks_unloaded_with_extras_installed():
    r = _run(f'import sys\nimport {_core_imports()}\n'
             'print(sorted(m for m in ("torch", "jax") if m in sys.modules))')
    assert r.returncode == 0 and r.stdout.strip() == "[]", (r.stdout, r.stderr)
