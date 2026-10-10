import ast, importlib.util, os, pathlib, subprocess, sys
import pytest
SRC = pathlib.Path(__file__).resolve().parents[1] / "src"

def _under(m, prefixes):
    """m is one of the packages in prefixes or a submodule of one (dot boundary: bolzakit.certifiedX is not under
    bolzakit.certified)."""
    return any(m == p or m.startswith(p + ".") for p in prefixes)

def _imports(pkgdir):
    out = set()
    files = pkgdir.rglob("*.py") if pkgdir.is_dir() else [pkgdir]
    for f in files:
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
           if _under(m, ("bolzakit",)) and not _under(m, ("bolzakit.sensitivity",))]
    assert bad == []

def test_checker_never_imports_bolzakit():
    bad = [(f, m) for f, m in _imports(SRC / "bolzakit_check")
           if m == "bolzakit" or m.startswith("bolzakit.")]
    assert bad == []

def test_certified_layering():
    allowed = ("bolzakit.certified", "bolzakit.exact", "bolzakit.geometry", "bolzakit.provenance", "bolzakit.errors")
    bad = [(f, m) for f, m in _imports(SRC / "bolzakit" / "certified")
           if _under(m, ("bolzakit",)) and not _under(m, allowed)]
    assert bad == []

def test_v01_modules_do_not_import_v02():
    # v0.1 subpackages and the top-level modules (__init__, errors, provenance)
    v01 = [p for p in (SRC / "bolzakit").iterdir()
           if (p.is_dir() and p.name not in ("certified", "sensitivity", "__pycache__", "data")) or p.suffix == ".py"]
    assert SRC / "bolzakit" / "__init__.py" in v01
    bad = []
    for d in v01:
        bad += [(f, m) for f, m in _imports(d) if _under(m, ("bolzakit.certified", "bolzakit.sensitivity"))]
    assert bad == []

def test_version_and_error():
    import bolzakit
    from bolzakit.errors import CertificationFailed, BolzakitError
    assert bolzakit.__version__ == "0.3.0"
    assert issubclass(CertificationFailed, BolzakitError)


# ------------------------------------------------------------------ sub-project F (spec 3.1)
SENS = SRC / "bolzakit" / "sensitivity"
NO_FRAMEWORK = [SENS / "core.py", SENS / "errors.py", SENS / "optimize.py", SENS / "trim.py", SENS / "__init__.py",
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


def _new_layer_refs(f):
    """Every way file f can reach bolzakit.sensitivity.optimize / .trim: module imports, names imported from
    bolzakit.sensitivity or bolzakit (`from bolzakit.sensitivity import optimize`, `from bolzakit import sensitivity`
    followed by `sensitivity.trim`), and attribute chains `<...>.sensitivity.optimize` / `.trim`."""
    new = ("bolzakit.sensitivity.optimize", "bolzakit.sensitivity.trim")
    tree = ast.parse(f.read_text(encoding="utf-8"))
    out = [m for _, m in _imports_file(f) if any(_is(m, n) for n in new)]
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and not node.level and node.module in ("bolzakit", "bolzakit.sensitivity"):
            out += [f"{node.module}.{a.name}" for a in node.names
                    if any(_is(f"{node.module}.{a.name}", n) for n in new)]
        if isinstance(node, ast.Attribute) and node.attr in ("optimize", "trim"):
            v = node.value
            if (isinstance(v, ast.Name) and v.id == "sensitivity") or (
                    isinstance(v, ast.Attribute) and v.attr == "sensitivity"):
                out.append(f"<attribute>.sensitivity.{node.attr}")
    return out


def test_nothing_outside_sensitivity_imports_the_new_layers():
    bad = []
    for f in list((SRC / "bolzakit").rglob("*.py")) + list((SRC / "bolzakit_check").rglob("*.py")):
        if SENS in f.parents:
            continue
        bad += [(f, m) for m in _new_layer_refs(f)]
    assert bad == []


def test_new_layer_check_sees_package_level_imports(tmp_path):
    cases = {"a.py": "from bolzakit.sensitivity import optimize\n", "b.py": "from bolzakit import sensitivity\n"
             "x = sensitivity.trim\n", "c.py": "import bolzakit.sensitivity.trim\n",
             "d.py": "import bolzakit\ny = bolzakit.sensitivity.optimize\n", "ok.py": "from bolzakit import exact\n"}
    found = {}
    for name, code in cases.items():
        (tmp_path / name).write_text(code, encoding="utf-8")
        found[name] = bool(_new_layer_refs(tmp_path / name))
    assert found == {"a.py": True, "b.py": True, "c.py": True, "d.py": True, "ok.py": False}


def _core_imports():
    """The core modules that exist so far (optimize.py arrives in Task 6)."""
    mods = ["bolzakit.sensitivity"] + [f"bolzakit.sensitivity.{m}" for m in ("optimize", "trim")
                                       if (SENS / f"{m}.py").exists()]
    return ", ".join(mods)


def _run(code):
    return subprocess.run([sys.executable, "-c", code], capture_output=True, text=True,
                          env={**os.environ, "PYTHONPATH": str(SRC)})


def test_importable_without_torch_and_jax():
    r = _run('import sys; sys.modules["torch"] = None; sys.modules["jax"] = None\n'
             f'import {_core_imports()}\nprint("ok")')
    assert r.returncode == 0 and "ok" in r.stdout, r.stderr


@pytest.mark.skipif(importlib.util.find_spec("torch") is None and importlib.util.find_spec("jax") is None,
                    reason="neither torch nor jax installed (optional check; no CI job installs them)")
def test_core_import_leaves_frameworks_unloaded_with_extras_installed():
    r = _run(f'import sys\nimport {_core_imports()}\n'
             'print(sorted(m for m in ("torch", "jax") if m in sys.modules))')
    assert r.returncode == 0 and r.stdout.strip() == "[]", (r.stdout, r.stderr)


def test_trim_imports_only_core_numpy_scipy():
    ok = ("bolzakit.sensitivity", "bolzakit.sensitivity.core", "numpy", "scipy", "__future__", "dataclasses")
    bad = [m for m in _mods(SENS / "trim.py") if not any(_is(m, o) for o in ok)]
    assert bad == []
