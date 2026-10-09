import ast, pathlib
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
    assert bolzakit.__version__ == "0.2.0"
    assert issubclass(CertificationFailed, BolzakitError)
