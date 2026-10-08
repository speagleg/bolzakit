import bolzakit as bk
from bolzakit.provenance import registry, LABELS


def test_top_level_names():
    for name in ("group", "automorphisms", "lengths", "systole", "spectral", "green", "vortex", "arith", "geometry"):
        assert hasattr(bk, name), name


def test_spectral_reexports():
    assert callable(bk.spectral.tau_bounds)
    assert callable(bk.spectral.selberg_check)
    assert callable(bk.spectral.eigenvalues)


def test_quickstart_runs():
    G = bk.group()
    assert len(G.GENS) == 8
    assert bk.automorphisms().order == 48
    assert bk.lengths(max_half_trace=5)[0]["n_oriented"] == 24
    assert abs(bk.systole() - 3.0571418389619960) < 1e-13


def test_green_callable_module(monkeypatch):
    import bolzakit.green.surface as S
    monkeypatch.setattr(S, "_green", lambda T, D: ("stub", T, D))
    assert callable(bk.green)
    assert bk.green(T=2.6, D_ball=13.0) == ("stub", 2.6, 13.0)
    assert callable(bk.green.constants) and callable(bk.green.bergman_at_weierstrass)
    assert bk.green.VALIDITY_RADIUS == 0.45


def test_every_public_function_labelled():
    import inspect
    import bolzakit.arith.curve as curve
    exempt_names = {"SurfaceGreen", "Ring"}      # classes; constants (G1, AREA, VALIDITY_RADIUS, SU_URL, ...) are non-callable
    exempt_other = {"cache_dir"}                 # path helper, not a mathematical result
    mods = [bk.geometry, bk.spectral, bk.green, bk.vortex, bk.arith, curve]
    funcs = {}
    for m in mods:
        for n, o in vars(m).items():
            if n.startswith("_") or not callable(o) or inspect.ismodule(o):
                continue
            if inspect.isclass(o) or n in exempt_names | exempt_other:
                continue
            if n == "provenance" or getattr(o, "__module__", "").startswith("bolzakit.provenance"):
                continue
            if getattr(o, "__module__", "").split(".")[0] != "bolzakit":
                continue
            funcs[f"{m.__name__}.{n}"] = o
    funcs.update({"systole": bk.systole, "lengths": bk.lengths, "automorphisms": bk.automorphisms})
    assert len(funcs) >= 15
    missing = [k for k, f in funcs.items() if not hasattr(f, "__provenance__")]
    assert not missing, missing
    reg = registry()
    assert len(reg) >= 15
    assert all(i.label in LABELS for i in reg.values())
