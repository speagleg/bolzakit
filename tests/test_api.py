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


def test_every_public_function_labelled():
    reg = registry()
    assert len(reg) >= 15
    assert all(i.label in LABELS for i in reg.values())
