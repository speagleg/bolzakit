import hashlib

from importlib import resources

from bolzakit import data
from bolzakit.exact import automorphisms as AU, group as G


def test_group_json_matches_code():
    d = data.load("group")
    assert d["relator"] == list(G.RELATOR)
    assert [[list(a), list(b)] for a, b in G.GENS] == d["generators"]


def test_manifest_checksums():
    man = data.load("MANIFEST")
    for fname, meta in man.items():
        raw = resources.files("bolzakit.data").joinpath(fname).read_bytes()
        assert hashlib.sha256(raw).hexdigest() == meta["sha256"], fname


def test_manifest_covers_all_data_files():
    man = data.load("MANIFEST")
    names = {p.name for p in resources.files("bolzakit.data").iterdir() if p.name.endswith(".json")}
    assert names - {"MANIFEST.json"} == set(man)


def test_automorphisms_json_matches_code():
    d = data.load("automorphisms")
    A = AU.automorphisms()
    assert d["order"] == A.order
    assert d["full_order"] == A.full_order
    assert {int(k): v for k, v in d["element_order_counts"].items()} == A.element_order_counts()
    W = A.weierstrass_points()
    Wj = [complex(a, b) for a, b in d["weierstrass_points"]]
    assert len(Wj) == len(W) == 6
    for w in Wj:
        assert any(AU.same_point_mod_gamma(w, z, tol=1e-9) for z in W)
    for z in W:
        assert any(AU.same_point_mod_gamma(w, z, tol=1e-9) for w in Wj)
