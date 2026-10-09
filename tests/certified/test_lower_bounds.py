"""Certified lambda_1 lower bound (Lehmann-Goerisch): flat-torus validation, Bolza n = 24 / 48, negative controls and
checker tamper tests."""
import hashlib
import json

import numpy as np
import pytest

flint = pytest.importorskip("flint")
from flint import arb

from bolzakit.certified import GapInput
from bolzakit.certified import lower_bounds as LB
from bolzakit.certified._arb import ball_str, prec
from bolzakit_check.check import check_file

GAP = GapInput(lower="4.0", source="trace-formula LP count (certified)")
P2_UPPER = 3.8388872830                       # rr_p2_n96 (certified P2 upper bound)
SRC_N24 = 3.8387028952042492                  # item 23 lb_bolza_24.json, [3.83870289520424920936766646644 +/- 4.33e-38]
SRC_N48 = 3.8388757010527246                  # item 23 lb_bolza_48.json certified_lambda1_lower
VALID_MSG = "VALID given stored Gram data (conditional on lambda_4 >= 4.0)"


# ------------------------------------------------------------------ flat torus (port of item 23 lb_torus.py)
def _torus(N):
    lab = lambda i, j: (i % N) * N + (j % N)
    tris, geo = [], []
    h = arb(1) / N
    for i in range(N):
        for j in range(N):
            # vertices: local r=0 <-> (1,0) in (a,b), r=1 <-> (0,1), r=2 <-> (0,0)
            for (P0, P1, P2) in [((i + 1, j), (i + 1, j + 1), (i, j)), ((i + 1, j + 1), (i, j + 1), (i, j))]:
                tris.append([lab(*P0), lab(*P1), lab(*P2)])
                Jm = [[(P0[0] - P2[0]) * h, (P1[0] - P2[0]) * h], [(P0[1] - P2[1]) * h, (P1[1] - P2[1]) * h]]
                g11 = Jm[0][0] ** 2 + Jm[1][0] ** 2; g12 = Jm[0][0] * Jm[0][1] + Jm[1][0] * Jm[1][1]
                g22 = Jm[0][1] ** 2 + Jm[1][1] ** 2
                D = abs(Jm[0][0] * Jm[1][1] - Jm[0][1] * Jm[1][0])
                geo.append(({(0, 0): arb(1)}, ({(0, 0): g11}, {(0, 0): g12}, {(0, 0): g22}), D))
    return np.array(tris), N * N, geo


def _torus_lb(N, rho_m1=78.9, ntrial=4):
    """lambda_1(R^2/Z^2) = 4 pi^2 (multiplicity 4), lambda_5 = 8 pi^2; rho = 1 + 78.9 <= 1 + 8 pi^2, m = 5, k = 4."""
    with prec(LB.PREC):
        tris, nv, geo = _torus(N)
        topo = LB.Topo(tris, nv)
        assert topo.chi == 0
        cache, elems = {}, []
        for g in geo:      # two congruent triangle shapes only
            key = (str(g[1][0][(0, 0)]), str(g[1][1][(0, 0)]), str(g[1][2][(0, 0)]))
            if key not in cache:
                cache[key] = LB.element(*g)
            elems.append(cache[key])
        res, *_, R, V, w = LB.lg_bound(topo, elems, ntrial, arb(1 + rho_m1), np.random.RandomState(0).rand(topo.np2))
        assert res["ok"] and LB.green_ok(R)
        return res["lower_bound_k_n"] - 1, 4 * arb.pi() ** 2


def test_flat_torus_validation():
    published = {4: 39.2041077545, 8: 39.4582495177, 16: 39.4770992124}   # item 23 lb_torus.out
    prev = None
    for N, pub in published.items():
        lb, exact = _torus_lb(N)
        assert lb < exact                                      # certified: whole ball below 4 pi^2
        assert abs(float(lb.mid()) - pub) < 1e-9
        if prev is not None:                                   # converging upward, O(h^4)
            assert prev < lb and (exact - lb) < (exact - prev) / 8
        prev = lb


# ------------------------------------------------------------------ Bolza
@pytest.fixture(scope="module")
def n24(tmp_path_factory):
    # sym=True (sectors 2-7 reuse 0-1) keeps the default suite fast; it reproduces SRC's SYM=0 value to 1e-12
    r = LB.certify_lam1_lower(n=24, rho=5, gap=GAP, sym=True)
    p = tmp_path_factory.mktemp("lg") / "lg_lower_n24.json"
    r.certificate.save(p)
    return r, p


def test_bolza_n24(n24):
    r, p = n24
    assert r.certified and r.lam1_lower is not None
    assert r.lam1_lower <= r.lam1_lower_inprocess <= P2_UPPER
    assert abs(r.lam1_lower_inprocess - SRC_N24) < 1e-12
    c = r.certificate
    assert c.kind == "lg_lower" and c.data["rho"] == "5" and c.data["shift"] == 1 and c.data["k"] == 3
    assert {d.id for d in c.dependencies} == {"D15", "gap"}
    assert c.claims[0].lhs == "lambda_1" and c.claims[0].op == ">="
    assert float(c.claims[0].rhs) <= r.lam1_lower
    assert set(c.inputs["trial"]) == {"arpack_v0", "sha256"}
    assert check_file(p) == (True, VALID_MSG)


@pytest.mark.slow
def test_bolza_n48(tmp_path):
    r = LB.certify_lam1_lower(n=48, rho=5, gap=GAP)
    assert r.certified
    assert abs(r.lam1_lower_inprocess - SRC_N48) < 1e-12
    assert r.lam1_lower <= r.lam1_lower_inprocess <= P2_UPPER
    p = tmp_path / "lg_lower_n48.json"
    r.certificate.save(p)
    assert check_file(p) == (True, VALID_MSG)


# ------------------------------------------------------------------ negative controls
def test_gap_hypothesis_enforced():
    with pytest.raises(ValueError, match="hypothesis"):
        LB.certify_lam1_lower(n=24, rho=5, gap=GapInput(lower="3.9", source="too weak"))
    with pytest.raises(ValueError, match="hypothesis"):
        LB.certify_lam1_lower(n=24, rho=7.5, gap=GAP)
    with pytest.raises(ValueError):
        LB.certify_lam1_lower(n=24, rho=5, gap=None)


def test_flipped_glued_edge_is_not_certified(monkeypatch):
    """Flip the normal component of one edge glued by a side pairing (second owner's flux DOFs): sigma is no longer
    in H(div) of the closed surface, A4 fails, and the result must not certify (B indefinite / Green residual)."""
    base = LB.Topo
    flipped = []

    class Flipped(base):
        def __init__(self, tri_labels, nverts):
            super().__init__(tri_labels, nverts)
            per = self.Nt // 8
            for key, (t1, t2) in self.owner.items():
                if abs(t1 // per - t2 // per) == 4:            # sectors k and k+4: a side-pairing edge
                    r = list(self.p2loc[t2][3:6]).index(nverts + self.ekey[key])
                    self.rtsign[t2, 3 * r:3 * r + 3] *= -1
                    flipped.append(key)
                    break
    monkeypatch.setattr(LB, "Topo", Flipped)
    r = LB.certify_lam1_lower(n=24, rho=5, gap=GAP, sym=True)
    assert flipped
    assert not r.certified and r.lam1_lower is None
    d = r.certificate.data["diagnostics"]
    assert d["reason"] in ("B not positive definite", "Green residual not negligible")
    assert d["green_residual_max"] > 1e-6


def test_nan_element_is_not_certified(monkeypatch, tmp_path):
    base = LB.element
    calls = []

    def poisoned(*a):
        K, Mm, S = base(*a)
        if not calls:
            K[0, 0] = arb("nan")
        calls.append(1)
        return K, Mm, S
    monkeypatch.setattr(LB, "element", poisoned)
    r = LB.certify_lam1_lower(n=24, rho=5, gap=GAP, sym=True)
    assert not r.certified and r.lam1_lower is None and r.lam1_lower_inprocess is None
    r.certificate.save(tmp_path / "nan.json")                  # strict JSON
    ok, why = check_file(tmp_path / "nan.json")
    assert not ok and "not certified" in why


# ------------------------------------------------------------------ checker tamper tests
def _tamper(p, fn, tmp_path, name):
    raw = json.loads(p.read_text())
    fn(raw)
    core = {k: raw[k] for k in ("schema", "kind", "inputs", "data")}
    raw["sha256"] = hashlib.sha256(json.dumps(core, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    q = tmp_path / f"{name}.json"
    q.write_text(json.dumps(raw))
    return check_file(q)


def _scale_diag(key, f):
    def fn(raw):
        M = raw["data"][key]
        with prec(256):
            for i in range(3):
                M[i][i] = ball_str(arb(M[i][i]) * arb(f))
    return fn


def _set_rho(v):
    return lambda raw: raw["data"].__setitem__("rho", v)


def _raise_claim(raw):            # claim (and enclosure) 1e-9 stronger than certified
    lo = float(arb(raw["data"]["lam1_lower_enclosure"]).mid()) + 1e-9
    raw["data"]["lam1_lower_enclosure"] = repr(lo)
    raw["claims"][0]["rhs"] = repr(lo)


@pytest.mark.parametrize("name,fn", [
    ("A0_up", _scale_diag("A0", "1.00000001")),      # larger A0: A less negative, nu_max closer to 0
    ("A1_down", _scale_diag("A1", "0.99999999")),    # smaller A1: same direction
    ("A2_up", _scale_diag("A2", "1.00000001")),      # larger A2: B larger, |nu| smaller
    ("rho_6", _set_rho("6")),                        # stronger rho, violates rho <= lambda_4 + 1 = 5
    ("rho_49_10", _set_rho("49/10")),                # admissible but weaker rho: the stored claim is not reached
    ("claim_up", _raise_claim),
])
def test_tamper_strengthening_is_invalid(n24, tmp_path, name, fn):
    ok, why = _tamper(n24[1], fn, tmp_path, name)
    assert not ok and why.startswith("INVALID"), (name, why)


def test_tamper_weakening_stays_valid(n24, tmp_path):
    def weaken(raw):              # a weaker claim (and enclosure) from the same Gram data is still certified
        raw["data"]["lam1_lower_enclosure"] = "3.8"
        raw["claims"][0]["rhs"] = "3.8"
    assert _tamper(n24[1], weaken, tmp_path, "weaker") == (True, VALID_MSG)
