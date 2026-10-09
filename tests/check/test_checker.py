import hashlib
import json
import subprocess
import sys

import pytest

flint = pytest.importorskip("flint")
from bolzakit.certified import Surface, certify_spectrum, count_eigenvalues
from bolzakit.certified.chart import certify_rank
from bolzakit.certified.cluster import GapInput, Lam1Input, certify_cluster
from bolzakit_check.check import check_file

LEGACY = Lam1Input(lower="3.8388862587999999", upper="3.8388882588000001", source="legacy SU13 interval")
GAP = GapInput(lower="4.0", source="trace-formula LP count (certified)")
# Synthetic lam1 input just below the n=24 P1 trial Rayleigh quotient mu = 3.8516016 (fixed ARPACK v0): the gap
# transfer error is then small and the n=24 first-order certificate certifies.  This exercises the checker's VALID
# path for cluster_first_order cheaply; it is not a statement about Bolza (the certificate is conditional on D1).
NEAR_MU = Lam1Input(lower="3.85150", upper="3.85159", source="synthetic (checker test)")


@pytest.fixture(scope="module")
def certs(tmp_path_factory):
    d = tmp_path_factory.mktemp("c")
    out = {
        "rr": certify_spectrum(Surface.bolza(), k=3, degree=2, n=12).certificate,
        "count": count_eigenvalues("bolza", window=(3.3, 4.0)).certificate,      # mode B (conditional)
        "countA": count_eigenvalues("bolza", window=(0, 3.3)).certificate,      # mode A
        "rank": certify_rank().certificate,
        "first_order": certify_cluster("F_A", degree=1, n=24, lam1=NEAR_MU, gap=GAP).certificates[0],
        "kyfan": certify_cluster("F_B", degree=2, n=24, lam1=LEGACY, gap=GAP, second_order=True).certificates[0],
    }
    paths = {}
    for k, c in out.items():
        assert c.certified, k
        p = d / f"{k}.json"
        c.save(p)
        paths[k] = p
    return paths


def _mutate(p, fn, tmp_path, resha=False):
    raw = json.loads(p.read_text())
    fn(raw)
    if resha:
        core = {k: raw[k] for k in ("schema", "kind", "inputs", "data")}
        raw["sha256"] = hashlib.sha256(json.dumps(core, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    q = tmp_path / ("m_" + p.name)
    q.write_text(json.dumps(raw))
    return q


def test_valid(certs):
    for k, p in certs.items():
        ok, why = check_file(p)
        assert ok, (k, why)
        assert why.startswith("VALID")


def test_mode_b_is_only_conditionally_valid(certs):
    ok, why = check_file(certs["count"])
    assert ok and why.startswith("VALID (conditional on lambda_1 >= 33/10)") and "assumes: lambda_1 >= 33/10" in why
    ok, why = check_file(certs["countA"])
    assert ok and why == "VALID"


def test_mode_b_without_hypothesis_is_invalid(certs, tmp_path):
    q = _mutate(certs["count"], lambda r: r["inputs"].pop("assumes"), tmp_path, resha=True)
    ok, why = check_file(q)
    assert not ok and "assumes" in why


def test_tamper_detected(certs, tmp_path):
    q = _mutate(certs["rr"], lambda r: r["data"].__setitem__("mu", "3.0"), tmp_path)
    ok, why = check_file(q)
    assert not ok and "sha256 mismatch" in why


def _bump_ball(s, by):
    return str(flint.arb(s) + by)


def _kyfan_bump(r):
    G = r["data"]["raw"]["G"]["K2"]
    G[0][0] = _bump_ball(G[0][0], 100)


def _bump_J(r):
    r["data"]["J"][0][0] = "[0 +/- 0]"
    r["data"]["J"][1] = r["data"]["J"][0]


TAMPERS = {
    "rr": [lambda r: r["data"].__setitem__("mu", "3.0"),
           lambda r: r["data"]["K"][1].__setitem__(1, _bump_ball(r["data"]["K"][1][1], 10))],
    "count": [lambda r: r["data"].__setitem__("target", 3),
              lambda r: r["data"]["c_coeffs"].__setitem__(0, "1/3")],
    "countA": [lambda r: r["data"].__setitem__("a", "4")],
    "rank": [_bump_J],
    "first_order": [lambda r: r["data"]["transfer_terms"].__setitem__(
        "mu", _bump_ball(r["data"]["transfer_terms"]["mu"], 0.1))],
    "kyfan": [_kyfan_bump,
              lambda r: r["data"]["raw"].__setitem__("lam", ["3.0", "3.1"])],
}


@pytest.mark.parametrize("kind,i", [(k, i) for k, fs in TAMPERS.items() for i in range(len(fs))])
def test_resha_after_tamper_still_invalid(certs, tmp_path, kind, i):
    q = _mutate(certs[kind], TAMPERS[kind][i], tmp_path, resha=True)
    ok, why = check_file(q)
    assert not ok and "sha" not in why.lower(), why


def test_flipped_claim_is_invalid(certs, tmp_path):
    # claims are outside the hashed core: a stored holds must still match the re-derivation
    q = _mutate(certs["rank"], lambda r: r["claims"][0].__setitem__("holds", False), tmp_path)
    ok, why = check_file(q)
    assert not ok and "re-derived True" in why
    q = _mutate(certs["kyfan"], lambda r: r["claims"][0].__setitem__("rhs", "-100"), tmp_path)
    assert not check_file(q)[0]


def test_cli_exit_codes(certs, tmp_path):
    r = subprocess.run([sys.executable, "-m", "bolzakit_check", str(certs["rr"])], capture_output=True, text=True)
    assert r.returncode == 0 and "VALID" in r.stdout
    q = _mutate(certs["rr"], lambda r: r["data"].__setitem__("mu", "3.0"), tmp_path)
    r = subprocess.run([sys.executable, "-m", "bolzakit_check", str(certs["rr"]), str(q)], capture_output=True,
                       text=True)
    assert r.returncode == 1 and f"{q}: INVALID" in r.stdout and f"{certs['rr']}: VALID" in r.stdout


def test_deep_mode_labels_consistency_only(certs, tmp_path):
    ok, why = check_file(certs["rr"], deep=True)
    assert ok and "consistency only" in why
    sample = json.loads(certs["rr"].read_text())["data"]["element_sample"]
    assert sample["n_elements"] == 1152 and len(sample["items"]) == 11

    def bad(r):
        K = r["data"]["element_sample"]["items"][0]["K"]
        K[0][0] = _bump_ball(K[0][0], 1e-3)
    q = _mutate(certs["rr"], bad, tmp_path, resha=True)
    assert check_file(q)[0]                      # the certified claims do not depend on the sample
    ok, why = check_file(q, deep=True)
    assert not ok and "consistency only" in why
    # a hand-picked sample (not the documented deterministic selection) is rejected
    q = _mutate(certs["rr"], lambda r: r["data"]["element_sample"]["items"].pop(), tmp_path, resha=True)
    ok, why = check_file(q, deep=True)
    assert not ok and "selection" in why


def test_trial_provenance_recorded(certs):
    for k in ("rr", "first_order", "kyfan"):
        t = json.loads(certs[k].read_text())["inputs"]["trial"]
        assert t["arpack_v0"].startswith("RandomState(0).rand(") and len(t["sha256"]) == 64


def _strict(text):
    def no_const(name):
        raise ValueError(f"non-strict JSON constant {name}")
    return json.loads(text, parse_constant=no_const)


def _fail_closed(monkeypatch):
    """Fail-closed certificates of every kind (NaN / infinite intermediate data)."""
    import bolzakit.certified.chart as chart
    import bolzakit.certified.cluster as C
    import bolzakit.certified.forms as F
    import bolzakit.certified.rayleigh_ritz as rr
    out = []
    with monkeypatch.context() as m:
        def no_upper(*a, **k):
            raise RuntimeError("could not certify an upper bound")
        m.setattr(rr, "certified_upper", no_upper)
        out.append(certify_spectrum(Surface.bolza(), k=3, degree=2, n=12).certificate)
    with monkeypatch.context() as m:
        real = C.element_dual

        def poisoned(*a):
            Ke, Me, dKe, dMe, ba, br = real(*a)
            return Ke, Me, dKe, dMe, flint.arb("inf"), br
        m.setattr(C, "element_dual", poisoned)
        out.append(certify_cluster("F_A", degree=1, n=12, lam1=LEGACY, gap=GAP).certificates[0])
    with monkeypatch.context() as m:
        real2 = F.element_p2

        def poisoned2(*a, **k):
            Ks, Ms, bet = real2(*a, **k)
            Ks = [[row[:] for row in Mx] for Mx in Ks]
            Ks[1][0][1] = Ks[1][1][0] = flint.arb("nan")
            return Ks, Ms, bet
        m.setattr(F, "element_p2", poisoned2)
        out.append(certify_cluster("F_B", degree=2, n=12, lam1=LEGACY, gap=GAP, second_order=True).certificates[0])
    with monkeypatch.context() as m:
        real3 = chart.jacobians

        def fake(x0):
            Jc, Jf, cval, fval = real3(x0)
            Jc = [list(row) for row in Jc]
            Jc[0][1] = flint.arb("nan")
            return Jc, Jf, cval, fval
        m.setattr(chart, "jacobians", fake)
        out.append(certify_rank().certificate)
    # NaN square roots in the main paths (lam1 ball reaching above the trial Rayleigh quotients)
    wide = Lam1Input(lower="3.5", upper="3.99", source="test")
    out.append(certify_cluster("F_A", degree=1, n=12, lam1=wide, gap=GapInput("4.5", "test")).certificates[0])
    out.append(certify_cluster("F_B", degree=2, n=12, lam1=wide, gap=GapInput("4.5", "test"),
                               second_order=True).certificates[0])
    return out


def test_every_certificate_is_strict_json(certs, monkeypatch, tmp_path):
    for p in certs.values():
        _strict(p.read_text())
    for i, c in enumerate(_fail_closed(monkeypatch)):
        assert not c.certified
        text = c.to_json()
        raw = _strict(text)
        assert raw["sha256"] == c.sha256()
        p = tmp_path / f"fail{i}.json"
        c.save(p)
        ok, why = check_file(p)
        assert not ok and why.startswith("INVALID") and "sha" not in why.lower(), why
