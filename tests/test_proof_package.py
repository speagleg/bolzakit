import json
import pathlib
import subprocess
import sys

import pytest

pytest.importorskip("flint")
ROOT = pathlib.Path(__file__).resolve().parents[1] / "proofs" / "bolza-local-max"
GOLDEN = {"lp_A_3.3.json", "lp_B_3.3_4.0.json", "rr_p1_n192.json", "rr_p2_n96.json",
          "first_order_FA_n192_su13.json", "first_order_FA_n192_lg.json",
          "kyfan_FB_n96_su13.json", "kyfan_FB_n96_lg.json", "rank.json", "lg_lower_n48.json"}


def _check(*files):
    return subprocess.run([sys.executable, "-m", "bolzakit_check", *map(str, files)], capture_output=True, text=True)


def test_files_present():
    for f in ["README.md", "PROOF.md", "NOVELTY.md", "LEAN.md", "LB_RESULTS.md", "LB_PLAN.md", "reproduce.py",
              "paper/main.tex", "paper/refs.bib", "paper/NOTES.md", "paper/main.pdf"]:
        assert (ROOT / f).exists(), f
    names = {p.name for p in (ROOT / "certificates").glob("*.json")}
    assert GOLDEN <= names
    assert len(list((ROOT / "reviews").glob("review_*.md"))) >= 6


def test_golden_certificates_valid():
    # exit 0 = all unconditionally VALID; 2 = some only conditionally VALID (the mode-B count); 1 = some INVALID
    files = sorted((ROOT / "certificates").glob("*.json"))
    assert GOLDEN <= {f.name for f in files}
    r = _check(*files)
    assert r.returncode in (0, 2), r.stdout + r.stderr
    assert "INVALID" not in r.stdout


def test_lp_b_is_conditional():
    r = _check(ROOT / "certificates" / "lp_B_3.3_4.0.json")
    assert r.returncode == 2 and "conditional" in r.stdout, r.stdout


def test_lg_lower_is_conditional():
    r = _check(ROOT / "certificates" / "lg_lower_n48.json")
    assert r.returncode == 2, r.stdout
    assert "VALID given stored Gram data (conditional on lambda_4 >= 4.0)" in r.stdout, r.stdout
    c = json.loads((ROOT / "certificates" / "lg_lower_n48.json").read_text(encoding="utf-8"))
    assert c["certified"] and float(c["claims"][0]["rhs"]) >= 3.8388757010   # covers the _lg certificates' lam1 input


def test_golden_lam1_sources():
    for name, src in [("first_order_FA_n192_lg.json", "Lehmann–Goerisch"), ("kyfan_FB_n96_lg.json", "Lehmann–Goerisch"),
                      ("first_order_FA_n192_su13.json", "SU13"), ("kyfan_FB_n96_su13.json", "SU13")]:
        c = json.loads((ROOT / "certificates" / name).read_text(encoding="utf-8"))
        assert c["certified"] and src in c["inputs"]["lam1"]["source"], name


HEADLINE = ("PROVED (computer-assisted), conditional on the cited inputs in work/PROOF.md §Dependencies: FBP Lemmas "
            "7.2–7.3 and §7.4 geodesic data, the Selberg trace formula, the Lehmann–Goerisch theorem (quoted from "
            "secondary sources), and standard facts. No longer depends on Strohmaier–Uski (SU13): λ₁ ≥ 3.8388757010 is "
            "now certified in-house (review_lb_1: ACCEPT). Every component second-read; integrated proof read: "
            "ACCEPT-WITH-FIXES, fixes applied.")


def test_headline_is_conditional():
    t = (ROOT / "README.md").read_text(encoding="utf-8")
    assert HEADLINE in t                      # verbatim from item 23's README Status row (SRC 53bc014)
    assert "PROVED (computer-assisted)" in t and "conditional" in t and "Strohmaier" in t
