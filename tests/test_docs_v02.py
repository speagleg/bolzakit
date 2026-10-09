# tests/test_docs_v02.py
import pathlib, re, subprocess, sys, pytest
ROOT = pathlib.Path(__file__).resolve().parents[1]


def _read(p):
    return (ROOT / p).read_text(encoding="utf-8")


def test_readme_sections():
    t = _read("README.md")
    for s in ["For software engineers", "For physicists", "For engineers", "bolzakit_check", "proofs/bolza-local-max"]:
        assert s in t, s
    # the honest headline status of the proof
    for s in ["conditional on the cited inputs", "No human expert has read", "optional"]:
        assert s in t, s


def test_docs_exist_with_sections():
    c = _read("docs/certified.md")
    for s in ["Trust model", "Validity limits", "not checked"]:
        assert s in c
    s = _read("docs/sensitivity.md")
    assert "discrete" in s and "Ky Fan" in s


def test_certified_md_states_known_limitations():
    c = _read("docs/certified.md")
    for s in ["|t| ≤ 0.05", "Bolza-only", "minimum-norm", "v0.3", "slow", "--deep", "consistency",
              "derivation-independent", "λ₁ ≥ a", "λ₄ ≥ 4.0", "hand", "wrong ρ"]:
        assert s in c, s


def test_sensitivity_md_near_tie_caveat():
    s = _read("docs/sensitivity.md")
    assert "diag(1, 1+δ, 3)" in s and "membrane_repeated_frequency.py" in s
    core = _read("src/bolzakit/sensitivity/core.py")
    assert "O(spread of the cluster)" not in core      # the old overclaim (R2)


def test_versions_consistent():
    cff = _read("CITATION.cff")
    assert "version: 0.2.0" in cff or 'version: "0.2.0"' in cff
    assert "2026-10-09" in cff
    assert "## [0.2.0]" in _read("CHANGELOG.md") or "## 0.2.0" in _read("CHANGELOG.md")
    import bolzakit
    assert bolzakit.__version__ == "0.2.0"


def test_verification_rows_cover_golden_certificates():
    v = _read("VERIFICATION.md")
    for f in (ROOT / "proofs/bolza-local-max/certificates").glob("*.json"):
        assert f.name in v, f.name


def test_no_stale_plate_example_references():
    assert not (ROOT / "examples/plate_repeated_frequency.py").exists()
    for p in ["README.md", "CHANGELOG.md", "docs/sensitivity.md", "docs/certified.md", ".github/workflows/ci.yml"]:
        assert "plate_repeated_frequency" not in _read(p), p


def test_ci_tiers():
    ci = _read(".github/workflows/ci.yml")
    for s in [".[dev,cert,skfem]", "python -m bolzakit_check proofs/bolza-local-max/certificates/*.json",
              "--runslow", "reproduce.py --quick"]:
        assert s in ci, s
    lean = _read(".github/workflows/lean.yml")
    gate = lean[lean.index("Axiom gate"):]
    assert "shell: bash" in gate.split("- name:")[0]


@pytest.mark.parametrize("ex", ["membrane_repeated_frequency.py", "bolza_deform.py", "bolza_strain_rates.py"])
def test_examples_run(ex):
    pytest.importorskip("flint"); pytest.importorskip("skfem")
    r = subprocess.run([sys.executable, str(ROOT / "examples" / ex)], capture_output=True, text=True, timeout=600)
    assert r.returncode == 0, r.stderr


def _slug(heading):
    # GitHub's heading anchors: lower case, punctuation dropped (letters, digits, '-', '_' and spaces kept), spaces -> '-'
    return re.sub(r"[^\w\- ]", "", heading.strip().lower()).replace(" ", "-")


def test_readme_in_page_anchors_resolve():
    t = (ROOT / "README.md").read_text(encoding="utf-8")
    slugs = {_slug(h) for h in re.findall(r"^#{1,6} (.+)$", t, re.M)}
    links = re.findall(r"\]\(#([^)]+)\)", t)
    assert links and [a for a in links if a not in slugs] == []


def test_theorem_headings_are_conditional():
    t = (ROOT / "README.md").read_text(encoding="utf-8")
    assert "## The Bolza surface is a strict local maximum of λ₁ (computer-assisted, conditional on cited inputs)" in t
    first = (ROOT / "proofs" / "bolza-local-max" / "README.md").read_text(encoding="utf-8").splitlines()[0]
    assert "computer-assisted, conditional" in first
