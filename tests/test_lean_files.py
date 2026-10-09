import pathlib, re, subprocess
ROOT = pathlib.Path(__file__).resolve().parents[1]

def test_localmax_files_and_no_sorry():
    files = list((ROOT / "lean/BolzaKit/LocalMax").glob("*.lean"))
    assert len(files) == 8
    for f in files + [ROOT / "lean/BolzaKit/LocalMax.lean"]:
        t = f.read_text(encoding="utf-8")
        assert not re.search(r"\bsorry\b|\badmit\b|native_decide|^axiom ", t, re.M), f

def test_axiom_gate_script(tmp_path):
    ok = tmp_path / "ok.log"; ok.write_text("'Bolza.x' depends on axioms: [propext, Classical.choice, Quot.sound]\n")
    bad = tmp_path / "bad.log"; bad.write_text("'Bolza.x' depends on axioms: [propext, sorryAx]\n")
    empty = tmp_path / "e.log"; empty.write_text("")
    run = lambda p: subprocess.run(["bash", str(ROOT / "scripts/check_axioms.sh"), str(p)]).returncode
    assert run(ok) == 0 and run(bad) == 1 and run(empty) == 1
