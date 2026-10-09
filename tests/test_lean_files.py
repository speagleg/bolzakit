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


def test_axiom_gate_joins_wrapped_lines(tmp_path):
    # Lean wraps long axiom lists; a disallowed axiom on a continuation line must still fail the gate.
    run = lambda p: subprocess.run(["bash", str(ROOT / "scripts/check_axioms.sh"), str(p)]).returncode
    wrapped_bad = tmp_path / "wb.log"
    wrapped_bad.write_text("'Bolza.bolza_strict_local_max' depends on axioms: [propext,\n Classical.choice,\n"
                           " sorryAx,\n Quot.sound]\n")
    wrapped_ok = tmp_path / "wo.log"
    wrapped_ok.write_text("'Bolza.x' depends on axioms: [propext,\n Classical.choice,\n Quot.sound]\n"
                          "'Bolza.y' does not depend on any axioms\n")
    unclosed = tmp_path / "u.log"; unclosed.write_text("'Bolza.x' depends on axioms: [propext,\n")
    assert run(wrapped_bad) == 1 and run(wrapped_ok) == 0 and run(unclosed) == 1
