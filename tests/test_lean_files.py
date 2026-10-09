import pathlib, re, subprocess
ROOT = pathlib.Path(__file__).resolve().parents[1]

def test_localmax_files_and_no_sorry():
    files = list((ROOT / "lean/BolzaKit/LocalMax").glob("*.lean"))
    assert len(files) == 8
    for f in files + [ROOT / "lean/BolzaKit/LocalMax.lean"]:
        t = f.read_text(encoding="utf-8")
        assert not re.search(r"\bsorry\b|\badmit\b|native_decide|^\s*((private|protected|noncomputable|unsafe|scoped)\s+)*axiom\b",
                             t, re.M), f

MAIN = "'Bolza.bolza_strict_local_max' depends on axioms: [propext, Classical.choice, Quot.sound]\n"


def test_axiom_gate_script(tmp_path):
    ok = tmp_path / "ok.log"
    ok.write_text(MAIN + "'Bolza.x' depends on axioms: [propext, Classical.choice, Quot.sound]\n")
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
    wrapped_ok.write_text(MAIN + "'Bolza.x' depends on axioms: [propext,\n Classical.choice,\n Quot.sound]\n"
                          "'Bolza.y' does not depend on any axioms\n")
    unclosed = tmp_path / "u.log"; unclosed.write_text("'Bolza.x' depends on axioms: [propext,\n")
    assert run(wrapped_bad) == 1 and run(wrapped_ok) == 0 and run(unclosed) == 1


def test_axiom_gate_needs_main_theorem_and_no_errors(tmp_path):
    run = lambda p: subprocess.run(["bash", str(ROOT / "scripts/check_axioms.sh"), str(p)]).returncode
    other = tmp_path / "o.log"
    other.write_text("'Bolza.x' depends on axioms: [propext, Classical.choice, Quot.sound]\n")
    err = tmp_path / "err.log"
    err.write_text(MAIN + "BolzaKit/LocalMax/KyFan.lean:3:0: error: unknown identifier 'foo'\n")
    lookalike = tmp_path / "l.log"
    lookalike.write_text(MAIN.replace("bolza_strict_local_max", "bolza_strict_local_maxX"))
    assert run(other) == 1 and run(err) == 1 and run(lookalike) == 1
