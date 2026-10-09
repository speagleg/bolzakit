#!/usr/bin/env bash
# Usage: check_axioms.sh LOGFILE
# Fail (exit 1) if LOGFILE has no `depends on axioms: [...]` list, or if any such list names an axiom outside
# {propext, Classical.choice, Quot.sound}.  Lean prints "does not depend on any axioms" for axiom-free theorems;
# those lines are fine but do not count toward the "at least one" requirement.
# Lean wraps long axiom lists across lines, so the log is joined into one line before it is parsed. A list whose
# closing bracket is missing (a truncated log) fails the gate.
set -euo pipefail
log="${1:?usage: check_axioms.sh LOGFILE}"
joined="$(tr '\r\n' '  ' < "$log")"
opened="$(grep -o 'depends on axioms: \[' <<< "$joined" | wc -l || true)"
n=0; bad=0
while IFS= read -r item; do
  [ -n "$item" ] || continue
  n=$((n+1))
  inner="${item#*depends on axioms: [}"; inner="${inner%]}"
  IFS=',' read -ra names <<< "$inner"
  for a in "${names[@]}"; do
    a="${a//[[:space:]]/}"
    case "$a" in
      propext|Classical.choice|Quot.sound) ;;
      *) echo "disallowed axiom '$a' in: $item" >&2; bad=1 ;;
    esac
  done
done < <(grep -o 'depends on axioms: \[[^]]*\]' <<< "$joined" || true)
if [ "$n" -eq 0 ]; then echo "no 'depends on axioms' lists in $log" >&2; exit 1; fi
if [ "$n" -ne "$opened" ]; then echo "unterminated 'depends on axioms: [' list in $log" >&2; exit 1; fi
[ "$bad" -eq 0 ] || exit 1
echo "axiom gate OK ($n theorems checked)"
