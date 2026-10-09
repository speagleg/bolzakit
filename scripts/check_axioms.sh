#!/usr/bin/env bash
# Usage: check_axioms.sh LOGFILE
# Fail (exit 1) if LOGFILE has no `depends on axioms: [...]` line, or if any such line names an axiom outside
# {propext, Classical.choice, Quot.sound}.  Lean prints "does not depend on any axioms" for axiom-free theorems;
# those lines are fine but do not count toward the "at least one" requirement.
set -euo pipefail
log="${1:?usage: check_axioms.sh LOGFILE}"
n=0; bad=0
while IFS= read -r line; do
  n=$((n+1))
  inner="${line#*depends on axioms: [}"; inner="${inner%%]*}"
  IFS=',' read -ra names <<< "$inner"
  for a in "${names[@]}"; do
    a="${a//[[:space:]]/}"
    case "$a" in
      propext|Classical.choice|Quot.sound) ;;
      *) echo "disallowed axiom '$a' in: $line" >&2; bad=1 ;;
    esac
  done
done < <(grep 'depends on axioms: \[' "$log" || true)
if [ "$n" -eq 0 ]; then echo "no 'depends on axioms' lines in $log" >&2; exit 1; fi
[ "$bad" -eq 0 ] || exit 1
echo "axiom gate OK ($n theorems checked)"
