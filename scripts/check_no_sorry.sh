#!/usr/bin/env bash
# Fail if any Lean source uses `sorry`, `admit` or `native_decide` (whole words, comments included), or declares an
# `axiom` (a line starting with `axiom`, possibly after modifiers such as `private`, so the word inside comments is
# ignored).
set -euo pipefail
cd "$(dirname "$0")/.."
if grep -rn --include=*.lean -w -e sorry -e admit -e native_decide lean/BolzaKit lean/BolzaKit.lean; then
  echo "sorry/admit/native_decide found"; exit 1
fi
if grep -rnE --include=*.lean '^\s*((private|protected|noncomputable|unsafe|scoped)\s+)*axiom\b' lean/BolzaKit lean/BolzaKit.lean; then
  echo "axiom declaration found"; exit 1
fi
echo "no sorry, no axiom"
