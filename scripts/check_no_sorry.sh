#!/usr/bin/env bash
# Fail if any Lean source uses `sorry`, `admit` or `native_decide` (whole words, comments included).
set -euo pipefail
cd "$(dirname "$0")/.."
if grep -rn --include=*.lean -w -e sorry -e admit -e native_decide lean/BolzaKit lean/BolzaKit.lean; then
  echo "sorry/admit/native_decide found"; exit 1
fi
echo "no sorry"
