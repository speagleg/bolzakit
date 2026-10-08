#!/usr/bin/env bash
set -euo pipefail
if grep -rn --include=*.lean -w sorry lean/BolzaKit lean/BolzaKit.lean; then echo "sorry found"; exit 1; fi
echo "no sorry"
