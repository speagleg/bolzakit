#!/usr/bin/env bash
# Local (no-lake, no-download) build of BolzaKit.LocalMax, for machines where `lake exe cache get` is not wanted.
# CI still uses `lake build`.  This compiles the 8 LocalMax files with plain `lean` against prebuilt Mathlib oleans.
#
#   PKGS=/path/to/copied/packages  ./build_local.sh     # PKGS/<pkg>/lean must hold each package's lib/lean oleans
#   (default: read <repo>/lean/.lake/packages/*/.lake/build/lib/lean directly; slow over /mnt/c, ~2 h)
# Copying the oleans to a local disk first makes the build take a few minutes.
# OUT (default ./.local-build) receives our own oleans.  Logs the #print axioms output to $OUT/axioms.log and runs
# the axiom gate on it.
set -euo pipefail
cd "$(dirname "$0")"
export PATH=$HOME/.elan/bin:$PATH
OUT=${OUT:-$(pwd)/.local-build}
if [ -n "${PKGS:-}" ]; then
  LP=$(for d in "$PKGS"/*; do printf "%s:" "$d/lean"; done)
else
  P=${PKGDIR:-$(pwd)/.lake/packages}
  LP=$(for d in "$P"/*; do printf "%s:" "$d/.lake/build/lib/lean"; done)
fi
mkdir -p "$OUT/BolzaKit/LocalMax"
export LEAN_PATH="${LP}${OUT}"
for f in RepTheory LocalMax Taylor RepReal Certificates KyFan Assembly Model; do
  echo "== $f"
  nice lean -o "$OUT/BolzaKit/LocalMax/$f.olean" -i "$OUT/BolzaKit/LocalMax/$f.ilean" "BolzaKit/LocalMax/$f.lean"
done
echo "== BolzaKit.LocalMax (root + #print axioms)"
nice lean BolzaKit/LocalMax.lean 2>&1 | tee "$OUT/axioms.log"
bash ../scripts/check_axioms.sh "$OUT/axioms.log"
