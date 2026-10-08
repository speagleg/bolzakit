"""Regenerate src/bolzakit/data/automorphisms.json. Usage: python scripts/regen_automorphisms.py [out_path]"""
import json
import sys

from bolzakit.exact.automorphisms import automorphisms

A = automorphisms()
out = {
    "schema": 1,
    "order": A.order,
    "full_order": A.full_order,
    "element_order_counts": {str(k): v for k, v in A.element_order_counts().items()},
    "weierstrass_points": [[z.real, z.imag] for z in A.weierstrass_points()],
    "gl23": {str(i): [list(r) for r in m] for i, m in A.isomorphism_to_GL23().items()},
}
path = sys.argv[1] if len(sys.argv) > 1 else "src/bolzakit/data/automorphisms.json"
with open(path, "w", encoding="utf-8", newline="\n") as fh:
    json.dump(out, fh, indent=1)
    fh.write("\n")
print(len(out["gl23"]), "elements ->", path)
