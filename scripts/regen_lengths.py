"""Regenerate src/bolzakit/data/lengths.json. Usage: python scripts/regen_lengths.py [T_max]  (default 120; ~minutes)."""
import json, sys, datetime
from bolzakit.exact import lengths as LS

T = float(sys.argv[1]) if len(sys.argv) > 1 else 120.0
rows = LS.length_spectrum(T)
out = {"schema": 1, "max_half_trace": T, "generated": datetime.date.today().isoformat(),
       "convention": "n_oriented = oriented primitive closed geodesics with this half-trace; n_unoriented = n_oriented/2",
       "entries": [{"halftrace": list(r["halftrace"]), "length": r["length"], "n_oriented": r["n_oriented"],
                    "n_unoriented": r["n_oriented"] // 2, "residual": r["residual"]} for r in rows]}
out["max_residual"] = max(e["residual"] for e in out["entries"])
assert out["max_residual"] < 1e-7
path = "src/bolzakit/data/lengths.json" if len(sys.argv) < 3 else sys.argv[2]
with open(path, "w", encoding="utf-8", newline="\n") as fh:
    json.dump(out, fh, indent=1)
    fh.write("\n")
print(len(out["entries"]), "entries ->", path)
