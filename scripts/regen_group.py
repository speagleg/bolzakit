"""Regenerate src/bolzakit/data/group.json. Usage: python scripts/regen_group.py"""
import json, datetime
from bolzakit.exact import group as G

out = {"schema": 1, "generated": datetime.date.today().isoformat(),
       "convention": "element (A,B) in Z[zeta8]^2 <-> [[A, mu B],[mu conj(B), conj(A)]], mu = sqrt(2+2sqrt2); "
                     "Z[zeta8] element = (a0,a1,a2,a3) meaning sum a_k zeta8^k",
       "generators": [[list(a), list(b)] for a, b in G.GENS], "relator": list(G.RELATOR),
       "generator_half_trace": [1, 1]}
with open("src/bolzakit/data/group.json", "w", encoding="utf-8", newline="\n") as fh:
    json.dump(out, fh, indent=1)
    fh.write("\n")
