"""Write src/bolzakit/data/green_constants.json: closed-form Weierstrass-point constants as 30-digit strings.

Usage (from the repo root):  python scripts/regen_green_constants.py
"""
import json
import os

import mpmath as mp

from bolzakit import green

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "src", "bolzakit", "data", "green_constants.json")

c = green.constants(dps=30)
out = {k: {"value": mp.nstr(v["value"], 30), "closed_form": v["closed_form"], "label": v["label"]}
       for k, v in c.items()}
with open(OUT, "w", newline="\n", encoding="utf-8") as fh:
    json.dump(out, fh, indent=1)
    fh.write("\n")
print(json.dumps(out, indent=1))
