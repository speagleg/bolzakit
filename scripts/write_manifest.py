"""Write src/bolzakit/data/MANIFEST.json (sha256 of every data file). Run LAST, after all data files are final."""
import hashlib, json, pathlib, datetime
import bolzakit

D = pathlib.Path("src/bolzakit/data")
SCRIPTS = {"group.json": "scripts/regen_group.py", "lengths.json": "scripts/regen_lengths.py",
           "automorphisms.json": "scripts/regen_automorphisms.py", "eigenvalues_fem.json": "scripts/regen_fem.py",
           "tau.json": "scripts/regen_tau.py", "green_constants.json": "scripts/regen_green_constants.py"}
man = {}
for f in sorted(D.glob("*.json")):
    if f.name == "MANIFEST.json":
        continue
    man[f.name] = {"sha256": hashlib.sha256(f.read_bytes()).hexdigest(), "script": SCRIPTS.get(f.name, "?"),
                   "command": f"python {SCRIPTS.get(f.name, '?')}", "bolzakit_version": bolzakit.__version__,
                   "date": datetime.date.today().isoformat()}
with open(D / "MANIFEST.json", "w", encoding="utf-8", newline="\n") as fh:
    json.dump(man, fh, indent=1)
    fh.write("\n")
print(list(man))
