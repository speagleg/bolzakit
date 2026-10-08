# Regenerate src/bolzakit/data/eigenvalues_fem.json.
# Adapted from math-pursuits item 02 checks/fem_extrapolate.py: runs fem.py at n=40 and n=80, groups eigenvalues
# into clusters, and does h^2 Richardson extrapolation on the cluster means: lam_ext = (4*lam_80 - lam_40)/3.
# Usage: python scripts/regen_fem.py   (takes minutes)
import json, subprocess, sys, os, re, numpy as np
here = os.path.dirname(os.path.abspath(__file__))
def run(n, k=30):
    out = subprocess.run([sys.executable, os.path.join(here, 'fem.py'), str(n), str(k)], capture_output=True, text=True, check=True).stdout
    nums = re.findall(r'[-+]?\d+\.\d+e[-+]\d+', out.split('\n', 1)[1])
    return np.array(sorted(float(x) for x in nums))[1:]   # drop the zero mode
def clusters(v, rel=0.008):
    cl = [[v[0]]]
    for x in v[1:]:
        (cl[-1].append(x) if abs(x - cl[-1][-1]) < rel * x else cl.append([x]))
    return cl
c40 = clusters(run(40)); c80 = clusters(run(80))
rows = []
for k, (a, b) in enumerate(zip(c40, c80), 1):
    if len(a) != len(b): print('cluster mismatch', k, a, b); break
    m40, m80 = np.mean(a), np.mean(b)
    rows.append({"value": float((4*m80-m40)/3), "multiplicity": len(b)})
    print(f"{k:>2} {len(a):>6} {len(b):>6} {m40:11.6f} {m80:11.6f} {(4*m80-m40)/3:11.6f}")
# The last cluster may be cut off by the k=30 eigenvalue limit (observed: x2 instead of x4), so drop it.
rows = rows[:-1]
out = {"schema": 1,
       "method": "P1 FEM on meshes n=40,80, h^2 Richardson on cluster means (item 02 checks/fem_extrapolate.py)",
       "rel_error_vs_SU2013": 2e-4, "clusters": rows}
path = os.path.join(here, "..", "src", "bolzakit", "data", "eigenvalues_fem.json")
with open(path, "w", encoding="utf-8", newline="\n") as fh:
    json.dump(out, fh, indent=1)
    fh.write("\n")
print(len(rows), "clusters ->", path)
