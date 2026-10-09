"""Certificate object, schema bolzakit-cert/1."""
from __future__ import annotations

import hashlib
import importlib
import json
import platform
import subprocess
from dataclasses import asdict, dataclass, field
from pathlib import Path


@dataclass
class Claim:
    lhs: str
    op: str  # one of "<", "<=", ">", ">=", "!=0", "encloses"
    rhs: str
    holds: bool


@dataclass
class Dependency:
    id: str
    statement: str
    source: str


@dataclass
class Certificate:
    kind: str  # "rr_upper" | "count" | "cluster_first_order" | "kyfan_second_order" | "rank" | "lg_lower"
    inputs: dict
    claims: list
    data: dict
    dependencies: list
    producer: dict = field(default_factory=dict)
    schema: str = "bolzakit-cert/1"

    @property
    def certified(self) -> bool:
        return len(self.claims) > 0 and all(c.holds for c in self.claims)

    def canonical_json(self) -> str:
        return json.dumps(
            {"schema": self.schema, "kind": self.kind, "inputs": self.inputs, "data": self.data},
            sort_keys=True, separators=(",", ":"), allow_nan=False)

    def sha256(self) -> str:
        return hashlib.sha256(self.canonical_json().encode("utf-8")).hexdigest()

    def to_json(self) -> str:
        d = {
            "schema": self.schema, "kind": self.kind, "inputs": self.inputs,
            "claims": [asdict(c) for c in self.claims], "data": self.data,
            "dependencies": [asdict(x) for x in self.dependencies],
            "producer": self.producer, "certified": self.certified, "sha256": self.sha256(),
        }
        return json.dumps(d, sort_keys=True, indent=2, allow_nan=False)  # strict JSON: no NaN / Infinity

    def save(self, path) -> None:
        Path(path).write_text(self.to_json() + "\n", encoding="utf-8")

    @classmethod
    def from_json(cls, s: str) -> "Certificate":
        d = json.loads(s)
        return cls(kind=d["kind"], inputs=d["inputs"],
                   claims=[Claim(**c) for c in d["claims"]], data=d["data"],
                   dependencies=[Dependency(**x) for x in d["dependencies"]],
                   producer=d.get("producer", {}), schema=d.get("schema", "bolzakit-cert/1"))

    @classmethod
    def load(cls, path) -> "Certificate":
        return cls.from_json(Path(path).read_text(encoding="utf-8"))


def trial_provenance(vectors, v0_rule="RandomState(0).rand(Nf)") -> dict:
    """Provenance of the (untrusted) float trial vectors: the ARPACK start vector rule and the sha256 of the vectors,
    stacked row-wise (one row per vector) as little-endian float64 bytes."""
    import numpy as np
    a = np.ascontiguousarray(np.array([np.asarray(v, dtype=float) for v in vectors]), dtype="<f8")
    return {"arpack_v0": v0_rule, "sha256": hashlib.sha256(a.tobytes()).hexdigest()}


def make_producer() -> dict:
    import flint
    # Not a static import: the certified layer may not import the top-level bolzakit package (layering test).
    __version__ = importlib.import_module("bolzakit").__version__
    here = Path(__file__).resolve()
    git = "unknown"
    try:
        # Only describe a checkout that actually tracks this file: a non-editable install inside an unrelated git
        # checkout (e.g. a venv in some other repository) must not record that repository's describe.
        tracked = subprocess.run(["git", "ls-files", "--error-unmatch", here.name], cwd=here.parent,
                                 capture_output=True, text=True, timeout=10)
        if tracked.returncode == 0:
            r = subprocess.run(["git", "describe", "--always", "--dirty"], cwd=here.parent, capture_output=True,
                               text=True, timeout=10)
            if r.returncode == 0 and r.stdout.strip():
                git = r.stdout.strip()
    except (OSError, subprocess.TimeoutExpired):
        git = "unknown"
    return {"bolzakit": __version__, "git": git, "flint": flint.__version__,
            "python": platform.python_version()}
