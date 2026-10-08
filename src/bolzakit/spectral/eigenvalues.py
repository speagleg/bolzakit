"""Laplace eigenvalues of the Bolza surface (area 4π, curvature −1)."""
from __future__ import annotations
import hashlib, os, pathlib, urllib.request
import mpmath as mp
from .. import data as _data
from ..errors import DataUnavailable
from ..provenance import provenance

SU_URL = "https://arxiv.org/src/1110.2150v4/anc/eig-bolza-refined0-1000.txt"
SU_SHA256 = "869f48bcc1403cd2ae722b9c7d8ae94e016a4030a353d49c68ce5e6eb4520345"
SU_NAME = "eig-bolza-refined0-1000.txt"
AREA = 4 * mp.pi


def cache_dir() -> pathlib.Path:
    d = pathlib.Path(os.environ.get("BOLZAKIT_CACHE", pathlib.Path.home() / ".cache" / "bolzakit"))
    d.mkdir(parents=True, exist_ok=True)
    return d


def _download(dest: pathlib.Path) -> None:
    with urllib.request.urlopen(SU_URL, timeout=60) as r:
        dest.write_bytes(r.read())


def _verified(path: pathlib.Path) -> bool:
    return hashlib.sha256(path.read_bytes()).hexdigest() == SU_SHA256


@provenance("QUOTED", "Strohmaier & Uski, Commun. Math. Phys. 317 (2013), arXiv:1110.2150v4 ancillary file",
            "as stated by the authors (see docs/spectral.md)")
def su2013(path: str | None = None, dps: int = 50) -> list:
    """The 1000 non-zero eigenvalues of Strohmaier–Uski (2013), downloaded on first use and sha256-checked."""
    p = pathlib.Path(path) if path else cache_dir() / SU_NAME
    if path is None and (not p.exists() or not _verified(p)):
        if p.exists():
            p.unlink()
        try:
            _download(p)
        except Exception as e:  # network or filesystem
            raise DataUnavailable(f"could not download {SU_URL} ({e}); download it manually and pass path=") from e
    if not _verified(p):
        if path is None:
            p.unlink()
        raise DataUnavailable(f"checksum mismatch for {p} (expected sha256 {SU_SHA256[:12]}…)")
    with mp.workdps(dps):
        return [mp.mpf(s) for s in p.read_text().split()]


def _clusters(vals, rel=1e-9):
    out = []
    for v in vals:
        if out and abs(v - out[-1][0]) <= rel * abs(v):
            out[-1][1] += 1
        else:
            out.append([v, 1])
    return [(float(v), m) for v, m in out]


@provenance("NUMERICAL", "FEM (bundled, item 02) or Strohmaier–Uski 2013 (QUOTED, downloaded)",
            "FEM: relative error <= 2e-4 vs Strohmaier–Uski")
def eigenvalues(source: str = "fem", n: int | None = None) -> list[tuple[float, int]]:
    """Non-zero Laplace eigenvalue clusters (value, multiplicity), ascending.

    The decorator labels this function NUMERICAL as a whole; ``source="su2013"`` values are QUOTED
    (see :func:`su2013`).
    """
    if source == "fem":
        cl = [(c["value"], c["multiplicity"]) for c in _data.load("eigenvalues_fem")["clusters"]]
        if n is not None and n > len(cl):
            raise DataUnavailable(f"FEM bundle has {len(cl)} clusters; use source='su2013' for more")
    elif source == "su2013":
        cl = _clusters(su2013())
        if n is not None and n > len(cl):
            raise DataUnavailable(f"Strohmaier–Uski data has {len(cl)} clusters")
    else:
        raise ValueError("source must be 'fem' or 'su2013'")
    return cl[:n] if n is not None else cl


def weyl_counting(lam: float) -> float:
    """Weyl law leading term N(λ) ≈ A λ / 4π = λ for the Bolza surface (A = 4π)."""
    return float(AREA * lam / (4 * mp.pi))
