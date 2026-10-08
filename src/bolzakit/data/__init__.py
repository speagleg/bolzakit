"""Bundled data files (versioned JSON). See MANIFEST.json for checksums and generating commands."""
from __future__ import annotations
import json
from importlib import resources


def load(name: str) -> dict:
    with resources.files(__package__).joinpath(f"{name}.json").open("r", encoding="utf-8") as fh:
        return json.load(fh)
