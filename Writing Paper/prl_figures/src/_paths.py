"""_paths.py -- resolve the repository's source trees once, by name.

The figure modules pull their data straight from the gate code, which lives in
two sibling trees ("No-go theorem" and "New no-go theory"). Each module used to
reach them with its own `os.path.join(HERE, "..", "..", "..", ...)` chain, so
the imports depended on the figure files staying exactly three levels down and
on two directory names that contain spaces. Moving or renaming anything broke
them silently at import time, and the test module additionally relied on
`import fig1_classes` having already mutated sys.path as a side effect.

repo_root() instead walks upward until it finds the directory that holds both
trees, and add_sources() takes the trees by name.
"""
from __future__ import annotations

import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))

# marker directories that identify the repository root
_MARKERS = ("No-go theorem", "New no-go theory")


def repo_root(start: str = _HERE) -> str:
    d = os.path.abspath(start)
    while True:
        if all(os.path.isdir(os.path.join(d, m)) for m in _MARKERS):
            return d
        parent = os.path.dirname(d)
        if parent == d:
            raise RuntimeError(
                "could not locate the repository root above "
                f"{start!r}: expected a directory containing {_MARKERS}"
            )
        d = parent


def _tree(*parts: str) -> str:
    return os.path.join(repo_root(), *parts)


SOURCES = {
    "figures": _HERE,
    "nogo": lambda: _tree("No-go theorem", "src"),
    "phase": lambda: _tree("New no-go theory", "src"),
    "gate_a": lambda: _tree("New no-go theory", "PhaseO_observable_inheritance", "src"),
    "gate_b": lambda: _tree("New no-go theory", "GateB_superconducting_witness", "src"),
    "gate_c": lambda: _tree("New no-go theory", "GateC_material_independence", "src"),
    "gate_d": lambda: _tree("New no-go theory", "GateD_robustness_discriminability", "src"),
    "roomt": lambda: _tree("New no-go theory", "RoomT", "src"),
}


def add_sources(*keys: str) -> list:
    """Put the named source trees on sys.path and return the resolved paths."""
    out = []
    for key in ("figures",) + keys:
        entry = SOURCES[key]
        path = entry if isinstance(entry, str) else entry()
        if not os.path.isdir(path):
            raise RuntimeError(f"source tree {key!r} not found at {path}")
        if path not in sys.path:
            sys.path.insert(0, path)
        out.append(path)
    return out
