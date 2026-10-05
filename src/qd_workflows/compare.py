"""
Compare two qdprops / Orchestr.AI PROPS results (props/properties.json) value by value.

The Phase 0 gate: Cd16 and Cd68 rerun on the cluster must match the Mac
results. The two runs differ in device and precision (Apple GPU in float32 vs
CUDA in float64), so values are compared with tolerances, not for equality.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

# (key fragment, absolute tolerance); the first fragment found in a key's path wins.
TOLERANCES = (
    ("cm1", 2.0),            # frequencies, cm-1
    ("meV", 0.5),            # meV and meV/atom
    ("eV_A", 2e-3),          # forces, eV/Å
    ("_A", 2e-3),            # lengths, Å
    ("eV", 2e-3),            # energies, eV
    ("_K", 1.0),
)
DEFAULT_RTOL = 1e-3
SKIP = ("time", "seconds", "wall", "steps", "timestamp", "date", "sha256", "commit", "dirty", "device", "dtype",
        "torch", "mace", "xtb", "binary", "path", "file")


def _leaves(obj, prefix=""):
    if isinstance(obj, dict):
        for k, v in obj.items():
            yield from _leaves(v, f"{prefix}.{k}" if prefix else str(k))
    elif isinstance(obj, list) and obj and all(isinstance(x, (int, float)) for x in obj):
        for i, v in enumerate(obj):
            yield f"{prefix}[{i}]", v
    elif isinstance(obj, (int, float)) and not isinstance(obj, bool):
        yield prefix, obj


def _atol(path: str) -> float | None:
    for frag, tol in TOLERANCES:
        if frag in path:
            return tol
    return None


def compare(ref, new, rtol: float = DEFAULT_RTOL) -> dict:
    """Numeric leaves of the two `summary` blocks; each outside its tolerance is listed."""
    a = json.loads(Path(ref).read_text())["summary"]
    b = json.loads(Path(new).read_text())["summary"]
    va, vb = dict(_leaves(a)), dict(_leaves(b))
    failures, compared = [], 0
    for key in sorted(set(va) & set(vb)):
        if any(s in key.lower() for s in SKIP):
            continue
        x, y = float(va[key]), float(vb[key])
        if math.isnan(x) and math.isnan(y):
            continue
        compared += 1
        atol = _atol(key)
        ok = abs(x - y) <= (atol if atol is not None else 0.0) + rtol * max(abs(x), abs(y))
        if not ok:
            failures.append({"key": key, "ref": x, "new": y, "diff": y - x, "atol": atol})
    return {"compared": compared, "failed": len(failures), "only_in_ref": sorted(set(va) - set(vb))[:20],
            "only_in_new": sorted(set(vb) - set(va))[:20],
            "failures": sorted(failures, key=lambda f: -abs(f["diff"]))}
