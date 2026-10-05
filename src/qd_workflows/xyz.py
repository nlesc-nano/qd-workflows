"""Plain XYZ reading and writing."""
from __future__ import annotations

from pathlib import Path

import numpy as np


def read_frames(path) -> list[tuple[list[str], np.ndarray, str]]:
    """All frames of an XYZ file as (symbols, coordinates in Å, comment)."""
    lines = Path(path).read_text().splitlines()
    frames, i = [], 0
    while i < len(lines) and lines[i].strip():
        n = int(lines[i].split()[0])
        comment = lines[i + 1] if i + 1 < len(lines) else ""
        rows = [lines[i + 2 + k].split() for k in range(n)]
        frames.append(([r[0] for r in rows], np.array([[float(x) for x in r[1:4]] for r in rows]), comment))
        i += n + 2
    if not frames:
        raise ValueError(f"no XYZ frame in {path}")
    return frames


def read_xyz(path) -> tuple[list[str], np.ndarray]:
    symbols, coords, _ = read_frames(path)[0]
    return symbols, coords


def write_xyz(path, symbols, coords, comment: str = "") -> None:
    rows = [f"{s:<3s} {x:16.8f} {y:16.8f} {z:16.8f}" for s, (x, y, z) in zip(symbols, np.asarray(coords))]
    Path(path).write_text(f"{len(symbols)}\n{comment}\n" + "\n".join(rows) + "\n")
