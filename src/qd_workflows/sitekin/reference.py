"""Step 2, reference (GPU): relax each built dot, no anneal (annealing moves atoms, and the sites could no longer
be identified). L-BFGS from the builder geometry, then rattled by rattle_A and relaxed again while that lowers it by
more than polish_gain_eV and the dot stays in one piece (at most polish_rounds); float64 energy. All sizes in one
batch. Writes reference/<tag>.json {positions, E, steps}."""
from __future__ import annotations

import json
import zlib
from pathlib import Path

import numpy as np
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

from . import recipe as rc


def one_piece(x, cut=3.2) -> bool:
    return connected_components(cKDTree(x).sparse_distance_matrix(cKDTree(x), cut), directed=False)[0] == 1


def relax_dots(engine, dots: dict, p: dict, log=print) -> dict:
    """dots: {tag: (symbols, positions)}. Returns {tag: {positions, E, steps, polish_rounds}}."""
    names = list(dots)
    # one rattle stream per dot (seed, tag): the same reference whichever dots share the batch
    rng = {k: np.random.default_rng([int(p.get("seed", 13)), zlib.crc32(k.encode())]) for k in names}
    res = engine.relax_many([(dots[k][0], np.asarray(dots[k][1], float)) for k in names], fmax=float(p["fmax"]),
                            max_steps=int(p["max_steps"]))
    cur = {k: (np.asarray(r.positions, float), float(r.energy), int(r.steps)) for k, r in zip(names, res)}
    rounds = {k: 0 for k in names}
    active = list(names)
    for _ in range(int(p.get("polish_rounds", 3))):
        if not active:
            break
        res = engine.relax_many([(dots[k][0], cur[k][0] + rng[k].normal(0, float(p.get("rattle_A", 0.05)), cur[k][0].shape))
                                 for k in active], fmax=float(p["fmax"]), max_steps=1500)
        nxt = []
        for k, r in zip(active, res):
            if r.energy < cur[k][1] - float(p.get("polish_gain_eV", 0.02)) and one_piece(r.positions):
                cur[k] = (np.asarray(r.positions, float), float(r.energy), cur[k][2] + int(r.steps))
                rounds[k] += 1
                nxt.append(k)
        active = nxt
    e64 = engine.ev64([(dots[k][0], cur[k][0]) for k in names])
    out = {}
    for k, (e, _f) in zip(names, e64):
        out[k] = {"positions": np.round(cur[k][0], 5).tolist(), "E": float(e), "steps": cur[k][2], "polish_rounds": rounds[k]}
        log(f"{k}: relaxed, E64 {float(e):.4f} eV, {cur[k][2]} steps, {rounds[k]} polish rounds kept")
    return out


def run(r: dict, out: Path, engine, sizes=None, log=print) -> list[Path]:
    out = Path(out)
    (out / "reference").mkdir(parents=True, exist_ok=True)
    dots = {}
    for e in r["sizes"]:
        t = rc.tag(e["size"])
        if (sizes and t not in sizes) or (out / "reference" / f"{t}.json").exists():
            continue
        d = json.loads((out / "dots" / f"{t}.json").read_text())
        dots[t] = (d["symbols"], d["positions"])
    paths = []
    if dots:
        for t, v in relax_dots(engine, dots, r["reference"], log).items():
            p = out / "reference" / f"{t}.json"
            p.write_text(json.dumps({"tag": t, **v}))
            paths.append(p)
    return paths
