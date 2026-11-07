"""The campaign's dot set: same libraries and parameters give the same set; groups follow the material's facet
labels (synthetic material with families F, PA, PB, where PA/PB are the two terminations of one polar family)."""
import json

import numpy as np

from qd_workflows.shapes import select_set


def _lib(tmp_path, name, seed, n=40):
    rng = np.random.default_rng(seed)
    S = []
    for i in range(n):
        kind = ["F", "PA", "PB", "PA+PB", "F+PA+PB"][i % 5]
        fams = kind.split("+")
        areas = {f: float(rng.uniform(50, 500)) for f in fams}
        edges = {"|".join(sorted((a, b))): {"n": 4, "length_A": float(rng.uniform(10, 80))}
                 for a in fams for b in fams if a <= b and (a != b or len(fams) == 1)}
        S.append({"id": f"{name}-{i:03d}", "fingerprint": f"{name}{i}", "n_atoms": int(rng.integers(40, 580)),
                  "m": i + 10, "k": 4 + i % 7,
                  "geometry": {"areas_A2": areas, "edges": edges, "n_vertices": int(rng.integers(4, 30))}})
    d = tmp_path / name
    d.mkdir()
    (d / "summary.json").write_text(json.dumps({"structures": S}))
    return str(d)


PARAMS = {"fit": {"tiers": [{"max_atoms": 300, "n": 8}, {"max_atoms": 600, "n": 4}], "per_shape": 3},
          "groups": {"max_atoms": 600, "n_pure": 3, "min_pure": 2, "dominated_share": 0.6, "n_dominated": 2},
          "strip": {"per_group": 1, "mixed": 1, "atoms": [40, 600], "mixed_min_families": 3}}


def test_select_is_reproducible_and_generic(tmp_path):
    libs = [_lib(tmp_path, "a", 1), _lib(tmp_path, "b", 2)]
    r1 = select_set(libs, PARAMS)
    r2 = select_set(libs, PARAMS)
    assert r1["ids"] == r2["ids"] and r1["strip"] == r2["strip"]
    assert set(r1["groups"]) == {"pure F", "pure PA", "pure PB", "pure PA+PB"}      # polar pair found from the labels
    assert len(r1["strip"]) == 5                                                     # one per group + one mixed
    assert set(r1["libraries"]) == {str((tmp_path / n).resolve()) for n in ("a", "b")}


def test_dominated_group_when_pure_dots_are_missing(tmp_path):
    lib = _lib(tmp_path, "c", 3)
    S = json.loads((tmp_path / "c" / "summary.json").read_text())["structures"]
    S = [s for s in S if list(s["geometry"]["areas_A2"]) != ["F"]]                  # no pure F dots (like CdSe cubes)
    (tmp_path / "c" / "summary.json").write_text(json.dumps({"structures": S}))
    r = select_set([lib], PARAMS)
    assert "pure F" not in r["groups"] and "dominated F" in r["groups"]
