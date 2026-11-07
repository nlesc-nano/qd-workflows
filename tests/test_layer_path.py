"""Layer path on a toy lineage edge: the 14-plane seed cut (the scan's Cd16Se13Cl6 planes) and a child with one {111}A
and one {111}B plane one layer out. The builder tests (a few layer states and one island) need QD_builder (`builder`,
the qd_builder env); the rest is lattice bookkeeping."""
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

from qd_workflows.layer_path import islands, lattice as L, treepaths

CIF = Path(__file__).parent / "data" / "cd16" / "CdSe_zb.cif"
SEED = {n: (4 if f in ("100", "111B") else 5) for f, n in L.PLANES}       # levels of the scan seed, central Se at 0
MOVES = [("111A", (1, 1, -1)), ("111B", (1, 1, 1))]
needs_builder = pytest.mark.skipif(importlib.util.find_spec("builder") is None, reason="QD_builder not installed")


def _cut(levels):
    LG, Ls = L.lattice_box(np.zeros((1, 3), int), pad=12)
    m = L.inside(LG, list(levels.items()))
    return Ls[m], LG[m]


def _struct(levels):
    s, G = _cut(levels)
    sym = np.where(s == "A", "Se", "Cd")
    order = np.argsort(sym != "Se", kind="stable")                      # a Se first: the sites' reference
    return {"xyz": [[str(e), "x", *map(float, g * L.Q)] for e, g in zip(sym[order], G[order])]}, int((s == "A").sum())


def toy_data(path):
    child = dict(SEED)
    for f, n in MOVES:
        child[n] += L.SPACING[f]
    (s0, m0), (s1, m1) = _struct(SEED), _struct(child)
    nodes = [{"label": "seed", "dot": 0, "parents": [], "children": [1], "size": m0},
             {"label": "child", "dot": 1, "parents": [0], "children": [], "size": m1}]
    path.write_text(json.dumps({"structs": [s0, s1], "lineage": {"nodes": nodes}}))
    return path


def test_class_and_dot_keys_are_symmetry_invariant():
    s, G = _cut(SEED)
    M, t = L.TD_OPS[7], np.array([4, 0, -4])
    assert L.canon(s, G) == L.canon(s, G @ M.T + t)
    sym = np.where(s == "A", "Se", "Cd")
    assert L.dot_key(sym, G * L.Q) == L.dot_key(sym, (G @ M.T + t) * L.Q)
    assert L.canon(s, G) != L.canon(*_cut({**SEED, (1, 1, 1): 8}))


def test_vicinal_planes():
    assert sorted(islands.vicinal("100", (1, 0, 0))) == [(2, -1, 0), (2, 0, -1), (2, 0, 1), (2, 1, 0)]
    assert sorted(islands.vicinal("111B", (1, 1, 1))) == [(1, 1, 2), (1, 2, 1), (2, 1, 1)]


def test_toy_edge_layer_states(tmp_path):
    data = toy_data(tmp_path / "data.json")
    treepaths._init({"data": str(data), "smax": 10})
    se = treepaths.size_edge((0, 1))
    assert (se["n_moves"], se["layers"], se["states"], se["neg"]) == (2, 2, 4, 0)
    r = treepaths.classes((0, 1, 14))
    assert r["mode"] == "all sites" and r["start_ok"] and r["end_ok"]
    assert sorted(r["moves"]) == ["111A(1, 1, -1)", "111B(1, 1, 1)"] and r["nlay"] == [1, 1]
    assert len(r["classes"]) == 4                         # seed, +111A, +111B, child: four distinct site sets
    assert len(r["edges"]) == 4 and r["first"] != r["last"]
    sizes = sorted(v["sites"] for v in r["classes"].values())
    assert sizes[0] == 35 and sizes[-1] == len(_cut({**SEED, (1, 1, -1): 9, (1, 1, 1): 8})[0])
    for v in r["classes"].values():                        # every plane in the builder recipe, far out if faceless
        assert len(v["planes"]) == 14 and min(h for _n, h in v["planes"]) > 0


@needs_builder
def test_toy_treepaths_and_one_island(tmp_path):
    data = toy_data(tmp_path / "data.json")
    treepaths.size_edges(data, tmp_path / "edges_size.json", nproc=1)
    net = tmp_path / "net"
    reg = treepaths.treepaths(data, net, mmax=100, smax=10, cif=CIF, edges_size=tmp_path / "edges_size.json", nproc=1)
    assert len(reg["classes"]) == 4 and all(v is not None for v in reg["class_dot"].values())
    e = json.loads((net / "edges" / "0_1.json").read_text())
    assert e["start_is_parent_cd"] and e["end_is_child_cd"] and len(e["dots"]) >= 3
    byid = {v["id"]: v for v in reg["dots"].values()}
    for d in e["dots"]:                                    # charge-balanced builder dots
        sy, _x = L.read_xyz(net / "builds" / byid[d]["dir"] / "out.xyz")
        assert 2 * (sy == "Cd").sum() - 2 * (sy == "Se").sum() - (sy == "Cl").sum() == 0
    # one island on a {111}B facet of each path dot: the largest one short of the full layer (on dots this small the
    # prepass strips most islands back to the parent; at least one must keep its new Cd)
    islands._init({"out": str(tmp_path / "isl"), "cif": str(CIF), "per": 1})
    n_cd = lambda fml: int(fml.split("Cd")[1].split("Se")[0])
    grown = []
    for parent in e["dots"]:
        rows = islands.states((parent, "111B", str(net / "builds" / byid[parent]["dir"])))
        assert rows and all(r["parent_cut_ok"] for r in rows)
        sub = [r for r in rows if r["s_sites"] < r["full_layer_Cd"]]
        if not sub:
            continue
        b = islands.build_one(max(sub, key=lambda r: r["s_sites"]))
        assert b["ok"] and b["Q"] == 0
        if n_cd(b["formula"]) > n_cd(byid[parent]["formula"]):
            grown.append(b)
    assert grown and all(b["key"] not in reg["dots"] for b in grown)
