"""
Small islands on a facet as builder cuts (the scripts islands.py and stage_islands.py).

An island state is the parent's own recipe with one facet plane n0 moved one layer out, and vicinal planes that cut the
new layer back to an island bounded by steps: on {100} (n0 = e_x) v = (2, +-1, 0), (2, 0, +-1) ({210}); on {111}
v = n0 + n0_i e_i ({211}). A vicinal plane is never put inside the parent's lattice cut, so the parent stays whole and
the new layer keeps only the sites the vicinal planes leave. Every distinct island (Td ops + translation) of each
new-Cd count s (s <= scap, or the full layer) is built with QD_builder, at most `per` per size and normal.
Parents are registry dots (their builds/<dir>/recipe.yaml and out.xyz); parents.json: [[dot_id, family], ...].

islands.json rows: parent, family, n0, s_sites, vicinal, levels, full_layer_Cd, parent_cut_ok, tag, ok, formula, Q,
key (dot key), in_registry, parent_formula, dup_of.

stage: the unique built islands with Q = 0 that differ from the parent and contain its lattice (a pure addition,
steps.embed), 1 <= dCd <= 10, not yet relaxed, into LUMI batches (wave/batchJ/iNNNN.xyz, balanced by atoms^1.5);
writes islands_dir/states_index.json (status todo / relaxed / dup / parent, dSe, dCd, dCl, pure, name).
"""
from __future__ import annotations

import collections
import itertools
import json
import shutil
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import yaml

from . import lattice as L
from .build import build

SCAP = 14             # lattice island sizes kept (the prepass strips the smallest; the full layer is always kept)
_CTX = {}


def _init(ctx):
    _CTX.update(ctx)


def vicinal(f, n0):
    n0 = np.array(n0)
    if f == "100":
        i = int(np.flatnonzero(n0)[0])
        out = []
        for j in range(3):
            if j != i:
                for s in (1, -1):
                    v = 2 * n0.copy()
                    v[j] = s
                    out.append(tuple(int(x) for x in v))
        return out
    return [tuple(int(x) for x in n0 + n0[i] * np.eye(3, dtype=int)[i]) for i in range(3)]


def level(n, h):
    return int(np.floor(h * np.linalg.norm(n) / L.Q + 1e-6))


def to_h(n, d):
    return d * L.Q / np.linalg.norm(n) + L.MARGIN


def states(task):
    """Island states of one parent dot (id, family, build folder): rows with the builder planes."""
    pid, fam, wd0 = task
    per, scap = _CTX.get("per", 2), _CTX.get("scap", SCAP)
    wd0 = Path(wd0)
    rec = yaml.safe_load((wd0 / "recipe.yaml").read_text())
    planes = [(tuple(int(x) for x in f["hkl"]), float(f["gamma"])) for f in rec["facets"]]
    hs, hG = L.arr(L.sites(*L.read_xyz(wd0 / "out.xyz"), np.zeros(3)))
    LG, Ls = L.lattice_box(hG, pad=8)
    lv = {n: level(n, h) for n, h in planes}
    cut0 = L.inside(LG, list(lv.items()))
    ok0 = {tuple(g) for g in LG[cut0 & (Ls == "Cd")]} == {tuple(g) for g in hG[hs == "Cd"]}
    rows = []
    for n0 in L.NORMALS[fam]:
        if n0 not in lv or lv[n0] > 10 ** 4:
            continue
        lv1 = dict(lv)
        lv1[n0] += L.SPACING[fam]
        full = L.inside(LG, list(lv1.items()))
        newC = full & ~cut0 & (Ls == "Cd")
        if newC.sum() == 0:
            continue
        V = vicinal(fam, n0)
        rng = []
        for v in V:
            dmin = int((LG[cut0] @ np.array(v)).max())
            dfull = int((LG[full] @ np.array(v)).max())
            rng.append(range(dmin, max(dmin, dfull) + 1))
        seen = {}
        for lev in itertools.product(*rng):
            m = L.inside(LG, list(lv1.items()) + list(zip(V, lev)))
            s = int((m & ~cut0 & (Ls == "Cd")).sum())
            if s == 0:
                continue
            k = L.canon(Ls[m], LG[m])
            if k not in seen:
                seen[k] = (s, lev)
        bys = {}
        for _k, (s, lev) in sorted(seen.items(), key=lambda kv: (kv[1][0], sum(kv[1][1]))):
            if (s <= scap or s == int(newC.sum())) and len(bys.setdefault(s, [])) < per:
                bys[s].append(lev)
        for s, levs in bys.items():
            for lev in levs:
                pA = [(n, (to_h(n, lv1[n]) if n == n0 else h)) for n, h in planes] + \
                     [(v, to_h(v, d)) for v, d in zip(V, lev) if d < int((LG[full] @ np.array(v)).max())]
                have_n = {n for n, _h in pA}
                for v in V:      # the builder wants every symmetry-equivalent plane of a family: the rest far out
                    for p_ in set(itertools.permutations(v)):
                        for sg in itertools.product((1, -1), repeat=3):
                            w = tuple(int(a * c) for a, c in zip(p_, sg))
                            if w not in have_n:
                                have_n.add(w)
                                pA.append((w, L.FAR))
                rows.append({"parent": pid, "family": fam, "n0": list(n0), "s_sites": s, "vicinal": [list(v) for v in V],
                             "levels": [int(d) for d in lev], "full_layer_Cd": int(newC.sum()), "parent_cut_ok": bool(ok0),
                             "planes": [[list(n), float(h)] for n, h in pA]})
    return rows


def tag_of(r):
    return (f"p{r['parent']}_{r['family']}_{''.join('m' if c < 0 else 'p' for c in r['n0'])}_s{r['s_sites']}_"
            f"{'_'.join(map(str, r['levels']))}")


def build_one(r):
    tag = tag_of(r)
    wd = Path(_CTX["out"]) / "builds" / tag
    ok = build([(tuple(n), h) for n, h in r["planes"]], wd, _CTX["cif"])
    r = dict(r, tag=tag, ok=bool(ok))
    if ok:
        sy, x = L.read_xyz(wd / "out.xyz")
        n = {e: int((sy == e).sum()) for e in ("Cd", "Se", "Cl")}
        r.update(formula=f"Cd{n['Cd']}Se{n['Se']}Cl{n['Cl']}", Q=2 * n["Cd"] - 2 * n["Se"] - n["Cl"],
                 key=repr(L.dot_key(sy, x)))
    return r


def islands(net, parents, out, cif, per=2, nproc=4, scap=SCAP, limit=0):
    """Build the islands of the parents ([(dot_id, family), ...]) into out/builds; writes out/islands.json.
    limit > 0: build only the first `limit` island states (a smoke test)."""
    net, out = Path(net), Path(out)
    reg = json.loads((net / "registry.json").read_text())
    byid = {v["id"]: (k, v) for k, v in reg["dots"].items()}
    out.mkdir(parents=True, exist_ok=True)
    tasks = [(int(pid), fam, str(net / "builds" / byid[int(pid)][1]["dir"])) for pid, fam in parents]
    ctx = {"out": str(out), "cif": str(cif), "per": per, "scap": scap}
    with Pool(nproc, initializer=_init, initargs=(ctx,)) as pool:
        rows = [r for rr in pool.map(states, tasks) for r in rr]
        if limit:
            rows = rows[:limit]
        print(f"{len(rows)} island states from {len(tasks)} parents; building", flush=True)
        built = pool.map(build_one, rows)
    for r in built:
        k = r.get("key")
        r["in_registry"] = reg["dots"][k]["id"] if k in reg["dots"] else None
        r.pop("planes")
        r["parent_formula"] = byid[r["parent"]][1]["formula"]
    keys = {}
    for r in built:
        if r.get("key"):
            keys.setdefault(r["key"], r["tag"])
            r["dup_of"] = None if keys[r["key"]] == r["tag"] else keys[r["key"]]
    (out / "islands.json").write_text(json.dumps(built, indent=1))
    for r in built:
        print(f"{r['tag']:40s} {r.get('parent_formula')} -> {r.get('formula')} Q={r.get('Q')} s={r['s_sites']}"
              f"/{r['full_layer_Cd']} cut_ok={r['parent_cut_ok']} reg={r['in_registry']} dup={r.get('dup_of')}")
    return built


def stage(net, isd, results, wave, n_batches):
    """Unique built island states -> LUMI batches; writes isd/states_index.json."""
    from .steps import embed, parse
    net, isd = Path(net), Path(isd)
    reg = json.loads((net / "registry.json").read_text())
    byid = {v["id"]: k for k, v in reg["dots"].items()}
    res = json.loads(Path(results).read_text()) if Path(results).exists() else {}
    rows = json.loads((isd / "islands.json").read_text())
    uniq, out = {}, []
    for r in rows:
        if not r.get("ok") or r.get("Q") != 0:
            continue
        pk = byid[r["parent"]]
        if r["key"] == pk:
            r["status"] = "parent"
            continue
        pa, pb = parse(pk), parse(r["key"])
        r.update(dSe=pb[2]["Se"] - pa[2]["Se"], dCd=pb[2]["Cd"] - pa[2]["Cd"], dCl=pb[2]["Cl"] - pa[2]["Cl"],
                 pure=embed(pa[:2], pb[:2]) is not None)
        if r["key"] in uniq:
            r["status"] = "dup"
            r["same_as"] = uniq[r["key"]]
            continue
        uniq[r["key"]] = r["tag"]
        r["status"] = "relaxed" if r["key"] in res else "todo"
        out.append(r)

    def natoms(r):
        return int((isd / "builds" / r["tag"] / "out.xyz").read_text().split("\n", 1)[0])
    todo = [r for r in out if r["status"] == "todo" and r["pure"] and 1 <= r["dCd"] <= 10]
    todo.sort(key=lambda r: -natoms(r))
    load = [0] * n_batches
    for i, r in enumerate(todo):
        j = load.index(min(load))
        load[j] += natoms(r) ** 1.5
        d = Path(wave) / f"batch{j}"
        d.mkdir(parents=True, exist_ok=True)
        r["name"] = f"i{i:04d}"
        shutil.copy(isd / "builds" / r["tag"] / "out.xyz", d / f"{r['name']}.xyz")
    (isd / "states_index.json").write_text(json.dumps(rows, indent=1))
    print(f"{len(rows)} builds; {len(out)} distinct island states ({sum(r['pure'] for r in out)} pure additions), "
          f"{sum(r['status'] == 'relaxed' for r in out)} already relaxed, {len(todo)} staged in {n_batches} batches")
    print("added (dSe, dCd, dCl):", collections.Counter((r["dSe"], r["dCd"], r["dCl"]) for r in out).most_common(12))
    return rows
