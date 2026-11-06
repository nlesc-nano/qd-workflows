"""
One lineage edge in detail (the seed Cd16Se13Cl6 -> a child): overlay, layer moves, symmetry-distinct states, layer
graph and the builder states (the scripts overlay.py, layers.py, lattice_check.py, distinct.py, graph.py, build2.py
and dedupe.py).

- overlay: every embedding of the parent's Cd/Se in the child (Td op + lattice translation), the planes each moves;
- layers: one embedding as whole-layer moves of single {100}/{111} planes; checks that the parent is the child inside
  the parent's 14 planes and that moving all listed planes gives the child; writes layers_<child>.json with the Cd/Se
  framework of every subset of moves;
- lattice_check: is each layer state's framework the full lattice cut by its 14 planes (and with the child's {110})?
- distinct: layer states -> symmetry-distinct frameworks;
- path: the same on lattice sites (anions species-blind), every distinct state built with QD_builder (states.json);
- graph: class of every layer state and the edges between classes one layer apart (graph_<child>.json);
- unique: unique built dots over several path folders (dot key), unique.json.

data.json: the Wulff scan's (QD_builder runs/wulff_scan_cdse): structs[i].xyz rows [symbol, role, x, y, z];
lineage.nodes[j]: label, dot, parents, children, size.
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np

from . import lattice as L
from .build import build


def _nodes(data, child, parent=None):
    S, N = data["structs"], data["lineage"]["nodes"]
    p = N[0] if parent is None else next(n for n in N if n["label"] == parent)
    c = next(n for n in N if n["label"] == child)
    return S, p, c


def _parse_planes(dd):
    return {tuple(map(int, k.strip("()").split(","))): v for k, v in dd.items()}


def _levels(dp, moves, lev):
    pl = [(n, dp[n]) for _f, n in L.PLANES]
    for (f, n), l in zip(moves, lev):
        pl = [(nn, d + l * L.SPACING[f] if nn == n else d) for nn, d in pl]
    return pl


def overlay(data, child_labels=(), parent=None):
    """Print, per child of the parent, its distinct embeddings and the planes moved (in layers)."""
    S, N = data["structs"], data["lineage"]["nodes"]
    seed = N[0] if parent is None else next(n for n in N if n["label"] == parent)
    ps, pG = L.grid(S[seed["dot"]]["xyz"])
    want = set(child_labels)
    for ci in seed["children"]:
        node = N[ci]
        if want and node["label"] not in want:
            continue
        cs, cG = L.grid(S[node["dot"]]["xyz"])
        emb = L.embeddings(ps, pG, cs, cG)
        seen = {}                     # distinct embeddings up to the child's own symmetry: the sets of added atoms
        for M, t, PG in emb:
            seen.setdefault(frozenset(L.keyset(cs, cG) - L.keyset(ps, PG)), (M, t, PG))
        print(f"\n=== {node['label']}  {S[node['dot']].get('unit', '')}  centre {S[node['dot']].get('centre', '')}: "
              f"{len(emb)} embeddings, {len(seen)} distinct added sets")
        for k, (_M, _t, PG) in enumerate(seen.values()):
            rows = []
            for fam, ns in L.NORMALS.items():
                for n in ns:
                    dp, dc = int((PG @ np.array(n)).max()), int((cG @ np.array(n)).max())
                    if dc != dp:
                        rows.append(f"{fam}{n}: +{(dc - dp) / L.SPACING[fam]:g}")
            add = L.keyset(cs, cG) - L.keyset(ps, PG)
            print(f"  embedding {k}: +{sum(1 for a in add if a[0] == 'Cd')} Cd +{sum(1 for a in add if a[0] == 'Se')} Se;"
                  "  planes moved (in layers):")
            for fam in L.NORMALS:
                r = [x for x in rows if x.startswith(fam + "(")]
                if r:
                    print("    " + "  ".join(r))


def layers(data, child, out=None, k_emb=0, parent=None):
    """Parent -> child as whole-layer moves of single {100}/{111} planes (one embedding); writes layers_<child>.json."""
    S, seed, node = _nodes(data, child, parent)
    ps, pG = L.grid(S[seed["dot"]]["xyz"])
    cs, cG = L.grid(S[node["dot"]]["xyz"])
    seen = {}
    for _M, _t, PG in L.embeddings(ps, pG, cs, cG):
        seen.setdefault(frozenset(L.keyset(cs, cG) - L.keyset(ps, PG)), PG)
    PG = list(seen.values())[k_emb]
    dp = {n: int((PG @ np.array(n)).max()) for _f, n in L.PLANES}
    dc = {n: int((cG @ np.array(n)).max()) for _f, n in L.PLANES}
    moves = [(f, n) for f, n in L.PLANES if dc[n] > dp[n]]
    seed_cut = L.inside(cG, [(n, dp[n]) for _f, n in L.PLANES])
    ok_seed = L.keyset(cs[seed_cut], cG[seed_cut]) == L.keyset(ps, PG)
    full = L.inside(cG, [(n, dc[n]) for _f, n in L.PLANES])
    print(f"{child} embedding {k_emb}: seed = child inside seed planes: {ok_seed}; "
          f"child inside its own 14 planes: {bool(full.all())}")
    print("moves:", ", ".join(f"{f}{n} +{(dc[n] - dp[n]) / L.SPACING[f]:g} layer" for f, n in moves))
    nlay = [(dc[n] - dp[n]) // L.SPACING[f] for f, n in moves]
    n_cd0, n_se0 = int((ps == "Cd").sum()), int((ps == "Se").sum())
    rows = []
    for lev in itertools.product(*[range(k + 1) for k in nlay]):
        m = L.inside(cG, _levels(dp, moves, lev))
        rows.append((lev, int((cs[m] == "Cd").sum()), int((cs[m] == "Se").sum())))
    for lev, ncd, nse in rows:
        if sum(lev) <= 1:
            tag = "seed" if not sum(lev) else "{}{}".format(*moves[lev.index(1)])
            print(f"  {tag:22s} Cd{ncd} Se{nse}  (+{ncd - n_cd0} Cd, +{nse - n_se0} Se)")
    print(f"  {len(rows)} layer states (layers per plane {nlay}); last: Cd{rows[-1][1]} Se{rows[-1][2]}")
    rec = {"child": child, "embedding": k_emb, "dp": {str(k): v for k, v in dp.items()},
           "dc": {str(k): v for k, v in dc.items()}, "moves": [[f, list(n)] for f, n in moves],
           "states": [[list(s), a, b] for s, a, b in rows]}
    Path(out or f"layers_{child}.json").write_text(json.dumps(rec))
    return rec


def _child_frame(data, child, layers_file, parent=None):
    S, _seed, node = _nodes(data, child, parent)
    cs, cG = L.grid(S[node["dot"]]["xyz"])
    lay = json.loads(Path(layers_file or f"layers_{child}.json").read_text())
    return cs, cG, lay, _parse_planes(lay["dp"]), [(f, tuple(n)) for f, n in lay["moves"]]


def lattice_check(data, child, layers_file=None, parent=None):
    """Count the layer states whose framework is not the lattice cut by their 14 planes (and with the child's {110})."""
    cs, cG, lay, dp, moves = _child_frame(data, child, layers_file, parent)
    LG, _Ls = L.lattice_box(cG, pad=8)
    D110 = {n: int((cG @ np.array(n)).max()) for n in L.NORMALS["110"]}
    bad = bad110 = 0
    for lev, ncd, nse in lay["states"]:
        pl = _levels(dp, moves, lev)
        bad += int(L.inside(LG, pl).sum() != ncd + nse)
        bad110 += int(L.inside(LG, pl + [(n, D110[n]) for n in L.NORMALS["110"]]).sum() != ncd + nse)
    print(f"{child}: {len(lay['states'])} states; lattice cut by 14 planes differs from intended in {bad}; "
          f"with the child's 110 planes added: {bad110}")
    return bad, bad110


def distinct(data, child, layers_file=None, parent=None):
    """Symmetry-distinct Cd/Se frameworks of the layer states (Td ops + lattice translation)."""
    cs, cG, lay, dp, moves = _child_frame(data, child, layers_file, parent)
    by = {}
    for lev, ncd, nse in lay["states"]:
        m = L.inside(cG, _levels(dp, moves, lev))
        by.setdefault(L.canon_cdse(cs[m], cG[m]), []).append((tuple(lev), ncd, nse))
    print(f"{child}: {len(lay['states'])} layer states -> {len(by)} symmetry-distinct frameworks")
    lv = {}
    for v in by.values():
        lv.setdefault(sum(v[0][0]), []).append(f"Cd{v[0][1]}Se{v[0][2]}")
    for nl in sorted(lv):
        print(f"  {nl} layers: {len(lv[nl])} distinct  {sorted(set(lv[nl]))}")
    return by


def _site_edge(data, child, parent=None):
    """Parent and child on lattice sites, the first embedding, the plane moves and the lattice box."""
    S, seed, node = _nodes(data, child, parent)
    seedS, childS = L.from_struct(S[seed["dot"]]), L.from_struct(S[node["dot"]])
    ps, pG = L.arr(seedS)
    cs, cG = L.arr(childS)
    emb = [PG for _M, _t, PG in L.embeddings(ps, pG, cs, cG, anion="A")]
    PG = emb[0]
    dp = {n: int((PG @ np.array(n)).max()) for _f, n in L.PLANES}
    dc = {n: int((cG @ np.array(n)).max()) for _f, n in L.PLANES}
    moves = [(f, n) for f, n in L.PLANES if dc[n] > dp[n]]
    nlay = [(dc[n] - dp[n]) // L.SPACING[f] for f, n in moves]
    LG, Ls = L.lattice_box(cG)
    return dict(ps=ps, pG=pG, cs=cs, cG=cG, emb=emb, PG=PG, dp=dp, dc=dc, moves=moves, nlay=nlay, LG=LG, Ls=Ls)


def _state_classes(e):
    classes = {}
    for lev in itertools.product(*[range(k + 1) for k in e["nlay"]]):
        pl = {n: e["dp"][n] for _f, n in L.PLANES}
        for (f, n), l in zip(e["moves"], lev):
            pl[n] += l * L.SPACING[f]
        m = L.inside(e["LG"], list(pl.items()))
        classes.setdefault(L.canon(e["Ls"][m], e["LG"][m]), (lev, pl, m))
    return classes


def path(data, child, out, cif, parent=None):
    """Parent -> child layer path on lattice sites; every symmetry-distinct state built with QD_builder (states.json)."""
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    e = _site_edge(data, child, parent)
    ps, cs, PG, LG, Ls, dp, moves = e["ps"], e["cs"], e["PG"], e["LG"], e["Ls"], e["dp"], e["moves"]
    print(f"{child}: seed sites Cd{(ps == 'Cd').sum()} A{(ps == 'A').sum()}, child sites Cd{(cs == 'Cd').sum()} "
          f"A{(cs == 'A').sum()}; {len(e['emb'])} embeddings")
    print("moves:", ", ".join(f"{f}{n} +{k}" for (f, n), k in zip(moves, e["nlay"])))
    O = L.central_anion(PG, ps)                 # the seed's central Se: the builder's origin
    for nm, G0, s0, lev in (("seed", PG, ps, e["dp"]), ("child", e["cG"], cs, e["dc"])):
        m = L.inside(LG, list(lev.items()))
        same = L.keyset(Ls[m], LG[m]) == L.keyset(s0, G0)
        print(f"  lattice cut by the {nm}'s 14 planes == {nm} sites: {same}")
    classes = _state_classes(e)
    print(f"  {int(np.prod([k + 1 for k in e['nlay']]))} layer states -> {len(classes)} symmetry-distinct")
    rows = []
    for i, (lev, pl, m) in enumerate(sorted(classes.values(), key=lambda v: (sum(v[0]), v[0]))):
        want = {(s, *map(int, g)) for s, g in zip(Ls[m], LG[m] - O)}
        pA, active = L.planes_A(LG, pl, O)
        wd = out / f"s{i:02d}"
        ok = build(pA, wd, cif)
        r = {"i": i, "levels": [int(x) for x in lev], "n_layers": int(sum(lev)),
             "moved": ["{}{}".format(*moves[j]) for j, l in enumerate(lev) if l],
             "active_planes": ["{}{}".format(f, n) for f, n in L.PLANES if active[n] < L.OFF],
             "sites": f"Cd{sum(1 for k in want if k[0] == 'Cd')}A{sum(1 for k in want if k[0] == 'A')}"}
        if not ok:
            r["error"] = (wd / "builder.log").read_text()[-300:]
            rows.append(r)
            print(f"s{i:02d} build failed")
            continue
        sy, x = L.read_xyz(wd / "out.xyz")
        ncd, nse, ncl = (int((sy == el).sum()) for el in ("Cd", "Se", "Cl"))
        have = L.sites(sy, x, np.zeros(3))
        cut = set(L.sites(*L.read_xyz(wd / "out_cut.xyz"), np.zeros(3)))
        r.update(formula=f"Cd{ncd}Se{nse}Cl{ncl}", unit=f"[CdSe]{nse}(CdCl2){ncd - nse}",
                 Q=2 * ncd - 2 * nse - ncl, cut_ok=(cut == want),
                 lattice_cl=sum(1 for v in have.values() if v == "Cl"),
                 offlattice_cl=ncl - sum(1 for v in have.values() if v == "Cl"))
        stripped = []
        for k in sorted(want - set(have)):
            v = np.array(k[1:]) + O
            outs = ["{}{}".format(*mv) for mv in moves if v @ np.array(mv[1]) > dp[mv[1]]]
            stripped.append(f"{'Cd' if k[0] == 'Cd' else 'anion'}@" + ("+".join(outs) if outs else "seed"))
        r["stripped"] = stripped
        r["extra"] = len(set(have) - want)
        rows.append(r)
        print(f"s{i:02d} {r['n_layers']}L sites {r['sites']:9s} -> {r['formula']:14s} {r['unit']:19s} Q={r['Q']} "
              f"cut_ok={r['cut_ok']} lattCl={r['lattice_cl']} Z-Cl={r['offlattice_cl']} "
              f"stripped={len(stripped)} {sorted(set(stripped)) if stripped else ''}"
              + (f" extra={r['extra']}" if r["extra"] else ""), flush=True)
    res = {"child": child, "origin_grid": O.tolist(), "moves": [[f, list(n)] for f, n in moves], "states": rows}
    (out / "states.json").write_text(json.dumps(res, indent=1))
    return res


def graph(data, child, out=None, parent=None):
    """Class of every layer state (same order as path's states.json) and the edges between classes one layer apart
    (with the moves that realise each edge); writes graph_<child>.json."""
    e = _site_edge(data, child, parent)
    LG, Ls, nlay = e["LG"], e["Ls"], e["nlay"]
    cls, first = {}, {}
    for lev in itertools.product(*[range(k + 1) for k in nlay]):
        pl = {n: e["dp"][n] for _f, n in L.PLANES}
        for (f, n), l in zip(e["moves"], lev):
            pl[n] += l * L.SPACING[f]
        m = L.inside(LG, list(pl.items()))
        k = L.canon(Ls[m], LG[m])
        cls[lev] = k
        first.setdefault(k, lev)
    order = sorted(first, key=lambda k: (sum(first[k]), first[k]))
    idx = {k: i for i, k in enumerate(order)}
    edges = {}
    for lev, k in cls.items():
        for j in range(len(lev)):
            if lev[j] < nlay[j]:
                l2 = list(lev)
                l2[j] += 1
                a, c = idx[k], idx[cls[tuple(l2)]]
                edges.setdefault((a, c), set()).add(j)
    rec = {"n_states": len(order), "moves": ["{}{}".format(*m) for m in e["moves"]],
           "class_of_levels": {"".join(map(str, l)): idx[k] for l, k in cls.items()},
           "edges": [[a, c, sorted(js)] for (a, c), js in sorted(edges.items()) if a != c],
           "null_moves": [[a, sorted(js)] for (a, c), js in sorted(edges.items()) if a == c]}
    Path(out or f"graph_{child}.json").write_text(json.dumps(rec, indent=0))
    print(len(order), "classes,", len(rec["edges"]), "edges; moves that change nothing:", rec["null_moves"])
    return rec


def unique(folders, out="unique.json"):
    """Unique built dots over several path folders (dot key); unique.json: first (folder, state) and every occurrence."""
    uniq, occ = {}, {}
    for folder in folders:
        for r in json.loads((Path(folder) / "states.json").read_text())["states"]:
            if "formula" not in r:
                continue
            p = Path(folder) / (f"s{r['i']:02d}" if r["i"] < 100 else f"s{r['i']}") / "out.xyz"
            k = L.dot_key(*L.read_xyz(p))
            uniq.setdefault(k, (str(folder), r["i"], r["formula"]))
            occ.setdefault(k, []).append(f"{Path(folder).name}:{r['i']}")
    print(sum(len(v) for v in occ.values()), "states ->", len(uniq), "unique dots")
    for f in folders:
        print(f"  {Path(f).name}: {sum(1 for v in uniq.values() if v[0] == str(f))} new")
    shared = [(uniq[k][2], v) for k, v in occ.items() if len({o.split(':')[0] for o in v}) > 1]
    print(len(shared), "dots shared between paths:")
    for f, v in sorted(shared, key=lambda t: len(t[0])):
        print("  ", f, v)
    rows = [{"folder": v[0], "i": v[1], "formula": v[2], "occurs": occ[k]} for k, v in uniq.items()]
    Path(out).write_text(json.dumps(rows, indent=0))
    return rows
