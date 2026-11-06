"""
Layer paths of every lineage parent -> child edge, all states in one registry (the scripts size_edges.py and
treepaths.py).

size_edges: per lineage edge, the parent's lattice sites embedded in the child's (Td op + lattice translation), the
{100}/{111} plane moves and the number of layer states (product of L + 1 over the moved planes): edges_size.json.

treepaths: for each edge with child m <= mmax and at most smax layer states, the planes go from the parent's outermost
sites to the child's one layer at a time; every layer state is the lattice cut by its planes. States are identified
by their site set up to Td ops + translation (class key); each new class is built once with QD_builder (build.build).
The built dot gets an id in net/registry.json (two classes can give the same dot after the prepass). Edge files
net/edges/<p>_<c>.json: the dots of the path and the dot-level layer edges. An edge whose path end misses the child's
Cd sites is retried with the 12 {110} planes as well. Run in waves of growing mmax; a rerun skips the edges done.

registry.json: classes (repr class key -> class id), dots (repr dot key -> id, formula, se, dir), class_dot
(class id -> dot id, None when the build failed), lineage (node -> its own scan dot, built under builds/L<node>).
"""
from __future__ import annotations

import itertools
import json
from multiprocessing import Pool
from pathlib import Path

import numpy as np

from . import lattice as L
from .build import build

_CTX = {}


def _init(ctx):
    """Pool initializer: the lineage (loaded once per worker), the net folder, the CIF and the state cap."""
    _CTX.update(ctx)
    if "data" in ctx and "S" not in ctx:
        d = json.loads(Path(ctx["data"]).read_text())
        _CTX.update(S=d["structs"], N=d["lineage"]["nodes"])


def site_arr(i):
    st = _CTX["S"][_CTX["N"][i]["dot"]]
    if not any(a[0] == "Se" for a in st["xyz"]):
        return None
    return L.arr(L.from_struct(st))


def size_edge(pc):
    p, c = pc
    N = _CTX["N"]
    (ps, pG), (cs, cG) = site_arr(p), site_arr(c)
    PG = L.embed_sites(ps, pG, cs, cG)
    if PG is None:
        return {"p": p, "c": c, "parent": N[p]["label"], "child": N[c]["label"], "error": "no embedding on lattice sites"}
    dp = {n: int((PG @ np.array(n)).max()) for _f, n in L.PLANES}
    dc = {n: int((cG @ np.array(n)).max()) for _f, n in L.PLANES}
    lay = [(dc[n] - dp[n]) // L.SPACING[f] for f, n in L.PLANES if dc[n] > dp[n]]
    return {"p": p, "c": c, "parent": N[p]["label"], "child": N[c]["label"], "mp": N[p]["size"], "mc": N[c]["size"],
            "n_moves": len(lay), "layers": int(sum(lay)), "states": int(np.prod([x + 1 for x in lay])),
            "neg": sum(1 for _f, n in L.PLANES if dc[n] < dp[n])}


def size_edges(data, out="edges_size.json", nproc=4):
    """Every lineage edge's embedding, plane moves and number of layer states."""
    N = json.loads(Path(data).read_text())["lineage"]["nodes"]
    pcs = [(p, c) for c, n in enumerate(N) for p in n["parents"]]
    with Pool(nproc, initializer=_init, initargs=({"data": str(data)},)) as pool:
        rows = pool.map(size_edge, pcs, chunksize=8)
    Path(out).write_text(json.dumps(rows, indent=0))
    ok = [r for r in rows if "error" not in r]
    print(len(rows), "edges,", len(rows) - len(ok), "without lattice-site embedding")
    for lim in (200, 1000, 5000, 20000, 100000):
        sel = [r for r in ok if r["states"] <= lim]
        print(f"  states <= {lim}: {len(sel)} edges, sum states {sum(r['states'] for r in sel)}, "
              f"max child m {max([r['mc'] for r in sel] or [0])}")
    print("negative moves (child plane inside parent):", sum(1 for r in ok if r["neg"]))
    return rows


def classes(task):
    """Layer states of one edge: class key -> builder planes and site count, the class edges one layer apart, the
    classes of the path's two ends, and whether those ends have the parent's and the child's Cd sites."""
    p, c, npl = task
    S, N, smax = _CTX["S"], _CTX["N"], _CTX["smax"]
    PL = L.PLANES if npl == 14 else L.PLANES26
    ps, pG = site_arr(p)
    cs, cG = site_arr(c)
    PG = L.embed_sites(ps, pG, cs, cG)
    mode = "all sites"
    if PG is None:                     # fallback: the parent's Cd + Se only (its lattice Cl may sit where the child has Cd)
        st = S[N[p]["dot"]]
        sym = [a[0] for a in st["xyz"]]
        P = np.array([a[2:5] for a in st["xyz"]], float)
        keep = {k for k, v in L.sites(np.array(sym), P, P[sym.index("Se")]).items() if v != "Cl"}
        use = np.array([(s, *g) in keep for s, g in zip(ps, map(tuple, pG))])
        PG = L.embed_sites(ps, pG, cs, cG, use)
        mode = "Cd+Se"
        if PG is None:
            return {"p": p, "c": c, "npl": npl, "error": "no embedding"}
    dp = {n: int((PG @ np.array(n)).max()) for _f, n in PL}
    dc = {n: int((cG @ np.array(n)).max()) for _f, n in PL}
    if any(dc[n] < dp[n] for _f, n in PL):
        return {"p": p, "c": c, "npl": npl, "error": "a child plane lies inside the parent"}
    moves = [(f, n) for f, n in PL if dc[n] > dp[n]]
    nlay = [(dc[n] - dp[n]) // L.SPACING[f] for f, n in moves]
    nst = int(np.prod([k + 1 for k in nlay]))
    if nst > smax:
        return {"p": p, "c": c, "npl": npl, "error": f"{nst} layer states > {smax}"}
    O = L.central_anion(PG, ps)
    LG, Ls = L.lattice_box(cG)
    cls, lev_key = {}, {}
    for lev in itertools.product(*[range(k + 1) for k in nlay]):
        pl = {n: dp[n] for _f, n in PL}
        for (f, n), l in zip(moves, lev):
            pl[n] += l * L.SPACING[f]
        m = L.inside(LG, list(pl.items()))
        if not m.any() or not (Ls[m] == "A").any():
            lev_key[lev] = None
            continue
        k = L.canon(Ls[m], LG[m])
        lev_key[lev] = k
        if k not in cls:
            pA, _active = L.planes_A(LG, pl, O, PL)
            cls[k] = {"planes": [(list(n), h) for n, h in pA], "sites": int(m.sum())}
    edges = set()
    for lev, k in lev_key.items():
        for j in range(len(lev)):
            if lev[j] < nlay[j]:
                l2 = list(lev)
                l2[j] += 1
                k2 = lev_key[tuple(l2)]
                if k is not None and k2 is not None and k != k2:
                    edges.add((k, k2))

    def cd_set(mask):
        return {tuple(g) for g in LG[mask & (Ls == "Cd")]}
    cd_first = cd_set(L.inside(LG, [(n, dp[n]) for _f, n in PL]))
    cd_last = cd_set(L.inside(LG, [(n, dc[n]) for _f, n in PL]))
    return {"p": p, "c": c, "npl": npl, "start_ok": cd_first == {tuple(g) for g in PG[ps == "Cd"]},
            "end_ok": cd_last == {tuple(g) for g in cG[cs == "Cd"]}, "mode": mode,
            "moves": ["{}{}".format(*m_) for m_ in moves], "nlay": [int(x) for x in nlay], "n_states": nst,
            "classes": cls, "edges": list(edges), "first": lev_key[tuple([0] * len(nlay))], "last": lev_key[tuple(nlay)]}


def build_one(args):
    key_id, planes = args
    wd = Path(_CTX["net"]) / "builds" / f"c{key_id:06d}"
    if not build([(tuple(n), h) for n, h in planes], wd, _CTX["cif"]):
        return key_id, None
    sy, x = L.read_xyz(wd / "out.xyz")
    return key_id, (repr(L.dot_key(sy, x)), L.formula(sy), repr(L.se_key(sy, x)))


def new_registry():
    return {"classes": {}, "dots": {}, "class_dot": {}, "lineage": {}}


def lineage_dot(reg, i, net, S, N):
    """The lineage node's own (scan) dot as a registry dot; its xyz under builds/L<i>/out.xyz."""
    if str(i) in reg["lineage"]:
        return reg["lineage"][str(i)]
    st = S[N[i]["dot"]]
    sy = np.array([a[0] for a in st["xyz"]])
    x = np.array([a[2:5] for a in st["xyz"]], float)
    dk = repr(L.dot_key(sy, x))
    if dk not in reg["dots"]:
        wd = Path(net) / "builds" / f"L{i:04d}"
        wd.mkdir(parents=True, exist_ok=True)
        L.write_xyz(wd / "out.xyz", sy, x, f"lineage node {i} {N[i]['label']}")
        reg["dots"][dk] = {"id": len(reg["dots"]), "formula": L.formula(sy), "se": repr(L.se_key(sy, x)), "dir": wd.name}
    reg["lineage"][str(i)] = {"dot": reg["dots"][dk]["id"], "se": reg["dots"][dk]["se"]}
    return reg["lineage"][str(i)]


def treepaths(data, net, mmax, smax, cif, edges_size="edges_size.json", nproc=4, chunk=8):
    """Build the layer paths of the lineage edges with child m <= mmax and <= smax layer states into net/."""
    net = Path(net)
    d = json.loads(Path(data).read_text())
    S, N = d["structs"], d["lineage"]["nodes"]
    (net / "edges").mkdir(parents=True, exist_ok=True)
    (net / "builds").mkdir(exist_ok=True)
    regp = net / "registry.json"
    reg = json.loads(regp.read_text()) if regp.exists() else new_registry()
    sz = {(r["p"], r["c"]): r for r in json.loads(Path(edges_size).read_text())}
    todo = [(p, c) for c, n in enumerate(N) for p in n["parents"]
            if N[c]["size"] <= mmax and not (net / "edges" / f"{p}_{c}.json").exists()
            and (p, c) in sz and ("error" in sz[(p, c)] or sz[(p, c)]["states"] <= smax)]
    todo.sort(key=lambda pc: (N[pc[1]]["size"], sz[pc].get("states", 0)))
    todo = [(p, c, 14) for p, c in todo]
    print(f"{len(todo)} edges to do (child m <= {mmax}, layer states <= {smax})", flush=True)
    ctx = {"data": str(data), "net": str(net), "cif": str(cif), "smax": smax}
    with Pool(nproc, initializer=_init, initargs=(ctx,)) as pool:
        while todo:
            batch, todo = todo[:chunk], todo[chunk:]
            res = pool.map(classes, batch)
            new = []
            for r in res:
                if "error" in r:
                    continue
                for k, v in r["classes"].items():
                    rk = repr(k)
                    if rk not in reg["classes"]:
                        reg["classes"][rk] = len(reg["classes"])
                        new.append((reg["classes"][rk], v["planes"]))
            for cid, out in pool.map(build_one, new):
                if out is None:
                    reg["class_dot"][str(cid)] = None
                    continue
                dk, fml, sk = out
                if dk not in reg["dots"]:
                    reg["dots"][dk] = {"id": len(reg["dots"]), "formula": fml, "se": sk, "dir": f"c{cid:06d}"}
                reg["class_dot"][str(cid)] = reg["dots"][dk]["id"]
            for r in res:
                p, c = r["p"], r["c"]
                if "error" in r:
                    (net / "edges" / f"{p}_{c}.json").write_text(json.dumps(r))
                    print(f"  [{r['npl']}] {N[p]['label']} -> {N[c]['label']}: {r['error']}", flush=True)
                    continue

                def dot(k):
                    return reg["class_dot"].get(str(reg["classes"][repr(k)]))
                de = sorted({(dot(a), dot(c_)) for a, c_ in r["edges"]
                             if dot(a) is not None and dot(c_) is not None and dot(a) != dot(c_)})
                dots = sorted({dot(k) for k in r["classes"] if dot(k) is not None})
                lp, lc = lineage_dot(reg, p, net, S, N), lineage_dot(reg, c, net, S, N)
                d0, d1 = dot(r["first"]), dot(r["last"])
                okp, okc = r["start_ok"], r["end_ok"]
                if not okc and r["npl"] == 14:
                    todo.append((p, c, 26))
                    print(f"  [14] {N[p]['label']} -> {N[c]['label']}: path end misses the child's Cd sites; "
                          f"retry with the {{110}} planes", flush=True)
                    continue
                # the lineage dots join the path ends they share the Cd framework with (isomer / decoration links)
                if okp and d0 is not None and lp["dot"] != d0:
                    de.append((lp["dot"], d0))
                if okc and d1 is not None and lc["dot"] != d1:
                    de.append((d1, lc["dot"]))
                rec = {"p": p, "c": c, "npl": r["npl"], "parent": N[p]["label"], "child": N[c]["label"],
                       "mode": r["mode"], "moves": r["moves"], "nlay": r["nlay"], "n_states": r["n_states"],
                       "n_classes": len(r["classes"]), "dots": dots, "edges": de, "path_start": d0, "path_end": d1,
                       "lineage_parent_dot": lp["dot"], "lineage_child_dot": lc["dot"],
                       "start_is_parent_cd": okp, "end_is_child_cd": okc}
                (net / "edges" / f"{p}_{c}.json").write_text(json.dumps(rec))
                print(f"  [{r['npl']}] {rec['parent']} -> {rec['child']}: Cd frameworks match {okp}/{okc}, "
                      f"{r['n_states']} states, {len(r['classes'])} classes, {len(dots)} dots", flush=True)
            regp.write_text(json.dumps(reg))
            print(f"[registry] {len(reg['classes'])} classes, {len(reg['dots'])} dots", flush=True)
    return reg
