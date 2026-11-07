"""
The m 20-50 fill: one-layer neighbours of the network's low-lying states, built and put in the registry (the scripts
fill2050/enumerate.py, links.py, build_all.py, register.py).

1. enumerate (no builds): parents = per m in [m_lo, m_hi] the states within dw eV of that m's lowest Omega at the
   given pools (at most cap per m). Each parent is a lattice cut by the 14 {100}/{111}A/{111}B planes; every plane
   moved one spacing out (a new sheet) or in (a sheet removed) gives a neighbour, kept when its class key is not in
   the registry yet. candidates.json: parents, and new class key -> parents, moves, builder planes, Cd and A counts.
2. links: per candidate, the registry dots and other candidates one layer away (links.json).
3. build_all (builds, no relaxation): the candidates (builder's own k0) and the shells k0 +- 1, +- 2 (z_type_target) of
   every new framework and of the low-lying parents with one k only; deduplicated against the registry and among
   themselves by dot key. OUTDIR/builds/<name>/ and OUTDIR/fill.json.
4. register: the kept dots (a list of fill names, after the relaxation wave) into the registry (builds copied to
   net/builds/f<name>), their lattice classes, and one edges file net/edges/fill2050.json with the dot-level steps:
   each candidate's parents, the registry dots and other candidates one layer away.
"""
from __future__ import annotations

import json
import re
import shutil
import time
from collections import Counter, defaultdict
from multiprocessing import Pool
from pathlib import Path

import numpy as np

from . import lattice as L
from .build import build

SHELLS = (+1, -1, +2, -2)
_CTX = {}


def _init(ctx):
    _CTX.update(ctx)


def _neighbours(G):
    """The lattice box of a class's sites, its 14 levels and cut mask, and every one-plane move (sign, family,
    normal, levels, mask) that changes the cut and keeps an anion."""
    LG, Ls = L.lattice_box(G)
    lev = {n: int((G @ np.array(n)).max()) for _f, n in L.PLANES}
    m0 = L.inside(LG, list(lev.items()))
    moves = []
    for f, n in L.PLANES:
        for sgn in (1, -1):
            pl = dict(lev)
            pl[n] += sgn * L.SPACING[f]
            msk = L.inside(LG, list(pl.items()))
            if (Ls[msk] == "A").any() and msk.sum() != m0.sum():
                moves.append((sgn, f, n, pl, msk))
    return LG, Ls, lev, m0, moves


def low_lying(network, thermo, params, m_lo=14, m_hi=55, dw=0.4, cap=6, T=500.0, eps=2.0, c_m=(7.6e-4, 3e-3),
              c_x=0.03 - 0.0038):
    """Registry id -> (m, k, label) of the states within dw of their m's lowest Omega at the pools c_m."""
    from .. import network_growth as ng
    net = ng.load(network, thermo, params)
    g, th, _ = net.G(T, eps)
    par = {}
    for cm in c_m:
        pl = ng._pool(th, cm, c_x)
        om = g - net.m * pl["mu_MA"] - net.k * pl["mu_MX"]
        om = om - om[net.seed]
        for m in range(m_lo, m_hi + 1):
            s = np.where((net.m == m) & (net.reg_id >= 0))[0]
            if not len(s):
                continue
            s = s[np.argsort(om[s])]
            for i in s[:cap]:
                if om[i] <= om[s[0]] + dw:
                    par[int(net.reg_id[i])] = (int(net.m[i]), int(net.k[i]), net.label[i])
    return par


def enumerate_candidates(registry, parents, out):
    """One-layer neighbours (new class keys) of the parents {registry id: (m, k, label)}; writes candidates.json."""
    reg = json.loads(Path(registry).read_text())
    print(f"{len(parents)} parents", flush=True)
    cid_key = {v: k for k, v in reg["classes"].items()}
    cls_of_dot = defaultdict(list)
    for cid, did in reg["class_dot"].items():
        if did is not None:
            cls_of_dot[int(did)].append(cid_key[int(cid)])
    known = set(reg["classes"])
    cand, nopar, notcut = {}, 0, 0
    for did, (m, k, lab) in sorted(parents.items(), key=lambda kv: kv[1]):
        keys = cls_of_dot.get(did)
        if not keys:
            nopar += 1                                  # a lineage node's own dot: no lattice class
            continue
        key = eval(keys[0], {"np": np})
        G = np.array([x[1:] for x in key], int)
        LG, Ls, _lev, m0, moves = _neighbours(G)
        if L.canon(Ls[m0], LG[m0]) != key:
            notcut += 1                                 # not a 14-plane cut (should not happen)
            continue
        O = L.central_anion(LG[m0], Ls[m0])
        for sgn, f, n, pl, msk in moves:
            rk = repr(L.canon(Ls[msk], LG[msk]))
            if rk in known:
                continue
            c = cand.setdefault(rk, {"parents": [], "moves": [], "n_Cd": int((Ls[msk] == "Cd").sum()),
                                     "n_A": int((Ls[msk] == "A").sum())})
            c["parents"].append([did, m, k, lab])
            c["moves"].append(f"{'+' if sgn > 0 else '-'}{f}{tuple(n)}")
            if "planes" not in c:
                c["planes"] = [[list(n2), h] for n2, h in L.planes_A(LG, pl, O)[0]]
    print(f"skipped: {nopar} parents without a class, {notcut} not 14-plane cuts")
    add = sum(1 for c in cand.values() if any(mv[0] == "+" for mv in c["moves"]))
    print(f"{len(cand)} new classes ({add} by adding a sheet, {len(cand) - add} by removing one)")
    hist = defaultdict(int)
    for c in cand.values():
        hist[(c["n_A"] // 10) * 10] += 1
    print("anion sites (Se + lattice Cl):", dict(sorted(hist.items())))
    res = {"parents": {str(k): list(v) for k, v in parents.items()}, "candidates": cand}
    Path(out).write_text(json.dumps(res))
    return res


def links(candidates, registry, out):
    d = json.loads(Path(candidates).read_text())
    reg = json.loads(Path(registry).read_text())
    cdot = {int(k): v for k, v in reg["class_dot"].items()}
    known = {k: cdot.get(v) for k, v in reg["classes"].items()}
    cand = d["candidates"]
    par = {int(k) for k in d["parents"]}
    nreg, npar, ncand = Counter(), Counter(), Counter()
    res = {}
    for rk in cand:
        _sym, G = L.key_sites(rk)
        LG, Ls, _lev, _m0, moves = _neighbours(G)
        r_, c_ = set(), set()
        for _sgn, _f, _n, _pl, msk in moves:
            k2 = repr(L.canon(Ls[msk], LG[msk]))
            if k2 in known and known[k2] is not None:
                r_.add(int(known[k2]))
            elif k2 in cand and k2 != rk:
                c_.add(k2)
        res[rk] = {"registry": sorted(r_), "candidates": sorted(c_)}
        nreg[min(len(r_), 6)] += 1
        npar[len(r_ & par)] += 1
        ncand[min(len(c_), 6)] += 1
    print("registry dots one layer away (6 = 6+):", dict(sorted(nreg.items())))
    print("  of them low-lying parents:", dict(sorted(npar.items())))
    print("other candidates one layer away:", dict(sorted(ncand.items())))
    print("candidates linking >= 2 low-lying parents:", sum(1 for v in res.values() if len(set(v["registry"]) & par) >= 2))
    Path(out).write_text(json.dumps(res))
    return res


def _one(task):
    name, planes, k_total = task
    wd = Path(_CTX["out"]) / "builds" / name
    try:
        ok = build(planes, wd, _CTX["cif"], k_total)
    except Exception as exc:                                  # the builder refuses some cuts
        return name, None, str(exc)[:80]
    if not ok:
        return name, None, "build failed"
    sy, x = L.read_xyz(wd / "out.xyz")
    return name, (repr(L.dot_key(sy, x)), L.formula(sy), repr(L.se_key(sy, x)), len(sy)), None


def _k_of(formula):
    return int(re.search(r"Cl(\d+)", formula).group(1)) // 2


def build_all(candidates, registry, out, cif, nproc=4):
    t0 = time.time()
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    CAND = json.loads(Path(candidates).read_text())
    reg = json.loads(Path(registry).read_text())
    se_k = {}
    for v in reg["dots"].values():
        se_k.setdefault(v["se"], set()).add(_k_of(v["formula"]))
    cand = CAND["candidates"]
    names = {f"n{j:04d}": rk for j, rk in enumerate(sorted(cand))}
    with Pool(nproc, initializer=_init, initargs=({"out": str(out), "cif": str(cif)},)) as pool:
        res = pool.map(_one, [(nm, [(tuple(n), h) for n, h in cand[rk]["planes"]], None) for nm, rk in names.items()],
                       chunksize=4)
        seen, dots, fail = set(reg["dots"]), {}, {}
        for nm, r, err in res:
            if r is None:
                fail[nm] = err
                continue
            if r[0] in seen:
                continue
            seen.add(r[0])
            dots[nm] = {"key": r[0], "formula": r[1], "se": r[2], "atoms": r[3], "kind": "layer", "class": names[nm],
                        "parents": cand[names[nm]]["parents"], "moves": cand[names[nm]]["moves"]}
        print(f"bases: {len(names)} built, {len(dots)} new dots, {len(fail)} failed, {time.time() - t0:.0f} s", flush=True)
        # shells: every new framework, and the low-lying parents whose framework has one k only
        sh = []
        for nm, v in list(dots.items()):
            k0 = _k_of(v["formula"])
            sh += [(f"{nm}_k{dk:+d}", [(tuple(n), h) for n, h in cand[v["class"]]["planes"]], k0 + dk)
                   for dk in SHELLS if k0 + dk >= 0]
        reg_by_id = {v["id"]: v for v in reg["dots"].values()}
        cdot = {int(c): d for c, d in reg["class_dot"].items()}
        cls_of = {}
        for rk, cid in reg["classes"].items():
            cls_of.setdefault(cdot.get(cid), rk)
        for did in CAND["parents"]:
            v = reg_by_id[int(did)]
            if len(se_k[v["se"]]) > 1 or int(did) not in cls_of:
                continue
            k0 = _k_of(v["formula"])
            pl = L.class_planes(cls_of[int(did)])
            sh += [(f"p{int(did):05d}_k{dk:+d}", pl, k0 + dk) for dk in SHELLS if k0 + dk >= 0]
        res = pool.map(_one, sh, chunksize=4)
    ns = 0
    for nm, r, err in res:
        if r is None:
            fail[nm] = err
            continue
        if r[0] in seen:
            continue
        seen.add(r[0])
        ns += 1
        dots[nm] = {"key": r[0], "formula": r[1], "se": r[2], "atoms": r[3], "kind": "shell"}
    print(f"shells: {len(sh)} built, {ns} new dots; {len(dots)} new dots in all; {len(fail)} failed; "
          f"{time.time() - t0:.0f} s", flush=True)
    rec = {"dots": dots, "failed": fail}
    (out / "fill.json").write_text(json.dumps(rec))
    return rec


def register(net, fill_dir, candidates, links_file, keep, edges_name="fill2050", label="m 20-50 one-layer neighbours"):
    """Put the kept fill dots (names) in the registry and write net/edges/<edges_name>.json."""
    net, fill_dir = Path(net), Path(fill_dir)
    reg = json.loads((net / "registry.json").read_text())
    cand = json.loads(Path(candidates).read_text())["candidates"]
    fill = json.loads((fill_dir / "fill.json").read_text())["dots"]
    lk = json.loads(Path(links_file).read_text())
    names = {f"n{j:04d}": rk for j, rk in enumerate(sorted(cand))}
    ef = net / "edges" / f"{edges_name}.json"
    if ef.exists():
        raise FileExistsError(f"{ef} exists: this fill is registered already")
    n0 = len(reg["dots"])
    for nm in sorted(keep):
        v = fill[nm]
        if v["key"] in reg["dots"]:
            continue
        dst = net / "builds" / f"f{nm}"
        shutil.copytree(fill_dir / "builds" / nm, dst, dirs_exist_ok=True)
        reg["dots"][v["key"]] = {"id": len(reg["dots"]), "formula": v["formula"], "se": v["se"], "dir": dst.name}
    # class -> dot for every candidate (also those built into an existing dot, or below m 14 and not kept)
    dot_of_class = {}
    for nm, rk in names.items():
        p = fill_dir / "builds" / nm / "out.xyz"
        if not p.exists():
            continue
        dk = repr(L.dot_key(*L.read_xyz(p)))
        if dk in reg["dots"]:
            dot_of_class[rk] = reg["dots"][dk]["id"]
            if rk not in reg["classes"]:
                reg["classes"][rk] = len(reg["classes"])
                reg["class_dot"][str(reg["classes"][rk])] = dot_of_class[rk]
    edges = set()
    for rk, c in cand.items():
        d = dot_of_class.get(rk)
        if d is None:
            continue
        nb = {pp[0] for pp in c["parents"]} | set(lk[rk]["registry"])
        nb |= {dot_of_class[k2] for k2 in lk[rk]["candidates"] if k2 in dot_of_class}
        for o in nb:
            if o != d:
                edges.add((min(o, d), max(o, d)))
    ef.parent.mkdir(parents=True, exist_ok=True)
    ef.write_text(json.dumps({"p": -1, "c": -1, "fill": label, "edges": sorted(map(list, edges))}))
    (net / "registry.json").write_text(json.dumps(reg))
    print(f"registry: {n0} -> {len(reg['dots'])} dots; {len(dot_of_class)} candidate classes on dots "
          f"({sum(1 for d in dot_of_class.values() if d < n0)} already in the registry); {len(edges)} layer edges")
    return reg
