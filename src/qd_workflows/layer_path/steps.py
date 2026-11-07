"""
What the registry's layer steps add (the scripts layer_steps.py, attach_dg.py, ladders.py, ladder_fit.py, chains.py).

layer_steps: for every dot-level layer edge a -> b of every edges file, a's lattice sites (Cd, and the anion sites, Se
or lattice Cl, species-blind) are embedded into b by a Td op + lattice translation. If they embed, the step is a pure
addition D = b - a. Its direction is the normal (14 {100}/{111}A/{111}B and 12 {110}) along which the new Cd lie
entirely beyond a's outermost Cd, the one whose heights span the least. Sheet = the new Cd in the lowest layer;
perimeter = sheet Cd with fewer in-plane Cd neighbours than a full sheet. Chains a -> b -> c in the same edge are
tested for b - a and c - b in the same slab (two rungs of one layer: same_slab_next).

attach_dg: dG = G_b - G_a from the current g_cache.json, and the parent's Se count m_a.

ladders: one-sheet pure additions from the same a along the same family with different sheet sizes s are rungs of one
island ladder; Delta G(s) = G_b - G_a, Delta Omega(s) = Delta G - dm mu_CdSe - dk mu_CdCl2. Also the per-family
regression Delta G = dm g_m + dk g_k + P eps over every one-sheet step (P the perimeter Cd).

ladder_fit: per ladder and pool, Delta Omega(s) = s delta + gamma_e sqrt(s) (no intercept); barrier -gamma_e^2 /
(4 delta) at s* = (gamma_e / (2 delta))^2, and the ladder's own top (highest rung before the first downhill one).

chains: same-sheet chains, island growth paths inside one layer: from every step followed by a second rung in the same
sheet, the pure steps b -> c -> d ... while every Cd added since a stays in one sheet beyond a and the sheet grows.
"""
from __future__ import annotations

import collections
import json
import re
from multiprocessing import Pool
from pathlib import Path

import numpy as np

from .lattice import NORMALS, SPACING, TD_OPS

# (mu_CdSe, mu_CdCl2) in eV of the three reference pools (alpha = the pool's CdCl2 share), T = 500 K
POOLS = {0.05: (-5.482, -8.220), 0.2: (-5.408, -8.234), 0.5: (-5.308, -8.294)}
KT = 8.617333e-5 * 500
FAMILIES = ("100", "111A", "111B", "110")

_PAT = re.compile(r"\(np\.str_\('(\w+)'\), (-?\d+), (-?\d+), (-?\d+)\)")
_FCC = np.array([v for v in {tuple(s * np.array(p)) for p in ((2, 2, 0), (2, 0, 2), (0, 2, 2)) for s in (1, -1)}] +
                [v for v in {(2, -2, 0), (-2, 2, 0), (2, 0, -2), (-2, 0, 2), (0, 2, -2), (0, -2, 2)}])
DIRS = [(f, np.array(n)) for f in FAMILIES for n in NORMALS[f]]
_CTX = {}


def _init(ctx):
    _CTX.update(ctx)
    if "registry" in ctx and "ID" not in ctx:
        reg = json.loads(Path(ctx["registry"]).read_text())
        _CTX["ID"] = {v["id"]: k for k, v in reg["dots"].items()}


def parse(dk, cl=True):
    """Lattice sites ('C' Cd, 'A' Se or lattice Cl; cl=False: Se only) of a dot key in a/4 units, and the counts of
    Se, Cd, Cl (all)."""
    sp, X = [], []
    cnt = {"Cd": 0, "Se": 0, "Cl": 0}
    for s, x, y, z in _PAT.findall(dk):
        cnt[s] += 1
        v = np.array([int(x), int(y), int(z)])
        if np.all(v % 4 == 0) and (cl or s != "Cl"):
            sp.append("C" if s == "Cd" else "A")
            X.append(v // 4)
    return np.array(sp), np.array(X), cnt


def _code(sp, X):
    return ((np.where(sp == "C", 1, 0).astype(np.int64) * 4096 + X[:, 0] + 2048) * 4096 + X[:, 1] + 2048) * 4096 \
        + X[:, 2] + 2048


def embed(a, b):
    """a's sites put on b's (Td op + translation): (a in b's frame, M, t), or None if a is not a subset of b."""
    (sa, Xa), (sb, Xb) = a, b
    cb = np.sort(_code(sb, Xb))
    anc_b = Xb[sb == "A"]
    for M in TD_OPS:
        P = Xa @ M.T
        T = anc_b - P[sa == "A"][0]
        for i in range(0, len(T), 64):
            Y = P[None] + T[i:i + 64, None]
            cc = _code(np.broadcast_to(sa, Y.shape[:2]).ravel(), Y.reshape(-1, 3)).reshape(Y.shape[:2])
            hit = np.isin(cc, cb).all(axis=1)
            if hit.any():
                j = int(np.argmax(hit))
                return Y[j], M, T[i + j]
    return None


def step(sa, Ya, sb, Xb):
    """Direction, sheet and perimeter of D = b - a (a already in b's frame)."""
    ka = set(map(tuple, np.c_[np.where(sa == "C", 1, 0), Ya]))
    kb = np.c_[np.where(sb == "C", 1, 0), Xb]
    inD = np.array([tuple(r) not in ka for r in kb])
    D, sD = Xb[inD], sb[inD]
    if not len(D):
        return {"nD": 0}
    out = {"nD": int(inD.sum()), "nD_A": int((sD == "A").sum()), "nD_C": int((sD == "C").sum())}
    C, Ca = D[sD == "C"], Ya[sa == "C"]
    if not len(C):
        out["dir"] = None                   # no new Cd: anion / Cl changes only
        return out
    # direction from the new Cd sites (the new Se often takes an anion site that a's Cl cap already held)
    best = None
    for f, n in DIRS:
        h = C @ n
        if h.min() > (Ca @ n).max():
            span = (h.max() - h.min()) / SPACING[f]
            if best is None or span < best[2]:
                best = (f, n, span)
    if best is None:
        out["dir"] = "mixed"                # new Cd on more than one side
        return out
    f, n, span = best
    h = C @ n
    lay = h == h.min()
    sheet = {tuple(v) for v in C[lay]}
    full = 6 if f.startswith("111") else (4 if f == "100" else 2)
    inpl = [sum(tuple(v + d) in sheet for d in _FCC if int(d @ n) == 0) for v in C[lay]]
    out.update({"dir": f, "normal": n.tolist(), "span": float(span), "sheet": int(lay.sum()),
                "perimeter": int(sum(c < full for c in inpl)), "n_sheets": int(len(set(h.tolist())))})
    return out


def _work(fn):
    e = json.loads(Path(fn).read_text())
    if "error" in e:
        return []
    ID = _CTX["ID"]
    cache, rows = {}, []

    def get(i):
        return cache.setdefault(i, parse(ID[i]))
    by_a = {}
    for a, b in e["edges"]:
        if a not in ID or b not in ID:
            continue
        pa, pb = get(a), get(b)
        if pb[2]["Se"] < pa[2]["Se"] or (pb[2]["Se"] == pa[2]["Se"] and pb[2]["Cd"] < pa[2]["Cd"]):
            a, b, pa, pb = b, a, pb, pa
        em = embed(pa[:2], pb[:2])
        Ya = None if em is None else em[0]
        fw = em is not None or embed(parse(ID[a], False)[:2], parse(ID[b], False)[:2]) is not None
        r = {"framework": fw, "edge": Path(fn).stem, "a": a, "b": b, "dSe": pb[2]["Se"] - pa[2]["Se"],
             "dCd": pb[2]["Cd"] - pa[2]["Cd"], "dCl": pb[2]["Cl"] - pa[2]["Cl"], "pure": Ya is not None}
        if Ya is not None:
            r.update(step(pa[0], Ya, pb[0], pb[1]))
            r["_Ya"] = Ya
        rows.append(r)
        by_a.setdefault(a, []).append(r)
    # two consecutive steps a -> b -> c that fill one slab (one sheet beyond a, along some normal) = two rungs of a layer
    for r in rows:
        if not r["pure"] or r.get("dir") in (None, "mixed"):
            continue
        for r2 in by_a.get(r["b"], []):
            if not r2["pure"] or r2.get("dir") in (None, "mixed"):
                continue
            pa, pb, pc = get(r["a"]), get(r["b"]), get(r2["b"])
            em = embed(pb[:2], pc[:2])
            if em is None:
                continue
            Yb, M, t = em
            Ya = r["_Ya"] @ M.T + t
            kb = set(map(tuple, np.c_[np.where(pb[0] == "C", 1, 0), Yb]))
            ka = set(map(tuple, np.c_[np.where(pa[0] == "C", 1, 0), Ya]))
            A1 = np.array([v for s_, v in zip(pb[0], Yb) if s_ == "C" and (1, *v) not in ka])
            A2 = np.array([v for s_, v in zip(pc[0], pc[1]) if s_ == "C" and (1, *v) not in kb])
            if not len(A1) or not len(A2):
                continue
            Aa = Ya[pa[0] == "C"]
            for f, n in DIRS:
                h = np.r_[A1 @ n, A2 @ n]
                if h.min() > (Aa @ n).max() and h.max() - h.min() < SPACING[f]:
                    r.setdefault("same_slab_next", []).append([r2["b"], f, len(A1), len(A2)])
                    break
    for r in rows:
        r.pop("_Ya", None)
    return rows


def g_by_id(reg, results, g_cache):
    """G at T per registry id, from the g_cache entries ("T|key|E") whose E is the current result's."""
    G = {}
    for k, v in g_cache.items():
        _t, dk, E = k.split("|", 2)
        if dk in results and abs(float(E) - results[dk]["E_eV"]) < 1e-9:
            G[reg["dots"][dk]["id"]] = v
    return G


def layer_steps(net, results, g_cache, out, nproc=4):
    net = Path(net)
    reg = json.loads((net / "registry.json").read_text())
    G = g_by_id(reg, json.loads(Path(results).read_text()), json.loads(Path(g_cache).read_text()))
    files = sorted((net / "edges").glob("*.json"))
    with Pool(nproc, initializer=_init, initargs=({"registry": str(net / "registry.json")},)) as pool:
        rows = [r for rr in pool.imap_unordered(_work, files, chunksize=2) for r in rr]
    for r in rows:
        if r["a"] in G and r["b"] in G:
            r["dG"] = G[r["b"]] - G[r["a"]]
    Path(out).write_text(json.dumps(rows))
    pure = [r for r in rows if r["pure"]]
    print(f"{len(rows)} layer steps, {len(pure)} pure additions, {sum(r.get('dir') is not None for r in pure)} with "
          f"one direction, {sum(bool(r.get('same_slab_next')) for r in pure)} with a second rung in the same slab, "
          f"{sum('dG' in r for r in rows)} with both G")
    return rows


def attach_dg(net, results, g_cache, steps_file):
    """Re-attach dG = G_b - G_a (and m_a) to a layer_steps.json in place."""
    reg = json.loads((Path(net) / "registry.json").read_text())
    G = g_by_id(reg, json.loads(Path(results).read_text()), json.loads(Path(g_cache).read_text()))
    mSe = {v["id"]: k.count("'Se'") for k, v in reg["dots"].items()}
    rows = json.loads(Path(steps_file).read_text())
    for r in rows:
        r.pop("dG", None)
        r["m_a"] = mSe.get(r["a"])
        if r["a"] in G and r["b"] in G:
            r["dG"] = G[r["b"]] - G[r["a"]]
    Path(steps_file).write_text(json.dumps(rows))
    print(f"{len(G)} dots with G; {sum('dG' in r for r in rows)} of {len(rows)} steps with both G")
    return rows


def ladders(steps_file, out):
    rows = json.loads(Path(steps_file).read_text())
    pure = [r for r in rows if r["pure"]]
    one = [r for r in pure if r.get("dir") in FAMILIES and r.get("n_sheets") == 1]
    print(f"{len(rows)} layer steps; {len(pure)} pure additions ({sum(r['framework'] for r in rows)} with the Cd+Se "
          f"frameworks overlaying); {len(one)} add one Cd sheet")
    print("direction of the pure additions:", dict(collections.Counter(r.get("dir") for r in pure)))
    for f in FAMILIES:
        S = [r for r in one if r["dir"] == f]
        if S:
            sz = collections.Counter(r["sheet"] for r in S)
            print(f"  {f}: {len(S)} one-sheet steps, sheet sizes {dict(sorted(sz.items()))}, "
                  f"{sum(bool(r.get('same_slab_next')) for r in S)} followed by a second rung in the same sheet")
    lad = collections.defaultdict(dict)          # same a, same family, several sheet sizes
    for r in one:
        if "dG" not in r:
            continue
        k, s = (r["a"], r["dir"]), r["sheet"]
        cur = lad[k].get(s)
        if cur is None or r["dG"] < cur["dG"]:
            lad[k][s] = {"b": r["b"], "dG": r["dG"], "dm": r["dSe"], "dk": r["dCd"] - r["dSe"], "P": r["perimeter"],
                         "m_a": r.get("m_a")}
    lads = {f"{a}|{f}": v for (a, f), v in lad.items() if len(v) >= 3}
    print(f"{len(lads)} ladders with >= 3 relaxed sheet sizes from one parent state:",
          dict(collections.Counter(k.split("|")[1] for k in lads)))
    fit = {}
    for f in FAMILIES:
        S = [r for r in one if r["dir"] == f and "dG" in r]
        if len(S) < 10:
            continue
        X = np.array([[r["dSe"], r["dCd"] - r["dSe"], r["perimeter"]] for r in S], float)
        y = np.array([r["dG"] for r in S])
        c, *_ = np.linalg.lstsq(X, y, rcond=None)
        fit[f] = {"n": len(S), "g_m": c[0], "g_k": c[1], "eps": c[2], "rmse": float(np.sqrt(((y - X @ c) ** 2).mean()))}
        print(f"  fit {f}: n = {len(S)}, g_m = {c[0]:.3f}, g_k = {c[1]:.3f}, eps = {c[2]:.3f} eV/perimeter Cd, "
              f"rmse {fit[f]['rmse']:.3f} eV")
    res = {"ladders": {k: {str(s): v for s, v in d.items()} for k, d in lads.items()}, "fit": fit}
    Path(out).write_text(json.dumps(res))
    return res


def ladder_fit(ladders_file, out, pools=POOLS):
    lad = json.loads(Path(ladders_file).read_text())["ladders"]
    res = {}
    for al, (mm, mk) in pools.items():
        fam = {}
        for key, d in lad.items():
            f = key.split("|")[1]
            s = np.array(sorted(int(x) for x in d), float)
            om = np.array([d[str(int(x))]["dG"] - d[str(int(x))]["dm"] * mm - d[str(int(x))]["dk"] * mk for x in s])
            X = np.c_[s, np.sqrt(s)]
            (dl, ge), *_ = np.linalg.lstsq(X, om, rcond=None)
            rms = float(np.sqrt(((om - X @ [dl, ge]) ** 2).mean()))
            dn = np.nonzero(om < 0)[0]
            top = float(max(0.0, om[:dn[0]].max())) if len(dn) and dn[0] > 0 else (0.0 if len(dn) else float(om.max()))
            bar = float(-ge ** 2 / (4 * dl)) if dl < 0 and ge > 0 else (0.0 if dl < 0 else np.inf)
            fam.setdefault(f, []).append({"a": key.split("|")[0], "m_a": next(iter(d.values())).get("m_a"), "n": len(s),
                                          "delta": float(dl), "gamma_e": float(ge), "rms": rms, "barrier_fit": bar,
                                          "top": top,
                                          "s_star": float((ge / (2 * dl)) ** 2) if dl < 0 and ge > 0 else None})
        res[al] = fam
        print(f"alpha = {al}: mu_CdSe {mm}, mu_CdCl2 {mk}")
        for f, Lf in fam.items():
            if len(Lf) < 5:
                print(f"  {f}: {len(Lf)} ladders (too few)")
                continue

            def g(k):
                return np.array([x[k] for x in Lf])
            dl, ge, bar, top, rms = g("delta"), g("gamma_e"), g("barrier_fit"), g("top"), g("rms")
            print(f"  {f}: {len(Lf)} ladders, downhill per Cd {(dl < 0).mean():.0%}; median delta {np.median(dl):.3f} "
                  f"eV/Cd, gamma_e {np.median(ge):.3f} eV/sqrt(Cd), fit rms {np.median(rms):.2f} eV")
            fb = bar[np.isfinite(bar)]
            print(f"      barrier (fit) median {np.median(fb):.2f} eV [25-75%: {np.percentile(fb, 25):.2f}-"
                  f"{np.percentile(fb, 75):.2f}], ladder top median {np.median(top):.2f} eV; top < 3 kT in "
                  f"{(top < 3 * KT).mean():.0%}")
            ma = np.array([x["m_a"] or 0 for x in Lf])
            for lo, hi in ((0, 30), (30, 60), (60, 100), (100, 150), (150, 400)):
                sel = (ma >= lo) & (ma < hi) & np.isfinite(bar)
                if sel.sum() >= 5:
                    print(f"      m_a {lo}-{hi}: {sel.sum()} ladders, delta {np.median(dl[sel]):.3f}, gamma_e "
                          f"{np.median(ge[sel]):.3f}, barrier {np.median(bar[sel]):.2f} eV, top {np.median(top[sel]):.2f} eV")
    out_ = {str(k): v for k, v in res.items()}
    Path(out).write_text(json.dumps(out_))
    return out_


def _trace(job):
    edge, starts, steps, maxd = job
    ID = _CTX["ID"]
    cache, emc = {}, {}

    def get(i):
        return cache.setdefault(i, parse(ID[i]))
    succ = {}
    for a, b in steps:
        succ.setdefault(a, []).append(b)

    def emb(u, v):
        if (u, v) not in emc:
            emc[(u, v)] = embed(get(u)[:2], get(v)[:2])
        return emc[(u, v)]

    def sheet(sa, Ya, sx, Xx, fam):
        ka = set(map(tuple, Ya[sa == "C"]))
        C = np.array([v for s_, v in zip(sx, Xx) if s_ == "C" and tuple(v) not in ka])
        if not len(C):
            return None
        Aa = Ya[sa == "C"]
        for f, n in DIRS:
            if f == fam:
                h = C @ n
                if h.min() > (Aa @ n).max() and h.max() - h.min() < SPACING[f]:
                    return len(C)
        return None

    out = []
    for a, fam in starts:
        pa = get(a)
        nodes, edges, stack = {a: 0}, set(), []
        for b in succ.get(a, []):
            em = emb(a, b)
            if em is None:
                continue
            s = sheet(pa[0], em[0], *get(b)[:2], fam)
            if s:
                stack.append((b, em[0], s, 1))
                edges.add((a, b))
                nodes[b] = s
        seen = set()
        while stack:
            x, Ya, s, d = stack.pop()
            if (x, d) in seen or d >= maxd:
                continue
            seen.add((x, d))
            for y in succ.get(x, []):
                em = emb(x, y)
                if em is None:
                    continue
                _, M, t = em
                Yy = Ya @ M.T + t
                s2 = sheet(pa[0], Yy, *get(y)[:2], fam)
                if s2 and s2 > s:
                    edges.add((x, y))
                    nodes[y] = max(nodes.get(y, 0), s2)
                    stack.append((y, Yy, s2, d + 1))
        out.append({"edge": edge, "a": a, "fam": fam, "nodes": nodes, "edges": sorted(edges)})
    return out


def chains(net, steps_file, out, nproc=4, max_depth=12):
    rows = json.loads(Path(steps_file).read_text())
    by_edge = {}
    for r in rows:
        if r["pure"] and r.get("dir") not in (None, "mixed"):
            by_edge.setdefault(r["edge"], {"steps": [], "starts": set()})["steps"].append((r["a"], r["b"]))
        if r.get("same_slab_next"):
            by_edge[r["edge"]]["starts"].add((r["a"], r["dir"]))
    jobs = [(e, sorted(v["starts"]), v["steps"], max_depth) for e, v in by_edge.items() if v["starts"]]
    with Pool(nproc, initializer=_init, initargs=({"registry": str(Path(net) / "registry.json")},)) as pool:
        res = [c for cc in pool.imap_unordered(_trace, jobs) for c in cc]
    for c in res:
        c["nodes"] = {str(k): v for k, v in c["nodes"].items()}
    Path(out).write_text(json.dumps(res))
    depth = [max(c["nodes"].values()) for c in res]
    print(f"{len(res)} chain graphs from {len(jobs)} lineage edges; {sum(len(c['nodes']) for c in res)} nodes, "
          f"{sum(len(c['edges']) for c in res)} steps; largest sheet median {np.median(depth) if depth else 0:.0f}, "
          f"max {max(depth) if depth else 0}")
    return res
