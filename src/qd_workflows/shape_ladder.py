"""
Sampled CdCl2 ladders on shape-library dots (qd_workflows.shapes): binding energy by site class, and full stripping
to the bare core on a few dots.

Site roles. Each dot carries the outward normals and family labels of its facets (from the builder's construction
planes, `planes_for`). On a given geometry every plane is fitted to the atoms: d_i = max over core atoms (Cd, Se) of
n_i . r, and the polyhedron of the fitted planes gives the exposed facets. A removable unit is one Cd with two Cl
(its own Cl first, completed by the nearest Cl): one neutral CdCl2, so every structure stays neutral and the bare
end point is Cd_m Se_m. The unit's role comes from the exposed planes (facet area >= MIN_AREA) its Cd lies within
TOL of: one plane = facet <family>, two = edge <family|family>, three or more = vertex, none = subsurface.
A removal must leave no Cl without a Cd and must keep the structure in one piece.

Steps (all dots batched together in each relaxation call, MACE-MH-1 through orchestr_ai's RelaxCache):
  relax     every dot, polished;
  binding   from each relaxed dot, `sites` units of every role removed one at a time: Delta E_bind(unit) =
            E(k-1) + E(CdCl2) - E(k)  (positive = bound), all polished;
  adsorb    `sites` Z-type CdCl2 additions of every role (orchestr_ai's adsorption-ladder sites on free Se):
            Delta E_bind = E(k) + E(CdCl2) - E(k+1); reaches facets that carry no CdCl2 to remove (e.g. {111}A);
  strip     dots flagged `strip`: greedy desorption to k = 0, each step `candidates` units sampled across roles
            (one per role in turn), relaxed, the lowest polished and kept.
State is saved after every step (state.json); a rerun continues where it stopped (relaxations come from the cache).

    python -m qd_workflows.shape_ladder input  <library dir> ... --ids set.json --strip set.json --out input.json
    python -m qd_workflows.shape_ladder run    <input.json> --out DIR --model M [--head omat_pbe] [--sites 3] [--candidates 4]
    python -m qd_workflows.shape_ladder prepare <input.json> --out DIR --cluster nhr --code SRC [--partition P] [--pieces]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from . import sxl_engine

BOND = 3.05        # A: Cd-Se and Cd-Cl bonds
LINK = 3.3         # A: atoms closer than this belong to one piece (connectivity check)
TOL = 1.0          # A: a unit's Cd this close to a fitted facet plane lies on that facet (covers the Cd just below
                   # Se-terminated {111}, 0.9 A under the outer Se layer)
MIN_AREA = 15.0    # A^2: a facet of the fitted polyhedron smaller than about one surface site (16.7 A^2 on {111})
                   # does not count as exposed
CORE = ("Cd", "Se")
CATION, LIGAND = "Cd", "Cl"


# ---------------- geometry and roles (numpy only) ----------------
def fitted_planes(symbols, positions, planes):
    """[(n, d_fit, label)]: each facet plane moved to the outermost core atom along its normal. Positions must be
    centred on the core (the polyhedron is built around the origin)."""
    s, p = np.asarray(symbols), np.asarray(positions, float)
    core = p[np.isin(s, CORE)]
    return [(np.asarray(n, float), float((core @ np.asarray(n, float)).max()), lab) for n, lab in planes]


def exposed(fplanes):
    """Indices of the fitted planes that keep a facet of non-zero area on their polyhedron."""
    from .shapes import polyhedron
    g = polyhedron([(n, d, str(i)) for i, (n, d, _lab) in enumerate(fplanes)])
    return sorted(int(k) for k, a in g["areas_A2"].items() if a >= MIN_AREA)


def centred(symbols, positions):
    s, p = np.asarray(symbols), np.asarray(positions, float)
    return p - p[np.isin(s, CORE)].mean(0)


def neighbours(positions, cutoff):
    from scipy.spatial import cKDTree
    t = cKDTree(np.asarray(positions, float))
    nb = [set() for _ in range(len(positions))]
    for i, j in t.query_pairs(cutoff, output_type="ndarray"):
        nb[i].add(int(j))
        nb[j].add(int(i))
    return nb


def connected(nb, removed):
    alive = [i for i in range(len(nb)) if i not in removed]
    if not alive:
        return True
    seen, todo = {alive[0]}, [alive[0]]
    while todo:
        for j in nb[todo.pop()]:
            if j not in removed and j not in seen:
                seen.add(j)
                todo.append(j)
    return len(seen) == len(alive)


def units(symbols, positions, planes):
    """Removable CdCl2 units with their site role: [{"atoms": [cd, cl, cl], "role": str}]."""
    s, p = np.asarray(symbols), centred(symbols, positions)
    cl_idx = np.flatnonzero(s == LIGAND)
    if len(cl_idx) < 2:
        return []
    nb = neighbours(p, BOND)
    link = neighbours(p, LINK)
    fp = fitted_planes(s, p, planes)
    exp = exposed(fp)
    out = []
    for cd in np.flatnonzero(s == CATION):
        own = sorted((j for j in nb[cd] if s[j] == LIGAND), key=lambda j: (sum(s[x] == CATION for x in nb[j]), np.linalg.norm(p[j] - p[cd])))
        if not own:
            continue
        pick = own[:2]
        if len(pick) < 2:
            rest = [j for j in cl_idx if j not in pick]
            pick.append(min(rest, key=lambda j: np.linalg.norm(p[j] - p[cd])))
        removed = {int(cd), *map(int, pick)}
        if any(not ({x for x in nb[j] if s[x] == CATION} - removed) for j in cl_idx if j not in removed):
            continue                                         # a Cl would be left without a Cd
        if not connected(link, removed):
            continue
        out.append({"atoms": [int(cd), *map(int, pick)], "role": role_of(p[cd], fp, exp)})
    return out


def role_of(point, fp, exp):
    """Site role of a Cd at `point` (core-centred frame) from the exposed fitted planes it lies within TOL of."""
    touch = [i for i in exp if fp[i][1] - np.asarray(point) @ fp[i][0] < TOL]
    return ("subsurface" if not touch else f"facet:{fp[touch[0]][2]}" if len(touch) == 1
            else f"edge:{'|'.join(sorted(fp[i][2] for i in touch))}" if len(touch) == 2 else "vertex")


def add_sites(symbols, positions, planes, bulk_bond):
    """Z-type CdCl2 additions on under-coordinated Se (orchestr_ai's adsorption-ladder site enumeration), each with
    the site role of the new Cd: [{"add": (symbols, positions), "role": str}]. Needs orchestr_ai.qd (cluster)."""
    from orchestr_ai.qd.steps import detachment as dt
    s, p = list(symbols), centred(symbols, positions)
    state = dt.State(s, p, list(range(len(s))), 0.0, [])
    charges, native, units_ = {"Cd": 2, "Se": -2, "Cl": -1}, ["Cd", "Se"], {"M": CATION, "X": LIGAND, "q": 2}
    d_mx = dt._surface_mx_bond(s, p, charges, native, bulk_bond, units_, 2.45)
    fp = fitted_planes(s, p, planes)
    exp = exposed(fp)
    return [{"add": site["add"], "role": role_of(site["center"], fp, exp)}
            for site in dt.enumerate_sites(state, charges, native, bulk_bond, units_, d_mx)], p


def remove(symbols, positions, atoms):
    keep = [i for i in range(len(symbols)) if i not in set(atoms)]
    return [symbols[i] for i in keep], np.asarray(positions, float)[keep]


def site_key(point, positions, n=24):
    """Environment fingerprint of a site: its sorted distances to the n nearest atoms (0.05 A). Symmetry-equivalent
    sites of a dot share it, so sampling them again would only repeat the same energy."""
    d = np.sort(np.linalg.norm(np.asarray(positions, float) - np.asarray(point, float), axis=1))[:n]
    return tuple(np.round(d / 0.05).astype(int))


def distinct(items, key):
    seen, out = set(), []
    for it in items:
        k = key(it)
        if k not in seen:
            seen.add(k)
            out.append(it)
    return out


def sample(us, per_role, rng, positions):
    """Up to per_role symmetry-distinct units of every role, spread over the dot (farthest-point on the Cd
    positions)."""
    p0 = np.asarray(positions, float)
    us = distinct(us, lambda u: (u["role"], site_key(p0[u["atoms"][0]], p0)))
    by = defaultdict(list)
    for u in us:
        by[u["role"]].append(u)
    out = []
    p = np.asarray(positions, float)
    for role in sorted(by):
        group = by[role]
        first = group[int(rng.integers(len(group)))]
        chosen = [first]
        while len(chosen) < min(per_role, len(group)):
            d = [min(np.linalg.norm(p[u["atoms"][0]] - p[c["atoms"][0]]) for c in chosen) for u in group]
            chosen.append(group[int(np.argmax(d))])
        out += chosen
    return out


def round_robin(us, n, rng):
    """n units sampled across roles: one per role in turn, random within a role."""
    by = defaultdict(list)
    for u in us:
        by[u["role"]].append(u)
    for g in by.values():
        rng.shuffle(g)
    out, roles = [], sorted(by)
    while len(out) < n and any(by.values()):
        for r in roles:
            if by[r] and len(out) < n:
                out.append(by[r].pop())
    return out


# ---------------- input ----------------
def planes_for(cif: str, recipe: dict):
    """(outward normal, family label) of every facet of the builder recipe (needs QD_builder)."""
    from .shapes import construction_planes
    return [(n.round(6).tolist(), lab) for n, _d, lab in construction_planes(cif, recipe)]


def make_input(libs: list[str], ids: list[str], strip: list[str], out: str) -> dict:
    """Dots of the shape libraries with their facet normals, for `run` (the cluster needs no builder)."""
    dots, want = [], set(ids) | set(strip)
    for lib in libs:
        lib = Path(lib)
        summ = json.loads((lib / "summary.json").read_text())
        cfg = __import__("yaml").safe_load(Path(summ["config"]).read_text())
        cif = str((Path(summ["config"]).parent / cfg["cif"]).resolve())
        for s in summ["structures"]:
            if s["id"] not in want:
                continue
            L = (lib / "structures" / s["id"] / "start.xyz").read_text().splitlines()
            n = int(L[0])
            recipe = {**cfg["preset"], "facets": s["facets"], "size_unit_cells": [s["size"] * a for a in s["aspect"]],
                      "construction_origin": ({"center_on_species": s["centre"]} if s["centre"] in CORE else
                                              {"fractional_shift": [float(x) for x in s["centre"].split("(")[1].rstrip(")").split()]})}
            dots.append({"id": s["id"], "library": lib.name, "m": s["m"], "k": s["k"], "n_atoms": n,
                         "symbols": [l.split()[0] for l in L[2:2 + n]],
                         "positions": [[float(x) for x in l.split()[1:4]] for l in L[2:2 + n]],
                         "planes": planes_for(cif, recipe), "strip": s["id"] in set(strip)})
    missing = want - {d["id"] for d in dots}
    if missing:
        raise SystemExit(f"not found in the libraries: {sorted(missing)[:5]}...")
    Path(out).write_text(json.dumps({"dots": dots}))
    print(f"{len(dots)} dots ({sum(d['strip'] for d in dots)} to strip) -> {out}")
    return {"dots": len(dots)}


# ---------------- run (MACE, on the cluster) ----------------
class Ladder:
    def __init__(self, inp, out, settings, sites=3, candidates=4, seed=0, log=print, bulk_bond=2.686):
        from orchestr_ai.qd.steps import detachment as dt
        self.dt, self.s, self.log = dt, settings, log
        self.bulk_bond = bulk_bond          # A, bulk Cd-Se bond (zinc-blende, a = 6.203 A: a sqrt(3) / 4)
        self.dots = {d["id"]: d for d in json.loads(Path(inp).read_text())["dots"]}
        self.out = Path(out)
        self.out.mkdir(parents=True, exist_ok=True)
        self.sites, self.candidates = sites, candidates
        self.rng = np.random.default_rng(seed)
        sxl_engine.configure()
        self.cache = dt.RelaxCache(SimpleNamespace(props=self.out, settings=settings, record_dir=self.out))
        self.E_MX = self._mx_energy()
        sp = self.out / "state.json"
        self.state = json.loads(sp.read_text()) if sp.is_file() else {"relaxed": {}, "binding": {}, "strip": {}}

    def _mx_energy(self):
        """MACE energy of the CdCl2 monomer, relaxed here (same settings as the dots)."""
        mx = self.cache.relax_many([(["Cd", "Cl", "Cl"], np.array([[0, 0, 0], [2.3, 0, 0], [-2.3, 0, 0]], float))],
                                   self.s, polish=True)[0]
        return float(mx[1])

    def save(self):
        tmp = self.out / "state.tmp"
        tmp.write_text(json.dumps(self.state))
        os.replace(tmp, self.out / "state.json")
        self.cache._save()

    def relax(self, items, polish=True):
        """items: [(symbols, positions)] -> [(positions, energy)], one batch, polished."""
        if not items:
            return []
        res = self.cache.relax_many(items, self.s)
        if polish:
            res = self.cache.relax_many([(r[0].get_chemical_symbols(), r[0].get_positions()) for r in res], self.s, polish=True)
        return [(r[0].get_positions(), float(r[1])) for r in res]

    def run(self):
        t0 = time.time()
        todo = [i for i in self.dots if i not in self.state["relaxed"]]
        for i, (p, e) in zip(todo, self.relax([(self.dots[i]["symbols"], self.dots[i]["positions"]) for i in todo])):
            self.state["relaxed"][i] = {"positions": p.round(6).tolist(), "E": e}
        self.save()
        self.log(f"[ladder] relaxed {len(todo)} dots ({time.time() - t0:.0f} s); E(CdCl2) = {self.E_MX:.4f} eV")
        # binding by site class at full coverage
        jobs = []
        for i, d in self.dots.items():
            if i in self.state["binding"]:
                continue
            r = self.state["relaxed"][i]
            us = sample(units(d["symbols"], r["positions"], d["planes"]), self.sites, self.rng, r["positions"])
            for u in us:
                jobs.append((i, u, *remove(d["symbols"], r["positions"], u["atoms"])))
        for (i, u, _s, _p), (p, e) in zip(jobs, self.relax([(s, p) for _i, _u, s, p in jobs])):
            E0 = self.state["relaxed"][i]["E"]
            self.state["binding"].setdefault(i, []).append({"role": u["role"], "atoms": u["atoms"], "E": e,
                                                            "dE_bind": e + self.E_MX - E0})
        self.save()
        self.log(f"[ladder] binding: {len(jobs)} removals ({time.time() - t0:.0f} s)")
        # adsorption by site class: Z-type CdCl2 added on free Se sites (reaches facets that carry no CdCl2)
        jobs = []
        for i, d in self.dots.items():
            if i in self.state.setdefault("adsorb", {}):
                continue
            r = self.state["relaxed"][i]
            sites, pc = add_sites(d["symbols"], r["positions"], d["planes"], self.bulk_bond)
            sites = distinct(sites, lambda st: (st["role"], site_key(np.asarray(st["add"][1], float)[0], pc)))
            by = defaultdict(list)
            for st in sites:
                by[st["role"]].append(st)
            for role in sorted(by):
                for st in [by[role][j] for j in self.rng.permutation(len(by[role]))[:self.sites]]:
                    syms, pts = st["add"]
                    jobs.append((i, role, list(d["symbols"]) + list(syms), np.vstack([pc, np.asarray(pts, float)])))
        for (i, role, _s, _p), (p, e) in zip(jobs, self.relax([(s, p) for _i, _r, s, p in jobs])):
            E0 = self.state["relaxed"][i]["E"]
            self.state["adsorb"].setdefault(i, []).append({"role": role, "E": e, "dE_bind": E0 + self.E_MX - e})
        for i in self.dots:
            self.state["adsorb"].setdefault(i, [])
        self.save()
        self.log(f"[ladder] adsorption: {len(jobs)} additions ({time.time() - t0:.0f} s)")
        # full stripping, all strip dots in lockstep
        while True:
            jobs = []
            for i, d in self.dots.items():
                if not d["strip"]:
                    continue
                st = self.state["strip"].setdefault(i, {"steps": [], "symbols": d["symbols"],
                                                        "positions": self.state["relaxed"][i]["positions"],
                                                        "E": self.state["relaxed"][i]["E"], "done": False})
                if st["done"]:
                    continue
                us = units(st["symbols"], st["positions"], d["planes"])
                if not us:
                    st["done"] = True
                    st["stopped"] = "no removable unit" if st["symbols"].count(LIGAND) else "bare"
                    continue
                p0 = centred(st["symbols"], st["positions"])
                us = distinct(us, lambda u: (u["role"], site_key(p0[u["atoms"][0]], p0)))
                for u in round_robin(us, self.candidates, self.rng):
                    jobs.append((i, u, *remove(st["symbols"], st["positions"], u["atoms"])))
            if not jobs:
                break
            res = self.relax([(s, p) for _i, _u, s, p in jobs], polish=False)
            best = {}
            for (i, u, s, _p), (p, e) in zip(jobs, res):
                if i not in best or e < best[i][3]:
                    best[i] = (u, s, p, e)
            pol = self.relax([(b[1], b[2]) for b in best.values()])
            for (i, (u, s, _p, _e)), (p, e) in zip(best.items(), pol):
                st = self.state["strip"][i]
                cands = [{"role": uu["role"], "E": ee} for (ii, uu, _s, _pp), (_x, ee) in zip(jobs, res) if ii == i]
                st["steps"].append({"k": s.count(LIGAND) // 2, "role": u["role"], "E": e, "dE": e + self.E_MX - st["E"],
                                    "candidates": cands})
                st["symbols"], st["positions"], st["E"] = s, p.round(6).tolist(), e
                if not s.count(LIGAND):
                    st["done"], st["stopped"] = True, "bare"
            self.save()
            left = sum(1 for i, d in self.dots.items() if d["strip"] and not self.state["strip"][i]["done"])
            self.log(f"[ladder] strip step: {len(jobs)} candidates, {left} dots still stripping ({time.time() - t0:.0f} s)")
        self.write()
        self.log("[ladder] done")

    def write(self):
        rows = []
        for i, d in self.dots.items():
            by, ad = defaultdict(list), defaultdict(list)
            for x in self.state["binding"].get(i, []):
                by[x["role"]].append(x["dE_bind"])
            for x in self.state.get("adsorb", {}).get(i, []):
                ad[x["role"]].append(x["dE_bind"])
            st = self.state["strip"].get(i)
            rows.append({"id": i, "library": d["library"], "m": d["m"], "k": d["k"], "n_atoms": d["n_atoms"],
                         "E": self.state["relaxed"][i]["E"],
                         "binding": {r: {"n": len(v), "mean": float(np.mean(v)), "min": float(min(v)), "max": float(max(v))} for r, v in by.items()},
                         "adsorption": {r: {"n": len(v), "mean": float(np.mean(v)), "min": float(min(v)), "max": float(max(v))} for r, v in ad.items()},
                         "strip": None if not st else {"steps": len(st["steps"]), "stopped": st.get("stopped"), "E_end": st["E"],
                                                       "dE": [x["dE"] for x in st["steps"]], "roles": [x["role"] for x in st["steps"]]}})
            if st and st["done"]:
                s, p = st["symbols"], st["positions"]
                (self.out / f"{i}_end.xyz").write_text(f"{len(s)}\n{i} end of stripping E={st['E']:.6f}\n" +
                                                       "\n".join(f"{a} {x:.5f} {y:.5f} {z:.5f}" for a, (x, y, z) in zip(s, p)) + "\n")
        (self.out / "ladder.json").write_text(json.dumps({"E_CdCl2": self.E_MX, "dots": rows}, indent=1))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python -m qd_workflows.shape_ladder")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("input", help="dots of shape libraries with their facet normals")
    p.add_argument("libs", nargs="+")
    p.add_argument("--ids", required=True, help="JSON file with a list of ids, or with keys 'ids' / 'ids_tier1', 'ids_tier2'")
    p.add_argument("--strip", default="", help="comma-separated ids to strip to the bare core (or a JSON list file)")
    p.add_argument("--out", required=True)
    p = sub.add_parser("run", help="relax, binding by site class, stripping (MACE)")
    p.add_argument("input")
    p.add_argument("--out", required=True)
    p.add_argument("--model", default=os.environ.get("QDPROPS_MACE_MODEL"))
    p.add_argument("--head", default="omat_pbe")
    p.add_argument("--sites", type=int, default=3)
    p.add_argument("--candidates", type=int, default=4)
    p = sub.add_parser("prepare", help="GPU job script (self-restarting pieces)")
    p.add_argument("input")
    p.add_argument("--out", required=True)
    p.add_argument("--cluster", default="default")
    p.add_argument("--code", required=True)
    p.add_argument("--partition")
    p.add_argument("--time", default="01:00:00")
    p.add_argument("--pieces", action="store_true")
    p.add_argument("--extra", default="")
    a = ap.parse_args(argv)
    if a.cmd == "input":
        raw = json.loads(Path(a.ids).read_text())
        ids = raw if isinstance(raw, list) else raw.get("ids", []) + raw.get("ids_tier1", []) + raw.get("ids_tier2", [])
        if a.strip.endswith(".json"):
            rs = json.loads(Path(a.strip).read_text())
            strip = rs["strip"] if isinstance(rs, dict) else rs          # a set file from `shapes select`, or a list
        else:
            strip = [x for x in a.strip.split(",") if x]
        make_input(a.libs, ids, strip, a.out)
    elif a.cmd == "prepare":
        from .cluster import load_cluster
        from .grid import prepare
        print(json.dumps(prepare(a.input, a.out, load_cluster(a.cluster), a.code, time_=a.time, extra=a.extra,
                                 partition=a.partition, pieces=a.pieces, module="qd_workflows.shape_ladder"), indent=1))
    else:
        from orchestr_ai.qd.run import Settings
        s = Settings(model=a.model, head=a.head, detach_quick=False, detach_screen=False, detach_rigid_refresh=False,
                     detach_frames=False, hessian="analytic")
        lad = Ladder(a.input, a.out, s, a.sites, a.candidates, log=lambda x: print(x, flush=True))

        def stop(signum, _frame):
            lad.save()
            lad.log(f"[ladder] signal {signum}: state saved, stopping")
            sys.exit(128 + signum)

        import signal
        signal.signal(signal.SIGTERM, stop)
        lad.run()
    return 0


if __name__ == "__main__":
    sys.exit(main())
