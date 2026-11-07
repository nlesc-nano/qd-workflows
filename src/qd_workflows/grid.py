"""
MA detachment grid: [MA]_m(MX_q)_k for every m below the dot's and every k that fits (MACE-MH-1).

The dot's MX_q ladder (props/relaxed.xyz, detach_<j>.xyz, adsorb_<j>.xyz from the detachment step)
is the top row, m = n_MA. The grid is built one row at a time, m -> m - 1:

    descend   every kept structure of (m, k) loses one bonded M-A pair (every pair tried, grouped
              by product fingerprint, relaxed), giving (m - 1, k);
    drop      products whose bond graph splits (a unit left the dot);
    sweep     every node of the row adds one MX_q (enumerate_sites) and loses one
              (enumerate_units); a lower product replaces its neighbour, until no node changes.
              The top of the row grows until no site is left: k_max(m - 1).

Every step is Orchestr.AI's detachment search (expand) with all candidates relaxed (no screening,
no local updates: in a small dot a removal relaxes the whole cluster), candidates within
POLISH_WINDOW of the best polished. A node keeps its structures within `delta` of its best, at
most `beam`. The candidates of a whole row are relaxed together in one batch before the
per-node searches run (which then find them in the cache), so the GPU sees large batches.

Energies are raw MACE totals; references (CdSe(CdCl2)_n complexes, chemical potentials) are left
to the analysis. Outputs in --out: grid.json (nodes, ensembles, timings), nodes/m<m>_k<k>.xyz,
detach_cache.json (relaxations, checkpointed: a rerun continues where the last one stopped).

    python -m qd_workflows.grid run <record dir> --out <dir> --model <mace model> [--head omat_pbe]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from . import sxl_engine

INF = 1e9
BEAM = 3              # structures kept per node
DELTA = 0.05          # eV above the node's best (~2 kT at 300 K)
EXTRA = 4             # children asked from each search beyond `beam`, so fragmented ones can be dropped
MAX_PASSES = 60       # row sweeps (the top of a row grows by one per pass)
SAVE_EVERY = 60.0     # s between checkpoints of the relaxation cache
CHANGE_TOL = 0.003    # eV: a node changes only if its best drops by more (relaxation noise of one minimum
                      # reached from two parents is 0.1-1 meV after polishing; below this the incumbent stays)


def _dt():
    from orchestr_ai.qd.steps import detachment
    return detachment


# ---------------- candidates ----------------
def enumerate_pairs(symbols, pts, charges, native, bulk_bond):
    """Bonded M-A pairs: every anion A with each cation M it is bonded to, in the unit dicts of
    detachment.enumerate_units (cation = M, ligands = [A]) so the search uses them unchanged."""
    from orchestr_ai.qd.steps.structure import _bonds, _cutoffs
    p = np.asarray(pts, float)
    nb = _bonds(list(symbols), p, _cutoffs(symbols, charges, native, bulk_bond))
    cat = next(e for e in native if charges.get(e, 0) > 0)
    an = next(e for e in native if charges.get(e, 0) < 0)
    lig = {i for i, s in enumerate(symbols) if s not in set(native)}
    units = []
    for a in (i for i, s in enumerate(symbols) if s == an):
        for c in sorted(j for j in nb[a] if symbols[j] == cat):
            units.append({"cation": c, "ligands": [a], "atoms": frozenset((c, a)), "molecular": True,
                          "cation_native_cn": sum(1 for j in nb[c] if symbols[j] == an),
                          "cation_ligands": sum(1 for j in nb[c] if j in lig),
                          "ligand_mu": [sum(1 for j in nb[a] if symbols[j] == cat)],
                          "M_X_A": [round(float(np.linalg.norm(p[a] - p[c])), 3)]})
    return units


def fragments(symbols, pts, charges, native, bulk_bond) -> int:
    """Connected pieces of the cation-anion bond graph (1: one dot)."""
    from orchestr_ai.qd.steps.structure import _bonds, _cutoffs
    nb = _bonds(list(symbols), np.asarray(pts, float), _cutoffs(symbols, charges, native, bulk_bond))
    seen, n = set(), 0
    for i in range(len(symbols)):
        if i in seen:
            continue
        n += 1
        todo = [i]
        while todo:
            j = todo.pop()
            if j not in seen:
                seen.add(j)
                todo += list(nb[j] - seen)
    return n


def enumerate_sites(state, charges, native, bulk_bond, units, d_mx, cn_bulk):
    """detachment.enumerate_sites with the anion's bulk coordination given (`cn_bulk`) instead of taken
    from the largest coordination in the structure: a small dot loses its only fully coordinated anion
    along the grid (Cd16Se13Cl6 has one Se of CN 4), after which every CN-3 Se would count as saturated
    and the missing-bond directions would be those of the wrong coordination (12 sites -> 3).
    Otherwise the same rules (to be moved into Orchestr.AI as an optional argument)."""
    from scipy.spatial import cKDTree
    from orchestr_ai.qd.steps.structure import _bonds, _cutoffs
    dt = _dt()
    sym, p = list(state.sym), state.pts
    nb = _bonds(sym, p, _cutoffs(sym, charges, native, bulk_bond))
    an_n = [e for e in native if charges.get(e, 0) < 0][0]
    m, x, q = units["M"], units["X"], units["q"]
    anions = [i for i, s in enumerate(sym) if s == an_n]
    cn_of = {i: [j for j in nb[i] if charges.get(sym[j], 0) > 0] for i in anions}
    centroid = p.mean(axis=0)
    tree = cKDTree(p)
    d_ma = float(bulk_bond)
    out = []
    for a in anions:
        if len(cn_of[a]) >= cn_bulk:
            continue
        outward = dt._unit_vec(p[a] - centroid)
        for d in dt._missing_directions([p[j] - p[a] for j in cn_of[a]], cn_bulk):
            if float(np.dot(d, outward)) < 0.1:
                continue
            mpos = p[a] + d_ma * d
            if any(j != a for j in tree.query_ball_point(mpos, dt.ADD_CLEAR * d_ma)):
                continue
            best = None
            for phi in np.deg2rad(np.arange(0.0, 120.0, 10.0)):
                sites = []
                for e in dt._ligand_directions(d, cn_bulk, phi):
                    xpos = mpos + d_mx * e
                    if tree.query_ball_point(xpos, dt.ADD_CLEAR * d_mx):
                        continue
                    partners = [j for j in tree.query_ball_point(xpos, dt.ADD_PARTNER * d_mx)
                                if charges.get(sym[j], 0) > 0]
                    sites.append((len(partners), float(np.dot(e, outward)), xpos, partners))
                sites.sort(key=lambda s: (-s[0], -s[1]))
                if len(sites) < q:
                    continue
                score = (sum(s[0] for s in sites[:q]), sum(s[1] for s in sites[:q]))
                if best is None or score > best[0]:
                    best = (score, sites[:q])
            if best is None:
                continue
            xs = best[1]
            partners = sorted({state.ids[j] for s in xs for j in s[3]})
            out.append({"add": ([m] + [x] * q, np.vstack([mpos] + [s[2] for s in xs])), "center": mpos,
                        "host": a, "key": ("+", state.ids[a], tuple(partners)),
                        "bridges": [s[0] for s in xs], "molecular": True})
    return out


def cations_without_anion(symbols, pts, charges, native, bulk_bond) -> int:
    """Cations bonded to no anion of the core (held by ligands only)."""
    from orchestr_ai.qd.steps.structure import _bonds, _cutoffs
    nb = _bonds(list(symbols), np.asarray(pts, float), _cutoffs(symbols, charges, native, bulk_bond))
    cat = next(e for e in native if charges.get(e, 0) > 0)
    an = next(e for e in native if charges.get(e, 0) < 0)
    return sum(1 for i, s in enumerate(symbols) if s == cat and not any(symbols[j] == an for j in nb[i]))


# ---------------- the grid ----------------
class Grid:
    def __init__(self, record_dir, out_dir, settings, beam=BEAM, delta=DELTA, m_min=1, log=print):
        from orchestr_ai.qd.references import binary_units
        from orchestr_ai.qd.run import load_context
        dt = _dt()
        self.dt, self.s, self.beam, self.delta, self.m_min, self.log = dt, settings, beam, delta, m_min, log
        self.record_dir, self.out = Path(record_dir).resolve(), Path(out_dir).resolve()
        (self.out / "nodes").mkdir(parents=True, exist_ok=True)
        ctx = load_context(self.record_dir, settings)
        self.charges, self.native = ctx.charges, ctx.native
        props = self.record_dir / "props"
        self.bulk_bond = json.loads((props / "structure.json").read_text())["summary"]["bulk_bond_A"]
        self.e_ref = json.loads((props / "relax.json").read_text())["summary"]["energy_eV"]
        sym0, pts0 = ctx.relaxed()
        self.u = binary_units(sym0, self.charges, self.native)
        if self.u is None or not self.u["X"]:
            raise ValueError("not an [MA]_n(MX_q)_m dot with Z-type ligands")
        mono = 1.08 * 2.4                                   # fallback M-X bond (no ligand on the dot)
        self.d_mx = dt._surface_mx_bond(sym0, np.asarray(pts0, float), self.charges, self.native,
                                        self.bulk_bond, self.u, mono)
        # bulk coordination of the anion: the largest in the intact dot (its core), kept for every node
        from orchestr_ai.qd.steps.structure import _bonds, _cutoffs
        nb0 = _bonds(list(sym0), np.asarray(pts0, float), _cutoffs(sym0, self.charges, self.native, self.bulk_bond))
        self.cn_bulk = max(sum(1 for j in nb0[i] if self.charges.get(sym0[j], 0) > 0)
                           for i, s in enumerate(sym0) if s == self.u["A"])
        sxl_engine.configure()
        self.cache = dt.RelaxCache(SimpleNamespace(props=self.out, settings=settings, record_dir=self.record_dir))
        # Symmetry numbers are not tracked during the search (each node's structures get theirs at the
        # end): a one-element rotation list makes every configuration_sigma 1, without point-group work.
        self.cache.perms = [np.arange(100000)]
        # The cache rewrites its whole JSON after every batch; the grid makes thousands of small batches
        # (tens of MB by the middle rows), so it is written at most every SAVE_EVERY s and after each phase.
        self._save_cache, self._saved, self._saved_n = self.cache._save, 0.0, self.cache.n_new
        self.cache._save = self._save_cache_throttled
        self.stats = {"relaxations": 0, "lazy": 0}
        self.nodes = {}                                      # (m, k) -> node dict
        self.timings, self.warnings = [], []
        prev = self.out / "grid.json"                        # earlier job pieces (each replays from the cache)
        self.pieces = json.loads(prev.read_text()).get("pieces", []) if prev.is_file() else []
        self.pieces.append({"job": os.environ.get("SLURM_JOB_ID"), "start": time.strftime("%Y-%m-%d %H:%M:%S"),
                            "timings": self.timings})
        self.t0 = time.time()

    # -- helpers --
    def _save_cache_throttled(self, force=False):
        """Checkpoint the relaxation cache: at most every SAVE_EVERY s (force: now), and only if it grew."""
        if self.cache.n_new == self._saved_n:
            return
        if force or time.time() - self._saved > SAVE_EVERY:
            self._save_cache()
            self._saved, self._saved_n = time.time(), self.cache.n_new

    def mk(self, symbols):
        a = sum(1 for s in symbols if s == self.u["A"])
        x = sum(1 for s in symbols if s == self.u["X"])
        return a, x // self.u["q"]

    def state(self, symbols, pts, energy_abs, src):
        st = self.dt.State(symbols, pts, range(len(symbols)), energy_abs - self.e_ref, [])
        st.src = src
        return st

    def enum_pairs(self, st):
        return enumerate_pairs(st.sym, st.pts, self.charges, self.native, self.bulk_bond)

    def enum_add(self, st):
        return enumerate_sites(st, self.charges, self.native, self.bulk_bond, self.u, self.d_mx, self.cn_bulk)

    def enum_remove(self, st):
        return self.dt.enumerate_units(st.sym, st.pts, self.charges, self.native, self.bulk_bond, self.u["q"])

    def phase(self, name, n0, t):
        rec = {"phase": name, "seconds": round(time.time() - t, 1), "relaxations_new": self.cache.n_new - n0,
               "elapsed_s": round(time.time() - self.t0, 1)}
        self._save_cache_throttled(force=True)
        self.timings.append(rec)
        self.log(f"[grid] {name}: {rec['seconds']} s, {rec['relaxations_new']} new relaxations "
                 f"(total {rec['elapsed_s']} s)")
        self.save()

    # -- batching: relax (and polish) every search's candidates at once --
    def prewarm(self, jobs):
        """jobs: [(parents, enum)]. Relaxes the representative of every product group of every job in one
        batch, then polishes, again in one batch, those within POLISH_WINDOW of each job's best (what the
        search polishes). The searches then find them in the cache. Returns each job's number of products
        in the window: the search takes no more children than that (a child beyond the window cannot be
        kept, delta < POLISH_WINDOW, and the search would polish it alone, one structure per GPU call)."""
        dt = self.dt
        reps, owner = [], []
        for ji, (parents, enum) in enumerate(jobs):
            for st in parents:
                for members in dt.plan(st, enum(st), INF).values():
                    reps.append(dt._apply(st, members[0][1]))
                    owner.append(ji)
        widths = [0] * len(jobs)
        if not reps:
            return widths
        res = self.cache.relax_many(reps, self.s)
        by_job = {}
        for i, (ji, (_a, e, _ok)) in enumerate(zip(owner, res)):
            by_job.setdefault(ji, []).append((e, i))
        pick = set()
        for ji, items in by_job.items():
            items.sort()
            window = [i for e, i in items if e <= items[0][0] + dt.POLISH_WINDOW]
            widths[ji] = len(window)
            pick |= set(window)
        if not getattr(self.s, "detach_quick", False):
            self.cache.relax_many([(res[i][0].get_chemical_symbols(), res[i][0].get_positions())
                                   for i in sorted(pick)], self.s, polish=True)
        return widths

    def search(self, parents, enum, move, origin, target, width=None):
        """One step from `parents` (states of one node) to the node `target`: the best distinct products,
        not fragmented; the node's ensemble entries {fingerprint: E} (every relaxed product, its final
        (polished) energy, fragmented ones left out); the fingerprints of the fragmented products."""
        from orchestr_ai.qd.records import fingerprint
        for st in parents:              # a state's candidate table belongs to one search (one kind of move)
            st.table, st.parent_table = {}, {}
        width = self.beam + EXTRA if width is None else max(1, min(width, self.beam + EXTRA))
        kids = self.dt.expand(parents, self.cache, self.s, self.e_ref, {}, self.stats, enum, INF, width, 0.0)
        out = []
        for ch in kids:
            if self.mk(ch.sym) != tuple(target):
                raise RuntimeError(f"{move} from {origin} gave {self.mk(ch.sym)}, not {target}")
            if fragments(ch.sym, ch.pts, self.charges, self.native, self.bulk_bond) > 1:
                continue
            st = self.state(ch.sym, ch.pts, ch.energy + self.e_ref,
                            {"move": move, "from": list(origin), **_step(getattr(ch, "step", {}))})
            out.append(st)
        # The ensemble from the final tables: one entry per relaxed group representative, with its polished
        # energy and geometry when it was polished (expand's own ensemble keeps the unpolished one too).
        ens, frag = {}, set()
        for st in parents:
            for key, entry in st.table.items():
                if not entry.get("exact") or entry.get("rep") != key or "product" not in entry:
                    continue
                a = entry["product"]
                sym, pos = a.get_chemical_symbols(), a.get_positions()
                fp = fingerprint(sym, pos)
                if fragments(sym, pos, self.charges, self.native, self.bulk_bond) > 1:
                    frag.add(fp)
                    continue
                e = st.energy + entry["J"] + self.e_ref
                ens[fp] = min(ens.get(fp, e), e)
        return out, frag, ens

    def merge(self, key, states, ensemble, fragmented=()):
        """Add states to node `key`; keep those within delta of the best, at most beam. A structure met
        again (same fingerprint) replaces the kept one only if lower by more than CHANGE_TOL. True if the
        node is new, its best dropped by more than CHANGE_TOL, or a new structure entered its kept set
        (its neighbours are then searched from it too)."""
        node = self.nodes.get(key)
        new = node is None
        if new:
            node = self.nodes[key] = {"states": [], "ensemble": {}, "fragmented": set()}
        old = list(node["states"])
        pool = {}
        for st in old + list(states):                       # one per structure; the incumbent unless clearly lower
            if st.fp not in pool or st.energy < pool[st.fp].energy - CHANGE_TOL:
                pool[st.fp] = st
        pool = sorted(pool.values(), key=lambda st: st.energy)
        if pool:
            node["states"] = [st for st in pool if st.energy <= pool[0].energy + self.delta][:self.beam]
        for fp, e in ensemble.items():
            node["ensemble"][fp] = min(node["ensemble"].get(fp, e), e)
        node["fragmented"] |= set(fragmented)
        if not node["states"]:
            del self.nodes[key]
            return False
        if new:
            return True
        best_drop = old[0].energy - node["states"][0].energy > CHANGE_TOL
        entered = any(st not in old and all(abs(st.energy - o.energy) > CHANGE_TOL for o in old)
                      for st in node["states"])
        return best_drop or entered

    # -- rows --
    def seed(self):
        """Top row: the dot's MX_q ladder, re-polished on this run's settings."""
        t, n0 = time.time(), self.cache.n_new
        props = self.record_dir / "props"
        from .xyz import read_xyz
        files = [props / "relaxed.xyz"] + sorted(props.glob("detach_*.xyz")) + sorted(props.glob("adsorb_*.xyz"))
        structs = [read_xyz(f) for f in files]
        res = self.cache.relax_many([(list(sy), p) for sy, p in structs], self.s, polish=True)
        for f, (sy, _p), (atoms, e, _ok) in zip(files, structs, res):
            st = self.state(atoms.get_chemical_symbols(), atoms.get_positions(), e, {"move": "seed", "file": f.name})
            self.merge(self.mk(sy), [st], {st.fp: e})
        m0 = self.mk(structs[0][0])[0]
        self.phase(f"seed row m={m0} ({len(files)} structures)", n0, t)
        return m0

    def descend(self, m):
        t, n0 = time.time(), self.cache.n_new
        row = sorted(k for (mm, k) in self.nodes if mm == m)
        widths = self.prewarm([(self.nodes[(m, k)]["states"], self.enum_pairs) for k in row])
        for k, w in zip(row, widths):
            kids, frag, ens = self.search(self.nodes[(m, k)]["states"], self.enum_pairs,
                                          f"-{self.u['M']}{self.u['A']}", (m, k), (m - 1, k), w)
            self.merge((m - 1, k), kids, ens, frag)
        self.phase(f"descend m={m} -> {m - 1} ({len(row)} nodes)", n0, t)

    def sweep(self, m):
        """Neighbours along k at fixed m, until no node changes; the top grows while sites are left."""
        t, n0 = time.time(), self.cache.n_new
        changed = {k for (mm, k) in self.nodes if mm == m}
        passes = 0
        while changed and passes < MAX_PASSES:
            passes += 1
            jobs = []
            for k in sorted(changed):
                if (m, k) not in self.nodes:
                    continue
                st = self.nodes[(m, k)]["states"]
                jobs.append((st, self.enum_add, f"+{self.mx}", (m, k), (m, k + 1)))
                if k > 0:
                    jobs.append((st, self.enum_remove, f"-{self.mx}", (m, k), (m, k - 1)))
            widths = self.prewarm([(p, e) for p, e, *_ in jobs])
            changed = set()
            for (parents, enum, move, origin, target), w in zip(jobs, widths):
                kids, frag, ens = self.search(parents, enum, move, origin, target, w)
                if self.merge(target, kids, ens, frag):
                    changed.add(target[1])
        if changed:
            msg = f"sweep m={m} stopped at MAX_PASSES={MAX_PASSES} with k = {sorted(changed)} still changing"
            self.warnings.append(msg)
            self.log(f"[grid] WARNING {msg}")
        row = sorted(k for (mm, k) in self.nodes if mm == m)
        self.phase(f"sweep m={m} ({passes} passes, k = {row[0]}..{row[-1]})", n0, t)

    @property
    def mx(self):
        return f"{self.u['M']}{self.u['X']}{self.u['q'] if self.u['q'] > 1 else ''}"

    def run(self, hessian=True):
        m = self.seed()
        self.sweep(m)
        while m > self.m_min and any(mm == m for (mm, _k) in self.nodes):
            self.descend(m)
            m -= 1
            self.sweep(m)
        if hessian:
            self.vibrations()
        self.save(final=True)

    def vibrations(self):
        """Hessian of every node's best structure (analytic: the grid's dots are small)."""
        from dataclasses import replace
        from orchestr_ai.qd.references import rotational_symmetry_number
        from orchestr_ai.qd.steps.hessian import thermo
        t, n0 = time.time(), self.cache.n_new
        s_h = replace(self.s, hessian="analytic")
        for key, node in sorted(self.nodes.items()):
            st = node["states"][0]
            fr = np.asarray(self.cache.frequencies(st.sym, st.pts, s_h))
            node["vib"] = {"n_imaginary": int((fr < -10.0).sum()), "lowest_cm1": round(float(fr.min()), 2),
                           "F_vib_eV": np.asarray(thermo(fr, self.s.temperatures)["F_vib_eV"]).round(6).tolist(),
                           "frequencies_cm1": fr.round(3).tolist()}
            for x in node["states"]:
                x.sigma = rotational_symmetry_number(x.sym, x.pts)
        self.phase(f"Hessians ({len(self.nodes)} nodes)", n0, t)

    # -- output --
    def save(self, final=False):
        nodes = []
        for (m, k), node in sorted(self.nodes.items()):
            best = node["states"][0]
            flags = {"cations_without_anion": cations_without_anion(best.sym, best.pts, self.charges, self.native,
                                                                    self.bulk_bond)}
            nodes.append({"m": m, "k": k, "n_atoms": len(best.sym), "formula": _formula(best.sym),
                          "E_eV": round(best.energy + self.e_ref, 6),
                          "states": [{"E_eV": round(x.energy + self.e_ref, 6), "fingerprint": x.fp,
                                      "sigma": getattr(x, "sigma", None), **x.src} for x in node["states"]],
                          "n_configurations": len(node["ensemble"]),
                          "ensemble_E_eV": sorted(round(e, 6) for e in node["ensemble"].values()),
                          "ensemble_fingerprints": [fp for fp, _e in sorted(node["ensemble"].items(), key=lambda x: x[1])],
                          "fragmented_products": len(node["fragmented"]), "flags": flags,
                          **({"vib": node["vib"]} if "vib" in node else {})})
            for i, x in enumerate(node["states"]):
                from .xyz import write_xyz
                write_xyz(self.out / "nodes" / f"m{m:02d}_k{k:02d}{'' if i == 0 else f'_{i}'}.xyz", x.sym, x.pts,
                          f"[{self.u['M']}{self.u['A']}]{m}({self.mx}){k} E={x.energy + self.e_ref:.6f} eV "
                          f"({x.src.get('move')} from {x.src.get('from', x.src.get('file'))})")
        rows = {}
        for n in nodes:
            rows.setdefault(n["m"], []).append(n["k"])
        out = {"record": self.record_dir.name, "unit_removed": f"{self.u['M']}{self.u['A']}", "ligand_unit": self.mx,
               "complete": final, "settings": {"beam": self.beam, "delta_eV": self.delta, "m_min": self.m_min,
                                                "head": self.s.head, "model": Path(self.s.model).name,
                                                "relaxer": self.s.relaxer(), "polish_fmax": self.dt.POLISH_FMAX,
                                                "polish_window_eV": self.dt.POLISH_WINDOW,
                                                "temperatures": self.s.temperatures, "cn_bulk": self.cn_bulk,
                                                "change_tol_eV": CHANGE_TOL},
               "rows": {str(m): {"k_min": min(ks), "k_max": max(ks), "n": len(ks)} for m, ks in sorted(rows.items())},
               "search": {**self.stats, "relaxations_cached_new": self.cache.n_new},
               "warnings": self.warnings, "timings": self.timings, "pieces": self.pieces, "nodes": nodes}
        tmp = self.out / "grid.json.tmp"
        tmp.write_text(json.dumps(out, indent=1))
        os.replace(tmp, self.out / "grid.json")


def _step(step: dict) -> dict:
    keep = ("unit", "J_eV", "cation_native_cn", "cation_ligands", "ligand_mu", "M_X_A", "host", "bridges",
            "n_candidates")
    return {k: (list(v) if isinstance(v, tuple) else v) for k, v in step.items() if k in keep}


def _formula(symbols) -> str:
    return "".join(f"{e}{symbols.count(e)}" for e in sorted(set(symbols), key=list(symbols).index))


def prepare(record: str, out_dir: str, cluster: dict, code: str, *, test: bool = True, time_: str = "01:00:00",
            extra: str = "", partition: str | None = None, pieces: bool = False,
            module: str = "qd_workflows.grid") -> dict:
    """GPU job script (the props job of the cluster profile) that runs the grid of one record; `code` is the
    qd-workflows src directory on the cluster that has this module. pieces: the job queues its successor and
    stops a minute before its limit; the successor continues from the relaxation cache, and the last piece
    (grid finished, or failed) cancels the one it queued."""
    import copy
    from .props import RUN_PIECES, _job, _seconds
    cl = copy.deepcopy(cluster)
    pp = cl["props"]
    if partition:
        pp["partition"] = partition
    elif test:
        pp["partition"] = pp.get("test_partition", pp["partition"])
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    run = (f"env PYTHONPATH={code}:$PYTHONPATH python -u -m {module} run {record} --out {out_dir} "
           f"--model {pp['model']} --head {pp.get('head', 'omat_pbe')} {extra}").rstrip()
    body = RUN_PIECES.format(seconds=max(_seconds(time_) - 60, 60), run=run, job="job.sh", then="") if pieces else run
    (out / "job.sh").write_text(_job(cl, "grid-" + Path(record).name, time_, None, body))
    return {"job": str(out / "job.sh"), "partition": pp["partition"], "pieces": pieces, "run": run}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python -m qd_workflows.grid")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("prepare", help="GPU job script for the grid of one record")
    p.add_argument("record", help="record directory on the cluster")
    p.add_argument("--out", required=True, help="job and output directory on the cluster")
    p.add_argument("--cluster", default="default", help="cluster profile")
    p.add_argument("--code", required=True, help="qd-workflows src directory on the cluster")
    p.add_argument("--regular", action="store_true", help="the regular GPU partition instead of the test one")
    p.add_argument("--time", default="01:00:00")
    p.add_argument("--extra", default="", help="more `run` options, e.g. '--beam 1'")
    p.add_argument("--partition", help="partition(s), e.g. 'gpu-a100:test,gpu-a100:shared' (overrides --regular)")
    p.add_argument("--pieces", action="store_true", help="resubmit itself until the grid is finished")
    p = sub.add_parser("run", help="build the grid of one record (props/ of a finished detachment step)")
    p.add_argument("record", help="record directory (record.json, start.xyz, props/)")
    p.add_argument("--out", required=True, help="output directory (grid.json, nodes/, relaxation cache)")
    p.add_argument("--model", default=os.environ.get("QDPROPS_MACE_MODEL"), help="MACE model")
    p.add_argument("--head", default="omat_pbe")
    p.add_argument("--beam", type=int, default=BEAM, help="structures kept per node")
    p.add_argument("--delta", type=float, default=DELTA, help="eV above a node's best that are kept")
    p.add_argument("--m-min", type=int, default=1, help="lowest number of MA units")
    p.add_argument("--no-hessian", action="store_true", help="skip the Hessians of the nodes")
    a = ap.parse_args(argv)
    if a.cmd == "prepare":
        from .cluster import load_cluster
        print(json.dumps(prepare(a.record, a.out, load_cluster(a.cluster), a.code, test=not a.regular,
                                 time_=a.time, extra=a.extra, partition=a.partition, pieces=a.pieces), indent=1))
        return 0
    from orchestr_ai.qd.run import Settings
    s = Settings(model=a.model, head=a.head, detach_quick=False, detach_screen=False, detach_rigid_refresh=False,
                 detach_frames=False, hessian="analytic")
    g = Grid(a.record, a.out, s, beam=a.beam, delta=a.delta, m_min=a.m_min, log=lambda x: print(x, flush=True))

    def stop(signum, _frame):           # `timeout` at the end of a job piece: keep what was relaxed since the last save
        g._save_cache_throttled(force=True)
        g.log(f"[grid] signal {signum}: relaxation cache saved, stopping")
        sys.exit(128 + signum)

    import signal
    signal.signal(signal.SIGTERM, stop)
    g.run(hessian=not a.no_hessian)
    return 0


if __name__ == "__main__":
    sys.exit(main())
