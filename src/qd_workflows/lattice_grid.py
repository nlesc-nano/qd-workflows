"""
MA detachment grid on the lattice: [MA]_m(MX_q)_k, every structure a configuration of the bulk lattice
(QD_builder rules), relaxed with MACE-MH-1 for its energy.

A node holds lattice configurations: atoms on the sites of the CIF lattice the builder cut the dot from
(M on cation sites, A and X on anion sites, as in the builder's start.xyz). Every move is decided on the
lattice positions, never on a relaxed geometry, so a relaxation cannot change what counts as a site.
The energy of a configuration is the MACE relaxation started from its lattice positions; the same
configuration reached by two paths is relaxed once.

Moves (QD_builder's definitions):
    +MX_q   a Z-type unit on a free virtual cation site above an under-coordinated surface anion
            (builder z_type_ligand_host_sites: anion below its bulk coordination, no ligand atom within
            3.5 A, site clear), its q X on the free anion sites around it (builder _anion_site_candidates,
            bridging and outward sites first);
    -MX_q   a surface M with q of its bonded X;
    -MA     a surface A with a surface M bonded to it.
A removal is allowed only if no X is left without a bonded M and the bond graph stays in one piece
(an isolated MX_q unit, for example, is not). Surface: the builder's shell, within surf_tol of the
outermost plane of each Wulff facet orientation (detect_facets_from_nc on the configuration).

Rows:
    top row     the builder's dot (start.xyz, k0): ladder up (+MX_q until no site) and down (-MX_q to k = 0);
    m -> m - 1  (--descent column, the default) every node (m, k) loses one MA: every allowed pair of each of its
                kept structures, the whole row relaxed in one batch, gives (m - 1, k); one batched neighbour pass
                (+MX_q from (m - 1, k - 1), -MX_q from (m - 1, k + 1)) keeps the lower of each node; then one ladder
                up from the top of row m - 1 fills the sites the removals opened (k_max(m - 1)).
                (--descent top) only the fully passivated node (m, k_max(m)) loses one MA; from there one ladder
                up (k_max(m - 1)) and one ladder down (to k = 0) fill the row.
Each step relaxes all candidates of the step together, polishes those within POLISH_WINDOW of the best
and keeps the configurations within `delta` of the best, at most `beam`, as the parents of the next step.

    python -m qd_workflows.lattice_grid run <record dir> --out <dir> --model <mace model> [--cif CdSe_zb.cif]
"""
from __future__ import annotations

import argparse
import itertools
import json
import os
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from .grid import BEAM, DELTA, SAVE_EVERY, _formula
from .xyz import read_xyz, write_xyz

BOND_TOL = 1.1        # x the bulk nearest-neighbour distance: bonds of a lattice configuration


class Lattice:
    """The builder's view of one record: CIF lattice, facets, surface shell, Z-type sites."""

    def __init__(self, record_dir, cif):
        from pymatgen.core import Structure
        from builder.config import _parse_hkl
        from builder.facets import expand_facets
        from builder.nc_types import Config, Facet, PassivationSpec
        rec = json.loads((Path(record_dir) / "record.json").read_text())
        recipe = rec["origin"]["recipe"]
        self.charges = {k: int(v) for k, v in recipe["charges"].items()}
        self.bulk = Structure.from_file(str(cif))
        self.native = {str(s.specie.symbol) for s in self.bulk.sites}
        self.cation = next(e for e in self.native if self.charges[e] > 0)
        self.anion = next(e for e in self.native if self.charges[e] < 0)
        self.ligand = recipe["passivation"]["ligand"]
        self.q = self.charges[self.cation] // -self.charges[self.ligand]
        self.surf_tol = float(recipe["passivation"].get("surf_tol", 2.0))
        seeds = [Facet(*_parse_hkl(f["hkl"]), float(f.get("gamma", 1.0)), f.get("termination"),
                       f.get("scope", "family")) for f in recipe["facets"]]
        proper = bool((recipe.get("symmetry") or {}).get("proper_rotations_only", True))
        self.facets = expand_facets(self.bulk, seeds, proper_only=proper)
        self.cfg = Config(mode="single", seeds=seeds, aspect=(1.0, 1.0, 1.0), shape_mode="wulff", sphere_planes=0,
                          size_unit_cells=None, proper_only=proper, pair_opposites=True,
                          passivation=PassivationSpec(ligand=self.ligand, surf_tol=self.surf_tol),
                          charges=dict(self.charges), materials=[])
        d = [n.nn_distance for n in self.bulk.get_neighbors(self.bulk.sites[0], r=4.0)]
        self.d_nn = float(min(d))

    # -- geometry of a configuration --
    def planes(self, syms, pts):
        from builder.facets import detect_facets_from_nc
        core = [i for i, s in enumerate(syms) if s in self.native]
        return detect_facets_from_nc([syms[i] for i in core], np.asarray(pts)[core], self.bulk.lattice,
                                     self.charges, self.facets, self.surf_tol)[1]

    def surface(self, syms, pts):
        pts = np.asarray(pts, float)
        mask = np.zeros(len(pts), bool)
        for n, d in self.planes(syms, pts):
            mask |= (float(d) - pts @ np.asarray(n, float)) < self.surf_tol
        return mask

    def bonds(self, syms, pts):
        from scipy.spatial import cKDTree
        nb = [set() for _ in syms]
        for i, j in cKDTree(np.asarray(pts)).query_pairs(BOND_TOL * self.d_nn):
            if self.charges.get(syms[i], 0) * self.charges.get(syms[j], 0) < 0:
                nb[i].add(j)
                nb[j].add(i)
        return nb

    @staticmethod
    def connected(nb, removed):
        alive = [i for i in range(len(nb)) if i not in removed]
        if not alive:
            return False
        seen, todo = set(), [alive[0]]
        while todo:
            a = todo.pop()
            if a in seen:
                continue
            seen.add(a)
            todo += [b for b in nb[a] if b not in removed and b not in seen]
        return len(seen) == len(alive)

    def _allowed(self, syms, nb, removed):
        """No X left without a bonded M, no M left without a bonded A (a Z-type MX_q stays on its anion: an MX_q
        held only through X bridges is not a passivating unit, and the shell would outgrow what the anions of the
        smaller core can host), and the bond graph in one piece."""
        for i in removed:
            for x in nb[i]:
                if x in removed:
                    continue
                if syms[x] == self.ligand and not (nb[x] - removed):
                    return False
                if syms[x] == self.cation and not any(syms[a] == self.anion for a in nb[x] - removed):
                    return False
        return self.connected(nb, removed)

    # -- moves: lists of (description, product symbols, product lattice positions) --
    def z_sites(self, syms, pts):
        """The builder's Z-type host sites of a configuration and, for each, its free anion (X) sites (dicts with
        pos and an id shared by the sites of neighbouring hosts)."""
        from builder.analysis import nc_bulk_frac_shift
        from builder.z_type_ligand_posttreat import (_anion_site_candidates, _free_sites, _site_users,
                                                     z_type_ligand_host_sites)
        pts = np.asarray(pts, float)
        f_shift = nc_bulk_frac_shift(list(syms), pts, self.bulk)
        if f_shift is None:
            raise RuntimeError("a lattice configuration does not map onto the CIF lattice")
        planes = self.planes(syms, pts)
        sites = z_type_ligand_host_sites(list(syms), pts, self.cfg, self.bulk, planes, set(self.native))
        if not sites:
            return [], []
        cands = _anion_site_candidates(sites, list(syms), pts, self.cfg, self.bulk, planes, f_shift)
        users = _site_users(cands)
        return sites, [_free_sites(j, cands, users, set(), set()) for j in range(len(sites))]

    def fillable(self, free):
        """Can every host get q distinct free X sites at once? (bipartite matching, each host taken q times)."""
        slots = [j for j, f in enumerate(free) for _ in range(self.q)]
        match = {}                                         # X-site id -> slot

        def augment(s, seen):
            for c in free[slots[s]]:
                if c["id"] in seen:
                    continue
                seen.add(c["id"])
                if c["id"] not in match or augment(match[c["id"]], seen):
                    match[c["id"]] = s
                    return True
            return False

        return all(augment(s, set()) for s in range(len(slots)))

    def add_mx(self, syms, pts):
        """+MX_q on every builder host site, with every choice of q of its free X sites, keeping only the additions
        after which every remaining host can still be filled (the Se deficits are all passivated at the top of the
        ladder: an X taking a site another unit needs would end the ladder early). When none keeps them fillable
        (the configuration already cannot fill all its hosts), every addition is kept."""
        pts = np.asarray(pts, float)
        sites, free = self.z_sites(syms, pts)
        out, blocked = [], []
        for site, f in zip(sites, free):
            for take in itertools.combinations(f, self.q):
                new = np.vstack([np.asarray(site["pos"], float)] + [np.asarray(c["pos"], float) for c in take])
                prod = (list(syms) + [self.cation] + [self.ligand] * self.q, np.vstack([pts, new]))
                _s2, free2 = self.z_sites(*prod)
                (out if self.fillable(free2) else blocked).append(({"host": int(site["hosts"][0])}, *prod))
        return out or blocked

    def surface_cations(self, syms, pts, nb):
        """Cations of the builder's surface shell that are not fully coordinated by the core anion (in a dot as
        small as Cd16Se13 the shell reaches the M bonded to the central A, which keep their bulk coordination)."""
        surf = self.surface(syms, pts)
        cn_bulk = max(sum(1 for j in nb[i] if syms[j] == self.anion) for i, s in enumerate(syms) if s == self.cation)
        return [i for i, s in enumerate(syms) if s == self.cation and surf[i]
                and sum(1 for j in nb[i] if syms[j] == self.anion) < max(cn_bulk, 4)], surf

    def remove_mx(self, syms, pts):
        """A surface M with q X: its own X, completed by the closest other X when it has fewer than q (on the
        builder's dot every Cl bridges two Cd, so no Cd owns two)."""
        pts = np.asarray(pts, float)
        nb = self.bonds(syms, pts)
        cats, _surf = self.surface_cations(syms, pts, nb)
        xs_all = [i for i, s in enumerate(syms) if s == self.ligand]
        out = []
        for c in cats:
            own = sorted(x for x in nb[c] if syms[x] == self.ligand)
            if not own:
                continue
            if len(own) >= self.q:
                combos = list(itertools.combinations(own, self.q))
            else:
                need = self.q - len(own)
                d = {x: float(np.linalg.norm(pts[x] - pts[c])) for x in xs_all if x not in own}
                if len(d) < need:
                    continue
                others = sorted(d, key=d.get)
                pool = [x for x in others if d[x] <= d[others[need - 1]] + 0.05]
                combos = [tuple(sorted(own + list(e))) for e in itertools.combinations(pool, need)]
            for combo in combos:
                removed = {c, *combo}
                if self._allowed(syms, nb, removed):
                    out.append(({"M": c, "X": list(combo), "own_X": len(own)}, *_without(syms, pts, removed)))
        return out

    def remove_ma(self, syms, pts):
        nb = self.bonds(syms, pts)
        cats, surf = self.surface_cations(syms, pts, nb)
        cats = set(cats)
        out = []
        for a in (i for i, s in enumerate(syms) if s == self.anion and surf[i]):
            for c in sorted(x for x in nb[a] if x in cats):
                removed = {a, c}
                if self._allowed(syms, nb, removed):
                    out.append(({"M": c, "A": a}, *_without(syms, pts, removed)))
        return out


def _without(syms, pts, removed):
    keep = [i for i in range(len(syms)) if i not in removed]
    return [syms[i] for i in keep], np.asarray(pts, float)[keep]


class Config_:
    """A lattice configuration with its relaxed structure and energy."""

    def __init__(self, syms, lat, relaxed, energy, fp, src):
        self.syms, self.lat, self.relaxed, self.energy, self.fp, self.src = syms, lat, relaxed, energy, fp, src


class LatticeGrid:
    def __init__(self, record_dir, out_dir, settings, cif, beam=BEAM, delta=DELTA, m_min=1, log=print, start="lattice"):
        from orchestr_ai.qd.steps import detachment as dt
        self.dt, self.s, self.beam, self.delta, self.m_min, self.log = dt, settings, beam, delta, m_min, log
        if start not in ("lattice", "relaxed"):
            raise ValueError(f"start: lattice or relaxed, not {start!r}")
        self.start = start
        self.record_dir, self.out = Path(record_dir).resolve(), Path(out_dir).resolve()
        (self.out / "nodes").mkdir(parents=True, exist_ok=True)
        self.lat = Lattice(self.record_dir, cif)
        self.cache = dt.RelaxCache(SimpleNamespace(props=self.out, settings=settings, record_dir=self.record_dir))
        self._save_cache, self._saved, self._saved_n = self.cache._save, 0.0, self.cache.n_new
        self.cache._save = self._save_cache_throttled
        self.nodes = {}                          # (m, k) -> {"states": [Config_], "ensemble": {fp: E}, "steps": {...}}
        self.timings = []
        prev = self.out / "grid.json"
        self.pieces = json.loads(prev.read_text()).get("pieces", []) if prev.is_file() else []
        self.pieces.append({"job": os.environ.get("SLURM_JOB_ID"), "start": time.strftime("%Y-%m-%d %H:%M:%S"),
                            "timings": self.timings})
        self.t0 = time.time()

    def _save_cache_throttled(self, force=False):
        if self.cache.n_new == self._saved_n:
            return
        if force or time.time() - self._saved > SAVE_EVERY:
            self._save_cache()
            self._saved, self._saved_n = time.time(), self.cache.n_new

    def mk(self, syms):
        return syms.count(self.lat.anion), syms.count(self.lat.ligand) // self.lat.q

    CLASH_A = 1.4        # Å: a new atom this close to a relaxed atom of the parent falls back to the lattice start

    def _start(self, st, info, syms, pts):
        """(start geometry, label) of a product. start='lattice': its lattice positions (label None). start='relaxed':
        the parent's relaxed coordinates with the removed atoms deleted, or with the new atoms at their lattice sites
        moved by the displacement of their host anion (QD_builder's growth from relaxed coordinates); the lattice
        configuration stays the bookkeeping of the moves."""
        if self.start == "lattice":
            return pts, None
        rel = np.asarray(st.relaxed, float)
        if "host" in info:                                  # +MX_q: appended at the end of the parent's atoms
            n = len(st.syms)
            shift = rel[info["host"]] - np.asarray(st.lat, float)[info["host"]]
            new = np.asarray(pts, float)[n:] + shift
            if np.linalg.norm(new[:, None] - rel[None], axis=-1).min() < self.CLASH_A:
                return pts, "lattice (clash)"
            return np.vstack([rel, new]), "relaxed"
        removed = {info["M"], *([info["A"]] if "A" in info else []), *info.get("X", [])}
        keep = [i for i in range(len(st.syms)) if i not in removed]
        return rel[keep], "relaxed"

    # -- steps: every move of every parent, relaxed together --
    def step(self, parents, moves, name, target):
        return self.step_many([(parents, moves, name, target)])[0]

    def step_many(self, jobs):
        """jobs: [(parents, moves, name, target node)]. The products of every job are relaxed in one batch and
        polished (POLISH_WINDOW of each job's best) in another; each target keeps its configurations within delta
        of its best, at most beam. Returns each target's kept states ([] when the job had no product)."""
        from orchestr_ai.qd.records import fingerprint
        t, n0 = time.time(), self.cache.n_new
        per_job = []
        for parents, moves, name, target in jobs:
            products = {}
            for st in parents:
                for info, syms, pts in moves(st.syms, st.lat):
                    if self.mk(syms) != tuple(target):
                        raise RuntimeError(f"{name} from {self.mk(st.syms)} gave {self.mk(syms)}, not {target}")
                    fp = fingerprint(syms, pts)
                    if fp in products:
                        continue
                    src = {"move": name, "from": [*self.mk(st.syms)], "parent": st.fp, **info}
                    start, how = self._start(st, info, syms, pts)
                    if how is not None:
                        src["start"] = how
                    products[fp] = (syms, pts, src, start)
            per_job.append(products)
        flat = [(ji, f) for ji, products in enumerate(per_job) for f in products]
        if not flat:
            return [[] for _ in jobs]
        res = self.cache.relax_many([(per_job[ji][f][0], per_job[ji][f][3]) for ji, f in flat], self.s)
        e = {jf: r[1] for jf, r in zip(flat, res)}
        relaxed = {jf: r[0] for jf, r in zip(flat, res)}
        best = {}
        for (ji, _f), v in e.items():
            best[ji] = min(best.get(ji, v), v)
        window = [jf for jf in flat if e[jf] <= best[jf[0]] + self.dt.POLISH_WINDOW]
        pol = self.cache.relax_many([(relaxed[jf].get_chemical_symbols(), relaxed[jf].get_positions()) for jf in window],
                                    self.s, polish=True)
        for jf, r in zip(window, pol):
            e[jf], relaxed[jf] = r[1], r[0]
        seconds, n_new = round(time.time() - t, 1), self.cache.n_new - n0
        out = []
        for ji, ((parents, _moves, name, target), products) in enumerate(zip(jobs, per_job)):
            if not products:
                out.append([])
                continue
            order = sorted(products, key=lambda f: e[(ji, f)])
            keep = [f for f in order if e[(ji, f)] <= e[(ji, order[0])] + self.delta][:self.beam]
            node = self.nodes.setdefault(tuple(target), {"states": [], "ensemble": {}, "steps": []})
            for f in products:
                node["ensemble"][f] = min(node["ensemble"].get(f, e[(ji, f)]), e[(ji, f)])
            pool = {st.fp: st for st in node["states"]}
            for f in keep:
                if f not in pool or e[(ji, f)] < pool[f].energy:
                    pool[f] = Config_(products[f][0], products[f][1], relaxed[(ji, f)].get_positions(), e[(ji, f)], f,
                                      products[f][2])
            node["states"] = sorted(pool.values(), key=lambda x: x.energy)
            node["states"] = [x for x in node["states"] if x.energy <= node["states"][0].energy + self.delta][:self.beam]
            rec = {"step": f"{name} -> m={target[0]} k={target[1]}", "parents": len(parents),
                   "candidates": len(products), "polished": sum(1 for j, _f in window if j == ji),
                   "batch": {"jobs": len(jobs), "relaxations_new": n_new, "seconds": seconds},
                   "elapsed_s": round(time.time() - self.t0, 1)}
            node["steps"].append(rec)
            self.timings.append(rec)
            self.log(f"[lattice-grid] {rec['step']}: {len(products)} candidates from {len(parents)} parents, best "
                     f"{e[(ji, order[0])]:.4f} eV")
            out.append(node["states"])
        self.log(f"[lattice-grid] batch of {len(jobs)} step(s): {len(flat)} candidates, {n_new} new relaxations, "
                 f"{seconds} s")
        self._save_cache_throttled(force=True)
        self.save()
        return out

    def ladder_up(self, m, k_top):
        """From the node (m, k_top): +MX_q one at a time until no site is left. Returns k_max(m)."""
        cur, k = self.nodes[(m, k_top)]["states"], k_top
        while True:
            kids = self.step(cur, self.lat.add_mx, f"+{self.mx}", (m, k + 1))
            if not kids:
                return k
            cur, k = kids, k + 1

    def neighbour_pass(self, m, ks):
        """One batched pass along a freshly descended row: every node (m, k) also gets +MX_q from the best structure
        of (m, k - 1) and -MX_q from the best of (m, k + 1), as left by the column descents, in one batch; the lower
        configurations stay. One pass only: the products are not searched further. (From the best structure only:
        with all kept ones the pass had ~3.5x the candidates of the column descents.)"""
        ks = set(ks)
        before = {k: self.nodes[(m, k)]["states"][0].energy for k in ks}
        jobs = []
        for k in range(min(ks), max(ks) + 1):            # a k no column descent reached is filled from its neighbours
            if k - 1 in ks:
                jobs.append((self.nodes[(m, k - 1)]["states"][:1], self.lat.add_mx, f"+{self.mx}", (m, k)))
            if k + 1 in ks:
                jobs.append((self.nodes[(m, k + 1)]["states"][:1], self.lat.remove_mx, f"-{self.mx}", (m, k)))
        if not jobs:
            return
        self.step_many(jobs)
        gained = {k: before[k] - self.nodes[(m, k)]["states"][0].energy for k in sorted(ks)}
        better = {k: round(v, 3) for k, v in gained.items() if v > 1e-6}
        filled = sorted(k for (mm, k) in self.nodes if mm == m and k not in ks and self.nodes[(mm, k)]["states"])
        self.log(f"[lattice-grid] neighbour pass m={m}: lower at k = {sorted(better)} (eV gained: {better})"
                 + (f"; filled k = {filled}" if filled else ""))

    def ladders(self, m, k_top):
        """From the node (m, k_top): one ladder up until no site is left, one ladder down to k = 0."""
        k_max = self.ladder_up(m, k_top)
        cur, k = self.nodes[(m, k_top)]["states"], k_top
        while k > 0:
            kids = self.step(cur, self.lat.remove_mx, f"-{self.mx}", (m, k - 1))
            if not kids:
                break
            cur, k = kids, k - 1
        self.log(f"[lattice-grid] row m={m}: k = {k}..{k_max} (top {k_top})")
        return k_max

    @property
    def mx(self):
        q = self.lat.q
        return f"{self.lat.cation}{self.lat.ligand}{q if q > 1 else ''}"

    def run(self, m_stop=None, descent="column"):
        """descent='column': every node (m, k) of a row loses one MA (all allowed pairs of its kept structures, the
        whole row in one batch), giving (m - 1, k) for every k that has an allowed removal; then one ladder up from
        the top of row m - 1 fills the sites the removals opened (k_max can grow; it shrinks where no removal is
        allowed). descent='top': only the fully passivated node loses one MA, and the row is filled by one ladder up
        and one ladder down from there."""
        from orchestr_ai.qd.records import fingerprint
        syms, pts = read_xyz(self.record_dir / "start.xyz")
        syms = list(syms)
        (atoms, e, _ok), = self.cache.relax_many([(syms, pts)], self.s)
        (atoms, e, _ok), = self.cache.relax_many([(atoms.get_chemical_symbols(), atoms.get_positions())], self.s,
                                                 polish=True)
        m, k = self.mk(syms)
        fp = fingerprint(syms, pts)
        self.nodes[(m, k)] = {"states": [Config_(syms, np.asarray(pts), atoms.get_positions(), e, fp,
                                                 {"move": "builder", "file": "start.xyz"})],
                              "ensemble": {fp: e}, "steps": []}
        self.log(f"[lattice-grid] builder dot m={m} k={k}: E = {e:.4f} eV")
        k_max = self.ladders(m, k)
        m_stop = self.m_min if m_stop is None else m_stop
        ma = f"-{self.lat.cation}{self.lat.anion}"
        while m > m_stop:
            if descent == "top":
                kids = self.step(self.nodes[(m, k_max)]["states"], self.lat.remove_ma, ma, (m - 1, k_max))
                if not kids:
                    self.log(f"[lattice-grid] no allowed {ma[1:]} removal from m={m} k={k_max}")
                    break
                m -= 1
                k_max = self.ladders(m, k_max)
                continue
            row = sorted(k for (mm, k) in self.nodes if mm == m and self.nodes[(mm, k)]["states"])
            kids = self.step_many([(self.nodes[(m, k)]["states"], self.lat.remove_ma, ma, (m - 1, k)) for k in row])
            got = [k for k, st in zip(row, kids) if st]
            if not got:
                self.log(f"[lattice-grid] no allowed {ma[1:]} removal from row m={m}")
                break
            m -= 1
            self.neighbour_pass(m, got)
            k_max = self.ladder_up(m, max(got))
            self.log(f"[lattice-grid] row m={m}: k = {min(got)}..{k_max} (descended: {got})")
        self.save(final=True)

    def save(self, final=False):
        nodes = []
        for (m, k), node in sorted(self.nodes.items()):
            if not node["states"]:
                continue
            best = node["states"][0]
            nodes.append({"m": m, "k": k, "n_atoms": len(best.syms), "formula": _formula(best.syms),
                          "E_eV": round(best.energy, 6),
                          "states": [{"E_eV": round(x.energy, 6), "fingerprint": x.fp, **x.src} for x in node["states"]],
                          "n_configurations": len(node["ensemble"]),
                          "ensemble_E_eV": sorted(round(v, 6) for v in node["ensemble"].values()),
                          "steps": node["steps"]})
            for i, x in enumerate(node["states"]):
                tag = f"m{m:02d}_k{k:02d}{'' if i == 0 else f'_{i}'}"
                write_xyz(self.out / "nodes" / f"{tag}.xyz", x.syms, x.relaxed, f"m={m} k={k} E={x.energy:.6f} eV relaxed")
                write_xyz(self.out / "nodes" / f"{tag}_lattice.xyz", x.syms, x.lat, f"m={m} k={k} lattice configuration")
        rows = {}
        for n in nodes:
            rows.setdefault(n["m"], []).append(n["k"])
        out = {"record": self.record_dir.name, "method": "lattice configurations (QD_builder rules), MACE relaxed",
               "complete": final,
               "settings": {"beam": self.beam, "delta_eV": self.delta, "m_min": self.m_min, "head": self.s.head,
                            "model": Path(self.s.model).name, "relaxer": self.s.relaxer(),
                            "polish_fmax": self.dt.POLISH_FMAX, "polish_window_eV": self.dt.POLISH_WINDOW,
                            "surf_tol_A": self.lat.surf_tol, "bond_tol": BOND_TOL, "d_nn_A": self.lat.d_nn,
                            "start": self.start},
               "rows": {str(m): {"k_min": min(ks), "k_max": max(ks), "n": len(ks)} for m, ks in sorted(rows.items())},
               "relaxations_new": self.cache.n_new, "pieces": self.pieces, "nodes": nodes}
        tmp = self.out / "grid.json.tmp"
        tmp.write_text(json.dumps(out, indent=1))
        os.replace(tmp, self.out / "grid.json")


def _find_cif(record_dir, cif):
    if cif:
        return cif
    name = json.loads((Path(record_dir) / "record.json").read_text())["origin"]["cif"]
    for d in os.environ.get("QDPROPS_CIF_DIRS", "").split(os.pathsep):
        if d and (Path(d) / name).is_file():
            return str(Path(d) / name)
    raise FileNotFoundError(f"{name}: give --cif or set QDPROPS_CIF_DIRS")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python -m qd_workflows.lattice_grid")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("prepare", help="GPU job script for the lattice grid of one record")
    p.add_argument("record")
    p.add_argument("--out", required=True)
    p.add_argument("--cluster", default="default")
    p.add_argument("--code", required=True, help="PYTHONPATH entries on the cluster (qd-workflows src, QD_builder src)")
    p.add_argument("--partition")
    p.add_argument("--time", default="01:00:00")
    p.add_argument("--pieces", action="store_true")
    p.add_argument("--extra", default="")
    p = sub.add_parser("run", help="build the grid of one record")
    p.add_argument("record")
    p.add_argument("--out", required=True)
    p.add_argument("--model", default=os.environ.get("QDPROPS_MACE_MODEL"))
    p.add_argument("--head", default="omat_pbe")
    p.add_argument("--cif", help="the record's CIF (default: origin.cif in $QDPROPS_CIF_DIRS)")
    p.add_argument("--beam", type=int, default=BEAM)
    p.add_argument("--delta", type=float, default=DELTA)
    p.add_argument("--m-min", type=int, default=1, help="lowest number of MA units")
    p.add_argument("--rows", type=int, help="stop after this many MA removals (test runs)")
    p.add_argument("--descent", choices=("column", "top"), default="column",
                   help="column: every (m, k) loses one MA; top: only the fully passivated node, then ladders")
    p.add_argument("--start", choices=("lattice", "relaxed"), default="lattice",
                   help="relax each product from its lattice positions, or from its parent's relaxed coordinates")
    a = ap.parse_args(argv)
    if a.cmd == "prepare":
        from .cluster import load_cluster
        from .grid import prepare
        print(json.dumps(prepare(a.record, a.out, load_cluster(a.cluster), a.code, time_=a.time, extra=a.extra,
                                 partition=a.partition, pieces=a.pieces, module="qd_workflows.lattice_grid"), indent=1))
        return 0
    from orchestr_ai.qd.run import Settings
    s = Settings(model=a.model, head=a.head, detach_quick=False, detach_screen=False, detach_rigid_refresh=False,
                 detach_frames=False, hessian="analytic")
    g = LatticeGrid(a.record, a.out, s, _find_cif(a.record, a.cif), beam=a.beam, delta=a.delta, m_min=a.m_min,
                    log=lambda x: print(x, flush=True), start=a.start)

    def stop(signum, _frame):
        g._save_cache_throttled(force=True)
        g.log(f"[lattice-grid] signal {signum}: relaxation cache saved, stopping")
        sys.exit(128 + signum)

    import signal
    signal.signal(signal.SIGTERM, stop)
    m0 = g.mk(list(read_xyz(Path(a.record) / "start.xyz")[0]))[0]
    g.run(m_stop=max(a.m_min, m0 - a.rows) if a.rows is not None else None, descent=a.descent)
    return 0


if __name__ == "__main__":
    sys.exit(main())
