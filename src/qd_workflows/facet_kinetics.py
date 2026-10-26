"""
Facet-resolved attachment kinetics of a monomer complex on the shape-library dots of the surface-energy fit: per facet
family, the adsorbed precursor, the incorporated state on lattice sites, the barrier between them (ASE NEB on
MACE-MH-1) and the harmonic free energies from analytic MACE Hessians. The units are added to the builder-built
(ideal) structures, and the dot, the precursor and the incorporated state are all relaxed from there, as the fit's dots
were. CdSe/CdCl2-specific for now (bond length, species, the complex, site densities): to be made general once the
growth model is shown to work. The growth model's facet rate constant follows as

    k_f = rho_f K_f^ads nu_f / (N_A 1e-27)   (A/s),  K_f^ads = e^(-dG_ads/kT) / c0,  nu_f = (kT/h) e^(-dG^/kT)

(attachment flux per area J_f = rho_f (1 - theta_f) K_f^ads C nu_f, v_f = Omega J_f).

Geometry (numpy only, from the ideal builder structure of each dot, whose atoms the relaxed one keeps in order):
  lattice   the four Cd->Se bond vectors b_k of the zinc-blende lattice; anion sites = cation sites + b_k
  sites     next-layer sites on the largest facet of a family: empty cation sites above surface Se, anion sites
            above surface Cd that are empty or hold an X-type Cl (Cl sit on anion lattice sites); each with the
            number of core atoms it bonds to (the anchoring bonds: 2 on a {100} terrace, 3 on {111}A, ...)
  unit      the complex CdSe.CdCl2 (Se, Cd_a, Cd_b, Cl, Cl) placed on lattice sites around the anchor: on an anion
            anchor Se goes there (a Cl on that site is displaced to a free anion site of the new Cd), the two Cd on
            cation sites bonded to it, the Cl on anion sites bonded to the Cd; on a cation anchor Cd_a goes there,
            Se on an anion site bonded to it, Cd_b on a cation site bonded to Se
  precursor the same five atoms moved outward along the facet normal (the facet unchanged): its relaxation is the
            adsorbed state, or falls into the incorporated one (no barrier for that site)
With `--on-relaxed` the units go on the relaxed dots instead, each new atom moved by the local displacement (relaxed -
ideal) of the atoms around it.

    python -m qd_workflows.facet_kinetics input  <ladder_input.json> <state.json> --dots ids.json --out input.json
    python -m qd_workflows.facet_kinetics run    <input.json> --out DIR --model M [--head omat_pbe]
    python -m qd_workflows.facet_kinetics report <DIR>/kinetics.json [--T 500] --out report.json
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from . import shape_ladder as sl

BOND = 2.686          # A, ideal zinc-blende Cd-Se bond (a sqrt(3)/4, a = 6.203 A)
CL_BOND = 2.45        # A, Cd-Cl
CLASH = 2.15          # A: a new atom closer than this to another atom does not fit
OCC = 0.6             # A: a lattice site with an atom this close is occupied
SHIFT = (2.5, 3.0, 3.5)   # A: outward displacement of the precursor (the first that fits)
LOCAL = 7.0           # A: atoms this close to the complex enter the partial Hessian (and stay free in the NEB)
NEB_FREE = 9.0        # A: atoms farther than this from the complex are fixed along the band
TEMPS = (450.0, 500.0, 550.0)
KB = 8.617333262e-5   # eV/K
H_EV = 4.135667696e-15    # eV s
NA = 6.02214076e23


# ---------------- lattice and sites ----------------
def bond_vectors(symbols, positions):
    """The four Cd->Se bond vectors of the ideal structure (A)."""
    s, p = np.asarray(symbols), np.asarray(positions, float)
    cd, se = p[s == "Cd"], p[s == "Se"]
    v = (se[None] - cd[:, None]).reshape(-1, 3)
    v = v[np.abs(np.linalg.norm(v, axis=1) - BOND) < 0.2]
    u = v / np.linalg.norm(v, axis=1)[:, None]
    out = []
    for x in u:
        if not any(np.dot(x, y) > 0.95 for y in out):
            out.append(x)
    if len(out) != 4:
        raise ValueError(f"expected 4 Cd->Se bond directions, found {len(out)}")
    return np.array(out) * BOND


class Surface:
    """Sites around the largest facet of one family on a dot (ideal, core-centred frame)."""

    def __init__(self, symbols, start, relaxed, planes, family):
        self.s = np.asarray(symbols)
        self.p0 = sl.centred(symbols, start)                       # ideal positions (builder)
        self.p1 = sl.centred(symbols, relaxed if relaxed is not None else start)   # where the units go
        self.B = bond_vectors(symbols, start)
        self.fp = sl.fitted_planes(self.s, self.p0, planes)
        self.exp = sl.exposed(self.fp)
        from .shapes import polyhedron
        areas = polyhedron([(n, d, str(i)) for i, (n, d, _l) in enumerate(self.fp)])["areas_A2"]
        mine = [i for i in self.exp if self.fp[i][2] == family]
        if not mine:
            raise ValueError(f"no exposed {family} facet")
        self.plane = max(mine, key=lambda i: areas.get(str(i), 0.0))
        self.n, self.d = self.fp[self.plane][0], self.fp[self.plane][1]
        self.family = family
        self.area = float(areas.get(str(self.plane), 0.0))
        from scipy.spatial import cKDTree
        self.tree = cKDTree(self.p0)
        self.core = np.isin(self.s, ("Cd", "Se"))

    def occupant(self, x, r=OCC):
        hit = self.tree.query_ball_point(x, r)
        return hit[0] if hit else None

    def bonded(self, x, species):
        """Atoms of `species` bonded to a lattice site x."""
        return [j for j in self.tree.query_ball_point(x, BOND + 0.3)
                if self.s[j] in species and abs(np.linalg.norm(self.p0[j] - x) - BOND) < 0.3]

    def inset(self, x):
        """How far a point lies inside the other exposed planes (A): small near the facet's edges."""
        return min((self.fp[i][1] - x @ self.fp[i][0] for i in self.exp if i != self.plane), default=99.0)

    def candidates(self):
        """Next-layer sites over the facet: [{"x", "type", "anchor", "cl", "role", "inset"}]."""
        top = self.core & (self.p0 @ self.n > self.d - 2.0)
        out, seen = [], []
        for j in np.flatnonzero(top):
            vecs = self.B if self.s[j] == "Cd" else -self.B           # Cd -> anion site, Se -> cation site
            for b in vecs:
                if b @ self.n < 0.3:
                    continue
                x = self.p0[j] + b
                if x @ self.n <= self.d + 0.1 or any(np.linalg.norm(x - y) < 0.3 for y in seen):
                    continue
                seen.append(x)
                typ = "anion" if self.s[j] == "Cd" else "cation"
                occ = self.occupant(x)
                if occ is not None and not (typ == "anion" and self.s[occ] == "Cl"):
                    continue
                partner = ("Cd",) if typ == "anion" else ("Se",)
                anchor = len([k for k in self.bonded(x, partner) if self.core[k]])
                out.append({"x": x, "type": typ, "anchor": anchor, "cl": occ,
                            "role": sl.role_of(x, self.fp, self.exp), "inset": float(self.inset(x))})
        return out

    def pick(self, n_terrace=2, n_edge=1):
        """[(site, states)]: terrace sites (facet role, most anchoring bonds, farthest from the edges, 6 A apart) and
        edge sites, tried in that order until enough units fit."""
        c = self.candidates()
        fam = f"facet:{self.family}"
        terr = [x for x in c if x["role"] == fam]
        edge = [x for x in c if x["role"].startswith("edge:") and self.family in x["role"].split(":")[1].split("|")]
        out = []
        for cls, group, n, key in (("terrace", terr, n_terrace, lambda x: (-x["anchor"], -x["inset"])),
                                   ("edge", edge, n_edge, lambda x: (-x["anchor"], x["x"] @ self.n))):
            got = []
            for x in sorted(group, key=key):
                if len(got) >= n:
                    break
                if any(np.linalg.norm(x["x"] - y["x"]) < 6.0 for y, _st in got):
                    continue
                st = self.states(x)
                if st is not None:
                    got.append(({**x, "class": cls}, st))
            out += got
        return out

    # -- the unit on lattice sites --
    def _free(self, x, taken, ignore=()):
        """x is empty in the ideal frame and clear of every atom and of the new atoms."""
        for j in self.tree.query_ball_point(x, CLASH):
            if j not in ignore:
                return False
        return all(np.linalg.norm(x - y) >= CLASH for y in taken)

    def _bonds(self, x, species, new, ign):
        """Bonds a lattice site x would make: atoms of `species` already there (ideal frame) plus the new ones."""
        old = [j for j in self.bonded(x, species) if j not in ign]
        return len(old) + sum(1 for sp, y in new if sp in species and abs(np.linalg.norm(y - x) - BOND) < 0.3)

    def _best(self, cands, species, new, ign):
        """The free candidate site with the most bonds (epitaxial), then the most outward."""
        taken = [y for _sp, y in new]
        ok = [v for v in cands if self._free(v, taken, ign)]
        if not ok:
            return None
        return max(ok, key=lambda v: (self._bonds(v, species, new, ign), float(v @ self.n)))

    def unit(self, site):
        """Positions (ideal frame) of [Se, Cd_a, Cd_b, Cl, Cl] and, when the anchor held a Cl, (index, new position)
        of that Cl; None when the unit does not fit. Every new atom takes the free lattice site with the most bonds to
        the atoms already there (dot and unit), so the unit grows the crystal epitaxially."""
        x, ign = site["x"], ({site["cl"]} if site["cl"] is not None else set())
        anions, cations = ("Se", "Cl"), ("Cd",)
        new = []
        if site["type"] == "anion":
            new.append(("Se", x))
            cd_a = self._best([x - b for b in self.B], anions, new, ign)
            if cd_a is None:
                return None
            new.append(("Cd", cd_a))
        else:
            new.append(("Cd", x))
            se = self._best([x + b for b in self.B], cations, new, ign)
            if se is None:
                return None
            new.insert(0, ("Se", se))
        se = new[0][1]
        # Cd_b: a cation site bonded to the new Se, else one next to Cd_a bonded to an anion already there
        cd_b = self._best([se - b for b in self.B], anions, new, ign)
        if cd_b is None:
            cd_a = new[1][1]
            cd_b = self._best([cd_a + b - b2 for b in self.B for b2 in self.B if np.linalg.norm(b - b2) > 0.1],
                              anions, new, ign)
        if cd_b is None:
            return None
        new.append(("Cd", cd_b))
        need = 3 if site["cl"] is not None else 2
        cls = []
        for _ in range(need):                                          # Cl on anion sites bonded to the new Cd
            c = self._best([cd + b for _sp, cd in new[1:3] for b in self.B], cations, new + cls, ign)
            if c is None or self._bonds(c, cations, new[1:3], set(range(len(self.s)))) == 0:
                break
            cls.append(("Cl", c))
        while len(cls) < need:                                         # along the normal from a new Cd
            cd = new[1 + len(cls) % 2][1]
            v = cd + CL_BOND * self.n
            if not self._free(v, [y for _sp, y in new + cls], ign):
                return None
            cls.append(("Cl", v))
        moved = (site["cl"], cls[2][1]) if site["cl"] is not None else None
        return np.array([new[0][1], new[1][1], new[2][1], cls[0][1], cls[1][1]]), moved

    def extended(self, symbols_add, positions_add, moved=None):
        """A copy with atoms added (ideal frame; builder geometry only), the same facet plane and lattice."""
        import copy
        from scipy.spatial import cKDTree
        o = copy.copy(self)
        o.s = np.concatenate([self.s, np.asarray(symbols_add)])
        o.p0 = np.vstack([self.p0, np.asarray(positions_add, float)])
        if moved is not None:
            o.p0[moved[0]] = moved[1]
        o.p1 = o.p0
        o.tree = cKDTree(o.p0)
        o.core = np.isin(o.s, ("Cd", "Se"))
        return o

    def island(self, sizes):
        """{N: (symbols, positions)}: compact islands of N units in the first new layer of the facet, grown one unit
        at a time from the best terrace site, each next anchor in the same layer next to the island with the most
        bonds (nested: the island of N + 1 contains that of N). Builder geometry only."""
        start = self.pick(1, 0)
        if not start:
            return {}
        site0 = start[0][0]
        h = float(site0["x"] @ self.n)
        cur, atoms, out = self, [], {}
        for k in range(1, max(sizes) + 1):
            if k == 1:
                tries = [site0]
            else:
                isl = np.array(atoms)
                cands = [c for c in cur.candidates() if abs(c["x"] @ self.n - h) < 0.5
                         and c["role"] == f"facet:{self.family}" and np.linalg.norm(isl - c["x"], axis=1).min() < 5.0]
                cen = isl.mean(0)
                tries = sorted(cands, key=lambda c: (-c["anchor"], float(np.linalg.norm(c["x"] - cen))))
            for site in tries:
                u = cur.unit(site)
                if u is not None:
                    break
            else:
                break
            five, moved = u
            cur = cur.extended(UNIT_SYMBOLS, five, moved)
            atoms += list(five[:3])
            if k in sizes:
                out[k] = ([str(x) for x in cur.s], cur.p0.copy())
        return out

    def to_relaxed(self, y):
        """Ideal-frame point -> relaxed frame, by the mean displacement of the atoms within 4 A."""
        y = np.atleast_2d(y)
        out = []
        for v in y:
            nb = self.tree.query_ball_point(v, 4.5)
            out.append(v + (self.p1[nb] - self.p0[nb]).mean(0) if nb else v)
        return np.array(out)

    def states(self, site):
        """(precursor positions, incorporated positions) of dot + [Se, Cd, Cd, Cl, Cl], relaxed frame, same order;
        None when the unit does not fit."""
        u = self.unit(site)
        if u is None:
            return None
        five, moved = u
        five_r = self.to_relaxed(five)
        inc = np.vstack([self.p1, five_r])
        if moved is not None:
            inc[moved[0]] = self.to_relaxed(moved[1])[0]
        from scipy.spatial import cKDTree
        if cKDTree(inc).query(inc, k=2)[0][:, 1].min() < CLASH - 0.05:  # the relaxed frame squeezed a contact
            return None
        t = cKDTree(self.p1)
        for sh in SHIFT:
            pre5 = five_r + sh * self.n
            if all(not t.query_ball_point(v, CLASH) for v in pre5):
                return np.vstack([self.p1, pre5]), inc
        return None


UNIT_SYMBOLS = ["Se", "Cd", "Cd", "Cl", "Cl"]


def make_input(ladder_input, state, dots, out, n_terrace=2, n_edge=1, on_relaxed=False, islands=()):
    """Precursor and incorporated geometries of every picked site on the dots {id: family}: on the builder-built
    structures (default), or on the relaxed ones from the ladder run's state.json; and islands of `islands` units."""
    inp = {d["id"]: d for d in json.loads(Path(ladder_input).read_text())["dots"]}
    rel = json.loads(Path(state).read_text())["relaxed"] if on_relaxed else {}
    rows, isl = [], []
    for i, fam in dots.items():
        d = inp[i]
        surf = Surface(d["symbols"], d["positions"], rel[i]["positions"] if on_relaxed else None, d["planes"], fam)
        picked = surf.pick(n_terrace, n_edge)
        for k, (site, (pre, inc)) in enumerate(picked):
            rows.append({"id": f"{i}|{fam}|{site['class']}{k}", "dot": i, "family": fam, "class": site["class"],
                         "role": site["role"], "anchor_type": site["type"], "anchor_bonds": site["anchor"],
                         "displaced_cl": site["cl"], "facet_area_A2": surf.area, "normal": surf.n.tolist(),
                         "symbols": list(d["symbols"]) + UNIT_SYMBOLS,
                         "precursor": pre.round(5).tolist(), "incorporated": inc.round(5).tolist(),
                         "dot_symbols": list(d["symbols"]), "dot_positions": np.asarray(surf.p1).round(5).tolist(),
                         "unit": list(range(len(d["symbols"]), len(d["symbols"]) + 5))})
        print(f"{i} {fam}: facet {surf.area:.0f} A^2, {len(surf.candidates())} next-layer sites, "
              f"{sum(1 for r in rows if r['dot'] == i)} states")
        if islands and not on_relaxed:
            for n, (sy, ps) in surf.island(islands).items():
                isl.append({"id": f"{i}|{fam}|island{n}", "dot": i, "family": fam, "n_units": n, "symbols": sy,
                            "positions": np.asarray(ps).round(5).tolist(), "dot_symbols": list(d["symbols"]),
                            "dot_positions": np.asarray(surf.p1).round(5).tolist()})
            print(f"   islands: {sorted(x['n_units'] for x in isl if x['dot'] == i)}")
    Path(out).write_text(json.dumps({"sites": rows, "islands": isl}))
    return rows


def cluster_input(thermo_path, complex_positions, out, sizes=range(1, 13), n_dirs=2, mu_mx=None, gap=7.0):
    """Attachment of the monomer complex to the grid's clusters, for the Becker-Doring barrier E_att(m): for each m,
    the lowest state of the shell k that minimises G(m,k) - k mu_MX (mu_mx; default: the best k at 500 K, mu_MX =
    the CdCl2 reference), the complex placed at `n_dirs` under-coordinated Cd, far (every contact beyond `gap` A, past
    the MACE cutoff) and bound (its Se 2.7 A from the Cd, pointing in). Rows in the format of make_input."""
    from .xyz import read_xyz
    th = json.loads(Path(thermo_path).read_text())
    base = Path(thermo_path).resolve().parent
    cx = np.asarray(complex_positions, float)
    cx = cx - cx[0]                                                    # Se of the complex at the origin
    u = cx[1:].mean(0)
    u /= np.linalg.norm(u)                                             # Se -> rest of the complex
    nodes = {(n["m"], n["k"]): n for n in th["nodes"]}
    ref = th.get("references", {}).get("CdCl2", {}).get("G_eV")
    mmx = mu_mx if mu_mx is not None else (ref[len(ref) // 2] if ref else 0.0)
    rows = []
    for m in sizes:
        ks = [k for (mm, k) in nodes if mm == m]
        if not ks:
            continue
        g = lambda k: min(s["E_eV"] for s in nodes[(m, k)]["states"]) - k * mmx    # noqa: E731
        k = min(ks, key=g)
        st = min(nodes[(m, k)]["states"], key=lambda s: s["E_eV"])
        sym, pts = read_xyz(base / st["file"])
        sym, pts = list(sym), np.asarray(pts, float) - np.asarray(pts, float).mean(0)
        from scipy.spatial import cKDTree
        t = cKDTree(pts)
        cds = [i for i, s in enumerate(sym) if s == "Cd"]
        cn = {i: len(t.query_ball_point(pts[i], 3.0)) - 1 for i in cds}
        order = sorted(cds, key=lambda i: (cn[i], -np.linalg.norm(pts[i])))
        picked = []
        for i in order:
            if len(picked) >= n_dirs:
                break
            if all(np.linalg.norm(pts[i] - pts[j]) > 4.0 for j in picked):
                picked.append(i)
        for j, i in enumerate(picked):
            d = pts[i] / max(np.linalg.norm(pts[i]), 1e-9)
            # rotate the complex so that Se -> rest points along d (outward)
            v = np.cross(u, d)
            s_, c_ = np.linalg.norm(v), float(u @ d)
            if s_ < 1e-8:
                rot = np.eye(3) if c_ > 0 else -np.eye(3)
            else:
                vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
                rot = np.eye(3) + vx + vx @ vx * ((1 - c_) / s_ ** 2)
            five = cx @ rot.T
            bound = five + pts[i] + 2.7 * d
            if cKDTree(pts).query(bound)[0].min() < 2.1:
                continue
            far = bound.copy()
            while cKDTree(pts).query(far)[0].min() < gap:
                far = far + 0.5 * d
            rows.append({"id": f"cluster|m{m}k{k}|dir{j}", "dot": f"cluster_m{m}_k{k}", "family": f"m{m}",
                         "class": "attach", "role": "cluster", "anchor_type": "cation", "anchor_bonds": cn[i],
                         "displaced_cl": None, "facet_area_A2": 0.0, "normal": d.tolist(),
                         "symbols": sym + UNIT_SYMBOLS, "precursor": np.vstack([pts, far]).round(5).tolist(),
                         "incorporated": np.vstack([pts, bound]).round(5).tolist(), "dot_symbols": sym,
                         "dot_positions": pts.round(5).tolist(), "unit": list(range(len(sym), len(sym) + 5)),
                         "m": m, "k": k})
    Path(out).write_text(json.dumps({"sites": rows, "islands": []}))
    print(f"{len(rows)} attachment paths over m = {sorted(set(r['m'] for r in rows))}")
    return rows


# ---------------- run (MACE, on the cluster) ----------------
class Kinetics:
    def __init__(self, inp, out, settings, log=print, neb_images=7, fmax=0.05, steps=600):
        from orchestr_ai.qd.steps import detachment as dt
        self.dt, self.s, self.log = dt, settings, log
        data = json.loads(Path(inp).read_text())
        self.sites, self.islands = data["sites"], data.get("islands", [])
        self.out = Path(out)
        self.out.mkdir(parents=True, exist_ok=True)
        self.cache = dt.RelaxCache(SimpleNamespace(props=self.out, settings=settings, record_dir=self.out))
        self.nimg, self.fmax, self.steps = neb_images, fmax, steps
        sp = self.out / "state.json"
        self.state = json.loads(sp.read_text()) if sp.is_file() else {}

    def save(self):
        tmp = self.out / "state.tmp"
        tmp.write_text(json.dumps(self.state))
        os.replace(tmp, self.out / "state.json")
        self.cache._save()

    def relax(self, items):
        if not items:
            return []
        res = self.cache.relax_many(items, self.s)
        res = self.cache.relax_many([(r[0].get_chemical_symbols(), r[0].get_positions()) for r in res], self.s, polish=True)
        return [(r[0].get_positions(), float(r[1])) for r in res]

    def calc(self):
        if not hasattr(self, "_calc"):
            from orchestr_ai.qd.engines import mace_calculator
            self._calc = mace_calculator(self.s.head, self.s.model, self.s.device, self.s.dtype)
        return self._calc

    def complex_energy(self):
        """The free complex CdSe.CdCl2, relaxed from a chain and a bent start (the lower kept)."""
        if "complex" in self.state:
            return self.state["complex"]
        starts = [np.array([[0, 0, 0], [2.5, 0, 0], [-2.5, 0, 0], [4.9, 0, 0], [-4.9, 0, 0]], float),
                  np.array([[0, 0, 0], [2.0, 1.6, 0], [-2.0, 1.6, 0], [4.1, 2.8, 0], [-4.1, 2.8, 0]], float)]
        res = self.relax([(UNIT_SYMBOLS, p) for p in starts])
        p, e = min(res, key=lambda r: r[1])
        self.state["complex"] = {"positions": p.round(6).tolist(), "E": e}
        self.save()
        return self.state["complex"]

    def neb(self, symbols, a, b, free):
        """Climbing-image NEB (ASE, IDPP start, FIRE) between relaxed end points; atoms outside `free` fixed.
        Orchestr.AI's qd NEB step when it is there, else the same band here."""
        fixed = [i for i in range(len(symbols)) if i not in set(free)]
        try:
            from orchestr_ai.qd.steps.neb import neb_path
        except ImportError:
            neb_path = None
        if neb_path is not None:
            p = neb_path(symbols, a, b, self.s, images=self.nimg, fixed=fixed, fmax=self.fmax, steps=self.steps,
                         calc=self.calc())
            return {"energies": p["energies"], "ts": np.round(p["saddle_positions"], 6).tolist(),
                    "ts_index": p["saddle_index"], "converged": p["converged"], "max_force": p["max_force"]}
        from ase import Atoms
        from ase.constraints import FixAtoms
        from ase.mep import NEB
        from ase.optimize import FIRE
        b = np.array(b, float)
        b[fixed] = np.asarray(a, float)[fixed]                         # ASE refuses fixed atoms that move
        images = []
        for k in range(self.nimg + 2):
            at = Atoms(symbols, positions=a if k < self.nimg + 1 else b)
            at.set_constraint(FixAtoms(indices=fixed))
            at.calc = self.calc()
            images.append(at)
        band = NEB(images, climb=False, allow_shared_calculator=True, method="improvedtangent")
        band.interpolate(method="idpp")
        FIRE(band, logfile=None).run(fmax=0.15, steps=self.steps // 2)
        band.climb = True
        conv = FIRE(band, logfile=None).run(fmax=self.fmax, steps=self.steps)
        es = [float(im.get_potential_energy()) for im in images]
        top = int(np.argmax(es))
        return {"energies": es, "ts": images[top].get_positions().round(6).tolist(), "ts_index": top,
                "converged": bool(conv)}

    def hessian_modes(self, symbols, positions, local):
        """Partial-Hessian frequencies (cm^-1, negative = imaginary) of the atoms `local`, the rest held fixed."""
        from ase import Atoms
        from orchestr_ai.qd.engines import mace_hessian
        at = Atoms(symbols, positions=positions)
        h = mace_hessian(at, self.calc())
        if not len(local):
            return []
        idx = np.array([3 * i + k for i in local for k in range(3)], dtype=int)
        m = np.repeat(at.get_masses()[local], 3)
        hl = h[np.ix_(idx, idx)] / np.sqrt(np.outer(m, m))
        w = np.linalg.eigvalsh(0.5 * (hl + hl.T))
        conv = 521.47090                                   # sqrt(eV/A^2/amu) -> cm^-1
        return (np.sign(w) * np.sqrt(np.abs(w)) * conv).tolist()

    @staticmethod
    def _far_moved(a, b, u):
        """Atoms away from the unit (beyond NEB_FREE in both end points) that moved more than 1 A between them."""
        da = np.min(np.linalg.norm(a[:, None] - a[u][None], axis=2), axis=1)
        db = np.min(np.linalg.norm(b[:, None] - b[u][None], axis=2), axis=1)
        shift = np.linalg.norm(a - b, axis=1)
        return [int(i) for i in np.flatnonzero((shift > 1.0) & (da >= NEB_FREE) & (db >= NEB_FREE))]

    def match_shells(self):
        """Where the dot's CdCl2 shell relaxed into different arrangements in the two end points (atoms far from the
        unit moved), also relax a precursor built on the incorporated state's shell (unit pulled out along the
        normal) and an incorporated state built on the precursor's shell; the NEB then uses the pair whose shells
        agree best (st["pair"]). The lowest of each kind still sets the binding energies."""
        from scipy.spatial import cKDTree
        jobs = []
        for r in self.sites:
            st = self.state[r["id"]]
            if "neb" in st or "pair" in st:
                continue
            a, b = np.array(st["precursor"]["positions"]), np.array(st["incorporated"]["positions"])
            u, n, dcl = r["unit"], np.asarray(r["normal"], float), r["displaced_cl"]
            if not self._far_moved(a, b, u):
                st["pair"] = ["precursor", "incorporated"]
                continue
            inc_alt = a.copy()
            inc_alt[u] = b[u]
            if dcl is not None:
                inc_alt[dcl] = b[dcl]
            jobs.append((r, "incorporated_alt", inc_alt))
            rest = [i for i in range(len(b)) if i not in set(u)]
            for sh in SHIFT:
                pre_alt = b.copy()
                pre_alt[u] = b[u] + sh * n
                if dcl is not None:
                    pre_alt[dcl] = a[dcl]
                if cKDTree(pre_alt[rest]).query(pre_alt[u])[0].min() >= CLASH:
                    jobs.append((r, "precursor_alt", pre_alt))
                    break
        if jobs:
            for (r, key, _x), (p, e) in zip(jobs, self.relax([(r["symbols"], x) for r, _k, x in jobs])):
                self.state[r["id"]][key] = {"positions": p.round(6).tolist(), "E": e}
        for r in self.sites:
            st = self.state[r["id"]]
            if "pair" in st:
                continue
            pres = [k for k in ("precursor", "precursor_alt") if k in st]
            incs = [k for k in ("incorporated", "incorporated_alt") if k in st]
            u = r["unit"]
            pairs = [(len(self._far_moved(np.array(st[p]["positions"]), np.array(st[i]["positions"]), u)),
                      st[p]["E"] + st[i]["E"], p, i) for p in pres for i in incs]
            best = min(pairs)
            st["pair"] = [best[2], best[3]]
            st["pair_far_moved"] = best[0]
            self.log(f"[kin] {r['id']}: shells matched with {best[2]} + {best[3]} ({best[0]} far atoms still moved)")
        self.save()

    def run(self):
        t0 = time.time()
        cx = self.complex_energy()
        self.log(f"[kin] complex E = {cx['E']:.4f} eV")
        # dots and end points, one batch
        dots = {r["dot"]: r for r in self.sites + self.islands}
        todo = [i for i in dots if f"dot|{i}" not in self.state]
        for i, (p, e) in zip(todo, self.relax([(dots[i]["dot_symbols"], dots[i]["dot_positions"]) for i in todo])):
            self.state[f"dot|{i}"] = {"positions": p.round(6).tolist(), "E": e}
        jobs = []
        for r in self.sites:
            st = self.state.setdefault(r["id"], {})
            for key in ("precursor", "incorporated"):
                if key not in st:
                    jobs.append((r, key, r[key]))
        for r in self.islands:
            if "island" not in self.state.setdefault(r["id"], {}):
                jobs.append((r, "island", r["positions"]))
        for (r, key, _x), (p, e) in zip(jobs, self.relax([(r["symbols"], x) for r, _k, x in jobs])):
            self.state[r["id"]][key] = {"positions": p.round(6).tolist(), "E": e}
        self.save()
        self.log(f"[kin] relaxed {len(todo)} dots and {len(jobs)} end points and islands ({time.time() - t0:.0f} s)")
        self.match_shells()
        # NEBs: every pending band at once when Orchestr.AI has neb_many (one batched MACE call per step), else in turn
        pend = []
        for r in self.sites:
            st = self.state[r["id"]]
            if "neb" in st:
                continue
            kp, ki = st.get("pair", ["precursor", "incorporated"])
            a, b = np.array(st[kp]["positions"]), np.array(st[ki]["positions"])
            u = r["unit"]
            moved = float(np.sqrt(((a[u] - b[u]) ** 2).sum(1).mean()))
            st["unit_rmsd"] = moved
            if moved < 0.5 or st[kp]["E"] <= st[ki]["E"] - 2.0:
                st["neb"] = {"skipped": f"end points alike (unit RMSD {moved:.2f} A)" if moved < 0.5 else "precursor far below"}
                continue
            # free: atoms near the unit, and any atom that moved between the end points (a shell that rearranged on
            # relaxation far from the site must move along the band too; fixing it would distort the end image)
            shift = np.linalg.norm(a - b, axis=1)
            free = set(i for i in range(len(a)) if shift[i] > 0.3 or min(np.linalg.norm(a[u] - a[i], axis=1).min(),
                                                                         np.linalg.norm(b[u] - b[i], axis=1).min()) < NEB_FREE)
            far = [i for i in range(len(a)) if shift[i] > 1.0 and min(np.linalg.norm(a[u] - a[i], axis=1).min(),
                                                                      np.linalg.norm(b[u] - b[i], axis=1).min()) >= NEB_FREE]
            st["shell_rearranged"] = len(far)
            part = st.get("neb_partial") or {}
            pend.append((r, {"symbols": r["symbols"], "initial": a, "final": b,
                             "fixed": [i for i in range(len(a)) if i not in free],
                             "band": part.get("band"), "stage": part.get("stage")}))
        try:
            from orchestr_ai.qd.steps.neb import neb_many
        except ImportError:
            neb_many = None
        if pend and neb_many is not None:
            self.log(f"[kin] {len(pend)} NEBs in lockstep ({sum(1 for _r, j in pend if j['band'])} resumed)")

            def checkpoint(parts):
                for (r, _j), part in zip(pend, parts):
                    self.state[r["id"]]["neb_partial"] = part
                self.save()

            try:
                paths = neb_many([j for _r, j in pend], self.s, images=self.nimg, fmax=self.fmax, steps=self.steps,
                                 max_atoms=4000, checkpoint=checkpoint)
            except TypeError:                               # an Orchestr.AI without checkpoints
                paths = neb_many([j for _r, j in pend], self.s, images=self.nimg, fmax=self.fmax, steps=self.steps,
                                 max_atoms=4000)
            for (r, _j), p in zip(pend, paths):
                self.state[r["id"]]["neb"] = {"energies": p["energies"], "ts": np.round(p["saddle_positions"], 6).tolist(),
                                              "ts_index": p["saddle_index"], "converged": p["converged"],
                                              "max_force": p["max_force"]}
                self.state[r["id"]].pop("neb_partial", None)
        else:
            for r, j in pend:
                self.state[r["id"]]["neb"] = self.neb(r["symbols"], j["initial"], j["final"],
                                                      [i for i in range(len(j["initial"])) if i not in set(j["fixed"])])
                self.save()
        self.save()
        for r in self.sites:
            st = self.state[r["id"]]
            self.log(f"[kin] {r['id']}: E_ads {st['precursor']['E'] - self.state['dot|' + r['dot']]['E'] - cx['E']:+.3f}, "
                     f"E_inc {st['incorporated']['E'] - self.state['dot|' + r['dot']]['E'] - cx['E']:+.3f} eV, "
                     f"barrier {max(st['neb'].get('energies', [0])) - st['precursor']['E'] if 'energies' in st['neb'] else 0:.3f} eV "
                     f"({time.time() - t0:.0f} s)")
        # partial Hessians: precursor, saddle, incorporated, the dot's same atoms, the free complex
        if "complex_modes" not in self.state:
            self.state["complex_modes"] = self.hessian_modes(UNIT_SYMBOLS, cx["positions"], list(range(5)))
        for r in self.sites:
            st = self.state[r["id"]]
            if "modes" in st:
                continue
            kp, ki = st.get("pair", ["precursor", "incorporated"])
            a = np.array(st[kp]["positions"])
            u = r["unit"]
            b_ = np.array(st[ki]["positions"])                 # near the unit in either end point (a far precursor
            local = [i for i in range(len(a)) if min(np.linalg.norm(a[u] - a[i], axis=1).min(),     # has none)
                                                     np.linalg.norm(b_[u] - b_[i], axis=1).min()) < LOCAL]
            dot_local = [i for i in local if i < u[0]]
            modes = {"local": local, "pair": [kp, ki],
                     "precursor": self.hessian_modes(r["symbols"], a, local),
                     "incorporated": self.hessian_modes(r["symbols"], st[ki]["positions"], local),
                     "dot": self.hessian_modes(r["dot_symbols"], self.state["dot|" + r["dot"]]["positions"], dot_local)}
            if "ts" in st["neb"]:
                modes["ts"] = self.hessian_modes(r["symbols"], st["neb"]["ts"], local)
            st["modes"] = modes
            self.save()
            self.log(f"[kin] {r['id']}: Hessians ({time.time() - t0:.0f} s)")
        # islands: partial Hessians of the island and its surroundings, and of the same atoms of the bare dot
        for r in self.islands:
            st = self.state.get(r["id"], {})
            if "island" not in st or "modes" in st:
                continue
            p = np.array(st["island"]["positions"])
            n0 = len(r["dot_symbols"])
            new = list(range(n0, len(p)))
            local = [i for i in range(len(p)) if np.linalg.norm(p[new] - p[i], axis=1).min() < LOCAL]
            st["modes"] = {"local": local,
                           "island": self.hessian_modes(r["symbols"], p, local),
                           "dot": self.hessian_modes(r["dot_symbols"], self.state["dot|" + r["dot"]]["positions"],
                                                     [i for i in local if i < n0])}
            self.save()
            self.log(f"[kin] {r['id']}: island Hessians ({time.time() - t0:.0f} s)")
        self.write()
        self.log("[kin] done")

    def write(self):
        cx = self.state["complex"]
        rows = []
        for r in self.sites:
            st = self.state[r["id"]]
            e_dot = self.state["dot|" + r["dot"]]["E"]
            row = {k: r[k] for k in ("id", "dot", "family", "class", "role", "anchor_type", "anchor_bonds",
                                     "displaced_cl", "facet_area_A2")}
            kp, ki = st.get("pair", ["precursor", "incorporated"])
            row.update(E_dot=e_dot, E_complex=cx["E"],
                       E_pre=min(st[k]["E"] for k in ("precursor", "precursor_alt") if k in st),
                       E_inc=min(st[k]["E"] for k in ("incorporated", "incorporated_alt") if k in st),
                       pair=[kp, ki], E_pair_pre=st[kp]["E"], E_pair_inc=st[ki]["E"],
                       pair_far_moved=st.get("pair_far_moved", 0),
                       unit_rmsd=st.get("unit_rmsd"), neb={k: v for k, v in st.get("neb", {}).items() if k != "ts"},
                       modes=st.get("modes"))
            rows.append(row)
        isl = [{"id": r["id"], "dot": r["dot"], "family": r["family"], "n_units": r["n_units"],
                "E": self.state[r["id"]]["island"]["E"], "E_dot": self.state["dot|" + r["dot"]]["E"],
                "dE": self.state[r["id"]]["island"]["E"] - self.state["dot|" + r["dot"]]["E"] - r["n_units"] * cx["E"],
                "modes": self.state[r["id"]].get("modes")}
               for r in self.islands if "island" in self.state.get(r["id"], {})]
        (self.out / "kinetics.json").write_text(json.dumps({"complex": cx, "complex_modes": self.state.get("complex_modes"),
                                                            "sites": rows, "islands": isl}, indent=1))


# ---------------- implicit solvation (no CP2K) ----------------
def solvation(state_path, input_path, eps: float = 2.0) -> dict:
    """Implicit-solvent correction of every site and island from the geometry alone (grid_solvation's model):
    a dot (with or without the unit) gets per_area(eps) x S, S its solvent-exposed area (vdW spheres x 1.4); the free
    complex its campaign COSMO data (the grid's m = 1, k = 1 row). Returns {id: {key: ddG_solv (eV)}} with
    ddG_solv(state) = dG_solv(dot + units) - dG_solv(dot) - n_units dG_solv(complex), keys precursor / incorporated /
    ts for sites and island for islands."""
    from .grid_solvation import SOLVATION, cosmo, per_area, surface
    st = json.loads(Path(state_path).read_text())
    inp = json.loads(Path(input_path).read_text())
    pa = per_area(eps)
    g_cx = cosmo(SOLVATION["complexes"]["1"], eps)
    s_dot, out = {}, {}

    def dot_s(r):
        if r["dot"] not in s_dot:
            s_dot[r["dot"]] = surface(r["dot_symbols"], st["dot|" + r["dot"]]["positions"])
        return s_dot[r["dot"]]

    for r in inp.get("sites", []):
        s = st.get(r["id"], {})
        if not s:
            continue
        kp, ki = s.get("pair", ["precursor", "incorporated"])
        res = {}
        for key, pos in (("precursor", s.get(kp, {}).get("positions")), ("incorporated", s.get(ki, {}).get("positions")),
                         ("ts", s.get("neb", {}).get("ts"))):
            if pos is not None:
                res[key] = pa * (surface(r["symbols"], pos) - dot_s(r)) - g_cx
        out[r["id"]] = res
    for r in inp.get("islands", []):
        s = st.get(r["id"], {})
        if "island" in s:
            out[r["id"]] = {"island": pa * (surface(r["symbols"], s["island"]["positions"]) - dot_s(r)) - r["n_units"] * g_cx}
    return {"eps": eps, "per_area_eV_A2": pa, "complex_eV": g_cx, "terms": out}


# ---------------- report ----------------
def vib_free(freqs_cm1, T, n_free=0):
    """Harmonic vibrational free energy (eV) of real modes above 50 cm^-1 (lower ones held at 50: quasi-harmonic).
    n_free: the molecule is free (gas or solution), its n_free smallest |modes| are translations and rotations
    (6, or 5 if linear) and are left out: gas_trans_rot counts them."""
    fr = sorted(freqs_cm1, key=abs)[int(n_free):]
    f = np.asarray([x for x in fr if x > 0], float)
    f = np.maximum(f, 50.0) * 1.239841984e-4
    kt = KB * T
    return float(np.sum(0.5 * f + kt * np.log1p(-np.exp(-f / kt))))


def gas_trans_rot(positions, symbols, T, c0=1.0):
    """Translational (1 M standard state) and rotational free energy (eV) of the free complex (rigid rotor)."""
    from ase import Atoms
    at = Atoms(symbols, positions=positions)
    m = at.get_masses().sum() * 1.66053907e-27
    kt_j = 1.380649e-23 * T
    lam = 6.62607015e-34 / math.sqrt(2 * math.pi * m * kt_j)
    v = 1.0 / (c0 * NA * 1e3)                                 # m^3 per molecule at c0 (M)
    g_trans = -KB * T * math.log(v / lam ** 3)
    I = np.asarray(at.get_moments_of_inertia()) * 1.66053907e-47     # principal moments, amu A^2 -> kg m^2
    I = I[I > 1e-50]
    hb = 1.054571817e-34
    if len(I) == 2:                                           # linear
        q = 2 * I.mean() * kt_j / hb ** 2 / 2
    else:
        q = math.sqrt(math.pi * np.prod(2 * I * kt_j / hb ** 2)) / 2
    return g_trans - KB * T * math.log(q)


def report(path, T=500.0, site_density=None):
    """Per site: dE_ads, dE_inc, the barrier, the harmonic dG_ads and dG^ at T and the rate constant k_f (A/s)."""
    k = json.loads(Path(path).read_text())
    cx = k["complex"]
    kt = KB * T
    g_cx = vib_free(k["complex_modes"], T, n_free=6) + gas_trans_rot(cx["positions"], UNIT_SYMBOLS, T)
    rho = site_density or {"100": 1 / 19.24, "110": 1 / 27.21, "111A": 1 / 16.66, "111B": 1 / 16.66}   # sites/A^2
    out = []
    for r in k["sites"]:
        dE_ads = r["E_pre"] - r["E_dot"] - r["E_complex"]              # the lowest precursor found
        dE_inc = r["E_inc"] - r["E_dot"] - r["E_complex"]              # the lowest incorporated state found
        dE_pair = r.get("E_pair_pre", r["E_pre"]) - r["E_dot"] - r["E_complex"]   # the band's own start
        neb = r["neb"]
        dE_ts = (max(neb["energies"]) - neb["energies"][0]) if "energies" in neb else 0.0
        row = {**{x: r[x] for x in ("id", "family", "class", "role", "anchor_bonds")},
               "dE_ads": dE_ads, "dE_inc": dE_inc, "dE_ads_band": dE_pair, "dE_barrier": dE_ts,
               "pair": r.get("pair"), "pair_far_moved": r.get("pair_far_moved", 0),
               "neb": neb.get("skipped") or ("converged" if neb.get("converged") else "not converged")}
        m = r.get("modes")
        if m and r["family"] in rho:                     # facet sites (cluster-attachment rows: energies only)
            f_pre, f_dot = vib_free(m["precursor"], T), vib_free(m["dot"], T)
            dG_ads = dE_pair + f_pre - f_dot - g_cx
            f_ts = vib_free(m["ts"], T) if "ts" in m else f_pre
            dG_ts = dE_ts + f_ts - f_pre
            n_imag = sum(1 for x in m.get("ts", []) if x < -50)
            nu = kt / H_EV * math.exp(-max(dG_ts, 0.0) / kt)
            K_ads = math.exp(-dG_ads / kt)                    # M^-1 (c0 = 1 M)
            k_tst = rho[r["family"]] * K_ads * nu / (NA * 1e-27)
            # no faster than the complexes hit the surface: the collision limit sqrt(kT / 2 pi m) (A/s)
            m_cx = sum({"Se": 78.971, "Cd": 112.414, "Cl": 35.453}[s] for s in UNIT_SYMBOLS) * 1.66053907e-27
            k_coll = math.sqrt(1.380649e-23 * T / (2 * math.pi * m_cx)) * 1e10
            row.update(dG_ads=dG_ads, dG_barrier=dG_ts, ts_imaginary_modes=n_imag, nu=nu, K_ads_per_M=K_ads,
                       k_tst=k_tst, k_collision=k_coll, k_f=min(k_tst, k_coll))
        out.append(row)
    # islands: energy of N units in the first new layer against N free complexes; the step to N from N - 1
    isl = sorted(k.get("islands", []), key=lambda r: (r["family"], r["n_units"]))
    prev = {}
    for r in isl:
        p = prev.get(r["family"], (0, 0.0))
        r["dE_step"] = (r["dE"] - p[1]) / max(r["n_units"] - p[0], 1)       # per unit added since the last size
        prev[r["family"]] = (r["n_units"], r["dE"])
        m = r.pop("modes", None)
        if m and m.get("island"):
            # island free energy against N free complexes at 1 M: Delta G_f(N), whose maximum is the layer barrier
            r["dG"] = r["dE"] + vib_free(m["island"], T) - vib_free(m["dot"], T) - r["n_units"] * g_cx
    for fam in sorted({r["family"] for r in isl}):
        rows = [r for r in isl if r["family"] == fam and "dG" in r]
        if rows:
            top = max(rows, key=lambda r: r["dG"])
            for r in rows:
                r["dG_max_family"] = top["dG"]
    return {"T": T, "G_complex_trans_rot_vib": g_cx, "sites": out, "islands": isl}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python -m qd_workflows.facet_kinetics")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("input", help="precursor and incorporated geometries on the fit dots")
    p.add_argument("ladder_input")
    p.add_argument("state", help="the ladder run's state.json (relaxed dots)")
    p.add_argument("--dots", required=True, help="JSON {dot id: facet family}")
    p.add_argument("--terrace", type=int, default=2)
    p.add_argument("--edge", type=int, default=1)
    p.add_argument("--on-relaxed", action="store_true", help="add the units to the relaxed dots, not the builder ones")
    p.add_argument("--islands", default="", help="island sizes in units, e.g. 1,2,3,4,6 (builder geometry)")
    p.add_argument("--out", required=True)
    p = sub.add_parser("run", help="relax, NEB, Hessians (MACE)")
    p.add_argument("input")
    p.add_argument("--out", required=True)
    p.add_argument("--model", default=os.environ.get("QDPROPS_MACE_MODEL"))
    p.add_argument("--head", default="omat_pbe")
    p.add_argument("--images", type=int, default=7)
    p = sub.add_parser("report", help="free energies and rate constants per site")
    p.add_argument("kinetics")
    p.add_argument("--T", type=float, default=500.0)
    p.add_argument("--out")
    p = sub.add_parser("clusters", help="attachment paths of the complex onto the grid's clusters (E_att(m))")
    p.add_argument("thermo")
    p.add_argument("complex", help="JSON with the relaxed complex positions [Se, Cd, Cd, Cl, Cl] (key 'positions')")
    p.add_argument("--m", default="1-12")
    p.add_argument("--dirs", type=int, default=2)
    p.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    if a.cmd == "input":
        make_input(a.ladder_input, a.state, json.loads(Path(a.dots).read_text()), a.out, a.terrace, a.edge, a.on_relaxed,
                   tuple(int(x) for x in a.islands.split(",") if x))
    elif a.cmd == "clusters":
        lo, hi = (int(x) for x in a.m.split("-"))
        cluster_input(a.thermo, json.loads(Path(a.complex).read_text())["positions"], a.out, range(lo, hi + 1), a.dirs)
    elif a.cmd == "report":
        rep = report(a.kinetics, a.T)
        if a.out:
            Path(a.out).write_text(json.dumps(rep, indent=1))
        for r in rep["sites"]:
            print(f"{r['id']:52s} dE_ads {r['dE_ads']:+.3f} dE_inc {r['dE_inc']:+.3f} barrier {r['dE_barrier']:.3f} eV"
                  + (f" | dG_ads {r['dG_ads']:+.3f} dG^ {r['dG_barrier']:.3f} eV k_f {r['k_f']:.2e} A/s" if "k_f" in r else "")
                  + f" ({r['neb']})")
        for r in rep["islands"]:
            print(f"{r['id']:52s} N {r['n_units']:2d}  dE {r['dE']:+.3f} eV  per added unit {r['dE_step']:+.3f} eV"
                  + (f"  dG {r['dG']:+.3f} eV" if "dG" in r else ""))
    else:
        from orchestr_ai.qd.run import Settings
        s = Settings(model=a.model, head=a.head, detach_quick=False, detach_screen=False, detach_rigid_refresh=False,
                     detach_frames=False, hessian="analytic")
        kin = Kinetics(a.input, a.out, s, log=lambda x: print(x, flush=True), neb_images=a.images)

        def stop(signum, _frame):
            kin.save()
            kin.log(f"[kin] signal {signum}: state saved, stopping")
            sys.exit(128 + signum)

        import signal
        signal.signal(signal.SIGTERM, stop)
        kin.run()
    return 0


if __name__ == "__main__":
    sys.exit(main())
