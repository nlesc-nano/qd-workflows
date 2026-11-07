"""
Nucleation workflow of one material: from the smallest Wulff dot with discernible facets down to the monomer,
[MA]_m(MX_q)_k for every m <= m0 and k, with free energies and the grand potential.

Stages, each resumable (a stage is done when its output says so):
  start    the smallest neutral Wulff dot of the material config whose atom-fitted polyhedron exposes every facet
           family of the recipe with at least `min_sites` surface sites each (laptop, builder only); written as a
           record folder (record.json, start.xyz, the CIF) for lattice_grid, with the two seeds of its row (seeds)
  grid     search.scheme seeded (the default): lattice_grid from the two seeds, relaxed start (optionally every
           candidate screened by a short relaxation and only the best few per node fully relaxed: screen_steps > 0);
           after each row the row anneal
           (MD + quench + polish of the nodes within anneal_window of the row's lowest Omega for some reachable
           mu_MX), whose minima start the descent to the next row.
           search.scheme legacy: lattice_grid twice from the start dot (lattice start and relaxed start, column
           descent, every candidate fully relaxed), then
  anneal   (legacy) MD anneal + quench from the best structure of every node of each start
  thermo   Hessians of the grid states and of the new annealed minima (seeded: only the nodes the row anneal chose)
  merge    one thermochemistry file with every structure's source, the Omega page, per-m trajectory files; (seeded)
           its solvation block: every state's surface and the m = 1 complexes' COSMO from `complexes`
  complexes (CPU, CP2K) ESP runs of the m = 1 complexes from this grid's structures, then COSMO on their isodensity:
           `complexes prepare` (inputs, job.sh), the job, `complexes collect` (cosmo.json, into thermo_merged.json)

Solvation: the grid's free energies are gas phase (ideal, 1 M). The growth model and the simulator add the implicit
solvation of every species (grid_solvation: COSMO + cavity, the campaign's model) from each state's geometric
surface; run `python -m qd_workflows.grid_solvation surfaces thermo_merged.json` after the merge. TO DO: explicit
solvent and ligands (oleic acid, amines, phosphines) on the monomers, on every (m, k) and on the facets.

Seeds (seeded): the start dot's faceted core at both ends of its MX_q range: fully Z-passivated (QD_builder,
post_treatment.z_type_target above any reachable total), laddered down to k = 0; and X-type only, the same dot with
every Z-type unit (a cation bonded to one anion) stripped with q X (k_X = k - their number), laddered up to k_max.
After each row the row anneal always takes the seed nodes and k_max besides the nodes the Omega window picks.

    python -m qd_workflows.nucleation start   <material.yaml> --out DIR
    python -m qd_workflows.nucleation seeds   DIR                         (laptop: the two seeds of an existing start)
    python -m qd_workflows.nucleation run     DIR --out DIR --model M      (cluster GPU: the unfinished stages in turn)
    python -m qd_workflows.nucleation prepare DIR --cluster nhr --code SRC [--pieces]

Material config: the shapes format (qd_workflows/templates/shapes/*.yaml): cif, native_elements, positive_q_mode,
preset (charges, passivation, facets with gamma, symmetry), plus a `nucleation:` block (see NUCLEATION).
"""
from __future__ import annotations

import argparse
import copy
import json
import os
import sys
import tempfile
from pathlib import Path

import numpy as np
import yaml

NUCLEATION = {   # defaults of the `nucleation:` block
    "sizes": {"start": 0.5, "stop": 3.0, "step": 0.125},      # unit cells, searched upwards
    "centres": "auto",            # construction centres: auto = every native species and the empty interstitial
                                  # sites of the CIF; or a list like the shapes config
    "min_sites": 1.0,             # an exposed facet family must hold at least this many surface sites of area
    "min_facet_atoms": 3,         # ... and one of its facet planes at least this many core atoms (a triangle)
    "min_core_atoms": 20,         # smallest core (cation + anion atoms), as the QDSpace library
    "anneal": {"ps": 4.0, "hot": [700.0, 900.0], "small_m": 6, "small_ps": 8.0, "small_hot": [700.0, 900.0, 1100.0]},
    "max_md_atoms": 6000,         # atoms of one batched MD run (one checkpoint)
    "ligand": None,               # carboxylate X instead of the Cl placeholders (qd_workflows.ligands): {"name": "oleate",
                                  # "seed": 7, "ff": "mmff", "d3": true, "md_timestep_fs": 1.0, "inner": "cbeta"}
    "m_min": 1,                   # lowest m of the grid (a larger value gives a short test run)
    "search": {
        "scheme": "seeded",       # seeded | legacy (lattice and relaxed starts from the start dot, no screening,
                                  # one anneal of every node at the end)
        "screen_steps": 0,        # optimiser steps of the screening relaxation of every candidate (0: all fully relaxed).
                                  # Off: on the CdSe m = 13 row 20 steps / best 3 lost 0.21 eV per ladder step on average
                                  # (up to 1.1 eV), 0.5-1.8 eV per node after 15 steps; even after 100 steps the true
                                  # minimum ranked up to 21st of ~30 (it comes from a rearrangement late in the relaxation)
        "screen_keep": 3,         # candidates per node fully relaxed and polished after screening
        "capacity": True,         # every +MX_q / -MA step also keeps its most Z-capable product (lattice_grid): without
                                  # it the CdSe rows m <= 9 stopped 1-3 CdCl2 below the grid with the lattice start
        "anneal_all_m": 10,       # rows m <= this anneal every node (small, cheap): with the Omega window alone 17 of the
                                  # 19 skipped nodes the merged CdSe grid marks relevant sat in these rows
        "anneal_starts": 2,       # MD starts per annealed node: its best state, then the best of each other seed path
                                  # (else the best again, with new velocities); the anneal minima differed by 0.1-0.4 eV
                                  # between two single-start runs of CdSe m = 13
        "anneal_window": 0.15,    # eV: a node is annealed when its Omega comes this close to the row's lowest ...
                                  # (G without vibrations, before the Hessians: 0.15 eV picks the nodes 0.1 eV picks
                                  # with them, 96 against 94 of the 160 CdSe nodes; 0.1 eV here missed 10 of the 94)
        "mu_T": [300.0, 700.0],   # ... for some temperature (K) in this range ...
        "mu_c": [1e-6, 1.0],      # ... and free MX_q concentration (M) in this range
    },
}


def _block(cfg: dict) -> dict:
    P = json.loads(json.dumps(NUCLEATION))
    for k, v in (cfg.get("nucleation") or {}).items():
        P[k] = P[k] | v if isinstance(v, dict) and isinstance(P.get(k), dict) else v
    return P


# ---------------- facet site areas and interstitial centres (from the CIF, any material) ----------------
def site_area(struct, normal, species: str, reach: int = 4) -> float:
    """Area per surface site of `species` on planes normal to `normal`: 1 / (n_species d), with n the species'
    number density and d the spacing of its planes along the normal (lattice centring included)."""
    lat = struct.lattice
    frac = np.array([s.frac_coords for s in struct.sites if s.specie.symbol == species])
    n_density = len(frac) / lat.volume
    rng = range(-reach, reach + 1)
    shifts = np.array([[i, j, k] for i in rng for j in rng for k in rng], float)
    pts = ((frac[None, :, :] + shifts[:, None, :]).reshape(-1, 3)) @ lat.matrix
    proj = np.unique(np.round(pts @ (np.asarray(normal, float) / np.linalg.norm(normal)), 4))
    gaps = np.diff(proj)
    d = float(gaps[gaps > 1e-3].min())
    return 1.0 / (n_density * d)


def interstitial_centres(struct, n_grid: int = 8, tol: float = 0.05) -> list[list[float]]:
    """Fractional positions of the empty high-symmetry sites of the cell, one per kind: grid points (n_grid per
    axis) with at least four equidistant nearest atoms, told apart by their first two neighbour shells (species and
    distance); e.g. zinc-blende's empty tetrahedral hole (4 cations nearest) and octahedral hole (4 anions nearest)."""
    lat = struct.lattice
    g = np.array([[i, j, k] for i in range(n_grid) for j in range(n_grid) for k in range(n_grid)], float) / n_grid
    frac = np.array([s.frac_coords for s in struct.sites])
    sym = np.array([s.specie.symbol for s in struct.sites])
    sh = np.array([[i, j, k] for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)], float)
    dist = np.linalg.norm(((g[:, None, None, :] - frac[None, None, :, :] - sh[None, :, None, :]) @ lat.matrix), axis=3)
    dist = dist.reshape(len(g), -1)                                     # every atom image
    species = np.tile(sym, len(sh))
    far = dist.min(1)
    out, seen = [], set()
    for i in range(len(g)):
        d = dist[i]
        dmin = d.min()
        if dmin < 0.9 * far.max():                                      # not a hole (on or near atoms and bonds)
            continue
        first = np.abs(d - dmin) < tol
        if first.sum() < 4 or len(set(species[first])) > 1:           # holes: >= 4 equidistant atoms of one kind
            continue
        rest = d[~first]
        d2 = rest.min()
        second = np.abs(d - d2) < tol
        key = (tuple(sorted(species[first])), round(float(dmin), 1), tuple(sorted(species[second])), round(float(d2), 1))
        if key not in seen:
            seen.add(key)
            out.append([round(float(x), 4) for x in g[i]])
    return out


# ---------------- the start: smallest Wulff dot with discernible facets ----------------
FACET_TOL = 0.6     # A: a core atom this close to a fitted facet plane lies on that facet


def exposed_sites(symbols, positions, planes, struct, cation: str) -> dict:
    """Facet families of the atom-fitted polyhedron (the builder's facet normals, each plane through the
    outermost core atom): {family: (surface sites of area, most core atoms on one of its planes)}."""
    from .shapes import polyhedron
    s, p = np.asarray(symbols), np.asarray(positions, float)
    core = np.isin(s, list({x.specie.symbol for x in struct.sites}))
    p = p - p[core].mean(0)
    fitted = [(np.asarray(n, float), float((p[core] @ n).max()), lab) for n, _d, lab in planes]
    g = polyhedron([(n, d, str(i)) for i, (n, d, _l) in enumerate(fitted)])
    area_by, out = {}, {}
    for i, (n, d, lab) in enumerate(fitted):
        area_by.setdefault(lab, site_area(struct, n, cation))
        a = g["areas_A2"].get(str(i), 0.0)
        if a <= 0:
            continue
        on = int((d - p[core] @ n < FACET_TOL).sum())
        sites, atoms = out.get(lab, (0.0, 0))
        out[lab] = (round(sites + a / area_by[lab], 2), max(atoms, on))
    return out


def find_start(config: str, out: str | None = None, log=print) -> dict:
    """Search sizes upwards (all construction centres at each size) for the first neutral dot whose atom-fitted
    polyhedron exposes every facet family of the recipe with >= min_sites sites; write it as a record folder."""
    from pymatgen.core import Structure
    from builder.library_record import describe_structure
    from builder.scripts.generate_library import _build
    from .shapes import construction_planes, family_label
    from builder.config import _parse_hkl
    cfg_path = Path(config).resolve()
    cfg = yaml.safe_load(cfg_path.read_text())
    P = _block(cfg)
    cif = str((cfg_path.parent / cfg["cif"]).resolve())
    struct = Structure.from_file(cif)
    native = list(cfg["native_elements"])
    cation = native[0]
    preset = cfg["preset"]
    want = sorted({family_label(*_parse_hkl(str(f["hkl"])), f.get("termination")) for f in preset["facets"]})
    if P["centres"] == "auto":
        centres = [{"center_on_species": e} for e in native] + [{"fractional_shift": c} for c in interstitial_centres(struct)]
    else:
        centres = [c if isinstance(c, dict) else {"center_on_species": c} for c in P["centres"]]
    sz = P["sizes"]
    sizes = [round(sz["start"] + i * sz["step"], 4) for i in range(int(round((sz["stop"] - sz["start"]) / sz["step"])) + 1)]
    tried = []
    with tempfile.TemporaryDirectory() as tmp:
        for size in sizes:
            found = []
            for c in centres:
                recipe = copy.deepcopy(preset)
                recipe["size_unit_cells"] = [size] * 3
                recipe["construction_origin"] = c
                res = _build(cif, recipe, cfg.get("positive_q_mode", "add"), Path(tmp) / f"b{len(tried)}")
                row = {"size": size, "centre": c}
                if "error" in res or not len(res.get("pts", [])):
                    tried.append(row | {"status": "build failed"})
                    continue
                desc = describe_structure(res["symbols"], res["pts"], native_order=native, charges=recipe["charges"],
                                          centre_index=res["centre_index"])
                row.update(formula=desc["formula"], n_atoms=desc["n_atoms"], charge=desc["total_charge"])
                if desc["total_charge"] != 0:
                    tried.append(row | {"status": "charged"})
                    continue
                if desc["size"]["n_core"] < P["min_core_atoms"]:
                    tried.append(row | {"status": f"core {desc['size']['n_core']} atoms"})
                    continue
                sites = exposed_sites(res["symbols"], res["pts"], construction_planes(cif, recipe), struct, cation)
                row["sites"] = sites
                ok = all(f in sites and sites[f][0] >= P["min_sites"] and sites[f][1] >= P["min_facet_atoms"] for f in want)
                tried.append(row | {"status": "start" if ok else "facets missing"})
                if ok:
                    found.append((list(res["symbols"]).count(native[1]), row, res, recipe, desc))
            if found:
                m, row, res, recipe, desc = min(found, key=lambda x: (x[0], x[3]["construction_origin"].__repr__()))
                log(f"[nucleation] start: {desc['formula']} (m = {m}) at {size} unit cells, centre {row['centre']}; "
                    f"sites per family {row['sites']}")
                break
        else:
            raise SystemExit(f"no dot up to {sizes[-1]} unit cells exposes {want} with >= {P['min_sites']} sites each")
    result = {"formula": desc["formula"], "m": m, "size_unit_cells": row["size"], "centre": row["centre"],
              "sites": row["sites"], "families": want, "tried": tried}
    if out:
        d = Path(out) / "start"
        d.mkdir(parents=True, exist_ok=True)
        (d / "start.xyz").write_text(f"{len(res['symbols'])}\n{desc['formula']} nucleation start\n" + "\n".join(
            f"{s} {x:.6f} {y:.6f} {z:.6f}" for s, (x, y, z) in zip(res["symbols"], np.asarray(res["pts"]))) + "\n")
        (d / Path(cif).name).write_text(Path(cif).read_text())
        record = {"id": f"{cfg['material']}-{desc['formula']}-nucleation-start", "material": cfg["material"],
                  "phase": cfg.get("phase"), "formula": desc["formula"], "n_atoms": desc["n_atoms"],
                  "origin": {"generator": "qd_workflows.nucleation start", "cif": Path(cif).name, "recipe": recipe,
                             "positive_q_mode": cfg.get("positive_q_mode", "add"), "config": str(cfg_path)}}
        (d / "record.json").write_text(json.dumps(record, indent=1))
        (Path(out) / "start.json").write_text(json.dumps(result, indent=1, default=float))
    return result


# ---------------- the seeds of the top row: both ends of the MX_q range (builder z_type_target) ----------------
SEED_X, SEED_FULL = "seed_x.xyz", "seed_full.xyz"
BOND_TOL = 1.1      # x the bulk nearest-neighbour distance (as lattice_grid)


def single_bond_cations(symbols, positions, cation: str, anion: str, d_nn: float) -> int:
    """Cations bonded to exactly one anion: Z-type MX_q units on a host anion (QD_builder's definition)."""
    from scipy.spatial import cKDTree
    p = np.asarray(positions, float)
    an = [i for i, s in enumerate(symbols) if s == anion]
    tree = cKDTree(p[an])
    return sum(1 for i, s in enumerate(symbols) if s == cation and len(tree.query_ball_point(p[i], BOND_TOL * d_nn)) == 1)


def strip_z_units(symbols, positions, cation: str, anion: str, ligand: str, q: int, d_nn: float):
    """The X-type-only dot of a Z-passivated one: every cation bonded to a single anion (a Z-type MX_q unit) removed,
    then q X per removed cation: first every X left without a cation, then, one at a time, the most crowded X (the
    most X on its cations; then the most bonded, then the outermost), so the X left spread evenly over the kept
    cations. (The units need not own q X each: on Cd28Se13Cl30 the 12 Z-type Cd hold 20 Cl; the 4 others are the
    Cl bridging three edge Cd, and the 6 Cl left bridge the 12 edge Cd pairwise.) None when more X are left without
    a cation than the units take, or the rest falls apart. Returns (symbols, positions, number of units removed)."""
    from scipy.spatial import cKDTree
    s, p = list(symbols), np.asarray(positions, float)
    nb = [set() for _ in s]
    for i, j in cKDTree(p).query_pairs(BOND_TOL * d_nn):
        if (s[i] == cation) != (s[j] == cation) and {s[i], s[j]} & {anion, ligand}:
            nb[i].add(j)
            nb[j].add(i)
    R = {i for i, e in enumerate(s) if e == cation and sum(s[j] == anion for j in nb[i]) == 1}
    if not R:
        return s, p, 0
    xs = {x for x in range(len(s)) if s[x] == ligand}
    gone = {x for x in xs if nb[x] <= R}
    if len(gone) > q * len(R):
        return None
    left = xs - gone
    n_x = lambda c: len(nb[c] & left)
    centre = p.mean(axis=0)
    while len(gone) < q * len(R):
        x = max(sorted(left), key=lambda x: (sum(n_x(c) for c in nb[x] - R), len(nb[x] - R),
                                              float(np.linalg.norm(p[x] - centre))))
        gone.add(x)
        left.discard(x)
    removed = R | gone
    keep = [i for i in range(len(s)) if i not in removed]
    alive, seen, todo = set(keep), set(), [keep[0]]
    while todo:
        a = todo.pop()
        if a not in seen:
            seen.add(a)
            todo += [b for b in nb[a] if b in alive and b not in seen]
    if len(seen) != len(keep):
        return None
    return [s[i] for i in keep], p[keep], len(R)


def _z_target(recipe: dict, total: int) -> dict:
    r = copy.deepcopy(recipe)
    pt = {k: v for k, v in (r.get("post_treatment") or {}).items() if k not in ("z_type_ligands", "z_type_displacement")}
    r["post_treatment"] = pt | {"z_type_target": {"target_total": int(total)}}
    return r


def build_seeds(rec_dir, log=print) -> dict:
    """The two seeds of the start record's row: fully Z-passivated, built by QD_builder with post_treatment.z_type_target
    on the start's recipe (a target no dot reaches: every free host site taken; laddered down), and X-type only, the
    same dot with its Z-type units stripped (strip_z_units: k_X = k_full - the single-bond cations; laddered up). When
    stripping fails, the builder's own displacement to k_X (z_type_target, D'Hondt over the families) stands in: on
    CdSe m = 13 that one kept a single-bond Cd and removed edge Cd, 3.4 eV above the stripped dot. Writes
    seed_full.xyz, seed_x.xyz and seeds.json into the record folder."""
    from pymatgen.core import Structure
    from builder.scripts.generate_library import _build
    from .grid import _formula
    rec_dir = Path(rec_dir)
    rec = json.loads((rec_dir / "record.json").read_text())
    o = rec["origin"]
    recipe, mode, cif = o["recipe"], o.get("positive_q_mode", "add"), str(rec_dir / o["cif"])
    struct = Structure.from_file(cif)
    charges = recipe["charges"]
    native = sorted({s.specie.symbol for s in struct.sites}, key=lambda e: -charges[e])
    cation, anion, lig = native[0], native[1], recipe["passivation"]["ligand"]
    q = charges[cation] // -charges[lig]
    d_nn = float(min(n.nn_distance for n in struct.get_neighbors(struct.sites[0], r=4.0)))
    mk = lambda s: (s.count(anion), s.count(lig) // q)
    s0, _p0 = _read_xyz(rec_dir / "start.xyz")
    out = []
    with tempfile.TemporaryDirectory() as tmp:
        def build(total, tag):
            res = _build(cif, _z_target(recipe, total), mode, Path(tmp) / tag)
            if "error" in res:
                raise RuntimeError(f"seed {tag}: {res['error']}")
            return (list(res["symbols"]), np.asarray(res["pts"], float),
                    (res.get("manifest") or {}).get("z_type_target_ledger"))
        sf, pf, lf = build(mk(s0)[1] + s0.count(anion), "full")
        m, k_full = mk(sf)
        k_x = k_full - single_bond_cations(sf, pf, cation, anion, d_nn)
        strip = strip_z_units(sf, pf, cation, anion, lig, q, d_nn)
        if strip is not None:
            sx, px, n_z = strip
            lx = {"method": "strip_z_units", "units_removed": n_z, "k_full": k_full}
        else:
            sx, px, lx = build(k_x, "x")
            lx = (lx or {}) | {"method": "z_type_target"}
        if mk(sx) != (m, k_x):
            raise RuntimeError(f"X-type seed: m, k = {mk(sx)}, not {(m, k_x)}")
        for f, s, p, led, way in ((SEED_FULL, sf, pf, lf, "down"), (SEED_X, sx, px, lx, "up")):
            _write_xyz(rec_dir / f, s, p, f"m={m} k={mk(s)[1]} nucleation seed")
            out.append({"file": f, "m": m, "k": mk(s)[1], "ladder": way, "formula": _formula(s),
                        "ledger": {x: v for x, v in (led or {}).items() if x != "order"}})
    if k_x == k_full:                                       # no Z-type unit: one seed, laddered both ways
        out = [out[0] | {"ladder": "both"}]
    res = {"seeds": out}
    (rec_dir / "seeds.json").write_text(json.dumps(res, indent=1, default=float))
    log("[nucleation] seeds: " + ", ".join(f"{x['formula']} (k = {x['k']}, ladder {x['ladder']})" for x in out))
    return res


# ---------------- cluster stages (MACE; orchestr_ai) ----------------
LINK = 3.3          # A: atoms closer than this belong to one piece (annealed quenches that split are dropped)
STARTS = ("lattice", "relaxed")
WINDOW, TOL_E = 0.3, 1e-3      # eV: annealed minima within WINDOW of their node's best get a Hessian; energies within
                               # TOL_E are one minimum


def _read_xyz(path):
    L = Path(path).read_text().splitlines()
    n = int(L[0])
    return [l.split()[0] for l in L[2:2 + n]], np.array([[float(x) for x in l.split()[1:4]] for l in L[2:2 + n]])


def _write_xyz(path, symbols, pts, comment):
    Path(path).write_text(f"{len(symbols)}\n{comment}\n" + "\n".join(
        f"{s:<2s} {x:12.6f} {y:12.6f} {z:12.6f}" for s, (x, y, z) in zip(symbols, np.asarray(pts))) + "\n")


def _connected(pos):
    from scipy.spatial import cKDTree
    nb = [[] for _ in range(len(pos))]
    for i, j in cKDTree(pos).query_pairs(LINK, output_type="ndarray"):
        nb[i].append(j)
        nb[j].append(i)
    seen, todo = {0}, [0]
    while todo:
        for j in nb[todo.pop()]:
            if j not in seen:
                seen.add(j)
                todo.append(j)
    return len(seen) == len(pos)


def _done(path) -> bool:
    return Path(path).is_file() and json.loads(Path(path).read_text()).get("complete", False)


def _settings(model, head, **md):
    from orchestr_ai.qd.run import Settings
    return Settings(model=model, head=head, detach_quick=False, detach_screen=False, detach_rigid_refresh=False,
                    detach_frames=False, hessian="analytic", **md)


def stage_grid(wd: Path, settings, P: dict, log=print, model=None, head="omat_pbe", refs=None):
    """seeded: lattice_grid from the seeds (relaxed start, optional screening, row anneal) into grid/. legacy: lattice_grid
    from the start dot, once per start (lattice positions; the parent's relaxed coordinates)."""
    from .lattice_grid import LatticeGrid
    rec = wd / "start"
    cif = rec / json.loads((rec / "record.json").read_text())["origin"]["cif"]
    S = P["search"]
    if S["scheme"] == "seeded":
        if _done(wd / "grid" / "grid.json"):
            return
        if not (rec / "seeds.json").is_file():
            raise SystemExit(f"{rec}/seeds.json missing: run `python -m qd_workflows.nucleation seeds {wd}` (builder)")
        seeds = [(x["file"], x["ladder"]) for x in json.loads((rec / "seeds.json").read_text())["seeds"]]
        log(f"[nucleation] grid, seeded: seeds {seeds}, {S['screen_steps']} screening steps, best {S['screen_keep']}")
        LatticeGrid(rec, wd / "grid", settings, str(cif), m_min=int(P["m_min"]), log=log, start="relaxed",
                    screen=int(S["screen_steps"]), screen_keep=int(S["screen_keep"]),
                    capacity=bool(S.get("capacity", True))).run_seeded(
            seeds, row_hook=row_anneal(wd, P, model, head, refs, log))
        return
    if S["scheme"] != "legacy":
        raise ValueError(f"search.scheme: seeded or legacy, not {S['scheme']!r}")
    for start in STARTS:
        out = wd / f"grid_{start}"
        if _done(out / "grid.json"):
            continue
        log(f"[nucleation] grid, {start} start")
        LatticeGrid(rec, out, settings, str(cif), m_min=int(P["m_min"]), log=log, start=start).run(descent="column")


def anneal_starts(wd: Path) -> list[dict]:
    """The best structure of every node of each grid, in a fixed order (m, k, start)."""
    out = []
    for start in STARTS:
        g = json.loads((wd / f"grid_{start}" / "grid.json").read_text())
        for n in g["nodes"]:
            s, p = _read_xyz(wd / f"grid_{start}" / "nodes" / f"m{n['m']:02d}_k{n['k']:02d}.xyz")
            out.append({"m": n["m"], "k": n["k"], "source": f"{start}_start", "E_eV": n["E_eV"], "symbols": s,
                        "positions": p.round(6).tolist()})
    return sorted(out, key=lambda x: (x["m"], x["k"], x["source"]))


def anneal_chunks(starts: list[dict], P: dict) -> list[dict]:
    """Starts grouped by schedule (m <= small_m: longer, hotter, more replicas) and cut, in order, into MD batches of
    at most max_md_atoms atoms (all replicas counted): each chunk is one checkpoint."""
    A = P["anneal"]
    chunks = []
    for small in (True, False):
        hot = A["small_hot"] if small else A["hot"]
        ps = A["small_ps"] if small else A["ps"]
        cur, atoms = [], 0
        for i, s in enumerate(starts):
            if (s["m"] <= A["small_m"]) != small:
                continue
            a = len(s["symbols"]) * len(hot)
            if cur and atoms + a > P["max_md_atoms"]:
                chunks.append({"starts": cur, "hot": hot, "ps": ps})
                cur, atoms = [], 0
            cur.append(i)
            atoms += a
        if cur:
            chunks.append({"starts": cur, "hot": hot, "ps": ps})
    return chunks


def anneal_chunk(starts, chunk, out: Path, settings, seed=0, cold=300.0, skip_ps=0.5, stride_ps=0.5, log=print):
    """MD replicas (linear ramps hot -> cold, CSVR, one batched MACE call per step) from each start; frames every
    stride_ps quenched, those within POLISH_WINDOW of their node's best polished, split quenches dropped.
    Writes anneal.json (per start: its best quench) and anneal_best.xyz; the MD frames are kept, so a rerun goes
    straight to the cached quenches."""
    import pickle
    from types import SimpleNamespace
    from orchestr_ai.qd.steps import detachment as dt
    from orchestr_ai.qd.steps import ensemble as ens
    out.mkdir(parents=True, exist_ok=True)
    if _done(out / "anneal.json"):
        return
    sel = [starts[i] for i in chunk["starts"]]
    md_file = out / "md_frames.pkl"
    h_mass = (getattr(settings, "ligand", None) or {}).get("h_mass")
    masses0 = ens._masses
    if h_mass:                          # heavy H (amu) for a longer MD step; quench minima do not depend on masses
        ens._masses = lambda sy: np.where(np.asarray(sy) == "H", float(h_mass), masses0(sy))
    if md_file.is_file():
        frames, reps = pickle.loads(md_file.read_bytes())
    else:
        rng = np.random.default_rng(seed)
        systems, vels, scheds, names, reps = [], [], [], [], []
        for i, st in enumerate(sel):
            masses = ens._masses(st["symbols"])
            for hot in chunk["hot"]:
                systems.append((st["symbols"], np.asarray(st["positions"], float)))
                vels.append(ens.maxwell_velocities(masses, hot, rng))
                scheds.append([hot, cold])
                names.append(f"{st['source']}_m{st['m']}_k{st['k']}_{int(hot)}K")
                reps.append({"start": i, "hot": hot})
        n_steps = int(round(chunk["ps"] * 1000 / settings.md_timestep_fs))
        fr, summ = ens.md_replicas(ens._evaluator(settings), systems, vels, scheds, names, settings, n_steps,
                                   int(round(skip_ps * 1000 / settings.md_timestep_fs)),
                                   int(round(stride_ps * 1000 / settings.md_timestep_fs)), rng)
        frames = [[(f["positions"], f["info"]["time_fs"]) for f in fl] for fl in fr]
        md_file.write_bytes(pickle.dumps((frames, reps)))
        log(f"[nucleation] anneal {out.name}: MD {len(systems)} replicas, {sum(x['stopped'] is not None for x in summ)} stopped")
    ens._masses = masses0
    cache = dt.RelaxCache(SimpleNamespace(props=out, settings=settings, record_dir=out))
    jobs = [(r, j) for r, fl in enumerate(frames) for j in range(len(fl))]
    key = lambda r: (sel[reps[r]["start"]]["m"], sel[reps[r]["start"]]["k"])
    lig = getattr(settings, "ligand", None) or {}
    screen = int(lig.get("quench_screen_steps", 0))
    if screen:
        # two-stage quench (long-chain X: most of a quench is the chains settling): every frame relaxed screen steps,
        # then per node the quench_keep lowest (and each start's lowest) relaxed to the end from there; the others
        # dropped. Oleate m = 13 traces: the final best was among the 3 lowest after 10-100 steps, top 8 at 50 steps
        # kept it at k = 3 and 7 and saved ~75% of the steps.
        import pickle as _pk
        sf = out / "quench_screen.pkl"
        if sf.is_file():
            scr = _pk.loads(sf.read_bytes())
        else:
            from orchestr_ai.qd import batch_relax as br
            from orchestr_ai.qd.engines import mace_calculator
            ev = br.BatchMACE(mace_calculator(settings.head, settings.model, settings.device, settings.dtype),
                              settings.relax_batch_atoms)
            rr = br.relax_many(cache._perturbed([(sel[reps[r]["start"]]["symbols"], frames[r][j][0]) for r, j in jobs]),
                               ev, max_steps=screen)
            scr = [(np.asarray(x.positions), float(x.energy)) for x in rr]
            sf.write_bytes(_pk.dumps(scr))
        groups = {}
        for i, (r, _j) in enumerate(jobs):
            groups.setdefault(key(r), []).append(i)
        kept = set()
        for idx in groups.values():
            idx = sorted(idx, key=lambda i: scr[i][1])
            kept |= set(idx[:int(lig.get("quench_keep", 8))])
            firsts = {}
            for i in idx:
                firsts.setdefault(reps[jobs[i][0]]["start"], i)
            kept |= set(firsts.values())
        kept = sorted(kept)
        full = cache.relax_many([(sel[reps[jobs[i][0]]["start"]]["symbols"], scr[i][0]) for i in kept], settings)
        res, ok = [None] * len(jobs), [False] * len(jobs)
        for i, x in zip(kept, full):
            res[i], ok[i] = x, _connected(x[0].get_positions())
        log(f"[nucleation] anneal {out.name}: {len(jobs)} quenches screened {screen} steps, {len(kept)} relaxed to the end")
    else:
        res = list(cache.relax_many([(sel[reps[r]["start"]]["symbols"], frames[r][j][0]) for r, j in jobs], settings))
        ok = [_connected(x[0].get_positions()) for x in res]
    best = {}
    for (r, _j), x, good in zip(jobs, res, ok):
        if good:
            best[key(r)] = min(best.get(key(r), x[1]), x[1])
    window = [i for i, ((r, _j), x) in enumerate(zip(jobs, res)) if ok[i] and x[1] <= best[key(r)] + dt.POLISH_WINDOW]
    pol = cache.relax_many([(res[i][0].get_chemical_symbols(), res[i][0].get_positions()) for i in window], settings, polish=True)
    for i, x in zip(window, pol):
        res[i] = x
        ok[i] = _connected(x[0].get_positions())
    cache._save()
    polished, per = set(window), {}
    for i, ((r, j), x) in enumerate(zip(jobs, res)):
        si = reps[r]["start"]
        if ok[i] and i in polished and (si not in per or x[1] < per[si][0]):
            per[si] = (x[1], x[0].get_positions(), reps[r]["hot"], frames[r][j][1])
    rows, text = [], []
    for si, st in enumerate(sel):
        e, p, hot, tfs = per.get(si, (None, None, None, None))
        rows.append({"m": st["m"], "k": st["k"], "source": st["source"], "E_start_eV": st["E_eV"], "E_anneal_eV": e,
                     "gain_eV": None if e is None else e - st["E_eV"], "hot_K": hot, "frame_fs": tfs})
        if p is not None:
            text.append(f"{len(p)}\nm={st['m']} k={st['k']} anneal from {st['source']} ({int(hot)} K ramp, frame {tfs:.0f} fs) "
                        f"E={e:.6f} start E={st['E_eV']:.6f}\n" +
                        "\n".join(f"{a} {x:.6f} {y:.6f} {z:.6f}" for a, (x, y, z) in zip(st["symbols"], p)))
    (out / "anneal_best.xyz").write_text("\n".join(text) + "\n")
    (out / "anneal.json").write_text(json.dumps({"complete": True, "chunk": chunk, "rows": rows}, indent=1))
    log(f"[nucleation] anneal {out.name}: {len(jobs)} quenches, {sum(1 for r in rows if r['gain_eV'] is not None and r['gain_eV'] < -0.05)} "
        f"of {len(rows)} starts lower by > 0.05 eV")


def stage_anneal(wd: Path, P: dict, model: str, head: str, log=print):
    s = md_settings(model, head)
    starts = anneal_starts(wd)
    chunks = anneal_chunks(starts, P)
    (wd / "anneal").mkdir(exist_ok=True)
    (wd / "anneal" / "chunks.json").write_text(json.dumps({"n_starts": len(starts), "chunks": chunks}))
    for c, ch in enumerate(chunks):
        anneal_chunk(starts, ch, wd / "anneal" / f"chunk{c:03d}", s, seed=c, log=log)


def chunk_dirs(wd: Path) -> list:
    """The anneal chunk folders (chunk000, chunk001, ...), not the chunks.json index beside them."""
    return [d for d in sorted((wd / "anneal").glob("chunk[0-9]*")) if d.is_dir()]


def stage_ensemble(wd: Path, log=print):
    """The annealed minima as a grid folder for grid_thermo: per node the best quench of each start that is lower
    than its start (a new minimum) and within WINDOW of the node's best; source anneal_from_<start source>."""
    from orchestr_ai.qd.records import fingerprint
    from .grid import _formula
    out = wd / "anneal" / "grid"
    if _done(out / "grid.json"):
        return
    (out / "nodes").mkdir(parents=True, exist_ok=True)
    found, start_best = {}, {}
    for d in chunk_dirs(wd):
        rows = json.loads((d / "anneal.json").read_text())["rows"]
        L = (d / "anneal_best.xyz").read_text().splitlines()
        frames, i = [], 0
        while i < len(L) and L[i].strip():
            n = int(L[i])
            frames.append(([l.split()[0] for l in L[i + 2:i + 2 + n]], [[float(x) for x in l.split()[1:4]] for l in L[i + 2:i + 2 + n]]))
            i += 2 + n
        fi = 0
        for r in rows:
            key = (r["m"], r["k"])
            start_best[key] = min(start_best.get(key, r["E_start_eV"]), r["E_start_eV"])
            if r["E_anneal_eV"] is None:
                continue
            s, p = frames[fi]
            fi += 1
            if r["gain_eV"] < -TOL_E:
                found.setdefault(key, []).append((r["E_anneal_eV"], f"anneal_from_{r['source']}", s, p))
    nodes = []
    for key in sorted(found, key=lambda x: (-x[0], x[1])):
        items = sorted(found[key], key=lambda x: x[0])
        top = min(items[0][0], start_best[key])
        states, ens_e = [], []
        for e, src, s, p in items:
            dup = [x for x in states if abs(x["E_eV"] - e) < TOL_E]
            if dup:
                if src not in dup[0]["source"]:
                    dup[0]["source"] += f" + {src}"
                continue
            if any(abs(e - x) < TOL_E for x in ens_e):
                continue
            ens_e.append(e)
            if e <= top + WINDOW:
                tag = f"m{key[0]:02d}_k{key[1]:02d}{'' if not states else f'_{len(states)}'}"
                _write_xyz(out / "nodes" / f"{tag}.xyz", s, p, f"m={key[0]} k={key[1]} {src} E={e:.6f} eV")
                states.append({"E_eV": e, "fingerprint": fingerprint(s, p), "source": src})
        if states:
            nodes.append({"m": key[0], "k": key[1], "formula": _formula(items[0][2]), "E_eV": states[0]["E_eV"],
                          "states": states, "ensemble_E_eV": ens_e})
    (out / "grid.json").write_text(json.dumps({"record": "annealed minima", "complete": True, "nodes": nodes}, indent=1))
    log(f"[nucleation] ensemble: {len(nodes)} nodes, {sum(len(n['states']) for n in nodes)} new minima")


# ---------------- seeded scheme: row anneal of the nodes the pool can reach ----------------
KB_EV = 8.617333262e-5


def omega_gaps(G: dict, g_mx, temps, t_range, c_range, n_c: int = 61) -> dict:
    """{k: smallest Omega(k) - min_k' Omega(k') over the conditions}: Omega(k) = G(k) - k mu_MX, mu_MX = G_MX(T) +
    kT ln(c / 1 M) for T in t_range (the temperatures of `temps` inside it) and c in c_range (log grid).
    G: {k: G(k; T) on temps}; g_mx: G_MX(T) on temps."""
    temps = np.asarray(temps, float)
    ts = temps[(temps >= t_range[0] - 1e-9) & (temps <= t_range[1] + 1e-9)]
    ks = sorted(G)
    gap = {k: np.inf for k in ks}
    lnc = np.log(np.logspace(np.log10(c_range[0]), np.log10(c_range[1]), n_c))
    for t in ts:
        gk = np.array([np.interp(t, temps, G[k]) for k in ks])
        mu = np.interp(t, temps, g_mx) + KB_EV * t * lnc                  # (n_c,)
        om = gk[:, None] - np.asarray(ks, float)[:, None] * mu[None, :]    # (k, n_c)
        d = (om - om.min(axis=0)).min(axis=1)
        for k, x in zip(ks, d):
            gap[k] = min(gap[k], float(x))
    return gap


def anneal_parents(states, n: int) -> list:
    """Up to n MD starts of a node, [(state, label)]: its best state, then the best state of each other seed path,
    then its other states, distinct structures only (fingerprint and energy); the best again (new velocities) for the
    starts still missing."""
    states = sorted(states, key=lambda x: x.energy)
    skip = (None, "anneal", "capacity", states[0].path)
    order = [states[0]] + [x for x in states[1:] if x.path not in skip] + [x for x in states[1:] if x.path != "capacity"]
    out = []
    for x in order:
        if len(out) == n:
            break
        if all(x.fp != y.fp and abs(x.energy - y.energy) > TOL_E for y, _l in out):
            out.append((x, "best" if not out else f"path:{x.path}" if x.path not in (None, "anneal", "capacity")
                        else "state"))
    while len(out) < n:
        out.append((states[0], f"best, repeat {len(out)}"))
    return out


def _anneal_results(d: Path) -> list:
    """[(row, (symbols, positions) or None)] of an anneal chunk folder, in start order."""
    rows = json.loads((d / "anneal.json").read_text())["rows"]
    L = (d / "anneal_best.xyz").read_text().splitlines()
    frames, i = [], 0
    while i < len(L) and L[i].strip():
        n = int(L[i])
        frames.append(([l.split()[0] for l in L[i + 2:i + 2 + n]],
                       np.array([[float(x) for x in l.split()[1:4]] for l in L[i + 2:i + 2 + n]])))
        i += 2 + n
    out, fi = [], 0
    for r in rows:
        if r["E_anneal_eV"] is None:
            out.append((r, None))
        else:
            out.append((r, frames[fi]))
            fi += 1
    return out


def md_settings(model, head):
    # QDW_MD_BATCH_ATOMS: atoms per batched call of the MD and its quenches (6000 on an 80 GB A100; 4000 on a 64 GB
    # MI250X GCD in float64)
    return _settings(model, head, md_thermostat="csvr", md_timestep_fs=2.0, md_tau_fs=100.0, md_max_force=50.0,
                     relax_batch_atoms=int(os.environ.get("QDW_MD_BATCH_ATOMS", 6000)))


def row_anneal(wd: Path, P: dict, model: str, head: str, refs, log=print):
    """The lattice_grid row hook: Omega gaps of the row's nodes (best states, G without vibrations for clusters and
    the MX_q monomer alike, over the search's mu_T x mu_c range), MD anneal + quench + polish of the nodes within
    anneal_window, and always of the seed nodes and k_max, whose ladder energies can be far off (k = 15 of CdSe m = 13:
    1.6 eV above its anneal minimum) (one batch: anneal_chunks), and every quench lower than its node's best added to the node (path
    'anneal'), so the descent to m - 1 starts from it. anneal_rows/m<m>/selection.json: gaps, chosen nodes, gains,
    wall time and force evaluations."""
    from .grid_thermo import TEMPERATURES, _gas, _monomers, unit_names
    from .lattice_grid import EVALS
    S = P["search"]

    def hook(g, m):
        import time
        t0, ev0 = time.time(), dict(EVALS)
        out = wd / "anneal_rows" / f"m{m:02d}"
        out.mkdir(parents=True, exist_ok=True)
        best = {k: g.nodes[(mm, k)]["states"][0] for (mm, k) in sorted(g.nodes) if mm == m and g.nodes[(mm, k)]["states"]}
        mx = unit_names([{"formula": _formula_of(st.syms)} for st in best.values() if st.syms.count(g.lat.ligand)])[0]
        mono = _monomers(refs, head, [mx], TEMPERATURES, vibrations=False) if refs else {}
        if mx in mono:
            G = {k: _gas(g.rsyms(st.syms), st.relaxed, st.energy, [], TEMPERATURES)[1] for k, st in best.items()}
            gaps = omega_gaps(G, mono[mx]["G_eV"], TEMPERATURES, S["mu_T"], S["mu_c"])
        else:
            log(f"[nucleation] row anneal m={m}: no {mx} monomer reference: every node annealed")
            gaps = {k: 0.0 for k in best}
        always = {k for (mm, k) in getattr(g, "seed_nodes", ()) if mm == m and k in best} | {max(best)}
        every = m <= int(S.get("anneal_all_m", 0))
        chosen = [k for k in best if every or gaps[k] <= S["anneal_window"] or k in always]
        parents, starts = [], []
        for k in chosen:
            for st, label in anneal_parents(g.nodes[(m, k)]["states"], int(S.get("anneal_starts", 1))):
                parents.append(st)
                starts.append({"m": m, "k": k, "source": f"row{len(starts)}", "label": label, "E_eV": st.energy,
                               "symbols": g.rsyms(st.syms), "positions": np.asarray(st.relaxed).round(6).tolist()})
        chunks = anneal_chunks(starts, P)
        s_md = md_settings(model, head)
        if P.get("ligand"):                      # H atoms: shorter MD step (D3, if asked, is in every BatchMACE)
            s_md.md_timestep_fs = float(P["ligand"].get("md_timestep_fs", 1.0))
            s_md.ligand = dict(P["ligand"])
        md_evals = 0
        for c, ch in enumerate(chunks):
            anneal_chunk(starts, ch, out / f"chunk{c:03d}", s_md, seed=1000 * m + c, log=log)
            n_steps = int(round(ch["ps"] * 1000 / s_md.md_timestep_fs))
            md_evals += sum(len(starts[i]["symbols"]) for i in ch["starts"]) * len(ch["hot"]) * n_steps
        gained = {}
        for c, ch in enumerate(chunks):
            for (r, fr), i in zip(_anneal_results(out / f"chunk{c:03d}"), ch["starts"]):
                k = starts[i]["k"]
                if fr is not None and r["E_anneal_eV"] < best[k].energy - TOL_E:
                    g.add_state(m, k, parents[i], fr[1], r["E_anneal_eV"],
                                {"move": "anneal", "source": "anneal", "hot_K": r["hot_K"], "frame_fs": r["frame_fs"],
                                 "anneal_start": starts[i]["label"]})
                    gained[k] = max(gained.get(k, 0.0), round(best[k].energy - r["E_anneal_eV"], 4))
        old = json.loads((out / "selection.json").read_text()) if (out / "selection.json").is_file() else {}
        sel = {"m": m, "window_eV": S["anneal_window"], "mu_T": S["mu_T"], "mu_c": S["mu_c"],
               "gaps_eV": {str(k): round(v, 4) for k, v in gaps.items()}, "annealed": chosen, "always": sorted(always),
               "starts": [{"k": x["k"], "label": x["label"], "E_eV": round(x["E_eV"], 6)} for x in starts],
               "gained_eV": {str(k): v for k, v in gained.items()},
               "seconds": old.get("seconds", round(time.time() - t0, 1)),
               "md_atom_force_evals": md_evals, "quench_force_evals": old.get("quench_force_evals", EVALS["structures"] - ev0["structures"]),
               "quench_atom_force_evals": old.get("quench_atom_force_evals", EVALS["atoms"] - ev0["atoms"])}
        (out / "selection.json").write_text(json.dumps(sel, indent=1))
        log(f"[nucleation] row anneal m={m}: {len(chosen)} of {len(best)} nodes annealed (k = {chosen}, "
            f"{len(starts)} starts); "
            f"lower at {gained}")
        g.save()
    return hook


def _formula_of(symbols) -> str:
    from .grid import _formula
    return _formula(symbols)


def annealed_nodes(wd: Path) -> list:
    """(m, k) of every node the row anneals chose (the Hessians of the seeded scheme)."""
    return [(json.loads(f.read_text())["m"], k) for f in sorted((wd / "anneal_rows").glob("m*/selection.json"))
            for k in json.loads(f.read_text())["annealed"]]


def stage_thermo(wd: Path, settings, refs, log=print, scheme="legacy"):
    from .grid_thermo import hessians
    if scheme == "seeded":
        if not _done(wd / "grid" / "thermo.json"):
            hessians(wd / "grid", settings, refs, log=log, only=annealed_nodes(wd))
        return
    for d in [wd / f"grid_{s}" for s in STARTS] + [wd / "anneal" / "grid"]:
        if not _done(d / "thermo.json"):
            hessians(d, settings, refs, log=log)


def merge_thermo(out: Path, runs: list, tol=TOL_E) -> dict:
    """Union of thermo.json files, node by node; one copy per minimum (energy within tol: grid fingerprints are
    those of lattice configurations, which the two starts share); every state keeps its source (its own, or the
    run's label; a minimum found by several runs lists them all). Structure files are written relative to the folder
    of `out`, so the merged thermo opens its structures wherever the runs sit."""
    import os
    merged, base = {}, None
    for label, f in runs:
        th = json.loads(Path(f).read_text())
        rel = os.path.relpath(Path(f).resolve().parent, Path(out).resolve().parent)
        base = base or th
        for n in th["nodes"]:
            node = merged.setdefault((n["m"], n["k"]), {"m": n["m"], "k": n["k"], "formula": n["formula"], "states": [],
                                                        "ensemble_E_eV": []})
            for s in n["states"]:
                src = s.get("source", label)
                same = [x for x in node["states"] if abs(s["E_eV"] - x["E_eV"]) < tol]
                if same:
                    if src not in same[0]["source"].split(" + "):
                        same[0]["source"] += f" + {src}"
                else:
                    node["states"].append(s | {"file": os.path.normpath(os.path.join(rel, s["file"])), "source": src})
            for e in n["ensemble_E_eV"]:
                if not any(abs(e - x) < tol for x in node["ensemble_E_eV"]):
                    node["ensemble_E_eV"].append(e)
    nodes = sorted(merged.values(), key=lambda n: (-n["m"], n["k"]))
    for n in nodes:
        n["states"].sort(key=lambda s: s["E_eV"])
        n["ensemble_E_eV"].sort()
    res = base | {"grid": " + ".join(l for l, _ in runs), "nodes": nodes, "complete": True}
    Path(out).write_text(json.dumps(res))
    return res


def write_trajectories(merged: dict, base: Path, out: Path):
    """One xyz per m, the lowest structure of every k as a frame, padded so every frame has the atoms of the largest
    k in one element order (trajectory viewers such as Morion then read every frame): an extra atom sits on a real
    atom of its element, or on an anion when the frame has none; n_pad in the comment line."""
    out.mkdir(parents=True, exist_ok=True)
    rows = {}
    for n in merged["nodes"]:
        rows.setdefault(n["m"], []).append(n)
    for m, nodes in sorted(rows.items()):
        frames = []
        for n in sorted(nodes, key=lambda x: x["k"]):
            st = n["states"][0]
            s, p = _read_xyz(Path(base) / st["file"])                   # relative to the merged thermo's folder
            frames.append((n, st, s, p))
        els = []
        for _n, _st, s, _p in frames:
            for e in s:
                if e not in els:
                    els.append(e)
        need = {e: max(s.count(e) for _n, _st, s, _p in frames) for e in els}
        anion = els[1] if len(els) > 1 else els[0]
        text = []
        for n, st, s, p in frames:
            atoms, pad = [], 0
            for e in els:
                pos = [x for a, x in zip(s, p) if a == e]
                host = pos or [x for a, x in zip(s, p) if a == anion]
                for j in range(need[e] - len(pos)):
                    pos.append(host[j % len(host)])
                    pad += 1
                atoms += [(e, x) for x in pos]
            text.append(f"{len(atoms)}\nm={m} k={n['k']} energy={st['E_eV']:.6f} formula={n['formula']} "
                        f"source={st['source'].replace(' ', '')} n_pad={pad}\n" +
                        "\n".join(f"{e:<2s} {x:12.6f} {y:12.6f} {z:12.6f}" for e, (x, y, z) in atoms))
        (out / f"m{m:02d}.xyz").write_text("\n".join(text) + "\n")


def stage_merge(wd: Path, log=print, scheme="legacy"):
    from .grid_thermo import page
    runs = ([("grid", wd / "grid" / "thermo.json")] if scheme == "seeded" else
            [(f"{s}_start", wd / f"grid_{s}" / "thermo.json") for s in STARTS] + [("anneal", wd / "anneal" / "grid" / "thermo.json")])
    merged = merge_thermo(wd / "thermo_merged.json", runs)
    rec = json.loads((wd / "start" / "record.json").read_text())
    s, _p = _read_xyz(wd / "start" / "start.xyz")
    lig = rec["origin"]["recipe"]["passivation"]["ligand"]
    q = rec["origin"]["recipe"]["charges"][next(e for e in s if e not in (lig,))] // -rec["origin"]["recipe"]["charges"][lig]
    start = json.loads((wd / "start.json").read_text())
    page(wd / "thermo_merged.json", wd / "omega.html", title=f"{rec['material']} Grid Grand Potential",
         start=(start["m"], s.count(lig) // q))
    write_trajectories(merged, wd, wd / "xyz")
    if scheme == "seeded":
        solvation_block(wd, log)
    (wd / "merge.json").write_text(json.dumps({"complete": True, "nodes": len(merged["nodes"]),
                                               "states": sum(len(n["states"]) for n in merged["nodes"])}))
    log(f"[nucleation] merged: {len(merged['nodes'])} nodes; omega.html, xyz/")


# ---------------- implicit solvation from the workflow's own structures ----------------
def solvation_block(wd: Path, log=print) -> dict:
    """The thermo_merged.json solvation block, from this workflow only: the geometric surface of every state
    (grid_solvation.add_surfaces) and the COSMO data of the m = 1 complexes (complexes/cosmo.json, the CP2K runs of
    `complexes`; {} until they are collected: the cluster model stands in, never another grid's values)."""
    from .grid_solvation import add_surfaces
    path = wd / "thermo_merged.json"
    th = json.loads(path.read_text())
    add_surfaces(th, wd, log=log)
    f = wd / "complexes" / "cosmo.json"
    cosmo = json.loads(f.read_text()) if f.is_file() else {}
    th["solvation"]["complexes"] = {k: {"E_inf_eV": v["E_inf_eV"], "S_A2": v["S_A2"]} for k, v in cosmo.items()}
    th["solvation"]["complexes_source"] = (f"CP2K COSMO of this grid's m = 1 states ({f.relative_to(wd)})" if cosmo
                                           else "none yet (cluster model)")
    path.write_text(json.dumps(th))
    return th["solvation"]


def complexes_dir(wd: Path) -> Path:
    return Path(wd) / "complexes"


def complexes_prepare(wd, cluster_name: str = "nhr", log=print) -> list:
    """CP2K ESP runs (the campaign monomer library's input, monomers.write_cp2k) of the best state of every m = 1
    node with k >= 1, the precursor complexes MA(MX_q)_k, and a CPU job script running them one after the other on a
    whole node (complexes/job.sh). k = 0 is the MA monomer of the library."""
    from .cluster import load_cluster, sbatch_account
    from .grid_thermo import unit_names
    from .monomers import _xyz, write_cp2k
    wd = Path(wd)
    g = json.loads((wd / "grid" / "grid.json").read_text())
    mx, ma = unit_names(g["nodes"])
    d = complexes_dir(wd)
    d.mkdir(parents=True, exist_ok=True)
    frames = []
    for n in sorted((n for n in g["nodes"] if n["m"] == 1 and n["k"] >= 1), key=lambda n: n["k"]):
        sym, pts = _read_xyz(wd / "grid" / "nodes" / f"m01_k{n['k']:02d}.xyz")
        frames.append(_xyz(sym, pts, f"MA {ma}_{mx}_{n['k']}"))
    (d / "complexes.xyz").write_text("".join(frames))
    cl = load_cluster(cluster_name)
    runs = write_cp2k(str(d / "lib"), str(d / "complexes.xyz"), cl)
    c = cl["cp2k"]
    n = int(cl["cores_per_node"])
    loop = " ".join(str(Path(r).relative_to(d)) for r in runs)
    (d / "job.sh").write_text(f"""#!/bin/bash
#SBATCH -J {ma}-complex-esp
#SBATCH -t 00:40:00
#SBATCH -p {c['test_partition']}
#SBATCH -N 1
#SBATCH --ntasks-per-node={n}
#SBATCH --mem={int(cl['node_memory_gb'])}G
{sbatch_account(cl)}# ESP cubes (COSMO input) of the grid's precursor complexes, one after the other, each on the whole node
source {c['env_script']}
export OMP_NUM_THREADS=1
for d in {loop}; do
    (cd "$d" && {c['launcher']} -n {n} --cpu-bind=cores {c['executable']} -i cp2k.inp -o cp2k.out > srun.log 2>&1)
    echo "$d: $(grep -c "PROGRAM ENDED" $d/cp2k.out) ended"
done
""")
    log(f"[nucleation] complexes: {len(runs)} CP2K ESP runs in {d / 'lib'}; submit {d / 'job.sh'}")
    return runs


def complexes_collect(wd, log=print) -> dict:
    """COSMO (eps -> inf) and isodensity surface of each complex from its ESP cubes (reaction.solvation, as the
    campaign did for its monomers): complexes/cosmo.json {k: ...}; then the thermo's solvation block again."""
    from .reaction import solvation as sv
    wd = Path(wd)
    d = complexes_dir(wd)
    out = {}
    for run in sorted((d / "lib").glob("*/esp")):
        k = run.parent.name.rsplit("_", 1)[-1]
        cubes = sv.esp_cubes(run)
        if not cubes:
            raise RuntimeError(f"{run}: no ESP cubes (CP2K not run or not finished)")
        r = sv.cosmo_isodensity(sv.read_cube(cubes[0]), sv.read_cube(cubes[1]))
        out[k] = {"E_inf_eV": r["E_inf_eV"], "S_A2": r["S_A2"], "V_A3": r["V_A3"], "n_points": r["n_points"],
                  "run": str(run.relative_to(wd))}
        log(f"[nucleation] complex k = {k}: E_inf {r['E_inf_eV']:.4f} eV, S {r['S_A2']:.1f} A2")
    (d / "cosmo.json").write_text(json.dumps(out, indent=1))
    if (wd / "thermo_merged.json").is_file():
        solvation_block(wd, log)
    return out


def run(wd: str, model: str, head: str = "omat_pbe", refs: str | None = None, log=print):
    """Every unfinished stage in turn (grid, anneal, ensemble, thermo, merge); a rerun continues where it stopped."""
    import time
    from .lattice_grid import EVALS, count_evals
    wd = Path(wd)
    P = _block(json.loads((wd / "nucleation.json").read_text()))
    scheme = P["search"]["scheme"]
    s = _settings(model, head)
    if P.get("ligand"):
        from .ligands import configure
        configure(s, P["ligand"])
        log(f"[nucleation] ligand {P['ligand']}: carboxylate X on the Cl placeholders, D3 {P['ligand'].get('d3', True)}")
    count_evals()
    tf = wd / "timings.json"
    timings = json.loads(tf.read_text()) if tf.is_file() else []

    def timed(name, f, *a, **kw):
        t0, ev0 = time.time(), dict(EVALS)
        f(*a, **kw)
        timings.append({"stage": name, "job": os.environ.get("SLURM_JOB_ID"), "seconds": round(time.time() - t0, 1),
                        "force_evals": EVALS["structures"] - ev0["structures"],
                        "atom_force_evals": EVALS["atoms"] - ev0["atoms"]})
        tf.write_text(json.dumps(timings, indent=1))

    timed("grid", stage_grid, wd, s, P, log, model, head, refs)
    if scheme == "legacy":
        timed("anneal", stage_anneal, wd, P, model, head, log)
        stage_ensemble(wd, log)
    timed("thermo", stage_thermo, wd, s, refs, log, scheme)
    if not _done(wd / "merge.json"):
        stage_merge(wd, log, scheme)
    log("[nucleation] done")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python -m qd_workflows.nucleation")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("start", help="find the smallest Wulff dot with discernible facets; write the workflow folder")
    p.add_argument("config")
    p.add_argument("--out", required=True)
    p = sub.add_parser("complexes", help="CP2K COSMO of the grid's m = 1 complexes: prepare (inputs + job.sh), collect")
    p.add_argument("action", choices=("prepare", "collect"))
    p.add_argument("wd")
    p.add_argument("--cluster", default="nhr")
    p = sub.add_parser("seeds", help="the two seeds of an existing start (builder): DIR/start/seed_*.xyz, seeds.json")
    p.add_argument("wd")
    p = sub.add_parser("run", help="the unfinished stages in turn (cluster GPU)")
    p.add_argument("wd")
    p.add_argument("--out", help="ignored (the workflow folder is wd); there for qd_workflows.grid.prepare")
    p.add_argument("--model", default=os.environ.get("QDPROPS_MACE_MODEL"))
    p.add_argument("--head", default="omat_pbe")
    p.add_argument("--refs", default=os.environ.get("QDPROPS_REFS"), help="references dir with the monomer library")
    p = sub.add_parser("prepare", help="GPU job script (self-restarting pieces)")
    p.add_argument("wd")
    p.add_argument("--cluster", default="default")
    p.add_argument("--code", required=True)
    p.add_argument("--partition")
    p.add_argument("--time", default="01:00:00")
    p.add_argument("--pieces", action="store_true")
    a = ap.parse_args(argv)
    if a.cmd == "start":
        out = Path(a.out)
        out.mkdir(parents=True, exist_ok=True)
        cfg = yaml.safe_load(Path(a.config).read_text())
        (out / "nucleation.json").write_text(json.dumps({"nucleation": cfg.get("nucleation", {}), "config": str(Path(a.config).resolve())}))
        r = find_start(a.config, a.out)
        print(json.dumps({k: v for k, v in r.items() if k != "tried"}, default=float, indent=1))
        if _block(cfg)["search"]["scheme"] == "seeded":
            build_seeds(out / "start")
    elif a.cmd == "seeds":
        build_seeds(Path(a.wd) / "start")
    elif a.cmd == "complexes":
        complexes_prepare(a.wd, a.cluster) if a.action == "prepare" else complexes_collect(a.wd)
    elif a.cmd == "prepare":
        from .cluster import load_cluster
        from .grid import prepare
        cl = load_cluster(a.cluster)
        extra = f"--refs {cl['props']['refs']}" if cl["props"].get("refs") else ""
        print(json.dumps(prepare(a.wd, a.wd, cl, a.code, time_=a.time, partition=a.partition, pieces=a.pieces,
                                 extra=extra, module="qd_workflows.nucleation"), indent=1))
    else:
        import signal

        def stop(signum, _frame):
            print(f"[nucleation] signal {signum}: stopping (stages resume from their caches)", flush=True)
            sys.exit(128 + signum)
        signal.signal(signal.SIGTERM, stop)
        run(a.wd, a.model, a.head, a.refs, log=lambda x: print(x, flush=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
