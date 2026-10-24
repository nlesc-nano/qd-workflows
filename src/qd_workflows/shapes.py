"""
Shape library of one material from QD_builder, for fitting facet energies, edge tensions and ligand binding.

Sweeps the builder recipe over facet surface-energy ratios (each facet family switched off or given one of a few
gamma values), cell sizes (isotropic and anisotropic: platelets, rods) and construction centres, keeps the unique
structures (builder fingerprint), and adds to each the geometry of its Wulff polyhedron, rebuilt from the builder's
own construction planes: facet areas by family and termination, edges by facet pair with their lengths, vertices.
`estimate` turns the library into a count of calculations and A100 hours for the surface-energy campaign.

    python -m qd_workflows.shapes sweep    <config.yaml> --out DIR [--jobs N] [--limit N]
    python -m qd_workflows.shapes select   <config.yaml> LIB [LIB ...] --out set.json
    python -m qd_workflows.shapes estimate LIB [LIB ...] [--set set.json] [--sites 3] [--candidates 4]

Config (see templates/shapes/cdse_zb.yaml):
    material, phase, cif (relative to the config), native_elements, centres, positive_q_mode,
    sizes: {start, stop, step}           unit cells, the isotropic size (smallest axis for anisotropic cells)
    aspects: [[1, 1, 1], [1, 1, 0.5]]    cell multipliers per axis (platelets < 1 along one axis, rods > 1)
    max_atoms, min_core_atoms, min_interatomic_distance
    preset: {charges, passivation, symmetry}
    facets: [{hkl, termination?, gamma: [null, 1.0, 1.25]}]   null = family absent
"""
from __future__ import annotations

import argparse
import copy
import csv
import itertools
import json
import os
import tempfile
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import yaml

SEC_PER_ATOM_RELAX = 0.022     # s per atom per relaxation, batched MACE-MH-1 on one A100 (Cd16 grid annealing quenches)


# ---------------- recipes ----------------
def gamma_sets(facets: list[dict]) -> list[list[dict]]:
    """Every combination of the families' gamma values (null = absent), one per distinct set of ratios."""
    seen, out = set(), []
    for values in itertools.product(*[f["gamma"] for f in facets]):
        on = [(f, g) for f, g in zip(facets, values) if g is not None]
        if not on:
            continue
        gmin = min(g for _f, g in on)
        key = tuple(sorted((f["hkl"], f.get("termination") or "", round(g / gmin, 4)) for f, g in on))
        if key in seen:
            continue
        seen.add(key)
        out.append([{"hkl": f["hkl"], "scope": "family", **({"termination": f["termination"]} if f.get("termination") else {}),
                     "gamma": round(g / gmin, 4)} for f, g in on])
    return out


def recipes(cfg: dict):
    """facets: a gamma grid (every combination, gamma_sets); facet_sets: explicit shapes, each a list of facets."""
    start, stop, step = (float(cfg["sizes"][k]) for k in ("start", "stop", "step"))
    sizes = [round(start + i * step, 4) for i in range(int(round((stop - start) / step)) + 1)]
    sets = gamma_sets(cfg["facets"]) if "facets" in cfg else []
    sets += [[{"scope": "family", **f} for f in s["facets"]] for s in cfg.get("facet_sets", [])]
    for facets in sets:
        for aspect in cfg.get("aspects", [[1, 1, 1]]):
            for size in sizes:
                for centre in cfg["centres"]:       # a species, or {fractional_shift: [...]} / {cartesian_shift: [...]}
                    r = copy.deepcopy(cfg["preset"])
                    r["facets"] = facets
                    r["size_unit_cells"] = [round(size * float(a), 4) for a in aspect]
                    if isinstance(centre, dict):
                        r["construction_origin"] = dict(centre)
                        label = ",".join(f"{k.split('_')[0]}({' '.join(f'{x:g}' for x in v)})" for k, v in centre.items())
                    else:
                        r["construction_origin"] = {"center_on_species": centre}
                        label = centre
                    yield r, {"size": size, "aspect": list(map(float, aspect)), "centre": label}


# ---------------- geometry of the construction polyhedron ----------------
def family_label(h: int, k: int, l: int, termination: str | None) -> str:
    fam = "".join(str(x) for x in sorted((abs(h), abs(k), abs(l)), reverse=True))
    if fam == "111":
        return "111A" if termination == "cation_rich" else "111B" if termination == "anion_rich" else "111"
    return fam


def construction_planes(cif: str, recipe: dict):
    """The builder's Wulff planes for this recipe (terminations resolved, symmetry expanded, auto-shifted off atom
    layers), with a family label per plane. Same calls as builder.main's first construction attempt."""
    from pymatgen.core import Structure
    from builder.config import _parse_hkl
    from builder.facets import expand_facets
    from builder.geometry import build_nanocrystal
    from builder.main import _resolve_facet_terminations, construction_origin_shift
    from builder.nc_types import Facet
    from builder.stack import size_unit_cells_to_radius_aspect

    struct = Structure.from_file(cif)
    seeds = [Facet(*_parse_hkl(str(f["hkl"])), float(f["gamma"]), termination=f.get("termination"),
                   scope=f.get("scope", "family")) for f in recipe["facets"]]
    seeds = _resolve_facet_terminations(struct, seeds, recipe["charges"])
    wulff = expand_facets(struct, seeds, proper_only=bool(recipe.get("symmetry", {}).get("proper_rotations_only", True)))
    radius, aspect = size_unit_cells_to_radius_aspect(struct, tuple(recipe["size_unit_cells"]))
    shift = construction_origin_shift(struct, recipe.get("construction_origin"))
    _s, _p, planes = build_nanocrystal(struct, wulff, radius, aspect=aspect, origin_shift=shift)
    return [(np.asarray(n, float), float(d), family_label(f.h, f.k, f.l, f.termination)) for (n, d), f in zip(planes, wulff)]


def polyhedron(planes, tol: float = 1e-6) -> dict:
    """Facet areas by family, edges by facet pair (count, total length), vertices, volume of the half-space
    intersection n.x <= d (origin inside)."""
    from scipy.spatial import ConvexHull, HalfspaceIntersection
    hs = np.array([np.append(n, -d) for n, d, _l in planes])
    verts = HalfspaceIntersection(hs, np.zeros(3)).intersections
    verts = np.unique(np.round(verts, 6), axis=0)
    on = [np.flatnonzero(np.abs(verts @ n - d) < 1e-4 * max(1.0, d)) for n, d, _l in planes]
    areas, faces = defaultdict(float), []
    for (n, d, lab), idx in zip(planes, on):
        if len(idx) < 3:
            continue
        pts = verts[idx]
        c = pts.mean(0)
        u = pts[0] - c
        u /= np.linalg.norm(u)
        w = np.cross(n, u)
        order = np.argsort(np.arctan2((pts - c) @ w, (pts - c) @ u))
        p = pts[order]
        a = 0.5 * np.linalg.norm(np.cross(p, np.roll(p, -1, axis=0)).sum(0))
        if a > tol:
            areas[lab] += a
            faces.append((lab, set(idx.tolist())))
    edges = defaultdict(lambda: [0, 0.0])
    for (la, ia), (lb, ib) in itertools.combinations(faces, 2):
        common = sorted(ia & ib)
        if len(common) >= 2:
            seg = verts[common]
            length = max(np.linalg.norm(x - y) for x, y in itertools.combinations(seg, 2))
            key = "|".join(sorted((la, lb)))
            edges[key][0] += 1
            edges[key][1] += float(length)
    used = sorted(set().union(*[s for _l, s in faces])) if faces else []
    vtypes = Counter("+".join(sorted({lab for lab, s in faces if v in s})) for v in used)
    return {"areas_A2": {k: round(v, 3) for k, v in sorted(areas.items())},
            "edges": {k: {"n": v[0], "length_A": round(v[1], 3)} for k, v in sorted(edges.items())},
            "n_vertices": len(used), "vertex_types": dict(vtypes),
            "volume_A3": round(float(ConvexHull(verts).volume), 3), "n_faces": len(faces)}


# ---------------- one build ----------------
def build_one(job):
    cif, recipe, meta, posq, native, limits = job
    from builder.library_record import bulk_volume_per_atom, describe_structure, min_distance
    from builder.scripts.generate_library import _build
    row = {**meta, "facets": recipe["facets"], "status": ""}
    with tempfile.TemporaryDirectory() as tmp:
        res = _build(cif, recipe, posq, Path(tmp) / "b")
    if "error" in res:
        return row | {"status": "rejected: " + res["error"].splitlines()[0][:160]}
    desc = describe_structure(res["symbols"], res["pts"], native_order=native, charges=recipe["charges"],
                              centre_index=res["centre_index"], bulk_volume=bulk_volume_per_atom(cif))
    syms = list(res["symbols"])
    cation, anion = native
    lig = [e for e in set(syms) if e not in native]
    m, k = syms.count(anion), syms.count(cation) - syms.count(anion)
    row.update(formula=desc["formula"], n_atoms=desc["n_atoms"], Q=desc["total_charge"], m=m, k=k,
               n_ligand=sum(syms.count(e) for e in lig), fingerprint=desc["fingerprint"], d_nm=desc["size"]["d_nm"])
    reasons = []
    if desc["total_charge"] != 0:
        reasons.append(f"charge {desc['total_charge']:+d}")
    if desc["size"]["n_core"] < limits["min_core_atoms"]:
        reasons.append(f"core {desc['size']['n_core']} atoms")
    if desc["n_atoms"] > limits["max_atoms"]:
        reasons.append(f"{desc['n_atoms']} atoms > {limits['max_atoms']}")
    if len(res["pts"]) and min_distance(res["pts"]) < limits["min_interatomic_distance"]:
        reasons.append("atom clash")
    if reasons:
        return row | {"status": "rejected: " + "; ".join(reasons)}
    try:
        geo = polyhedron(construction_planes(cif, recipe))
        vol_atoms = desc["size"]["n_core"] * next(iter(bulk_volume_per_atom(cif).values()))   # same for every element
        geo["volume_ratio"] = round(geo["volume_A3"] / vol_atoms, 3)    # construction polyhedron / volume of the core atoms
    except Exception as exc:  # keep the structure; flag the missing geometry
        geo = {"error": f"{type(exc).__name__}: {exc}"}
    return row | {"status": "built", "geometry": geo, "symbols": syms, "pts": np.round(res["pts"], 5).tolist()}


# ---------------- sweep ----------------
def sweep(config: str, out: str, jobs: int | None = None, limit: int | None = None, positive_q_mode: str | None = None) -> dict:
    cfg_path = Path(config).resolve()
    cfg = yaml.safe_load(cfg_path.read_text())
    if positive_q_mode:
        cfg["positive_q_mode"] = positive_q_mode
    cif = str((cfg_path.parent / cfg["cif"]).resolve())
    native = list(cfg["native_elements"])
    limits = {"max_atoms": int(cfg.get("max_atoms", 2000)), "min_core_atoms": int(cfg.get("min_core_atoms", 20)),
              "min_interatomic_distance": float(cfg.get("min_interatomic_distance", 1.8))}
    todo = [(cif, r, meta, cfg.get("positive_q_mode", "add"), native, limits) for r, meta in recipes(cfg)]
    if limit:
        todo = todo[:limit]
    out = Path(out)
    (out / "structures").mkdir(parents=True, exist_ok=True)
    rows, kept = [], {}
    with ProcessPoolExecutor(max_workers=jobs or os.cpu_count()) as ex:
        for i, row in enumerate(ex.map(build_one, todo, chunksize=4)):
            if row["status"] == "built":
                fp = row["fingerprint"]
                if fp in kept:
                    kept[fp]["also"].append({k: row[k] for k in ("size", "aspect", "centre", "facets")})
                    row = {k: v for k, v in row.items() if k not in ("symbols", "pts")} | {"status": f"duplicate of {kept[fp]['id']}"}
                else:
                    sid = f"{cfg['material']}-{row['formula']}-{len(kept):05d}"
                    kept[fp] = {"id": sid, "row": row, "also": []}
                    row = {k: v for k, v in row.items() if k not in ("symbols", "pts")} | {"status": "kept", "id": sid}
            rows.append(row)
            if (i + 1) % 200 == 0:
                print(f"[shapes] {i + 1}/{len(todo)} builds, {len(kept)} unique", flush=True)
    for fp, rec in kept.items():
        r = rec["row"]
        d = out / "structures" / rec["id"]
        d.mkdir(exist_ok=True)
        d.joinpath("start.xyz").write_text(f"{len(r['symbols'])}\n{rec['id']}\n" + "\n".join(
            f"{s} {x:.5f} {y:.5f} {z:.5f}" for s, (x, y, z) in zip(r["symbols"], r["pts"])) + "\n")
        rec_json = {k: v for k, v in r.items() if k not in ("symbols", "pts", "status")} | {
            "id": rec["id"], "material": cfg["material"], "phase": cfg.get("phase"),
            "positive_q_mode": cfg.get("positive_q_mode", "add"), "also_built_by": rec["also"]}
        d.joinpath("record.json").write_text(json.dumps(rec_json, indent=1))
    cols = ["status", "id", "formula", "m", "k", "n_atoms", "size", "aspect", "centre", "facets"]
    with open(out / "sweep.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r | {"facets": ";".join(f"{f['hkl']}{(':' + f['termination'][0]) if f.get('termination') else ''}={f['gamma']}"
                                                for f in r["facets"])})
    status = Counter(r["status"].split(":")[0].split(" of ")[0] for r in rows)
    summary = {"config": str(cfg_path), "builds": len(rows), "status": dict(status), "unique": len(kept),
               "unique_mk": len({(rec["row"]["m"], rec["row"]["k"]) for rec in kept.values()}),
               "structures": [json.loads((out / "structures" / rec["id"] / "record.json").read_text()) for rec in kept.values()]}
    (out / "summary.json").write_text(json.dumps(summary))
    print(f"[shapes] {len(rows)} builds: {dict(status)}; {len(kept)} unique structures, {summary['unique_mk']} distinct (m, k)")
    return summary


# ---------------- selection ----------------
def design(S: list[dict]):
    """Geometry matrix of the fit sum_i gamma_i A_i + sum_j tau_j L_j + eps_v N_v: one column per facet family,
    per edge type and the vertex count; rows: structures with a geometry."""
    S = [s for s in S if "areas_A2" in s.get("geometry", {})]
    fams = sorted({f for s in S for f in s["geometry"]["areas_A2"]})
    etypes = sorted({e for s in S for e in s["geometry"]["edges"]})
    cols = [f"A_{f}" for f in fams] + [f"L_{e}" for e in etypes] + ["N_v"]
    X = np.array([[s["geometry"]["areas_A2"].get(f, 0.0) for f in fams] +
                  [s["geometry"]["edges"].get(e, {}).get("length_A", 0.0) for e in etypes] +
                  [s["geometry"]["n_vertices"]] for s in S], float)
    return S, cols, X


SELECTION = {   # defaults of the `selection:` block of a material config
    "exposed_share": 0.02,                                       # a family counts as exposed above this area share
    "fit": {"tiers": [{"max_atoms": 300, "n": 60}, {"max_atoms": 600, "n": 30}], "per_shape": 6},
    "groups": {"max_atoms": 300, "n_pure": 5, "min_pure": 3, "dominated_share": 0.6, "n_dominated": 4},
    "strip": {"per_group": 2, "mixed": 2, "atoms": [60, 200], "mixed_min_families": 3},
}


def load_libraries(dirs: list[str]):
    """Structures of the libraries in the given order, duplicates (fingerprint) after the first dropped; and a
    sha256 of each library's summary.json, so a selection names exactly the libraries it was made from."""
    import hashlib
    S, seen, digests = [], set(), {}
    for d in dirs:
        raw = (Path(d) / "summary.json").read_bytes()
        digests[str(Path(d).resolve())] = hashlib.sha256(raw).hexdigest()
        for s in json.loads(raw)["structures"]:
            if s["fingerprint"] not in seen and "areas_A2" in s.get("geometry", {}):
                seen.add(s["fingerprint"])
                S.append(s | {"library": str(Path(d).resolve())})
    return S, digests


def exposed_families(s: dict, share: float) -> list[str]:
    A = s["geometry"]["areas_A2"]
    tot = sum(A.values()) or 1.0
    return sorted(f for f, a in A.items() if a / tot >= share)


def select_set(dirs: list[str], params: dict | None = None, out: str | None = None) -> dict:
    """The dots of a surface-energy campaign, by fixed rules (same libraries and parameters -> same set):
    fit      greedy D-optimal per size tier on the geometry matrix (facet areas, edge lengths by type, vertex count;
             columns normalised): each step adds the dot that most increases det(X^T X); ties go to the earlier dot;
             at most per_shape dots per exposed-facet combination;
    groups   facet-resolved binding: one group per single-family shape, and per pair of polar terminations of one
             family (e.g. 111A + 111B, octahedra); a family with fewer than min_pure pure dots gets a 'dominated'
             group instead (area share >= dominated_share); n_pure / n_dominated dots spread over size;
    strip    per_group dots of every group and `mixed` fit dots with >= mixed_min_families families, spread over size
             within the atom window: stripped to the bare core.
    Group names come from the facet labels of the library, so the rules apply to any material."""
    P = json.loads(json.dumps(SELECTION))
    for k, v in (params or {}).items():
        P[k] = P[k] | v if isinstance(v, dict) and isinstance(P.get(k), dict) else v
    S, digests = load_libraries(dirs)
    _S, cols, X = design(S)
    Xn = X / np.maximum(np.abs(X).max(0), 1e-12)
    fam = [exposed_families(s, P["exposed_share"]) for s in S]
    shape = ["+".join(f) for f in fam]
    atoms = np.array([s["n_atoms"] for s in S])

    def dopt(pool, n, start, per_shape):
        chosen, M = list(start), 1e-6 * np.eye(X.shape[1])
        for i in chosen:
            M += np.outer(Xn[i], Xn[i])
        count = Counter(shape[i] for i in chosen)
        new = []
        for _ in range(n):
            cand = [i for i in pool if i not in chosen and count[shape[i]] < per_shape]
            if not cand:
                break
            Minv = np.linalg.inv(M)
            gain = [Xn[i] @ Minv @ Xn[i] for i in cand]
            b = cand[int(np.argmax(gain))]                       # argmax: first of equal gains
            chosen.append(b)
            new.append(b)
            count[shape[b]] += 1
            M += np.outer(Xn[b], Xn[b])
        return new

    def spread(pool, n):
        pool = sorted(pool, key=lambda i: (atoms[i], S[i]["id"]))
        if len(pool) <= n:
            return pool
        return [pool[int(round(j * (len(pool) - 1) / (n - 1)))] for j in range(n)] if n > 1 else pool[len(pool) // 2: len(pool) // 2 + 1]

    reason, fit, lo = {}, [], 0
    for t, tier in enumerate(P["fit"]["tiers"]):
        pool = [i for i in range(len(S)) if lo < atoms[i] <= tier["max_atoms"]]
        new = dopt(pool, tier["n"], fit, P["fit"]["per_shape"])
        fit += new
        for i in new:
            reason[i] = [f"fit tier {t + 1}"]
        lo = tier["max_atoms"]
    G = P["groups"]
    small = [i for i in range(len(S)) if atoms[i] <= G["max_atoms"]]
    base = lambda f: f.rstrip("AB")
    families = sorted({f for fs in fam for f in fs})
    groups = {}
    for f in families:
        groups[f"pure {f}"] = [i for i in small if fam[i] == [f]]
    for a, b in __import__("itertools").combinations(families, 2):
        if base(a) == base(b) and a != b:
            groups[f"pure {a}+{b}"] = [i for i in small if fam[i] == [a, b]]
    for f in families:
        if len(groups[f"pure {f}"]) < G["min_pure"]:
            del groups[f"pure {f}"]
            groups[f"dominated {f}"] = [i for i in small if S[i]["geometry"]["areas_A2"].get(f, 0) /
                                        sum(S[i]["geometry"]["areas_A2"].values()) >= G["dominated_share"]]
    picked = {}
    for g in sorted(groups):
        picked[g] = spread(groups[g], G["n_pure"] if g.startswith("pure") else G["n_dominated"])
        for i in picked[g]:
            reason.setdefault(i, []).append(f"binding {g}")
    st = P["strip"]
    lo_a, hi_a = st["atoms"]
    win = lambda pool: [i for i in pool if lo_a <= atoms[i] <= hi_a]
    strip = []
    for g in sorted(groups):
        strip += [i for i in spread(win(groups[g]), st["per_group"]) if i not in strip]
    strip += [i for i in spread(win([i for i in fit if len(fam[i]) >= st["mixed_min_families"]]), st["mixed"]) if i not in strip]
    for i in strip:
        reason.setdefault(i, []).append("strip")
    chosen = sorted(reason, key=lambda i: (atoms[i], S[i]["id"]))
    res = {"parameters": P, "libraries": digests, "columns": cols, "n": len(chosen), "n_strip": len(strip),
           "condition": float(np.linalg.cond(Xn[chosen])), "groups": {g: len(v) for g, v in groups.items()},
           "ids": [S[i]["id"] for i in chosen], "strip": [S[i]["id"] for i in strip],
           "dots": [{"id": S[i]["id"], "library": S[i]["library"], "n_atoms": int(atoms[i]), "m": S[i]["m"], "k": S[i]["k"],
                     "shape": shape[i], "reason": reason[i]} for i in chosen]}
    if out:
        Path(out).write_text(json.dumps(res, indent=1))
    print(f"[shapes] set: {res['n']} dots ({len(fit)} fit, groups {dict((g, len(v)) for g, v in picked.items())}), "
          f"{len(strip)} to strip; condition {res['condition']:.1f}")
    return res


# ---------------- estimate ----------------
def estimate(dirs: list[str], set_path: str | None = None, sites: int = 3, candidates: int = 4, out: str | None = None) -> dict:
    """Calculations of the sampled ladders (qd_workflows.shape_ladder) for a set (or every dot of the libraries,
    each stripped):
    relax       one relaxation per dot, polished (counted twice);
    binding     `sites` CdCl2 removals per site class (facet families, edge types, vertices), polished;
    adsorption  `sites` Z-type CdCl2 additions per site class, polished;
    strip       strip dots only: k steps of `candidates` removals plus the polish of the step's best."""
    S, _d = load_libraries(dirs)
    if set_path:
        sel = json.loads(Path(set_path).read_text())
        ids, strip = set(sel["ids"]), set(sel["strip"])
    else:
        ids = strip = {s["id"] for s in S}
    tot, rows = defaultdict(float), []
    for s in S:
        if s["id"] not in ids:
            continue
        g = s["geometry"]
        n_cls = len(g.get("areas_A2", {})) + len(g.get("edges", {})) + (1 if g.get("n_vertices") else 0)
        n, k = s["n_atoms"], s["k"]
        relax, bind, ads = 2, 2 * n_cls * sites, 2 * n_cls * sites
        st = (candidates + 1) * k if s["id"] in strip else 0
        atoms_st = sum((candidates + 1) * (n - 3 * j) for j in range(k)) if s["id"] in strip else 0
        rows.append({"id": s["id"], "n_atoms": n, "k": k, "classes": n_cls, "strip": s["id"] in strip})
        for key, v in (("relax", relax), ("binding", bind), ("adsorption", ads), ("strip", st)):
            tot[key] += v
        tot["atom_relax"] += n * (relax + bind) + (n + 3) * ads + atoms_st
    hours = tot["atom_relax"] * SEC_PER_ATOM_RELAX / 3600
    res = {"dots": len(rows), "strip": sum(r["strip"] for r in rows),
           "relaxations": {k: int(v) for k, v in tot.items() if k != "atom_relax"},
           "a100_hours": round(hours, 1), "core_hours_shared": round(hours * 150),
           "assumptions": {"sites_per_class": sites, "candidates_per_step": candidates, "sec_per_atom_relax": SEC_PER_ATOM_RELAX,
                           "polish": "relax, binding and adsorption counted twice"}}
    if out:
        Path(out).write_text(json.dumps(res | {"per_dot": rows}, indent=1))
    print(json.dumps(res, indent=1))
    return res


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python -m qd_workflows.shapes")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("sweep", help="build the shape library of one material")
    p.add_argument("config")
    p.add_argument("--out", required=True)
    p.add_argument("--jobs", type=int)
    p.add_argument("--limit", type=int, help="build only the first N recipes (timing tests)")
    p.add_argument("--positive-q-mode", choices=("add", "remove"), help="override the config (add: Cl added; remove: Cd removed)")
    p = sub.add_parser("select", help="the campaign's dot set by fixed rules (config `selection:` block)")
    p.add_argument("config", help="material config; its `selection:` block overrides the defaults")
    p.add_argument("libs", nargs="+", help="library directories, in a fixed order (earlier wins on duplicates)")
    p.add_argument("--out", required=True)
    p = sub.add_parser("estimate", help="calculations and A100 hours of the sampled ladders")
    p.add_argument("libs", nargs="+")
    p.add_argument("--set", help="set file from `select` (default: every dot, each stripped)")
    p.add_argument("--sites", type=int, default=3)
    p.add_argument("--candidates", type=int, default=4)
    p.add_argument("--out")
    a = ap.parse_args(argv)
    if a.cmd == "sweep":
        sweep(a.config, a.out, a.jobs, a.limit, a.positive_q_mode)
    elif a.cmd == "select":
        select_set(a.libs, yaml.safe_load(Path(a.config).read_text()).get("selection"), a.out)
    else:
        estimate(a.libs, a.set, a.sites, a.candidates, a.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
