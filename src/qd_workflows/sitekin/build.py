"""Step 1, build: the dot of each size from QD_builder, with the builder's own site role of every atom.

The roles come from the build manifest (site_roles.roles, from builder.polyhedron on the as-built dot: "facet:<f>",
"edge:<f>|<g>", "vertex:<...>", "interior") and are carried by atom through the surface treatment; planes are never
re-fitted on a modified dot. Surface treatments:
  plain      the builder's Cl charge balance
  swap       the builder's Cl swap reconstruction
  half_swap  swap, then (as build_half.py of the 2026-10 runs)
               1. every Cl on {111}A removed
               2. {111}A Cd left three-coordinated paired by bridging Cl (off the pair's middle by 0.5 A toward the
                  hollow with no Se below, Cl-Cd 2.55 A; maximum matching, bridges >= 3.6 A apart); a bridge Cl takes
                  the role of the Cd it is nearest to (a tie between the
                  two Cd of the pair is broken by the k-d tree; only Cl roles depend on it, sites are picked on Cd)
               3. on {111}B as many swapped Cl kept as keep the dot neutral, spread evenly over the facets and,
                  within a facet, farthest-point apart from its rim; the others go back to Se
Writes dots/<tag>.json: symbols, positions (as built, centred on the builder dot's mean), roles, Cd per role, the recipe used."""
from __future__ import annotations

import contextlib
import io
import json
import tempfile
from collections import Counter
from pathlib import Path

import numpy as np
import yaml
from scipy.spatial import cKDTree

from . import recipe as rc

RECON = {"surface_reconstruction": {"enabled": True, "ligand": "Cl", "anion_pattern": "swap"}}


def builder_recipe(r: dict, size: float, surface: str) -> dict:
    b = r["builder"]
    out = {"charges": r["material"]["charges"], "passivation": b["passivation"], "symmetry": b["symmetry"],
           "construction_origin": {"center_on_species": b["centre"]},
           "facets": [{"scope": "family", **f} for f in b["facets"]], "size_unit_cells": [float(size)] * 3}
    if surface in ("swap", "half_swap"):
        out["post_treatment"] = RECON
    return out


def run_builder(cif: str, rec: dict, positive_q_mode: str = "add"):
    """QD_builder on one recipe: (symbols, positions, manifest)."""
    from builder.main import main as builder_main
    with tempfile.TemporaryDirectory() as tmp:
        yml, xyz = Path(tmp) / "r.yaml", Path(tmp) / "out.xyz"
        yml.write_text(yaml.safe_dump(rec, sort_keys=False))
        log = io.StringIO()
        with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
            rc_ = builder_main([cif, str(yml), "-o", str(xyz), "--positive-q-mode", positive_q_mode])
        if rc_ not in (0, None) or not xyz.exists():
            raise RuntimeError(f"QD_builder failed ({rc_}): {log.getvalue()[-500:]}")
        man = json.loads((Path(tmp) / "out.json").read_text())
        lines = xyz.read_text().splitlines()
    n = int(lines[0])
    s = np.array([ln.split()[0] for ln in lines[2:2 + n]])
    x = np.array([[float(v) for v in ln.split()[1:4]] for ln in lines[2:2 + n]])
    return s, x, man


def facet_planes(man: dict, s: np.ndarray, x: np.ndarray, core: np.ndarray):
    """Faces of the builder's atom-fitted polyhedron through the outermost core atom: (normals, labels) per face."""
    from builder.polyhedron import clip
    planes = man["facet_geometry"]["atom_fitted"]["planes"]
    N = np.array([p["normal"] for p in planes], float)
    lab = [p["label"] for p in planes]
    d = (x[core] @ N.T).max(0)
    poly = clip([(n, v) for n, v in zip(N, d)], lab)
    keys = sorted(poly.faces)
    return np.array([poly.normals[k] for k in keys], float), [poly.labels[k] for k in keys]


def half_swap(s: np.ndarray, x: np.ndarray, N: np.ndarray, lab: list[str], log=print):
    """The half-swap dot from the swap-reconstructed one (see the module text). Returns (symbols, positions, src):
    src[i] is the index in the swap dot of new atom i, or -1 for a bridging Cl."""
    import networkx as nx
    d_pl = (x[s != "Cl"] @ N.T).max(0)
    t = cKDTree(x)
    nb = [[j for j in t.query_ball_point(x[i], 3.2) if j != i] for i in range(len(s))]

    def plane(i):
        return int(np.argmin(d_pl - x[i] @ N.T))

    cl = np.flatnonzero(s == "Cl")
    ncd = {int(k): sum(s[j] == "Cd" for j in nb[k]) for k in cl}
    on_a = [int(k) for k in cl if lab[plane(k)] == "111A" and (ncd[k] == 1 or (d_pl[plane(k)] - x[k] @ N[plane(k)]) < -0.5)]
    on_b = [int(k) for k in cl if lab[plane(k)] == "111B" and ncd[k] >= 3]
    gone = set(on_a)
    shift = 0.5
    planes_a = [p for p in range(len(lab)) if lab[p] == "111A"]

    def a_plane(i):
        g = [p for p in planes_a if abs(d_pl[p] - x[i] @ N[p]) < 1.0]
        return g[0] if len(g) == 1 else None

    def coord_after(i):
        return sum(1 for j in nb[i] if s[j] == "Se" or (s[j] == "Cl" and j not in gone))

    elig = [int(i) for i in np.flatnonzero(s == "Cd") if a_plane(i) is not None and coord_after(i) == 3]
    keep_idx = [i for i in range(len(s)) if i not in gone]
    tk = cKDTree(x[keep_idx])
    se_x = x[s == "Se"]
    G = nx.Graph()
    site = {}
    for a in elig:
        for b in elig:
            if a < b and a_plane(a) == a_plane(b) and np.linalg.norm(x[a] - x[b]) < 4.7:
                n_ = N[a_plane(a)]
                m = (x[a] + x[b]) / 2
                u = np.cross(n_, x[b] - x[a])
                u /= np.linalg.norm(u)
                if np.linalg.norm(se_x - (m + u), axis=1).min() < np.linalg.norm(se_x - (m - u), axis=1).min():
                    u = -u
                lat = np.hypot(np.linalg.norm(x[b] - x[a]) / 2, shift)
                p = m + shift * u + np.sqrt(max(2.55 ** 2 - lat ** 2, 0.64)) * n_
                near = [keep_idx[j] for j in tk.query_ball_point(p, 2.7)]
                if set(near) <= {a, b}:
                    G.add_edge(a, b)
                    site[(a, b)] = p
    match = sorted(tuple(sorted(e)) for e in nx.max_weight_matching(G, maxcardinality=True))
    bridges = []
    for e in match:
        p = site[e]
        if all(np.linalg.norm(p - q) >= 3.6 for _e, q in bridges):
            bridges.append((e, p))
    n_cd, n_se, n_cl = int((s == "Cd").sum()), int((s == "Se").sum()), int((s == "Cl").sum())
    q_now = 2 * n_cd - 2 * n_se - (n_cl - len(on_a) + len(bridges))
    keep_sw = len(on_b) - q_now
    if not 0 <= keep_sw <= len(on_b):
        raise ValueError(f"half swap: cannot neutralise (charge {q_now:+d}, {len(on_b)} {{111}}B swaps)")
    byf = {}
    for k in on_b:
        byf.setdefault(plane(k), []).append(k)
    faces = sorted(byf)
    quota = {f: keep_sw // len(faces) + (1 if i < keep_sw % len(faces) else 0) for i, f in enumerate(faces)}
    kept = []
    for f in faces:
        cand = byf[f]
        c0 = x[cand].mean(0)
        sel = [max(cand, key=lambda k: np.linalg.norm(x[k] - c0))] if quota[f] else []
        while len(sel) < quota[f]:
            sel.append(max((k for k in cand if k not in sel), key=lambda k: min(np.linalg.norm(x[k] - x[j]) for j in sel)))
        kept += sel
    back = [k for k in on_b if k not in kept]
    s2 = s.copy()
    s2[back] = "Se"
    idx = [i for i in range(len(s)) if i not in gone]
    s3 = np.concatenate([s2[idx], ["Cl"] * len(bridges)])
    x3 = np.vstack([x[idx]] + [p for _e, p in bridges]) if bridges else x[idx]
    src = idx + [-1] * len(bridges)
    n = Counter(s3.tolist())
    q = 2 * n["Cd"] - 2 * n["Se"] - n["Cl"]
    log(f"half swap: {len(on_a)} Cl off {{111}}A, {len(bridges)} bridges, {len(kept)} of {len(on_b)} {{111}}B swaps kept, "
        f"{len(back)} back to Se; {len(s3)} atoms {dict(n)}, charge {q:+d}")
    if q != 0:
        raise ValueError(f"half swap: charge {q:+d}")
    return s3, x3, src


def build_size(r: dict, entry: dict, log=print) -> dict:
    """One dot of the recipe: symbols, positions (centred on the mean), role per atom, Cd per role."""
    size, surface = float(entry["size"]), entry["surface"]
    rec = builder_recipe(r, size, surface)
    pq = r["builder"].get("positive_q_mode", "add")
    s, x, man = run_builder(r["material"]["cif"], rec, pq)
    roles = list(man["site_roles"]["roles"])
    if r["builder"].get("rebuild_from_distances", True):
        # the dot again from its construction planes (facet_distances, as the 2026-10 layer dots were made): the same
        # atoms in the order those runs had, which fixes the tie-breaks of the half swap; roles mapped by position
        d = {}
        for pl in man["facet_geometry"]["construction"]["planes"]:
            d.setdefault(pl["label"], []).append(float(pl["d_A"]))
        rec = {k: v for k, v in rec.items() if k not in ("size_unit_cells", "facets")}
        rec["facet_distances"] = {k: float(round(float(np.mean(v)), 4)) for k, v in d.items()}
        s2, x2, man = run_builder(r["material"]["cif"], rec, pq)
        dist, j = cKDTree(x).query(x2)
        if len(s2) != len(s) or dist.max() > 1e-3 or any(s[j] != s2):
            raise RuntimeError(f"rebuild from facet distances differs from the size build (max {dist.max():.2e} A)")
        roles = [roles[int(k)] for k in j]
        s, x = s2, x2
    N, lab = facet_planes(man, s, x, s != "Cl") if surface == "half_swap" else (None, None)
    # centred on the as-built dot and rounded to 1e-5 A before the half swap, as the 2026-10 runs stored it: the
    # farthest-point choice of the kept {111}B swaps breaks ties on these coordinates
    x = np.round(x - x.mean(0), 5)
    if surface == "half_swap":
        s, x_new, src = half_swap(s, x, N, lab, log)
        roles_b = roles
        roles = [roles_b[i] if i >= 0 else None for i in src]
        cd = np.flatnonzero(s == "Cd")
        tcd = cKDTree(x_new[cd])
        for i, j in enumerate(src):
            if j < 0:                                    # bridging Cl: the role of the nearest Cd
                roles[i] = roles[int(cd[tcd.query(x_new[i])[1]])]
        x = x_new
    elif surface not in ("plain", "swap"):
        raise ValueError(f"unknown surface treatment {surface!r}")
    counts = Counter(roles[i] for i in np.flatnonzero(s == "Cd"))
    comp = Counter(s.tolist())
    return {"tag": rc.tag(size), "size": size, "surface": surface, "builder_recipe": rec,
            "symbols": s.tolist(), "positions": np.round(x, 5).tolist(), "roles": roles,
            "composition": dict(comp), "cd_per_role": dict(sorted(counts.items())),
            "site_roles_source": "QD_builder manifest site_roles.roles on the as-built dot, carried by atom"}


def run(r: dict, out: Path, sizes=None, log=print, recipe_path: Path | None = None) -> list[Path]:
    """Every size of the recipe (or `sizes`); the recipe file is copied to <run>/recipe.yaml, which every later step
    reads (an existing different one is an error: one run, one recipe)."""
    run_dir = Path(out)
    run_dir.mkdir(parents=True, exist_ok=True)
    own = run_dir / rc.RUN_RECIPE
    text = Path(recipe_path or rc.template_path()).read_text()
    if own.exists() and own.read_text() != text:
        raise ValueError(f"{own} exists and differs from {recipe_path or 'the template'}: one run, one recipe")
    own.write_text(text)
    out = run_dir / "dots"
    out.mkdir(parents=True, exist_ok=True)
    paths = []
    for e in r["sizes"]:
        if sizes and rc.tag(e["size"]) not in sizes:
            continue
        d = build_size(r, e, log)
        p = out / f"{d['tag']}.json"
        p.write_text(json.dumps(d))
        log(f"{d['tag']} ({d['surface']}): {len(d['symbols'])} atoms {d['composition']}; Cd per role "
            + ", ".join(f"{k} {v}" for k, v in d["cd_per_role"].items()))
        paths.append(p)
    return paths
