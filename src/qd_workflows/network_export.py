"""
Export of the relaxed layer-path network (the tree registry) for the network growth model (network_growth.py).

Reads the tree registry and its products and writes one self-contained folder:

    network.json   states (registry id, formula, m = N_Se, k = N_Cd - N_Se, framework, imaginary modes), the steps
                   between them (a -> b with b the larger: dm >= 0, and dk > 0 when dm = 0; kind, facet, sheet size),
                   and the island table (per parent and facet family: added dm, dk and the row of the structure)
    arrays.npz     per structure (the states, then the islands and ladder rungs): G_gas on the grid's temperatures
                   (ideal gas at 1 M, no solvation; the rotational symmetry number included) and the vdW surface S
                   (A^2), so that G(T, eps) = interp(T, G_gas) + per_area(eps) S, as for the grid's clusters
    geometry.npz   symbols (codes into network.json's elements) and positions of every state (the structure panel)

Inputs (formats of the tree scripts): registry.json (dots: key -> id, framework "se"; lineage "0": the seed dot),
results_by_key.json (key -> symbols, positions, E_eV, frequencies_cm1, n_imaginary; a state relaxed without a
Hessian has F_vib_model_eV instead, the learned F_vib (fvib_model) on the thermo temperatures), g_cache.json
("T|key|E" -> G at T with the solvation at eps, which gives back the symmetry number), layer_steps.json (a, b, dSe,
dCd, dCl, dir, sheet, same_slab_next), ladders.json ("parent|family" -> rungs b, dm, dk), the island index
(states_index.json: parent, family, dSe, dCd, status, pure, key, name) and the relaxed island states (states.json).

    python -m qd_workflows.network_export --registry registry.json --results results_by_key.json --g-cache g_cache.json \
        --steps layer_steps.json --ladders ladders.json --islands states_index.json --island-states batch*/states.json \
        --thermo thermo_merged.json --params templates/growth/cdse_zb.yaml --out DIR [--nproc 5]
"""
from __future__ import annotations

import argparse
import collections
import json
import re
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import yaml

from . import grid_solvation as gs
from . import grid_thermo as gt

FACETS = ("100", "111A", "111B", "110")
KEY_ELEMENTS = ("Cd", "Se", "Cl")                    # the registry key's sites (Cl: the X site, whatever the ligand)
ELEMENTS = ("Cd", "Se", "Cl", "C", "O", "H")         # geometry.npz symbol codes (C, O, H: carboxylate ligands)
_ATOM = re.compile(r"'(Cd|Se|Cl)'")


def counts(key: str) -> dict:
    c = collections.Counter(_ATOM.findall(key))
    return {e: c.get(e, 0) for e in KEY_ELEMENTS}


def formula(symbols) -> str:
    """Formula of the relaxed structure (Cd and Se always, other elements when present)."""
    c = collections.Counter(symbols)
    return "".join(f"{e}{c[e]}" for e in ELEMENTS if c[e] or e in ("Cd", "Se"))


def _vib(s):
    """(frequencies, None) of a state with a Hessian, or ([], learned F_vib on the temperatures) of one without."""
    if s.get("frequencies_cm1") is not None:
        return s["frequencies_cm1"], None
    return [], s["F_vib_model_eV"]


def _thermal(job):
    """G_gas with symmetry number 1 on the temperatures, and the vdW surface; sigma when asked (no cache). Without
    frequencies: translation, rotation and pV from ideal_gas_g plus the learned F_vib."""
    symbols, pos, energy, (freqs, fvib), temps, scale, want_sigma = job
    from ase import Atoms
    from orchestr_ai.qd.references import ideal_gas_g, rotational_symmetry_number
    pos = np.asarray(pos, float)
    sigma = rotational_symmetry_number(symbols, pos) if want_sigma else 1
    g = np.array(ideal_gas_g(Atoms(list(symbols), positions=pos), energy, np.asarray(freqs, float), temps, sigma), float)
    if fvib is not None:
        g += np.asarray(fvib, float)
    return g, gs.surface(list(symbols), pos, scale), int(sigma)


def export(a) -> dict:
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    reg = json.loads(Path(a.registry).read_text())
    key_of = {v["id"]: k for k, v in reg["dots"].items()}
    fw_of = {v["id"]: v["se"] for v in reg["dots"].values()}
    seed = int(reg["lineage"]["0"]["dot"])
    print(f"[export] registry: {len(key_of)} dots, seed {seed}", flush=True)
    res = json.loads(Path(a.results).read_text())
    have = {i for i, k in key_of.items() if k in res}
    thermo = json.loads(Path(a.thermo).read_text())
    temps = thermo["temperatures"]
    P = yaml.safe_load(Path(a.params).read_text())["solvation"]

    # ---- steps between relaxed dots, deduplicated; continuations of a sheet already started
    rows = json.loads(Path(a.steps).read_text())
    cont = set()
    for r in rows:
        for c, _f, _n1, _n2 in r.get("same_slab_next") or []:
            cont.add((r["b"], c))
    steps = {}
    for r in rows:
        x, y = r["a"], r["b"]
        if x == y or x not in have or y not in have:
            continue
        steps.setdefault((x, y), r)
    # ---- the shell moves: same CdSe framework, neighbouring k (lowest-G isomer to lowest-G isomer), and each isomer to
    # the lowest one of its (framework, k)
    gc = json.loads(Path(a.g_cache).read_text())
    t_cache = float(next(iter(gc)).split("|", 1)[0])
    g_at = {}
    for ck, v in gc.items():
        _t, dk, e = ck.split("|", 2)
        if dk in res and abs(float(e) - res[dk]["E_eV"]) < 1e-9:
            g_at[reg["dots"][dk]["id"]] = v
    comp = {i: counts(key_of[i]) for i in have}
    mk = {i: (c["Se"], c["Cd"] - c["Se"]) for i, c in comp.items()}
    by_fw = collections.defaultdict(lambda: collections.defaultdict(list))
    for i in have:
        if i in g_at:
            by_fw[fw_of[i]][mk[i][1]].append(i)
    shell = set()
    for levels in by_fw.values():
        low = {k: min(v, key=g_at.get) for k, v in levels.items()}
        for k, v in levels.items():
            shell.update((low[k], i) for i in v if i != low[k])
        ks = sorted(levels)
        shell.update((low[k1], low[k2]) for k1, k2 in zip(ks, ks[1:]))
    # ---- the component of the seed
    adj = collections.defaultdict(set)
    for x, y in list(steps) + list(shell):
        adj[x].add(y); adj[y].add(x)
    comp_ids, todo = {seed}, [seed]
    while todo:
        for y in adj[todo.pop()]:
            if y not in comp_ids and y in g_at:
                comp_ids.add(y); todo.append(y)
    ids = sorted(comp_ids)
    if a.limit:                                      # smoke test: the seed's first dots only
        ids = sorted(ids, key=lambda i: (mk[i][0], i))[:a.limit]
    row = {i: j for j, i in enumerate(ids)}
    print(f"[export] {len(steps)} steps between relaxed dots, {len(shell)} shell links; seed component {len(ids)} dots",
          flush=True)

    def orient(x, y):
        (m1, k1), (m2, k2) = mk[x], mk[y]
        return (x, y) if (m2, k2) >= (m1, k1) else (y, x)

    st_out, seen = [], set()
    for (x, y), r in steps.items():
        if x not in row or y not in row:
            continue
        x, y = orient(x, y)
        if (x, y) in seen:
            continue
        seen.add((x, y))
        dm, dk = mk[y][0] - mk[x][0], mk[y][1] - mk[x][1]
        f = r.get("dir") if r.get("dir") in FACETS else None
        if dm > 0 and f is not None and (r["a"], r["b"]) not in cont:
            kind = "layer"
        elif dm > 0:
            kind = "attach"
        elif dk != 0:
            kind = "shell"
        else:
            kind = "iso"
        st_out.append({"a": row[x], "b": row[y], "dm": dm, "dk": dk, "kind": kind, "facet": f,
                       "sheet": int(r.get("sheet") or 0)})
    for x, y in shell:
        if x not in row or y not in row:
            continue
        x, y = orient(x, y)
        if (x, y) in seen:
            continue
        seen.add((x, y))
        dm, dk = mk[y][0] - mk[x][0], mk[y][1] - mk[x][1]
        st_out.append({"a": row[x], "b": row[y], "dm": dm, "dk": dk, "kind": "shell" if dk else "iso", "facet": None,
                       "sheet": 0})

    # ---- islands (the relaxed builder islands, as island_barriers.py selects them) and the ladders' one-sheet rungs
    istates = {}
    for f in a.island_states:
        istates.update(json.loads(Path(f).read_text())["states"])
    extra = []                                       # (parent id, family, dm, dk, source, structure)
    for r in json.loads(Path(a.islands).read_text()):
        if r.get("status") not in ("todo", "relaxed") or not r.get("pure") or not 1 <= r.get("dCd", 0) <= 10:
            continue
        s = istates.get(r.get("name")) if r["status"] == "todo" else res.get(r["key"])
        if s is None or r["parent"] not in row:
            continue
        extra.append((r["parent"], r["family"], r["dSe"], r["dCd"] - r["dSe"], "island", s))
    parents = {(p, f) for p, f, *_ in extra}
    lad = json.loads(Path(a.ladders).read_text())["ladders"]
    for p, f in sorted(parents):
        for v in lad.get(f"{p}|{f}", {}).values():
            if key_of.get(v["b"]) in res:
                extra.append((p, f, v["dm"], v["dk"], "rung", res[key_of[v["b"]]]))
    print(f"[export] islands: {sum(e[4] == 'island' for e in extra)} relaxed island states and "
          f"{sum(e[4] == 'rung' for e in extra)} ladder rungs on {len(parents)} (parent, facet) pairs", flush=True)

    # ---- thermal data: every state and every island structure
    if a.limit:
        extra = extra[:a.limit]
    jobs = [(res[key_of[i]]["symbols"], res[key_of[i]]["positions"], res[key_of[i]]["E_eV"],
             _vib(res[key_of[i]]), temps, P["vdw_scale"], False) for i in ids]
    jobs += [(s["symbols"], s["positions"], s["E_eV"], _vib(s), temps, P["vdw_scale"], True)
             for *_x, s in extra]
    if any(fv is not None and len(fv) != len(temps) for *_x, (_f, fv), _t, _sc, _w in jobs):
        raise SystemExit("F_vib_model_eV is not on the thermo temperatures")
    n_model = sum(fv is not None for *_x, (_f, fv), _t, _sc, _w in jobs)
    if n_model:
        print(f"[export] {n_model} of {len(jobs)} structures without a Hessian: learned F_vib", flush=True)
    G = np.zeros((len(jobs), len(temps)))
    S = np.zeros(len(jobs))
    sig = np.ones(len(jobs), int)
    with Pool(a.nproc) as pool:
        for j, (g, s_, sg) in enumerate(pool.imap(_thermal, jobs, chunksize=4)):
            G[j], S[j], sig[j] = g, s_, sg
            if j % 500 == 0:
                print(f"[export] thermal {j}/{len(jobs)}", flush=True)
    # symmetry numbers of the states from the cached G at t_cache: G_cache = G_gas(sigma=1) + kT ln sigma + per_area S
    kt = gt.KB_EV * t_cache
    pa = gs.per_area(float(P["eps"]), P)
    bad = 0
    for j, i in enumerate(ids):
        resid = g_at[i] - (gt._at(t_cache, temps, G[j]) + pa * S[j])
        s_ = max(1, int(round(np.exp(resid / kt))))
        if abs(kt * np.log(s_) - resid) > 2e-3:
            bad += 1
        sig[j] = s_
        G[j] += gt.KB_EV * np.asarray(temps) * np.log(s_)       # rotation: + kT ln sigma
    print(f"[export] symmetry numbers: {collections.Counter(sig[:len(ids)].tolist())}; "
          f"{bad} states whose cached G is off the recomputed one by > 2 meV", flush=True)

    # ---- write
    states = []
    for i in ids:
        s = res[key_of[i]]
        states.append({"id": i, "formula": formula(s["symbols"]), "m": mk[i][0], "k": mk[i][1],
                       "fw": fw_of[i], "n_imag": int(s.get("n_imaginary") or 0)})
        if s.get("frequencies_cm1") is None:                   # learned F_vib; out of the model's domain flagged
            states[-1]["fvib"] = "model_ood" if s.get("fvib_out_of_domain") else "model"
    fw_index = {}
    for s in states:
        s["fw"] = fw_index.setdefault(s["fw"], len(fw_index))
    islands = [{"parent": row[p], "family": f, "dm": dm, "dk": dk, "src": src, "row": len(ids) + j}
               for j, (p, f, dm, dk, src, _s) in enumerate(extra)]
    meta = {"source": {k: str(v) for k, v in vars(a).items() if k not in ("island_states",)},
            "temperatures": temps, "seed": row[seed], "seed_id": seed, "vdw_scale": P["vdw_scale"], "elements": ELEMENTS,
            "states": states, "steps": st_out, "islands": islands}
    (out / "network.json").write_text(json.dumps(meta, separators=(",", ":")))
    np.savez_compressed(out / "arrays.npz", G_gas=G, S=S, sigma=sig)
    sym, pos, off = [], [], [0]
    for i in ids:
        s = res[key_of[i]]
        sym += [ELEMENTS.index(x) for x in s["symbols"]]
        pos.append(np.asarray(s["positions"], np.float32))
        off.append(len(sym))
    np.savez_compressed(out / "geometry.npz", symbols=np.array(sym, np.int8), positions=np.concatenate(pos),
                        offsets=np.array(off, np.int64))
    kinds = collections.Counter(s["kind"] for s in st_out)
    print(f"[export] {out}: {len(states)} states (m {min(s['m'] for s in states)}-{max(s['m'] for s in states)}), "
          f"{len(st_out)} steps {dict(kinds)}, {len(islands)} island rows", flush=True)
    return meta


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python -m qd_workflows.network_export")
    ap.add_argument("--registry", required=True)
    ap.add_argument("--results", required=True)
    ap.add_argument("--g-cache", required=True)
    ap.add_argument("--steps", required=True)
    ap.add_argument("--ladders", required=True)
    ap.add_argument("--islands", required=True, help="states_index.json of the island builds")
    ap.add_argument("--island-states", nargs="+", required=True, help="the relaxed island batches' states.json")
    ap.add_argument("--thermo", required=True, help="the grid's thermo_merged.json (temperatures)")
    ap.add_argument("--params", required=True, help="growth template with the solvation block (eps, vdw_scale)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--nproc", type=int, default=5)
    ap.add_argument("--limit", type=int, default=0, help="smoke test: only this many states and island rows")
    export(ap.parse_args(argv))
    return 0


if __name__ == "__main__":
    sys.exit(main())
