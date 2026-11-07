"""
Relaxation of the registry dots on LUMI and their free energies (the scripts stage.py, relax_states.py, collect.py and
gcache.py).

- stage: registry dots without results into batches wave/batchJ/d<id>.xyz, balanced by atoms^1.5; staged.json keeps
  the ids already staged, so the next wave skips them;
- relax_states: on the cluster, every *.xyz of a batch relaxed (no anneal) and polished with the (m, k) grid's relaxer
  (Orchestr.AI RelaxCache, MACE-MH-1 omat_pbe float64) and its FD Hessian; <name>_relaxed.xyz and states.json;
- collect: the batches' states.json into results_by_key.json (keyed by the registry's dot key: d<id> names, or older
  runs identified by their builder geometry);
- gcache: G per relaxed registry dot at the template's T into g_cache.json ("T|key|E" -> G: ideal gas on the thermo
  temperatures + per_area(eps) * vdW surface).

Job script: templates/slurm/layer_path_relax_lumi.sh.
"""
from __future__ import annotations

import json
import random
import shutil
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import yaml

from . import lattice as L


def stage(net, results, wave, n_batches, staged="staged.json"):
    net = Path(net)
    reg = json.loads((net / "registry.json").read_text())
    done = json.loads(Path(results).read_text()) if Path(results).exists() else {}
    staged_p = Path(staged)
    old = set(json.loads(staged_p.read_text())) if staged_p.exists() else set()
    todo = []
    for dk, v in reg["dots"].items():
        if dk in done or v["id"] in old:
            continue
        p = net / "builds" / v["dir"] / "out.xyz"
        todo.append((int(p.read_text().split("\n", 1)[0]), v["id"], p))
    todo.sort(reverse=True)
    load = [0] * n_batches
    for n, i, p in todo:
        j = load.index(min(load))
        load[j] += n ** 1.5
        d = Path(wave) / f"batch{j}"
        d.mkdir(parents=True, exist_ok=True)
        shutil.copy(p, d / f"d{i:05d}.xyz")
    staged_p.write_text(json.dumps(sorted(old | {i for _n, i, _p in todo})))
    print(f"{len(todo)} dots staged in {n_batches} batches; atoms max {todo[0][0] if todo else 0}, "
          f"load {[round(x) for x in load]}")
    return todo


def relax_states(root, model, out, pattern="*.xyz"):
    """Relax (no anneal), polish and FD-Hessian every structure of a batch folder (on the cluster, with Orchestr.AI)."""
    from orchestr_ai.qd.run import Settings
    from orchestr_ai.qd.steps import detachment as dt
    root, out = Path(root), Path(out)
    st = Settings(model=model, head="omat_pbe", dtype="float64", hessian="fd", detach_quick=False, detach_screen=False,
                  detach_rigid_refresh=False, detach_frames=False)
    cache = dt.RelaxCache(SimpleNamespace(props=root, settings=st, record_dir=root))
    files = sorted(p for p in root.glob(pattern) if not p.stem.endswith("_relaxed"))
    structs = [(list(s), x) for s, x in map(L.read_xyz, files)]
    t0 = time.time()
    res = cache.relax_many(structs, st)
    print(f"relaxed {len(res)} in {time.time() - t0:.0f} s", flush=True)
    pol = cache.relax_many([(a.get_chemical_symbols(), a.get_positions()) for a, _e, _ok in res], st, polish=True)
    print(f"polished in {time.time() - t0:.0f} s", flush=True)
    rows = {}
    for p, (_sy0, x0), (a, e, ok) in zip(files, structs, pol):
        sy, x = a.get_chemical_symbols(), a.get_positions()
        f = np.asarray(cache.frequencies(sy, x, st))
        rows[p.stem] = {"symbols": sy, "positions": np.round(x, 6).tolist(), "start_positions": np.round(x0, 6).tolist(),
                        "E_eV": float(e), "converged": bool(ok), "frequencies_cm1": f.round(3).tolist(),
                        "n_imaginary": int((f < -10.0).sum()), "lowest_cm1": float(f.min())}
        L.write_xyz(root / f"{p.stem}_relaxed.xyz", sy, x, f"{p.stem} E={e:.6f} eV converged={ok}")
        print(f"{p.stem}: E = {e:.5f} eV  converged {ok}  imaginary {rows[p.stem]['n_imaginary']}  "
              f"lowest {f.min():.1f} cm-1  ({time.time() - t0:.0f} s)", flush=True)
        out.write_text(json.dumps({"settings": {"model": model, "head": st.head, "dtype": st.dtype, "hessian": st.hessian,
                                                "relaxer": st.relaxer(), "polish_fmax": dt.POLISH_FMAX},
                                   "states": rows}))
    cache.close()
    print("done", flush=True)
    return rows


def collect(net, results, states_files):
    """Merge relaxed batches into results_by_key.json."""
    reg = json.loads((Path(net) / "registry.json").read_text())
    by_id = {v["id"]: k for k, v in reg["dots"].items()}
    resf = Path(results)
    res = json.loads(resf.read_text()) if resf.exists() else {}
    n0 = len(res)
    for f in states_files:
        for name, s in json.loads(Path(f).read_text())["states"].items():
            if name.startswith("d") and name[1:].isdigit():
                k = by_id.get(int(name[1:]))
            else:   # earlier runs: identify the dot by its builder geometry
                k = repr(L.dot_key(np.array(s["symbols"]), np.array(s["start_positions"])))
            if k is not None and k in reg["dots"]:
                res[k] = {kk: s[kk] for kk in ("symbols", "positions", "start_positions", "E_eV", "converged",
                                                 "frequencies_cm1", "n_imaginary", "lowest_cm1")}
    resf.write_text(json.dumps(res))
    print(f"results: {n0} -> {len(res)} of {len(reg['dots'])} registry dots")
    return res


def g_model(thermo, params):
    """(T, g_of): G of a relaxed state (symbols, positions, E_eV, frequencies_cm1) at the template's T, with the
    solvation (per_area(eps) * vdW surface)."""
    from .. import grid_solvation as gs
    from .. import grid_thermo as gt
    temps = json.loads(Path(thermo).read_text())["temperatures"]
    tmpl = yaml.safe_load(Path(params).read_text())
    T = float(tmpl["conditions"]["T"])
    P = tmpl["solvation"]
    eps = float(P["eps"])

    def g_of(s):
        pos = np.array(s["positions"])
        _sig, g = gt._gas(s["symbols"], pos, s["E_eV"], s["frequencies_cm1"], temps)
        return gt._at(T, temps, g) + gs.per_area(eps, P) * gs.surface(list(s["symbols"]), pos, P["vdw_scale"])
    return T, g_of


def gcache(net, results, g_cache, thermo, params, check=0):
    """G per relaxed registry dot at the template's T into g_cache.json (check > 0: recompute that many and print)."""
    reg = json.loads((Path(net) / "registry.json").read_text())
    res = json.loads(Path(results).read_text())
    gcf = Path(g_cache)
    gc = json.loads(gcf.read_text()) if gcf.exists() else {}
    T, g_of = g_model(thermo, params)
    if check:
        for ck in random.Random(0).sample([k for k in gc if k.split("|", 2)[1] in res], check):
            print(f"cached {gc[ck]:.6f}  recomputed {g_of(res[ck.split('|', 2)[1]]):.6f}")
        return gc
    new = 0
    for dk in reg["dots"]:
        s = res.get(dk)
        if s is None:
            continue
        ck = f"{T}|{dk}|{s['E_eV']}"
        if ck not in gc:
            gc[ck] = g_of(s)
            new += 1
    gcf.write_text(json.dumps(gc))
    print(f"g_cache: {new} new, {len(gc)} in all")
    return gc
