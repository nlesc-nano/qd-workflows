"""Step 4, states (GPU): every state of every site relaxed on its own, with a local free energy.

  1. free CdCl2: relaxed (float64), all 3N - 5 modes, translation at 1 M and rotation
  2. relax every state of every site in one batch (float32, L-BFGS to states.fmax), float64 energies
  3. partial Hessian of each state over the site's local set plus the state's own site atoms (+-step_A, float64);
     F_vib harmonic over the real modes, below floor_cm held at the floor; G = E64 + F_vib
  4. dG_bind,X = G(ref) - G(empty) - G(CdCl2, 1 M),  dG_M_minus_ref = G(M) - G(ref)
Check: imaginary modes (below states.imag_cm) per state, and how far the relaxation moved atoms.
Writes states/<tag>.json; resumable (relaxations and Hessians already done are kept)."""
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np

from . import engine as en
from . import recipe as rc


def _save(path: Path, S: dict):
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(S))
    os.replace(tmp, path)


def run_job(engine, job: dict, path: Path, r: dict, log=print) -> dict:
    T = float(r["T_K"])
    ph, st = r["phva"], r["states"]
    S = json.loads(path.read_text()) if path.exists() else {"dot": job["dot"], "reference": job["reference"], "sites": {}}
    if "CdCl2" not in S:
        cat, cl = r["material"]["cation"], r["ligand"]["anion"]
        S["CdCl2"] = en.free_molecule(engine, [cat, cl, cl], [[0, 0, 0], [2.27, 0, 0], [-2.27, 0, 0]], T,
                                      float(ph["step_A"]), float(ph["floor_cm"]))
        _save(path, S)
        log(f"free CdCl2: E {S['CdCl2']['E']:.4f}, G({T:.0f} K, 1 M) {S['CdCl2']['G']:.4f} eV")
    sites = job["sites"]
    todo = [(t, k) for t, s in sites.items() for k in s["states"] if k not in S["sites"].get(t, {}).get("relaxed", {})]
    if todo:
        items = [(sites[t]["states"][k]["symbols"], np.array(sites[t]["states"][k]["positions"], float)) for t, k in todo]
        log(f"relaxing {len(items)} structures ({sum(len(s) for s, _x in items)} atoms)")
        res = engine.relax_many(items, fmax=float(st["fmax"]), max_steps=int(st["max_steps"]))
        E = [float(e) for e, _f in engine.ev64([(s, np.asarray(q.positions, float)) for (s, _x), q in zip(items, res)])]
        for (t, k), (_s, x0), q, e in zip(todo, items, res, E):
            xr = np.asarray(q.positions, float)
            mv = np.linalg.norm(xr - x0, axis=1)
            S["sites"].setdefault(t, {"relaxed": {}})["relaxed"][k] = {
                "E": e, "steps": int(q.steps), "reason": q.reason, "max_move": float(mv.max()),
                "n_moved_1A": int((mv > 1).sum()), "positions": np.round(xr, 5).tolist()}
        _save(path, S)
        log("relaxed")
    for t, s in sites.items():
        R = S["sites"][t]["relaxed"]
        n0 = s["n_context"]
        local = list(s["local"])
        sets = {"ref": local + [n0, n0 + 1, n0 + 2], "empty": local, "M": local + [n0, n0 + 1]}
        V = S["sites"][t].setdefault("vib", {})
        for k, idx in sets.items():
            if k in V:
                continue
            fr = en.phva(engine, s["states"][k]["symbols"], np.array(R[k]["positions"], float), idx, float(ph["step_A"]))
            V[k] = {"F_vib": en.vib_free(fr, T, 0, float(ph["floor_cm"])), "n_imag": int(sum(f < float(st["imag_cm"]) for f in fr)),
                    "lowest": sorted(fr)[:6], "n_local": len(idx)}
            _save(path, S)
        g = {k: R[k]["E"] + V[k]["F_vib"] for k in sets}
        S["sites"][t].update({"role": s["role"], "env": s["env"], "cd": s["cd"], "G": g,
                              "dG_bind_X": g["ref"] - g["empty"] - S["CdCl2"]["G"], "dG_M_minus_ref": g["M"] - g["ref"]})
        _save(path, S)
        log(f"{t:44s} dG_bind(X) {S['sites'][t]['dG_bind_X']:+.3f}  G(M)-G(ref) {g['M'] - g['ref']:+.3f} eV  "
            f"imag {[V[k]['n_imag'] for k in sets]}  moved>1A {[R[k]['n_moved_1A'] for k in sets]}")
    S["done"] = True
    _save(path, S)
    return S


def run(r: dict, out: Path, engine, sizes=None, log=print) -> list[Path]:
    out = Path(out)
    (out / "states").mkdir(parents=True, exist_ok=True)
    paths = []
    for e in r["sizes"]:
        t = rc.tag(e["size"])
        if sizes and t not in sizes:
            continue
        job = json.loads((out / "sites" / f"{t}.json").read_text())
        p = out / "states" / f"{t}.json"
        run_job(engine, job, p, r, log)
        paths.append(p)
    return paths
