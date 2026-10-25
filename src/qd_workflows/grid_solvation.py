"""
Implicit solvation of the nucleation grid and of its precursor pool (the campaign's model,
qd_workflows.reaction.solvation): every solution species, monomer or cluster, gets

    dG_solv(eps) = (1 - 1/eps) E_inf + (alpha + gamma)(eps) S,     (alpha + gamma) = 2.16 ln(eps) - 14.79 mN/m

    monomers (CdSe, CdCl2)  E_inf (COSMO on the CP2K potential, eps -> inf) and S (density isosurface) of the
                            campaign's monomer library
    clusters (m, k)         S from the geometry (vdW spheres x 1.4, which gives the isodensity surface of both
                            monomers within 5 %), E_inf = e_el S with e_el the campaign dots' electrostatic energy
                            per area (196 CdSe dot variants: -1.06 +- 0.23 meV/A^2, no trend with size)

eps = 1 is the gas phase (no solvation); the cavity term is defined from eps = 1.8 (pentane) up.

IMPLICIT SOLVENT ONLY. Explicit solvent and ligands (oleic acid, amines, phosphines) binding to the monomers, the
clusters and the facets are not included; they are the next step of the workflow and change every (m, k).

    python -m qd_workflows.grid_solvation surfaces <thermo.json> [--scale 1.4]    (adds S_vdW_A2 to every state)
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

SOLVATION = {
    "eps": 2.0,                                   # solvent dielectric constant (about octadecene); 1 = gas phase
    "vdw_scale": 1.4,                             # cluster surface: vdW spheres x this (isodensity within 5 %)
    "e_el_eV_A2": -1.06e-3,                       # campaign CdSe dots: E_inf / S (COSMO, eps -> inf)
    "monomers": {                                 # campaign monomer library (CP2K PBE, COSMO on the isodensity)
        "CdSe": {"E_inf_eV": -0.327178, "S_A2": 117.799},
        "CdCl2": {"E_inf_eV": -0.198261, "S_A2": 135.843},
    },
    "complexes": {                                # the grid's m = 1 row, CdSe(CdCl2)_k: the same CP2K ESP + COSMO
        "1": {"E_inf_eV": -0.279721, "S_A2": 208.647},   # (NHR job 9619590, MACE geometries of the grid)
        "2": {"E_inf_eV": -0.461963, "S_A2": 287.318},
        "3": {"E_inf_eV": -0.511390, "S_A2": 351.680},
    },
    "note": "implicit solvent only; explicit solvent and ligands (oleic acid, amines, phosphines) come later",
}
EPS_MIN = 1.8                                     # the cavity term's lower bound (pentane 1.84)


def surface(symbols, positions, scale: float = SOLVATION["vdw_scale"], density: float = 4.0) -> float:
    """Solvent-exposed area (A^2) of the union of the scaled vdW spheres (points on each sphere outside all others)."""
    from .reaction.solvation import VDW, _sphere_points
    x = np.asarray(positions, float)
    r = np.array([VDW.get(s, 2.0) for s in symbols]) * scale
    total = 0.0
    for a in range(len(x)):
        n = max(12, int(4 * math.pi * r[a] ** 2 * density))
        p = x[a] + r[a] * _sphere_points(n)
        d = np.linalg.norm(p[:, None, :] - x[None, :, :], axis=2)
        d[:, a] = np.inf
        total += (d > r[None, :]).all(axis=1).sum() * 4 * math.pi * r[a] ** 2 / n
    return float(total)


def per_area(eps: float, P: dict = SOLVATION) -> float:
    """Solvation free energy per area of a neutral nanocrystal surface (eV/A^2): (1 - 1/eps) e_el + (alpha + gamma)."""
    from .reaction.solvation import cavity_term
    if eps <= 1.0:
        return 0.0
    return (1.0 - 1.0 / eps) * float(P["e_el_eV_A2"]) + cavity_term(max(eps, EPS_MIN), 1.0, 0.0)


def cosmo(entry: dict, eps: float) -> float:
    """(1 - 1/eps) E_inf + (alpha + gamma) S of a species with its own COSMO data (E_inf_eV, S_A2)."""
    from .reaction.solvation import cavity_term
    if eps <= 1.0:
        return 0.0
    return (1.0 - 1.0 / eps) * float(entry["E_inf_eV"]) + cavity_term(max(eps, EPS_MIN), float(entry["S_A2"]), 0.0)


def monomer(name: str, eps: float, P: dict = SOLVATION) -> float:
    """Solvation free energy (eV) of a library monomer at eps."""
    return cosmo(P["monomers"][name], eps)


def state_shift(node: dict, state: dict, eps: float, P: dict = SOLVATION, monomer_names=("CdSe", "CdCl2")):
    """dG_solv of one grid state: the library monomer for the bare MA monomer (m = 1, k = 0), the precursor
    complex's own COSMO data for the rest of the m = 1 row when given (complexes: {k: ...}), else the cluster
    model on the state's surface."""
    if eps <= 1.0:
        return 0.0
    if node["m"] == 1 and node["k"] == 0 and monomer_names[0] in P["monomers"]:
        return monomer(monomer_names[0], eps, P)
    if node["m"] == 1 and str(node["k"]) in P.get("complexes", {}):
        return cosmo(P["complexes"][str(node["k"])], eps)
    if "S_vdW_A2" not in state:
        raise ValueError(f"state {state.get('file')} has no surface: run python -m qd_workflows.grid_solvation surfaces")
    return per_area(eps, P) * float(state["S_vdW_A2"])


def add_surfaces(thermo: dict, base: Path, scale: float = SOLVATION["vdw_scale"], log=print) -> int:
    """S_vdW_A2 of every state from its structure file (relative to the thermo's folder); returns states done."""
    from .xyz import read_xyz
    done = 0
    for n in thermo["nodes"]:
        for s in n["states"]:
            if "S_vdW_A2" in s:
                continue
            sym, pts = read_xyz(base / s["file"])
            s["S_vdW_A2"] = round(surface(list(sym), pts, scale), 2)
            done += 1
    thermo.setdefault("solvation", {})["surface"] = f"vdW spheres x {scale} (qd_workflows.grid_solvation)"
    log(f"[grid-solvation] surfaces of {done} states")
    return done


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python -m qd_workflows.grid_solvation")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("surfaces", help="add the geometric surface S_vdW_A2 to every state of a thermo file")
    p.add_argument("thermo")
    p.add_argument("--scale", type=float, default=SOLVATION["vdw_scale"])
    a = ap.parse_args(argv)
    path = Path(a.thermo)
    th = json.loads(path.read_text())
    add_surfaces(th, path.resolve().parent, a.scale)
    path.write_text(json.dumps(th))
    return 0


if __name__ == "__main__":
    sys.exit(main())
