"""
MACE-MH-1 part of a reaction run: for every species, relax from the DFT structure to the MACE
minimum, take the Hessian there and store the vibrations (species/<name>/mace.json).

The DFT structure is not a stationary point of MACE; a Hessian taken there has a residual
gradient and spurious soft or imaginary modes. The electronic energy stays the DFT one;
only G^corr comes from these vibrations. The heavy-atom RMSD of the relaxation is reported
and flagged above settings.rmsd_warn.

    python -m qd_workflows.reaction.mace <run dir> [--model PATH] [--device auto] [--dtype float64]
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from ..xyz import read_xyz, write_xyz


def kabsch_rmsd(a, b, mask=None) -> float:
    """RMSD after optimal superposition, over the atoms in `mask` (all by default)."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    if mask is not None:
        a, b = a[mask], b[mask]
    a = a - a.mean(axis=0)
    b = b - b.mean(axis=0)
    u, _, vt = np.linalg.svd(a.T @ b)
    d = np.sign(np.linalg.det(u @ vt))
    r = u @ np.diag([1.0, 1.0, d]) @ vt
    return float(np.sqrt(((a @ r - b) ** 2).sum(axis=1).mean()))


def relax(atoms, fmax: float, max_steps: int) -> dict:
    from ase.optimize import BFGS, FIRE, LBFGS
    first = BFGS if len(atoms) <= 500 else LBFGS
    opt = first(atoms, logfile=None)
    converged = bool(opt.run(fmax=fmax, steps=max_steps))
    steps, name = opt.get_number_of_steps(), first.__name__
    if not converged:
        opt = FIRE(atoms, logfile=None)
        converged = bool(opt.run(fmax=fmax, steps=max_steps))
        steps += opt.get_number_of_steps()
        name += "+FIRE"
    return {"converged": converged, "optimizer": name, "steps": int(steps),
            "fmax_eV_A": float(np.linalg.norm(atoms.get_forces(), axis=1).max())}


def run_species(sp: dict, out: Path, settings: dict, model: str, device: str, dtype: str) -> dict:
    from ase import Atoms
    from orchestr_ai.qd.engines import mace_calculator, mace_hessian, mace_provenance
    from orchestr_ai.qd.references import rotational_symmetry_number
    from orchestr_ai.qd.steps.hessian import _fd_hessian, vibrations

    t0 = time.time()
    symbols, pos0 = read_xyz(out / "dft.xyz")
    atoms = Atoms(list(symbols), positions=pos0)
    atoms.calc = mace_calculator(settings["head"], model, device, dtype)
    e0 = float(atoms.get_potential_energy())
    f0 = float(np.linalg.norm(atoms.get_forces(), axis=1).max()) if len(atoms) > 1 else 0.0
    rel = {"converged": True, "optimizer": "none", "steps": 0, "fmax_eV_A": 0.0}
    if len(atoms) > 1:
        rel = relax(atoms, settings["fmax"], settings["max_steps"])
    e1 = float(atoms.get_potential_energy())
    pos = atoms.get_positions()
    heavy = np.array([s != "H" for s in symbols])
    rmsd = kabsch_rmsd(pos0, pos)
    rmsd_heavy = kabsch_rmsd(pos0, pos, heavy) if heavy.sum() >= 2 else 0.0
    write_xyz(out / "mace_relaxed.xyz", symbols, pos, f"{sp['name']} MACE-MH-1/{settings['head']} minimum E={e1:.6f} eV")

    freqs, tr, method = np.zeros(0), np.zeros(0), "none"
    if len(atoms) > 1:
        n = len(atoms)
        method = settings["hessian"]
        if method == "auto":
            method = "analytic" if n <= settings["analytic_max_atoms"] else "fd"
        h = mace_hessian(atoms, atoms.calc) if method == "analytic" else _fd_hessian(atoms, settings["fd_step"], 4000)
        freqs, _, tr = vibrations(h, pos, atoms.get_masses())
    sigma = sp.get("symmetry_number") or (rotational_symmetry_number(symbols, pos) if len(atoms) > 1 else 1)
    flags = []
    if rmsd_heavy > settings["rmsd_warn"]:
        flags.append(f"heavy-atom RMSD DFT -> MACE minimum {rmsd_heavy:.3f} Å > {settings['rmsd_warn']} Å")
    if not rel["converged"]:
        flags.append(f"MACE relaxation not converged (fmax {rel['fmax_eV_A']:.3g} eV/Å)")
    if (freqs < -50.0).any():
        flags.append(f"{int((freqs < -50).sum())} imaginary modes below -50 cm-1 (lowest {freqs.min():.0f})")
    res = {
        "name": sp["name"], "n_atoms": len(atoms), "symbols": list(symbols),
        "positions_A": pos.tolist(), "masses_amu": atoms.get_masses().tolist(),
        "energy_mace_dft_geometry_eV": e0, "fmax_dft_geometry_eV_A": f0, "energy_mace_min_eV": e1,
        "relaxation_energy_eV": e1 - e0, "rmsd_A": rmsd, "rmsd_heavy_A": rmsd_heavy, "relax": rel,
        "hessian": method, "frequencies_cm1": freqs.tolist(),
        "n_imaginary": int((freqs < -10).sum()), "lowest_cm1": float(freqs.min()) if freqs.size else None,
        "tr_residual_max_cm1": float(np.abs(tr).max()) if tr.size else None,
        "symmetry_number": int(sigma), "flags": flags, "seconds": round(time.time() - t0, 1),
        "provenance": mace_provenance(settings["head"], model, device, dtype),
    }
    (out / "mace.json").write_text(json.dumps(res, indent=1))
    return res


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir")
    ap.add_argument("--model", default=None, help="MACE-MH-1 model (default: $QDPROPS_MACE_MODEL)")
    ap.add_argument("--device", default="auto")
    ap.add_argument("--dtype", default="float64")
    ap.add_argument("--species", help="comma-separated subset")
    ap.add_argument("--force", action="store_true", help="recompute species that have mace.json")
    a = ap.parse_args(argv)
    from orchestr_ai.qd.engines import DEFAULT_MODEL
    run = Path(a.run_dir)
    rx = json.loads((run / "reaction.json").read_text())
    names = a.species.split(",") if a.species else sorted(rx["species"], key=lambda n: rx["species"][n]["n_atoms"])
    failed = []
    for name in names:
        out = run / "species" / name
        if (out / "mace.json").is_file() and not a.force:
            print(f"[mace] {name}: done", flush=True)
            continue
        try:
            r = run_species(rx["species"][name], out, rx["settings"], a.model or DEFAULT_MODEL, a.device, a.dtype)
            print(f"[mace] {name}: {r['n_atoms']} atoms, RMSD {r['rmsd_heavy_A']:.3f} Å, "
                  f"{r['n_imaginary']} imaginary, lowest {r['lowest_cm1']}, {r['seconds']} s"
                  + "".join(f"\n[mace]   ! {f}" for f in r["flags"]), flush=True)
        except Exception as exc:                              # one species must not stop the others
            failed.append(name)
            print(f"[mace] {name}: FAILED {type(exc).__name__}: {exc}", flush=True)
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
