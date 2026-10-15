"""
Collect a reaction run into results/results.json: per species E^DFT, the thermal terms on the
T grid (thermo.components from the MACE vibrations), GB on the Hirshfeld charges, the SCCS
knots (electrostatics, surface, volume), the ALPB check; flags. Then plots and reaction.html.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from ..xyz import read_xyz
from . import solvation, thermo
from .cp2k import eps_dir

HARTREE_EV = 27.211386245988


def _hirshfeld(sdir: Path):
    if (sdir / "charges.txt").is_file():
        rows = [l.split() for l in (sdir / "charges.txt").read_text().splitlines() if l.strip() and not l.startswith("#")]
        return np.array([float(r[-1]) for r in rows]), "user"
    d = sdir / "cp2k" / "charges"
    if not (d / ".done").is_file():
        return None, None
    from orchestr_ai.qd.field import read_hirshfeld
    return read_hirshfeld(d), "CP2K PBE (GAPW) Hirshfeld"


def _dipole(out: Path):
    """|dipole| (Debye) printed by &MOMENTS."""
    if not out.is_file():
        return None
    import re
    m = re.findall(r"Dipole moment \[Debye\]\s*\n.*Total=\s*(\S+)", out.read_text(errors="ignore"))
    return float(m[-1]) if m else None


def _dispersion(out: Path):
    """D3 dispersion energy (eV) printed in the energy summary of a CP2K run."""
    import re
    m = re.findall(r"Dispersion energy:\s+(-?\d+\.\d+)", out.read_text(errors="ignore"))
    return float(m[-1]) * HARTREE_EV if m else None


def _point_dipole(pos, q) -> float:
    x = np.asarray(pos, float)
    return float(np.linalg.norm((np.asarray(q)[:, None] * (x - x.mean(0))).sum(0)) * 4.80320)


def _cp2k_energy(out: Path):
    if not out.is_file():
        return None
    import re
    e = re.findall(r"ENERGY\| Total FORCE_EVAL \( QS \) energy[^\n-]*?(-?\d+\.\d+)\s*$", out.read_text(errors="ignore"), re.M)
    return float(e[-1]) * HARTREE_EV if e else None


def species_result(rx: dict, name: str, run: Path) -> dict:
    sp = dict(rx["species"][name])
    sdir = run / "species" / name
    flags = []
    sym, pos = read_xyz(sdir / "dft.xyz")
    mace = json.loads((sdir / "mace.json").read_text()) if (sdir / "mace.json").is_file() else None
    if mace is None:
        flags.append("no MACE vibrations (mace.json missing): G_corr unavailable")
        sp["thermo"] = {k: [float("nan")] * len(rx["T"]) for k in
                        ("zpe", "dU_vib", "H_trans", "H_rot", "S_trans", "S_rot", "S_vib", "S_el", "H_corr", "S", "G_corr")}
    else:
        flags += mace["flags"]
        sp["thermo"] = thermo.components(mace["symbols"], mace["positions_A"], mace["masses_amu"],
                                         mace["frequencies_cm1"], rx["T"], mace["symmetry_number"], sp["mult"])
        sp["energy_mace_min_eV"] = mace["energy_mace_min_eV"]
        sp["mace"] = {k: mace[k] for k in ("rmsd_A", "rmsd_heavy_A", "relaxation_energy_eV", "fmax_dft_geometry_eV_A",
                                           "n_imaginary", "lowest_cm1", "symmetry_number", "hessian", "relax")}
        sp["frequencies_cm1"] = mace["frequencies_cm1"]
    # charges for the charge-based solvation: RESP (fitted to the DFT potential, cp2k/esp) where
    # available; the Hirshfeld charges (CP2K default, Gaussian shape functions) only for comparison
    q_h, src = _hirshfeld(sdir)
    q_r = solvation.read_resp(sdir / "cp2k" / "esp" / "resp.dat")
    sp["gb_inf_eV"] = sp["cosmo_inf_eV"] = None
    dip = {"dft": _dipole(sdir / "cp2k" / "charges" / "cp2k.out")}
    if q_h is not None and len(q_h) == len(sym):
        sp["gb_hirshfeld_inf_eV"] = solvation.gb_inf(sym, pos, q_h)
        dip["hirshfeld"] = _point_dipole(pos, q_h)
    if q_r is not None and len(q_r) == len(sym):
        sp["charges_source"] = "RESP (CP2K, fitted to the potential of a gas-phase single point)"
        sp["gb_inf_eV"] = solvation.gb_inf(sym, pos, q_r)
        sp["cosmo_resp_inf_eV"] = solvation.cosmo_inf(sym, pos, q_r, solvation.COSMO_SCALE)["E_inf_eV"]
        dip["resp"] = _point_dipole(pos, q_r)
    sp["dipole_D"] = dip
    cubes = solvation.esp_cubes(sdir / "cp2k" / "esp")
    if cubes:
        # COSMO on the DFT potential, cavity = the density isosurface (the default; see solvation.py)
        iso = solvation.cosmo_isodensity(solvation.read_cube(cubes[0]), solvation.read_cube(cubes[1]))
        sp["cosmo_inf_eV"] = iso["E_inf_eV"]
        sp["cavity"] = {"S_A2": iso["S_A2"], "V_A3": iso["V_A3"], "source": f"density isosurface {iso['rho0']:.2e} e/bohr^3",
                        "n_points": iso["n_points"]}
    else:
        flags.append("no potential/density cubes (cp2k/esp): COSMO falls back to RESP charges in 1.2 x vdW spheres")
        sp["cosmo_inf_eV"] = sp.get("cosmo_resp_inf_eV")
        _, area, vol = solvation.cavity(sym, pos, solvation.CAVITY_SCALE)
        sp["cavity"] = {"S_A2": float(area.sum()), "V_A3": vol, "source": f"{solvation.CAVITY_SCALE} x vdW spheres"}
    if sp["cosmo_inf_eV"] is None:
        flags.append("no RESP charges or cubes (cp2k/esp): COSMO electrostatics unavailable")
    energies = {}
    sp["dispersion_eV"] = {}
    for key, d in (("pbe0_given", "energy_pbe0"), ("pbe0_mace", "energy_pbe0_mace"),
                   ("pbe0d3_given", "energy_pbe0_d3"), ("pbe0d3_mace", "energy_pbe0_mace_d3")):
        r = solvation.read_sccs(sdir / "cp2k" / d / "cp2k.out")
        if r:
            energies[key] = r["E_eV"]
            e_d = _dispersion(sdir / "cp2k" / d / "cp2k.out")
            if e_d is not None:
                sp["dispersion_eV"][key] = e_d
    if mace is not None:
        energies["mace"] = mace["energy_mace_min_eV"]
    sp["energies"] = energies
    e_cp2k = _cp2k_energy(sdir / "cp2k" / "charges" / "cp2k.out")
    sp_f = rx["settings"].get("sp_functional", "PBE")
    if e_cp2k is not None:
        sp["energy_cp2k_sp_eV"] = e_cp2k
        # same functional as the given energies: they must agree (to the settings); otherwise information only
        if rx.get("energy_method", "PBE").upper() == sp_f.upper() and abs(e_cp2k - sp["energy_eV"]) > 0.01:
            flags.append(f"CP2K {sp_f} single point differs from the given energy by {1000 * (e_cp2k - sp['energy_eV']):+.0f} meV "
                         "(different settings, or not the same structure)")
    vac = solvation.read_sccs(sdir / "cp2k" / "sccs_vac" / "cp2k.out")
    ef = rx["settings"].get("energy_check_functional")
    if ef:
        chk = solvation.read_sccs(sdir / "cp2k" / f"energy_{ef.lower()}" / "cp2k.out")
        if chk:
            sp["energy_cp2k_check_eV"] = chk["E_eV"]
    points = {}
    for e in rx["settings"]["sccs_eps"]:
        r = solvation.read_sccs(sdir / "cp2k" / eps_dir(e) / "cp2k.out")
        if vac and r and "S_A2" in r:
            points[str(float(e))] = {"dG_el_eV": r["E_eV"] - vac["E_eV"], "S_A2": r["S_A2"], "V_A3": r["V_A3"],
                                     "E_pol_eV": r.get("E_pol_eV")}
    if points:
        sp["sccs"] = {**solvation.sccs_knots(points), "points": points}
        if len(points) < len(rx["settings"]["sccs_eps"]):
            flags.append(f"SCCS: {len(points)} of {len(rx['settings']['sccs_eps'])} eps finished")
    elif rx["settings"]["sccs"] and sp["n_atoms"] <= rx["settings"]["sccs_max_atoms"]:
        flags.append("no SCCS results: SCCS electrostatics and the cavity/dispersion term unavailable")
    if (sdir / "alpb.json").is_file():
        sp["alpb"] = json.loads((sdir / "alpb.json").read_text())
    sp["flags"] = flags
    for k in ("symbols", "xyz", "charges_file"):
        sp.pop(k, None)
    return sp


def collect(run_dir: str) -> dict:
    run = Path(run_dir)
    rx = json.loads((run / "reaction.json").read_text())
    res = {k: rx.get(k) for k in ("name", "equation", "nu", "T", "conditions", "scan", "per", "settings", "energy_method",
                                  "ligand", "equilibria", "recipe")}
    res["equilibria"] = res["equilibria"] or []
    res["species"] = {n: species_result(rx, n, run) for n in rx["species"]}
    # the electronic reaction energy from the other sources at hand, to show how much rides on it
    st = rx["settings"]
    alt = {f"given ({rx.get('energy_method', 'PBE')})": "energy_eV", "MACE-MH-1 at the MACE minima": "energy_mace_min_eV",
           f"CP2K {st.get('sp_functional', 'PBE')} (GAPW) at the given geometries": "energy_cp2k_sp_eV"}
    # the electronic-energy sets the page can switch between: complete ones only
    method = rx.get("energy_method", "PBE")
    labels = {"pbe0_given": f"PBE0@{method} geometries", "mace": "MACE-MH1",
              "pbe0_mace": "PBE0@MACE-MH1 geometries", "pbe0d3_given": f"PBE0-D3(BJ)@{method} geometries",
              "pbe0d3_mace": "PBE0-D3(BJ)@MACE-MH1 geometries"}
    res["energy_sets"] = {"given": method}
    for k, lab in labels.items():
        if all(k in res["species"][n].get("energies", {}) for n in res["species"]):
            res["energy_sets"][k] = lab
    res["dE_alternatives_eV"] = {
        label: (sum(v * res["species"][n][key] for n, v in rx["nu"].items())
                if all(res["species"][n].get(key) is not None for n in rx["nu"]) else None)
        for label, key in alt.items()}
    res["nonel_model"] = {"alpha_gamma_mN_m": [solvation.HILLE_SLOPE, solvation.HILLE_OFFSET],
                          "beta_GPa": solvation.HILLE_BETA_GPA, "min_eps": solvation.NONEL_MIN_EPS,
                          "beta0_alpha_gamma_mN_m": [solvation.BETA0_SLOPE, solvation.BETA0_OFFSET],
                          "cavity_scale": solvation.CAVITY_SCALE, "cosmo_scale": solvation.COSMO_SCALE,
                          "rho_min": solvation.RHO_MIN, "rho_max": solvation.RHO_MAX,
                          "ref": "Hille et al., J. Chem. Phys. 150, 041710 (2019)"}
    return res


def analyse(run_dir: str) -> dict:
    from . import report
    run = Path(run_dir)
    res = collect(run_dir)
    out = run / "results"
    out.mkdir(exist_ok=True)
    (out / "results.json").write_text(json.dumps(res, indent=1, default=float))
    summary = report.summary(res)
    (out / "summary.json").write_text(json.dumps(summary, indent=1))
    report.plots(res, out)
    (out / "reaction.html").write_text(report.page(res))
    return summary
