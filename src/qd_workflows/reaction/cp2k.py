"""
CP2K part of a reaction run, per species (species/<name>/cp2k/):

    charges/      GAPW PBE single point at the DFT geometry, the settings of every other single
                  point (qdw cp2k prepare --step charges): Hirshfeld charges for GB
    sccs_vac/     GPW, no solvent: reference of the SCCS runs (restarts from charges/)
    sccs_e<eps>/  GPW SCCS at each eps of settings.sccs_eps (restart from sccs_vac/)

OT/DIIS, MAX_SCF 200, no outer loop; a run whose SCF has not converged is started again from
the .wfn it wrote (job_cp2k.sh, up to 4 times). Species run side by side on the job's cores.
"""
from __future__ import annotations

import math
import re
from collections import Counter
from pathlib import Path

import numpy as np

from .. import TEMPLATES, render
from ..cluster import cp2k_resources, sbatch_account
from ..cp2k import inputs
from ..cp2k.basis import kind_blocks, kinds_for
from ..xyz import read_xyz, write_xyz
from .solvation import RHO_MAX, RHO_MIN

SMALL_ATOMS = 20          # molecules this small run on SMALL_RANKS ranks
SMALL_RANKS = 2          # exact exchange on a few atoms: 8-24 ranks were 10-30x slower per SCF step
SCCS_VACUUM_A = 16.0      # added to the extent along each axis (periodic cell; neutral solutes)


ADMM_BASIS = "BASIS_ADMM_UZH"
ADMM_FAMILY = "admm-dzp"
# ADMM2 as in the CdSe PBE0 benchmark (CP2K/CdSe_dots/pbe0: 357 atoms, 48 ranks x 4 threads,
# 13.6 s per OT step once the integrals are stored)
ADMM = """    &AUXILIARY_DENSITY_MATRIX_METHOD
      ADMM_TYPE ADMM2
    &END
"""
HFX = """      &HF
        FRACTION 0.25
        &SCREENING
          EPS_SCHWARZ 1.0E-6
          SCREEN_ON_INITIAL_P TRUE
        &END SCREENING
        &MEMORY
          MAX_MEMORY {mem}          # MB per rank: the integrals are kept after the first step
          EPS_STORAGE_SCALING 0.1
        &END MEMORY
        &INTERACTION_POTENTIAL
{potential}
        &END INTERACTION_POTENTIAL
      &END HF"""
# isolated (non-periodic cell): the full Coulomb operator; periodic cell (SCCS): truncated at
# CUTOFF_RADIUS, below half the cell, so exchange sees no periodic image
COULOMB = "          POTENTIAL_TYPE COULOMB"
TRUNCATED = """          POTENTIAL_TYPE TRUNCATED
          CUTOFF_RADIUS {rc:.2f}
          T_C_G_DATA {tcg}"""


def xc_blocks(functional: str, *, periodic: bool, cell_min: float = 0.0, t_c_g: str = "",
              hfx_memory_mb: int = 16000) -> dict:
    """XC, ADMM, extra basis file and preconditioner of a GPW input for PBE, HLE17 or PBE0."""
    if functional != "PBE0":
        return {"XC": inputs.XC[functional], "ADMM": "", "AUX": None, "PRECONDITIONER": "FULL_KINETIC"}
    if periodic:
        if not t_c_g:
            raise KeyError("PBE0 in a periodic cell needs the truncated-Coulomb table: cluster cp2k.t_c_g")
        pot = TRUNCATED.format(rc=min(10.0, 0.5 * cell_min - 0.5), tcg=t_c_g)
    else:
        pot = COULOMB
    xc = ("    &XC\n      &XC_FUNCTIONAL\n        &PBE\n          SCALE_X 0.75\n          SCALE_C 1.0\n"
          "        &END PBE\n      &END XC_FUNCTIONAL\n" + HFX.format(mem=hfx_memory_mb, potential=pot) + "\n    &END XC")
    return {"XC": xc, "ADMM": ADMM, "AUX": ADMM_FAMILY, "PRECONDITIONER": "FULL_SINGLE_INVERSE"}


def sccs_cell(coords) -> str:
    """Periodic cell of the SCCS runs: the extent plus SCCS_VACUUM_A, at least 18 A."""
    ext = np.ptp(np.asarray(coords, float), axis=0) + SCCS_VACUUM_A
    return " ".join(f"{max(x, 18.0):.3f}" for x in ext)


def sccs_input(xyz: Path, out: Path, project: str, eps: float | None, cluster: dict, charge: int = 0,
               data_dir: str | None = None, functional: str = "PBE", periodic: bool = True,
               hfx_memory_mb: int = 16000) -> int:
    """
    <out>/{cp2k.inp, geom.xyz} for one GPW single point; returns n_bf. periodic: the SCCS setting
    (eps, or eps None for its vacuum reference); otherwise an isolated vacuum single point
    (MULTIPOLE, non-periodic cell, as the other single points), e.g. the PBE0 energy check.
    """
    cp = cluster["cp2k"]
    symbols, coords = read_xyz(xyz)
    local = inputs.local_data_dir(cluster, data_dir)
    kinds = kinds_for(symbols, local / cp["basis_file"], local / cp["potential_file"], cp["basis"], cp["potential"])
    cell = sccs_cell(coords) if periodic else inputs.cell_abc(coords)
    xc = xc_blocks(functional, periodic=periodic, cell_min=min(float(x) for x in cell.split()),
                   t_c_g=cp.get("t_c_g", ""), hfx_memory_mb=hfx_memory_mb)
    values = {
        "CHARGE": charge, "BASIS_FILE": f"{cp['data_dir']}/{cp['basis_file']}",
        "EXTRA_BASIS": f"    BASIS_SET_FILE_NAME {cp['data_dir']}/{ADMM_BASIS}\n" if xc["AUX"] else "",
        "POTENTIAL_FILE": f"{cp['data_dir']}/{cp['potential_file']}", "COORD_FILE": "geom.xyz",
        "XC": xc["XC"], "ADMM": xc["ADMM"], "XC_NAME": functional,
        # SCCS changes the Hamiltonian a lot once the solvent is switched on: FULL_SINGLE_INVERSE,
        # built from the starting one, made PBE0 SCCS runs at eps 10-80 diverge; FULL_KINETIC converges
        "PRECONDITIONER": "FULL_KINETIC" if periodic else xc["PRECONDITIONER"],
        # With the solvent on, OT/DIIS sloshes the density of large polarisable solutes (InAs dot
        # at eps 10: the energy jumped by ~1 Ha between steps, gradient stuck at 1e-3, although each
        # SCCS cycle converged): conjugate gradients with a two-point line search for SCCS runs
        "OT_MINIMIZER": ("        MINIMIZER  CG\n        LINESEARCH  2PNT\n" if eps is not None
                         else "        N_DIIS  7\n        MINIMIZER  DIIS\n"),
        "POISSON": "PERIODIC" if periodic else "MULTIPOLE", "PERIODIC": "XYZ" if periodic else "None",
        "CELL": cell, "KINDS": kind_blocks(kinds, xc["AUX"]), "PROJECT": project,
        "WFN_FILE": f"{project}-RESTART.wfn", "SCCS": "OFF" if eps is None else "ON",
        "EPS": f"{eps if eps is not None else 1.0:g}", "RHO_MAX": f"{RHO_MAX:.6e}", "RHO_MIN": f"{RHO_MIN:.6e}",
    }
    out.mkdir(parents=True, exist_ok=True)
    write_xyz(out / "geom.xyz", symbols, coords, f"{project}")
    (out / "cp2k.inp").write_text(render((TEMPLATES / "cp2k" / "sccs.inp").read_text(), values))
    counts = Counter(symbols)
    return sum(kinds[e].functions * n for e, n in counts.items())


# ESP of a species for charge- or potential-based solvation in post-processing: an isolated GPW
# single point (MT Poisson solver, cell twice the density's extent), RESP charges fitted to the
# potential on shells around the vdW surface, and the Hartree potential as a cube.
ESP_PROPERTIES = """  &PROPERTIES
    &RESP
      INTEGER_TOTAL_CHARGE T
      RESTRAIN_HEAVIES_TO_ZERO F
      STRIDE 2
      &SPHERE_SAMPLING
        AUTO_VDW_RADII_TABLE UFF
        AUTO_RMIN_SCALE 1.4
        AUTO_RMAX_SCALE 2.0
      &END SPHERE_SAMPLING
      &PRINT
        &RESP_CHARGES_TO_FILE
          FILENAME =resp.dat
        &END
      &END PRINT
    &END RESP
  &END PROPERTIES
"""
V_HARTREE = """    &PRINT
      &V_HARTREE_CUBE
        STRIDE 2 2 2
        FILENAME =v_hartree
      &END V_HARTREE_CUBE
      &E_DENSITY_CUBE
        STRIDE 2 2 2
        FILENAME =density
      &END E_DENSITY_CUBE
    &END PRINT
"""


def esp_from_single_point(sp_dir: Path, out: Path) -> None:
    """
    <out>/cp2k.inp: the species' single point (charges/, GAPW, MULTIPOLE, non-periodic) restarted from
    its converged .wfn, with RESP charges and the Hartree-potential cube added: one SCF cycle. (The
    MULTIPOLE potential gives the same RESP fit as an isolated MT cell: HCl +-0.174 against +-0.189 e,
    relative RMS 0.44 against 0.39, without MT's doubled cell.)
    """
    t = (sp_dir / "cp2k.inp").read_text()
    t = t.replace("      &HIRSHFELD", V_HARTREE.split("\n", 1)[1].rsplit("    &END PRINT", 1)[0] + "      &HIRSHFELD", 1)
    t = t.replace("  &SUBSYS", ESP_PROPERTIES + "  &SUBSYS", 1)
    out.mkdir(parents=True, exist_ok=True)
    (out / "cp2k.inp").write_text(t)
    (out / "geom.xyz").write_text((sp_dir / "geom.xyz").read_text())


D3BJ = """      &VDW_POTENTIAL
        DISPERSION_FUNCTIONAL PAIR_POTENTIAL
        &PAIR_POTENTIAL
          TYPE DFTD3(BJ)
          PARAMETER_FILE_NAME {param}
          REFERENCE_FUNCTIONAL {functional}
          R_CUTOFF [angstrom] 40.0        # beyond a whole dot (the default is 10.6 A)
        &END PAIR_POTENTIAL
      &END VDW_POTENTIAL
"""


def with_d3(src: Path, out: Path, param_file: str, functional: str = "PBE0") -> None:
    """
    <out>/{cp2k.inp, geom.xyz, .wfn}: the single point in <src> with Grimme's D3(BJ) for
    `functional` added, restarted from its converged .wfn. D3 is a pairwise energy that does not
    change the density: the SCF converges in its first cycles, and CP2K prints the dispersion
    energy separately (energy_<f>_d3 next to energy_<f>).
    """
    t = (src / "cp2k.inp").read_text()
    # into &XC, before its own &END (not into &XC_FUNCTIONAL, whose &END starts with the same text)
    t, n = re.subn(r"^(    &END XC\s*)$", D3BJ.format(param=param_file, functional=functional).rstrip("\n") + r"\n\1",
                   t, count=1, flags=re.M)
    if n != 1:
        raise ValueError(f"{src}/cp2k.inp: no '    &END XC' line to add D3 before")
    out.mkdir(parents=True, exist_ok=True)
    (out / "cp2k.inp").write_text(t)
    (out / "geom.xyz").write_text((src / "geom.xyz").read_text())
    for w in src.glob("*-RESTART.wfn"):
        (out / w.name).write_bytes(w.read_bytes())


def eps_dir(eps: float) -> str:
    return f"sccs_e{eps:g}"


def energy_dir(functional: str) -> str:
    return f"energy_{functional.lower()}"


def prepare_species(sp: dict, sdir: Path, cluster: dict, settings: dict, data_dir: str | None = None) -> dict:
    """All CP2K inputs of one species; returns {n_bf, runs, small} (ranks: assign_ranks)."""
    if sp["mult"] != 1:
        raise ValueError(f"{sp['name']}: open-shell species (multiplicity {sp['mult']}) are not set up for CP2K yet")
    xyz = sdir / "dft.xyz"
    base = sdir / "cp2k"
    runs = []
    if not sp.get("charges_file"):
        info = inputs.prepare(str(xyz), "charges", f"{sp['name']}-charges", str(base / "charges"), cluster,
                              charge=sp["charge"], data_dir=data_dir, xc=settings["sp_functional"])
        (base / "charges" / "job.sh").unlink()                  # run by job_cp2k.sh
        n_bf = info["basis_functions"]
        runs.append("charges")
    mem = int(settings.get("hfx_memory_mb", 16000))
    ef = settings.get("energy_check_functional")
    if ef:        # isolated vacuum single point at another level for the reaction energy (e.g. PBE0)
        n_bf = sccs_input(xyz, base / energy_dir(ef), f"{sp['name']}-{ef.lower()}", None, cluster, sp["charge"],
                          data_dir, ef, periodic=False, hfx_memory_mb=mem)
        runs.append(energy_dir(ef))
    # SCCS runs in a periodic cell (the non-periodic solvers are wrong with SCCS): neutral species only
    sccs = settings["sccs"] and sp["n_atoms"] <= settings["sccs_max_atoms"] and sp["charge"] == 0
    if sccs:
        f = settings["sccs_functional"]
        n_bf = sccs_input(xyz, base / "sccs_vac", f"{sp['name']}-vac", None, cluster, sp["charge"], data_dir, f,
                          hfx_memory_mb=mem)
        runs.append("sccs_vac")
        # the eps runs go single-threaded (job_cp2k.sh): stored integrals per rank shrink accordingly
        mem_1 = max(1000, mem // int(settings.get("omp_threads", 4)))
        for e in settings["sccs_eps"]:
            sccs_input(xyz, base / eps_dir(e), f"{sp['name']}-e{e:g}", float(e), cluster, sp["charge"], data_dir, f,
                       hfx_memory_mb=mem_1)
            runs.append(eps_dir(e))
    if not runs:
        return {"n_bf": 0, "runs": [], "small": True}
    return {"n_bf": int(n_bf), "runs": runs, "small": sp["n_atoms"] <= SMALL_ATOMS}


def assign_ranks(plan: dict, cluster: dict) -> int:
    """Ranks per species (in plan) and of the job: small molecules SMALL_RANKS each, the rest of
    whole nodes shared evenly among the larger species (SCCS on their MT cells is grid-bound)."""
    per_node = int(cluster["cores_per_node"])
    small = [n for n, p in plan.items() if p["runs"] and p["small"]]
    big = [n for n, p in plan.items() if p["runs"] and not p["small"]]
    for n in small:
        plan[n]["ranks"] = SMALL_RANKS
    if not big:
        return max(SMALL_RANKS * len(small), SMALL_RANKS)
    need = sum(cp2k_resources(cluster, plan[n]["n_bf"]).cores for n in big) + SMALL_RANKS * len(small)
    total = per_node * max(1, math.ceil(need / per_node))
    for n in big:
        plan[n]["ranks"] = (total - SMALL_RANKS * len(small)) // len(big)
    return total


def launch_threads(cp: dict) -> str:
    """MPI launch of $1 ranks with $2 OpenMP threads each."""
    if cp.get("launcher") == "srun":
        return "srun --exact -n $1 -c $2 --cpu-bind=cores"
    return "mpirun -np $1"


def job(run: Path, cluster: dict, plan: dict, *, time: str | None = None, qos: str | None = None,
        test_queue: bool = False, threads: int = 4) -> str:
    """job_cp2k.sh: every species' runs, species side by side on the job's ranks. test_queue: the
    profile's cp2k.test_partition, 1 h pieces that resubmit themselves (then the analysis)."""
    cp = cluster["cp2k"]
    per_node = int(cluster["cores_per_node"])
    total = assign_ranks(plan, cluster)
    n_bf = max(p["n_bf"] for p in plan.values())
    # SCCS cost is set by the MT grid, not by the basis: the profile's chain time (the job skips
    # finished runs when resubmitted)
    wall = time or ("01:00:00" if test_queue else cp["time"].get("chain", "12:00:00"))
    res = cp2k_resources(cluster, n_bf)
    # whole nodes: all their memory (stored HFX integrals); otherwise per allocated CPU
    mem = "#SBATCH --mem=0" if total >= per_node else f"#SBATCH --mem-per-cpu={res.mem_per_cpu}"
    order = sorted((n for n in plan if plan[n]["runs"]), key=lambda n: -plan[n]["n_bf"])
    values = {
        "JOB_NAME": f"rx-cp2k-{run.name}"[:64], "TIME": wall, "NTASKS": max(1, total // threads),
        "NODES": -(-total // per_node), "MEM": mem,
        "QOS": qos or cluster.get("qos", "regular"), "SETUP": "\n".join(cluster.get("setup", [])),
        "CP2K_SETUP": inputs.cp2k_setup(cp), "LAUNCH_N": launch_threads(cp),
        "THREADS": threads,
        "PARTITION": f"#SBATCH -p {cp['test_partition']}\n" if test_queue else inputs.cp2k_partition(cp),
        "PIECES": 1 if test_queue else 0, "ACCOUNT": sbatch_account(cluster), "OMP_NUM_THREADS": threads,
        "EXECUTABLE": cp["executable"],
        "SPECIES": "\n".join(f"launch {n} {plan[n]['ranks']}" for n in order),
    }
    return render((TEMPLATES / "slurm" / "cp2k_reaction.sh").read_text(), values)
