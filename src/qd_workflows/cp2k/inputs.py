"""CP2K inputs and job scripts for one structure and one step."""
from __future__ import annotations

import json
import os
import re
from collections import Counter
from pathlib import Path

import numpy as np

from .. import TEMPLATES, render
from ..cluster import cp2k_resources, node_tasks, sbatch_account
from ..xyz import read_xyz, write_xyz
from .basis import kind_blocks, kinds_for

STEPS = ("geo_opt", "pdos", "trexio", "sample", "charges")
VACUUM_A = 12.0          # added to the extent along each axis (6 Å on each side)
# auto = BFGS at every size (Ivan, 2026-10-22): above 1,000 atoms LBFGS stopped on its energy criterion
# (REL_REDUCTION_OF_F <= FACTR*EPSMCH) with the maximum force 1-15x above MAX_FORCE (15 dots on NHR, 2 on Hyperion).
LBFGS_ABOVE_ATOMS = None

PDOS_PRINT = """    &PRINT
      &DOS
        NLUMO -1
        &PDOS
        &END
      &END
    &END"""
NO_RESTART = """      &PRINT
        &RESTART OFF
        &END
      &END"""
FORCES_PRINT = """  &PRINT
    &FORCES ON
      FILENAME =forces.xyz
      NDIGITS 8
    &END
  &END"""
# charges: Hirshfeld charges (and the dipole) of a structure at its own geometry, from scratch
# (atomic guess unless its own .wfn exists). An SCF that has not converged in MAX_SCF cycles
# still writes its .wfn; the job then restarts from it (qdw reaction).
CHARGES_PRINT = """    &PRINT
      &MOMENTS
        PERIODIC F
      &END MOMENTS
      &HIRSHFELD
        FILENAME =hirshfeld.dat
      &END HIRSHFELD
    &END PRINT"""
IGNORE_SCF_FAILURE = """      IGNORE_CONVERGENCE_FAILURE
      &PRINT
        &RESTART              # every 10 cycles: a run cut off by the time limit continues from here
          BACKUP_COPIES 0
          &EACH
            QS_SCF 10
          &END
        &END
      &END"""
XC = {
    "PBE": """    &XC
      &XC_FUNCTIONAL PBE
      &END
    &END""",
    # Verma and Truhlar's meta-GGA (keeps the gap of narrow-gap dots open); grid as in Ivan's hle_sp.in
    "HLE17": """    &XC
      &XC_FUNCTIONAL
        &MGGA_XC_HLE17
        &END MGGA_XC_HLE17
      &END XC_FUNCTIONAL
      &XC_GRID
        XC_DERIV  spline3
        XC_SMOOTH_RHO  none
      &END
    &END""",
}
# Converged on the forces only (CP2K's force limits; the step limits loose): on the flat energy surface of
# Cl-rich dots BFGS kept stepping 0.04 bohr for 20+ steps with the forces 4x below the limit, then lost
# them again (campaign 2026-10-18, Ivan: a plateau with converged forces is accepted).
GEO_CONVERGENCE = """    MAX_FORCE 4.5E-4
    RMS_FORCE 3.0E-4
    MAX_DR 1.0
    RMS_DR 1.0
"""


def geo_motion(lbfgs: bool) -> str:
    return "&MOTION\n  &GEO_OPT\n" + ("    OPTIMIZER LBFGS\n" if lbfgs else "") + GEO_CONVERGENCE + "  &END\n&END\n"


def to_bfgs(inp: str) -> str:
    """An existing geo_opt input with BFGS (CP2K's default optimizer) instead of LBFGS; the limits are kept."""
    return re.sub(r"^\s*OPTIMIZER\s+LBFGS\s*\n", "", inp, flags=re.M | re.I)


def forces_only(inp: str) -> str:
    """An existing geo_opt input converged on the forces only (GEO_CONVERGENCE), for chains made before."""
    if "MAX_DR" in inp:
        return inp
    if re.search(r"^\s*&GEO_OPT\s*$", inp, flags=re.M):
        return re.sub(r"^(\s*&GEO_OPT\s*\n)", lambda m: m.group(1) + GEO_CONVERGENCE, inp, count=1, flags=re.M)
    return re.sub(r"^(&GLOBAL\s*)$", geo_motion(False).replace("\\", "\\\\") + r"\n\g<1>", inp, count=1, flags=re.M)


LBFGS_MOTION = geo_motion(True)


# One GAPW single point (CUTOFF 400 Ry, 4 grids, non-periodic MULTIPOLE cell) holds the cell's
# real-space grids, which outweigh the basis for small and medium dots: Cd16Se13Cl6 in a 23.6 A cube
# took 13.3 GB over 3 ranks (hyperion, CP2K 2026.2; 4.5 GB per rank), whatever the rank count.
GRID_MEM_REF_GB = 13.5
GRID_MEM_REF_VOL_A3 = 23.631 ** 3


def grid_memory_gb(cell: str) -> float:
    """Memory (GB) of one single point of the dot in its cell, scaled from the reference by volume."""
    a, b, c = (float(x) for x in str(cell).split())
    return GRID_MEM_REF_GB * a * b * c / GRID_MEM_REF_VOL_A3


# ---------------- the cell on a common grid ----------------
# CP2K (pw_grid_info.F) gives a cell edge of length L (bohr) n = 2 floor(L sqrt(Ecut[Ry]) / 2 pi) + 1
# points at the density cutoff, then the next length its FFT accepts: 2^a 3^b 5^c (b <= 3, c <= 2). The
# FFTW3 source lists factors 7 and 11 as well, but CP2K 2026.2 on NHR skips them (a 36.574 Å edge got
# 450 points, not 440 = 2^3 5 11; 15.71 Å got 192, not 189 = 3^3 7; checked 2026-10-17). The spacing L/n then differs slightly from edge to edge
# and from system to system. Every edge here is n h instead, with h just under the largest spacing the
# cutoff allows and n a length CP2K picks for n h: the same grid spacing in every direction of every
# dot, stripped dot and monomer, so energy differences between them are free of grid jitter.
CUTOFF_RY = 400                 # CUTOFF of every template (tests/test_cell.py checks)
BOHR_A = 0.529177210903
GRID_MARGIN = 1e-6              # h = (1 - margin) x the largest spacing


def fft_lengths(limit: int = 20000) -> list[int]:
    """The FFT lengths CP2K picks, up to `limit` (2^a 3^b 5^c, checked on NHR's CP2K 2026.2)."""
    out = set()
    for b in range(4):
        for c in range(3):
            n = 3 ** b * 5 ** c
            while n < limit:
                out.add(n)
                n *= 2
    return sorted(out)


_FFT = fft_lengths()


def grid_points(length_a: float, cutoff_ry: float = CUTOFF_RY) -> int:
    """The points CP2K puts along a cell edge of length_a (Å) at this density cutoff."""
    n_min = 2 * int(np.floor(length_a / BOHR_A * np.sqrt(cutoff_ry) / (2 * np.pi))) + 1
    return next(n for n in _FFT if n >= n_min)


def grid_spacing(cutoff_ry: float = CUTOFF_RY) -> float:
    """The common grid spacing (Å): just under pi / sqrt(Ecut[Ry]) bohr."""
    return (1.0 - GRID_MARGIN) * np.pi / np.sqrt(cutoff_ry) * BOHR_A


def aligned_length(min_length_a: float, cutoff_ry: float = CUTOFF_RY) -> float:
    """The shortest edge >= min_length_a that is n x grid_spacing with CP2K choosing exactly n points."""
    h = grid_spacing(cutoff_ry)
    for n in _FFT:
        if n * h >= min_length_a and grid_points(n * h, cutoff_ry) == n:
            return n * h
    raise ValueError(f"cell edge {min_length_a} Å beyond the FFT lengths")


def cell_abc(coords: np.ndarray, vacuum: float = VACUUM_A, min_length: float = 0.0,
             cutoff_ry: float | None = CUTOFF_RY) -> str:
    """ABC of a non-periodic cell: the extent plus `vacuum` (6 Å each side), at least `min_length`,
    each edge rounded up to the common grid (cutoff_ry None: not rounded)."""
    extent = np.maximum(np.ptp(np.asarray(coords, float), axis=0) + vacuum, min_length)
    if cutoff_ry is None:
        return " ".join(f"{x:.3f}" for x in extent)
    return " ".join(f"{aligned_length(x, cutoff_ry):.6f}" for x in extent)


def cp2k_setup(cp: dict) -> str:
    """Shell lines that make CP2K available: an environment script (own builds) or a module."""
    if cp.get("env_script"):
        return f"unset I_MPI_ROOT\nsource {cp['env_script']}"
    return f"module purge\nmodule load {cp['module']}"


def cp2k_launch(cp: dict, ranks: str = "${SLURM_NTASKS}", step_mem: bool = False) -> str:
    """
    MPI launch prefix: srun (SLURM-native builds) or mpirun (default). step_mem: srun asks for
    $STEP_MEM per node when it is set (runs side by side: without it the first step of a job that
    asks --mem per node takes all of it and the others wait), on $STEP_NODES nodes when it is set
    (a run that fits on one node stays there).
    """
    if cp.get("launcher") == "srun":
        mem = " ${STEP_NODES:+-N $STEP_NODES} ${STEP_MEM:+--mem=$STEP_MEM}" if step_mem else ""
        return f"srun --exact -n {ranks}{mem} --cpu-bind=cores"
    return f"mpirun -np {ranks}"


def cp2k_partition(cp: dict) -> str:
    """An #SBATCH -p line when the profile names the CPU partition CP2K was built for."""
    return f"#SBATCH -p {cp['partition']}\n" if cp.get("partition") else ""


def local_data_dir(cluster: dict, override: str | None = None) -> Path:
    """Where the basis and potential files are read from on this machine (inputs keep the cluster path)."""
    for cand in (override, os.environ.get("QDW_CP2K_DATA"), cluster["cp2k"]["data_dir"]):
        if cand and (Path(cand) / cluster["cp2k"]["basis_file"]).is_file():
            return Path(cand)
    raise FileNotFoundError("CP2K data files not found: pass --data-dir or set QDW_CP2K_DATA")


# GAPW (FORCE_PAW, GAPW_ACCURATE_XCINT, REL_CUTOFF 60) is needed mostly for the forces of the ML-frame labels; the
# chains still to run use the faster GPW with the default REL_CUTOFF (Ivan, 2026-10-21: first above 600 atoms, then for
# all sizes). A chain keeps one method for all its steps (the file METHOD in the chain directory).
GPW_ABOVE_ATOMS = 0
METHODS = ("gapw", "gpw")


def chain_method(n_atoms: int) -> str:
    return "gpw" if n_atoms > GPW_ABOVE_ATOMS else "gapw"


def to_gpw(inp: str) -> str:
    """The same input with GPW: METHOD gpw, no FORCE_PAW / GAPW_ACCURATE_XCINT, REL_CUTOFF at its default."""
    t = re.sub(r"^(\s*METHOD\s+)gapw\b", r"\g<1>gpw", inp, flags=re.M | re.I)
    t = re.sub(r"^\s*(FORCE_PAW|GAPW_ACCURATE_XCINT|REL_CUTOFF)\b.*\n", "", t, flags=re.M | re.I)
    return t


def prepare(
    xyz: str,
    step: str,
    project: str,
    out_dir: str,
    cluster: dict,
    *,
    wfn: str | None = None,
    window: dict | None = None,
    charge: int = 0,
    optimizer: str = "auto",
    data_dir: str | None = None,
    time: str | None = None,
    qos: str | None = None,
    xc: str = "PBE",
    method: str = "gapw",
) -> dict:
    """
    Write `<out_dir>/{cp2k.inp, geom.xyz, job.sh, job.json}` for one step.

    `wfn` is the geo_opt restart file every later step reads (a path relative to
    out_dir or absolute); geo_opt restarts from its own. `window` (from
    `outputs.mo_window`) sets ADDED_MOS for the trexio step. `time` and `qos`
    override the cluster profile (e.g. a short test queue). `method`: gapw (the
    templates) or gpw (to_gpw).
    """
    if method not in METHODS:
        raise ValueError(f"unknown method {method!r}; one of {', '.join(METHODS)}")
    if step not in STEPS:
        raise ValueError(f"unknown step {step!r}; one of {', '.join(STEPS)}")
    cp = cluster["cp2k"]
    symbols, coords = read_xyz(xyz)
    local = local_data_dir(cluster, data_dir)
    kinds = kinds_for(symbols, local / cp["basis_file"], local / cp["potential_file"], cp["basis"], cp["potential"])
    counts = Counter(symbols)
    n_bf = sum(kinds[e].functions * n for e, n in counts.items())
    n_elec = sum(kinds[e].valence * n for e, n in counts.items()) - charge

    values = {
        "CHARGE": charge,
        "BASIS_FILE": f"{cp['data_dir']}/{cp['basis_file']}",
        "POTENTIAL_FILE": f"{cp['data_dir']}/{cp['potential_file']}",
        "COORD_FILE": "geom.xyz",
        "CELL": cell_abc(coords),
        "KINDS": kind_blocks(kinds),
        "PROJECT": project,
        "WFN_FILE": wfn or f"{project}-RESTART.wfn",
        "XC": XC[xc],
    }
    if step not in ("geo_opt", "charges") and not wfn:
        raise ValueError(f"step {step} restarts from the geo_opt .wfn: pass wfn")
    if step == "geo_opt":
        template = "geo_opt.inp"
        use_lbfgs = optimizer.upper() == "LBFGS" or (optimizer == "auto" and LBFGS_ABOVE_ATOMS is not None
                                                     and len(symbols) > LBFGS_ABOVE_ATOMS)
        values["MOTION"] = geo_motion(use_lbfgs)
    elif step == "charges":
        template = "single_point_ot.inp"
        values.update({"DFT_PRINT": CHARGES_PRINT, "SCF_PRINT": IGNORE_SCF_FAILURE, "FORCE_EVAL_PRINT": "",
                       "RUN_TYPE": "ENERGY"})
    elif step in ("pdos", "sample"):
        template = "single_point_ot.inp"
        values.update({
            "DFT_PRINT": PDOS_PRINT if step == "pdos" else "",
            "SCF_PRINT": NO_RESTART,
            "FORCE_EVAL_PRINT": FORCES_PRINT if step == "sample" else "",
            "RUN_TYPE": "ENERGY" if step == "pdos" else "ENERGY_FORCE",
        })
    else:
        if not window:
            raise ValueError("step trexio needs the MO window from the pdos step")
        template = "trexio.inp"
        values["ADDED_MOS"] = int(window["added_mos"])

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    write_xyz(out / "geom.xyz", symbols, coords, f"{project} {step}")
    inp = render((TEMPLATES / "cp2k" / template).read_text(), values)
    (out / "cp2k.inp").write_text(to_gpw(inp) if method == "gpw" else inp)

    res = cp2k_resources(cluster, n_bf)
    job = {
        "JOB_NAME": f"{step}-{project}"[:64],
        "TIME": time or cp["time"].get(step, cp["time"]["sample"]),
        "NTASKS": res.cores,
        "NODES": -(-res.cores // int(cluster["cores_per_node"])),   # as few nodes as the cores need
        "NODE_TASKS": node_tasks(cluster, res.cores),   # spread evenly: --mem is per node
        "MEM": res.mem_per_node,
        "QOS": qos or cluster.get("qos", "regular"),
        # Shell setup the job needs whatever shell submitted it (e.g. Lmod for `module`).
        "SETUP": "\n".join(cluster.get("setup", [])),
        "CP2K_SETUP": cp2k_setup(cp), "LAUNCH": cp2k_launch(cp), "PARTITION": cp2k_partition(cp), "ACCOUNT": sbatch_account(cluster),
        "OMP_NUM_THREADS": cp.get("omp_num_threads", 1),
        "EXECUTABLE": cp["executable"],
        "INPUT": "cp2k.inp",
        "OUTPUT": "cp2k.out",
    }
    (out / "job.sh").write_text(render((TEMPLATES / "slurm" / "cp2k.sh").read_text(), job))
    info = {
        "project": project, "step": step, "cluster": cluster["name"], "n_atoms": len(symbols), "method": method,
        "composition": dict(sorted(counts.items())), "basis_functions": n_bf, "electrons": n_elec,
        "occupied_mos": n_elec // 2, "cell_A": values["CELL"], "wfn": values["WFN_FILE"],
        "kinds": {e: {"basis": k.basis, "potential": k.potential, "valence": k.valence, "functions": k.functions}
                  for e, k in kinds.items()},
        "resources": res.as_dict(),
    }
    if step == "trexio":
        info["window"] = window
    (out / "job.json").write_text(json.dumps(info, indent=2))
    return info


def chain_time(n_bf: int, cores: int, cap: str = "12:00:00") -> str:
    """
    Wall time of a whole chain from the cost of a PBE geometry step on CdSe dots (NHR and
    hyperion runs: 1.3 core-h at 5.9k basis functions, 3.6-4.8 at 20k, 9-11 at 31k, 22-38 at
    34-43k, i.e. ~4.1e-7 n_bf^1.72, conservative in the middle): ~160 step-equivalents
    (60 geo_opt steps, the field runs at ~5 steps each sharing the cores, PDOS, TREXIO),
    1.5 x for safety, at least three hours, at most `cap`. A flat 12 h kept a 35-atom chain
    waiting for hours for a backfill slot; one hour was too short for the whole chain with ESP and
    QDEX (campaign 2026-10-18: 92-357 atoms took 45 min to 1 h 52 min, 5 of 6 small chains timed out).
    """
    hours = 1.5 * 160 * 4.1e-7 * max(n_bf, 1) ** 1.72 / max(cores, 1)
    h, m, s = (int(x) for x in cap.split(":"))
    hours = min(max(hours, 3.0), h + m / 60 + s / 3600)
    total = int(round(hours * 60))
    return f"{total // 60:02d}:{total % 60:02d}:00"


def chain(xyz: str, project: str, out_dir: str, cluster: dict, *, cluster_ref: str = "default",
          record: str | None = None, charge: int = 0, optimizer: str = "auto", data_dir: str | None = None,
          time: str | None = None, qos: str | None = None, field: bool = True) -> dict:
    """
    One job for the whole CP2K track of a structure: `<out_dir>/geo_opt/` is
    prepared now; `<out_dir>/chain.sh` runs geo_opt, the finite-field single points
    (`field`: Born charges, Raman tensors, Hirshfeld charges, Orchestr.AI qd.field),
    then prepares and runs the PDOS and TREXIO steps and trims the TREXIO file,
    skipping finished steps. With a `record` (record.json), QDEX runs last on the
    trimmed orbitals, and the field step adds the record's monomers (props/monomers.xyz).
    Cores and memory are sized for the structure, the same for every step.
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    if record:
        (out / "record.json").write_text(Path(record).read_text())
    method = (out / "METHOD").read_text().strip() if (out / "METHOD").is_file() else chain_method(len(read_xyz(xyz)[0]))
    (out / "METHOD").write_text(method + "\n")      # every step of the chain (chain.sh: --method $(cat METHOD))
    info = prepare(xyz, "geo_opt", project, str(out / "geo_opt"), cluster, charge=charge,
                   optimizer=optimizer, data_dir=data_dir, time=time, qos=qos, method=method)
    cp = cluster["cp2k"]
    res = info["resources"]
    py = cluster.get("qdw_python")
    if not py:
        raise KeyError("cluster profile needs qdw_python (the python with qd_workflows installed)")
    total_time = time or chain_time(info["basis_functions"], res["cores"], cp["time"].get("chain", cp["time"]["geo_opt"]))
    job = {
        "JOB_NAME": f"cp2k-{project}"[:64], "TIME": total_time, "NTASKS": res["cores"],
        "NODES": -(-res["cores"] // int(cluster["cores_per_node"])), "NODE_TASKS": node_tasks(cluster, res["cores"]),
        "MEM": res["mem_per_node"],
        "QOS": qos or cluster.get("qos", "regular"), "SETUP": "\n".join(cluster.get("setup", [])),
        "CP2K_SETUP": cp2k_setup(cp), "LAUNCH": cp2k_launch(cp), "LAUNCH_N": cp2k_launch(cp, "$1", step_mem=True),
        "PARTITION": cp2k_partition(cp), "ACCOUNT": sbatch_account(cluster),
        "OMP_NUM_THREADS": cp.get("omp_num_threads", 1), "EXECUTABLE": cp["executable"],
        "PROJECT": project, "QDW_PYTHON": py, "QDW_PYTHONPATH": cluster.get("qdw_pythonpath", ""),
        "CLUSTER": cluster_ref, "QDEX": cluster.get("qdex", {}).get("command", "qdex"),
        "FIELD": "1" if field else "0", "FIELD_PARALLEL": int(cp.get("field_parallel", 13)),
        "ESP_PARALLEL": int(cp.get("esp_parallel", 3)),
        "GRID_MEM_GB": f"{grid_memory_gb(info['cell_A']):.0f}", "GRID_MEM_MB": int(1024 * grid_memory_gb(info["cell_A"])),
        "RECORD_DIR": str(Path(record).resolve().parent) if record else "",
        "MONOMER_LIBRARY": cp.get("monomer_library", ""),
        "QDW_LOG": os.environ.get("QDW_LOG", ""),      # the run's progress log, if it has one (qdw log)
    }
    (out / "chain.sh").write_text(render((TEMPLATES / "slurm" / "cp2k_chain.sh").read_text(), job))
    return {**info, "chain": str(out / "chain.sh")}
