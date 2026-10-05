"""CP2K inputs and job scripts for one structure and one step."""
from __future__ import annotations

import json
import os
from collections import Counter
from pathlib import Path

import numpy as np

from .. import TEMPLATES, render
from ..cluster import cp2k_resources
from ..xyz import read_xyz, write_xyz
from .basis import kind_blocks, kinds_for

STEPS = ("geo_opt", "pdos", "trexio", "sample")
VACUUM_A = 12.0          # added to the extent along each axis (6 Å on each side)
LBFGS_ABOVE_ATOMS = 1000

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
LBFGS_MOTION = """&MOTION
  &GEO_OPT
    OPTIMIZER LBFGS
  &END
&END
"""


def cell_abc(coords: np.ndarray, vacuum: float = VACUUM_A) -> str:
    extent = np.ptp(np.asarray(coords, float), axis=0) + vacuum
    return " ".join(f"{x:.3f}" for x in extent)


def local_data_dir(cluster: dict, override: str | None = None) -> Path:
    """Where the basis and potential files are read from on this machine (inputs keep the cluster path)."""
    for cand in (override, os.environ.get("QDW_CP2K_DATA"), cluster["cp2k"]["data_dir"]):
        if cand and (Path(cand) / cluster["cp2k"]["basis_file"]).is_file():
            return Path(cand)
    raise FileNotFoundError("CP2K data files not found: pass --data-dir or set QDW_CP2K_DATA")


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
) -> dict:
    """
    Write `<out_dir>/{cp2k.inp, geom.xyz, job.sh, job.json}` for one step.

    `wfn` is the geo_opt restart file every later step reads (a path relative to
    out_dir or absolute); geo_opt restarts from its own. `window` (from
    `outputs.mo_window`) sets ADDED_MOS for the trexio step. `time` and `qos`
    override the cluster profile (e.g. a short test queue).
    """
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
    }
    if step != "geo_opt" and not wfn:
        raise ValueError(f"step {step} restarts from the geo_opt .wfn: pass wfn")
    if step == "geo_opt":
        template = "geo_opt.inp"
        use_lbfgs = optimizer.upper() == "LBFGS" or (optimizer == "auto" and len(symbols) > LBFGS_ABOVE_ATOMS)
        values["MOTION"] = LBFGS_MOTION if use_lbfgs else ""
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
    (out / "cp2k.inp").write_text(render((TEMPLATES / "cp2k" / template).read_text(), values))

    res = cp2k_resources(cluster, n_bf)
    job = {
        "JOB_NAME": f"{step}-{project}"[:64],
        "TIME": time or cp["time"][step],
        "NTASKS": res.cores,
        "MEM_PER_CPU": res.mem_per_cpu,
        "QOS": qos or cluster.get("qos", "regular"),
        "MODULE": cp["module"],
        "OMP_NUM_THREADS": cp.get("omp_num_threads", 1),
        "EXECUTABLE": cp["executable"],
        "INPUT": "cp2k.inp",
        "OUTPUT": "cp2k.out",
    }
    (out / "job.sh").write_text(render((TEMPLATES / "slurm" / "cp2k.sh").read_text(), job))
    info = {
        "project": project, "step": step, "cluster": cluster["name"], "n_atoms": len(symbols),
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
