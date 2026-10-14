"""A reaction run directory and its jobs: MACE + ALPB (GPU), CP2K (CPU), then the analysis."""
from __future__ import annotations

import json
import shutil
from pathlib import Path

import yaml

from .. import props
from ..xyz import read_xyz, write_xyz
from . import cp2k, spec

SUBMIT = """#!/bin/bash
# MACE (GPU) and CP2K (CPU) run side by side; the analysis starts when both have ended.
cd "$(dirname "$0")"
gpu=$(sbatch --parsable job_gpu.sh) || exit 1
{cpu_line}
sbatch --dependency=afterany:$gpu{cpu_dep} job_analyse.sh
"""


SUBMIT_TEST = """#!/bin/bash
# Test queues: MACE (GPU), then the CP2K pieces, which submit the analysis when every run has ended.
cd "$(dirname "$0")"
gpu=$(sbatch --parsable job_gpu.sh) || exit 1
{cpu_line}
"""


def _to_test(script: str, partition: str, test_partition: str | None) -> str:
    """The same job in the test partition with the 1 h limit."""
    if not test_partition:
        raise KeyError("the cluster profile names no test partition")
    lines = []
    for line in script.splitlines():
        if line.startswith("#SBATCH -p "):
            line = f"#SBATCH -p {test_partition}"
        elif line.startswith("#SBATCH -t "):
            line = "#SBATCH -t 01:00:00"
        lines.append(line)
    return "\n".join(lines) + "\n"


def prepare(path: str, out_dir: str, cluster: dict, *, gpu_time: str = "01:00:00", cpu_time: str | None = None,
            qos: str | None = None, data_dir: str | None = None, test_queue: bool = False) -> dict:
    """
    The run directory. test_queue: every job in the profile's test partitions (1 h limit); the
    CP2K job then runs in pieces that resubmit themselves and submits the analysis at the end.
    """
    rx = spec.load(path)
    run = Path(out_dir).resolve()
    run.mkdir(parents=True, exist_ok=True)
    shutil.copy(path, run / "reaction.yaml")
    plan = {}
    for name, sp in rx["species"].items():
        sdir = run / "species" / name
        sdir.mkdir(parents=True, exist_ok=True)
        sym, pos = read_xyz(sp["xyz"])
        write_xyz(sdir / "dft.xyz", sym, pos, f"{name} DFT structure (from {sp['xyz']}), E = {sp['energy_input']}")
        if sp["charges_file"]:
            shutil.copy(sp["charges_file"], sdir / "charges.txt")
        if "cp2k" in cluster:
            plan[name] = cp2k.prepare_species(sp, sdir, cluster, rx["settings"], data_dir)
    (run / "reaction.json").write_text(json.dumps(rx, indent=1))

    py_path = cluster.get("qdw_pythonpath", "")
    # paths relative to the submit directory: the run can be prepared here and copied to the cluster
    run_gpu = (f"export PYTHONPATH={py_path}${{PYTHONPATH:+:$PYTHONPATH}}\ncd \"$SLURM_SUBMIT_DIR\"\n"
               "python -m qd_workflows.reaction.mace .; rc=$?\n"
               "python -m qd_workflows.reaction.solvation alpb .\nexit $rc")
    gpu_job = props._job(cluster, f"rx-mace-{run.name}", gpu_time, qos, run_gpu)
    if test_queue:
        gpu_job = _to_test(gpu_job, cluster["props"]["partition"], cluster["props"].get("test_partition"))
    (run / "job_gpu.sh").write_text(gpu_job)
    run_an = (f"export PYTHONPATH={py_path}${{PYTHONPATH:+:$PYTHONPATH}}\ncd \"$SLURM_SUBMIT_DIR\"\n"
              "python -m qd_workflows.cli reaction analyse .")
    an = props._job(cluster, f"rx-analyse-{run.name}", "00:30:00", qos, run_an, cpu=True)
    an = an.replace(f"#SBATCH -c {cluster['props']['cpu'].get('cpus', 8)}", "#SBATCH -c 4").replace(
        f"#SBATCH --mem={cluster['props']['cpu'].get('mem', '32G')}", "#SBATCH --mem=16G")
    if test_queue:
        an = _to_test(an, cluster["props"]["cpu"]["partition"], cluster["cp2k"].get("test_partition"))
    (run / "job_analyse.sh").write_text(an)
    has_cp2k = any(p["runs"] for p in plan.values())
    if has_cp2k:
        (run / "job_cp2k.sh").write_text(cp2k.job(run, cluster, plan, time=cpu_time, qos=qos, test_queue=test_queue,
                                                  threads=int(rx["settings"].get("omp_threads", 4))))   # sets the ranks
    if test_queue:      # one running job per user there: MACE, then CP2K pieces (they submit the analysis)
        (run / "submit.sh").write_text(SUBMIT_TEST.format(cpu_line="sbatch --dependency=afterany:$gpu job_cp2k.sh"
                                                          if has_cp2k else "sbatch --dependency=afterany:$gpu job_analyse.sh"))
    else:
        (run / "submit.sh").write_text(SUBMIT.format(
        cpu_line="cpu=$(sbatch --parsable job_cp2k.sh) || exit 1" if has_cp2k else "",
        cpu_dep=":$cpu" if has_cp2k else ""))
    (run / "plan.yaml").write_text(yaml.safe_dump({"cp2k": plan}, sort_keys=False))
    return {"run": str(run), "species": {n: {"n_atoms": s["n_atoms"], "nu": rx["nu"][n], **plan.get(n, {})}
                                         for n, s in rx["species"].items()}}
