"""Orchestr.AI PROPS runs (MACE-MH-1 / xtb properties of library records) on a GPU cluster."""
from __future__ import annotations

import os
from pathlib import Path

import yaml

from . import TEMPLATES, render

ALL_STEPS = ["relax", "structure", "hessian", "vibspec", "electronic", "stability", "detachment",
             "solvation", "sites", "report"]


RUN = "python -m orchestr_ai.postprocessing config.yaml"
# Short QoS (e.g. 10-min test jobs): stop a minute before the limit and submit the
# next piece. Finished steps, desorption relaxations and Raman modes are checkpointed,
# so every piece continues where the last one stopped.
RUN_PIECES = """timeout {seconds} {run}; rc=$?
if [ $rc -eq 124 ]; then echo "[props] time limit: next piece"; cd "$SLURM_SUBMIT_DIR" && sbatch job.sh && exit 0; fi
exit $rc"""


def _seconds(t: str) -> int:
    days, _, hms = t.rpartition("-")
    parts = [int(x) for x in hms.split(":")]
    while len(parts) < 3:
        parts.insert(0, 0)
    return int(days or 0) * 86400 + parts[0] * 3600 + parts[1] * 60 + parts[2]


def prepare(records, out_dir: str, cluster: dict, *, steps=None, time: str | None = None,
            qos: str | None = None, max_atoms: int | None = None, name: str = "props",
            pieces: bool = False) -> dict:
    """
    Write `<out_dir>/{config.yaml, job.sh}`: one GPU job running `run_type: PROPS`
    over `records` (directories with record.json and start.xyz). Paths to the
    MACE model, g-xTB, the bulk CIFs and the shared references come from the
    cluster profile's `props` section. With `pieces`, the job resubmits itself
    until the run finishes (for queues where only short jobs start soon).
    """
    pp = cluster["props"]
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    recs = [str(Path(r).resolve()) for r in records]
    missing = [r for r in recs if not (Path(r) / "record.json").is_file()]
    if missing:
        raise FileNotFoundError(f"no record.json in {', '.join(missing)}")
    config = {
        "run_type": "PROPS",
        "platform": "mace",
        "model_path": pp["model"],
        "mace_head": pp.get("head", "omat_pbe"),
        "props": {"records": recs, "steps": list(steps or ALL_STEPS), "force": False,
                  "settings": {"device": "auto", "dtype": "float64"}},
    }
    if max_atoms:
        config["props"]["max_atoms"] = int(max_atoms)
    (out / "config.yaml").write_text(yaml.safe_dump(config, sort_keys=False))
    wall = time or pp.get("time", "01:00:00")
    run = RUN_PIECES.format(seconds=max(_seconds(wall) - 60, 60), run=RUN) if pieces else RUN
    (out / "job.sh").write_text(_job(cluster, name, wall, qos, run))
    return {"records": recs, "config": str(out / "config.yaml"), "job": str(out / "job.sh")}


def bench(out_dir: str, cluster: dict, *, sizes: str = "100,300,1000,2000,5000", time: str = "01:00:00",
          qos: str | None = None, extra: str = "") -> dict:
    """GPU job for the size benchmark (orchestr_ai.qd.bench): force calls and Hessians vs atom count."""
    pp = cluster["props"]
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    run = (f"python -m orchestr_ai.qd.bench --model {pp['model']} --head {pp.get('head', 'omat_pbe')} "
           f"--sizes {sizes} -o bench.json {extra}".rstrip())
    (out / "job.sh").write_text(_job(cluster, "bench-mace", time, qos, run))
    return {"job": str(out / "job.sh"), "run": run}


def _job(cluster: dict, name: str, wall: str, qos: str | None, run: str) -> str:
    pp = cluster["props"]
    job = {
        "JOB_NAME": name[:64], "TIME": wall, "RUN": run, "PARTITION": pp["partition"],
        "QOS": qos or cluster.get("qos", "regular"), "GRES": pp["gres"],
        "CONSTRAINT": f"#SBATCH --constraint={pp['constraint']}\n" if pp.get("constraint") else "", "CPUS": pp.get("cpus", 8),
        "MEM": pp.get("mem", "32G"), "SETUP": "\n".join(cluster.get("setup", [])),
        "CONDA_INIT": pp.get("conda_init", ""), "CONDA_ENV": pp["conda_env"], "MODEL": pp["model"],
        "GXTB": pp["gxtb"], "CIF_DIRS": os.pathsep.join(pp.get("cif_dirs", [])), "REFS": pp["refs"],
    }
    return render((TEMPLATES / "slurm" / "props_gpu.sh").read_text(), job)
