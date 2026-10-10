"""Orchestr.AI PROPS runs (MACE-MH-1 / xtb properties of library records) on a GPU cluster."""
from __future__ import annotations

import os
from pathlib import Path

import yaml

from . import TEMPLATES, render

ALL_STEPS = ["relax", "structure", "hessian", "vibspec", "electronic", "stability", "detachment",
             "solvation", "sites", "report", "wigner", "md"]


RUN = "python -m orchestr_ai.postprocessing config.yaml"
# Short QoS (e.g. 10-min test jobs): every piece first queues its successor (it starts
# when this one ends, however it ends), stops a minute before the limit, and cancels
# the successor once the run is done. Queuing at the start keeps the job slot: a
# per-user submit limit cannot break the chain at the handover. Finished steps,
# desorption relaxations and Raman modes are checkpointed, so every piece continues
# where the last one stopped.
RUN_PIECES = """queue_next() {{ (cd "$SLURM_SUBMIT_DIR" && sbatch --parsable --dependency=afterany:$SLURM_JOB_ID job.sh 2>/dev/null) > .next_piece.$SLURM_JOB_ID; [ -s .next_piece.$SLURM_JOB_ID ]; }}
# the submit limit may be reached now (another job queued): keep trying in the background
(for i in $(seq 1 18); do queue_next && break; sleep 30; done) &
retry=$!
timeout {seconds} {run}; rc=$?
kill $retry 2>/dev/null; wait $retry 2>/dev/null
next=$(cat .next_piece.$SLURM_JOB_ID 2>/dev/null); rm -f .next_piece.$SLURM_JOB_ID
if [ $rc -eq 124 ]; then
    [ -n "$next" ] || {{ queue_next && next=$(cat .next_piece.$SLURM_JOB_ID); rm -f .next_piece.$SLURM_JOB_ID; }}
    echo "[props] time limit: next piece ${{next:-NOT QUEUED (submit limit)}}"; exit 0
fi
[ -n "$next" ] && scancel "$next" && echo "[props] run ended (exit $rc): cancelled $next"
exit $rc"""


def _seconds(t: str) -> int:
    days, _, hms = t.rpartition("-")
    parts = [int(x) for x in hms.split(":")]
    while len(parts) < 3:
        parts.insert(0, 0)
    return int(days or 0) * 86400 + parts[0] * 3600 + parts[1] * 60 + parts[2]


def prepare(records, out_dir: str, cluster: dict, *, steps=None, time: str | None = None,
            qos: str | None = None, max_atoms: int | None = None, name: str = "props",
            pieces: bool = False, md: bool = False) -> dict:
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
                  # md runs only for the MD subset (wigner always)
                  "settings": {"device": "auto", "dtype": "float64", "md_enabled": bool(md)}},
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
        "EXTRA_ENV": "\n".join(f"export {k}={v}" for k, v in (pp.get("env") or {}).items()),
    }
    return render((TEMPLATES / "slurm" / "props_gpu.sh").read_text(), job)
