"""Orchestr.AI PROPS runs (MACE-MH-1 / xtb properties of library records) on a GPU cluster."""
from __future__ import annotations

import os
from pathlib import Path

import yaml

from . import TEMPLATES, render

ALL_STEPS = ["relax", "structure", "hessian", "vibspec", "electronic", "stability", "detachment",
             "solvation", "sites", "report", "wigner", "md"]
# With a `props.cpu` section in the cluster profile the run is split in two jobs: MACE-MH-1
# steps on a GPU, then the xTB steps (and the report) on a CPU node, which the GPU job submits
# when it ends. g-xTB and GFN2-xTB are CPU programs: on the GPU node they left the A100 idle
# (54 min of g-xTB Raman for a 141-atom dot).
GPU_STEPS = ["relax", "structure", "hessian", "stability", "detachment", "wigner", "md"]
CPU_STEPS = ["electronic", "vibspec", "solvation", "sites", "report"]


RUN = "python -m orchestr_ai.postprocessing {config}"
# Short QoS (e.g. 10-min test jobs): every piece first queues its successor (it starts
# when this one ends, however it ends), stops a minute before the limit, and cancels
# the successor once the run is done. Queuing at the start keeps the job slot: a
# per-user submit limit cannot break the chain at the handover. Finished steps,
# desorption relaxations and Raman modes are checkpointed, so every piece continues
# where the last one stopped.
RUN_PIECES = """queue_next() {{ (cd "$SLURM_SUBMIT_DIR" && sbatch --parsable --dependency=afterany:$SLURM_JOB_ID {job} 2>/dev/null) > .next_piece.$SLURM_JOB_ID; [ -s .next_piece.$SLURM_JOB_ID ]; }}
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
{then}exit $rc"""
# The GPU stage hands over to the CPU stage when its run ends (not at a piece's time limit).
# Records whose GPU steps failed then fail in the CPU stage, which needs them cached.
THEN_CPU = """(cd "$SLURM_SUBMIT_DIR" && sbatch job_cpu.sh) || echo "[props] could not submit the CPU stage: sbatch job_cpu.sh"
"""


def _seconds(t: str) -> int:
    days, _, hms = t.rpartition("-")
    parts = [int(x) for x in hms.split(":")]
    while len(parts) < 3:
        parts.insert(0, 0)
    return int(days or 0) * 86400 + parts[0] * 3600 + parts[1] * 60 + parts[2]


def prepare(records, out_dir: str, cluster: dict, *, steps=None, time: str | None = None,
            qos: str | None = None, max_atoms: int | None = None, name: str = "props",
            pieces: bool = False, md: bool = False, cpu_time: str | None = None) -> dict:
    """
    Write `<out_dir>/{config.yaml, job.sh}`: a GPU job running `run_type: PROPS` over
    `records` (directories with record.json and start.xyz). If the cluster profile has a
    `props.cpu` section, the xTB steps go to `{config_cpu.yaml, job_cpu.sh}`, a CPU job
    that the GPU job submits when it ends and that only reads the GPU steps from the cache.
    Paths to the MACE model, g-xTB, the bulk CIFs and the shared references come from the
    profile's `props` section. With `pieces`, each job resubmits itself until its run
    finishes (for queues where only short jobs start soon).
    """
    pp = cluster["props"]
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    recs = [str(Path(r).resolve()) for r in records]
    missing = [r for r in recs if not (Path(r) / "record.json").is_file()]
    if missing:
        raise FileNotFoundError(f"no record.json in {', '.join(missing)}")
    steps = list(steps or ALL_STEPS)
    split = bool(pp.get("cpu"))
    gpu_steps = [s for s in steps if s in GPU_STEPS] if split else steps
    cpu_steps = [s for s in steps if s in CPU_STEPS] if split else []

    def config(stage_steps, upstream=None):
        c = {"run_type": "PROPS", "platform": "mace", "model_path": pp["model"],
             "mace_head": pp.get("head", "omat_pbe"),
             "props": {"records": recs, "steps": stage_steps, "force": False,
                       # md runs only for the MD subset (wigner always)
                       "settings": {"device": "auto", "dtype": "float64", "md_enabled": bool(md)}}}
        if upstream:
            c["props"]["upstream"] = upstream
        if max_atoms:
            c["props"]["max_atoms"] = int(max_atoms)
        return c

    def runner(cfg, job, wall, then=""):
        run = RUN.format(config=cfg)
        if pieces:
            return RUN_PIECES.format(seconds=max(_seconds(wall) - 60, 60), run=run, job=job, then=then)
        return f"{run}; rc=$?\n{then}exit $rc" if then else run

    result = {"records": recs}
    if gpu_steps:
        wall = time or pp.get("time", "01:00:00")
        (out / "config.yaml").write_text(yaml.safe_dump(config(gpu_steps), sort_keys=False))
        (out / "job.sh").write_text(_job(cluster, name, wall, qos,
                                         runner("config.yaml", "job.sh", wall, THEN_CPU if cpu_steps else "")))
        result.update({"config": str(out / "config.yaml"), "job": str(out / "job.sh"), "gpu_steps": gpu_steps})
    if cpu_steps:
        wall = cpu_time or pp["cpu"].get("time", "01:00:00")
        (out / "config_cpu.yaml").write_text(yaml.safe_dump(config(cpu_steps, "cached"), sort_keys=False))
        (out / "job_cpu.sh").write_text(_job(cluster, f"{name}-cpu", wall, qos,
                                             runner("config_cpu.yaml", "job_cpu.sh", wall), cpu=True))
        result.update({"config_cpu": str(out / "config_cpu.yaml"), "job_cpu": str(out / "job_cpu.sh"),
                       "cpu_steps": cpu_steps})
        if not gpu_steps:
            result["job"] = str(out / "job_cpu.sh")
    return result


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


def _job(cluster: dict, name: str, wall: str, qos: str | None, run: str, cpu: bool = False) -> str:
    pp = cluster["props"]
    part = pp["cpu"] if cpu else pp
    job = {
        "JOB_NAME": name[:64], "TIME": wall, "RUN": run, "PARTITION": part["partition"],
        "QOS": qos or cluster.get("qos", "regular"), "GRES": pp.get("gres", ""),
        "CONSTRAINT": f"#SBATCH --constraint={part['constraint']}\n" if part.get("constraint") else "",
        "CPUS": part.get("cpus", 8), "MEM": part.get("mem", "32G"), "SETUP": "\n".join(cluster.get("setup", [])),
        "CONDA_INIT": pp.get("conda_init", ""), "CONDA_ENV": pp["conda_env"], "MODEL": pp["model"],
        "GXTB": pp["gxtb"], "CIF_DIRS": os.pathsep.join(pp.get("cif_dirs", [])), "REFS": pp["refs"],
        "EXTRA_ENV": "\n".join(f"export {k}={v}" for k, v in (pp.get("env") or {}).items()),
    }
    template = "props_cpu.sh" if cpu else "props_gpu.sh"
    return render((TEMPLATES / "slurm" / template).read_text(), job)
