"""Orchestr.AI PROPS runs (MACE-MH-1 / xtb properties of library records) on a GPU cluster."""
from __future__ import annotations

import os
from pathlib import Path

import yaml

from . import TEMPLATES, render

ALL_STEPS = ["relax", "structure", "hessian", "vibspec", "electronic", "stability", "detachment",
             "solvation", "sites", "report"]


def prepare(records, out_dir: str, cluster: dict, *, steps=None, time: str | None = None,
            qos: str | None = None, max_atoms: int | None = None, name: str = "props") -> dict:
    """
    Write `<out_dir>/{config.yaml, job.sh}`: one GPU job running `run_type: PROPS`
    over `records` (directories with record.json and start.xyz). Paths to the
    MACE model, g-xTB, the bulk CIFs and the shared references come from the
    cluster profile's `props` section.
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
    job = {
        "JOB_NAME": name[:64], "TIME": time or pp.get("time", "01:00:00"), "PARTITION": pp["partition"],
        "QOS": qos or cluster.get("qos", "regular"), "GRES": pp["gres"], "CPUS": pp.get("cpus", 8),
        "MEM": pp.get("mem", "32G"), "SETUP": "\n".join(cluster.get("setup", [])),
        "CONDA_INIT": pp.get("conda_init", ""), "CONDA_ENV": pp["conda_env"], "MODEL": pp["model"],
        "GXTB": pp["gxtb"], "CIF_DIRS": os.pathsep.join(pp.get("cif_dirs", [])), "REFS": pp["refs"],
    }
    (out / "job.sh").write_text(render((TEMPLATES / "slurm" / "props_gpu.sh").read_text(), job))
    return {"records": recs, "config": str(out / "config.yaml"), "job": str(out / "job.sh")}
