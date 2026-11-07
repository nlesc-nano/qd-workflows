"""Batch scripts for the GPU steps (reference, states, neb) from the recipe's cluster block: one GPU, the run
directory as working directory, the step run through `python -m qd_workflows.cli sites <step>`."""
from __future__ import annotations

from pathlib import Path

GPU_STEPS = ("reference", "states", "neb")


def script(r: dict, run_dir: str, steps, recipe_path: str | None = None, time: str = "02:00:00", name: str = "sitekin",
           extra: str = "") -> str:
    c = r["cluster"]
    for s in steps:
        if s not in GPU_STEPS:
            raise ValueError(f"{s!r} is not a GPU step ({', '.join(GPU_STEPS)})")
    rec = f" --recipe {recipe_path}" if recipe_path else ""      # none: the packaged template on the cluster
    cmds = " && ".join(f"python -u -m qd_workflows.cli sites {s} {run_dir}{rec}{(' ' + extra) if extra else ''} > log_{s}.txt 2>&1"
                       for s in steps)
    run = c["run"].replace("{cmd}", cmds)                  # the template is written as it goes in the script
    return "\n".join([
        "#!/bin/bash",
        f"#SBATCH -A {c['account']}",
        f"#SBATCH -p {c['partition']}",
        f"#SBATCH -t {time}",
        "#SBATCH -N 1",
        f"#SBATCH --gpus-per-node={c.get('gpus', 1)}",
        "#SBATCH --ntasks-per-node=1",
        f"#SBATCH -c {c.get('cpus', 7)}",
        f"#SBATCH --mem={c.get('mem', '60G')}",
        f"#SBATCH -J {name}",
        f"#SBATCH -o {name}-%j.out",
        f"# site-kinetics protocol, steps {', '.join(steps)} (qdw sites)",
        f"cd {run_dir}",
        f"source {c['env']}",
        run, ""])


def write(r: dict, run_dir: str, steps, recipe_path: str | None, path: Path, **kw) -> Path:
    path = Path(path)
    path.write_text(script(r, run_dir, steps, recipe_path, **kw))
    return path
