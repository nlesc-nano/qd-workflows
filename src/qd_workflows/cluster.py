"""HPC cluster profiles (clusters/<name>.yaml) and the resources a job asks for."""
from __future__ import annotations

import math
import os
from dataclasses import dataclass
from pathlib import Path

import yaml

from . import REPO_DIR


def cluster_dirs() -> list[Path]:
    dirs = [Path(p) for p in os.environ.get("QDW_CLUSTERS", "").split(os.pathsep) if p]
    return dirs + [REPO_DIR / "clusters"]


def load_cluster(name_or_path: str = "default") -> dict:
    """A cluster profile by name (looked up in $QDW_CLUSTERS, then clusters/) or by path."""
    path = Path(name_or_path)
    if not path.is_file():
        path = next((d / f"{name_or_path}.yaml" for d in cluster_dirs() if (d / f"{name_or_path}.yaml").is_file()), None)
        if path is None:
            raise FileNotFoundError(f"no cluster profile {name_or_path!r} in {', '.join(map(str, cluster_dirs()))}")
    return yaml.safe_load(path.read_text())


def sbatch_account(cluster: dict) -> str:
    """An #SBATCH --account line when the profile names the account to charge (else the user's default)."""
    return f"#SBATCH --account={cluster['account']}\n" if cluster.get("account") else ""


@dataclass
class Resources:
    cores: int
    mem_per_cpu: str      # per core
    mem_per_node: str     # what the job asks for (#SBATCH --mem): the fullest node's cores x mem_per_cpu
    memory_gb: float      # estimated total

    def as_dict(self) -> dict:
        return {"cores": self.cores, "mem_per_cpu": self.mem_per_cpu, "mem_per_node": self.mem_per_node,
                "memory_gb": round(self.memory_gb, 1)}


def cp2k_resources(cluster: dict, basis_functions: int) -> Resources:
    """
    Cores and memory for a CP2K job. Both grow with the square of the number of
    basis functions (dense matrices), scaled from the profile's calibration run.
    """
    cal = cluster["cp2k"]["calibration"]
    scale = (basis_functions / cal["basis_functions"]) ** 2
    cores = max(int(cluster.get("min_cores", 1)), math.ceil(cal["cores"] * scale))
    if cluster.get("request", "cores") == "nodes":
        per_node = int(cluster["cores_per_node"])
        cores = per_node * math.ceil(cores / per_node)
    memory_gb = cal["memory_gb"] * scale
    floor_gb = _gb(cluster.get("mem_per_cpu", "2g"))
    mem = max(floor_gb, math.ceil(memory_gb / cores))
    if cluster.get("node_memory_gb"):
        tasks = node_tasks(cluster, cores)
        mem = min(mem, int(cluster["node_memory_gb"] // tasks))
        if mem < 1:
            raise ValueError(f"{tasks} ranks per node do not fit in {cluster['node_memory_gb']} GB")
    return Resources(cores=cores, mem_per_cpu=f"{mem}g", mem_per_node=node_memory(cluster, cores, f"{mem}g"),
                     memory_gb=memory_gb)


def node_tasks(cluster: dict, cores: int) -> int:
    """Ranks per node of a job with `cores` ranks spread evenly on as few nodes as they need
    (#SBATCH --ntasks-per-node: Slurm would otherwise fill the first node, e.g. 192 + 10)."""
    nodes = -(-cores // int(cluster["cores_per_node"]))
    return -(-cores // nodes)


def node_memory(cluster: dict, cores: int, mem_per_cpu: str) -> str:
    """
    The job's memory per node (#SBATCH --mem) for `cores` ranks of `mem_per_cpu` each, at most what a
    node holds (`node_memory_gb`). Jobs ask per node, not per CPU: on NHR's cpu-genoa a job with
    -n 48 --mem-per-cpu=7g is planned on an idle node every backfill cycle and never started,
    while --mem=336G starts at once (and Slurm counts --mem-per-cpu per hardware thread).
    """
    gb = node_tasks(cluster, cores) * _gb(mem_per_cpu)
    if cluster.get("node_memory_gb"):
        gb = min(gb, int(cluster["node_memory_gb"]))
    return f"{max(gb, 1)}G"


def _gb(text: str) -> int:
    text = str(text).strip().lower()
    if text.endswith("g"):
        return int(float(text[:-1]))
    if text.endswith("m"):
        return max(1, math.ceil(float(text[:-1]) / 1024))
    return int(float(text))
