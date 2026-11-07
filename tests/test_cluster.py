"""Cores and memory of CP2K jobs: memory asked per node (#SBATCH --mem), runs side by side per step."""
from pathlib import Path

from qd_workflows.cluster import cp2k_resources, load_cluster, node_memory
from qd_workflows.cp2k.inputs import cp2k_launch
from qd_workflows.reaction import cp2k as rx


def test_nhr_memory_per_node():
    nhr = load_cluster("nhr")                         # request: nodes, every CP2K job on whole nodes
    small = cp2k_resources(nhr, 3000)                 # a small dot: one whole node, all 192 cores as ranks
    assert (small.cores, small.mem_per_cpu, small.mem_per_node) == (192, "7g", "1344G")
    full = cp2k_resources(nhr, 21000)                 # 184 cores asked: rounded up to the node
    assert full.cores == 192 and full.mem_per_node == "1344G"
    two = cp2k_resources(nhr, 24000)                  # 240 cores asked: two whole nodes
    assert two.cores == 384 and two.mem_per_node == "1344G"
    assert node_memory(nhr, 192, "8g") == "1400G"     # capped by what a node holds


def test_side_by_side_steps_ask_for_their_memory():
    cp = {"launcher": "srun"}
    assert "--mem" not in cp2k_launch(cp)             # one step on all ranks: the whole job's memory
    assert cp2k_launch(cp, "$1", step_mem=True) == "srun --exact -n $1 ${STEP_NODES:+-N $STEP_NODES} ${STEP_MEM:+--mem=$STEP_MEM} --cpu-bind=cores"
    assert cp2k_launch({"launcher": "mpirun"}, "$1", step_mem=True) == "mpirun -np $1"


def test_reaction_job_memory_per_node():
    nhr = load_cluster("nhr")
    plan = {"HCl": {"n_bf": 40, "runs": ["charges"], "small": True},
            "AsH3": {"n_bf": 60, "runs": ["charges"], "small": True}}
    script = rx.job(Path("rx"), nhr, plan)
    assert "#SBATCH --mem=" in script and "--mem-per-cpu" not in script
    assert "${STEP_MEM:+--mem=$STEP_MEM}" in script and 'STEP_MEM=$(step_mem "$n")' in script
