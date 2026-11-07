#!/bin/bash
#SBATCH -A {{ACCOUNT}}
#SBATCH -p {{PARTITION}}
#SBATCH -t {{TIME}}
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 32
#SBATCH --mem=120G
#SBATCH -o steps-%j.out
# After a relaxation wave is collected (collect, gcache): layer steps of every edges file, dG, ladders, then the
# network export into OUT. CPU only.
# Fill in: ACCOUNT (project_465003437), PARTITION (small), TIME (04:00:00), WORKDIR (with net/, tree/ =
# results_by_key.json, g_cache.json and the island batches, islands/ = states_index.json, thermo_merged.json),
# ENV_SH, QDW_SRC, OAI_SRC, PARAMS (templates/growth/cdse_zb.yaml), OUT.
D={{WORKDIR}}
cd $D
source {{ENV_SH}}
srun singularity exec $SIF bash -c "\$WITH_CONDA; source $W/venv/bin/activate; \
  export PYTHONPATH={{QDW_SRC}}:{{OAI_SRC}} OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1; set -e; \
  LP='python -m qd_workflows.layer_path'; \
  \$LP layer-steps net tree/results_by_key.json tree/g_cache.json tree/layer_steps.json --nproc 32; \
  \$LP attach-dg net tree/results_by_key.json tree/g_cache.json tree/layer_steps.json; \
  \$LP ladders tree/layer_steps.json tree/ladders.json | tail -3; \
  mkdir -p {{OUT}}; \
  python -m qd_workflows.network_export --registry net/registry.json --results tree/results_by_key.json \
    --g-cache tree/g_cache.json --steps tree/layer_steps.json --ladders tree/ladders.json \
    --islands islands/states_index.json --island-states tree/isl_res/isl/batch*/states.json \
    --thermo thermo_merged.json --params {{PARAMS}} --out {{OUT}} --nproc 32 | grep -v 'thermal [0-9]'"
