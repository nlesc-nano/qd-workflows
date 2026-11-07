#!/bin/bash
#SBATCH -A {{ACCOUNT}}
#SBATCH -p {{PARTITION}}
#SBATCH -t {{TIME}}
#SBATCH -N 1
#SBATCH --gpus-per-node=1
#SBATCH --ntasks-per-node=1
#SBATCH -c 7
#SBATCH --mem=60G
#SBATCH -o relax-%x-%j.out
# Layer-path states (registry dots d<id>.xyz, islands i<n>.xyz): relax (no anneal) + polish + FD Hessian,
# MACE-MH-1 omat_pbe float64, the (m, k) grid's relaxer (python -m qd_workflows.layer_path relax).
# Fill in: ACCOUNT (project_465003437), PARTITION (small-g; dev-g for a short test), TIME (05:00:00), WORKDIR (the
# folder with the waves), ENV_SH (sets SIF, W = the MACE venv over the PyTorch ROCm container, WITH_CONDA, and M =
# the MACE model), QDW_SRC (qd-workflows src/ on the cluster), OAI_SRC (Orchestr.AI src/), SXL_JSON (the
# Symmetrix-XL export of the same model, symmetrix_extract_mace --head omat_pbe; the venv needs
# symmetrix-xl-rocm6-gfx90a). QDW_SXL=1 relaxes with Symmetrix-XL (same MACE numbers, ~3x faster per GCD);
# set it to 0 for mace-torch.
#   sbatch -J w1b0 layer_path_relax_lumi.sh wave1/batch0
cd {{WORKDIR}}
source {{ENV_SH}}
B=$1
srun singularity exec $SIF bash -c "\$WITH_CONDA; source $W/venv/bin/activate; source $W/env.sh; \
  export OMP_NUM_THREADS=7 PYTHONPATH={{QDW_SRC}}:{{OAI_SRC}} QDW_SXL=1 QDW_SXL_JSON={{SXL_JSON}}; \
  python -u -m qd_workflows.layer_path relax $B \$M $B/states.json"
