#!/bin/bash
#SBATCH -J {{JOB_NAME}}
#SBATCH -t {{TIME}}
#SBATCH -p {{PARTITION}}
#SBATCH --qos={{QOS}}
#SBATCH --gres={{GRES}}
{{CONSTRAINT}}#SBATCH -c {{CPUS}}
#SBATCH --mem={{MEM}}
{{SETUP}}
{{CONDA_INIT}}
conda activate {{CONDA_ENV}}
# torch's pip wheels otherwise load the system libstdc++, which numpy cannot use.
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
export ORCHESTRAI_SINGLE_ENV=1
export PYTORCH_ALLOC_CONF=expandable_segments:True   # less fragmentation for varying batch sizes
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export QDPROPS_MACE_MODEL={{MODEL}}
export QDPROPS_XTB=$CONDA_PREFIX/bin/xtb
export QDPROPS_GXTB={{GXTB}}
export QDPROPS_CIF_DIRS={{CIF_DIRS}}
export QDPROPS_REFS={{REFS}}
{{EXTRA_ENV}}

nvidia-smi -L
{{RUN}}
