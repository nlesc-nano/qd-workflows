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
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export QDPROPS_MACE_MODEL={{MODEL}}
export QDPROPS_XTB=$CONDA_PREFIX/bin/xtb
export QDPROPS_GXTB={{GXTB}}
export QDPROPS_CIF_DIRS={{CIF_DIRS}}
export QDPROPS_REFS={{REFS}}

nvidia-smi -L
python -m orchestr_ai.postprocessing config.yaml
