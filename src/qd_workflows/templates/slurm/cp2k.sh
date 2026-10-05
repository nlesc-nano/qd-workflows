#!/bin/bash
#SBATCH -J {{JOB_NAME}}
#SBATCH -t {{TIME}}
#SBATCH -n {{NTASKS}}
#SBATCH --mem-per-cpu={{MEM_PER_CPU}}
#SBATCH --qos={{QOS}}

module purge
module load {{MODULE}}

echo "${SLURM_NTASKS} processes"

export OMP_NUM_THREADS={{OMP_NUM_THREADS}}

mpirun -np ${SLURM_NTASKS} {{EXECUTABLE}} -i {{INPUT}} -o {{OUTPUT}}
