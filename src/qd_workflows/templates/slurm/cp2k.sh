#!/bin/bash
#SBATCH -J {{JOB_NAME}}
#SBATCH -t {{TIME}}
#SBATCH -n {{NTASKS}}
#SBATCH -N {{NODES}}
#SBATCH --mem-per-cpu={{MEM_PER_CPU}}
#SBATCH --qos={{QOS}}
{{PARTITION}}{{SETUP}}
{{CP2K_SETUP}}

echo "${SLURM_NTASKS} processes"

export OMP_NUM_THREADS={{OMP_NUM_THREADS}}

{{LAUNCH}} {{EXECUTABLE}} -i {{INPUT}} -o {{OUTPUT}}
