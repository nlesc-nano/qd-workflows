#!/bin/bash
#SBATCH -J {{JOB_NAME}}
#SBATCH -t {{TIME}}
#SBATCH -p {{PARTITION}}
#SBATCH --qos={{QOS}}
#SBATCH --gres={{GRES}}
{{CONSTRAINT}}{{ACCOUNT}}#SBATCH -c {{CPUS}}
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
# Progress log of the run (qdw log, $QDW_LOG): start and end of the job, and each step as it ends
# with its key result (a follower reads this job's output).
export QDW_LOG=${QDW_LOG:-{{QDW_LOG}}}
qdwlog() { [ -z "$QDW_LOG" ] || env -u PYTHONHOME PYTHONPATH={{QDW_PYTHONPATH}} python -m qd_workflows.cli log "$@" > /dev/null 2>&1 || true; }
T0=$(date +%s); STOP=.follow_stop.$SLURM_JOB_ID; rm -f $STOP
OUT=$(scontrol show job $SLURM_JOB_ID 2>/dev/null | sed -n 's/^ *StdOut=//p'); [ -n "$OUT" ] || OUT=$SLURM_SUBMIT_DIR/slurm-$SLURM_JOB_ID.out
qdwlog event gpu "start: {{JOB_NAME}}, $(nvidia-smi --query-gpu=name --format=csv,noheader | head -1), limit {{TIME}}"
[ -z "$QDW_LOG" ] || [ -z "{{RECORDS}}" ] || { qdwlog follow "$OUT" --stage gpu --records {{RECORDS}} --stop $STOP & FOLLOW=$!; }
on_exit() {
    local rc=$?
    touch $STOP; [ -z "${FOLLOW:-}" ] || wait $FOLLOW; rm -f $STOP
    if [ $rc -eq 0 ]; then qdwlog event gpu "end: {{JOB_NAME}}" --since $T0
    else qdwlog event gpu "ERROR {{JOB_NAME}} ended with exit $rc (read $OUT)" --since $T0; fi
}
trap on_exit EXIT

{{RUN}}
