#!/bin/bash
#SBATCH -p general
#SBATCH --qos=regular
#SBATCH --constraint=emerald
#SBATCH -N 1
#SBATCH --ntasks-per-node=112
#SBATCH --cpus-per-task=1
#SBATCH -t 24:00:00
#SBATCH -o hyp-%j.out
# One 24 h piece of a campaign CP2K chain on Hyperion (qdw_campaign/records/<id>/cp2k, chain.sh from
# `qdw props chains --cluster hyperion_campaign`): one emerald node (2 TB), all 112 ranks on every CP2K run, one run after
# the other. One node, not two: on 2026-10-21 Slurm planned a 1-node job 35 h before a 2-node one (emerald nodes are
# shared with many small multi-day jobs). Memory is sized for geo_opt (~800 GB at 2,000 atoms, ~ n^2, x 1.25), not for
# TREXIO (Ivan: one SCF cycle from the converged .wfn). Submit from the chain directory:
# sbatch -J hyp-<id> --mem=<max(400, 1000 (n/2000)^2) G> hyp_chain.sh
# TREXIO that fails in that memory gets one short piece of its own with TREXIO_MEM on an emerald node (the chain goes on
# from there: trim, QDEX); its partial orbitals file is renamed *.partial-<job>, else the chain would take it as finished.
# The piece first queues its successor (afterany): a piece stopped by the time limit or a node failure is continued by
# it, and chain.sh skips every finished step. When chain.sh ends (done or failed) the successor is cancelled.
# At most MAX_PIECES pieces per chain; the outcome of each piece is a line in logs/hyp_chains.tsv.
C=/scratch/iinfante/qdw_campaign
MAX_PIECES=6
TREXIO_MEM=1800G
cd "$SLURM_SUBMIT_DIR" || exit 1
ID=$(basename "$(dirname "$PWD")")
row() { printf '%s\t%s\t%s\t%s\t%s\n' "$(date -Is)" "$ID" "$SLURM_JOB_ID" "$1" "$2" >> $C/logs/hyp_chains.tsv; }
if [ -f .done ]; then row 0 "already done"; exit 0; fi
n=$(( $(cat .hyp_pieces 2>/dev/null || echo 0) + 1 )); echo $n > .hyp_pieces
if [ $n -gt $MAX_PIECES ]; then row $n "STOPPED: more than $MAX_PIECES pieces"; exit 1; fi
# a TREXIO run that did not end (killed with its job) leaves an orbitals file the chain would take as finished
if [ -f trexio/cp2k.out ] && ! grep -q "PROGRAM ENDED" trexio/cp2k.out; then
    for f in trexio/orbitals*.h5 trexio/orbitals*.trexio; do [ -e "$f" ] && mv "$f" "$f.partial-before-$SLURM_JOB_ID" && row $n "renamed unfinished $f"; done
fi
next=$(sbatch --parsable --dependency=afterany:$SLURM_JOB_ID -J "$SLURM_JOB_NAME" --mem=${SLURM_MEM_PER_NODE}M $C/bin/hyp_chain.sh)
row $n "start on $SLURM_JOB_NODELIST ($SLURM_NTASKS ranks, ${SLURM_MEM_PER_NODE} MB per node), successor ${next:-NOT QUEUED}"
export QDW_LOG=$C/logs/campaign.log
bash chain.sh > chain-$SLURM_JOB_ID.log 2>&1; rc=$?
[ -n "$next" ] && scancel "$next"
if [ $rc -eq 0 ] && [ -f .done ]; then row $n "DONE"; exit 0; fi
last=$(grep -o '^\[chain\] [a-z_]*' chain-$SLURM_JOB_ID.log | tail -1)
if [ "$last" = "[chain] trexio" ] && [ ! -f .trexio_mem_tried ]; then
    date -Is > .trexio_mem_tried
    for f in trexio/orbitals*.h5 trexio/orbitals*.trexio; do [ -e "$f" ] && mv "$f" "$f.partial-$SLURM_JOB_ID"; done
    j=$(sbatch --parsable -J "$SLURM_JOB_NAME" --mem=$TREXIO_MEM -t 08:00:00 $C/bin/hyp_chain.sh)
    row $n "TREXIO failed in ${SLURM_MEM_PER_NODE} MB (rc=$rc, read $PWD/trexio/cp2k.out): TREXIO piece with $TREXIO_MEM, job ${j:-NOT QUEUED}"
    exit $rc
fi
row $n "FAILED rc=$rc in ${last#\[chain\] } (read $PWD/chain-$SLURM_JOB_ID.log)"
exit $rc
