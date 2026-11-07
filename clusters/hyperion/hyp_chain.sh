#!/bin/bash
#SBATCH -p general
#SBATCH --qos=regular
#SBATCH -n 240
#SBATCH -c 1
#SBATCH -t 24:00:00
#SBATCH -o hyp-%j.out
# One 24 h piece of a campaign CP2K chain on Hyperion (qdw_campaign/records/<id>/cp2k, chain.sh from
# `qdw props chains --cluster hyperion_campaign`). Submitted by hyp_submit.sh, which writes the chain's job shape to
# .hyp_shape (e.g. "-n 240 --mem-per-cpu=4G"): ranks spread over as many nodes as Slurm finds free cores on (Ivan: a
# fixed node count never starts on Hyperion), all of one CPU type: with no --constraint Hyperion's submit plugin sets
# [cascadelake|icelake|emerald|zen3] (all nodes of one feature, so one CP2K build per job; general has no Zen3 nodes,
# which have no CP2K 2026.2). Every CP2K run uses all ranks (mpirun), one after
# the other.
# The piece first queues its successor (afterany, same shape): a piece stopped by the time limit or a node failure is
# continued by it, and chain.sh skips every finished step. When chain.sh ends (done or failed) the successor is cancelled.
# TREXIO or QDEX failing (their memory sits on few ranks / the first node) get one piece of their own on one node:
# TREXIO on an emerald node with 1800G, QDEX on 48 cores with 240G; the chain goes on from there.
# At most MAX_PIECES pieces per chain; the outcome of each piece is a line in logs/hyp_chains.tsv.
C=/scratch/iinfante/qdw_campaign
MAX_PIECES=8
cd "$SLURM_SUBMIT_DIR" || exit 1
ID=$(basename "$(dirname "$PWD")")
SHAPE=$(cat .hyp_shape 2>/dev/null || echo "-n 240 --mem-per-cpu=4G")
row() { printf '%s\t%s\t%s\t%s\t%s\n' "$(date -Is)" "$ID" "$SLURM_JOB_ID" "$1" "$2" >> $C/logs/hyp_chains.tsv; }
if [ -f .done ]; then row 0 "already done"; exit 0; fi
n=$(( $(cat .hyp_pieces 2>/dev/null || echo 0) + 1 )); echo $n > .hyp_pieces
if [ $n -gt $MAX_PIECES ]; then row $n "STOPPED: more than $MAX_PIECES pieces"; exit 1; fi
# a TREXIO run that did not end (killed with its job) leaves an orbitals file the chain would take as finished
if [ -f trexio/cp2k.out ] && ! grep -q "PROGRAM ENDED" trexio/cp2k.out; then
    for f in trexio/orbitals*.h5 trexio/orbitals*.trexio; do [ -e "$f" ] && mv "$f" "$f.partial-before-$SLURM_JOB_ID" && row $n "renamed unfinished $f"; done
fi
next=$(sbatch --parsable --dependency=afterany:$SLURM_JOB_ID -J "$SLURM_JOB_NAME" $SHAPE $C/bin/hyp_chain.sh)
row $n "start on $SLURM_JOB_NODELIST ($SLURM_NTASKS ranks on $SLURM_JOB_NUM_NODES nodes: $SLURM_TASKS_PER_NODE), successor ${next:-NOT QUEUED}"
export QDW_LOG=$C/logs/campaign.log
bash chain.sh > chain-$SLURM_JOB_ID.log 2>&1; rc=$?
[ -n "$next" ] && scancel "$next"
if [ $rc -eq 0 ] && [ -f .done ]; then row $n "DONE"; exit 0; fi
last=$(grep -o '^\[chain\] [a-z_]*' chain-$SLURM_JOB_ID.log | tail -1); last=${last#\[chain\] }
one_node() {   # $1 step, $2.. sbatch options: one piece of its own for that step, once per chain
    local step=$1; shift
    [ -f .$step.alone_tried ] && return 1
    date -Is > .$step.alone_tried
    j=$(sbatch --parsable -J "$SLURM_JOB_NAME" "$@" $C/bin/hyp_chain.sh)
    row $n "$step failed (rc=$rc, read $PWD/chain-$SLURM_JOB_ID.log): $step piece on one node ($*), job ${j:-NOT QUEUED}"
}
case $last in
    trexio)
        for f in trexio/orbitals*.h5 trexio/orbitals*.trexio; do [ -e "$f" ] && mv "$f" "$f.partial-$SLURM_JOB_ID"; done
        one_node trexio --constraint=emerald -N 1 -n 112 --mem=1800G -t 08:00:00 && exit $rc ;;
    qdex|qdex_trap)
        one_node qdex -N 1 -n 48 --mem=240G -t 12:00:00 && exit $rc ;;
esac
row $n "FAILED rc=$rc in ${last:-?} (read $PWD/chain-$SLURM_JOB_ID.log)"
exit $rc
