#!/bin/bash
# Bring finished Hyperion chains back into the NHR records. Runs ON NHR (ssh hyperion with ~/.ssh/id_ed25519_hyperion,
# Host hyperion in ~/.ssh/config); deployed as qdw_campaign/logs/hyperion_sync.sh, so the NHR watchdog can call it.
# Usage: hyperion_sync.sh [id...]   (default: every Hyperion dot with cp2k/.done and no cp2k/.synced_to_nhr)
# Per dot (record path from logs/hyperion_dots.tsv): Hyperion records/<id>/cp2k -> <record>/cp2k without the .wfn
# restart files (~20 GB per dot; they stay on Hyperion), Hyperion's chain.sh (NHR keeps its own) and the wrapper's
# counters. The esp/MA, esp/MX links to the Hyperion monomer library are pointed to NHR's library, a file HYPERION says
# where the chain ran, and .done is written last: the watchdog's stage 3 starts from it. Then the held NHR chain job of
# the dot (job name cp2k-<id>, pending) is cancelled. Nothing else is deleted.
H=/scratch/iinfante/qdw_campaign
N=/scratch/projects/bep00157/qdw_campaign
export QDW_LOG=$N/logs/campaign.log PYTHONPATH=$N/code/qd-workflows/src:$N/code/Orchestr.AI/src
PY=/home/becivinf/miniforge3/envs/qdw-compute/bin/python
ids="$*"
[ -n "$ids" ] || ids=$(ssh hyperion "cd $H/records && for d in */; do d=\${d%/}; [ -f \$d/cp2k/.done ] && [ ! -f \$d/cp2k/.synced_to_nhr ] && echo \$d; done; true")
for id in $ids; do (
    set -e
    rec=$(awk -F'\t' -v i="$id" '$1 == i {print $3}' $N/logs/hyperion_dots.tsv)
    [ -n "$rec" ] || { echo "$id: not in logs/hyperion_dots.tsv, skipped"; exit 0; }
    ssh hyperion "test -f $H/records/$id/cp2k/.done" || { echo "$id: no .done on Hyperion, skipped"; exit 0; }
    mkdir -p $rec/cp2k
    rsync -a --exclude=.done --exclude='*.wfn' --exclude='*.wfn.*' --exclude=/chain.sh --exclude='/chain.sh.*' \
        --exclude=.hyp_pieces --exclude=.hyp_shape --exclude='.*.alone_tried' --exclude=.synced_to_nhr \
        "hyperion:$H/records/$id/cp2k/" "$rec/cp2k/"
    cd $rec/cp2k
    for l in $(find esp -maxdepth 1 -type l 2>/dev/null); do
        t=$(readlink $l)
        case $t in $H/monomers/*) ln -sfn $N/monomers/${t#$H/monomers/} $l ;; esac
    done
    miss=""; for f in relaxed.xyz window.json orbitals.h5 qdex/webapp.json qdex_trap/webapp.json; do [ -s $f ] || miss="$miss $f"; done
    [ -z "$miss" ] || { echo "$id: missing$miss after the sync: .done not written"; exit 1; }
    jobs=$(ssh hyperion "awk -F'\t' -v i=$id '\$2 == i {print \$3}' $H/logs/hyp_chains.tsv | sort -u | tr '\n' ' '")
    cat > HYPERION <<EOT
Chain run on Hyperion (DIPC), CP2K/2026.2-foss-2025a cp2k.psmp, mpirun, all ranks of the job (1 thread each) on every
CP2K run, method $(cat METHOD 2>/dev/null); Slurm jobs $jobs(records/$id/cp2k/hyp-<job>.out, chain-<job>.log);
profile clusters/hyperion_campaign.yaml, wrapper clusters/hyperion/hyp_chain.sh (qd-workflows feat/campaign).
Synced to NHR $(date -u +%Y-%m-%dT%H:%MZ). The .wfn restart files are not here: $H/records/$id/cp2k on Hyperion.
EOT
    touch .done
    held=$(squeue -u becivinf -h -t PD -n cp2k-$id -o %i | tr '\n' ' ')
    [ -z "$held" ] || scancel $held
    $PY -m qd_workflows.cli log event chain "hyperion: $id chain done on Hyperion (jobs $jobs), synced to $rec/cp2k; held NHR job(s) cancelled: ${held:-none}" > /dev/null 2>&1 || true
    ssh hyperion "date -u +%Y-%m-%dT%H:%MZ > $H/records/$id/cp2k/.synced_to_nhr"
    echo "$id: synced to $rec/cp2k; NHR held job(s) cancelled: ${held:-none}"
) done
