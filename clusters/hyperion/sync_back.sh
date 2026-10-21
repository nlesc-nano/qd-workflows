#!/bin/bash
# Bring finished Hyperion chains back into the NHR records (run on the laptop: the two clusters do not see each other).
# Usage: sync_back.sh [id...]   (default: every Hyperion dot with cp2k/.done and no cp2k/.synced_to_nhr)
# Per dot: Hyperion records/<id>/cp2k -> laptop relay -> NHR <record>/cp2k (record path from NHR logs/hyperion_dots.tsv).
# Not copied: the .wfn restart files (~20 GB per dot; they stay on Hyperion), Hyperion's chain.sh (NHR keeps its own) and
# the wrapper's counters. The esp/MA, esp/MX links to the Hyperion monomer library are pointed to NHR's library, a file
# HYPERION says where the chain ran, and .done is written last: the NHR watchdog's stage 3 starts from it. Then the
# held NHR chain job of the dot (job name cp2k-<id>, pending) is cancelled. Nothing else is deleted.
set -e
RELAY=${RELAY:-$HOME/.cache/qdw_hyperion_relay/back}
H=/scratch/iinfante/qdw_campaign
N=/scratch/projects/bep00157/qdw_campaign
ids="$*"
[ -n "$ids" ] || ids=$(ssh -x hyperion "cd $H/records && for d in */; do d=\${d%/}; [ -f \$d/cp2k/.done ] && [ ! -f \$d/cp2k/.synced_to_nhr ] && echo \$d; done; true")
for id in $ids; do
    rec=$(ssh -x nhr "awk -F'\t' '\$1 == \"$id\" {print \$3}' $N/logs/hyperion_dots.tsv")
    [ -n "$rec" ] || { echo "$id: not in NHR logs/hyperion_dots.tsv, skipped"; continue; }
    ssh -x hyperion "test -f $H/records/$id/cp2k/.done" || { echo "$id: no .done on Hyperion, skipped"; continue; }
    mkdir -p "$RELAY/$id"
    rsync -a -e "ssh -x" --exclude=.done --exclude='*.wfn' --exclude='*.wfn.*' --exclude=/chain.sh --exclude='/chain.sh.*' \
        --exclude=.hyp_pieces --exclude=.synced_to_nhr "hyperion:$H/records/$id/cp2k/" "$RELAY/$id/"
    rsync -a -e "ssh -x" "$RELAY/$id/" "nhr:$rec/cp2k/"
    jobs=$(ssh -x hyperion "awk -F'\t' '\$2 == \"$id\" {print \$3}' $H/logs/hyp_chains.tsv | sort -u | tr '\n' ' '")
    ssh -x nhr "bash -s" <<EOF
set -e
cd $rec/cp2k
for l in \$(find esp -maxdepth 1 -type l 2>/dev/null); do
    t=\$(readlink \$l)
    case \$t in $H/monomers/*) ln -sfn $N/monomers/\${t#$H/monomers/} \$l ;; esac
done
miss=""; for f in relaxed.xyz window.json orbitals.h5 qdex/webapp.json qdex_trap/webapp.json; do [ -s \$f ] || miss="\$miss \$f"; done
[ -z "\$miss" ] || { echo "$id: missing\$miss on NHR after the sync: .done not written"; exit 1; }
cat > HYPERION <<EOT
Chain run on Hyperion (DIPC), emerald node(s), CP2K/2026.2-foss-2025a cp2k.psmp, mpirun, all cores as MPI ranks, 1 thread,
method \$(cat METHOD 2>/dev/null), Slurm jobs $jobs; profile clusters/hyperion_campaign.yaml.
Synced to NHR $(date -u +%Y-%m-%dT%H:%MZ). The .wfn restart files are not here: $H/records/$id/cp2k on Hyperion.
EOT
touch .done
held=\$(squeue -u becivinf -h -t PD -n cp2k-$id -o %i | tr '\n' ' ')
[ -z "\$held" ] || scancel \$held
QDW_LOG=$N/logs/campaign.log PYTHONPATH=$N/code/qd-workflows/src:$N/code/Orchestr.AI/src /home/becivinf/miniforge3/envs/qdw-compute/bin/python \
    -m qd_workflows.cli log event chain "hyperion: $id chain done on Hyperion (jobs $jobs), synced to $rec/cp2k; held NHR job(s) cancelled: \${held:-none}" || true
echo "$id: synced to $rec/cp2k; NHR held job(s) cancelled: \${held:-none}"
EOF
    ssh -x hyperion "date -u +%Y-%m-%dT%H:%MZ > $H/records/$id/cp2k/.synced_to_nhr"
done
