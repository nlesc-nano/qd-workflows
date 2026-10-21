#!/bin/bash
# Copy campaign dots from NHR to Hyperion and prepare their CP2K chains there (run on the laptop: NHR and Hyperion
# do not see each other). Usage: migrate.sh <ledger record path on NHR>...
# Per dot: start.xyz, record.json and the props files the chain reads (relaxed.xyz, detach_*.xyz, adsorb_*.xyz,
# monomers.xyz) go to qdw_campaign/records/<id>; the chain is written there by `qdw props chains --cluster
# hyperion_campaign`. NHR progress (a stopped geo_opt; finished field/esp/pdos steps) is copied to records/<id>/nhr_cp2k
# and laid over the new chain, its trajectory newer than geom.xyz so geo_opt continues from the last geometry. Never copied: chain
# scripts, slurm outputs, TREXIO_ON_LARGE and trexio/ (a run killed by OOM leaves an orbitals.h5 the chain would take
# as finished). Nothing is deleted on either side.
RELAY=${RELAY:-$HOME/.cache/qdw_hyperion_relay}
H=/scratch/iinfante/qdw_campaign
mkdir -p "$RELAY"
one() (   # one dot; all of it again after a failure (Hyperion closes ssh connections now and then)
    set -e
    rec=$1; id=$(basename "$rec")
    mkdir -p "$RELAY/$id"
    rsync -a -e "ssh -x" --include=/start.xyz --include=/record.json --include=/props/ --include=/props/relaxed.xyz \
        --include='/props/detach_*.xyz' --include='/props/adsorb_*.xyz' --include=/props/monomers.xyz --exclude='*' \
        "nhr:$rec/" "$RELAY/$id/"
    # NHR chain progress: geo_opt trajectory, output and .wfn; finished later steps
    if ssh -x nhr "test -d $rec/cp2k"; then rsync -a -e "ssh -x" --exclude=trexio/ --exclude=TREXIO_ON_LARGE --exclude='chain*.sh' --exclude='slurm-*.out' \
        --exclude='chain-*.log' --exclude='*.bak*' --exclude=record.json --exclude=geo_opt/job.sh --exclude=geo_opt/job.json \
        --exclude=geo_opt/geom.xyz --exclude=geo_opt/cp2k.inp --exclude='geo_opt/cp2k.inp.*' \
        "nhr:$rec/cp2k/" "$RELAY/$id/nhr_cp2k/"; fi
    # an NHR cp2k folder with only its prepared geo_opt input (nothing run) carries no progress
    if [ -d "$RELAY/$id/nhr_cp2k" ] && [ -z "$(find "$RELAY/$id/nhr_cp2k" -type f | head -1)" ]; then rm -r "$RELAY/$id/nhr_cp2k"; fi
    rsync -a -e "ssh -x" "$RELAY/$id/" "hyperion:$H/records/$id/"
    ssh -x hyperion "bash -s" <<EOF
set -e
cd $H/records/$id
export QDW_LOG=$H/logs/campaign.log
[ -e cp2k/chain.sh ] && { echo "$id: cp2k/chain.sh exists, left as it is"; exit 0; }
# NHR progress first (NHR times kept: the sync back skips the files that did not change): its METHOD (gapw or gpw)
# then holds for the chain written here
[ -d nhr_cp2k ] && mkdir -p cp2k && cp -rp nhr_cp2k/. cp2k/
env -u PYTHONHOME PYTHONPATH=$H/code/qd-workflows/src:$H/code/Orchestr.AI/src /dipc/iinfante/.conda/envs/qdw-compute/bin/python \
    -m qd_workflows.cli props chains . --cluster hyperion_campaign > /dev/null
[ -s cp2k/chain.sh ] || { echo "$id: no chain written (no props/relaxed.xyz?)"; exit 1; }
touch cp2k/geo_opt/*-pos-1.xyz* 2>/dev/null || true   # newer than the new geom.xyz: geo_opt goes on from it
echo "$id: chain (\$(cat cp2k/METHOD))\$([ -d nhr_cp2k ] && echo " + NHR progress: \$(cd nhr_cp2k && ls | tr '\n' ' ')")"
EOF
)
for rec in "$@"; do
    for try in 1 2 3 4 5; do one "$rec" && break; echo "$(basename "$rec"): try $try failed, again in 60 s" >&2; sleep 60; done
done
