#!/bin/bash
# migrate.sh for many dots at once: one file list on NHR, one rsync NHR -> laptop, one rsync laptop -> Hyperion, one ssh
# to write the chains (migrate.sh needs ~6 ssh connections per dot, ~4 min each over the laptop link).
# Usage: bulk_migrate.sh <file with one dot id per line>   (ids of NHR qdw_campaign/logs/hyperion_dots.tsv)
# Copied per dot: start.xyz, record.json, props/{relaxed.xyz, detach_*.xyz, adsorb_*.xyz, monomers.xyz}; NHR chain
# progress: cp2k/METHOD and a stopped geo_opt's cp2k.out, trajectory and .wfn (no *.bak*, cubes, ensembles, other props).
# Dots whose Hyperion cp2k/chain.sh exists are left as they are. Nothing is deleted on either side.
set -e
ids=$1
RELAY=${RELAY:-$HOME/.cache/qdw_hyperion_relay}/bulk
H=/scratch/iinfante/qdw_campaign
N=/scratch/projects/bep00157/qdw_campaign
mkdir -p "$RELAY/nhr" "$RELAY/staged"
scp -q -o ForwardX11=no "$ids" nhr:$N/logs/hyp_bulk_ids.txt     # shared scratch: login nodes differ between connections
ssh -x nhr "python3 - $N/logs/hyp_bulk_ids.txt $N/logs/hyperion_dots.tsv" > "$RELAY/files.txt" <<'EOF'
import glob, os, sys
want = set(open(sys.argv[1]).read().split())
for line in open(sys.argv[2]):
    i, n, rec = line.split("\t")[:3]
    if i not in want:
        continue
    pats = ["start.xyz", "record.json", "props/relaxed.xyz", "props/detach_*.xyz", "props/adsorb_*.xyz",
            "props/monomers.xyz", "cp2k/METHOD", "cp2k/geo_opt/cp2k.out", "cp2k/geo_opt/*-pos-1.xyz*",
            "cp2k/geo_opt/*-RESTART.wfn"]
    for p in pats:
        for f in sorted(glob.glob(f"{rec}/{p}")):
            if ".bak" not in f:
                print(f.lstrip("/"))
EOF
echo "$(wc -l < "$RELAY/files.txt") files"
rsync -a -e "ssh -x" --files-from="$RELAY/files.txt" nhr:/ "$RELAY/nhr/"
# <record>/x -> staged/<id>/x, <record>/cp2k/x -> staged/<id>/nhr_cp2k/x (hard links)
python3 - "$RELAY" <<'EOF'
import os, sys
r = sys.argv[1]
for f in open(f"{r}/files.txt").read().split():
    parts = f.split("/")
    k = next(j for j, p in enumerate(parts) if p in ("start.xyz", "record.json", "props", "cp2k"))
    rid, rest = parts[k - 1], parts[k:]
    if rest[0] == "cp2k":
        rest = ["nhr_cp2k"] + rest[1:]
    dst = os.path.join(r, "staged", rid, *rest)
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    if not os.path.exists(dst):
        os.link(os.path.join(r, "nhr", f), dst)
EOF
rsync -a -e "ssh -x" "$RELAY/staged/" "hyperion:$H/records/"
ssh -x hyperion "bash -s" < <(echo "IDS='$(tr '\n' ' ' < "$ids")'"; cat <<EOF
cd $H/records
export QDW_LOG=$H/logs/campaign.log
for id in \$IDS; do (
    cd \$id || exit 0
    [ -e cp2k/chain.sh ] && { echo "\$id: cp2k/chain.sh exists, left as it is"; exit 0; }
    [ -d nhr_cp2k ] && mkdir -p cp2k && cp -rp nhr_cp2k/. cp2k/
    env -u PYTHONHOME PYTHONPATH=$H/code/qd-workflows/src:$H/code/Orchestr.AI/src /dipc/iinfante/.conda/envs/qdw-compute/bin/python \
        -m qd_workflows.cli props chains . --cluster hyperion_campaign > /dev/null
    [ -s cp2k/chain.sh ] || { echo "\$id: NO CHAIN written"; exit 0; }
    touch cp2k/geo_opt/*-pos-1.xyz* 2>/dev/null || true
    echo "\$id: chain (\$(cat cp2k/METHOD))\$(ls nhr_cp2k/geo_opt/*-pos-1.xyz >/dev/null 2>&1 && echo ' + NHR geo_opt progress')"
) done
EOF
)
