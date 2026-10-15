#!/bin/bash
#SBATCH -J {{JOB_NAME}}
#SBATCH -t {{TIME}}
#SBATCH -n {{NTASKS}}
#SBATCH -c {{THREADS}}
#SBATCH -N {{NODES}}
{{MEM}}
#SBATCH --qos={{QOS}}
{{PARTITION}}{{ACCOUNT}}{{SETUP}}
# CP2K single points of a reaction run (qdw reaction): per species a GAPW single point at the
# given geometry (Hirshfeld charges, .wfn), then the GPW vacuum reference and one SCCS run per
# eps, restarted from it. Species run side by side.
# An SCF that has not converged in 200 cycles is started again from the .wfn it wrote (up to
# 8 runs, then .failed). Finished runs (.done) are skipped, so the job can be resubmitted.
# PIECES=1 (short queues): the job first queues its successor, which starts when this one
# ends however it ends; the last piece cancels it and submits the analysis.
{{CP2K_SETUP}}
export OMP_NUM_THREADS={{OMP_NUM_THREADS}} OMP_PLACES=cores OMP_PROC_BIND=close OMP_STACKSIZE=256M
ulimit -s unlimited
THREADS={{THREADS}}         # OpenMP threads per rank; launch gives cores, the chains run cores/THREADS ranks
cd "$SLURM_SUBMIT_DIR"
ROOT=$PWD
mkdir -p logs
PIECES={{PIECES}}
MAX_PIECES=40
next=""
piece=$(( $(cat .pieces 2>/dev/null || echo 0) + 1 )); echo $piece > .pieces
# Stop cleanly before the time limit: CP2K ends its SCF when it finds an EXIT file and writes its
# .wfn; such a run is put aside as cp2k.out.wall-* and does not count as an attempt.
rm -f .stop
left=$(squeue -h -j "$SLURM_JOB_ID" -o %L | awk -F'[-:]' '{s=0; for (i=1;i<=NF;i++) s=s*60+$i; if (NF==4) s=$1*86400+$2*3600+$3*60+$4; print s}')
if [ -n "$left" ] && [ "$left" -gt 300 ]; then
    ( sleep $(( left - 150 )); touch "$ROOT/.stop"
      for d in species/*/cp2k/*/; do [ -f "$d/cp2k.out" ] && [ ! -f "$d/.done" ] && touch "$d/EXIT"; done ) &
fi
if [ "$PIECES" = 1 ] && [ $piece -lt $MAX_PIECES ]; then
    next=$(sbatch --parsable --dependency=afterany:$SLURM_JOB_ID job_cp2k.sh) || next=""
    echo "[reaction] piece $piece, next: ${next:-not queued}"
fi
cp2k_n() { OMP_NUM_THREADS=$2 {{LAUNCH_N}} {{EXECUTABLE}} -i cp2k.inp -o cp2k.out; }   # $1 ranks, $2 threads each
# SCCS runs: one thread per rank. With OpenMP threads the SCCS polarisation loop of CP2K 2026.2
# goes wrong (HCl, eps 80: diverges with 4 threads, converges in 13 steps with 1; AsH3, eps 2-4:
# +3.8 and +8.3 eV with 4 threads). Other runs: THREADS per rank.
layout() { case "$1" in */sccs_e*) echo "$2 1" ;; *) echo "$(( $2 / THREADS )) $THREADS" ;; esac; }   # dir cores -> ranks threads

scf_ok() {   # the run ended and its last SCF converged
    grep -q "PROGRAM ENDED" cp2k.out 2>/dev/null || return 1
    [ "$(grep -E 'SCF run (NOT )?converged' cp2k.out | tail -1 | grep -c NOT)" = 0 ] && grep -q "SCF run converged" cp2k.out
}

sccs_ok() {   # an SCCS run really solvated: the cavity is printed once the SCCS cycle runs
    case "$1" in */sccs_e*) grep -q "SCCS| Volume of the solute cavity" cp2k.out ;; *) return 0 ;; esac
}

# An SCCS run whose SCF converged before the solvent was switched on is not finished: again.
for d in species/*/cp2k/sccs_e*/; do
    [ -f "$d/.done" ] || continue
    if ! grep -q "SCCS| Volume of the solute cavity" "$d/cp2k.out" 2>/dev/null; then
        rm -f "$d/.done" "$d"/*-RESTART.wfn; mv "$d/cp2k.out" "$d/cp2k.out.nosccs" 2>/dev/null
        echo "[reaction] $d: SCCS never ran, again"
    fi
done

run_one() {  # dir cores [parent .wfn]
    local d=$1 parent=$3 p k n th
    read -r n th <<< "$(layout "$1" "$2")"
    [ "$n" -ge 1 ] || n=1
    [ -d "$d" ] || return 0
    [ -f "$d/.done" ] || [ -f "$d/.failed" ] && return 0
    cd "$d" || return 1
    p=$(awk '$1 == "PROJECT" {print $2; exit}' cp2k.inp)
    if [ ! -s "$p-RESTART.wfn" ] && [ -n "$parent" ] && [ -s "$parent" ]; then cp "$parent" "$p-RESTART.wfn"; fi
    while true; do
        [ -f "$ROOT/.stop" ] && break       # the piece is ending: the next one continues
        k=$(ls cp2k.out.prev-* 2>/dev/null | wc -l)
        if [ "$k" -ge 8 ]; then touch .failed; echo "[reaction] $d: no converged SCF in $k runs, given up"; break; fi
        [ -s cp2k.out ] && mv cp2k.out "cp2k.out.prev-$((k + 1))"
        cp2k_n "$n" "$th"
        if [ -f "$ROOT/.stop" ]; then        # stopped for the time limit (CP2K deletes EXIT): not a failed attempt
            rm -f EXIT; mv cp2k.out "cp2k.out.wall-$(date +%s)"; echo "[reaction] $d: stopped before the time limit"; break
        fi
        if scf_ok && sccs_ok "$d"; then touch .done; echo "[reaction] $d done"; break; fi
        if scf_ok; then         # converged before SCCS started: again from the single point's .wfn ($w of the chain)
            local from=${w:-$parent}
            rm -f "$p-RESTART.wfn"; [ -n "$from" ] && [ -s "$from" ] && cp "$from" "$p-RESTART.wfn"
            echo "[reaction] $d: SCF converged before SCCS started, again from $from"; continue
        fi
        echo "[reaction] $d: SCF not converged, restarting from $p-RESTART.wfn"
    done
    cd "$ROOT"
    [ -f "$d/.done" ]
}

species_chain() {  # name cores
    local s=species/$1/cp2k n=$2 e
    run_one "$s/charges" "$n"
    local w="$ROOT/$s/charges/$1-charges-RESTART.wfn"
    for e in "$s"/energy_*; do run_one "$e" "$n" "$w"; done      # isolated single points (energy check)
    [ -d "$s/sccs_vac" ] || return 0
    run_one "$s/sccs_vac" "$n" "$w" || echo "[reaction] $1: vacuum reference failed"
    # a ladder in eps: the lowest from the single point's .wfn (not the converged vacuum one:
    # restarted from that, the SCF converges in its first cycle, before the SCCS cycle starts),
    # each next eps from the one below, so the dielectric is switched on in small steps
    local prev=$w
    for e in $(ls -d "$s"/sccs_e* | awk -F'sccs_e' '{print $NF, $0}' | sort -g | cut -d' ' -f2); do
        run_one "$e" "$n" "$prev"
        prev="$ROOT/$e/$1-e${e##*sccs_e}-RESTART.wfn"
    done
}

TOTAL=$(( SLURM_NTASKS * THREADS ))      # cores
used=0
declare -A RANKS
reap() { local p; for p in "${!RANKS[@]}"; do kill -0 "$p" 2>/dev/null || { used=$((used - RANKS[$p])); unset "RANKS[$p]"; }; done; }
launch() {  # name ranks: start when enough ranks are free (never more than the job has)
    local n=$2
    [ "$n" -gt "$TOTAL" ] && n=$TOTAL
    while [ $((used + n)) -gt "$TOTAL" ]; do wait -n || true; reap; done
    species_chain "$1" "$n" > "logs/cp2k_$1.log" 2>&1 &
    RANKS[$!]=$n
    used=$((used + n))
    echo "[reaction] $1 on $n ranks"
}

{{SPECIES}}
wait
all_done() { local d; for d in species/*/cp2k/*/; do [ -f "$d/.done" ] || [ -f "$d/.failed" ] || return 1; done; }
if all_done; then
    echo "[reaction] CP2K done ($(ls -d species/*/cp2k/*/.failed 2>/dev/null | wc -l) runs failed)"
    [ -n "$next" ] && scancel "$next" && echo "[reaction] cancelled $next"
    [ "$PIECES" = 1 ] && sbatch job_analyse.sh
fi
