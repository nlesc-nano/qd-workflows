#!/bin/bash
#SBATCH -J {{JOB_NAME}}
#SBATCH -t {{TIME}}
#SBATCH -n {{NTASKS}}
#SBATCH -N {{NODES}}
#SBATCH --ntasks-per-node={{NODE_TASKS}}
#SBATCH --mem={{MEM}}
#SBATCH --qos={{QOS}}
{{PARTITION}}{{ACCOUNT}}{{SETUP}}
# CP2K track of one structure in one job: geo_opt -> field single points -> PDOS -> MO window
# -> TREXIO -> trim -> QDEX, all transitions and without the surface traps (when record.json is here).
# Every step is skipped when its output is already there, so a job stopped by the
# time limit continues where it stopped when resubmitted.
{{CP2K_SETUP}}
export OMP_NUM_THREADS={{OMP_NUM_THREADS}}
set -e

P={{PROJECT}}
QDW_PY={{QDW_PYTHON}}
# qdw runs in its own environment: the CP2K module's PYTHONPATH (its own numpy) and
# libraries must not leak into it, nor its libraries into CP2K.
qdw() { env -u PYTHONHOME PYTHONPATH={{QDW_PYTHONPATH}} LD_LIBRARY_PATH=$(dirname "$(dirname "$QDW_PY")")/lib "$QDW_PY" -m qd_workflows.cli "$@"; }
# Orchestr.AI writes and checks the field inputs (same environment as qdw).
orch_field() { env -u PYTHONHOME PYTHONPATH={{QDW_PYTHONPATH}} LD_LIBRARY_PATH=$(dirname "$(dirname "$QDW_PY")")/lib "$QDW_PY" -m orchestr_ai.qd.field "$@"; }
# QDEX: its own environment too, on the cores of the first node.
QDEX={{QDEX}}
qdex_env() { env -u PYTHONPATH -u PYTHONHOME LD_LIBRARY_PATH=$(dirname "$(dirname "$QDEX")")/lib "$@"; }
qdex_run() { qdex_env OMP_NUM_THREADS=$1 "$QDEX" --config config.yaml; }
cp2k() { {{LAUNCH}} {{EXECUTABLE}} -i cp2k.inp -o cp2k.out; }
# $1 ranks, $STEP_MEM per node (see step_mem), on as few nodes as the ranks need: a run that fits on one
# node stays on one (spread over two it holds memory on both and talks across them); srun's and CP2K's
# messages in srun.log next to the run
PPN=${SLURM_NTASKS_PER_NODE:-{{NODE_TASKS}}}; PPN=${PPN%%(*}
cp2k_n() { local STEP_NODES=$(( ($1 + PPN - 1) / PPN )); {{LAUNCH_N}} {{EXECUTABLE}} -i cp2k.inp -o cp2k.out > srun.log 2>&1; }
# Ranks of each of $1 runs side by side: on several nodes, as many as let ceil($1 / nodes) runs share a
# node (13 runs on 2 x 96 cores: 7 per node, 13 ranks each); else an equal share of the job's ranks.
ranks_each() {
    local n=$1 nn=${SLURM_NNODES:-1} r
    if [ "$nn" -gt 1 ] && [ "$n" -ge "$nn" ]; then r=$(( PPN / ((n + nn - 1) / nn) )); else r=$(( SLURM_NTASKS / n )); fi
    [ $r -ge 1 ] || r=1
    echo $r
}
# How many single points of the dot fit side by side, at most $1: each holds the cell's real-space
# grids (~{{GRID_MEM_GB}} GB, whatever its rank count); 80 % of the job's memory is used.
side_by_side() {
    local mem_mb n nn=${SLURM_NNODES:-1}
    if [ "${SLURM_MEM_PER_NODE:-0}" -gt 0 ]; then
        n=$(( SLURM_MEM_PER_NODE * 8 / 10 / {{GRID_MEM_MB}} ))   # per node; on several nodes one
        [ "$nn" -le 1 ] || n=$(( (n - 1) * nn ))                  # place less (a run across two, step_mem)
    else
        if [ -n "${SLURM_MEM_PER_CPU:-}" ]; then
            mem_mb=$(( SLURM_MEM_PER_CPU * SLURM_NTASKS ))
        else                                 # the whole node(s)
            mem_mb=$(( $(awk '/MemTotal/ {print int($2 / 1024)}' /proc/meminfo) * nn ))
        fi
        n=$(( mem_mb * 8 / 10 / {{GRID_MEM_MB}} ))
    fi
    [ $n -ge 1 ] || n=1; [ $n -le $1 ] || n=$1
    echo $n
}
# Memory per node of each of $1 runs of $2 ranks side by side (srun --mem): an equal share of 90 % of a
# node among the runs on it (a run of more ranks than a node holds counts on each of its nodes). The job
# asks for memory per node, and a step that asks for none takes all of it: the other runs would wait for
# it. Empty without a per-node figure from Slurm (then srun asks for none).
step_mem() {
    [ "${SLURM_MEM_PER_NODE:-0}" -gt 0 ] || return 0
    local nn=${SLURM_NNODES:-1} span=$(( (${2:-1} + PPN - 1) / PPN ))
    local k=$(( ($1 * span + nn - 1) / nn ))
    [ $k -ge 1 ] || k=1; [ $k -le $1 ] || k=$1
    echo $(( SLURM_MEM_PER_NODE * 9 / 10 / k ))M
}
# The last lines srun and CP2K printed for the unfinished runs $2... of step directory $1.
why() { local d=$1 r; shift; for r in "$@"; do echo "[chain] $d/$r: $(tail -n 2 "$d/$r/srun.log" 2>/dev/null | tr '\n' ' ')"; done; }
EXTRA="--cluster {{CLUSTER}} --time {{TIME}} --qos {{QOS}}"
MONOLIB="{{MONOMER_LIBRARY}}"   # CP2K monomer library (qdw monomers), empty: every chain runs its monomers
# Progress log of the run (qdw log, $QDW_LOG): a line when each step starts and ends, with its key
# result; on exit, how the job ended and the step it was in. Never stops the chain.
export QDW_LOG=${QDW_LOG:-{{QDW_LOG}}}
plog() { [ -z "$QDW_LOG" ] || qdw log "$@" > /dev/null 2>&1 || true; }
T0=$(date +%s); STEP=start; CHAIN_DIR=$(pwd); t=
on_exit() {
    local rc=$?
    if [ $rc -eq 0 ]; then plog event chain "chain done: $CHAIN_DIR" --since $T0
    else plog event chain "ERROR chain stopped in $STEP (exit $rc): read $CHAIN_DIR/slurm-$SLURM_JOB_ID.out" --since $T0; fi
}
trap on_exit EXIT
plog event chain "start $P: $SLURM_NTASKS ranks on ${SLURM_NNODES:-1} node(s), ${SLURM_MEM_PER_NODE:-?} MB per node, limit {{TIME}}"

# 1. geo_opt (from the MACE-relaxed structure)
cd geo_opt
if ! grep -q "GEOMETRY OPTIMIZATION COMPLETED" cp2k.out 2>/dev/null; then
    if [ -s "$P-pos-1.xyz" ]; then      # keep the trajectory of a stopped run
        n=$(ls "$P-pos-1.xyz.prev-"* 2>/dev/null | wc -l)
        mv cp2k.out "cp2k.out.prev-$n" 2>/dev/null || true
        mv "$P-pos-1.xyz" "$P-pos-1.xyz.prev-$n"
    fi
    last=$(ls -t "$P-pos-1.xyz.prev-"* 2>/dev/null | head -1 || true)
    if [ -n "$last" ] && [ "$last" -nt geom.xyz ]; then   # continue from its last geometry
        qdw cp2k relaxed "$last" -o geom.xyz
    fi
    echo "[chain] geo_opt"; STEP=geo_opt
    plog event chain "geo_opt started$([ -n "$last" ] && echo " (continuing from $last)")"; cp2k
    plog cp2k chain geo_opt .
fi
grep -q "GEOMETRY OPTIMIZATION COMPLETED" cp2k.out || { echo "[chain] geo_opt did not converge"; plog event chain "ERROR geo_opt did not converge (read $CHAIN_DIR/geo_opt/cp2k.out)"; exit 1; }
cd ..
[ -s relaxed.xyz ] || qdw cp2k relaxed "geo_opt/$P-pos-1.xyz" -o relaxed.xyz

# 2. Finite-field single points: 13 restarted from the geo_opt .wfn ({{FIELD_PARALLEL}} at a time, the
#    job's ranks shared among them) -> Born charges, Raman tensors, Hirshfeld charges of the dot
#    (Orchestr.AI field step). Finished runs are skipped.
if [ "{{FIELD}}" = "1" ]; then
    [ -s field/base/cp2k.inp ] || qdw cp2k prepare relaxed.xyz --step sample --project "$P-field" --wfn "../../geo_opt/$P-RESTART.wfn" --out field/base $EXTRA
    [ -s field/runs.txt ] || orch_field prepare field/base/cp2k.inp field
    todo=$(orch_field status field)
    if [ "$todo" != "all done" ]; then
        echo "[chain] field: $todo"
        FP=$(side_by_side {{FIELD_PARALLEL}})
        R=$(ranks_each $FP)
        STEP_MEM=$(step_mem $FP $R); echo "[chain] field: $FP side by side, $R ranks${STEP_MEM:+ and $STEP_MEM} each"
        STEP=field; t=$(date +%s)
        plog event chain "field: $(echo $todo | wc -w) single points, $FP side by side, $R ranks${STEP_MEM:+, $STEP_MEM} each"
        k=0
        for r in $todo; do
            (cd "field/$r" && cp2k_n $R) &
            k=$((k + 1)); if [ $k -ge $FP ]; then wait; k=0; fi
        done
        wait
        # a run that stopped at MAX_SCF goes on once from its own .wfn (written every 50 cycles)
        for r in $(orch_field status field); do
            if [ "$(qdw cp2k scf-restart "field/$r" | sed 's/.*: //')" = "restarted" ]; then
                plog event chain "field/$r: SCF not converged at MAX_SCF, restarted from its own .wfn"
                (cd "field/$r" && cp2k_n $R) || true      # a second failure is reported below, not fatal here
            fi
        done
    fi
    todo=$(orch_field status field)
    [ "$todo" = "all done" ] || { echo "[chain] field runs incomplete"; why field $todo; plog event chain "ERROR field runs incomplete: $todo (read field/<run>/srun.log)"; exit 1; }
    [ -z "${t:-}" ] || plog event chain "field done: $(ls -d field/*/ | wc -l) runs" --since $t
fi

# 2b. Solvation: potential and density cubes of the intact dot (restart of the geo_opt .wfn), the
#     fully stripped dot (last state of the MACE desorption ladder), the fully added dot (last state of
#     the adsorption ladder, when there is one) and the MA and MX_q monomers ->
#     COSMO on the DFT potential at the density isosurface + its area (solvation.json, read by
#     Orchestr.AI's solvation step). Needs the field step's base input. Finished runs are skipped.
if [ "{{FIELD}}" = "1" ] && [ ! -s solvation.json ]; then (   # a failure here must not stop the chain
    PROPS="{{RECORD_DIR}}/props"
    STRIP=""; kmax=-1
    for f in "$PROPS"/detach_*.xyz; do
        [ -s "$f" ] || continue; k=${f##*detach_}; k=${k%.xyz}
        if [ "$k" -gt "$kmax" ]; then kmax=$k; STRIP=$f; fi
    done
    ADD=""; kmax=-1
    for f in "$PROPS"/adsorb_*.xyz; do
        [ -s "$f" ] || continue; k=${f##*adsorb_}; k=${k%.xyz}
        if [ "$k" -gt "$kmax" ]; then kmax=$k; ADD=$f; fi
    done
    [ -s esp/runs.txt ] || qdw cp2k esp prepare field/base/cp2k.inp esp $( [ -n "$STRIP" ] && echo --stripped "$STRIP" ) \
        $( [ -n "$ADD" ] && echo --added "$ADD" ) \
        $( [ -s "$PROPS/monomers.xyz" ] && echo --monomers "$PROPS/monomers.xyz" ) > /dev/null
    if [ -n "$MONOLIB" ] && [ -s "$PROPS/monomers.xyz" ]; then
        qdw monomers link --esp esp --monomers "$PROPS/monomers.xyz" --library "$MONOLIB" > esp/monomers_link.json 2>&1 || true
        plog event chain "esp monomers: $(grep -c '"status": "linked"' esp/monomers_link.json) of $(grep -c '"run":' esp/monomers_link.json) from the monomer library"
    fi
    todo=$(qdw cp2k esp status esp)
    if [ "$todo" != "all done" ]; then
        echo "[chain] esp: $todo"
        EP=$(side_by_side {{ESP_PARALLEL}})   # the dots (intact, stripped, added) side by side if they fit
        NR=$(ranks_each $EP)
        STEP_MEM=$(step_mem $EP $NR); echo "[chain] esp: $EP side by side, $NR ranks${STEP_MEM:+ and $STEP_MEM} each"
        t=$(date +%s); plog event chain "esp: $(echo $todo | wc -w) cube runs, $EP dot(s) side by side, $NR ranks each"
        k=0
        for r in $todo; do
            case $r in
                dot_*) (cd "esp/$r" && cp2k_n $NR) &
                       k=$((k + 1)); if [ $k -ge $EP ]; then wait; k=0; fi ;;
                *)     (cd "esp/$r" && cp2k_n $(( SLURM_NTASKS < 8 ? SLURM_NTASKS : 8 ))) ;;
            esac
        done
        wait
        for r in $(qdw cp2k esp status esp); do                  # stopped at MAX_SCF: once more from its own .wfn
            case $r in dot_*)
                if [ "$(qdw cp2k scf-restart "esp/$r" | sed 's/.*: //')" = "restarted" ]; then
                    plog event chain "esp/$r: SCF not converged at MAX_SCF, restarted from its own .wfn"
                    (cd "esp/$r" && cp2k_n $NR) || true
                fi ;;
            esac
        done
    fi
    todo=$(qdw cp2k esp status esp)
    if [ "$todo" = "all done" ]; then
        qdw cp2k esp collect esp -o solvation.json > /dev/null && echo "[chain] solvation.json"
        plog event chain "esp done: solvation.json ($(ls -d esp/*/ | wc -l) runs)" ${t:+--since $t}
    else
        echo "[chain] esp runs incomplete: no DFT solvation (vacuum)"; why esp $todo
        plog event chain "WARN esp runs incomplete ($todo): no DFT solvation, vacuum (read esp/<run>/srun.log)"
    fi
) || echo "[chain] solvation step failed: solvation from the field step's charges"; fi

# 3. PDOS single point (all MO energies)
[ -d pdos ] || qdw cp2k prepare relaxed.xyz --step pdos --project "$P-pdos" --wfn "../geo_opt/$P-RESTART.wfn" --out pdos $EXTRA
if ! ls pdos/*.pdos >/dev/null 2>&1; then echo "[chain] pdos"; STEP=pdos; plog event chain "pdos started"; (cd pdos && cp2k); plog cp2k chain pdos pdos; fi

# 4. MO window: 5 eV below the valence band edge to 5 eV above the conduction band edge
[ -s window.json ] || { STEP=window; qdw cp2k window pdos/*.pdos -o window.json > /dev/null; plog window chain window.json; }

# 4b. Large dots (file TREXIO_ON_LARGE here): everything up to the MO window runs on standard nodes; the TREXIO
#     single point needs a large-memory node, so a job elsewhere stops here and submits its own continuation there.
if [ -f TREXIO_ON_LARGE ] && [ "${SLURM_JOB_PARTITION:-}" != "cpu-clx:large" ] \
   && ! ls trexio/orbitals*.h5 trexio/orbitals*.trexio >/dev/null 2>&1; then
    nl=1; [ "$(head -1 relaxed.xyz)" -ge 1500 ] && nl=2
    j=$(sbatch --parsable -p cpu-clx:large -N $nl -n $((nl * 96)) --ntasks-per-node=96 --mem=700G -t 06:00:00 chain.sh)
    echo "[chain] TREXIO on a large-memory node: job $j"; plog event chain "$P: up to the MO window done; TREXIO, trim and QDEX continue on $nl large-memory node(s): job $j"
    exit 0
fi

# 5. TREXIO single point (one diagonalisation, ADDED_MOS = virtual MOs in the window)
[ -d trexio ] || qdw cp2k prepare relaxed.xyz --step trexio --project "$P-trexio" --wfn "../geo_opt/$P-RESTART.wfn" --window window.json --out trexio $EXTRA
T=$(ls trexio/orbitals*.h5 trexio/orbitals*.trexio 2>/dev/null | head -1 || true)
if [ -z "$T" ]; then echo "[chain] trexio"; STEP=trexio; plog event chain "trexio started"; (cd trexio && cp2k); plog cp2k chain trexio trexio; T=$(ls trexio/orbitals*.h5 trexio/orbitals*.trexio 2>/dev/null | head -1); fi

# 6. Keep only the window's MOs
[ -s orbitals.h5 ] || { STEP=trim; qdw cp2k trim "$T" orbitals.h5 --window window.json > /dev/null; plog event chain "trim done: orbitals.h5 $(du -h orbitals.h5 | cut -f1)"; }

# 7. QDEX: QP gap, diagonal sBSE (Resta), SOC, fuzzy bands, PDOS and COOP -> qdex_electronic.h5,
#    qdex_excitations.h5, HOMO-1..LUMO+1 cubes; webapp.json (summary, band-edge excitons, spectra);
#    the HTML dashboards (fuzzy bands and exciton analysis, SF and SOC) rebuilt from the HDF5 files.
#    qdex: every transition in the sBSE; qdex_trap: the same without the transitions from or to
#    surface states (Facet or Localized), on the same orbitals.
if [ -s record.json ]; then
    NT=${SLURM_CPUS_ON_NODE:-8}
    for Q in qdex qdex_trap; do
        t=
        [ -s $Q/config.yaml ] || qdw qdex prepare record.json --window window.json --out $Q --threads "$NT" --cluster {{CLUSTER}} $([ $Q = qdex_trap ] && echo --trap-filter)
        if [ ! -s $Q/qdex_excitations.h5 ]; then
            echo "[chain] $Q"; STEP=$Q; t=$(date +%s); plog event chain "$Q started ($NT threads)"
            (cd $Q && qdex_run "$NT" > qdex.out 2>&1) || { plog event chain "ERROR $Q failed (read $(pwd)/$Q/qdex.out)"; exit 1; }
        fi
        qdw qdex webapp $Q -o $Q/webapp.json > /dev/null
        [ -z "${t:-}" ] || plog qdex chain $Q $Q
        if [ ! -s $Q/exciton_analysis_sf.html ]; then
            echo "[chain] $Q dashboards"
            (cd $Q && qdex_env "$QDEX" dashboards . > dashboards.out 2>&1) && plog event chain "$Q dashboards done" \
                || { echo "[chain] $Q dashboards failed ($Q/dashboards.out)"; plog event chain "WARN $Q dashboards failed (read $(pwd)/$Q/dashboards.out)"; }
        fi
    done
fi
echo "[chain] done: $(pwd)"
touch .done
