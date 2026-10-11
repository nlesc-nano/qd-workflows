#!/bin/bash
#SBATCH -J {{JOB_NAME}}
#SBATCH -t {{TIME}}
#SBATCH -n {{NTASKS}}
#SBATCH -N {{NODES}}
#SBATCH --mem-per-cpu={{MEM_PER_CPU}}
#SBATCH --qos={{QOS}}
{{PARTITION}}{{SETUP}}
# CP2K track of one structure in one job: geo_opt -> field single points -> PDOS -> MO window
# -> TREXIO -> trim -> QDEX (when record.json is here).
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
qdex_run() { env -u PYTHONPATH -u PYTHONHOME LD_LIBRARY_PATH=$(dirname "$(dirname "$QDEX")")/lib OMP_NUM_THREADS=$1 "$QDEX" --config config.yaml; }
cp2k() { {{LAUNCH}} {{EXECUTABLE}} -i cp2k.inp -o cp2k.out; }
cp2k_n() { {{LAUNCH_N}} {{EXECUTABLE}} -i cp2k.inp -o cp2k.out > /dev/null 2>&1; }   # $1 ranks
EXTRA="--cluster {{CLUSTER}} --time {{TIME}} --qos {{QOS}}"

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
    echo "[chain] geo_opt"; cp2k
fi
grep -q "GEOMETRY OPTIMIZATION COMPLETED" cp2k.out || { echo "[chain] geo_opt did not converge"; exit 1; }
cd ..
[ -s relaxed.xyz ] || qdw cp2k relaxed "geo_opt/$P-pos-1.xyz" -o relaxed.xyz

# 2. Finite-field single points: 13 restarted from the geo_opt .wfn ({{FIELD_PARALLEL}} at a time, the
#    job's ranks shared among them) + the MA and MX_q monomers -> Born charges, Raman tensors,
#    Hirshfeld charges (Orchestr.AI field step). Finished runs are skipped.
if [ "{{FIELD}}" = "1" ]; then
    [ -s field/base/cp2k.inp ] || qdw cp2k prepare relaxed.xyz --step sample --project "$P-field" --wfn "../../geo_opt/$P-RESTART.wfn" --out field/base $EXTRA
    MONO="{{RECORD_DIR}}/props/monomers.xyz"
    [ -s field/runs.txt ] || orch_field prepare field/base/cp2k.inp field $( [ -s "$MONO" ] && echo --monomers "$MONO" )
    todo=$(orch_field status field)
    if [ "$todo" != "all done" ]; then
        echo "[chain] field: $todo"
        R=$(( SLURM_NTASKS / {{FIELD_PARALLEL}} )); [ $R -ge 1 ] || R=1
        k=0
        for r in $todo; do
            n=$R; case $r in mono_*) n=$(( R < 8 ? R : 8 ));; esac
            (cd "field/$r" && cp2k_n $n) &
            k=$((k + 1)); if [ $k -ge {{FIELD_PARALLEL}} ]; then wait; k=0; fi
        done
        wait
    fi
    [ "$(orch_field status field)" = "all done" ] || { echo "[chain] field runs incomplete"; exit 1; }
fi

# 3. PDOS single point (all MO energies)
[ -d pdos ] || qdw cp2k prepare relaxed.xyz --step pdos --project "$P-pdos" --wfn "../geo_opt/$P-RESTART.wfn" --out pdos $EXTRA
if ! ls pdos/*.pdos >/dev/null 2>&1; then echo "[chain] pdos"; (cd pdos && cp2k); fi

# 4. MO window: 5 eV below the valence band edge to 5 eV above the conduction band edge
[ -s window.json ] || qdw cp2k window pdos/*.pdos -o window.json

# 5. TREXIO single point (one diagonalisation, ADDED_MOS = virtual MOs in the window)
[ -d trexio ] || qdw cp2k prepare relaxed.xyz --step trexio --project "$P-trexio" --wfn "../geo_opt/$P-RESTART.wfn" --window window.json --out trexio $EXTRA
T=$(ls trexio/orbitals*.h5 trexio/orbitals*.trexio 2>/dev/null | head -1 || true)
if [ -z "$T" ]; then echo "[chain] trexio"; (cd trexio && cp2k); T=$(ls trexio/orbitals*.h5 trexio/orbitals*.trexio 2>/dev/null | head -1); fi

# 6. Keep only the window's MOs
[ -s orbitals.h5 ] || qdw cp2k trim "$T" orbitals.h5 --window window.json

# 7. QDEX: QP gap, diagonal sBSE (Resta), SOC, fuzzy bands, PDOS and COOP -> qdex_electronic.h5,
#    qdex_excitations.h5, HOMO-1..LUMO+1 cubes; webapp.json (summary, band-edge excitons, spectra)
if [ -s record.json ]; then
    NT=${SLURM_CPUS_ON_NODE:-8}
    [ -s qdex/config.yaml ] || qdw qdex prepare record.json --window window.json --out qdex --threads "$NT" --cluster {{CLUSTER}}
    if [ ! -s qdex/qdex_excitations.h5 ]; then echo "[chain] qdex"; (cd qdex && qdex_run "$NT" > qdex.out 2>&1); fi
    qdw qdex webapp qdex -o qdex/webapp.json > /dev/null
fi
echo "[chain] done: $(pwd)"
touch .done
