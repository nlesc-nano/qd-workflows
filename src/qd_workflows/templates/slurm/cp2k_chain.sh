#!/bin/bash
#SBATCH -J {{JOB_NAME}}
#SBATCH -t {{TIME}}
#SBATCH -n {{NTASKS}}
#SBATCH -N {{NODES}}
#SBATCH --mem-per-cpu={{MEM_PER_CPU}}
#SBATCH --qos={{QOS}}
{{SETUP}}
# CP2K track of one structure in one job: geo_opt -> PDOS -> MO window -> TREXIO -> trim.
# Every step is skipped when its output is already there, so a job stopped by the
# time limit continues where it stopped when resubmitted.
module purge
module load {{MODULE}}
export OMP_NUM_THREADS={{OMP_NUM_THREADS}}
set -e

P={{PROJECT}}
QDW_PY={{QDW_PYTHON}}
# qdw runs in its own environment; its libraries do not leak into CP2K's.
qdw() { LD_LIBRARY_PATH=$(dirname "$(dirname "$QDW_PY")")/lib "$QDW_PY" -m qd_workflows.cli "$@"; }
cp2k() { mpirun -np ${SLURM_NTASKS} {{EXECUTABLE}} -i cp2k.inp -o cp2k.out; }
EXTRA="--cluster {{CLUSTER}} --time {{TIME}} --qos {{QOS}}"

# 1. geo_opt (from the MACE-relaxed structure)
cd geo_opt
if ! grep -q "GEOMETRY OPTIMIZATION COMPLETED" cp2k.out 2>/dev/null; then
    if [ -s "$P-pos-1.xyz" ]; then      # continue from the last geometry of a stopped run
        n=$(ls cp2k.out.prev-* 2>/dev/null | wc -l)
        mv cp2k.out "cp2k.out.prev-$n" 2>/dev/null || true
        mv "$P-pos-1.xyz" "$P-pos-1.xyz.prev-$n"
        qdw cp2k relaxed "$P-pos-1.xyz.prev-$n" -o geom.xyz
    fi
    echo "[chain] geo_opt"; cp2k
fi
grep -q "GEOMETRY OPTIMIZATION COMPLETED" cp2k.out || { echo "[chain] geo_opt did not converge"; exit 1; }
cd ..
[ -s relaxed.xyz ] || qdw cp2k relaxed "geo_opt/$P-pos-1.xyz" -o relaxed.xyz

# 2. PDOS single point (all MO energies)
[ -d pdos ] || qdw cp2k prepare relaxed.xyz --step pdos --project "$P-pdos" --wfn "../geo_opt/$P-RESTART.wfn" --out pdos $EXTRA
if ! ls pdos/*.pdos >/dev/null 2>&1; then echo "[chain] pdos"; (cd pdos && cp2k); fi

# 3. MO window: 5 eV below the valence band edge to 5 eV above the conduction band edge
[ -s window.json ] || qdw cp2k window pdos/*.pdos -o window.json

# 4. TREXIO single point (one diagonalisation, ADDED_MOS = virtual MOs in the window)
[ -d trexio ] || qdw cp2k prepare relaxed.xyz --step trexio --project "$P-trexio" --wfn "../geo_opt/$P-RESTART.wfn" --window window.json --out trexio $EXTRA
T=$(ls trexio/orbitals*.h5 trexio/orbitals*.trexio 2>/dev/null | head -1 || true)
if [ -z "$T" ]; then echo "[chain] trexio"; (cd trexio && cp2k); T=$(ls trexio/orbitals*.h5 trexio/orbitals*.trexio | head -1); fi

# 5. Keep only the window's MOs
[ -s orbitals.h5 ] || qdw cp2k trim "$T" orbitals.h5 --window window.json
echo "[chain] done: $(pwd)/orbitals.h5"
