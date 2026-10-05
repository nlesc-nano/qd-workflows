# Phase 0 runbook: Cd16 and Cd68 on the cluster

Goal: the MACE-MH-1 / xtb properties of Cd16Se13Cl6 and Cd68Se55Cl26, rerun on
the cluster with Orchestr.AI (`run_type: PROPS`), match the Mac results; and the
CP2K PBE track produces a geo_opt, a PDOS, a window, a TREXIO file and QDEX
dashboards for Cd16.

## 1. Code

```bash
git clone -b feature/qd-properties git@github.com:nlesc-nano/Orchestr.AI.git
git clone git@github.com:nlesc-nano/QD_Builder.git          # bulk CIFs (examples/)
git clone git@github.com:nlesc-nano/qd-workflows.git        # private repository
```

## 2. Environment (GPU node)

Either the container (`containers/compute.def`; put the static g-xTB binary
from https://github.com/grimme-lab/g-xtb in `third_party/gxtb/bin/xtb` first):

```bash
apptainer build qdw-compute.sif containers/compute.def
apptainer exec --nv qdw-compute.sif pip install --no-deps -e Orchestr.AI -e qd-workflows
```

or a conda environment from `envs/compute-gpu.yml`, with g-xTB unpacked anywhere.

Then:

```bash
export QDPROPS_MACE_MODEL=/path/to/macemh1model           # same file as on the Mac (sha256 a522eb7f...)
export QDPROPS_XTB=$(which xtb)                           # xtb 6.7.1, GFN2-xTB
export QDPROPS_GXTB=/path/to/g-xtb/bin/xtb                # g-xTB
export QDPROPS_CIF_DIRS=$PWD/QD_Builder/examples/cifs:$PWD/QD_Builder/examples/library/cifs
export QDPROPS_REFS=/shared/path/qdprops_references       # shared by all jobs
export ORCHESTRAI_SINGLE_ENV=1                            # one environment: no re-dispatch
```

## 3. Properties (MACE track)

Copy `record.json` and `start.xyz` of both records (QDSpaceWebApp
`qd-frontend/public/II-VI/CdSe/builder/<id>/`) to a run directory, point
`Orchestr.AI/config_files/postprocessing/props/config_props.yaml` at them and set
`model_path`, then on a GPU node:

```bash
python -m orchestr_ai.postprocessing config_props.yaml
qdw compare <webapp>/.../CdSe-Se-Cd16Se13Cl6-clean/props/properties.json run/CdSe-Se-Cd16Se13Cl6-clean/props/properties.json
qdw compare <webapp>/.../CdSe-Se-Cd68Se55Cl26-clean/props/properties.json run/CdSe-Se-Cd68Se55Cl26-clean/props/properties.json
```

Gate: no value outside tolerance, or every difference explained by float32
(Mac) vs float64 (cluster).

## 4. CP2K track (CPU nodes), Cd16

```bash
R=run/CdSe-Se-Cd16Se13Cl6-clean
P=CdSe-Se-Cd16Se13Cl6-clean
qdw cp2k prepare $R/props/relaxed.xyz --step geo_opt --project $P --out $R/cp2k/geo_opt
(cd $R/cp2k/geo_opt && sbatch job.sh)
# after it converges:
qdw cp2k relaxed $R/cp2k/geo_opt/$P-pos-1.xyz -o $R/cp2k/relaxed.xyz
qdw cp2k prepare $R/cp2k/relaxed.xyz --step pdos --project $P-pdos --wfn ../geo_opt/$P-RESTART.wfn --out $R/cp2k/pdos
(cd $R/cp2k/pdos && sbatch job.sh)
qdw cp2k window $R/cp2k/pdos/*.pdos -o $R/cp2k/window.json
qdw cp2k prepare $R/cp2k/relaxed.xyz --step trexio --project $P-trexio --wfn ../geo_opt/$P-RESTART.wfn --window $R/cp2k/window.json --out $R/cp2k/trexio
(cd $R/cp2k/trexio && sbatch job.sh)
qdw cp2k trim $R/cp2k/trexio/orbitals.h5 $R/cp2k/orbitals.h5 --window $R/cp2k/window.json
```

Check on the way:

- `--data-dir` (or `QDW_CP2K_DATA`) points at the cluster's `cp2k_basis` when
  preparing on a login node; the inputs always use `clusters/default.yaml`'s path.
- The PDOS file names (`*-k1-1.pdos` before the new `&DOS`/`&PDOS` layout):
  `qdw cp2k window` reads MO index, eigenvalue and occupation from the first
  three columns.
- TREXIO's output name (`orbitals.h5` or `orbitals.trexio`, depending on the
  build) and that `ADDED_MOS` = the window's virtual MOs gives the same MOs as
  `ADDED_MOS -1` inside the window.
- A sample single point (`--step sample`) writes `forces.xyz` next to `cp2k.out`;
  `qdw cp2k label geom.xyz --out cp2k.out --forces forces.xyz` refuses it if any
  force is missing.
