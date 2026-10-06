# qd-workflows

Workflows that tie the quantum-dot packages together: structures from
QD_builder, MACE-MH-1 properties and Wigner ensembles from Orchestr.AI,
CP2K PBE and QDEX for the QDSpace webapp, and the fine-tuning loop.

The packages keep their own repositories; this one pins them
(`components.yaml`), describes each HPC cluster (`clusters/`), and holds the job
templates and glue code (`src/qd_workflows`).

## Layout

| Path | Holds |
| --- | --- |
| `components.yaml` | repository and commit of each package a result was made with |
| `clusters/<name>.yaml` | how one HPC cluster takes jobs: cores per node, core or node requests, memory, modules, data paths |
| `src/qd_workflows/templates/cp2k/` | CP2K inputs (geo_opt, PDOS, TREXIO, sample energy-force) with `{{PLACEHOLDERS}}` |
| `src/qd_workflows/templates/slurm/` | job scripts |
| `src/qd_workflows/cp2k/` | input generation, MO window, TREXIO trimming, force parsing |
| `src/qd_workflows/compare.py` | `qdw compare`: two `props/properties.json` value by value, with tolerances (Phase 0 gate) |
| `schema/` | `qd-schema`, the shared file formats (ensemble manifest, CP2K summary; later the library record) |
| `envs/compute-gpu.yml`, `containers/compute.def` | GPU compute environment and its Apptainer image (untested) |

## CP2K track (webapp structures up to 2,000 atoms)

```
geo_opt (from the MACE-relaxed structure)
  -> pdos      OT single point, DOS over all unoccupied MOs + PDOS
  -> window    MOs from 5 eV below the valence band edge to 5 eV above the conduction band edge
  -> trexio    one diagonalisation cycle, ADDED_MOS = virtual MOs in the window, TREXIO print
  -> trim      drop the occupied MOs below the window from the .h5
  -> QDEX      QP gap (bulk, scaled vertex), diagonal sBSE (Resta) balanced to the window,
               SOC, fuzzy bands, PDOS, COOP (with the cation-ligand pair); summary.json
```

One job runs it all and skips finished steps when resubmitted:

```bash
qdw cp2k chain props/relaxed.xyz --project <id> --record record.json --out cp2k
(cd cp2k && sbatch chain.sh)
```

QDEX writes `qdex_electronic.h5` (orbitals, projections, PDOS, COOP, fuzzy-band
weights), `qdex_excitations.h5` (every exciton up to 4.5 eV with its descriptors)
and the HOMO-1..LUMO+1 cubes on a 0.8 Å grid; no dashboards. `qdex/summary.json`
holds what the webapp lists: DFT and QP HOMO, LUMO and gap, lowest and first
bright exciton (f >= 0.05), spin-free and with SOC. QDEX runs
in its own environment (`envs/qdex.yml`, then `pip install --no-build-isolation -e QDEX`)
on the cores of the chain job's first node.

Every step after geo_opt restarts from the geo_opt `.wfn`. Wigner samples use
the `sample` input (energy and forces, forces written to `forces.xyz`, no
restart file written).

Single steps:

```bash
qdw cp2k prepare relaxed.xyz --step geo_opt --project CdSe-Se-Cd68Se55Cl26-clean --out run/geo_opt
qdw cp2k window run/pdos/*.pdos -o run/window.json
qdw cp2k prepare geo.xyz --step trexio --wfn ../geo_opt/X-RESTART.wfn --window run/window.json --out run/trexio
qdw cp2k trim run/trexio/orbitals.h5 run/orbitals.h5 --window run/window.json
```

Resources come from the cluster profile: the core count scales with the square
of the number of basis functions, calibrated on a 4.2 nm CdSe run (24k basis
functions, 240 cores, 340 GB).

## Phase 0 gate

Run the properties of Cd16Se13Cl6 and Cd68Se55Cl26 on the cluster with
Orchestr.AI (`run_type: PROPS`, branch `feature/qd-properties`) and compare with
the Mac results:

```bash
qdw compare <webapp>/.../CdSe-Se-Cd16Se13Cl6-clean/props/properties.json <cluster run>/props/properties.json
```

The Mac ran on the Apple GPU in float32 and the cluster runs CUDA in float64,
so values are compared with tolerances (2 cm-1, 2 meV, 0.5 meV/atom, 2e-3 Å,
relative 1e-3).

Where only short jobs start soon, `qdw props prepare --pieces` makes the job stop
a minute before its limit and submit itself again until the run is done.

## Size benchmark (Phase 1)

```bash
qdw props bench --out jobs/bench --sizes 100,300,1000,2000,5000
```

times one force call, batched force calls, the analytic Hessian and the batched
finite-difference Hessian on CdSe spheres of each size on one GPU, and fits
t = a N^b to each (`bench.json`).
