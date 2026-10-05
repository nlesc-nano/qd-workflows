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

## CP2K track (webapp structures up to 2,000 atoms)

```
geo_opt (from the MACE-relaxed structure)
  -> pdos      OT single point, DOS over all unoccupied MOs + PDOS
  -> window    MOs from 5 eV below the valence band edge to 5 eV above the conduction band edge
  -> trexio    one diagonalisation cycle, ADDED_MOS = virtual MOs in the window, TREXIO print
  -> trim      drop the occupied MOs below the window from the .h5
  -> QDEX
```

Every step after geo_opt restarts from the geo_opt `.wfn`. Wigner samples use
the `sample` input (energy and forces, forces written to `forces.xyz`, no
restart file written).

```bash
qdw cp2k prepare relaxed.xyz --step geo_opt --project CdSe-Se-Cd68Se55Cl26-clean --out run/geo_opt
qdw cp2k window run/pdos/*.pdos -o run/window.json
qdw cp2k prepare geo.xyz --step trexio --wfn ../geo_opt/X-RESTART.wfn --window run/window.json --out run/trexio
qdw cp2k trim run/trexio/orbitals.h5 run/orbitals.h5 --window run/window.json
```

Resources come from the cluster profile: the core count scales with the square
of the number of basis functions, calibrated on a 4.2 nm CdSe run (24k basis
functions, 240 cores, 340 GB).
