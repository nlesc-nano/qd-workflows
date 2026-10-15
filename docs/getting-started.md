# Getting started

## Code

```bash
git clone -b feature/qd-properties git@github.com:nlesc-nano/Orchestr.AI.git
git clone git@github.com:nlesc-nano/QD_Builder.git          # bulk CIFs (examples/)
git clone git@github.com:nlesc-nano/qd-workflows.git
```

`components.yaml` names the commit of each package that a result was made with; check out those
commits to reproduce a run.

## Environments

| Environment | Used by | From |
|---|---|---|
| `qdw-compute` | the GPU stage (MACE-MH-1), the CPU stage, `qdw` inside the CP2K chain | `envs/compute-gpu.yml`, or the Apptainer image `containers/compute.def` |
| `qdex` | QDEX at the end of the chain | `envs/qdex.yml`, then `pip install --no-build-isolation -e QDEX` |
| CP2K ≥ 2026.2 | the CP2K chain | a module or an own build (GAPW, TREXIO, MULTIPOLE Poisson solver) |

```bash
conda env create -f envs/compute-gpu.yml        # qdw-compute
conda activate qdw-compute
pip install --no-deps -e Orchestr.AI -e qd-workflows
pip install scikit-image                         # the solvation step (marching cubes)
```

The MACE-MH-1 model file, the bulk CIF directories and the shared reference cache are set in the
cluster profile (`props.model`, `props.cif_dirs`, `props.refs`).

## A cluster profile

Every cluster is one YAML file in `clusters/` (or in a directory listed in `$QDW_CLUSTERS`): how
it takes jobs (cores, memory, partitions, account), where CP2K and the environments are, and the
CP2K resource calibration. See [Cluster profiles](reference/clusters.md). `clusters/nhr.yaml` is
the production profile (NHR@ZIB); `clusters/default.yaml` is hyperion.

## A batch of records

```bash
qdw props prepare <record dirs> --cluster nhr --out jobs/batch
(cd jobs/batch && sbatch job.sh)
```

- `job.sh` is the GPU job (MACE-MH-1). When it ends, it submits one CP2K chain per record and then
  `job_cpu.sh`, the CPU stage, which waits for every chain.
- Records above 2,000 atoms are left out; `prepare` lists them under `skipped_too_large`.
- Every step is cached in `<record>/props/` and `<record>/cp2k/`: a job that stops is resubmitted
  and continues where it stopped.

Results per record: `props/properties.json` and the dashboards (report), `cp2k/qdex/webapp.json`
and the QDEX files. See [Output files](reference/outputs.md).

## One CP2K chain by hand

```bash
qdw cp2k chain props/relaxed.xyz --project <record id> --record record.json --out cp2k --cluster nhr
(cd cp2k && sbatch chain.sh)
```

## A reaction free energy

```bash
qdw reaction check   reaction.yaml
qdw reaction prepare reaction.yaml --cluster nhr --out runs/x --test-queue
(cd runs/x && bash submit.sh)
qdw reaction analyse runs/x
```

See [Reaction free energies](reactions.md).
