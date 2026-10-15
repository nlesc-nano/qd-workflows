# Validation

## Phase 0 gate

Cd16Se13Cl6 and Cd68Se55Cl26, rerun on the cluster with Orchestr.AI (`run_type: PROPS`), must
match the results made on the Mac; the CP2K track must produce a geo_opt, a PDOS, a window, a
TREXIO file and the QDEX outputs. Step by step: [Phase 0 runbook](../phase0_runbook.md).

```bash
qdw compare <webapp>/.../CdSe-Se-Cd16Se13Cl6-clean/props/properties.json <cluster run>/props/properties.json
```

The Mac runs on the Apple GPU in float32 and the cluster in CUDA float64, so values are compared
with tolerances; the first key fragment that matches decides:

| Key fragment | Absolute tolerance |
|---|---|
| `dipole_D` | 0.001 D |
| `cm1` | 2 cm⁻¹ |
| `meV` | 0.5 meV (also meV/atom) |
| `eV_A` | 0.002 eV/Å |
| `_A` | 0.002 Å |
| `eV` | 0.002 eV |
| `_K` | 1 K |
| anything else | relative 10⁻³ |

Bookkeeping (timings, cache counters, versions, paths, devices) is skipped. The desorption path
is judged on its cumulative energy: no worse than the reference by 0.05 eV.

## Phase 1: size benchmark

```bash
qdw props bench --out jobs/bench --sizes 100,300,1000,2000,5000
```

times one force call, batched force calls, the analytic Hessian and the batched finite-difference
Hessian on CdSe spheres of each size on one GPU, and fits $t = aN^b$ to each (`bench.json`).
Measured on an A100 80 GB (hyperion, float64):

| Atoms | One force call | Batched, per structure | Analytic Hessian | Peak memory (analytic) |
|---|---|---|---|---|
| 100 | 0.036 s | 0.014 s | 15 s | 2.4 GB |
| 300 | 0.060 s | 0.042 s | 86 s | 7.0 GB |
| 1,000 | 0.16 s | 0.15 s | 713 s | 23 GB |
| 2,000 | 0.32 s | | | |

MACE-MH-1 in float64 needs about 11.5 MB per atom for a force call, which sets the batch limit
of 4,000 atoms per call (46 GB).

## End-to-end check

One small record through all three stages (GPU stage, CP2K chain, CPU stage) on the production
cluster, checking that every chain step finishes, `solvation.json` holds the intact dot, the
fully stripped dot and both monomers, and the report is built.
