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

### Finite-difference against autograd Hessian

On the production GPU (NHR, A100 80 GB, float64, cuEquivariance), finite differences with graph
reuse and 4,000 atoms per call, compared with the stored autograd Hessian of the same record
(frequencies above 50 cm⁻¹):

| Atoms | Record | FD time | Force calls | Max Δν (cm⁻¹) | Median Δν (cm⁻¹) |
|---|---|---|---|---|---|
| 141 | CdSe Cd60Se39Cl42 | 6.7 s | 31 | 0.01 | 0.004 |
| 149 | CdSe Cd68 | 6.7 s | 35 | 0.01 | 0.003 |
| 316 | InP In119P80Cl117 | 19.3 s | 158 | 0.03 | 0.006 |
| 507 | PbS Pb236S201Cl70 | 51.7 s | 444 | 0.03 | 0.020 |
| 679 | ZnSe Zn301Se224Cl154 | 100.5 s | 891 | 0.02 | 0.004 |
| 957 | HgTe Hg453Te402Cl102 | 173.1 s | 1,436 | 0.02 | 0.001 |
| 1,946 | InAs In844As715Cl387 | 737.2 s | 5,838 | 0.02 | 0.003 |

The autograd Hessian takes about 3–3.5× longer at the same sizes (about 43 min at 2,000 atoms),
because cuEquivariance kernels have no `vmap` batching rule and MACE falls back to one row per
backward pass.

## End-to-end check

One small record through all three stages (GPU stage, CP2K chain, CPU stage) on the production
cluster, checking that every chain step finishes, `solvation.json` holds the intact dot, the
fully stripped dot and both monomers, and the report is built.
