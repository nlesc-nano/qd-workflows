# Cluster profiles

A profile is one YAML file in `clusters/` (or in a directory of `$QDW_CLUSTERS`), selected with
`--cluster <name>` or a path. It says how the cluster takes jobs and where the programs are; no
job script hard-codes a path. Shipped profiles: `nhr` (NHR@ZIB, production) and `default`
(hyperion).

## Top level

| Key | Example (nhr) | Meaning |
|---|---|---|
| `name` | `nhr` | the profile's name, passed on to jobs that call `qdw` again |
| `scheduler` | `slurm` | |
| `cores_per_node` | 192 | CPU cores of a node (CP2K jobs) |
| `request` | `cores` | `cores`: shared nodes, a job asks for a core count; `nodes`: whole nodes |
| `min_cores` | 48 | smallest CP2K job |
| `mem_per_cpu` | `7g` | memory per core, a floor for CP2K jobs |
| `node_memory_gb`, `threads_per_core` | 1400, 2 | usable node memory; Slurm counts `--mem-per-cpu` per hardware thread |
| `qos` | `normal` | default QoS |
| `account` | `bep00157` | `#SBATCH --account` of every job (omit for the user's default) |
| `setup` | `[]` | shell lines at the top of every job (e.g. `source /etc/profile.d/lmod.sh`) |
| `qdw_python` | path | the Python with qd_workflows, for `qdw` inside CP2K chains |
| `qdw_pythonpath` | `.../qd-workflows/src:.../Orchestr.AI/src` | prepended to `PYTHONPATH` for those calls |

## `cp2k`

| Key | Meaning |
|---|---|
| `module` or `env_script` | how CP2K is loaded (a module, or a script of an own build) |
| `executable` | `cp2k.psmp` / `cp2k.popt` |
| `partition`, `test_partition` | CPU partition; the 1 h partition for `--test-queue` |
| `launcher` | `srun` or `mpirun` |
| `field_parallel` | most finite-field single points side by side (13); the chain lowers it to what the memory holds |
| `omp_num_threads` | 1 (pure MPI; with threads every rank ended up on one core on NHR) |
| `data_dir`, `basis_file`, `potential_file`, `basis`, `potential` | basis and pseudopotential files and families |
| `t_c_g`, `dftd3` | truncated-Coulomb table (PBE0 in periodic cells), D3 parameters (`qdw reaction`) |
| `time` | wall time per step and the cap of a chain |
| `calibration` | `basis_functions`, `cores`, `memory_gb` of a reference run: resources scale with the square of the basis-function count |

## `qdex`

`command`, `gth_soc_file`, `max_states`, `eps_out`, `g_shell` (see [QDEX](../pipeline/qdex.md)).

## `props` (GPU and CPU stages)

| Key | Meaning |
|---|---|
| `dft_chain` | the GPU job submits one CP2K chain per record and the CPU job waits for them |
| `partition`, `test_partition`, `gres`, `constraint` | the GPU request (hyperion: the GPU model goes in `--constraint`) |
| `cpus`, `mem`, `time` | the GPU job's cores, memory and wall time |
| `conda_init`, `conda_env` | how the `qdw-compute` environment is activated |
| `env` | extra environment variables (e.g. `PYTHONPATH` to the Orchestr.AI checkout) |
| `model`, `head` | the MACE-MH-1 model file and head (`omat_pbe`) |
| `cif_dirs` | bulk CIFs (references, fuzzy bands) |
| `refs` | shared cache of the reference species (bulk, monomers) per material and head |
| `cpu.partition`, `cpu.cpus`, `cpu.mem`, `cpu.time` | the CPU stage's job |

## A new cluster

Copy `clusters/default.yaml`, change what differs, and check one record end to end:

```bash
qdw props prepare <one small record> --cluster <new> --out jobs/check
(cd jobs/check && sbatch job.sh)
```

Notes from the clusters in use:

- **NHR**: the test partitions (`cpu-genoa:test`, `gpu-a100:test`) share one running job per
  user; idle GPU nodes may be powered down and take a while to start.
- **hyperion**: `sbatch` from a non-interactive shell does not pass the `module` function on (the
  profile's `setup` sources Lmod); torch's pip wheels need `$CONDA_PREFIX/lib` first on
  `LD_LIBRARY_PATH` (the GPU job template does it).
