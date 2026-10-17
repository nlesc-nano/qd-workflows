# 2 · CP2K chain

One Slurm job per dot (up to 2,000 atoms) runs the whole DFT track. It is written by
`qdw cp2k chain` (or by `qdw props chains` for a batch) as `<record>/cp2k/chain.sh`. Every step is
skipped when its output is already there, so a job stopped by the time limit continues where it
stopped when resubmitted.

```mermaid
flowchart TB
    GEO["1 · geo_opt<br/>from the MACE structure"] --> FIELD["2 · field<br/>13 finite-field single points<br/>+ MA, MX_q monomers"]
    FIELD --> ESP["2b · esp<br/>potential + density cubes<br/>→ solvation.json"]
    ESP --> PDOS["3 · pdos<br/>all MO energies"] --> WIN["4 · MO window"]
    WIN --> TRX["5 · TREXIO<br/>one diagonalisation"] --> TRIM["6 · trim"] --> QDEX["7 · QDEX<br/>qdex + qdex_trap<br/>→ webapp.json, dashboards"]
    GEO -. ".wfn" .-> FIELD
    GEO -. ".wfn" .-> ESP
    GEO -. ".wfn" .-> PDOS
    GEO -. ".wfn" .-> TRX
```

## Common DFT settings

All steps share one setup (templates in `src/qd_workflows/templates/cp2k/`):

| Setting | Value |
|---|---|
| Functional | PBE |
| Method | GAPW (`FORCE_PAW T`, `GAPW_ACCURATE_XCINT T`, CP2K ≥ 2026.2) |
| Basis, pseudopotentials | DZVP-MOLOPT-GGA-GTH, GTH-GGA (`BASIS_MOLOPT_UZH`, `POTENTIAL_UZH`), per valence charge |
| Grid | `CUTOFF 400`, `REL_CUTOFF 60`, `NGRIDS 4` |
| Cell | non-periodic, the dot's extent + 12 Å (6 Å on each side); Poisson solver `MULTIPOLE` |
| SCF | OT, DIIS minimiser, `FULL_KINETIC` preconditioner, `EPS_SCF 1e-6`, `MAX_SCF 200`, no outer loop |

Every step after geo_opt restarts from the geo_opt wavefunction.

## 1 · geo_opt

Geometry optimisation from the MACE-relaxed structure (BFGS; L-BFGS above 1,000 atoms).
`IGNORE_CONVERGENCE_FAILURE` lets the first step use all 200 SCF cycles; later steps converge in
far fewer. A stopped run continues from the last frame of its trajectory.

## 2 · field (finite-field single points)

Thirteen single points in a uniform electric field (0.002 a.u. ≈ 0.1 V/Å), restarted from the
geo_opt wavefunction, written and read by Orchestr.AI `qd.field`:

| Run | Field |
|---|---|
| `f0` | none: reference forces, dipole, Hirshfeld charges |
| `fxp`, `fxm`, `fyp`, `fym`, `fzp`, `fzm` | $\pm E$ along one axis |
| `fxypp`, `fxymm`, `fxzpp`, `fxzmm`, `fyzpp`, `fyzmm` | $(+E,+E)$ and $(-E,-E)$ along two axes |

plus zero-field single points of the MA and MX$_q$ monomers. They give the Born charges, the Raman
tensors and the polarisability ([IR and Raman](../theory/spectra.md)). On Cd16 and Cd68: Born
charge sum rule to $10^{-4}\,e$, polarisability within 0.01 % of linear response.

The runs go side by side, as many as the job's memory holds (see [Memory](#memory)), the job's
ranks shared among them.

## 2b · esp (solvation cubes)

Hartree-potential and electron-density cubes (stride 2) of up to five structures, with the field
step's single-point settings:

| Run | Structure | Start |
|---|---|---|
| `dot_0` | the intact dot | restart from the geo_opt .wfn |
| `dot_K` | the fully stripped dot, last state of the desorption ladder (`props/detach_K.xyz`) | atomic guess, own cell |
| `dot_addK` | the fully added dot, last state of the adsorption ladder (`props/adsorb_K.xyz`), when there is one | atomic guess, own cell |
| `MA`, `MX` | the monomers (`props/monomers.xyz`) | atomic guess, cell ≥ 15 Å |

The dots run side by side when the memory holds them, the monomers after them. Orchestr.AI
interpolates the stripped states between `dot_0` and `dot_K`; the webapp interpolates the added
states (ladder up) between `dot_0` and `dot_addK` (`k_added` in `solvation.json`).

`qdw cp2k esp collect` turns each pair of cubes into the COSMO conductor energy on the density
isosurface and the isosurface area, written to `cp2k/solvation.json`
([Solvation](../theory/solvation.md)). A failure of this step does not stop the chain.

## 3–6 · pdos, window, TREXIO, trim

- **pdos**: an OT single point that prints the energies and occupations of all MOs, with PDOS.
- **window** (`qdw cp2k window`): the MOs from 5 eV below the valence band edge to 5 eV above the
  conduction band edge; `added_mos` is the number of virtual MOs in the window.
- **TREXIO**: one diagonalisation cycle reading the geo_opt .wfn, `ADDED_MOS` = the window's
  virtual MOs, the orbitals written in TREXIO format (spherical AOs).
- **trim** (`qdw cp2k trim`): keeps only the window's MOs in the HDF5 file; QDEX reads any MO
  subset with its energies and occupations.

## 7 · QDEX

With the record's `record.json` present, QDEX runs twice on the trimmed orbitals, on the first
node's cores ([QDEX](qdex.md)):

| Run | sBSE |
|---|---|
| `qdex/` | every transition of the balanced window |
| `qdex_trap/` | the same, without the transitions from or to surface states (`excitations: trap_filter: true`) |

In each, `qdw qdex webapp` writes `webapp.json`, and `qdex dashboards .` rebuilds the HTML pages
from the HDF5 files (log in `dashboards.out`; a failure there is reported and the chain goes on).
The trap-filtered run reads the same `../orbitals.h5` and takes about a minute for Cd68Se55Cl26.

## Resources

Cores and memory grow with the square of the number of basis functions, scaled from the
profile's calibration run (`cp2k.calibration`: 24,000 basis functions, 240 cores, 340 GB for a
4.2 nm CdSe dot), at least `min_cores`, at least `mem_per_cpu` per core. The wall time of the
whole chain comes from the cost of a PBE geometry step on CdSe dots: 1.5 × the estimate, at least
one hour, at most the profile's cap.

### Memory

Each GAPW single point holds the cell's real-space grids, which outweigh the basis for small and
medium dots: Cd16Se13Cl6 in a 23.6 Å cube used 13.3 GB over 3 ranks, whatever the rank count. The
chain scales this with the cell volume,

$$
M_\mathrm{run} \approx 13.5\ \mathrm{GB}\times\frac{V_\mathrm{cell}}{(23.631\ \text{Å})^3},
$$

reads at run time the memory Slurm granted, and runs side by side at most
$\lfloor 0.8\, M_\mathrm{job} / M_\mathrm{run}\rfloor$ field single points (at most 13) and esp dots
(at most 3). Examples: 5 at a time in 96 GB, all 13 in 336 GB (NHR, 48 × 7 GB).

The job asks for its memory per node, `#SBATCH --mem` = the ranks on a node × `mem_per_cpu`, at
most `node_memory_gb`, with the ranks spread evenly over the nodes (`--ntasks-per-node`). It does
not ask per core: on NHR's cpu-genoa a job with `-n 48 --mem-per-cpu=7g` was planned on an idle
node in every backfill cycle and never started (the Cd68 chain waited over an hour), while the same
total as `--mem=336G` starts at once. Asked per node, the memory belongs to the job, not to its
ranks, and an `srun` step that names none takes all of it: the next run side by side waits for it.
So every run side by side asks for its share,

$$
M_\mathrm{step} = \frac{0.9\, M_\mathrm{node}}{k},\qquad
k = \left\lceil \frac{n\, s}{N} \right\rceil,\qquad s = \left\lceil \frac{r}{p} \right\rceil
$$

per node for $n$ runs of $r$ ranks side by side on $N$ nodes of $p$ ranks (`srun --mem`). A run that
fits on one node ($s = 1$) stays on one (`srun -N 1`): spread over two it would hold its memory on
both and talk across them. On the Cd176 pilot (2 CLX nodes of 96 cores) the 13 field runs of 14
ranks each were spread over both nodes and only 8 fitted at once. Since then, on several nodes each
run gets as many ranks as let $\lceil n/N \rceil$ runs share a node (13 runs on 2 × 96 cores: 7 per
node, 13 ranks each), and all of them start together. The chain prints the runs, ranks and
memory of each step (`[chain] field: 13 side by side, 3 ranks and 23819M each`); what `srun` and
CP2K print goes to `srun.log` in each run's directory, and the chain prints its last lines for the
runs that did not finish. The reaction job (`qdw reaction`) does the same for its species, with a
share by cores.

## Running a chain

```bash
qdw cp2k chain relaxed.xyz --project <id> --record record.json --out cp2k --cluster nhr
(cd cp2k && sbatch chain.sh)
```

Single steps can also be prepared alone:

```bash
qdw cp2k prepare relaxed.xyz --step geo_opt --project <id> --out run/geo_opt
qdw cp2k window run/pdos/*.pdos -o run/window.json
qdw cp2k prepare geo.xyz --step trexio --wfn ../geo_opt/<id>-RESTART.wfn --window run/window.json --out run/trexio
qdw cp2k trim run/trexio/orbitals.h5 run/orbitals.h5 --window run/window.json
```
