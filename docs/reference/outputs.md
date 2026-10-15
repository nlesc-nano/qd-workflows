# Output files

Everything a record produces stays in its directory, next to `record.json` and `start.xyz`.

## `props/` (Orchestr.AI PROPS)

One `<step>.json` per step (summary, data and the hash of its inputs, which makes it a cache
entry), plus:

| File | Step | Holds |
|---|---|---|
| `relaxed.xyz` | relax | the MACE-MH-1 minimum, centred at the centre of mass |
| `structure.json` | structure | bonds, core / surface / ligand roles, ligand denticity, start-vs-relaxed changes |
| `hessian.json`, `modes.npz` | hessian | frequencies, thermochemistry on the T grid, VDOS; normal modes and masses |
| `stability.json`, `monomers.xyz` | stability | decomposition and binding free energies, references; the monomers for the CP2K chain |
| `detachment.json`, `detach_<k>.xyz`, `detach_cache.json` | detachment | the ladder, the structures along the best path, the relaxation checkpoint |
| `wigner.json`, `md.json`, `ensemble/` | wigner, md | ensembles with energies, forces and descriptors (`*.extxyz`, `embeddings.npz`, `manifest.json`) |
| `field.json`, `field.npz` | field | Born charges, Raman tensors, polarisability, dipole, Hirshfeld charges |
| `vibspec.json`, `spectra.json` | vibspec | IR and Raman per mode, symmetry labels, mode classes; broadened spectra |
| `solvation.json` | solvation | $E_\infty$, $S$ and their source per species |
| `sites.json` | sites | site classes, orbits, free energies, the faceted shape |
| `report.json`, `properties.json`, dashboards (`*.html`, `*.png`) | report | the summaries of every step and the report |

## `cp2k/` (the CP2K chain)

| Path | Holds |
|---|---|
| `chain.sh`, `slurm-<id>.out` | the job and its log (`[chain]` lines name each step) |
| `geo_opt/`, `relaxed.xyz` | the PBE optimisation, its restart `.wfn`, the last frame |
| `field/` | `base/` (the single-point input), the 13 field runs and `mono_MA`, `mono_MX`, `runs.txt` |
| `esp/`, `solvation.json` | the cube runs (`dot_0`, `dot_<K>`, `MA`, `MX`) and their COSMO energies and areas |
| `pdos/`, `window.json` | MO energies and PDOS; the MO window |
| `trexio/`, `orbitals.h5` | the TREXIO single point and the trimmed orbitals |
| `qdex/` | `qdex_electronic.h5`, `qdex_excitations.h5`, MO cubes, `webapp.json` |
| `.done` | written when the chain finished |

## `qdw reaction` run directory

| Path | Holds |
|---|---|
| `reaction.yaml`, `reaction.json` | the input and its normalised form |
| `species/<name>/` | `dft.xyz`, `mace.json` (relaxation, frequencies, $G^\mathrm{corr}$), `cp2k/<run>/` |
| `results/results.json` | everything the page needs (per species and per energy set) |
| `results/summary.json` | the numbers at 298.15 K for each model |
| `results/reaction.html` | the interactive page |
| `results/*.png` | static plots: ΔG(T), decomposition, maps, ε dependence, solvation, recipe |

## Schemas

`schema/` (`qd-schema`) holds the shared file formats: the ensemble manifest and the CP2K
summary; the library record will follow.
