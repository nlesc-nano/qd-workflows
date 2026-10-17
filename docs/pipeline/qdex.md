# QDEX

QDEX is the last step of the CP2K chain. It reads the trimmed TREXIO orbitals (the MO window
from 5 eV below the valence band edge to 5 eV above the conduction band edge) and computes:

- the **quasiparticle gap**, from a bulk-calibrated, scaled-vertex model;
- **excitons** from a diagonal screened Bethe–Salpeter equation (Resta screening), with the
  occupied and virtual MOs balanced in the window (at most `max_states` = 1,000);
- **spin–orbit coupling** (GTH SOC potentials, `qdex.gth_soc_file`);
- **fuzzy bands** (the MOs unfolded on the bulk Brillouin zone, `g_shell` = 2: 125 G replicas),
  **PDOS** and **COOP** for the cation–anion, like-element and cation–ligand pairs.

The QDEX equations are in the QDEX documentation; qd-workflows prepares its input
(`qdw qdex prepare`: the record, the window, the CIF for the fuzzy bands, the cores) and reads
its output.

## Settings from the cluster profile

| Key | Value | Meaning |
|---|---|---|
| `qdex.command` | path | the `qdex` executable of the `qdex` environment |
| `qdex.gth_soc_file` | `GTH_SOC_POTENTIALS_2024` | SOC pseudopotential parameters |
| `qdex.max_states` | 1000 | occupied + virtual MOs in the diagonal sBSE |
| `qdex.eps_out` | 2.24 | dielectric constant outside the dot (toluene) |
| `qdex.g_shell` | 2 | G shells folded into the first Brillouin zone |

## Outputs

| File | Holds |
|---|---|
| `qdex_electronic.h5` | orbitals, projections, PDOS, COOP, fuzzy-band weights |
| `qdex_excitations.h5` | every exciton up to 4.5 eV with its descriptors |
| `spatial_MO_*.cube` | HOMO−1 … LUMO+1 on a 0.8 Å grid |
| `webapp.json` (~140 kB) | what the webapp loads, below |
| `fuzzy_dashboard_sf.html`, `fuzzy_dashboard_soc.html` | fuzzy bands, PDOS, IPR, surface/core character, COOP, the MO cubes |
| `exciton_analysis_sf.html`, `exciton_analysis_soc.html` | the lowest 100 excitons: descriptors, NTOs, transition densities |
| `spectrum_{sf,soc}.dat`, `spectrum_*.png` | absorption spectra |
| `qdex.out`, `dashboards.out` | the logs of the run and of the dashboards |

`webapp.json` (`qdw qdex webapp`):

- **summary**: DFT and QP HOMO, LUMO and gap; the lowest exciton and the first bright one
  ($f \ge 0.05$); spin-free and with SOC;
- **band-edge excitons**: up to 0.5 eV above the lowest (20 to 300 states) with their descriptors;
- **absorption spectra** of every computed state on a fixed grid (0.5–6.0 eV), broadenings 0.03
  and 0.10 eV, and the sticks: energy and $f$ of every stored state with $f \ge 10^{-5}$, from
  which the webapp broadens the spectrum at any width.

`qdw qdex summary` prints the frontier levels, gaps and excitons as JSON.

QDEX itself writes the HDF5 files, the MO cubes and the spectra as data (`output.html: false`,
`plot: false`): `qdex dashboards <run dir>` rebuilds the HTML pages and the spectrum plots from `qdex_electronic.h5` and
`qdex_excitations.h5` without recomputing anything, so the pages can be remade whenever they are
needed. The chain does it after each run.

## Surface traps filtered out

The chain runs QDEX a second time in `cp2k/qdex_trap/` with the same settings plus
`excitations: trap_filter: true` (`qdw qdex prepare --trap-filter`): every MO of the active space
is classified by its angular coverage of the dot, and the transitions from or to a *Facet* or
*Localized* (surface) state are left out of the spin-free and SOC sBSE. It reads the
same `../orbitals.h5` and writes the same files, with its own `webapp.json` and dashboards; the
webapp import publishes them as `qdex/*_traps_filtered.{html,json}` next to the unfiltered ones.

On Cd68Se55Cl26 the lowest spin-free exciton moves from 2.69 to 2.77 eV and its oscillator
strength from 0.32 to 0.55; the SOC excitons stay the same.
