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
| `webapp.json` (~40 kB) | what the webapp loads, below |

`webapp.json` (`qdw qdex webapp`):

- **summary**: DFT and QP HOMO, LUMO and gap; the lowest exciton and the first bright one
  ($f \ge 0.05$); spin-free and with SOC;
- **band-edge excitons**: up to 0.5 eV above the lowest (20 to 300 states) with their descriptors;
- **absorption spectra** of every computed state on a fixed grid (0.5–6.0 eV), broadenings 0.03
  and 0.10 eV.

`qdw qdex summary` prints the frontier levels, gaps and excitons as JSON.
