# Site kinetics: K<sub>X</sub>, Δμ<sub>f</sub>, k<sub>dir</sub>, k<sub>disp</sub>

`qdw sites` is the standard protocol for the surface-site inputs of the per-facet size law:

$$
\frac{dm}{dt} = N_A \sum_f A_f(m)\,\frac{c^*_{M,\infty} - C_{\mathrm{eq},\infty}\,e^{\Delta\mu_f(m;c_X)/k_BT}}{\delta_f(m)/D + 1/k^{(f)}_\mathrm{eff}(c_X)},
\qquad
k^{(f)}_\mathrm{eff} = \frac{k^{(f,\mathrm{edge})}_\mathrm{dir} + k^{(f,\mathrm{edge})}_\mathrm{disp}\,K^{(f)}_{X,\mathrm{edge}}\,c_X}{1 + K^{(f)}_{X,\mathrm{eff}}\,c_X}
$$

It reproduces the CdSe/CdCl<sub>2</sub> runs of October 2026 (LUMI site jobs 23272989, 23279035, 23283088; NEB jobs
23286150, 23286888, 23288336–40): the same dots, roles, environments, sites and NEB endpoints. Computed values only:
nothing reads the growth template or the webapp's placeholders.

## The recipe

`qdw sites recipe > recipe.yaml` prints the CdSe/CdCl<sub>2</sub> template: material and builder settings, the sizes
with their surface treatment (`plain`, `swap`, `half_swap`), the model and every step's settings, and the cluster
block for job scripts. `$QD_BUILDER` and `$M` are expanded when the recipe is read. `build` copies the
recipe into the run directory, and every later step (on the laptop or on the cluster) reads `run/recipe.yaml`, so one
run has one recipe.

## Steps

One run directory holds every step's JSON. GPU steps run through `qdw sites job`.

| step | where | what | output |
|---|---|---|---|
| `build` | laptop | QD_builder dot per size; every atom's role from the build manifest (`site_roles.roles`), carried by atom through the half swap; planes are never re-fitted | `dots/<tag>.json` |
| `reference` | GPU | relax each dot (L-BFGS 0.03 eV/Å, rattle-polish while it pays), **no anneal** | `reference/<tag>.json` |
| `pick` | laptop | one site per symmetry-distinct environment of each role (element/distance list within 4.5 Å in the as-built dot), the CdCl<sub>2</sub> unit and the states ref (dot–X), empty (dot\*), M (dot–CdSe) | `sites/<tag>.json`, `sites/<tag>_summary.json` |
| `states` | GPU | relax every state, partial Hessian on the site's local set (6 Å), free CdCl<sub>2</sub> | `states/<tag>.json` |
| `table` | laptop | ΔG<sub>bind,X</sub>, K<sub>X</sub>, Δμ<sub>f</sub> per size, role and environment, with the Cd count n as weight | `role_table.json` |
| `neb-prepare` | laptop | k<sub>disp</sub> and k<sub>dir</sub> endpoints from the relaxed states, one site per role (largest environment), marked xyz for review | `neb/<tag>/<site>/<kind>/` |
| `neb` | GPU | loose bands (paths A and B for k<sub>disp</sub>), window climbing image with a drift guard, PHVA at the TS and at the precursor | `.../result.json` |
| `inputs` | laptop | role table and NEB results in one file | `kinetics_inputs.json` |

```bash
qdw sites build run --recipe recipe.yaml           # copies the recipe to run/recipe.yaml; every later step reads it
qdw sites job /scratch/.../run --steps reference -o job_ref.sh       # copy run/ to the cluster, sbatch
qdw sites pick run
qdw sites job /scratch/.../run --steps states -o job_states.sh
qdw sites table run
qdw sites neb-prepare run                          # look at the *_marked.xyz first
qdw sites job /scratch/.../run --steps neb -o job_neb.sh
qdw sites inputs run
```

## Definitions

$$
\Delta G_{\mathrm{bind},X} = G(\text{dot–X}) - G(\text{dot*}) - G^\circ_{\mathrm{CdCl_2}},\qquad
K_X = e^{-\Delta G_{\mathrm{bind},X}/k_BT},\qquad
\Delta\mu_f = G(\text{dot–CdSe}) - G(\text{dot–X}) + G^\circ_{\mathrm{CdCl_2}} - g_\mathrm{bulk}
$$

$$
k_1 = \frac{k_BT}{h}\,e^{-(G_\ddagger - G_P)/k_BT},\qquad
\Delta G_P = G_P - G_\mathrm{site} - G^\circ_\mathrm{CdSe},\qquad
\Delta G_\ddagger = G_\ddagger - G_\mathrm{site} - G^\circ_\mathrm{CdSe}
$$

G = E + F<sub>vib</sub> (partial Hessian, modes below 50 cm<sup>−1</sup> held at 50); free species at 1 M with
translation and rotation. G<sub>site</sub> is G(dot–X) for k<sub>disp</sub> and G(dot\*) for k<sub>dir</sub>.

**Reacting species.** Every monomer package in solution reacts: CdSe and CdSe(CdCl<sub>2</sub>)<sub>n</sub>, n = 1–3.
All of them are taken to react with the same kinetics, computed with bare CdSe to keep the reaction clean. The flux
therefore uses the total monomer supply (all packages), never the bare CdSe share. The packages enter only through
the free energies: μ<sub>CdSe</sub> of the pool sets the driving term and Δμ<sub>f</sub>.

The size law uses one Δμ per builder role: `dmu_mean` in `role_table.json`, the plain Cd-weighted mean over the role's
environments (the cost per unit of a whole layer). There is no channel per environment and no exponential average;
summing the environments as separate channels would let the most expensive one set the dot (Ivan, 2026-10-30).

## Checks

- `build`: Cd per role, printed and stored; the rebuild from facet distances must give the same atoms.
- `pick`: every environment is listed, those without a whole CdCl<sub>2</sub> unit as `NO UNIT` (no Δk = −1 channel).
- `states`: imaginary modes per state and how many atoms moved more than 1 Å.
- `neb`: `TS_ok` — the lowest mode below −5 cm<sup>−1</sup> and along the path (overlap ≥ 0.5), the next one above
  −10 cm<sup>−1</sup> (recipe `neb.ts_*`). A soft mode can sit above
  −30 cm<sup>−1</sup>, so the count alone is not the check. `converged_CI` false means the TS is the highest image of the
  refined window (about 0.1 eV). `no_barrier` (collision-limited) is set from the path's ends: look at the path when a
  deep minimum comes before the highest image.

## Known choice kept from the 2026-10 runs

The rotational partition function of every free molecule is halved (symmetry number 2). That is right for CdCl<sub>2</sub>
but CdSe has σ = 1, so G(CdSe) is k<sub>B</sub>T ln 2 = 0.030 eV too high at 500 K and ΔG<sub>P</sub>, ΔG<sub>‡</sub> against
free CdSe 0.030 eV too low. Kept so the protocol reproduces the computed values; to be decided.

## Open in the growth model

How k<sub>dir</sub> and k<sub>disp</sub> become a velocity (Å/s) is not part of the protocol: either k<sub>1</sub> × site
density with precursor saturation, or the second-order rate from ΔG<sub>‡</sub> against 1 M, both capped at the
collision rate.
