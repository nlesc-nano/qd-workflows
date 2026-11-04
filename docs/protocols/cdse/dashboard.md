# The dashboard

```
python -m qd_workflows.network_sim --network DIR --thermo thermo_merged.json --params TEMPLATE --port 8777
```

Run this from branch `claude/fervent-greider-938b1f`. Ivan: this is the page to use. Acetate:
`--network network_growth_data/cdse_acetate_p2 --thermo …/thermo_acetate_pool.json --params
src/qd_workflows/templates/growth/cdse_zb_acetate.yaml`.

| card | shows | read it as |
|---|---|---|
| Conditions | T, ε, $c_\mathrm{P,0}$, $c_\mathrm{MX,0}$, $k_\mathrm{conv}$, $E_a$, time, $s_\nu$, η, site area, wall limit, network-edge stop; island barriers on or off; Run and Stop; the time at which to look at Ω | the inputs of [the master equation](master-equation.md) §7 |
| 1 · Nucleation | the grid's ΔΩ(m) up to 13 with ln S, ΔΩ*, m*, k*; the P(m, k) map; barrier and pool along t; ΔΩ(m, t); the profile from the pool through the grid to the dots, with the run's −kT ln N(m) | [the grid](grid.md) §5; the profile as in [the master equation](master-equation.md) §8.3 |
| 2 · Growth | time bar with play; α, CdSe in solution, dots, mean diameter, ⟨m⟩ ± σ, share on the network's edge; pool, nucleation, diameter and size-distribution charts | the phases of [Results](results.md) |
| 3 · Facets | CdSe added per facet family; the island-barrier table per family and size bin | which facets grow or give material back |
| 4 · Map | one cell per (m, k) coloured by Ω; populations as circles; the steps with the most net flux | where the dots are, where they go (blue ahead), what dissolves (populated red) |
| 5 · Structure | the 8 most populated states at the bar's time; space-filling view (Cd orange, Se blue, ligands pale); .xyz downloads | the actual structures |
| 6 · Growth movie | the dots over time; a growth path (widest integrated net transfer from the nuclei); a dissolution path; .webm and .xyz export | the routes material took |

The diameter is $d=(6mV/\pi)^{1/3}$ with V = 59.667 Å³ per CdSe (zinc blende, a = 6.203 Å). Grid structures come from
the thermo folder; network structures come from `geometry.npz`.

Run the growth itself on a cluster for long cases. The laptop serves the page and runs short cases only.
