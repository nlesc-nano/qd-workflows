# Generalised protocol and gaps

## 1. Step by step for a new material and ligand

Each stage lists its inputs, outputs, code, where to run it and the checks before moving on.

**0. Decide the chemistry.** Cation M, anion A, the X-type ligand and its Z-type salt MX_q; the precursor P; solvent
ε; the CIF; the builder recipe (facet families, passivation, prepass). Use the simplest ligand for site energetics
(acetate for carboxylates); a long chain enters only as solubility and steric coverage.

**1. Start dot and seeds** (`nucleation start`, `nucleation seeds`; laptop, builder only).
Out: the smallest neutral dot exposing every facet family; seed_full and seed_x (X-type only). Check: the X-type seed
is the builder's first Wulff dot (CdSe: (13, 3)).

**2. References** (`ligands refs` for carboxylates; campaign monomers otherwise; one GPU).
Out: the MA and MX_q monomers with Hessians. Check: lowest conformer, imaginary modes.

**3. Grid** (`nucleation run WD --model M --refs WD/refs`; one LUMI node or one A100 in 1-h pieces).
Settings to start from: seeded scheme, screen_steps 0, capacity, anneal_all_m = top, 3 starts, 4 ps from 700/900 K
(8 ps and 1100 K for m ≤ 6), beam 3, δ 0.05 eV, polish 0.005 eV/Å within 0.3 eV. Out: `grid.json`, Hessians,
`thermo_merged.json`. Checks: every row reaches k = 0 and its $k_\mathrm{max}$; no skipped node; imaginary modes per
state < 3; row-anneal gains settle. Cost: CdSe/CdCl2 ~4 A100-h; CdSe/acetate ~41 LUMI node-h with D3.

**4. COSMO of free species** (CP2K, one whole CPU node, all cores as MPI ranks, 100-cycle SCF restarts).
Out: E_inf, S for the monomers and the m = 1 complexes, from this grid's own geometries.

**5. Simulator check** (`simulator`, `grid_thermo`): ΔΩ(m), ΔΩ*, m*, ln S along the conversion at the target
conditions. Check: a barrier exists at realistic supersaturation.

**6. Wulff scan and lineage** (QD_builder `scan.py` on a CPU cluster; lineage script). Out: unique dots, the lineage tree
from the first Wulff dot. Check: the first lineage node is the grid's hand-over dot.

**7. Layer states and registry** (`python -m qd_workflows.layer_path size-edges data.json`, then `python -m qd_workflows.layer_path treepaths data.json net --mmax MMAX
--smax SMAX --cif CdSe_zb.cif` in waves of growing MMAX).
Out: `registry.json`, edge files, builds. Checks: failed edges per reason; the seed in the registry.

**8. Fill and shells** where the low-Ω band is sparse (`layer_path fill-enumerate`, `fill-links`, `fill-build`, `fill-register`):
one-layer neighbours of low-lying parents, Cd shells k0 ± 1, ± 2 by `z_type_target`. Plan the count and cost first.

**9. Islands** (`layer_path islands` on ladder parents; `layer_path stage-islands`). Cover every facet family and parent sizes up to
the largest m the growth will reach (CdSe stopped at m 150, which left the m 150–240 range without islands).

**10. Relax** (`layer_path stage`, `layer_path relax` with `templates/slurm/layer_path_relax_lumi.sh`, then
`layer_path collect` and `gcache`; or `p1.py` for a dressed ligand; GPU batches, 8 per LUMI node or 4 per A100 node).
No anneal above the grid. FD Hessians on a training set spanning the m range (CdSe: m ≤ 60 plus part of the large
sizes), learned F_vib for the rest. Checks: dressing failures, convergence, the F_vib model's domain flags.

**11. Steps, ΔG, ladders, export** (`layer_path layer-steps`, or reuse the steps of another ligand on the same
registry; `layer_path attach-dg`; `layer_path ladders`; `network_export`; one CPU node,
`templates/slurm/layer_path_steps_lumi.sh`). Checks: state count of the seed's component, σ values,
G-cache consistency (< 2 meV), island curves per family.

**12. Growth** (`network_growth run --network OUT --thermo thermo_merged.json --params TEMPLATE`, or the page).
Checks: detailed-balance residual (`network_growth check`); stop reason; the share of dots on the network's edge; the
nucleation time, burst length and conversion against experiments for that system. Tune k_conv, activities and s_ν to
experiments; do not rerun atomistic stages for exactness.

**13. Extend the map where the run stops** (the edge share > 0.5 or the dots piling where states thin out), then
re-export and rerun.

Before step 1 for a new system: commit the grid and ligand code listed in the [overview](index.md) as one module, so
every stage runs from one tagged version (the layer-path code is the module `qd_workflows.layer_path`). The layer
path is CdSe zinc-blende only: the lattice constant, the anion species and the builder preset (Cl ligand, Cd²⁺,
Se²⁻) are fixed in `layer_path/lattice.py` and `build.py`.


## 2. Current state and gaps (2026-11-04)

- **Grids**: CdCl2 `cdse_full2`; acetate `wd_full_acetate_pool` (definitive).
- **Networks**: CdCl2 tree + fill (15,160 states); acetate P1 + P2 (15,131 states, m 7–223).
- **On hold (Ivan, 2026-11-04)**: the acetate fill, A = wave 4 up to m 239 + 31 fill dots + islands (~470 LUMI
  GPU-h), B = m 240–405 (~145 GPU-h). The FD Hessian check C was dropped. Plan page:
  [Acetate Network Fill Plan](https://claude.ai/artifact/VXmjtdWMhq9NT8NgeeiqXw).
- **Gaps in the map**: wave 4 (11,311 dots, m 79–405) built but not relaxed for acetate; no islands on parents above
  m 150; the registry stops at m 405.
- **Model gaps**:
  - (fixed 2026-11-04, commit 5af3e3f) the conversion now runs on the salt left in solution; the open-loop form had
    produced a spurious second burst and a stall at t ≈ 356 s in the acetate default case;
  - layer steps that also change k carry only CdSe in their rate;
  - the acetate complexes have no COSMO of their own (they take the CdSe monomer radius).
- **Hand-over**: the network's relaxed (13, 3) is 0.32 eV above the grid's annealed node, and network dots with m ≤ 13
  are 1.24 eV above the grid on average; the model replaces them by the grid nodes. Whether to seed the network from
  the grid's (13, 3) geometry was never decided.
- **Not found in the sources**: the script that wrote the lineage; the rule that chose the 32 island parents; the
  script that cut the P1/P2 lists and chose their islands.
- **MACE distillation**: a plan to make large-m relaxations cheaper is being discussed in its own session.
