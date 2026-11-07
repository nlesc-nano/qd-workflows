# Sources

Where every piece of code, data and run lives.

| what | where |
|---|---|
| grid construction code (CdCl2) | `src/qd_workflows/nucleation.py`, `lattice_grid.py`, `grid.py`, `templates/nucleation/cdse_zb.yaml` (`python -m qd_workflows.nucleation`; QD_builder main 330fc4b or later) |
| carboxylate grid code | `src/qd_workflows/ligands.py` (ligand block of `nucleation.json`), `grid_thermo.py` (partial Hessian); scripts `pool_grids.py`, `complete_grid.py`, `plan.html` in worktree `cdse-carboxylate-plan-501d58` `docs/plans/cdse_carboxylate/pilot/` |
| thermo, solvation, pool | `src/qd_workflows/grid_thermo.py`, `grid_solvation.py`, `reaction/solvation.py`; docs `theory/thermochemistry.md`, `theory/solvation.md` |
| CdCl2 grid run | `cdse-grid/cd16_lattice_grid/nucleation_builder_main/full_run2/` (`thermo_merged.json`, `grid/grid.json`, `complexes/`) |
| acetate grid run | LUMI `/scratch/project_465003437/qdw/wd_full_acetate_pool/`; local `thermo_acetate_pool.json` in `network_growth_data/cdse_acetate_p2/` |
| QD_builder | `~/Documents/University/escience/QD_builder`, pinned at main 330fc4b (`components.yaml`), which contains 6561a07 (fix/followups), the builder of the layer and island states; `runs/wulff_scan_cdse/` (scan, `data.json`) |
| layer-path code | `src/qd_workflows/layer_path/` (`python -m qd_workflows.layer_path <command>`: lattice, builds, treepaths, relax staging, layer steps, ladders, islands, barriers, fill); LUMI job templates `templates/slurm/layer_path_relax_lumi.sh`, `layer_path_steps_lumi.sh`; tests `tests/test_layer_path.py`. The original scripts stay in `network_growth_data/pipeline_scripts_layer_path/` for reference (one-off diagnostics there were not moved) |
| registry, builds | growth-models session folder `overlay/net/` (registry 711 MB, 40,048 build folders); NHR `$WORK/qdw_tests/acetate_net_export/net/registry.json`; LUMI `layer_path_tree/` (waves, `isl/`) |
| island builds | `network_growth_data/island_builds_cdcl2/` (`states_index.json`, `islands.json`, `builds/`) |
| acetate network runs | NHR `$WORK/qdw_tests/acetate_net_p{0,1,2}_nhr/`, LUMI `qdw/acetate_net_p2/` |
| export | NHR `$WORK/qdw_tests/acetate_net_export/` (`ac_build.py`, `attach_dg.py`, `ladders.py`, `job_net_p2.sh`, `network_p2/`); `src/qd_workflows/network_export.py` |
| learned F_vib | branch `claude/ml-dot-energy-2af8ae` (fvib_model); NHR `$WORK/qdw_tests/acetate_fvib/` |
| growth model and page | branch `claude/fervent-greider-938b1f`: `network_growth.py`, `network_sim.py`, `templates/network/index.html`, `templates/growth/cdse_zb*.yaml` |
| older size law | branch `feat/growth-unified`: `growth.py`, `facet_kinetics.py`, `builder_geometry.py`, `sitekin/` |
| exported networks | `network_growth_data/cdse_cdcl2_tree_fill/`, `cdse_acetate_tree/`, `cdse_acetate_p2/` |
| plans | [Acetate Layer-Path Network](https://claude.ai/artifact/6VCgDSB4A5az6HN1YuqBkc), [Acetate Network Fill Plan](https://claude.ai/artifact/VXmjtdWMhq9NT8NgeeiqXw), [qdw Nucleation Recipe](https://claude.ai/artifact/6AjayTE3CL3mzJHvbdFu1a) |
| figures of this section | `docs/protocols/cdse/figure_scripts/` (`render.py`: structures; `run_fig.py`: run phases; `map_fig.py`: Ω maps, run from the `claude/fervent-greider-938b1f` worktree with `PYTHONPATH=src`) |
