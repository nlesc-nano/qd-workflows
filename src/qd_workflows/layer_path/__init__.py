"""Layer path of the CdSe growth network: the layer states and islands between the Wulff dots of the lineage, built
with QD_builder, and the registry that network_export turns into the network of network_growth.

    python -m qd_workflows.layer_path <command> ...

Pipeline (docs/protocols/cdse/network.md, islands.md, protocol.md):
  size-edges        every lineage edge: embedding, plane moves, number of layer states       -> edges_size.json
  treepaths         layer states of the edges (waves of growing --mmax), one build per class -> net/registry.json,
                                                                                               net/edges, net/builds
  stage / relax / collect   LUMI batches, relax + FD Hessian, results by dot key           -> results_by_key.json
  gcache            G at the template's T per relaxed dot                                    -> g_cache.json
  layer-steps, attach-dg    what each layer step adds, dG                                    -> layer_steps.json
  ladders, ladder-fit, chains   island ladders already in the registry and their fits
  islands, stage-islands    step-bounded islands on the facets of chosen parents             -> islands.json,
                                                                                               states_index.json
  island-barriers, island-cdse  island curves and barriers
  fill-enumerate, fill-links, fill-build, fill-register   the m 20-50 fill of one-layer neighbours
then python -m qd_workflows.network_export. One edge in detail: overlay, layers, lattice-check, distinct, path, graph,
unique (seed_child).

Modules: lattice (sites, planes, class and dot keys), build (per-plane QD_builder builds), seed_child, treepaths,
relax, steps, islands, barriers, fill. Builds need QD_builder with fix/followups (6561a07). Job scripts:
templates/slurm/layer_path_relax_lumi.sh, layer_path_steps_lumi.sh.
"""
