"""Command line of the layer path: python -m qd_workflows.layer_path <command> ... (see __init__ for the steps)."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def _data(p):
    return json.loads(Path(p).read_text())


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python -m qd_workflows.layer_path")
    sub = ap.add_subparsers(dest="cmd", required=True)

    def cmd(name, help_):
        return sub.add_parser(name, help=help_)

    def nproc(p, default=4):
        p.add_argument("--nproc", type=int, default=default)

    # one lineage edge in detail (seed_child)
    p = cmd("overlay", "embeddings of the parent in each child and the planes they move")
    p.add_argument("data")
    p.add_argument("children", nargs="*")
    p.add_argument("--parent", help="parent label (default: lineage node 0, the seed)")
    p = cmd("layers", "one edge as whole-layer moves of {100}/{111} planes (layers_<child>.json)")
    p.add_argument("data")
    p.add_argument("child")
    p.add_argument("--embedding", type=int, default=0)
    p.add_argument("--out")
    p.add_argument("--parent")
    for name, help_ in (("lattice-check", "layer frameworks vs the lattice cut by their planes"),
                        ("distinct", "symmetry-distinct Cd/Se frameworks of the layer states")):
        p = cmd(name, help_)
        p.add_argument("data")
        p.add_argument("child")
        p.add_argument("--layers", help="layers_<child>.json")
        p.add_argument("--parent")
    p = cmd("path", "every distinct layer state of one edge built with QD_builder (states.json)")
    p.add_argument("data")
    p.add_argument("child")
    p.add_argument("out")
    p.add_argument("--cif", required=True)
    p.add_argument("--parent")
    p = cmd("graph", "classes of the layer states and the edges one layer apart (graph_<child>.json)")
    p.add_argument("data")
    p.add_argument("child")
    p.add_argument("--out")
    p.add_argument("--parent")
    p = cmd("unique", "unique built dots over several path folders")
    p.add_argument("folders", nargs="+")
    p.add_argument("--out", default="unique.json")

    # every lineage edge: the registry (treepaths)
    p = cmd("size-edges", "embedding and number of layer states of every lineage edge (edges_size.json)")
    p.add_argument("data")
    p.add_argument("--out", default="edges_size.json")
    nproc(p)
    p = cmd("treepaths", "layer paths of the lineage edges into net/registry.json, net/edges, net/builds")
    p.add_argument("data")
    p.add_argument("net")
    p.add_argument("--mmax", type=int, required=True, help="largest child m (N_Se) of this wave")
    p.add_argument("--smax", type=int, required=True, help="largest number of layer states of an edge")
    p.add_argument("--cif", required=True)
    p.add_argument("--edges-size", default="edges_size.json")
    nproc(p)

    # relaxation and free energies (relax)
    p = cmd("stage", "registry dots without results into balanced batches WAVE/batchJ/d<id>.xyz")
    p.add_argument("net")
    p.add_argument("results")
    p.add_argument("wave")
    p.add_argument("n_batches", type=int)
    p.add_argument("--staged", default="staged.json")
    p = cmd("relax", "relax + polish + FD Hessian of a batch folder (on the cluster, Orchestr.AI)")
    p.add_argument("batch")
    p.add_argument("model")
    p.add_argument("out")
    p.add_argument("--pattern", default="*.xyz")
    p = cmd("collect", "relaxed batches into results_by_key.json")
    p.add_argument("net")
    p.add_argument("results")
    p.add_argument("states", nargs="+")
    p = cmd("gcache", "G per relaxed registry dot at the template's T into g_cache.json")
    p.add_argument("net")
    p.add_argument("results")
    p.add_argument("g_cache")
    p.add_argument("--thermo", required=True, help="the grid's thermo_merged.json (temperatures)")
    p.add_argument("--params", required=True, help="growth template (conditions.T, solvation)")
    p.add_argument("--check", type=int, default=0)

    # layer steps and ladders (steps)
    p = cmd("layer-steps", "what each registry layer step adds (layer_steps.json)")
    p.add_argument("net")
    p.add_argument("results")
    p.add_argument("g_cache")
    p.add_argument("out")
    nproc(p)
    p = cmd("attach-dg", "re-attach dG and m_a to a layer_steps.json")
    p.add_argument("net")
    p.add_argument("results")
    p.add_argument("g_cache")
    p.add_argument("steps")
    p = cmd("ladders", "island ladders of the one-sheet steps and the per-family regression")
    p.add_argument("steps")
    p.add_argument("out")
    p = cmd("ladder-fit", "classical island fit of every ladder at the reference pools")
    p.add_argument("ladders")
    p.add_argument("out")
    p = cmd("chains", "same-sheet chains of steps")
    p.add_argument("net")
    p.add_argument("steps")
    p.add_argument("out")
    p.add_argument("--max-depth", type=int, default=12)
    nproc(p)

    # islands
    p = cmd("islands", "island states of the parents built with QD_builder (OUT/islands.json)")
    p.add_argument("net")
    p.add_argument("parents", help="JSON [[dot_id, family], ...]")
    p.add_argument("out")
    p.add_argument("--cif", required=True)
    p.add_argument("--per", type=int, default=2, help="islands per size and normal")
    p.add_argument("--limit", type=int, default=0, help="build only this many island states (smoke test)")
    nproc(p)
    p = cmd("stage-islands", "unique pure Q = 0 islands into batches; ISLANDS/states_index.json")
    p.add_argument("net")
    p.add_argument("islands")
    p.add_argument("results")
    p.add_argument("wave")
    p.add_argument("n_batches", type=int)
    p = cmd("island-barriers", "island curves and barriers per parent and facet (and island_states.json)")
    p.add_argument("net")
    p.add_argument("islands")
    p.add_argument("results")
    p.add_argument("g_cache")
    p.add_argument("ladders")
    p.add_argument("ladder_fit")
    p.add_argument("out")
    p.add_argument("--states", nargs="+", required=True, help="the relaxed island batches' states.json")
    p.add_argument("--thermo", required=True)
    p.add_argument("--params", required=True)
    p.add_argument("--types", help="also write the island states split by type (island_states.json)")
    p = cmd("island-cdse", "CdSe-growth barrier with the shell equilibrated at each size")
    p.add_argument("island_states")
    p.add_argument("ladders")
    p.add_argument("ladder_fit")
    p.add_argument("out")

    # m 20-50 fill
    p = cmd("fill-enumerate", "one-layer neighbours of the network's low-lying states (candidates.json)")
    p.add_argument("out")
    p.add_argument("--network", required=True, help="network_export folder")
    p.add_argument("--registry", required=True)
    p.add_argument("--thermo", required=True)
    p.add_argument("--params", required=True)
    p.add_argument("--m-lo", type=int, default=14)
    p.add_argument("--m-hi", type=int, default=55)
    p.add_argument("--dw", type=float, default=0.4)
    p.add_argument("--cap", type=int, default=6)
    p = cmd("fill-links", "registry dots and candidates one layer from each candidate")
    p.add_argument("candidates")
    p.add_argument("registry")
    p.add_argument("out")
    p = cmd("fill-build", "the candidates and their shells built with QD_builder (OUT/fill.json)")
    p.add_argument("candidates")
    p.add_argument("registry")
    p.add_argument("out")
    p.add_argument("--cif", required=True)
    nproc(p)
    p = cmd("fill-register", "the kept fill dots into the registry and net/edges/fill2050.json")
    p.add_argument("net")
    p.add_argument("fill_dir")
    p.add_argument("candidates")
    p.add_argument("links")
    p.add_argument("keep", help="JSON list of the fill names to keep")
    p.add_argument("--edges-name", default="fill2050")

    a = ap.parse_args(argv)
    from . import barriers, fill, islands, relax, seed_child, steps, treepaths
    c = a.cmd
    if c == "overlay":
        seed_child.overlay(_data(a.data), a.children, a.parent)
    elif c == "layers":
        seed_child.layers(_data(a.data), a.child, a.out, a.embedding, a.parent)
    elif c == "lattice-check":
        seed_child.lattice_check(_data(a.data), a.child, a.layers, a.parent)
    elif c == "distinct":
        seed_child.distinct(_data(a.data), a.child, a.layers, a.parent)
    elif c == "path":
        seed_child.path(_data(a.data), a.child, a.out, a.cif, a.parent)
    elif c == "graph":
        seed_child.graph(_data(a.data), a.child, a.out, a.parent)
    elif c == "unique":
        seed_child.unique(a.folders, a.out)
    elif c == "size-edges":
        treepaths.size_edges(a.data, a.out, a.nproc)
    elif c == "treepaths":
        treepaths.treepaths(a.data, a.net, a.mmax, a.smax, a.cif, a.edges_size, a.nproc)
    elif c == "stage":
        relax.stage(a.net, a.results, a.wave, a.n_batches, a.staged)
    elif c == "relax":
        relax.relax_states(a.batch, a.model, a.out, a.pattern)
    elif c == "collect":
        relax.collect(a.net, a.results, a.states)
    elif c == "gcache":
        relax.gcache(a.net, a.results, a.g_cache, a.thermo, a.params, a.check)
    elif c == "layer-steps":
        steps.layer_steps(a.net, a.results, a.g_cache, a.out, a.nproc)
    elif c == "attach-dg":
        steps.attach_dg(a.net, a.results, a.g_cache, a.steps)
    elif c == "ladders":
        steps.ladders(a.steps, a.out)
    elif c == "ladder-fit":
        steps.ladder_fit(a.ladders, a.out)
    elif c == "chains":
        steps.chains(a.net, a.steps, a.out, a.nproc, a.max_depth)
    elif c == "islands":
        islands.islands(a.net, _data(a.parents), a.out, a.cif, a.per, a.nproc, limit=a.limit)
    elif c == "stage-islands":
        islands.stage(a.net, a.islands, a.results, a.wave, a.n_batches)
    elif c == "island-barriers":
        isl = barriers.island_states(a.net, a.islands, a.results, a.g_cache, a.thermo, a.params, a.states)
        barriers.island_barriers(isl, a.ladders, a.ladder_fit, a.out)
        if a.types:
            barriers.island_types(isl, a.types)
    elif c == "island-cdse":
        barriers.island_cdse(a.island_states, a.ladders, a.ladder_fit, a.out)
    elif c == "fill-enumerate":
        par = fill.low_lying(a.network, a.thermo, a.params, a.m_lo, a.m_hi, a.dw, a.cap)
        fill.enumerate_candidates(a.registry, par, a.out)
    elif c == "fill-links":
        fill.links(a.candidates, a.registry, a.out)
    elif c == "fill-build":
        fill.build_all(a.candidates, a.registry, a.out, a.cif, a.nproc)
    elif c == "fill-register":
        fill.register(a.net, a.fill_dir, a.candidates, a.links, set(_data(a.keep)), a.edges_name)
    return 0


if __name__ == "__main__":
    sys.exit(main())
