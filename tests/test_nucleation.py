"""Nucleation workflow, the parts without MACE: facet site areas and interstitial centres from the CIF, the MD
chunks of the annealing stage, the merge of thermochemistry files and the padded per-m trajectories."""
import json
from pathlib import Path

import numpy as np
import pytest

from qd_workflows import nucleation as nu

CIF = Path(__file__).parent / "data" / "cd16" / "CdSe_zb.cif"


def test_site_areas_and_holes_of_zinc_blende():
    from pymatgen.core import Structure
    s = Structure.from_file(str(CIF))
    a = s.lattice.a
    cat = next(x.specie.symbol for x in s.sites if x.specie.symbol == "Cd")
    assert np.isclose(nu.site_area(s, [1, 1, 1], cat), np.sqrt(3) * a * a / 4, rtol=1e-3)
    assert np.isclose(nu.site_area(s, [1, 0, 0], cat), a * a / 2, rtol=1e-3)
    assert np.isclose(nu.site_area(s, [1, 1, 0], cat), a * a / np.sqrt(2), rtol=1e-3)
    holes = nu.interstitial_centres(s)
    assert len(holes) == 2                                # empty tetrahedral and octahedral holes


def test_chunk_dirs_skip_the_index(tmp_path):
    for c in ("chunk001", "chunk000"):
        (tmp_path / "anneal" / c).mkdir(parents=True)
    (tmp_path / "anneal" / "chunks.json").write_text("{}")
    assert [d.name for d in nu.chunk_dirs(tmp_path)] == ["chunk000", "chunk001"]


def test_anneal_chunks_respect_schedule_and_size():
    P = nu._block({"nucleation": {"max_md_atoms": 300}})
    starts = [{"m": m, "k": k, "symbols": ["Cd"] * (10 + 3 * m)} for m in range(1, 13) for k in range(2)]
    chunks = nu.anneal_chunks(starts, P)
    seen = sorted(i for c in chunks for i in c["starts"])
    assert seen == list(range(len(starts)))               # every start once
    for c in chunks:
        small = {starts[i]["m"] <= P["anneal"]["small_m"] for i in c["starts"]}
        assert len(small) == 1                            # one schedule per chunk
        atoms = sum(len(starts[i]["symbols"]) * len(c["hot"]) for i in c["starts"])
        assert atoms <= 300 or len(c["starts"]) == 1


def _thermo(path, label_states):
    nodes = [{"m": m, "k": k, "formula": f, "ensemble_E_eV": [e for e, _s in st], "states":
              [{"E_eV": e, "file": f"nodes/m{m:02d}_k{k:02d}{'' if i == 0 else f'_{i}'}.xyz", "fingerprint": "x", **({"source": s} if s else {})}
               for i, (e, s) in enumerate(st)]} for (m, k, f), st in label_states.items()]
    path.write_text(json.dumps({"temperatures": [300], "vibrations": True, "references": {}, "nodes": nodes, "complete": True}))
    return path


def test_merge_keeps_sources_and_one_copy_per_minimum(tmp_path):
    for d in ("grid_lattice", "grid_relaxed", "anneal/grid"):
        (tmp_path / d).mkdir(parents=True)
    a = _thermo(tmp_path / "grid_lattice" / "thermo.json", {(2, 1, "Cd3Se2Cl2"): [(-10.0, None), (-9.5, None)]})
    b = _thermo(tmp_path / "grid_relaxed" / "thermo.json", {(2, 1, "Cd3Se2Cl2"): [(-10.0004, None), (-9.8, None)]})
    c = _thermo(tmp_path / "anneal" / "grid" / "thermo.json", {(2, 1, "Cd3Se2Cl2"): [(-10.3, "anneal_from_lattice_start")]})
    res = nu.merge_thermo(tmp_path / "m.json", [("lattice_start", a), ("relaxed_start", b), ("anneal", c)])
    st = res["nodes"][0]["states"]
    assert [round(s["E_eV"], 2) for s in st] == [-10.3, -10.0, -9.8, -9.5]
    assert st[0]["source"] == "anneal_from_lattice_start"
    assert st[0]["file"] == "anneal/grid/nodes/m02_k01.xyz"                     # where the file is, from m.json
    assert st[3]["file"] == "grid_lattice/nodes/m02_k01_1.xyz"
    assert st[1]["source"] == "lattice_start + relaxed_start"          # found by both starts: one minimum


def test_trajectories_have_one_atom_list(tmp_path):
    (tmp_path / "g" / "nodes").mkdir(parents=True)
    rng = np.random.default_rng(0)
    nodes = []
    for k in range(3):
        s = ["Cd"] * (2 + k) + ["Se"] * 2 + ["Cl"] * (2 * k)
        nu._write_xyz(tmp_path / "g" / "nodes" / f"m02_k{k:02d}.xyz", s, rng.normal(size=(len(s), 3)), "x")
        nodes.append({"m": 2, "k": k, "formula": "x", "states": [{"E_eV": -k, "source": "lattice_start", "file": f"g/nodes/m02_k{k:02d}.xyz"}]})
    nu.write_trajectories({"nodes": nodes}, tmp_path, tmp_path / "xyz")
    L = (tmp_path / "xyz" / "m02.xyz").read_text().splitlines()
    frames, i = [], 0
    while i < len(L) and L[i].strip():
        n = int(L[i])
        frames.append(([l.split()[0] for l in L[i + 2:i + 2 + n]], L[i + 1]))
        i += 2 + n
    assert len(frames) == 3 and len({tuple(s) for s, _c in frames}) == 1   # same elements in the same order
    assert "n_pad=6" in frames[0][1] and "n_pad=0" in frames[2][1]


def test_start_finder_picks_a_faceted_dot(tmp_path):
    pytest.importorskip("builder.scripts.generate_library")
    cfg = Path(__file__).parents[1] / "src" / "qd_workflows" / "templates" / "nucleation" / "cdse_zb.yaml"
    if not Path(__import__("yaml").safe_load(cfg.read_text())["cif"]).is_file():
        pytest.skip("QD_builder CIF not on this machine")
    r = nu.find_start(str(cfg), str(tmp_path), log=lambda _x: None)
    assert all(r["sites"][f][1] >= 3 for f in r["families"])
    assert (tmp_path / "start" / "record.json").is_file() and (tmp_path / "start" / "start.xyz").is_file()


# ---------------- screened search: screening, seeds, row-anneal selection ----------------
def test_omega_gaps_follow_the_reachable_mu_range():
    """Omega(k) = G(k) - k mu: with G linear in k (slope s) every node ties at mu = s; a node bent above the line by
    0.3 eV never comes within 0.3 eV of the lowest, and only the ends win outside mu = s."""
    temps = [300.0, 500.0, 700.0]
    s = -2.0
    G = {k: [s * k + (0.3 if k == 2 else 0.0)] * 3 for k in range(5)}
    at_s = nu.omega_gaps(G, [s] * 3, temps, (300, 700), (1.0, 1.0))           # mu = s exactly
    assert at_s[0] == pytest.approx(0) and at_s[4] == pytest.approx(0) and at_s[2] == pytest.approx(0.3)
    lo = nu.omega_gaps(G, [s - 1.0] * 3, temps, (300, 700), (1.0, 1.0))      # mu far below s: k = 0 wins
    assert lo[0] == pytest.approx(0) and lo[4] == pytest.approx(4.0)
    wide = nu.omega_gaps(G, [s] * 3, temps, (300, 700), (1e-6, 1.0))        # lower c only lowers mu
    assert wide[4] == pytest.approx(0) and wide[2] == pytest.approx(0.3, abs=1e-6)


def test_single_bond_cations_count_z_type_units():
    d = 2.6
    syms = ["Se", "Cd", "Cd", "Se", "Cd"]
    pts = [[0, 0, 0], [d, 0, 0], [-d, 0, 0], [2 * d, 0, 0], [0, d, 0]]     # Cd 1 bridges two Se; Cd 2, Cd 4: one Se each
    assert nu.single_bond_cations(syms, pts, "Cd", "Se", d) == 2


def _record(tmp_path):
    """The CdSe test start (Cd28Se13Cl30, 1.25 cells, Se centre) as a record folder."""
    pytest.importorskip("builder.scripts.generate_library")
    import copy
    import yaml
    from builder.scripts.generate_library import _build
    cfg = yaml.safe_load((Path(__file__).parents[1] / "src" / "qd_workflows" / "templates" / "nucleation" /
                          "cdse_zb.yaml").read_text())
    if not Path(cfg["cif"]).is_file():
        pytest.skip("QD_builder CIF not on this machine")
    recipe = copy.deepcopy(cfg["preset"]) | {"size_unit_cells": [1.25] * 3, "construction_origin": {"center_on_species": "Se"}}
    rec = tmp_path / "start"
    rec.mkdir()
    res = _build(cfg["cif"], recipe, "add", tmp_path / "b")
    nu._write_xyz(rec / "start.xyz", list(res["symbols"]), res["pts"], "start")
    (rec / "CdSe_zb.cif").write_text(Path(cfg["cif"]).read_text())
    (rec / "record.json").write_text(json.dumps({"material": "CdSe", "origin": {"cif": "CdSe_zb.cif", "recipe": recipe,
                                                                                "positive_q_mode": "add"}}))
    return rec


def test_seeds_are_both_ends_of_the_mx_range(tmp_path):
    """Invariants of the two seeds for whatever dot the builder cuts: the full seed is the start's composition (no free
    host site on it), the X-type seed has k_X = k_full - the single-bond cations and none of them left, same core. On the
    CdSe start this pilot used (Cd28Se13Cl30) the X-type seed is the symmetric Cd16Se13Cl6."""
    from collections import Counter
    from scipy.spatial import cKDTree
    from builder.z_type_target import run_z_type_target  # noqa: F401  (a QD_builder with z_type_target)
    rec = _record(tmp_path)
    s0, p0 = nu._read_xyz(rec / "start.xyz")
    res = nu.build_seeds(rec, log=lambda _x: None)
    by = {x["ladder"]: x for x in res["seeds"]}
    sf, pf = nu._read_xyz(rec / nu.SEED_FULL)
    sx, px = nu._read_xyz(rec / nu.SEED_X)
    d = 2.632
    n1 = nu.single_bond_cations(sf, pf, "Cd", "Se", d)
    assert (by["down"]["m"], by["down"]["k"]) == (13, s0.count("Cl") // 2)                 # the start's own k
    assert (by["up"]["m"], by["up"]["k"]) == (13, by["down"]["k"] - n1)
    assert by["up"]["ledger"] == {"method": "strip_z_units", "units_removed": n1, "k_full": by["down"]["k"]}
    assert nu.single_bond_cations(sx, px, "Cd", "Se", d) == 0
    assert 2 * sx.count("Cd") == 2 * sx.count("Se") + sx.count("Cl")                       # neutral
    if Counter(s0) == Counter({"Cd": 28, "Se": 13, "Cl": 30}):
        t = cKDTree(px)
        cn = Counter((sum(sx[j] == "Se" for j in t.query_ball_point(px[i], 2.9)),
                      sum(sx[j] == "Cl" for j in t.query_ball_point(px[i], 2.9))) for i, e in enumerate(sx) if e == "Cd")
        assert cn == {(4, 0): 4, (2, 1): 12}              # the symmetric X-type dot: one bridging Cl per edge Cd


def test_strip_z_units_takes_the_orphans_then_the_most_crowded_x():
    """Two Se bridged by a Cd; on each Se a Z-type Cd (one Se) with its terminal Cl; an extra Cl bridges the core Cd and
    both Z-type Cd, and one more Cl sits on the core Cd only. Stripping both Z units (q = 2) takes the 3 Cl bonded only
    to them and one of the two equally crowded Cl on the core Cd, leaving one Cl there. A Z-type Cd holding more terminal
    Cl than q cannot be stripped (None)."""
    d = 2.6
    syms = ["Se", "Cd", "Se", "Cd", "Cl", "Cl", "Cd", "Cl", "Cl", "Cl"]
    pts = [[0, 0, 0], [d, 0, 0], [2 * d, 0, 0],               # Se - Cd - Se
           [0, d, 0], [0, 2 * d, 0], [-d, d, 0],             # Z Cd on Se 0, two terminal Cl
           [2 * d, d, 0], [2 * d, 2 * d, 0],                 # Z Cd on Se 2, one terminal Cl
           [d, d, 0],                                        # Cl bridging the core Cd and both Z Cd
           [d, -d, 0]]                                       # Cl on the core Cd only
    s, p, n = nu.strip_z_units(syms, pts, "Cd", "Se", "Cl", 2, d)
    assert n == 2 and (s.count("Cd"), s.count("Se"), s.count("Cl")) == (1, 2, 1)
    three = ["Se", "Cd", "Se", "Cd", "Cl", "Cl", "Cl"]   # Z Cd on Se 0 with three terminal Cl
    p3 = pts[:6] + [[0, d, d]]
    assert nu.strip_z_units(three, p3, "Cd", "Se", "Cl", 2, d) is None


class _Atoms:
    def __init__(self, syms, pos):
        self.syms, self.pos = syms, np.asarray(pos, float)

    def get_positions(self):
        return self.pos

    def get_chemical_symbols(self):
        return self.syms


def _fake_grid(monkeypatch, screen):
    """A LatticeGrid without MACE: energy = first coordinate of the first atom (the screened energy the same, plus a
    per-candidate offset), so which candidates get fully relaxed can be read off the relaxation calls."""
    import sys
    import types
    from qd_workflows import lattice_grid as lg
    rec = types.ModuleType("orchestr_ai.qd.records")
    rec.fingerprint = lambda s, p: f"{len(s)}:{np.asarray(p, float)[0, 0]:.3f}"
    monkeypatch.setitem(sys.modules, "orchestr_ai.qd.records", rec)
    calls = []

    class Cache:
        n_new = 0

        def relax_many(self, structures, settings, polish=False):
            calls.append((polish, len(structures)))
            return [(_Atoms(list(s), p), float(np.asarray(p)[0, 0]), True) for s, p in structures]

    g = object.__new__(lg.LatticeGrid)
    g.dt, g.s, g.beam, g.delta, g.log = types.SimpleNamespace(POLISH_WINDOW=0.3), None, 3, 0.05, lambda _x: None
    g.start, g.screen, g.screen_keep, g.n_screened = "lattice", screen, 3, 0
    g.capacity, g._sites, g.seed_nodes = False, {}, set()
    g.cache, g.nodes, g.timings, g.t0 = Cache(), {}, [], 0.0
    g.lat = types.SimpleNamespace(anion="Se", ligand="Cl", q=2, cation="Cd")
    g._save_cache_throttled = lambda force=False: None
    g.save = lambda final=False: None
    g._screen = lambda items: [(np.asarray(p, float), float(np.asarray(p)[0, 0])) for _s, p in items]
    return g, calls


def _moves(n):
    syms = ["Cd", "Se", "Cl", "Cl"]                       # m = 1, k = 1

    def moves(_syms, _pts):
        return [({"i": i}, syms, np.array([[float(i), 0, 0], [3, 0, 0], [0, 3, 0], [0, 0, 3]])) for i in range(n)]
    return moves


def test_screening_relaxes_only_the_best_few_per_node(monkeypatch):
    from qd_workflows import lattice_grid as lg
    g, calls = _fake_grid(monkeypatch, screen=20)
    parent = lg.Config_(["Cd"], np.zeros((1, 3)), np.zeros((1, 3)), 0.0, "p", {}, None)
    kept = g.step([parent], _moves(10), "+CdCl2", (1, 1), path="seed_x.xyz")
    assert calls[0] == (False, 3)                         # 10 candidates screened, the best 3 fully relaxed
    assert [x.energy for x in kept] == [0.0]              # within delta of the best: one
    assert all(x.path == "seed_x.xyz" for x in kept)
    rec = g.timings[-1]
    assert (rec["candidates"], rec["relaxed"]) == (10, 3)
    assert sorted(g.nodes[(1, 1)]["ensemble"].values()) == [0.0, 1.0, 2.0]   # screened-out ones are not minima
    g2, calls2 = _fake_grid(monkeypatch, screen=0)
    g2.step([parent], _moves(10), "+CdCl2", (1, 1))
    assert calls2[0] == (False, 10)                       # exhaustive (legacy): every candidate fully relaxed


def test_each_seed_path_keeps_its_best_state(monkeypatch):
    """Two seeds reaching one node: the node keeps both paths' bests even when one is far above the other (beyond
    delta), and each ladder goes on from its own products."""
    from qd_workflows import lattice_grid as lg
    g, _calls = _fake_grid(monkeypatch, screen=0)
    parent = lg.Config_(["Cd"], np.zeros((1, 3)), np.zeros((1, 3)), 0.0, "p", {}, None)
    syms = ["Cd", "Se", "Cl", "Cl"]
    hi = lambda _s, _p: [({"i": 0}, syms, np.array([[5.0, 0, 0], [3, 0, 0], [0, 3, 0], [0, 0, 3]]))]
    a = g.step([parent], _moves(2), "+CdCl2", (1, 1), path="seed_x.xyz")
    b = g.step([parent], hi, "-CdCl2", (1, 1), path="seed_full.xyz")
    assert [x.energy for x in b] == [5.0] and [x.energy for x in a] == [0.0]
    assert sorted((x.energy, x.path) for x in g.nodes[(1, 1)]["states"]) == [(0.0, "seed_x.xyz"), (5.0, "seed_full.xyz")]
    g_old, _c = _fake_grid(monkeypatch, screen=0)          # legacy (no path): the 5 eV state is dropped
    g_old.step([parent], _moves(2), "+CdCl2", (1, 1))
    assert [x.energy for x in g_old.step([parent], hi, "-CdCl2", (1, 1))] == [0.0]


def test_prune_and_shortlist():
    from qd_workflows import lattice_grid as lg
    st = [lg.Config_([], None, None, e, str(e), {}, p) for e, p in [(0.0, None), (0.01, None), (0.02, None),
                                                                     (0.03, None), (0.9, "anneal")]]
    assert [x.energy for x in lg.prune(st, 0.05, 3)] == [0.0, 0.01, 0.9]
    assert lg.shortlist({(1, 1): [("a", 3.0), ("b", 1.0), ("c", 2.0)], (1, 2): [("d", 0.0)]}, 2) == \
        {(1, 1): ["b", "c"], (1, 2): ["d"]}


def test_anneal_parents_take_each_seed_path_then_repeat():
    from qd_workflows import lattice_grid as lg
    st = lambda e, fp, path: lg.Config_([], None, None, e, fp, {}, path)
    a, b, c = st(-10.0, "a", "seed_x.xyz"), st(-9.8, "b", "seed_full.xyz"), st(-9.99, "c", "seed_x.xyz")
    got = nu.anneal_parents([c, b, a], 2)
    assert [(x.fp, l) for x, l in got] == [("a", "best"), ("b", "path:seed_full.xyz")]   # the other path, not c
    got = nu.anneal_parents([a, c], 3)
    assert [(x.fp, l) for x, l in got] == [("a", "best"), ("c", "path:seed_x.xyz"), ("a", "best, repeat 2")]
    assert [x.fp for x, _l in nu.anneal_parents([a], 2)] == ["a", "a"]                   # one state: new velocities


def test_row_checkpoint_resumes_without_the_seeds(monkeypatch, tmp_path):
    """After a finished row the grid is checkpointed; a rerun restores the nodes (states, paths, seed nodes) and goes
    straight to the descent from that row, without relaxing or laddering the seeds again."""
    from qd_workflows import lattice_grid as lg
    g, _calls = _fake_grid(monkeypatch, screen=0)
    g.out = tmp_path
    rel = np.arange(12.0).reshape(4, 3)
    g.nodes = {(13, 3): {"states": [lg.Config_(["Cd", "Se", "Cl", "Cl"], rel + 0.5, rel, -89.3, "fx",
                                               {"move": "builder", "file": "seed_x.xyz"}, "seed_x.xyz")],
                         "ensemble": {"fx": -89.3}, "steps": [{"step": "s"}]}}
    g.seed_nodes = {(13, 3)}
    g.save_rows(13, 15)
    g2, _c = _fake_grid(monkeypatch, screen=0)
    g2.out, g2.seed_nodes = tmp_path, set()
    g2.seed = lambda *a, **k: (_ for _ in ()).throw(AssertionError("seed relaxed again"))
    called = []
    g2.descend = lambda m, k_max, m_stop=None, descent="column", row_hook=None: called.append((m, k_max, descent))
    g2.run_seeded([("seed_x.xyz", "up")], row_hook=lambda *_a: None)
    assert called == [(13, 15, "column")]
    x = g2.nodes[(13, 3)]["states"][0]
    assert (x.energy, x.fp, x.path, x.src["file"]) == (-89.3, "fx", "seed_x.xyz", "seed_x.xyz")
    assert np.allclose(x.relaxed, rel) and np.allclose(x.lat, rel + 0.5) and g2.seed_nodes == {(13, 3)}


def test_capacity_state_keeps_the_most_z_capable_product(monkeypatch):
    """With capacity on, a +MX_q step keeps, besides its lowest products, the one with the most Z-type host sites
    (protected path 'capacity', returned to the ladder); a later product with fewer sites does not replace it;
    remove_mx steps keep none; grid.json shows it only within delta of the best."""
    from qd_workflows import lattice_grid as lg
    g, _calls = _fake_grid(monkeypatch, screen=0)
    g.capacity = True
    sites = {0.0: 1, 1.0: 5, 2.0: 2}                    # by first coordinate (= energy of the fake relaxation)
    g.lat.z_sites = lambda syms, lat: ([None] * sites.get(round(float(lat[0, 0]), 3), 0), [])
    g.lat.remove_mx = lambda s, p: []
    parent = lg.Config_(["Cd"], np.zeros((1, 3)), np.zeros((1, 3)), 0.0, "p", {}, None)
    kept = g.step([parent], _moves(3), "+CdCl2", (1, 1), path="seed_x.xyz")
    by = sorted((x.energy, x.path) for x in g.nodes[(1, 1)]["states"])
    assert by == [(0.0, "seed_x.xyz"), (1.0, "capacity")]
    assert sorted(x.path for x in kept) == ["capacity", "seed_x.xyz"]
    assert next(x for x in kept if x.path == "capacity").src["z_sites"] == 5
    two = lambda _s, _p: [({"i": 2}, ["Cd", "Se", "Cl", "Cl"], np.array([[2.0, 0, 0], [3, 0, 0], [0, 3, 0], [0, 0, 3]]))]
    g.step([parent], two, "+CdCl2", (1, 1))                # 2 sites < 5: the capacity state stays
    assert [x.energy for x in g.nodes[(1, 1)]["states"] if x.path == "capacity"] == [1.0]
    g.step([parent], g.lat.remove_mx, "-CdCl2", (1, 2))    # no products and no capacity state for -MX_q
    assert (1, 2) not in g.nodes or not any(x.path == "capacity" for x in g.nodes[(1, 2)]["states"])


def test_solvation_block_uses_only_this_grid(tmp_path):
    """thermo_merged.json gets the surface of every state and the m = 1 complexes' COSMO from complexes/cosmo.json;
    without it the block says none ({}), so the parameters' complexes of another grid are not used."""
    pytest.importorskip("qd_workflows.reaction.solvation")
    from qd_workflows.grid_thermo import thermo_at
    (tmp_path / "grid" / "nodes").mkdir(parents=True)
    nu._write_xyz(tmp_path / "grid" / "nodes" / "m01_k01.xyz", ["Cd", "Se", "Cd", "Cl", "Cl"],
                  [[0, 0, 0], [2.6, 0, 0], [5.2, 0, 0], [7.6, 0, 0], [-2.4, 0, 0]], "x")
    g = lambda t: [-11.0 - 0.001 * t for t in (300.0, 500.0)]
    th = {"temperatures": [300.0, 500.0], "references": {"CdCl2": {"G_eV": [-7.0, -7.1]}, "CdSe": {"G_eV": [-2.3, -2.4]}},
          "nodes": [{"m": 1, "k": 1, "formula": "Cd2Se1Cl2", "ensemble_E_eV": [-11.0],
                     "states": [{"E_eV": -11.0, "G_eV": g(0), "file": "grid/nodes/m01_k01.xyz", "source": "grid"}]}]}
    (tmp_path / "thermo_merged.json").write_text(json.dumps(th))
    blk = nu.solvation_block(tmp_path, log=lambda _x: None)
    assert blk["complexes"] == {} and "none" in blk["complexes_source"]
    saved = json.loads((tmp_path / "thermo_merged.json").read_text())
    assert saved["nodes"][0]["states"][0]["S_vdW_A2"] > 0
    other = {"eps": 2.0, "monomers": {"CdSe": {"E_inf_eV": -0.3, "S_A2": 118.0}, "CdCl2": {"E_inf_eV": -0.2, "S_A2": 136.0}},
             "e_el_eV_A2": -1.06e-3, "complexes": {"1": {"E_inf_eV": -9.0, "S_A2": 1.0}}}       # another grid's value
    cluster_model = thermo_at(saved, 300.0, other)["G"][(1, 1)]
    (tmp_path / "complexes").mkdir()
    (tmp_path / "complexes" / "cosmo.json").write_text(json.dumps({"1": {"E_inf_eV": -0.28, "S_A2": 208.0}}))
    blk = nu.solvation_block(tmp_path, log=lambda _x: None)
    assert blk["complexes"] == {"1": {"E_inf_eV": -0.28, "S_A2": 208.0}}
    own_th = json.loads((tmp_path / "thermo_merged.json").read_text())
    own = thermo_at(own_th, 300.0, other)["G"][(1, 1)]
    bare = {k: v for k, v in own_th.items() if k != "solvation"}  # a thermo without its own block: the parameters'
    params_based = thermo_at(bare, 300.0, other)["G"][(1, 1)]
    from qd_workflows.grid_solvation import cosmo
    assert abs(own - cluster_model) > 1e-6
    assert abs((params_based - own) - (cosmo(other["complexes"]["1"], 2.0) - cosmo({"E_inf_eV": -0.28, "S_A2": 208.0}, 2.0))) < 1e-9


def test_climb_goes_past_the_top_from_the_capacity_states_below(monkeypatch):
    """The row's top node has no free site; a capacity state two nodes below has: the climb from below reaches past the
    top (a ladder from the top node alone would stop there)."""
    from qd_workflows import lattice_grid as lg
    g, _calls = _fake_grid(monkeypatch, screen=0)
    g.capacity = True
    g.lat.remove_mx = lambda s, p: []
    g.lat.z_sites = lambda syms, lat: ([None] * int(lat[0, 1]), [])
    comp = lambda k: ["Cd"] * (1 + k) + ["Se"] + ["Cl"] * (2 * k)

    def config(k, sites, e):
        p = np.zeros((len(comp(k)), 3))
        p[0] = [e, sites, 0]
        p[1:, 2] = np.arange(1, len(p))
        return comp(k), p

    def add(syms, lat):                                # one product per free site left, up to k = 5
        k = syms.count("Cl") // 2
        n = int(lat[0, 1])
        if n == 0 or k >= 5:
            return []
        s2, p2 = config(k + 1, n - 1, float(lat[0, 0]) - 1.0)
        return [({"host": 0}, s2, p2)]
    g.lat.add_mx = add
    for k, sites in ((1, 4), (2, 0), (3, 0)):          # k = 3 is the top; only k = 1 has capacity left
        s, p = config(k, sites, -float(k))
        g.nodes[(1, k)] = {"states": [lg.Config_(s, p, p, -float(k), f"s{k}", {}, None)], "ensemble": {}, "steps": []}
    assert g.ladder_up(1, 3) == 3                      # from the top node alone: no site
    assert g.climb(1, 1) == 5                          # from k = 1 through 2 and 3 to the cap of the fake rule
