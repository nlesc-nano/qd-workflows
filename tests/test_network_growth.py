"""Network growth model (network_growth) on a synthetic network: detailed balance, conservation, island barriers."""
import json

import numpy as np

from qd_workflows import network_growth as ng
from qd_workflows.grid_thermo import pool

T = [100.0, 300.0, 500.0]


def _node(m, k, energies, g_shift=0.1):
    states = [{"fingerprint": f"{m}-{k}-{i}", "E_eV": e, "G_eV": [e + g_shift * t / 300 for t in T],
               "file": f"nodes/m{m:02d}_k{k:02d}.xyz", "S_vdW_A2": 150.0 * m + 100.0 * k}
              for i, e in enumerate(energies[:1])]
    return {"m": m, "k": k, "formula": f"Cd{m + k}Se{m}" + (f"Cl{2 * k}" if k else ""),
            "ensemble_E_eV": sorted(energies), "states": states}


def _thermo():
    nodes = [_node(1, 0, [-2.2]), _node(1, 1, [-9.5]), _node(1, 2, [-16.6, -16.5]), _node(1, 3, [-23.4]),
             _node(2, 0, [-6.0]), _node(2, 1, [-13.8, -13.75]), _node(2, 2, [-21.0])]
    return {"temperatures": T, "nodes": nodes,
            "references": {"CdCl2": {"G_eV": [-6.76, -6.86, -6.97]}, "CdSe": {"G_eV": [-2.3, -2.4, -2.5]}}}


def _network(tmp_path, island_dg=-1.9):
    """Rows: (2,1) (the grid node), (3,1), (4,2), (4,1); a {111}B and a {100} layer step from (3,1), one shell
    step; islands on (3,1) {111}B with dm = 1 (uphill) and dm = 2 (downhill)."""
    states = [{"id": i, "formula": f, "m": m, "k": k, "fw": i, "n_imag": 0}
              for i, (f, m, k) in enumerate([("Cd3Se2Cl2", 2, 1), ("Cd4Se3Cl2", 3, 1), ("Cd6Se4Cl4", 4, 2),
                                             ("Cd5Se4Cl2", 4, 1)])]
    steps = [{"a": 0, "b": 1, "dm": 1, "dk": 0, "kind": "attach", "facet": None, "sheet": 0},
             {"a": 1, "b": 2, "dm": 1, "dk": 1, "kind": "layer", "facet": "111B", "sheet": 3},
             {"a": 1, "b": 3, "dm": 1, "dk": 0, "kind": "layer", "facet": "100", "sheet": 2},
             {"a": 3, "b": 2, "dm": 0, "dk": 1, "kind": "shell", "facet": None, "sheet": 0}]
    islands = [{"parent": 1, "family": "111B", "dm": 1, "dk": 0, "src": "island", "row": 4},
               {"parent": 1, "family": "111B", "dm": 2, "dk": 0, "src": "island", "row": 5}]
    g500 = np.array([-13.6, -16.6, -27.0, -19.6, -16.6 + island_dg, -16.6 - 5.6])
    G = np.array([g + 0.05 * (np.array(T) - 500) / 300 for g in g500])
    S = np.array([400.0, 550.0, 800.0, 700.0, 600.0, 650.0])
    (tmp_path / "network.json").write_text(json.dumps({"temperatures": T, "seed": 1, "seed_id": 1, "states": states,
                                                        "steps": steps, "islands": islands}))
    np.savez(tmp_path / "arrays.npz", G_gas=G, S=S, sigma=np.ones(len(S), int))
    return ng.Network(tmp_path, _thermo(), {"eps": 1.0, "vdw_scale": 1.4, "e_el_eV_A2": -1.06e-3})


def test_grid_node_takes_the_small_network_state(tmp_path):
    net = _network(tmp_path)
    assert net.n_grid == 3 and net.n == 6                       # (2,0) (2,1) (2,2), then rows 1-3
    assert net.state_of_row[0] == net.grid_keys.index((2, 1))
    kinds = {int(k) for k in net.kind}
    assert ng.KIND["layer"] in kinds and ng.KIND["grid"] in kinds and len(net.src) > 0


def test_detailed_balance_at_a_fixed_pool(tmp_path):
    net = _network(tmp_path)
    _g, th, _pa = net.G(500.0, 1.0)
    for c_m in (1e-4, 3e-3):
        pl = ng._pool(th, c_m, 0.02)
        assert ng.stationary_residual(net, 500.0, 1.0, pl) < 1e-10


def test_pool_species_are_the_sources(tmp_path):
    """Omega of a pool complex is -kT ln c_n (the source state's population)."""
    net = _network(tmp_path)
    _g, th, _pa = net.G(500.0, 1.0)
    pl = pool(th, 2e-3, 0.02)
    kt = ng.KB_EV * 500
    for n, x in pl["speciation"].items():
        om = th["G"][(1, n)] - pl["mu_MA"] - n * pl["mu_MX"]
        assert abs(om + kt * np.log(2e-3 * x)) < 1e-9


def test_barrier_rule():
    assert ng.Network.barrier({1: 0.7, 2: -0.4, 3: 1.5}) == 0.7          # before the first downhill point only
    assert ng.Network.barrier({1: -0.2, 2: 0.5}) == 0.0
    assert ng.Network.barrier({1: 0.3, 2: 0.9}) == 0.9                   # never downhill: the highest point


def test_island_barrier_on_its_layer_step_only(tmp_path):
    net = _network(tmp_path)
    _g, th, _pa = net.G(500.0, 1.0)
    pl = ng._pool(th, 2e-3, 0.02)
    B, table = net.barriers(500.0, 1.0, pl["mu_MA"], pl["mu_MX"])
    i111 = np.where((net.kind == ng.KIND["layer"]) & (net.fam == ng.FAMILIES.index("111B")))[0][0]
    i100 = np.where((net.kind == ng.KIND["layer"]) & (net.fam == ng.FAMILIES.index("100")))[0][0]
    expect = -1.9 - pl["mu_MA"]
    assert abs(B[i111] - max(expect, 0.0)) < 1e-9 and table["own"]["1|111B"] == B[i111]
    assert B[i100] == table["all"]                              # {100} has no islands: the median over all
    assert np.all(B[net.kind != ng.KIND["layer"]] == 0)


def test_conservation_and_island_slowdown(tmp_path):
    cond = {"T": 500.0, "eps": 1.0, "c_P0": 0.01, "c_MX0": 0.02, "k_conv_500K": 10.0, "t_start": 1e-3,
            "t_end": 5.0, "n_out": 6}
    runs = {}
    for dg in (-1.9, -0.5):                     # island Delta G: Omega*(1) = dg - mu_MA, ~0.7 and ~2.1 eV
        net = _network(tmp_path, island_dg=dg)
        r = ng.simulate(net, cond)
        x = np.array(r["alpha"]) * cond["c_P0"]
        held_m = np.array([net.m @ n for n in r["N"]])
        held_k = np.array([net.k @ n for n in r["N"]])
        assert np.allclose(np.array(r["c_MA"]) + held_m, x, rtol=1e-6, atol=1e-12)
        assert np.allclose(np.array(r["c_X_left"]) + held_k, cond["c_MX0"] - x, rtol=1e-6, atol=1e-12)
        runs[dg] = r
    lo, hi = runs[-1.9]["cdse_flux"]["111B"], runs[-0.5]["cdse_flux"]["111B"]
    assert max(lo[1:]) > max(hi[1:])                            # the higher barrier grows {111}B slower


def test_run_stops_with_its_outputs(tmp_path):
    cond = {"T": 500.0, "eps": 1.0, "c_P0": 0.01, "c_MX0": 0.02, "k_conv_500K": 10.0, "t_start": 1e-3,
            "t_end": 5.0, "n_out": 6, "dense_out": False}
    net = _network(tmp_path)
    full = ng.simulate(net, cond)
    assert full["stopped"] is None and len(full["t"]) == 6
    calls = []
    r = ng.simulate(net, cond, stop=lambda: len(calls.append(0) or calls) > full["steps"] // 2)   # from the page
    assert r["stopped"].startswith("stopped from the page") and 2 <= len(r["t"]) <= 6
    n = len(r["t"]) - 1                                       # the outputs before the stop, then its last state
    assert np.allclose(r["c_MA"][:n], full["c_MA"][:n], rtol=1e-9) and full["t"][n - 1] < r["t"][-1] < full["t"][n]
    r = ng.simulate(net, {**cond, "max_wall_s": 1.0}, stop=lambda: __import__("time").sleep(0.3) or False)
    assert r["stopped"].startswith("wall-time limit") and len(r["t"]) < 6


def test_dense_outputs_keep_the_log_grid(tmp_path):
    cond = {"T": 500.0, "eps": 1.0, "c_P0": 0.01, "c_MX0": 0.02, "k_conv_500K": 10.0, "t_start": 1e-3,
            "t_end": 5.0, "n_out": 6}
    net = _network(tmp_path)
    log_only = ng.simulate(net, {**cond, "dense_out": False})
    dense = ng.simulate(net, cond)
    t = np.array(dense["t"])
    assert len(t) > len(log_only["t"]) and np.all(np.diff(t) > 0)
    for t_log, c in zip(log_only["t"], log_only["c_MA"]):        # the log-spaced outputs are still there, unchanged
        j = int(np.argmin(np.abs(t - t_log)))
        assert abs(t[j] - t_log) < 1e-9 * t_log and abs(dense["c_MA"][j] - c) < 1e-9 * c
