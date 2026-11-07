"""Network growth model (network_growth) on a synthetic network: detailed balance, conservation, island barriers."""
import json
import math

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
             {"a": 1, "b": 2, "dm": 1, "dk": 1, "kind": "layer", "facet": "111B", "sheet": 3, "area_A2": 60.0},
             {"a": 1, "b": 3, "dm": 1, "dk": 0, "kind": "layer", "facet": "100", "sheet": 2, "area_A2": 25.0},
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


def test_pruned_run(tmp_path):
    cond = {"T": 500.0, "eps": 1.0, "c_P0": 0.01, "c_MX0": 0.02, "k_conv_500K": 10.0, "t_start": 1e-3,
            "t_end": 5.0, "n_out": 6}
    net = _network(tmp_path)
    exact = ng.simulate(net, {**cond, "prune_eV": 0.0})       # on the whole network (sets its matrix pattern)
    wide = ng.simulate(net, {**cond, "prune_eV": 15.0})       # keeps every state: the same run
    assert wide["states"] == net.n
    assert np.allclose(np.array(wide["N"]), np.array(exact["N"]), rtol=1e-6, atol=1e-14)
    narrow = ng.simulate(net, {**cond, "prune_eV": 0.01})     # grid nodes, the seed and each size's best only
    assert narrow["states"] < net.n
    x = np.array(narrow["alpha"]) * cond["c_P0"]
    held = np.array([net.m @ n for n in narrow["N"]])
    assert np.allclose(np.array(narrow["c_MA"]) + held, x, rtol=1e-5, atol=1e-12)


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


def test_conversion_takes_the_salt_left_in_solution(tmp_path):
    """dx/dt = k_conv (c_P0 - x) C_X with C_X the salt not held by the states: slower than the open-loop closed form
    once the states hold salt, and the salt in solution never runs negative."""
    from qd_workflows.grid_thermo import converted
    cond = {"T": 500.0, "eps": 1.0, "c_P0": 0.01, "c_MX0": 0.0105, "k_conv_500K": 10.0, "t_start": 1e-3,
            "t_end": 50.0, "n_out": 12}
    net = _network(tmp_path)
    r = ng.simulate(net, cond)
    t, x, c_x = np.array(r["t"]), np.array(r["alpha"]) * cond["c_P0"], np.array(r["c_X_left"])
    held_k = np.array([net.k @ n for n in r["N"]])
    assert r["stopped"] is None and np.all(c_x > 0)
    assert np.allclose(c_x + held_k, cond["c_MX0"] - x, rtol=1e-6, atol=1e-12)
    assert x[-1] < 0.9 * converted(t[-1], 10.0, cond["c_P0"], cond["c_MX0"])
    q = 10.0 * (cond["c_P0"] - x) * c_x                        # the rate on the outputs, by trapezoids
    dx, qt = np.diff(x), 0.5 * (q[1:] + q[:-1]) * np.diff(t)
    assert np.median(np.abs(dx / qt - 1)) < 0.02


def test_collision_areas_come_from_the_structures(tmp_path):
    """Collision-limited (a tiny viscosity): nu / area is one constant per carrier, with the export's area_A2 for a
    layer step and the dot's S for an attachment; a network exported without area_A2 is refused."""
    net = _network(tmp_path)
    nu, src = net.prefactors(500.0, 1e-6, 1.0)                  # k_d ~ 1e4 k_r
    lay = np.where(net.kind == ng.KIND["layer"])[0]
    att = np.where(net.kind == ng.KIND["attach"])[0][0]
    assert sorted(net.area[lay]) == [25.0, 60.0]
    per = nu[lay] / net.area[lay]
    assert abs(per[0] / per[1] - 1) < 1e-3 and abs(nu[att] / net.S[net.a[att]] / per[0] - 1) < 1e-3
    assert np.all(src > 0)
    meta = json.loads((tmp_path / "network.json").read_text())
    del meta["steps"][1]["area_A2"]
    (tmp_path / "network.json").write_text(json.dumps(meta))
    try:
        ng.Network(tmp_path, _thermo(), {"eps": 1.0, "vdw_scale": 1.4, "e_el_eV_A2": -1.06e-3})
    except ValueError as e:
        assert "area_A2" in str(e)
    else:
        raise AssertionError("a layer step without area_A2 was accepted")


def test_species_names_from_the_seed():
    """X from the seed's formula (Cl, or the carboxylate with its short name); the template's block overrides."""
    cl = {"seed": 0, "states": [{"formula": "Cd16Se13Cl6", "m": 13, "k": 3}]}
    ac = {"seed": 0, "states": [{"formula": "Cd19Se16C12O12H18", "m": 16, "k": 3}]}
    assert ng.species(cl) == {"M": "Cd", "A": "Se", "X": "Cl", "MA": "CdSe", "MX": "CdCl2"}
    assert ng.species(ac)["X"] == "OAc" and ng.species(ac)["MX"] == "Cd(OAc)2"
    assert ng.species(ac, {"X": "acetate", "MX": "Cd(OAc)2"})["X"] == "acetate"


def test_solvent_sets_eps_and_viscosity():
    from qd_workflows import solvents
    hot, cool = ng.conditions({"T": 550.0}), ng.conditions({"T": 420.0})
    assert hot["solvent"] == "1-octadecene" and hot["eps"] == solvents.properties("1-octadecene", 550.0)[0]
    assert hot["eta_mPas"] < cool["eta_mPas"]                   # thinner when hotter
    assert ng.conditions({"eps": 1.0})["eps"] == 1.0 and ng.conditions({"eps": 1.0})["solvent"] == "custom"
    for row in solvents.table():                                  # every entry: through its points, thinner when hotter
        assert 1.0 < row["eps"] < 50 and row["basis"] in ("measured", "correlation", "estimate")
        for T, eta in row["eta"]:
            assert abs(solvents.viscosity(row["name"], T) / eta - 1) < 1e-9
        etas = [solvents.viscosity(row["name"], T) for T in (350.0, 450.0, 550.0, 650.0)]
        assert all(a > b > 0 for a, b in zip(etas, etas[1:]))


def test_syringe_injection_conserves_amounts(tmp_path):
    """The syringe adds precursor and volume: CdSe and salt balance in amounts per litre of the starting flask
    (concentration / phi), the flask ends at V0 + rate x time, and the pool without clusters matches the closed form
    once the syringe has nothing to add (c_syringe_M -> its own amount in the flask)."""
    cond = {"T": 500.0, "eps": 1.0, "c_P0": 0.0, "c_MX0": 0.05, "k_conv_500K": 10.0, "t_start": 1e-3,
            "t_end": 50.0, "n_out": 10, "inject": True, "c_syringe_M": 0.2, "rate_mL_h": 360.0, "V0_mL": 10.0,
            "inj_duration_s": 20.0}
    net = _network(tmp_path)
    r = ng.simulate(net, cond)
    assert r["stopped"] is None and 20.0 in [round(t, 9) for t in r["t"]]
    V = np.array(r["V_mL"]); phi = 10.0 / V
    assert abs(V[-1] - 12.0) < 1e-9 and np.all(np.diff(V) >= 0)
    total = 0.2 * 2.0 / 10.0                                     # mol per litre of the starting flask
    feed = ng.Feed(ng.conditions(cond))
    assert abs(feed.total - total) < 1e-15
    x = np.array(r["alpha"]) * np.array([feed.precursor(t) for t in r["t"]])   # alpha: of the precursor delivered
    held_m = np.array([net.m @ n for n in r["N"]]); held_k = np.array([net.k @ n for n in r["N"]])
    assert np.allclose((np.array(r["c_MA"]) + held_m) / phi, x, rtol=1e-6, atol=1e-12)
    assert np.allclose((np.array(r["c_X_left"]) + held_k) / phi, cond["c_MX0"] - x, rtol=1e-6, atol=1e-12)
    cm, cx, a = ng.precursor_only(cond, 1e4)                    # long after: everything delivered and converted
    assert abs(a - 1.0) < 1e-6 and abs(cm - total * 10.0 / 12.0) < 1e-9
    flask = {**cond, "inject": False, "c_P0": 0.03}               # no syringe: the closed form
    from qd_workflows.grid_thermo import converted
    assert abs(ng.precursor_only(flask, 5.0)[0] - converted(5.0, 10.0, 0.03, 0.05)) < 1e-15
