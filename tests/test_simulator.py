"""Simulator backend: a run returns the page's series and cohort summaries; the builder button refuses dots above
the atom limit and, with QD_builder installed, builds the most numerous dot."""
import json

import pytest

pytest.importorskip("scipy")
from qd_workflows import simulator as sm  # noqa: E402
from qd_workflows.grid_thermo import converted  # noqa: E402
from test_growth import _params, _thermo  # noqa: E402

PRESET = {"material": "CdSe", "cif": "", "native_elements": ["Cd", "Se"], "positive_q_mode": "add",
          "preset": {"charges": {"Cd": 2, "Se": -2, "Cl": -1}}}


def _sim(max_atoms=sm.MAX_ATOMS):
    return sm.Simulator(_thermo(), _params(nucleation={"J0": 1e-4}), PRESET, max_atoms)


def test_run_returns_series_and_cohorts():
    sim = _sim()
    out = sim.run({"T": 500.0, "t_end": 30.0})
    assert out["run_id"] in sim.runs and out["series"][-1]["N_dots"] > 0
    json.dumps(out, allow_nan=False)                         # valid JSON for the page: an infinite barrier is null
    co = out["snapshots"][-1]["cohorts"][0]
    assert {"i", "N", "m", "k", "d_nm", "n_atoms", "facets"} <= set(co)
    assert abs(sum(co["facets"].values()) - 1.0) < 1e-9
    with pytest.raises(ValueError):
        sim.run({"c_P0": 0.3, "c_MX0": 0.2})


def test_omega_at_the_pool_of_a_time():
    sim = _sim()
    out = sim.run({"t_end": 30.0})
    t_peak = max(out["series"], key=lambda r: r["C_M"])["t"]                         # monomer builds up, then the dots
    early, peak = sim.omega(out["run_id"], 1e-3), sim.omega(out["run_id"], t_peak)  # take it
    assert early["mu_MA"] < peak["mu_MA"] and early["barrier"] >= peak["barrier"]
    late = sim.omega(out["run_id"], 30.0)
    assert len(late["Omega"]) == len(late["ms"]) and len(late["Omega"][0]) == len(late["ks"])
    with pytest.raises(FileNotFoundError):
        sim.structure(late["m0"], 0)                                                 # the synthetic grid has no files


def test_live_omega_is_the_pre_nucleation_pool():
    sim = _sim()
    o = sim.omega_live({"T": 500.0, "c_P0": 0.1, "c_MX0": 0.2, "k_conv_500K": 1.0, "Ea_conv": 0.0, "t": 10.0})
    assert abs(o["C_M"] - converted(10.0, 1.0, 0.1, 0.2)) < 1e-12 and abs(o["C_X"] - (0.2 - o["C_M"])) < 1e-12
    assert o["m_star"] <= o["m0"] and o["barrier"] >= 0 and len(o["curve"]) > 10
    h = o["heat"]                                                                    # rows m = 1..m0, columns t
    assert len(h["dOmega"]) == len(h["m"]) == o["m0"] and len(h["dOmega"][0]) == len(h["t"]) == len(o["curve"])
    assert all(abs(x) < 1e-12 for x in h["dOmega"][0])                               # Delta Omega(1) = 0
    for j, c in enumerate(o["curve"]):
        assert abs(max(row[j] for row in h["dOmega"]) - c["barrier"]) < 1e-9          # the barrier is the column max
    run = sim.run({"T": 500.0, "c_P0": 0.1, "c_MX0": 0.2, "k_conv_500K": 1.0, "Ea_conv": 0.0, "t_end": 30.0})
    early = sim.omega(run["run_id"], 1e-3)                                           # no nuclei yet: the same pool
    live = sim.omega_live({"T": 500.0, "c_P0": 0.1, "c_MX0": 0.2, "k_conv_500K": 1.0, "Ea_conv": 0.0, "t": early["t"]})
    assert abs(early["mu_MA"] - live["mu_MA"]) < 0.02


def test_facet_overrides_are_checked_and_marked():
    sim = _sim()
    d = sim.defaults()["surface"]
    assert set(d["editable"]) == set(sm.EDITABLE)
    plain = sim.run({"t_end": 30.0})
    assert plain["modified"] is False
    f = next(iter(d["families"]))
    edited = sim.run({"t_end": 30.0, "families": {f: {"gamma0": 60.0}}})
    assert edited["modified"] is True and edited["surface_used"][f]["gamma0"] == 60.0     # used as typed (no rescale)
    same = sim.run({"t_end": 30.0, "families": {f: {"gamma0": d["families"][f]["gamma0"]}}})
    assert same["modified"] is False                                                     # unchanged value: not edited
    for bad in ({f: {"gamma0": 1e4}}, {f: {"tau": 1.0}}, {"999": {"gamma0": 30.0}}):
        with pytest.raises(ValueError):
            sim.run({"t_end": 30.0, "families": bad})


def test_build_refuses_large_dots():
    sim = _sim(max_atoms=1)
    out = sim.run({"t_end": 30.0})
    with pytest.raises(ValueError, match="atom limit"):
        sim.build(out["run_id"], 30.0, None)
    with pytest.raises(KeyError):
        sim.build("nope", 30.0, None)


def test_injection_temperature_program_and_activities_reach_the_model():
    sim = _sim()
    cond = {"t_end": 600.0, "inject": True, "c_syringe_M": 0.3, "rate_mL_h": 120.0, "V0_mL": 10.0,
            "inj_duration_s": 60.0, "two_T": True, "T_growth": 480.0, "t_switch_s": 30.0, "ramp_s": 0.0,
            "gamma_MA": 2.0, "gamma_MX": 0.5, "beta_scale": 0.8}
    P = sim.params_for(cond)
    c = P["conditions"]
    assert c["c_P0"] == 0.0 and abs(c["injection"]["rate_L_h"] - 0.12) < 1e-12 and c["injection"]["V0_L"] == 0.01
    assert c["T_growth"] == 480.0 and c["gamma_MA"] == 2.0 and c["gamma_MX"] == 0.5 and P["kinetics"]["beta_scale"] == 0.8
    out = sim.run(cond)
    last = out["series"][-1]
    assert abs(last["V_L"] - 0.012) < 1e-12 and abs(last["T"] - 480.0) < 1.0
    assert abs(out["conservation"]["MA"]) < 1e-9 and abs(out["conservation"]["MX"]) < 1e-9
    om = out["omega_run"]                                          # the grid's Omega along the run, for the live panel
    assert len(om["t"]) == len(om["barrier"]) == len(om["m_star"]) > 0
    assert len(om["dOmega"]) == len(om["m"]) and len(om["dOmega"][0]) == len(om["t"])
    json.dumps(om, allow_nan=False)
    live = sim.omega_live({**cond, "t": 45.0})
    assert live["injection"] and abs(live["T_now"] - 480.0) < 1e-9 and 0 < live["alpha"] <= 1.0
    with pytest.raises(ValueError):                                                    # more Se than CdCl2
        sim.params_for({**cond, "c_syringe_M": 1.0, "rate_mL_h": 1000.0, "inj_duration_s": 3600.0})


def test_classical_twin_from_the_atomistic_model():
    sim = _sim()
    m = sim.from_atomistic({"T": 500.0})
    mat = m["material"]
    g = sm.gr.Growth(sim.thermo, sim.params_for({"T": 500.0}))
    g._arrays()
    cm, cx = g.ref_pool()
    assert abs(g.pool(mat["C_inf_M"], cx)["mu_MA"] - g.gb) < 1e-6                 # C_inf: mu_MA = g_bulk
    assert 0.01 <= mat["gamma_J_m2"] <= 2.0 and mat["Vm_A3"] == g.Om and mat["D_cm2_s"] > 0 and mat["k_r_A_s"] > 0
    if m["raw"]["gamma_J_m2"] < 0.01:                     # the synthetic grid's full shell makes gamma < 0: clipped, noted
        assert any("gamma_J_m2" in n for n in m["notes"])
    out = sim.classical_run({"T": 500.0, "t_end": 60.0}, mat)
    assert abs(out["conservation"]["MA"]) < 1e-9 and out["series"][-1]["N_dots"] > 0
    json.dumps(out, allow_nan=False)
    inj = sim.classical_run({"t_end": 300.0, "inject": True, "c_syringe_M": 0.3, "rate_mL_h": 120.0, "V0_mL": 10.0,
                             "inj_duration_s": 60.0}, mat)
    assert abs(inj["series"][-1]["V_L"] - 0.012) < 1e-12 and abs(inj["conservation"]["MA"]) < 1e-9


def test_pre_pool_closed_form_and_injection_balance():
    sim = _sim()
    P = sim.params_for({"T": 500.0})
    _T, cm, _cx = sim.pre_pool(P, [10.0])[0]
    c = P["conditions"]
    assert abs(cm - converted(10.0, c["k_conv_500K"], c["c_P0"], c["c_MX0"])) < 1e-12
    Pi = sim.params_for({"inject": True, "c_syringe_M": 0.1, "rate_mL_h": 36.0, "V0_mL": 10.0, "inj_duration_s": 100.0})
    (_T, cm1, cx1), = sim.pre_pool(Pi, [1e4])                    # long after: all converted, everything diluted
    v_end = 0.01 + 0.036 / 3600 * 100
    assert abs(cm1 * v_end - 0.1 * 0.036 / 3600 * 100) < 1e-9
    assert abs((cx1 + cm1) * v_end - Pi["conditions"]["c_MX0"] * 0.01) < 1e-9
