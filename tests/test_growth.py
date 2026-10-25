"""Growth model: dot geometry, mass conservation, the no-nucleation limit (the precursor model alone) and the
builder recipe. A small synthetic grid thermochemistry stands in for a real one."""
import copy
import math

import numpy as np
import pytest

pytest.importorskip("scipy")
from qd_workflows import growth as gr  # noqa: E402
from qd_workflows.grid_thermo import converted  # noqa: E402

TEMPS = [300.0, 500.0, 700.0]


def _thermo():
    """Nodes [CdSe]m(CdCl2)k for m = 1..13, k = 0..3, a free energy that falls with size (a barrier near m = 2),
    and the CdCl2 and CdSe monomer references."""
    nodes = []
    for m in range(1, 14):
        for k in range(4):
            g = -5.2 * m - 7.0 * k - 0.6 * math.sqrt(m) * 0 + 1.2 * m ** (2 / 3) - 0.3 * k
            nodes.append({"m": m, "k": k, "formula": f"Cd{m + k}Se{m}Cl{2 * k}" if k else f"Cd{m}Se{m}",
                          "ensemble_E_eV": [g], "states": [{"E_eV": g, "G_eV": [g] * 3, "fingerprint": f"{m}-{k}", "file": "x"}]})
    return {"temperatures": TEMPS, "nodes": nodes,
            "references": {"CdCl2": {"G_eV": [-6.9] * 3}, "CdSe": {"G_eV": [-4.0] * 3}}}


def _params(**over):
    P = copy.deepcopy(gr.PARAMS)
    P["conditions"].update(t_end=60.0, n_out=12)
    P["material"]["g_bulk"] = -5.2                 # the large-size limit of the synthetic free energy
    P["solvation"]["eps"] = 1.0                    # gas phase: the synthetic grid has no surfaces
    for k, v in over.items():
        P[k].update(v)
    return P


def test_wulff_shape_has_the_nucleus_volume():
    sh = gr.Shape(gr.ZB_NORMALS)
    h = sh.wulff({"100": 0.035, "111A": 0.025, "111B": 0.055, "110": 0.045}, 13 * gr.PARAMS["material"]["Omega_m"])
    _a, vol, _s = sh(h)
    assert math.isclose(vol, 13 * gr.PARAMS["material"]["Omega_m"], rel_tol=1e-6)
    cube = sh({"100": 5.0, "111A": 50.0, "111B": 50.0, "110": 50.0})
    assert math.isclose(cube[0]["100"], 6 * 100.0, rel_tol=1e-6) and cube[0]["111A"] == 0.0


def test_coverage_law_obeys_gibbs_adsorption():
    kt = gr.KB_EV * 500.0
    p = {"gamma0": 20.0, "dG_bind": -0.5, "Gamma_max": 2.0, "w": 0.4}
    for c in (1e-4, 1e-2, 0.5):
        s0, s1 = gr.facet_state(p, c, kt), gr.facet_state(p, c * 1.001, kt)
        dg = (s1["gamma"] - s0["gamma"]) / math.log(1.001)
        assert math.isclose(dg, -s0["Gamma_max"] * s0["theta"] * kt, rel_tol=2e-3)       # d gamma / d ln c = -Gamma kT
    lang = gr.facet_state({**p, "w": 0.0}, 0.01, kt)                                   # w = 0: the Langmuir form
    K = math.exp(0.5 / kt)
    assert math.isclose(lang["theta"], K * 0.01 / (1 + K * 0.01), rel_tol=1e-9)
    assert math.isclose(lang["gamma"], 0.020 - 0.02 * kt * math.log1p(K * 0.01), rel_tol=1e-9)
    assert gr.facet_state(p, 0.01, kt)["theta"] < lang["theta"]                       # repulsion lowers coverage


def test_binding_free_energy_from_eps0():
    P = _params()
    g = gr.Growth(_thermo(), P)
    for f, p in g.F.items():
        assert math.isclose(p["dG_bind"], -p["eps0"] + g.dg_th)
    assert g.dg_th == 0.0                                      # the synthetic thermo has no monomer energy


def test_automatic_handover_picks_a_grid_size():
    P = _params()
    P["material"].update(m0="auto", m0_min=6)
    g = gr.Growth(_thermo(), P)
    ho = g.handover
    assert ho is not None and 6 <= ho["m0"] <= 13 and g.m0 == ho["m0"]
    best = min(abs(r["diff_eV"]) for r in ho["rows"])
    assert math.isclose(abs(next(r for r in ho["rows"] if r["m"] == ho["m0"])["diff_eV"]), best, abs_tol=0.005)
    assert all(r["m"] >= 6 for r in ho["rows"])
    P2 = _params()
    P2["material"]["m0"] = 13
    assert gr.Growth(_thermo(), P2).handover is None and gr.Growth(_thermo(), P2).m0 == 13


def test_implicit_solvation_raises_clusters_against_monomers():
    from qd_workflows import grid_solvation as gs
    from qd_workflows.grid_thermo import thermo_at
    th = _thermo()
    for n in th["nodes"]:
        for s in n["states"]:
            s["S_vdW_A2"] = 120.0 * (n["m"] + n["k"]) ** (2 / 3)                     # a surface growing as N^(2/3)
    sol = {**gs.SOLVATION, "eps": 2.0}
    gas, liq = thermo_at(th, 500.0), thermo_at(th, 500.0, sol)
    assert math.isclose(liq["dG_solv_mx"], gs.monomer("CdCl2", 2.0)) and liq["dG_solv_mx"] < 0
    assert math.isclose(liq["G"][(1, 0)] - gas["G"][(1, 0)], gs.monomer("CdSe", 2.0), abs_tol=1e-9)
    s13 = 120.0 * 16 ** (2 / 3)
    assert math.isclose(liq["G"][(13, 3)] - gas["G"][(13, 3)], gs.per_area(2.0) * s13, rel_tol=1e-9)
    # Omega(13, 3) relative to the monomers rises: the cluster keeps less solvation than its 13 + 3 units had
    rel = lambda t: t["G"][(13, 3)] - 13 * t["G"][(1, 0)] - 3 * t["g_mx"]          # noqa: E731
    assert rel(liq) > rel(gas)
    P = _params()
    P["solvation"]["eps"] = 2.0
    g = gr.Growth(th, P)
    for f, p in g.F.items():
        assert math.isclose(p["dG_bind"], -p["eps0"] + g.dg_th - liq["dG_solv_mx"])  # binding gives up solvation


def test_becker_doring_flux():
    kt, om_m = gr.KB_EV * 500.0, gr.PARAMS["material"]["Omega_m"]
    P = _params()
    ms = np.arange(1, 14)
    c1 = 0.01
    om = np.full(13, -kt * math.log(c1))                       # flat: every step as easy as the first
    om[5] += 0.5                                               # one rare cluster at m = 6
    J = gr.nucleation_flux(om, ms, 13, kt, om_m, P)
    d1 = P["kinetics"]["D"] * 1e16
    r = (3 * ms * om_m / (4 * math.pi)) ** (1 / 3)
    f6 = 4 * math.pi * d1 * (1 + 6 ** (-1 / 3)) * (r[0] + r[5]) * 6.02214076e-4 * c1
    assert math.isclose(J, f6 * c1 * math.exp(-0.5 / kt), rel_tol=0.05)     # the rarest cluster sets the flux
    assert J < gr.nucleation_flux(om - kt, ms, 13, kt, om_m, P)               # more monomer, more flux
    slow = _params(nucleation={"attach_barrier_eV": 0.1})
    assert math.isclose(gr.nucleation_flux(om, ms, 13, kt, om_m, slow), J * math.exp(-0.1 / kt), rel_tol=1e-9)
    const = _params(nucleation={"flux": "constant", "J0": 2.0})
    assert math.isclose(gr.nucleation_flux(om, ms, 13, kt, om_m, const, barrier=0.2), 2.0 * math.exp(-0.2 / kt))


def test_no_nucleation_is_the_precursor_model():
    P = _params(nucleation={"flux": "constant", "J0": 0.0})
    res = gr.Growth(_thermo(), P).run()
    c = P["conditions"]
    k = gr.k_conv_at(c["T"], c["k_conv_500K"], c["Ea_conv"])
    for row in res["series"]:
        assert row["N_dots"] == 0.0
        assert math.isclose(row["C_M"], converted(row["t"], k, c["c_P0"], c["c_MX0"]), rel_tol=2e-2, abs_tol=1e-6)


def test_mass_is_conserved_with_nucleation_and_growth():
    res = gr.Growth(_thermo(), _params(nucleation={"J0": 1e-4})).run()
    assert res["series"][-1]["N_dots"] > 0
    for c in gr.conservation(res):
        assert abs(c["MA"]) < 1e-6 and abs(c["MX"]) < 1e-6


def test_builder_recipe_uses_the_builder_keys():
    res = gr.Growth(_thermo(), _params(nucleation={"J0": 1e-4})).run()
    r = gr.builder_recipe(res, 60.0, {"charges": {"Cd": 2, "Se": -2, "Cl": -1}, "facets": [{"hkl": "100"}]})
    rec = r["recipe"]
    assert set(rec["facet_distances"]) == {"100", "111A", "111B", "110"} and "facets" not in rec
    assert rec["target"]["m"] == round(r["m"]) and rec["post_treatment"]["z_type_target"]["target_total"] == round(r["k"])


def test_size_law_limits():
    g = gr.Growth(_thermo(), _params(kinetics={"k_r": None}))
    g._arrays()
    A, r, th = np.array([400.0, 1600.0]), np.array([6.0, 12.0]), np.array([0.3, 0.3])
    c_eq, n = np.array([1e-4, 1e-4]), 6.02214076e-4
    dif = g.size_rates(A, r, c_eq, th, 1e-3)                                          # k_r None: diffusion limit
    assert np.allclose(dif, 4 * math.pi * g.D * r * n * (1e-3 - c_eq), rtol=1e-12)
    g.P["kinetics"]["k_r"], g.D = 50.0, 1e40                                          # fast diffusion: the reaction
    assert np.allclose(g.size_rates(A, r, c_eq, th, 1e-3), A * n * 50.0 * (1 - th) * (1e-3 - c_eq), rtol=1e-9)
    assert np.allclose(g.size_rates(A, r, c_eq, th, 1e-5), A * n * 50.0 * (1 - th) * (1e-5 - c_eq), rtol=1e-9)   # dissolving


def test_dot_chemical_potential_is_the_wulff_gibbs_thomson():
    g = gr.Growth(_thermo(), _params())
    g._arrays()
    _pl, fam = g.chem(1e-3, 0.1)
    gam = {f: fam[f]["gamma"] for f in g.F}
    h = g.shape.wulff(gam, 300 * g.Om)
    m, _k, area, _s = g.dot(h, fam)
    _A, r, dmu, _th = g.size_state([{"m": m, "area": area}], fam, {f: 0.0 for f in g.F})
    for f in g.F:
        if area[f] > 0:                                                               # every exposed facet: 2 gamma Om / h
            assert math.isclose(dmu[0], 2 * gam[f] * g.Om / h[f], rel_tol=1e-6)
    assert math.isclose(r[0], (3 * 300 * g.Om / (4 * math.pi)) ** (1 / 3), rel_tol=1e-6)


def test_ripening_without_supply_conserves_and_coarsens():
    P = _params(conditions={"c_P0": 0.0, "t_end": 600.0}, nucleation={"flux": "constant", "J0": 0.0})
    for f in P["families"].values():
        f["eps0"] = 0.3                     # weak binding: every gamma > 0 (a full shell makes gamma(100) < 0 and the
                                            # sizes converge instead, the inverse ripening of a negative surface energy)
    res = gr.Growth(_thermo(), P).run(seed=[{"N": 1e-5, "m": 150}, {"N": 1e-5, "m": 400}])
    end = res["snapshots"][-1]["cohorts"]
    assert len(end) == 1 or end[0]["m"] < 150 < 400 < end[1]["m"]                     # the small dot feeds the large
    assert max(co["m"] for co in end) > 400
    for c in gr.conservation(res):
        assert abs(c["MA"]) < 1e-9 and abs(c["MX"]) < 1e-9


def test_grid_mu_is_the_two_step_free_energy():
    g = gr.Growth(_thermo(), _params())
    g._arrays()
    kt, mu_mx = g.kt, -7.2
    gt = {}
    for m in range(1, 14):
        e = [g.th["G"][(m, k)] - k * mu_mx for k in range(4)]
        gt[m] = min(e) - kt * math.log(sum(math.exp(-(x - min(e)) / kt) for x in e))
    ms, mu = g.grid_mu(mu_mx)
    assert ms.tolist() == list(range(3, 14))
    assert np.allclose(mu, [0.5 * (gt[m] - gt[m - 2]) for m in range(3, 14)], atol=1e-9)
    assert g.grid_mu_on and g.below.sum() == 13                                       # the chain runs to the grid's top
    cont = gr.Growth(_thermo(), _params(kinetics={"mu_small": "continuum"}))
    cont._arrays()
    assert not cont.grid_mu_on and cont.below.sum() == cont.m0


def test_grid_dots_dissolve_or_grow_by_the_grid_mu():
    """A dot within the grid sizes sits in equilibrium with the monomer level that matches the grid's step there."""
    P = _params(conditions={"c_P0": 0.0, "t_end": 1e-3, "n_out": 2, "t_start": 1e-4, "c_M0": 1e-3},
                nucleation={"flux": "constant", "J0": 0.0})
    g = gr.Growth(_thermo(), P)
    res = g.run(seed=[{"N": 1e-9, "m": 10.0}])
    pl, _fam = g.chem(1e-3, P["conditions"]["c_MX0"])
    ms, mu = g.grid_mu(pl["mu_MX"])
    m_end = res["snapshots"][-1]["cohorts"][0]["m"] if res["snapshots"][-1]["cohorts"] else 0.0
    grows = pl["mu_MA"] > float(np.interp(10.0, ms, mu))
    assert (m_end > 10.0) == grows
    for c in gr.conservation(res):
        assert abs(c["MA"]) < 1e-9 and abs(c["MX"]) < 1e-9
