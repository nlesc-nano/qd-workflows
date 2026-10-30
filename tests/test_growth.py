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
            assert math.isclose(dmu[0] - g.junction(m), 2 * gam[f] * g.Om / h[f], rel_tol=1e-6)   # + grid junction
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


def test_facet_law_reduces_to_the_size_law_on_a_wulff_dot():
    facets = {f: {"k_f": 50.0} for f in ("100", "111A", "111B", "110")}
    P = _params(kinetics={"growth_law": "facet", "facets": facets, "k_r": 50.0})
    g = gr.Growth(_thermo(), P)
    g._arrays()
    pl, fam = g.chem(1e-3, 0.1)
    for f in g.F:                                                     # one coverage on every facet
        fam[f]["theta"] = 0.3
    th = {f: 0.3 for f in g.F}
    gam = {f: fam[f]["gamma"] for f in g.F}
    h = g.shape.wulff(gam, 400 * g.Om)
    m, k, area, _s = g.dot(h, fam)
    co = {"m": m, "k": k, "h": h, "area": area}
    A, r, dmu, thb = g.size_state([co], fam, th)
    C, C_M = 2e-3, 1e-3
    c_eq = C_M * np.exp((g.gb + dmu - pl["mu_MA"]) / g.kt)
    size = g.size_rates(A, r, c_eq, thb, C)[0]
    facet, v = g._facet_rate(co, C, C_M, pl["mu_MA"], fam, th)
    assert math.isclose(facet, size, rel_tol=1e-6)
    exposed = [f for f in g.F if area[f] > 0]
    assert np.allclose([v[f] for f in exposed], v[exposed[0]], rtol=1e-6)   # Wulff: every facet at one velocity


def test_island_limited_facet_grows_large():
    cond = {"c_P0": 0.05, "t_end": 30.0, "n_out": 4, "t_start": 1.0, "c_M0": 2e-2}

    def weak(P):                                                       # weak binding: every gamma > 0, several facets
        for f in P["families"].values():
            f["eps0"] = 0.3
        return P

    probe = gr.Growth(_thermo(), weak(_params(conditions=cond)))
    probe._arrays()
    fam0 = probe.chem(1e-4, 0.2)[1]
    area0 = probe.dot(probe.shape.wulff({f: fam0[f]["gamma"] for f in probe.F}, 2000 * probe.Om), fam0)[2]
    slow = min((f for f in area0 if area0[f] > 0), key=lambda f: area0[f])      # the smallest exposed facet
    facets = {f: ({"island_beta_eV": 0.25} if f == slow else {"k_f": 1e6}) for f in probe.F}
    P = weak(_params(conditions=cond, nucleation={"flux": "constant", "J0": 0.0},
                     kinetics={"growth_law": "facet", "facets": facets, "k_r": 1e6}))
    g = gr.Growth(_thermo(), P)
    res = g.run(seed=[{"N": 1e-6, "m": 2000}])
    first, last = res["snapshots"][0]["cohorts"][0], res["snapshots"][-1]["cohorts"][0]
    share = lambda co: co["area"][slow] / sum(co["area"].values())   # noqa: E731
    assert last["m"] > first["m"]
    assert share(last) > share(first)                                 # the slow facet takes over the surface
    for c in gr.conservation(res):
        assert abs(c["MA"]) < 1e-9 and abs(c["MX"]) < 1e-9


def test_shape_rule_keeps_the_size_law_flux_and_follows_the_slow_facets():
    cond = {"c_P0": 0.05, "t_end": 30.0, "n_out": 4, "t_start": 1.0, "c_M0": 2e-2}

    def weak(P):                                                       # weak binding: every gamma > 0, several facets
        for f in P["families"].values():
            f["eps0"] = 0.3
        return P

    probe = gr.Growth(_thermo(), weak(_params(conditions=cond)))
    probe._arrays()
    fam0 = probe.chem(1e-4, 0.2)[1]
    area0 = probe.dot(probe.shape.wulff({f: fam0[f]["gamma"] for f in probe.F}, 2000 * probe.Om), fam0)[2]
    main = max(area0, key=lambda f: area0[f])                          # the Wulff dot's dominant family
    runs = {}
    for law, ratio in (("size", None), ("slow", 0.01), ("fast", 100.0)):
        facets = None if ratio is None else {f: {"k_f": 1.0 if f == main else 1.0 / ratio} for f in probe.F}
        P = weak(_params(conditions=cond, nucleation={"flux": "constant", "J0": 0.0},
                         kinetics={"growth_law": "size" if ratio is None else "shape", "facets": facets}))
        g = gr.Growth(_thermo(), P)
        runs[law] = g.run(seed=[{"N": 1e-6, "m": 2000}])
        for c in gr.conservation(runs[law]):
            assert abs(c["MA"]) < 1e-9 and abs(c["MX"]) < 1e-9
    share = lambda co: co["area"][main] / sum(co["area"].values())   # noqa: E731
    last = {k: r["snapshots"][-1]["cohorts"][0] for k, r in runs.items()}
    for k in ("slow", "fast"):                                         # the size law sets how much the dot grows
        assert math.isclose(last[k]["m"], last["size"]["m"], rel_tol=0.2)
    assert share(last["slow"]) > share(last["size"]) - 1e-6            # slow: it keeps the surface
    assert share(last["fast"]) < 0.5 * share(last["size"])             # fast: it grows out, the others take over


def test_activity_coefficients_shift_the_potentials_and_keep_the_balance():
    g = gr.Growth(_thermo(), _params())
    g._arrays()
    kt = g.kt
    ref = g.pool(2e-3, 0.1)
    g.P["conditions"]["gamma_MA"] = 2.0
    pl = g.pool(2e-3, 0.1)
    assert math.isclose(pl["mu_MA"] - ref["mu_MA"], kt * math.log(2.0), abs_tol=1e-9)    # MA: kT ln gamma exactly
    assert math.isclose(pl["c_MX_free"], ref["c_MX_free"], rel_tol=1e-9)
    g.P["conditions"].update(gamma_MA=1.0, gamma_MX=0.5)
    pl = g.pool(2e-3, 0.1)
    assert math.isclose(pl["a_MX_free"], 0.5 * pl["c_MX_free"], rel_tol=1e-12)
    mu_from_activity = g.th["g_mx"] + kt * math.log(pl["a_MX_free"])
    assert math.isclose(pl["mu_MX"], mu_from_activity, abs_tol=1e-9)
    bound = 2e-3 * pl["mean_n"]                                       # the MX_q balance stays in concentrations
    assert math.isclose(pl["c_MX_free"] + bound, 0.1, rel_tol=1e-9)
    assert pl["mu_MX"] < ref["mu_MX"]
    _pl, fam = g.chem(2e-3, 0.1)                                      # the facets see the activity: less coverage
    _pl, fam_ref = gr.Growth(_thermo(), _params()).chem(2e-3, 0.1)
    assert all(fam[f]["theta"] <= fam_ref[f]["theta"] + 1e-12 for f in fam)


def test_surface_reaction_from_the_facet_rates():
    facets = {"100": {"k_f": 2e4}, "111A": {"k_f": 3e5}, "111B": {"k_f": 1e3}, "110": {"k_f": 7e4}}
    g = gr.Growth(_thermo(), _params(kinetics={"k_r": "facets", "facets": facets}))
    g._arrays()
    _pl, fam = g.chem(1e-3, 0.1)
    h = g.shape.wulff({f: fam[f]["gamma"] for f in g.F}, 500 * g.Om)
    m, k, area, _s = g.dot(h, fam)
    th = {f: fam[f]["theta"] for f in g.F}
    A, r, _dmu, thb = g.size_state([{"m": m, "area": area}], fam, th)
    want = sum(facets[f]["k_f"] * (1 - th[f]) * area[f] for f in g.F) / A[0]
    assert math.isclose(g._keff[0], want, rel_tol=1e-12)
    res_d, res_r = g.resistances(A, r, thb)
    assert math.isclose(res_r[0], 1.0 / want, rel_tol=1e-12)
    res = g.run()
    assert res["series"][-1]["N_dots"] > 0 and res["series"][-1]["k_eff"] > 0
    for c in gr.conservation(res):
        assert abs(c["MA"]) < 1e-9 and abs(c["MX"]) < 1e-9


SITES = {"centre": {"n": 30, "per_area": 0.03, "dG_bind": -0.6, "dW_dir": 0.05, "dW_disp": 0.9},
         "edge": {"n": 10, "per_length": 0.04, "dG_bind": -0.3, "dW_dir": 0.03}}


def _sited(**kin):
    """The synthetic params with {100} and {111}B given site classes ({111}B a single class)."""
    P = _params(kinetics={"k_r": "facets", "facets": {"100": {"k_f": 2e4}, "111A": {"k_f": 3e5}, "111B": {"k_f": 1e3},
                                                      "110": {"k_f": 7e4}}, **kin})
    P["families"]["100"].update(sites=copy.deepcopy(SITES), perimeter_per_area_ref=0.25)
    P["families"]["111B"].update(sites={"centre": {"n": 11, "per_area": 0.01, "dG_bind": -1.0, "dW_dir": 0.07,
                                                   "dW_disp": 1.2}})
    return P


def test_site_classes_are_langmuir_per_class():
    kt = gr.KB_EV * 500.0
    p = {"Gamma_max": 2.0, "gamma0": 20.0, "sites": copy.deepcopy(SITES)}
    for c in (1e-6, 1e-3, 0.1):
        s = gr.facet_state(p, c, kt)
        K = {k: math.exp(-q["dG_bind"] / kt) for k, q in SITES.items()}
        want = {k: K[k] * c / (1 + K[k] * c) for k in K}
        assert all(math.isclose(s["theta_s"][k], want[k], rel_tol=1e-12) for k in K)
        assert math.isclose(s["theta"], (30 * want["centre"] + 10 * want["edge"]) / 40, rel_tol=1e-12)
        dg = (gr.facet_state(p, c * 1.001, kt)["gamma"] - s["gamma"]) / math.log(1.001)   # Gibbs adsorption
        assert math.isclose(dg, -0.02 * kt * s["theta"], rel_tol=2e-3)
    q = {**p, "gamma_eff": 12.0, "Gamma_eff": 1.5}                                         # gamma from the dots
    assert math.isclose(gr.facet_state(q, 0.03, kt)["gamma"], 0.012 - 0.015 * kt * math.log(0.03), rel_tol=1e-12)


def test_site_rate_limits_and_inverse():
    g = gr.Growth(_thermo(), _sited())
    kt, f = g.kt, "100"
    v_c = g._speed(374.7)
    v_cx = g._speed(374.7 * 183.3 / (374.7 + 183.3))
    for th in (0.1, 0.5, 0.99):
        c = g._c_of_theta(f, th)
        assert math.isclose(gr.facet_state(g.F[f], c, kt)["theta"], th, rel_tol=1e-9)
    open_ = (0.03 * 14.9 * v_c * math.exp(-0.05 / kt) + 0.04 * 0.25 * 14.9 * v_c * math.exp(-0.03 / kt))
    assert math.isclose(g.site_rate(f, 1e-12), open_, rel_tol=1e-6)               # bare: attachment only
    full = 0.03 * 42.6 * v_cx * math.exp(-0.9 / kt)                                 # capped: the centre's exchange
    assert math.isclose(g.site_rate(f, 1.0), full, rel_tol=1e-6)
    assert math.isclose(g.site_rate(f, 1.0, p_over_a=1.0), full, rel_tol=1e-6)      # edges: no exchange, blocked
    assert g.site_rate(f, 0.5, p_over_a=0.5) > g.site_rate(f, 0.5, p_over_a=0.1)    # more perimeter, faster
    assert g.site_rate(f, 1e-12, p_over_a=1e6) == v_c                                # never above the collisions
    assert g.k_block("111A", 0.4, 3e5) == 3e5 * 0.6                                  # no classes: k_f (1 - theta)


def test_site_rates_follow_temperature_and_each_dot():
    g = gr.Growth(_thermo(), _sited())
    th = {"100": 0.1, "111A": 0.3, "111B": 0.99999, "110": 0.3}
    g.kinetic_survivors(th)
    g.set_temperature(350.0)
    fresh = gr.Growth(_thermo(), _sited())
    fresh.set_temperature(350.0)
    assert set(g.kinetic_survivors(th)) == set(fresh.kinetic_survivors(th))       # no survivors from the old T
    g.set_temperature(500.0)
    g._arrays()
    _pl, fam = g.chem(1e-3, 0.1)
    for m in (50, 5000):                                  # the reported k_eff is the one the dot grows with
        h = g.shape.wulff({f: fam[f]["gamma"] for f in g.F}, m * g.Om)
        co = {"N": 1.0, "m": m, "h": h, "area": g.shape(h)[0]}
        g.size_state([co], fam, th)
        assert math.isclose(g.mean_keff([co], th), g._keff[0], rel_tol=1e-12)


def test_gamma_eff_only_with_site_classes():
    P = _params()
    P["families"]["100"]["gamma_eff"] = 12.0              # no site classes: gamma0 stays the one in use
    P["material"]["match_m0"] = False
    g = gr.Growth(_thermo(), P)
    assert g.F["100"]["gamma0"] == P["families"]["100"]["gamma0"] and g.F["100"]["gamma_eff"] == 12.0
    S = _sited()
    S["solvation"]["eps"] = 2.0
    S["families"]["100"].update(gamma_eff=12.0, Gamma_eff=2.0)
    th = _thermo()
    for n in th["nodes"]:
        for s in n["states"]:
            s["S_vdW_A2"] = 120.0 * (n["m"] + n["k"]) ** (2 / 3)
    g = gr.Growth(th, S)
    want = 12.0 + 1000.0 * g.solv_area - 10.0 * 2.0 * g.th["dG_solv_mx"]
    assert math.isclose(g.F["100"]["gamma_eff"], want, rel_tol=1e-12)


def test_perimeter_of_a_cube():
    s = gr.Shape({"100": gr.ZB_NORMALS["100"]})
    assert math.isclose(s.perimeter({"100": 3.0})["100"], 6 * 4 * 6.0, rel_tol=1e-9)


def test_growth_with_site_classes_conserves():
    g = gr.Growth(_thermo(), _sited(growth_law="shape"))
    res = g.run()
    assert res["series"][-1]["N_dots"] > 0 and res["series"][-1]["k_eff"] > 0
    for c in gr.conservation(res):
        assert abs(c["MA"]) < 1e-9 and abs(c["MX"]) < 1e-9


def test_facet_law_with_site_classes_reduces_to_the_size_law():
    one = {"c": {"n": 1, "per_area": 0.002, "dG_bind": -0.5, "dW_dir": 0.3}}

    def growth(k_r):
        P = _params(kinetics={"growth_law": "facet", "facets": {f: {} for f in ("100", "111A", "111B", "110")},
                              "k_r": k_r})
        for p in P["families"].values():
            p["sites"] = copy.deepcopy(one)
        return gr.Growth(_thermo(), P)
    k_site = growth(1.0).site_rate("100", 0.3)
    g = growth(k_site / 0.7)                              # the size law's k_r (1 - theta) at theta = 0.3
    g._arrays()
    pl, fam = g.chem(1e-3, 0.1)
    th = {f: 0.3 for f in g.F}
    h = g.shape.wulff({f: fam[f]["gamma"] for f in g.F}, 400 * g.Om)
    m, k, area, _s = g.dot(h, fam)
    co = {"m": m, "k": k, "h": h, "area": area}
    A, r, dmu, thb = g.size_state([co], fam, th)
    C, C_M = 2e-3, 1e-3
    c_eq = C_M * np.exp((g.gb + dmu - pl["mu_MA"]) / g.kt)
    facet, v = g._facet_rate(co, C, C_M, pl["mu_MA"], fam, th)
    assert math.isclose(facet, g.size_rates(A, r, c_eq, np.array([0.3]), C)[0], rel_tol=1e-6)


def test_each_dot_has_its_own_capacitance():
    gr.SHAPE_CAPS.clear()
    g = gr.Growth(_thermo(), _params())
    g._arrays()
    _pl, fam = g.chem(1e-3, 0.1)
    g._R1 = 3.0
    h = {"100": 20.0, "111A": 19.0, "111B": 21.0, "110": 23.0}                 # all four families: near round
    m, _k, area, _s = g.dot(h, fam)
    cube = {f: (1.0 if f == "100" else 10.0) * 30.0 for f in g.F}             # only {100}: a cube
    mc, _kc, area_c, _sc = g.dot(cube, fam)
    c, con = g.dot_capacitance({"m": mc, "h": cube, "area": area_c})
    assert abs(c / 60.0 - 0.6607) < 0.01                                      # cube of edge 60 A
    assert con > c + 3.0                                                      # a cube's grown body beats C + R_1
    c, con = g.dot_capacitance({"m": m, "h": h, "area": area})
    c2, _ = g.dot_capacitance({"m": 8 * m, "h": {f: 2 * h[f] for f in g.F}, "area": area})
    assert len(gr.SHAPE_CAPS) == 2 and math.isclose(c2, 2 * c)               # one entry per shape, C ~ size
    r = (3 * m * g.Om / (4 * math.pi)) ** (1 / 3)
    assert 1.0 < c / r < 1.1 and abs(con - (c + 3.0)) / con < 0.02              # near-round: C + R_1 within 2 %
    # the shell: R_h grows by a + b k / A, the contact surface by R_1 + a' + b' k / A
    g.P["kinetics"]["shell_offset"] = {"R_h_A": [0.5, 80.0], "R_se_A": [0.2, 80.0]}
    k = 0.03 * sum(area.values())
    ch, cc = g.dot_capacitance({"m": m, "h": h, "area": area, "k": k})
    assert math.isclose(ch, c + 0.5 + 80.0 * 0.03) and abs(cc - (con + 0.2 + 2.4)) / cc < 0.02


def test_shape_capacitances_are_shared_by_runs(monkeypatch):
    gr.SHAPE_CAPS.clear()
    g = gr.Growth(_thermo(), _params())
    g._arrays()
    _pl, fam = g.chem(1e-3, 0.1)
    g._R1 = 3.0
    h = {"100": 20.0, "111A": 19.0, "111B": 21.0, "110": 23.0}
    m, _k, area, _s = g.dot(h, fam)
    want = g.dot_capacitance({"m": m, "h": h, "area": area})
    # a later run (or the same run after a temperature change, which redoes _arrays) takes them from the cache
    from qd_workflows import transport
    monkeypatch.setattr(transport, "capacitance_polyhedron", lambda *a, **k: pytest.fail("recomputed"))
    g2 = gr.Growth(_thermo(), _params())
    g2._arrays()
    g2._R1 = 3.0
    assert g2.dot_capacitance({"m": m, "h": h, "area": area}) == want


def test_shapes_within_the_key_share_a_capacitance():
    gr.SHAPE_CAPS.clear()
    g = gr.Growth(_thermo(), _params())
    g._arrays()
    _pl, fam = g.chem(1e-3, 0.1)
    g._R1 = 3.0
    h = {"100": 20.0, "111A": 19.0, "111B": 21.0, "110": 23.0}
    m, _k, area, _s = g.dot(h, fam)
    r = (3 * m * g.Om / (4 * math.pi)) ** (1 / 3)
    c, _ = g.dot_capacitance({"m": m, "h": h, "area": area})
    # {110} moved out by 0.01 r_eq: within the 0.05 key, the same C(P) per r_eq (only r_eq itself differs)
    h2 = dict(h, **{"110": h["110"] + 0.01 * r})
    c2, _ = g.dot_capacitance({"m": m, "h": h2, "area": area})
    assert len(gr.SHAPE_CAPS) == 1 and c2 == c
    # moved by 0.1 r_eq: a new shape
    g.dot_capacitance({"m": m, "h": dict(h, **{"110": h["110"] + 0.1 * r}), "area": area})
    assert len(gr.SHAPE_CAPS) == 2


def test_bulk_free_energy_table_replaces_the_einstein_estimate():
    g = gr.Growth(_thermo(), _params(material={"g_bulk_T": {"T": [300.0, 700.0], "G_eV": [-5.3, -5.7]}}))
    assert math.isclose(g.gb, -5.3 + (g.T - 300.0) / 400.0 * -0.4)
    g.set_temperature(500.0)
    assert math.isclose(g.gb, -5.5)


def test_ideal_growth_keeps_the_nucleation_and_grows_spheres():
    full = gr.Growth(_thermo(), _params())
    ideal = gr.Growth(_thermo(), _params(kinetics={"growth_law": "sphere", "sphere": {"k_r_A_s": 5e3}}))
    for g in (full, ideal):
        g._arrays()
        pl, _fam = g.chem(1e-3, 0.1)                                        # sets sigma at this pool
    J_full, b_full = full.nucleation(pl["mu_MA"], pl["mu_MX"])[:2]
    J_ideal, b_ideal = ideal.nucleation(pl["mu_MA"], pl["mu_MX"])[:2]
    assert math.isclose(J_full, J_ideal, rel_tol=1e-12) and math.isclose(b_full, b_ideal, rel_tol=1e-12)
    full.chem(*full.ref_pool())                                             # default gamma: the Wulff dot's at the reference pool
    sigma_w = full._sigma / ((36 * math.pi) ** (1 / 3) * full.Om ** (2 / 3)) / 0.0624150907
    assert math.isclose(ideal.sphere_used["gamma_J_m2"], sigma_w, rel_tol=1e-9) and ideal.sphere_used["k_r_A_s"] == 5e3
    ideal = gr.Growth(_thermo(), _params(kinetics={"growth_law": "sphere", "sphere": {"k_r_A_s": 5e3, "gamma_J_m2": 0.2}}))
    ideal._arrays()
    used = ideal.sphere_used
    assert math.isclose(used["gamma_J_m2"], 0.2)
    _pl, fam = ideal.chem(1e-3, 0.1)
    h = ideal.shape.wulff({f: fam[f]["gamma"] for f in ideal.F}, 2000 * ideal.Om)
    m, k, area, _s = ideal.dot(h, fam)
    assert k == 0.0                                                         # no MX_q shell on an ideal dot
    A, r, dmu, _th = ideal.size_state([{"m": m, "area": area}], fam, {f: 0.0 for f in ideal.F})
    assert math.isclose(r[0], (3 * m * ideal.Om / (4 * math.pi)) ** (1 / 3)) and math.isclose(A[0], 4 * math.pi * r[0] ** 2)
    gam = used["gamma_J_m2"] * 0.0624150907
    assert math.isclose(dmu[0] - ideal.junction(m), 2 * gam * ideal.Om / r[0], rel_tol=1e-12)
    assert math.isclose(ideal.resistances(A, r, _th)[1][0], 1 / 5e3)
    res = ideal.run()
    assert res["sphere"] == used and res["series"][-1]["N_dots"] > 0
    for c in gr.conservation(res):
        assert abs(c["MA"]) < 1e-9 and abs(c["MX"]) < 1e-9


def test_growth_control_and_nucleation_window():
    g = gr.Growth(_thermo(), _params())
    A, r, th = np.array([400.0]), np.array([8.0]), np.array([0.3])
    rd, rr = g.resistances(A, r, th)
    g.P["kinetics"]["growth_control"] = "diffusion"
    d, q = g.resistances(A, r, th)
    assert np.allclose(d, rd) and np.allclose(q, 0.0)                      # the surface reaction is fast
    g.P["kinetics"]["growth_control"] = "reaction"
    d, q = g.resistances(A, r, th)
    assert np.allclose(d, 0.0) and np.allclose(q, rr)                      # transport is fast
    res = gr.Growth(_thermo(), _params(conditions={"t_end": 30.0})).run()
    w = res["nucleation_window"]
    assert w and w["t10"] <= w["t50"] <= w["t90"] <= 30.0 and w["nuclei_mol"] > 0
    assert abs(res["series"][-1]["nucleated_M"] * res["series"][-1]["V_L"] - w["nuclei_mol"]) < 1e-12


def test_a_fast_facet_grows_out_of_the_surface_reaction():
    facets = {"100": {"k_f": 1e4}, "111A": {"k_f": 4e11}, "111B": {"k_f": 3e3}, "110": {"k_f": 2e11}}
    g = gr.Growth(_thermo(), _params(kinetics={"k_r": "facets", "facets": facets}))
    th = {f: 0.0 for f in g.F}
    assert g.kinetic_survivors(th) == {"111B"}                               # the slowest family alone
    # a cube with a thin {111}A sliver: the sliver adds nothing, the slowest exposed facet ({100}) sets the rate
    area = {"100": 990.0, "111A": 10.0, "111B": 0.0, "110": 0.0}
    assert math.isclose(g.facet_keff(area, th), 1e4)
    # {100} and {111}B together: {100} grows out, {111}B alone counts
    area = {"100": 800.0, "111A": 10.0, "111B": 200.0, "110": 0.0}
    assert math.isclose(g.facet_keff(area, th), 3e3)
    # comparable rates stay together: area-weighted over both
    g2 = gr.Growth(_thermo(), _params(kinetics={"k_r": "facets", "facets": {**facets, "100": {"k_f": 4e3}}}))
    assert g2.kinetic_survivors(th) == {"100", "111B"}
    assert math.isclose(g2.facet_keff(area, th), (4e3 * 800 + 3e3 * 200) / 1000)


def test_growth_does_not_change_its_params():
    P = _params()
    before = copy.deepcopy(P)
    gr.Growth(_thermo(), P)
    gr.Growth(_thermo(), P)
    assert P == before


def test_injection_dilutes_and_conserves_amounts():
    inj = {"c_syringe_M": 0.3, "rate_L_h": 0.36, "V0_L": 0.01, "duration_s": 20.0}       # 2 mL over 20 s into 10 mL
    P = _params(conditions={"c_P0": 0.0, "t_end": 40.0, "n_out": 8, "injection": inj})
    res = gr.Growth(_thermo(), P).run()
    last = res["series"][-1]
    assert math.isclose(last["V_L"], 0.012, rel_tol=1e-9)
    assert math.isclose(last["n_inj_mol"], 0.3 * 0.002, rel_tol=1e-9)
    for c in gr.conservation(res):
        assert abs(c["MA"]) < 1e-9 and abs(c["MX"]) < 1e-9
    assert res["series"][0]["P"] < 0.3 * 0.002 / 0.012                                   # the precursor arrives slowly


def test_temperature_program_switches_the_free_energies():
    P = _params(conditions={"t_end": 30.0, "n_out": 6, "T_growth": 450.0, "t_switch_s": 5.0, "ramp_s": 10.0})
    g = gr.Growth(_thermo(), P)
    assert g.temperature_at(1.0) == P["conditions"]["T"] and g.temperature_at(30.0) == 450.0
    assert math.isclose(g.temperature_at(10.0), 0.5 * (P["conditions"]["T"] + 450.0))
    res = g.run()
    Ts = [r["T"] for r in res["series"]]
    assert Ts[0] == P["conditions"]["T"] and abs(Ts[-1] - 450.0) < 1.0 and g.T == Ts[-1]
    ref = gr.Growth(_thermo(), _params(conditions={"T": 450.0}))
    assert math.isclose(g.gb, ref.gb, abs_tol=1e-3) and math.isclose(g.k_conv, ref.k_conv, rel_tol=1e-2)
    for c in gr.conservation(res):
        assert abs(c["MA"]) < 1e-9 and abs(c["MX"]) < 1e-9


def test_results_do_not_depend_on_the_output_times():
    sizes = []
    for n_out in (8, 40):
        res = gr.Growth(_thermo(), _params(conditions={"t_end": 30.0, "n_out": n_out})).run()
        sizes.append((res["series"][-1]["N_dots"], res["series"][-1]["m_mean"]))
    assert math.isclose(sizes[0][0], sizes[1][0], rel_tol=0.05) and math.isclose(sizes[0][1], sizes[1][1], rel_tol=0.05)


def test_atomistic_transport_enters_flux_and_dots():
    th = _thermo()
    for nd in th["nodes"]:                                            # geometry of every node: radii grow with size
        for s in nd["states"]:
            s["R_se_A"] = 2.4 * nd["m"] ** (1 / 3) + 0.3 * nd["k"]
            s["R_h_A"] = 2.6 * nd["m"] ** (1 / 3) + 0.3 * nd["k"]
    g = gr.Growth(th, _params(kinetics={"transport": "stokes", "viscosity_mPa_s": 0.4}))
    g._arrays()
    assert g.stokes
    pl, fam = g.chem(2e-3, 0.1)
    J1 = g.nucleation(pl["mu_MA"] + 0.3, pl["mu_MX"])[0]
    assert g._R1 > 2.4 and g._D1 > 0
    g2 = gr.Growth(th, _params(kinetics={"transport": "stokes", "viscosity_mPa_s": 0.8}))
    g2._arrays()
    g2.chem(2e-3, 0.1)
    J2 = g2.nucleation(pl["mu_MA"] + 0.3, pl["mu_MX"])[0]
    assert abs(J1 / J2 - 2.0) < 1e-6                                   # diffusion-limited: J ~ D ~ 1/eta
    g._capf = g._cap_factor(fam)                                       # a dot: its Wulff polyhedron's shape factor
    assert 1.0 < g._capf < 1.2
    plain = gr.Growth(_thermo(), _params(kinetics={"transport": "stokes"}))
    plain._arrays()
    assert not plain.stokes                                             # no geometry: the fixed D


def test_steps_resolve_the_pool_potentials():
    # no step moves the pool's mu_MA or mu_MX by more than dmu_step_kT kT (beyond the floor dt): the record of the
    # step-end pools resolves them; without the limit the same run takes larger jumps
    def run(dmu):
        g = gr.Growth(_thermo(), _params(material={"dmu_step_kT": dmu, "dt_rel_max": 0.05}))
        g.run()
        mus = [g.pool(cm, cx) for _t, _T, cm, cx in g.pools]
        return g, max(max(abs(b["mu_MA"] - a["mu_MA"]), abs(b["mu_MX"] - a["mu_MX"])) / g.kt for a, b in zip(mus, mus[1:]))
    g, jump = run(0.1)
    _g0, jump0 = run(0.0)
    assert len(g.pools) == len(g.steps) and all(a[0] < b[0] for a, b in zip(g.pools, g.pools[1:]))
    assert jump0 > 0.1 and jump < jump0


def test_pool_caps_the_free_mx_at_its_solubility():
    from qd_workflows.grid_thermo import pool, thermo_at
    th = thermo_at(_thermo(), 500.0)
    kt = gr.KB_EV * 500.0
    free = pool(th, 0.01, 0.2)
    assert free["c_MX_solid"] == 0.0
    same = pool(th, 0.01, 0.2, c_sat=10 * free["c_MX_free"])                 # above the balance: no solid
    assert math.isclose(same["mu_MX"], free["mu_MX"], rel_tol=1e-12) and same["c_MX_solid"] == 0.0
    sat = 0.1 * free["c_MX_free"]
    cap = pool(th, 0.01, 0.2, c_sat=sat)
    assert math.isclose(cap["c_MX_free"], sat) and cap["c_MX_solid"] > 0
    assert math.isclose(cap["mu_MX"], th["g_mx"] + kt * math.log(sat), rel_tol=1e-12)
    assert math.isclose(cap["c_MX_free"] + 0.01 * cap["mean_n"] + cap["c_MX_solid"], 0.2, rel_tol=1e-9)
    bare = pool(th, 0.01, 0.2, bare=True, c_sat=1e-3)
    assert math.isclose(bare["c_MX_free"], 1e-3) and math.isclose(bare["c_MX_solid"], 0.2 - 1e-3)


def test_solubility_from_the_solid_and_mass_balance():
    g = gr.Growth(_thermo(), _params())
    assert g.c_sat is None
    g_s = g.th["g_mx"] + g.kt * math.log(1e-3)                               # the solid at which 1 mM saturates
    tab = {"T": [300.0, 700.0], "G_eV": [g_s, g_s]}
    assert math.isclose(gr.Growth(_thermo(), _params(material={"g_MX_solid_T": tab})).c_sat, 1e-3, rel_tol=1e-9)
    both = gr.Growth(_thermo(), _params(material={"g_MX_solid_T": tab}, conditions={"c_MX_sat": 2e-3}))
    assert both.c_sat == 2e-3                                                  # a measured solubility wins
    res = gr.Growth(_thermo(), _params(nucleation={"J0": 1e-4}, conditions={"c_MX_sat": 1e-3})).run()
    assert res["series"][-1]["N_dots"] > 0
    assert all(r["c_MX_free"] <= 1e-3 * (1 + 1e-9) for r in res["series"])
    assert any(r["c_MX_solid"] > 0 for r in res["series"])
    for c in gr.conservation(res):
        assert abs(c["MA"]) < 1e-6 and abs(c["MX"]) < 1e-6


ROLES = {"rate_form": "second_order",
         "envs": {"facet:100": {20: [[8, 0.5, -0.9]], 60: [[20, 0.6, -0.8]]},
                  "edge:100|111B": {20: [[6, 0.4, -1.2], [2, 0.7, -0.6]], 60: [[12, 0.45, -1.1]]},
                  "vertex:100+111A+111B": {20: [[4, 0.3, -0.3]], 60: [[4, 0.35, -0.2]]},
                  "edge:111A|111B": {60: [[3, 0.2, -0.5]]}},
         "neb": {"facet:100": {"dir": {"no_barrier": True}, "disp": {"dG_P": -0.3, "dG_TS": 0.45}},
                 "edge:100|111B": {"dir": {"dG_P": -1.3, "dG_TS": 0.5}, "disp": {"dG_P": -0.6, "dG_TS": 0.55}},
                 "vertex:100+111A+111B": {"dir": {"no_barrier": True}, "disp": {"dG_P": -1.9, "dG_TS": 0.4}}}}


def _per_facet(growth_law="shape", **roles):
    return _params(kinetics={"size_law": "facets", "growth_law": growth_law, "roles": {**copy.deepcopy(ROLES), **roles}})


def test_role_environments_follow_the_computed_dots():
    g = gr.Growth(_thermo(), _per_facet())
    at = g.role_envs(60)
    assert sorted((r, n) for r, n, _d, _b in at) == sorted(
        [("facet:100", 20), ("edge:100|111B", 12), ("vertex:100+111A+111B", 4), ("edge:111A|111B", 3)])
    # between the sizes: both sizes' environments, weighted linearly in m^(-1/3)
    t = (20 ** (-1 / 3) - 35 ** (-1 / 3)) / (20 ** (-1 / 3) - 60 ** (-1 / 3))
    mid = g.role_envs(35)
    assert math.isclose(sum(n for r, n, *_ in mid if r == "facet:100"), (1 - t) * 8 + t * 20, rel_tol=1e-12)
    assert math.isclose(sum(n for r, n, *_ in mid if r == "edge:111A|111B"), t * 3, rel_tol=1e-12)
    # beyond the largest: facets as m^(2/3), edges m^(1/3), vertices fixed
    big = {r: n for r, n, *_ in g.role_envs(480)}
    assert math.isclose(big["facet:100"], 20 * 4.0, rel_tol=1e-12)
    assert math.isclose(big["edge:100|111B"], 12 * 2.0, rel_tol=1e-12)
    assert big["vertex:100+111A+111B"] == 4


def test_per_facet_law_is_the_channel_sum():
    """dm/dt = n sum_r g_r (C - C_eq,inf e^(Delta mu_r/kT)), g_r = N_r w_rf / (res_d / a_0 + 1/kappa_r), with the role
    channels by hand: N_r the role's Cd count, Delta mu_r = <dmu>_r + kT ln a_X (plain Cd-weighted mean), kappa_r =
    (kappa_dir + kappa_disp <K_e>_r a_X) / (1 + K_eff,f a_X)."""
    g = gr.Growth(_thermo(), _per_facet())
    g._arrays()
    C_M, C_X = 1e-3, 0.05
    pl, fam = g.chem(C_M, C_X)
    th = {f: fam[f]["theta"] for f in g.F}
    h = {"100": 8.0, "111A": 9.0, "111B": 10.0, "110": 50.0}
    h = {f: x * (35 * g.Om / g.shape(h)[1]) ** (1 / 3) for f, x in h.items()}        # a dot of 35 units
    m, k, area, _s = g.dot(h, fam)
    co = {"m": m, "k": k, "h": h, "area": area}
    A, r, dmu, thb = g.size_state([co], fam, th)
    res_d, _ = g.resistances(A, r, thb)
    kt, n1, ax = g.kt, 6.02214076e-4, pl["a_MX_free"]
    envs = g.role_envs(m)
    a0 = A[0] / sum(e[1] for e in envs)
    coll = a0 * g._speed(191.385)

    def kap(q):
        return coll if q.get("no_barrier") else min((kt / 4.135667696e-15) * math.exp(-q["dG_TS"] / kt) / n1, coll)
    num, den = {}, {}
    for role, n, _d, b in envs:
        for f, w in gr.Growth.role_facets(role).items():
            num[f] = num.get(f, 0.0) + w * n * math.exp(-b / kt)
            den[f] = den.get(f, 0.0) + w * n
    per = {}
    for role, n, d, b in envs:
        per.setdefault(role, []).append((n, d, math.exp(-b / kt)))
    C, want = 2e-3, 0.0
    for role, rows in per.items():
        if role not in ROLES["neb"]:
            continue                                          # no NEB: no channel
        q = ROLES["neb"][role]
        n = sum(x[0] for x in rows)
        d = sum(x[0] * x[1] for x in rows) / n
        K = sum(x[0] * x[2] for x in rows) / n
        for f, w in gr.Growth.role_facets(role).items():
            kk = (kap(q["dir"]) + kap(q["disp"]) * K * ax) / (1 + num[f] / den[f] * ax)
            ceq = C_M * math.exp((g.gb + d + kt * math.log(ax) - pl["mu_MA"]) / kt)
            want += n1 * w * n / (res_d[0] / a0 + 1 / kk) * (C - ceq)
    c_eq = C_M * np.exp((g.gb + dmu - pl["mu_MA"]) / kt)
    assert math.isclose(g.size_rates(A, r, c_eq, thb, C)[0], want, rel_tol=1e-9)
    assert "edge:111A|111B: no NEB, no channel" in g.role_flags


def test_role_channel_takes_the_plain_mean_delta_mu():
    """One channel per role and facet; Delta mu_r is the Cd-weighted mean of the role's environments, not a Boltzmann
    mean (edge:100|111B at size 20: 6 Cd at 0.4 and 2 at 0.7 -> 0.475)."""
    g = gr.Growth(_thermo(), _per_facet())
    pl = dict(g.chem(1e-3, 0.05)[0], a_MX_free=1.0)
    ch = g.role_channels(20, 1000.0, pl)
    i = [j for j, r in enumerate(ch["role"]) if r == "edge:100|111B"]
    assert len(i) == 2 and sorted(ch["f"][j] for j in i) == ["100", "111B"]
    assert ch["N"][i].sum() == pytest.approx(8.0)
    assert ch["dmu"][i] == pytest.approx(0.475 + g.th["dG_solv_mx"])


def test_displacement_keeps_a_capped_site_open():
    """At high a_X the empty sites close (1 / (1 + K_eff a_X)), but the displacement route keeps kappa finite."""
    g = gr.Growth(_thermo(), _per_facet())
    pl = dict(g.chem(1e-3, 0.05)[0])

    def kappa(ax, disp=True):
        p = dict(pl, a_MX_free=ax)
        if not disp:
            g.P["kinetics"]["roles"]["neb"]["edge:100|111B"]["disp"] = {"dG_P": 0.0, "dG_TS": 50.0}
        ch = g.role_channels(60, 3000.0, p)
        return ch["kappa"][[i for i, f in enumerate(ch["f"]) if f == "111B"][0]]
    hi, lo = kappa(1e3), kappa(1e2)
    assert hi / lo == pytest.approx(1.0, rel=1e-3)
    assert kappa(1e3, disp=False) < 1e-6 * hi


def test_per_facet_law_runs_and_conserves():
    for law in ("size", "shape"):
        P = _per_facet(growth_law=law)
        P["conditions"]["t_end"] = 5.0
        res = gr.Growth(_thermo(), P).run()
        assert res["series"][-1]["N_dots"] > 0
        for c in gr.conservation(res):
            assert abs(c["MA"]) < 1e-9 and abs(c["MX"]) < 1e-9


def test_per_facet_law_needs_a_rate_form():
    P = _per_facet(rate_form=None)
    g = gr.Growth(_thermo(), P)
    with pytest.raises(ValueError, match="rate_form"):
        g.role_kappa(g.chem(1e-3, 0.05)[0], 10.0)


def test_template_roles_are_the_computed_dots():
    import yaml
    from importlib import resources
    P = yaml.safe_load(resources.files("qd_workflows").joinpath("templates/growth/cdse_zb.yaml").read_text())
    R = P["kinetics"]["roles"]
    # the per-role block stays as the computed data; the template's law is the facet-area one built from it
    assert R["rate_form"] is None and P["kinetics"]["size_law"] == "facet_area"
    assert {int(s) for r in R["envs"].values() for s in r} == {13, 45, 127, 279, 441}
    assert set(R["neb"]) <= set(R["envs"])
    fams = set(P["families"])
    assert all(set(gr.Growth.role_facets(r)) <= fams for r in R["envs"])
