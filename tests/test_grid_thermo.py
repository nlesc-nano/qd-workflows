"""Grand-potential analysis of a grid (grid_thermo.omega) on synthetic thermochemistry."""
import numpy as np

from qd_workflows.grid_thermo import KB_EV, converted, k_conv_at, node_g, omega, unit_names

T = [100.0, 300.0, 500.0]


def _node(m, k, energies, g_shift=0.1):
    states = [{"fingerprint": f"{m}-{k}-{i}", "E_eV": e, "G_eV": [e + g_shift * t / 300 for t in T],
               "file": f"nodes/m{m:02d}_k{k:02d}.xyz"} for i, e in enumerate(energies[:1])]
    return {"m": m, "k": k, "formula": f"Cd{m + k}Se{m}" + (f"Cl{2 * k}" if k else ""),
            "ensemble_E_eV": sorted(energies), "states": states}


def _thermo():
    nodes = [_node(1, 0, [-2.2]), _node(1, 1, [-9.5]), _node(1, 2, [-16.6, -16.5]), _node(1, 3, [-23.4]),
             _node(2, 0, [-6.0]), _node(2, 1, [-13.8, -13.75]), _node(2, 2, [-21.0])]
    return {"temperatures": T, "nodes": nodes,
            "references": {"CdCl2": {"G_eV": [-6.76, -6.86, -6.97]}, "CdSe": {"G_eV": [-2.3, -2.4, -2.5]}}}


def test_unit_names():
    assert unit_names(_thermo()["nodes"]) == ["CdCl2", "CdSe"]


def test_node_free_energy_sums_the_isomers():
    n = _node(2, 1, [-13.8, -13.75])
    kt = KB_EV * 300
    g_best = -13.8 + 0.1
    expect = g_best - kt * np.log(1 + np.exp(-0.05 / kt))
    assert abs(node_g(n, 300.0, T) - expect) < 1e-9


def test_precursors_are_the_m1_row_and_add_up():
    th = _thermo()
    r = omega(th, 300.0, c_p0=1e-2, c_mx0=5e-2, alpha=0.1)
    kt = KB_EV * 300
    assert abs(sum(r["precursor_speciation"].values()) - 1) < 1e-12
    c1 = {k: np.exp(-r["Omega_eV"][f"1,{k}"] / kt) for k in range(4)}        # c(1,n)/c°
    assert abs(sum(c1.values()) - 1e-3) < 1e-12                              # the MA supply alpha c_P0
    for n, x in r["precursor_speciation"].items():
        assert abs(c1[n] / 1e-3 - x) < 1e-9
    assert set(r["per_m"]) == {1, 2} and r["per_m"][2]["k_star"] in (0, 1, 2)


def test_mx_balance_after_monomer_formation():
    """Free MX + MX in the complexes = initial MX - MX consumed by the conversion."""
    th = _thermo()
    r = omega(th, 300.0, c_p0=1e-2, c_mx0=5e-2, alpha=0.3)
    in_complexes = 3e-3 * sum(n * x for n, x in r["precursor_speciation"].items())
    assert abs(r["c_MX_free_M"] + in_complexes - (5e-2 - 3e-3)) < 1e-9
    r2 = omega(th, 300.0, c_p0=1e-2, c_mx0=5e-2, alpha=0.9)                 # more conversion: less MX left
    assert r2["mu_MX_eV"] < r["mu_MX_eV"] and r2["c_MA_M"] > r["c_MA_M"]


def _integrate(k, a0, b0, t, n=200000):
    """dx/dt = k (a0 - x)(b0 - x) by RK4, the reference for the closed form."""
    h, x = t / n, 0.0
    f = lambda y: k * (a0 - y) * (b0 - y)
    for _ in range(n):
        k1 = f(x); k2 = f(x + h * k1 / 2); k3 = f(x + h * k2 / 2); k4 = f(x + h * k3)
        x += h * (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return x


def test_conversion_closed_form_matches_the_rate_law():
    for a0, b0 in ((5e-2, 1e-2), (1e-2, 5e-2), (2e-2, 2e-2)):
        for t in (10.0, 1e3, 1e5):
            assert abs(converted(t, 0.5, b0, a0) - _integrate(0.5, a0, b0, t)) < 1e-9
    assert abs(converted(1e9, 1.0, 1e-2, 5e-2) - 1e-2) < 1e-12        # the limiting reagent is used up
    assert converted(0.0, 1.0, 1e-2, 5e-2) == 0.0


def test_arrhenius_scaling_of_k_conv():
    assert k_conv_at(500.0, 2.0, 0.5) == 2.0
    assert abs(k_conv_at(600.0, 2.0, 0.5) / 2.0 - np.exp(-0.5 / KB_EV * (1 / 600 - 1 / 500))) < 1e-12


def test_bare_monomer_is_the_pre_nucleation_formula():
    th = _thermo()
    r = omega(th, 300.0, 1e-2, 5e-2, 0.2, bare=True)
    kt = KB_EV * 300
    g10 = node_g(th["nodes"][0], 300.0, T)                              # the m = 1, k = 0 node: bare CdSe
    g_mx = -6.86
    for key in ("2,1", "1,2"):
        m, k = map(int, key.split(","))
        expect = r["G_eV"][key] - m * (g10 + kt * np.log(2e-3)) - k * (g_mx + kt * np.log(5e-2 - 2e-3))
        assert abs(r["Omega_eV"][key] - expect) < 1e-9


def test_page_embeds_the_thermochemistry(tmp_path):
    import json
    from qd_workflows.grid_thermo import page
    th = _thermo()
    th["vibrations"] = False
    src = tmp_path / "thermo.json"
    src.write_text(json.dumps(th))
    html = (tmp_path / "p.html").read_text() if page(src, tmp_path / "p.html", "Test Grid") else ""
    assert "<title>Test Grid</title>" in html and "{{DATA}}" not in html
    data = json.loads(html.split("const D = ", 1)[1].split(";\nconst KB", 1)[0])
    assert data["units"] == {"MX": "CdCl2", "MA": "CdSe"} and data["vibrations"] is False and len(data["nodes"]) == 7


def test_ligand_binding_lowers_mu_mx():
    th = _thermo()
    free = omega(th, 300.0, 1e-2, 5e-2, 0.1)
    bound = omega(th, 300.0, 1e-2, 5e-2, 0.1, k_l=1e4, c_l=1.0)
    assert bound["mu_MX_eV"] < free["mu_MX_eV"] - 0.1
    assert abs(bound["c_MX_free_M"] * (1 + 1e4) + bound["c_MX_in_complexes_M"] - (5e-2 - 1e-3)) < 1e-9
