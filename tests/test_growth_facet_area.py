"""Facet-area size law (kinetics.size_law = facet_area): one channel per facet family on its area, terrace K_X and
Delta mu, the step band along the edges, signed facet moves, no surface energies; the dot's geometry from the model's
polyhedron here (QD_builder's on-the-fly geometry is checked when a builder with lattice_fitted is importable)."""
import copy
import math
from pathlib import Path

import numpy as np
import pytest
import yaml

pytest.importorskip("scipy")
from qd_workflows import growth as gr  # noqa: E402
from test_growth import _params, _thermo  # noqa: E402

TEMPLATE = Path(gr.__file__).parent / "templates" / "growth" / "cdse_zb.yaml"


def _fa_params(**over):
    P = _params(**over)
    T = yaml.safe_load(TEMPLATE.read_text())
    P["kinetics"]["size_law"] = "facet_area"
    P["kinetics"]["facet_area"] = copy.deepcopy(T["kinetics"]["facet_area"])
    P["kinetics"]["facet_area"]["geometry"] = "shape"
    P["material"]["m0"] = 13
    return P


def test_template_block_is_the_terrace_and_edge_data():
    fa = yaml.safe_load(TEMPLATE.read_text())["kinetics"]["facet_area"]
    assert set(fa["facets"]) == set(gr.FA_FAMILIES)
    assert fa["facets"]["100"]["dG_bind_eV"] == pytest.approx(-0.938, abs=1e-3)      # Cd-weighted over sizes
    assert fa["facets"]["111B"]["dG_bind_eV"] == pytest.approx(-1.800, abs=1e-3)
    assert set(fa["facets"]["100"]["dmu_eV"]) == {127, 279, 441}                   # terrace sizes only
    assert fa["edges"]["111A|111B"]["neb"] is None                                  # no borrowed value


def test_edge_lengths_of_a_cube():
    sh = gr.Shape({f: gr.ZB_NORMALS[f] for f in gr.FA_FAMILIES})
    e = sh.edge_lengths({"100": 1.0, "111A": 10.0, "111B": 10.0})
    assert e == {"100|100": pytest.approx(24.0)}
    e = sh.edge_lengths({"100": 1.0, "111A": 1.6, "111B": 1.6})                    # corners cut by {111}
    assert set(e) == {"100|100", "100|111A", "100|111B"}


def test_no_gamma_in_the_facet_area_families():
    g = gr.Growth(_thermo(), _fa_params())
    assert set(g.F) == set(gr.FA_FAMILIES) and g.handover is None and g.m0 == 13
    for f, p in g.F.items():
        assert p["w"] == 0.0 and p["Gamma_max"] * 1e-2 == pytest.approx(1.0 / g.fa_site_area(f))
    g.chem(1e-3, 0.03)
    assert g._sigma is None


def test_dmu_is_held_outside_the_computed_sizes():
    g = gr.Growth(_thermo(), _fa_params())
    tab = g._fa()["facets"]["100"]["dmu_eV"]
    assert g.fa_dmu("100", 20) == pytest.approx(tab[127])
    assert g.fa_dmu("100", 2000) == pytest.approx(tab[441])
    x = lambda m: m ** (-1 / 3)                                                      # noqa: E731
    mid = g.fa_dmu("100", 200)
    t = (x(127) - x(200)) / (x(127) - x(279))
    assert mid == pytest.approx(tab[127] + t * (tab[279] - tab[127]))


def test_facet_rate_is_bounded_and_coverage_below_one():
    g = gr.Growth(_thermo(), _fa_params())
    g.chem(1e-3, 0.03)
    v1 = g._speed(191.385)
    h = g.fa_start_h(300.0)
    co = {"m": 300.0, "h": h, "area": g.shape(h)[0], "edges": g.shape.edge_lengths(h)}
    ch = g._fa_channels(co)
    assert np.all(ch["kappa"] > 0) and np.all(ch["kappa"] <= v1 * (1 + 1e-12))
    for ax in (1e-6, 1e-2, 1.0):
        for f in gr.FA_FAMILIES:
            th = gr.facet_state(g.F[f], ax, g.kt)["theta"]
            assert 0.0 <= th <= 1.0
            assert g.fa_keff(f, {f: 1e-6}, co["edges"], ax, 1e-18) <= v1 * (1 + 1e-12)   # tiny facet: band capped


def test_mass_is_conserved_and_no_gamma_is_read():
    P = _fa_params(nucleation={"J0": 1e-4}, conditions={"t_end": 5.0, "n_out": 6})
    for p in P["families"].values():                        # the facet-area law must not need these
        for k in ("gamma0", "eps0", "Gamma_max", "w"):
            p[k] = float("nan")
    res = gr.Growth(_thermo(), P).run()
    last = gr.conservation(res)[-1]
    assert abs(last["MA"]) < 1e-10 and abs(last["MX"]) < 1e-10
    assert res["geometry"]["source"] == "shape"


def test_signed_step_keeps_the_volume():
    g = gr.Growth(_thermo(), _fa_params())
    h = g.fa_start_h(200.0)
    co = {"m": 200.0, "h": h}
    area = g.shape(h)[0]
    v = {f: (1.0 if f == "111A" else -0.5) * (1.0 if area[f] > 0 else 0.0) for f in gr.FA_FAMILIES}
    hn = g._fa_step(co, 210.0, v)
    assert g.shape(hn)[1] == pytest.approx(210.0 * g.Om, rel=1e-8)
    for f in gr.FA_FAMILIES:
        if area[f] <= 0:
            assert hn[f] == h[f]                             # hidden planes held


def _builder_with_lattice_fitted():
    try:
        from builder import main as bm
    except Exception:
        return False
    return "_lattice_fitted" in dir(bm)


@pytest.mark.skipif(not _builder_with_lattice_fitted(), reason="QD_builder without facet_geometry.lattice_fitted")
def test_builder_geometry_reproduces_the_handover_dot():
    from qd_workflows.builder_geometry import BuilderGeometry
    fa = yaml.safe_load(TEMPLATE.read_text())["kinetics"]["facet_area"]
    bg = BuilderGeometry(fa["builder"], gr.FA_FAMILIES)
    area, edges, b = bg.get(fa["start_h_A"], 13)
    assert b["m"] == 13 and b["h_err_A"] < 0.45
    assert area["100"] == 0.0 and area["111A"] == pytest.approx(100.0, abs=1.0)
    assert edges["111A|111B"] == pytest.approx(91.2, abs=0.5)
