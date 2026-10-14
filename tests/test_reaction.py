"""qdw reaction: equation parsing, balancing, the eps interpolation and the reaction formulas."""
import math

import numpy as np
import pytest

from qd_workflows.reaction import model, solvation, spec


def test_parse_equation():
    assert spec.parse_equation("origin + 12 HCl -> QD_Cl + 4 AsH3 + 8 CH3NH2") == {
        "origin": -1, "HCl": -12, "QD_Cl": 1, "AsH3": 4, "CH3NH2": 8}
    with pytest.raises(ValueError):
        spec.parse_equation("a + b")


def test_suggest_coefficients():
    sp = {"o": {"composition": {"C": 21, "H": 105, "N": 21, "In": 52, "As": 35, "Cl": 51}, "charge": 0},
          "HCl": {"composition": {"H": 1, "Cl": 1}, "charge": 0},
          "q": {"composition": {"C": 13, "H": 65, "N": 13, "In": 52, "As": 31, "Cl": 63}, "charge": 0},
          "AsH3": {"composition": {"As": 1, "H": 3}, "charge": 0},
          "MA": {"composition": {"C": 1, "H": 5, "N": 1}, "charge": 0}}
    assert spec.suggest_coefficients(sp, ["o", "HCl"], ["q", "AsH3", "MA"]) == {"o": 1, "HCl": 12, "q": 1, "AsH3": 4, "MA": 8}


def test_interp_f_and_nonel():
    knots = [(0.0, 0.0), (0.5, -1.0), (0.9875, -2.0)]
    assert solvation.interp_f(knots, 2.0) == pytest.approx(-1.0)
    assert solvation.interp_f(knots, 80.0) == pytest.approx(-2.0)
    assert solvation.interp_f(knots, 1.0) == pytest.approx(0.0)
    assert solvation.nonel(1.0, 100.0, 100.0) == 0.0
    assert math.isnan(solvation.nonel(1.5, 100.0, 100.0))
    # 36 + 3 ln 2.38 mN/m on 100 A^2, -0.5 GPa on 100 A^3
    assert solvation.nonel(2.38, 100.0, 100.0) == pytest.approx(
        (36.0 + 3.0 * math.log(2.38)) * 6.241509e-5 * 100 - 0.5 * 6.241509e-3 * 100, rel=1e-6)


def _species(e, c, c0=None):
    T = [300.0, 400.0]
    th = {k: [0.0, 0.0] for k in ("zpe", "dU_vib", "H_trans", "H_rot", "S_trans", "S_rot", "S_vib", "S_el", "G_corr")}
    s = {"energy_eV": e, "thermo": th, "phase": "solute", "c": c, "gb_inf_eV": 0.0, "sccs": None}
    if c0 is not None:
        s["c0"] = c0
    return s


def test_reaction_and_equilibrium_a_to_b():
    res = {"T": [300.0, 400.0], "nu": {"A": -1, "B": 1}, "per": None,
           "species": {"A": _species(0.0, 1e-2), "B": _species(-0.05, 1e-3)}}
    cond = {"eps": 1.0, "gas": True}
    r = model.reaction(res, cond)
    kT = model.KB_EV * np.array(res["T"])
    assert np.allclose(r["dG0"], -0.05)
    assert np.allclose(r["dG"], -0.05 + kT * math.log(1e-3 / 1e-2))
    eq = model.equilibrium(res, cond)
    K = np.exp(0.05 / kT)                       # B / A = K, A + B = 1e-2
    assert np.allclose(eq["conversion"], K / (1 + K), rtol=1e-9)


def test_reaction_a_to_a_is_zero():
    res = {"T": [300.0, 400.0], "nu": {"A": -1, "A2": 1}, "per": None,
           "species": {"A": _species(-3.0, 1e-2), "A2": _species(-3.0, 1e-2)}}
    assert np.allclose(model.reaction(res, {"eps": 1.0, "gas": True})["dG"], 0.0)


def test_ligand_term_and_gas_phase():
    res = {"T": [300.0, 400.0], "nu": {"A": -1, "B": 1, "L": 2}, "per": None, "ligand": {"species": "L", "dG_eV": 0.1},
           "species": {"A": _species(0.0, 1e-2), "B": _species(-0.5, 1e-3), "L": {**_species(0.0, 1.0), "phase": "gas",
                                                                                   "p_bar": 1e-3, "gb_inf_eV": -9.0}}}
    T = np.array(res["T"])
    kT = model.KB_EV * T
    r = model.reaction(res, {"eps": 4.0, "electrostatics": "gb", "nonel": False})
    # gas species: no solvation term (GB would give (1 - 1/4) * -9 eV), activity p / (c0 R T)
    assert np.allclose(r["parts"]["ddG_solv_el"], 0.0)
    assert np.allclose(r["parts"]["ligand"], 0.2)
    lna_l = np.log(1e-3 * 1e5 / (model.C0 * model.R_GAS * T))
    assert np.allclose(r["dG"], -0.5 + 0.2 + kT * (math.log(1e-3 / 1e-2) + 2 * lna_l))
    assert np.allclose(model.reaction(res, {"eps": 1.0, "gas": True, "dG_lig": 0.0})["parts"]["ligand"], 0.0)
