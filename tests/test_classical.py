"""Classical twin: CNT of spheres, mass conservation with the syringe and the temperature program, and the LSW
(diffusion-limited) and Wagner (reaction-limited) ripening rate constants."""
import copy
import math

import numpy as np

from qd_workflows import classical as cl


def _p(**over):
    P = copy.deepcopy(cl.PARAMS)
    for k, v in over.items():
        P[k].update(v)
    return P


def test_cnt_of_spheres():
    m = cl.Classical(_p())
    C = 50 * m.c_inf
    nc = m.cnt(C)
    L = m.kt * math.log(50)
    assert math.isclose(nc["r_star"], 2 * m.gam * m.Vm / L, rel_tol=1e-12)
    assert math.isclose(nc["dG_star"], 16 * math.pi * m.gam ** 3 * m.Vm ** 2 / (3 * L * L), rel_tol=1e-12)
    assert math.isclose(float(m.radius(nc["n_star"])), nc["r_star"], rel_tol=1e-12)    # n* is the sphere of r*
    dG = lambda n: -n * L + m.gam * m.a * n ** (2 / 3)                                # noqa: E731
    assert math.isclose(dG(nc["n_star"]), nc["dG_star"], rel_tol=1e-12)
    assert math.isclose(dG(nc["n_birth"]), nc["dG_star"] - m.kt, rel_tol=1e-9) and nc["n_birth"] > nc["n_star"]
    assert m.cnt(0.5 * m.c_inf)["J"] == 0.0                                            # undersaturated: none
    assert m.cnt(100 * m.c_inf)["J"] > nc["J"]


def test_growth_rate_limits():
    m = cl.Classical(_p(material={"k_r_A_s": 1e14}))                                   # diffusion limit
    r = 30.0
    n = 4 * math.pi * r ** 3 / (3 * m.Vm)
    C = 10 * m.c_inf
    gt = math.exp(2 * m.gam * m.Vm / (r * m.kt))
    want = 4 * math.pi * r * m.D * cl.NA * (C - m.c_inf * gt)
    assert math.isclose(float(m.rate(n, C)), want, rel_tol=1e-4)                      # k_r finite: 3e-5


def test_growth_control_switches_the_limiting_step():
    r, C = 30.0, None
    rates = {}
    for ctl in ("mixed", "diffusion", "reaction"):
        m = cl.Classical(_p(material={"control": ctl, "k_r_A_s": 100.0}))
        n = 4 * math.pi * r ** 3 / (3 * m.Vm)
        C = 10 * m.c_inf
        gt = math.exp(2 * m.gam * m.Vm / (r * m.kt))
        rates[ctl] = float(m.rate(n, C))
        drive = cl.NA * (C - m.c_inf * gt) * 4 * math.pi * r * r
    assert math.isclose(rates["diffusion"], drive * m.D / r, rel_tol=1e-9)
    assert math.isclose(rates["reaction"], drive * m.k_r, rel_tol=1e-9)
    assert math.isclose(1 / rates["mixed"], 1 / rates["diffusion"] + 1 / rates["reaction"], rel_tol=1e-9)
    res = cl.Classical(_p()).run()
    w = res["nucleation_window"]
    assert w["t10"] <= w["t50"] <= w["t90"] and w["nuclei_mol"] > 0


def test_mass_is_conserved_with_syringe_and_temperature_program():
    for cond in ({"t_end": 600.0},
                 {"c_P0": 0.0, "t_end": 3000.0, "T_growth": 470.0, "t_switch_s": 100.0, "ramp_s": 50.0,
                  "injection": {"c_syringe_M": 0.3, "rate_L_h": 0.06, "V0_L": 0.01, "duration_s": 1200.0}}):
        res = cl.Classical(_p(conditions=cond)).run()
        assert res["series"][-1]["N_dots"] > 0
        for c in cl.conservation(res):
            assert abs(c["MA"]) < 1e-9
        if "injection" in cond:
            assert abs(res["series"][-1]["V_L"] - 0.01 - 0.06 / 3600 * 1200) < 1e-12
            assert abs(res["series"][-1]["T"] - 470.0) < 1.0


def _ripening(material, t_end, n_seed=2000, window=(0.03, 0.3)):
    """Mean radius against t while a fraction `window` of the seeded cohorts is left: past the transient of the
    uniform start, before the last few dots make the mean jumpy."""
    P = _p(conditions={"c_P0": 0.0, "t_end": t_end, "t_start": t_end * 1e-4, "n_out": 200, "c_M0": 1e-12},
           material={"J0_mode": "set", "J0_M_s": 0.0, **material},
           numerics={"merge_tol": 0.0, "max_cohorts": 100000})
    m = cl.Classical(P)
    rng = np.random.default_rng(1)
    r0 = rng.uniform(20.0, 40.0, n_seed)
    seed = [{"N": 1e-7, "n": 4 * math.pi * r ** 3 / (3 * m.Vm)} for r in r0]           # far more CdSe than C_inf
    res = m.run(seed=seed)
    ts, rbar = [], []
    for s in res["snapshots"]:
        if not window[0] * n_seed <= len(s["cohorts"]) <= window[1] * n_seed:
            continue
        w = np.array([x["N"] for x in s["cohorts"]])
        r = np.array([5.0 * x["d_nm"] for x in s["cohorts"]])                        # radius, A
        ts.append(s["t"])
        rbar.append((w * r).sum() / w.sum())
    assert len(ts) >= 5
    return m, np.array(ts), np.array(rbar)


def test_lsw_diffusion_limited_rate_constant():
    # the diffusion limit approaches LSW slowly from a uniform start (1.55, 1.25, 1.04 x K at 30, 10, 3 % left)
    m, t, r = _ripening({"k_r_A_s": 1e14, "C_inf_M": 1e-3, "gamma_J_m2": 0.2, "D_cm2_s": 1e-5}, 3.0,
                        n_seed=5000, window=(0.01, 0.03))
    K = 8 * m.gam * m.Vm ** 2 * cl.NA * m.c_inf * m.D / (9 * m.kt)                  # A^3/s
    slope = np.polyfit(t, r ** 3, 1)[0]
    assert 0.85 < slope / K < 1.15


def test_wagner_reaction_limited_rate_constant():
    m, t, r = _ripening({"k_r_A_s": 1e3, "C_inf_M": 1e-2, "gamma_J_m2": 0.2, "D_cm2_s": 10.0}, 1e4)
    K = 64 * m.gam * m.Vm ** 2 * cl.NA * m.c_inf * m.k_r / (81 * m.kt)              # A^2/s
    slope = np.polyfit(t, r ** 2, 1)[0]
    assert 0.85 < slope / K < 1.15
