"""
Reaction formulas on results.json (mirrored in JS by report.JS_MODEL for the live page).

    mu_i      = G_i^sol(T, eps) + kT ln a_i
    G_i^sol   = E_i^DFT + G_i^corr(T) + dG_solv,i(eps)             (gas-phase species: no dG_solv)
    a_i       = c_i / 1 M (solute) | c_liquid / 1 M (solvent) | p / (c0 R T) (gas)
    dG_r      = sum_i nu_i mu_i + nu_L dG_lig = dG_r^0 + kT ln Q
    dG_lig    per released ligand L (res["ligand"]): what the model leaves out (chain packing, conformations)
    dG(T, dmu_X) = dG_r^0 + nu_X dmu_X + sum_{i != X} nu_i kT ln a_i,   dmu_X = kT ln a_X
    equilibrium:  sum_solutes nu_i ln((c_i0 + nu_i xi) / 1 M) = -(dG_r^0 + sum_others nu_i kT ln a_i) / kT
    recipe (res["recipe"], res["equilibria"]): the reaction and the side equilibria (e.g. amine + HCl
        -> adduct) at equilibrium together, from what was put in the flask; see speciate()
"""
from __future__ import annotations

import math

import numpy as np

from .solvation import dg_solv

KB_EV = 8.617333262e-5
R_GAS = 8.314462618
C0 = 1000.0               # mol/m^3


def ln_activity(sp: dict, T, c: float | None = None) -> np.ndarray:
    """ln a_i on the T array (c in mol/L overrides the species' own for solutes)."""
    T = np.asarray(T, float)
    if sp["phase"] == "solvent":
        return np.log(float(sp["c_liquid"])) + 0 * T
    if sp["phase"] == "gas":
        return np.log(float(sp["p_bar"]) * 1e5 / (C0 * R_GAS * T))
    return np.log(float(sp["c"] if c is None else c)) + 0 * T


def cavity_model(cond: dict) -> str:
    """cond["cavity"]: beta0 | hille | none (older results: cond["nonel"] True/False)."""
    if "cavity" in cond:
        return cond["cavity"]
    return "beta0" if cond.get("nonel", True) else "none"


def solv(sp: dict, cond: dict) -> float:
    if sp["phase"] == "gas" or cond.get("gas"):
        return 0.0
    return dg_solv(sp, cond["eps"], cond.get("electrostatics", "cosmo"), cavity_model(cond))


def energy(sp: dict, cond: dict) -> float:
    """Electronic energy of the set cond["energy"] (given | pbe0_given | pbe0_mace | mace); nan if missing."""
    key = cond.get("energy", "given")
    if key == "given":
        return sp["energy_eV"]
    e = (sp.get("energies") or {}).get(key)
    return float("nan") if e is None else e


def g_std(sp: dict, cond: dict) -> np.ndarray:
    """E + G^corr(T) + dG_solv(eps): the 1 M standard free energy in solution (eV)."""
    return energy(sp, cond) + np.asarray(sp["thermo"]["G_corr"]) + solv(sp, cond)


def per_factor(res: dict) -> float:
    per = res.get("per")
    if per is None:
        return 1.0
    if isinstance(per, str):
        return abs(res["nu"][per])
    return float(per)


def reaction(res: dict, cond: dict, conc: dict | None = None) -> dict:
    """dG_r^0(T), dG_r(T) at the concentrations, and the decomposition (eV per reaction / per)."""
    conc = conc or {}
    T = np.asarray(res["T"], float)
    kT = KB_EV * T
    sps, nu = res["species"], res["nu"]
    g0 = sum(v * g_std(sps[n], cond) for n, v in nu.items())
    lnq = sum(v * ln_activity(sps[n], T, conc.get(n)) for n, v in nu.items())
    th = {k: sum(v * np.asarray(sps[n]["thermo"][k]) for n, v in nu.items())
          for k in ("zpe", "dU_vib", "H_trans", "H_rot", "S_trans", "S_rot", "S_vib", "S_el")}
    el_cond = {**cond, "cavity": "none"}
    solv_el = sum(v * solv(sps[n], el_cond) for n, v in nu.items())
    solv_all = sum(v * solv(sps[n], cond) for n, v in nu.items())
    k = per_factor(res)
    lig = res.get("ligand")
    lig_term = nu[lig["species"]] * cond.get("dG_lig", lig["dG_eV"]) if lig else 0.0
    g0 = g0 + lig_term
    parts = {
        "dE_DFT": sum(v * energy(sps[n], cond) for n, v in nu.items()) + 0 * T,
        "dZPE": th["zpe"], "dH_thermal": th["dU_vib"] + th["H_trans"] + th["H_rot"],
        "-TdS_trans": -T * th["S_trans"], "-TdS_rot": -T * th["S_rot"], "-TdS_vib": -T * th["S_vib"],
        "-TdS_el": -T * th["S_el"], "ddG_solv_el": solv_el + 0 * T, "ddG_solv_nonel": (solv_all - solv_el) + 0 * T,
        "ligand": lig_term + 0 * T,
    }
    return {"T": T, "dG0": g0 / k, "dG": (g0 + kT * lnq) / k, "kTlnQ": kT * lnq / k,
            "dH0": (parts["dE_DFT"] + parts["dZPE"] + parts["dH_thermal"] + solv_all + lig_term) / k,
            "dS0": (th["S_trans"] + th["S_rot"] + th["S_vib"] + th["S_el"]) / k,
            "parts": {p: v / k for p, v in parts.items()}, "per": k}


def scan_map(res: dict, cond: dict, species: str, log10_c, conc: dict | None = None) -> dict:
    """dG(T, c_X) on the T grid x log10 c_X grid; the dG = 0 line log10 c_X*(T)."""
    T = np.asarray(res["T"], float)
    kT = KB_EV * T
    r = reaction(res, cond, conc)
    nu_x = res["nu"][species]
    k = r["per"]
    lnx_now = ln_activity(res["species"][species], T, (conc or {}).get(species))
    base = r["dG"] - nu_x * kT * lnx_now / k                  # dG without the X term
    lc = np.asarray(log10_c, float)
    grid = base[:, None] + nu_x * kT[:, None] * np.log(10.0) * lc[None, :] / k
    star = -base * k / (nu_x * kT * np.log(10.0))
    return {"T": T, "log10_c": lc, "dG": grid, "log10_c_star": star}


def equilibrium(res: dict, cond: dict, c0: dict | None = None) -> dict:
    """Extent of reaction xi(T) from initial concentrations (reactants: their c, products: 0)."""
    T = np.asarray(res["T"], float)
    kT = KB_EV * T
    sps, nu = res["species"], res["nu"]
    k = per_factor(res)
    dg0 = reaction(res, cond)["dG0"] * k
    sol = [n for n in nu if sps[n]["phase"] == "solute"]
    other = [n for n in nu if n not in sol]
    init = {n: (c0 or {}).get(n, sps[n].get("c0", sps[n]["c"] if nu[n] < 0 else 0.0)) for n in sol}
    lo = max([-init[n] / nu[n] for n in sol if nu[n] > 0], default=0.0)
    hi = min([init[n] / -nu[n] for n in sol if nu[n] < 0] + [math.inf])
    if not math.isfinite(hi) or hi <= lo:
        return {"T": T, "xi": np.full_like(T, np.nan), "conversion": np.full_like(T, np.nan), "xi_max": hi}
    xi = np.empty_like(T)
    for j, t in enumerate(T):
        target = -(dg0[j] + sum(nu[n] * kT[j] * ln_activity(sps[n], [t])[0] for n in other)) / kT[j]
        a, b = lo, hi
        for _ in range(200):
            m = 0.5 * (a + b)
            f = sum(nu[n] * math.log(max(init[n] + nu[n] * m, 1e-300)) for n in sol) - target
            if f < 0:
                a = m
            else:
                b = m
        xi[j] = 0.5 * (a + b)
    return {"T": T, "xi": xi, "conversion": (xi - lo) / (hi - lo), "xi_max": hi, "init": init}


# --------------------------------------------------------------------------
# Recipe: the reaction and its side equilibria from the amounts put in the flask
# --------------------------------------------------------------------------
#
#   c_i = c_i0 + sum_r nu_ri xi_r                  (dissolved species: mass balance)
#   a_i fixed                                      (gas at its partial pressure, pure liquid)
#   g_r(xi) = dG_r^0 / kT + sum_i nu_ri ln a_i = 0 for every reaction r
# G is convex in the extents, so g_r rises with xi_r: the main extent by bisection, each side
# extent (for a given main extent) by bisection, swept until they settle. Fraction of dots
# converted: theta = xi_0 / c_dot.

def dg0_on_T(res: dict, cond: dict, nu: dict) -> np.ndarray:
    """sum_i nu_i (E + G^corr + dG_solv) on the T grid (eV, per reaction as written)."""
    return sum(v * g_std(res["species"][n], cond) for n, v in nu.items()) + 0 * np.asarray(res["T"], float)


def _bounds(nu, amounts, balanced):
    lo, hi = -math.inf, math.inf
    for n, v in nu.items():
        if n not in balanced:
            continue
        if v > 0:
            lo = max(lo, -amounts[n] / v)
        else:
            hi = min(hi, amounts[n] / -v)
    return lo, hi


def _bisect(f, lo, hi, n=80):
    for _ in range(n):
        m = 0.5 * (lo + hi)
        if f(m) < 0:
            lo = m
        else:
            hi = m
    return 0.5 * (lo + hi)


def speciate(res: dict, cond: dict, c_dot: float | None = None, equivalents: float | None = None,
             c_additive: float | None = None, dg_side: list | None = None, ti: int | None = None) -> dict:
    """Equilibrium of the recipe on the T grid: theta (fraction of dots converted) and every concentration (M)."""
    rc = res["recipe"]
    c_dot = rc["c_dot"] if c_dot is None else c_dot
    eqv = rc["equivalents"] if equivalents is None else equivalents
    c_add = rc["c_additive"] if c_additive is None else c_additive
    T = np.asarray(res["T"], float)
    sps = res["species"]
    reactions = [res["nu"]] + [e["nu"] for e in res.get("equilibria") or []]
    lig = res.get("ligand")
    lig_term = res["nu"][lig["species"]] * cond.get("dG_lig", lig["dG_eV"]) if lig else 0.0
    dg0 = [dg0_on_T(res, cond, reactions[0]) + lig_term]
    dg0 += [dg0_on_T(res, cond, nu) + ((dg_side or [])[k] if k < len(dg_side or []) else 0.0)
            for k, nu in enumerate(reactions[1:])]
    names = list(dict.fromkeys(n for nu in reactions for n in nu))
    balanced = [n for n in names if sps[n]["phase"] == "solute"]
    n0 = {n: 0.0 for n in balanced}
    n0[rc["dot"]] = c_dot
    n0[rc["reagent"]] = n0.get(rc["reagent"], 0.0) + eqv * c_dot
    if rc.get("additive") in n0:
        n0[rc["additive"]] += c_add
    out = {"T": T, "theta": np.empty_like(T), "c": {n: np.empty_like(T) for n in balanced}}
    for j, t in enumerate(T):
        if ti is not None and j != ti:
            out["theta"][j] = np.nan
            for n in balanced:
                out["c"][n][j] = np.nan
            continue
        kT = KB_EV * t
        fixed = [sum(v * ln_activity(sps[n], [t])[0] for n, v in nu.items() if n not in balanced) for nu in reactions]

        def amounts(xi):
            a = dict(n0)
            for nu, x in zip(reactions, xi):
                for n, v in nu.items():
                    if n in a:
                        a[n] += v * x
            return a

        def g(r, xi):
            a = amounts(xi)
            return dg0[r][j] / kT + fixed[r] + sum(v * math.log(max(a[n], 1e-300)) for n, v in reactions[r].items() if n in a)

        def sides(x0):
            xi = [x0] + [0.0] * (len(reactions) - 1)
            for _ in range(30 if len(reactions) > 2 else 1):
                old = list(xi)
                for r in range(1, len(reactions)):
                    rest = list(xi)
                    rest[r] = 0.0
                    lo, hi = _bounds(reactions[r], amounts(rest), balanced)
                    if not (math.isfinite(lo) and math.isfinite(hi)) or hi <= lo:
                        continue
                    xi[r] = _bisect(lambda x: g(r, xi[:r] + [x] + xi[r + 1:]), lo, hi)
                if max(abs(p - q) for p, q in zip(xi, old)) <= 1e-12 * max(c_dot, 1e-30):
                    break
            return xi

        lo, hi = _bounds(reactions[0], n0, balanced)
        if not (math.isfinite(lo) and math.isfinite(hi)) or hi <= lo or not np.isfinite(dg0[0][j]):
            out["theta"][j] = np.nan
            for n in balanced:
                out["c"][n][j] = np.nan
            continue
        x0 = _bisect(lambda x: g(0, sides(x)), lo, hi)
        a = amounts(sides(x0))
        out["theta"][j] = x0 / c_dot
        for n in balanced:
            out["c"][n][j] = max(a[n], 0.0)
    return out
