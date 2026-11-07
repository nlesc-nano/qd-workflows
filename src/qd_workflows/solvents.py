"""Solvents of the growth runs: the static dielectric constant (the implicit solvent's eps) and the viscosity at the
run's temperature (Stokes-Einstein diffusion of the collision rates), so a run names its solvent instead of setting
eps and eta.

eta(T) is interpolated linearly in ln eta against 1/T through each solvent's points (mPa s, K), and beyond them
extended with the slope of the two end points (an Andrade law, ln eta = A + B / T, on that pair); past the hottest
point it is an extrapolation (the page says so). eps is taken at one temperature: for these nonpolar or weakly polar
liquids it changes by a few per cent up to 600 K. "basis" says how good the numbers are: measured, a published
correlation (DIPPR, VDI-PPDS, ChemSep, Viswanath & Natarajan 1989), or an estimate from an analogue; every source is
in the entry. Collected 2026-11-06; TOPO and paraffin oil are left out (no dielectric constant found).
"""
from __future__ import annotations

import math

import numpy as np

# name: formula, eps (at eps_T_K), eta [(T K, mPa s)], normal boiling and melting points (K), basis, sources
SOLVENTS: dict[str, dict] = {
    "1-octadecene": {
        "formula": "C18H36", "eps": 2.15, "eps_T_K": 293.0, "Tb_K": 588.0, "Tm_K": 290.8, "basis": "correlation",
        "eta": [(290.76, 4.469), (324.94, 2.256), (359.12, 1.298), (393.31, 0.8215), (427.49, 0.5595),
                (461.67, 0.4034), (495.85, 0.3043), (530.04, 0.2380)],
        "source": "eta: ChemSep v8.3 correlation (Cheméo), checked against Sigma-Aldrich 822112, 3.6 mm2/s at 40 °C; "
                  "eps: estimated from the CRC 1-alkenes (1-decene 2.136, 1-dodecene 2.152); Tb, Tm: KDB, NIST. "
                  "No measured eta above 40 °C; ODE polymerises slowly at 120-320 °C (Dhaene et al. 2019)"},
    "squalane": {
        "formula": "C30H62", "eps": 1.9106, "eps_T_K": 373.2, "Tb_K": 623.0, "Tm_K": 235.0, "basis": "measured",
        "eta": [(293.0, 36.1), (298.15, 28.33), (323.0, 10.7), (373.0, 3.10), (413.2, 1.68), (453.15, 1.062),
                (473.05, 0.853)],
        "source": "eta: Comuñas et al., J. Phys. Chem. Ref. Data 42, 033101 (2013); Schmidt, Trusler et al., J. Chem. "
                  "Eng. Data 60, 137 (2015); eps: CRC (373 K); Tb, Tm: Stenutz"},
    "n-hexadecane": {
        "formula": "C16H34", "eps": 2.046, "eps_T_K": 293.2, "Tb_K": 560.0, "Tm_K": 291.3, "basis": "measured",
        "eta": [(293.15, 3.486), (303.15, 2.772), (323.15, 1.840), (373.15, 0.8927), (423.15, 0.5351),
                (473.15, 0.3440), (523.15, 0.2360), (573.15, 0.1705)],
        "source": "eta: Klein, Fröba et al., J. Chem. Eng. Data 64 (2019), doi:10.1021/acs.jced.9b00525 (saturated "
                  "liquid); eps: CRC; Tb: CRC; Tm: NIST"},
    "n-octadecane": {
        "formula": "C18H38", "eps": 2.05, "eps_T_K": 303.0, "Tb_K": 589.3, "Tm_K": 301.3, "basis": "measured",
        "eta": [(303.35, 3.67), (318.15, 2.83), (373.15, 1.136), (473.15, 0.408), (573.15, 0.222)],
        "source": "eta: Thermochim. Acta 668, 159 (2018) to 318 K, DIPPR (Perry's 2-313) above, consistent with Baled "
                  "et al., J. Chem. Thermodyn. 72, 108 (2014); eps: estimated from the CRC alkanes; Tb, Tm: NIST"},
    "oleic acid": {
        "formula": "C18H34O2", "eps": 2.336, "eps_T_K": 293.2, "Tb_K": 633.0, "Tm_K": 289.45, "basis": "correlation",
        "eta": [(293.15, 39.1), (325.03, 12.85), (363.52, 5.07), (402.02, 2.50), (440.52, 1.45), (479.01, 0.952),
                (517.51, 0.683), (556.01, 0.525)],
        "source": "eta: Thermo Scientific spec at 293 K, ChemSep v8.3 correlation (Cheméo) above; eps: CRC; Tb: "
                  "ChemSep (decomposes); Tm: NIST"},
    "oleylamine": {
        "formula": "C18H37N", "eps": 2.7, "eps_T_K": 298.0, "Tb_K": 622.0, "Tm_K": 294.0, "basis": "estimate",
        "eta": [(298.15, 11.0), (373.15, 2.0), (473.15, 0.63)],
        "source": "eta: no data found; estimated as ODE (ChemSep) x the n-octylamine / n-octane ratio; eps: estimated "
                  "from the CRC octadecyl- and hexadecylamine; technical grades are 70-90 % pure"},
    "1-hexadecylamine": {
        "formula": "C16H35N", "eps": 2.71, "eps_T_K": 328.35, "Tb_K": 603.0, "Tm_K": 319.0, "basis": "estimate",
        "eta": [(323.15, 4.5), (373.15, 1.7), (473.15, 0.58)],
        "source": "eta: no data found; estimated as n-hexadecane x the n-octylamine / n-octane ratio; eps: CRC (328 K); "
                  "Tb, Tm: NIST"},
    "trioctylamine": {
        "formula": "C24H51N", "eps": 2.25, "eps_T_K": 293.0, "Tb_K": 639.2, "Tm_K": 238.6, "basis": "estimate",
        "eta": [(283.15, 14.25), (293.15, 9.96), (298.15, 8.39), (303.15, 7.11), (373.15, 2.1), (423.15, 1.2),
                (473.15, 0.80)],
        "source": "eta: Fang et al., J. Chem. Thermodyn. 68 (2014) to 303 K; above, the 298 K value scaled by "
                  "tri-n-hexylamine's curve (Viswanath & Natarajan); eps: estimated from the CRC trialkylamines; "
                  "Tb, Tm: NIST"},
    "trioctylphosphine": {
        "formula": "C24H51P", "eps": 2.4, "eps_T_K": 298.0, "Tb_K": None, "Tm_K": None, "basis": "estimate",
        "eta": [(298.15, 8.4), (373.15, 2.1), (473.15, 0.8)],
        "source": "no measured eta or eps found: eta as trioctylamine (same size, no H-bonds), eps by analogy with the "
                  "trialkylamines (2.3-2.6); boils 557-564 K at 6.7 kPa (Sigma-Aldrich), no normal boiling point"},
    "dioctyl ether": {
        "formula": "C16H34O", "eps": 2.7, "eps_T_K": 295.0, "Tb_K": 559.7, "Tm_K": 265.6, "basis": "estimate",
        "eta": [(298.15, 3.0), (373.15, 0.89), (473.15, 0.34)],
        "source": "eta: no data found; as n-hexadecane (same mass and chain); eps: between the CRC dipentyl (2.80) and "
                  "didecyl ethers (2.64); Tb, Tm: NIST"},
    "dibenzyl ether": {
        "formula": "C14H14O", "eps": 3.821, "eps_T_K": 293.2, "Tb_K": 571.2, "Tm_K": 276.8, "basis": "correlation",
        "eta": [(293.15, 5.33), (298.15, 4.64), (323.15, 2.60), (353.15, 1.57), (373.15, 1.20), (423.15, 0.74),
                (473.15, 0.53)],
        "source": "eta: Viswanath & Natarajan (1989) equation, valid 270-380 K, extrapolated above; 5.33 mPa s at "
                  "293 K measured (SJTU ai4fuel); eps: CRC; Tb: NIST; Tm: 276.8 K (CRC 3.6 °C)"},
    "diphenyl ether": {
        "formula": "C12H10O", "eps": 3.726, "eps_T_K": 283.2, "Tb_K": 531.5, "Tm_K": 300.0, "basis": "measured",
        "eta": [(303.15, 3.33), (323.15, 2.17), (353.15, 1.33), (373.15, 1.03), (423.15, 0.498), (473.15, 0.301),
                (523.15, 0.200)],
        "source": "eta: Cabaleiro et al., J. Chem. Thermodyn. 93, 86 (2016) to 373 K, DIPPR (Perry's 2-313) above (8 % "
                  "below the measurement at 373 K); eps: CRC (283 K, supercooled); Tb, Tm: NIST"},
    "1-dodecanethiol": {
        "formula": "C12H26S", "eps": 3.4, "eps_T_K": 293.0, "Tb_K": 547.7, "Tm_K": 265.7, "basis": "estimate",
        "eta": [(293.15, 3.24), (373.15, 1.2), (473.15, 0.56)],
        "source": "eta: Sigma-Aldrich 443816 at 293 K, above scaled by 1-octanethiol's DIPPR curve; eps: estimated "
                  "from the CRC 1-alkanethiols by molar volume; Tb: NIST"},
    "toluene": {
        "formula": "C7H8", "eps": 2.379, "eps_T_K": 296.35, "Tb_K": 383.8, "Tm_K": 178.1, "basis": "correlation",
        "eta": [(298.15, 0.554), (323.15, 0.426), (353.15, 0.321), (373.15, 0.268), (423.15, 0.184), (473.15, 0.136)],
        "source": "eta: DIPPR (Perry's 2-313) to Tb, VDI-PPDS saturated liquid above (pressurised); eps: CRC; Tb, Tm: NIST"},
    "n-hexane": {
        "formula": "C6H14", "eps": 1.8865, "eps_T_K": 293.2, "Tb_K": 341.9, "Tm_K": 178.0, "basis": "measured",
        "eta": [(293.15, 0.3127), (303.15, 0.2809), (323.15, 0.2340), (373.15, 0.1550), (423.15, 0.1026),
                (473.15, 0.0655)],
        "source": "eta: Klein et al., J. Chem. Eng. Data 64 (2019) (saturated liquid, pressurised above Tb); eps: CRC; "
                  "Tb, Tm: NIST"},
    "chloroform": {
        "formula": "CHCl3", "eps": 4.8069, "eps_T_K": 293.2, "Tb_K": 334.3, "Tm_K": 209.6, "basis": "correlation",
        "eta": [(298.15, 0.539), (323.15, 0.429), (333.15, 0.395), (373.15, 0.300)],
        "source": "eta: DIPPR (Perry's 2-313) to 333 K, VDI-PPDS saturated liquid at 373 K; eps: CRC; Tb, Tm: NIST"},
}


def viscosity(name: str, T: float) -> float:
    """eta (mPa s) at T: linear in ln eta against 1/T through the points, the end pair's slope beyond them."""
    pts = sorted(SOLVENTS[name]["eta"], key=lambda p: -p[0])      # 1/T ascending
    x = np.array([1.0 / t for t, _ in pts])
    y = np.log([e for _, e in pts])
    u = 1.0 / T
    if u < x[0]:
        return float(math.exp(y[0] + (u - x[0]) * (y[1] - y[0]) / (x[1] - x[0])))
    if u > x[-1]:
        return float(math.exp(y[-1] + (u - x[-1]) * (y[-1] - y[-2]) / (x[-1] - x[-2])))
    return float(math.exp(np.interp(u, x, y)))


def properties(name: str, T: float) -> tuple[float, float]:
    """(eps, eta in mPa s at T)."""
    if name not in SOLVENTS:
        raise ValueError(f"unknown solvent {name!r}: one of {', '.join(SOLVENTS)}, or 'custom'")
    return float(SOLVENTS[name]["eps"]), viscosity(name, T)


def table() -> list[dict]:
    """The page's list: name, eps, the viscosity points, the boiling point, the basis and the sources."""
    return [{"name": name, "formula": s["formula"], "eps": s["eps"], "eta": s["eta"], "Tb_K": s["Tb_K"],
             "Tm_K": s["Tm_K"], "basis": s["basis"], "source": s["source"]} for name, s in SOLVENTS.items()]
