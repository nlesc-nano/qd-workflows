"""
Thermal correction G^corr(T) = G - E of one species: ideal gas at the 1 M standard state.

    H_trans = 5/2 kT,   S_trans = k [ln((2 pi m kT / h^2)^(3/2) / (N_A c0)) + 5/2]
    H_rot   = 3/2 kT,   S_rot   = k [ln(sqrt(pi)/sigma sqrt(T^3 / (Th_A Th_B Th_C))) + 3/2]   (linear: kT, k [ln(T/(sigma Th)) + 1])
    ZPE, U_vib, S_vib   harmonic U, quasi-RRHO S (Orchestr.AI qd.steps.hessian.thermo, the library's treatment)
    S_el    = k ln(2S + 1)
    G^corr  = ZPE + (U_vib - ZPE) + H_trans + H_rot - T (S_trans + S_rot + S_vib + S_el)

The same numbers as Orchestr.AI qd.references.ideal_gas_g (tested), split into the terms the
reaction decomposition reports.
"""
from __future__ import annotations

import numpy as np

KB_EV = 8.617333262e-5
H_SI, KB_SI, NA, AMU = 6.62607015e-34, 1.380649e-23, 6.02214076e23, 1.66053906660e-27
C0 = 1000.0                     # mol / m^3 (1 M)


def geometry(positions) -> str:
    p = np.asarray(positions, float)
    if len(p) == 1:
        return "monatomic"
    sv = np.linalg.svd(p - p.mean(axis=0), compute_uv=False)
    return "linear" if sv[1] < 1e-3 * sv[0] else "nonlinear"


def moments_si(symbols, positions) -> np.ndarray:
    """Principal moments of inertia (kg m^2), ascending."""
    from ase import Atoms
    return np.asarray(Atoms(list(symbols), positions=np.asarray(positions, float)).get_moments_of_inertia()) * AMU * 1e-20


def components(symbols, positions, masses, freqs_cm1, temperatures, sigma: int = 1, mult: int = 1) -> dict:
    """Every term of G^corr (eV, and eV/K for entropies) on `temperatures`; freqs are the 3N-6 (3N-5) vibrations."""
    from orchestr_ai.qd.steps.hessian import thermo
    T = np.asarray(temperatures, float)
    kT = KB_EV * T
    m = float(np.sum(masses)) * AMU
    geo = geometry(positions)
    s_trans = KB_EV * (np.log((2 * np.pi * m * KB_SI * T / H_SI ** 2) ** 1.5 / (NA * C0)) + 2.5)
    h_trans = 2.5 * kT
    if geo == "monatomic":
        h_rot, s_rot = 0 * T, 0 * T
    else:
        mom = moments_si(symbols, positions)
        if geo == "linear":
            theta = H_SI ** 2 / (8 * np.pi ** 2 * mom.max() * KB_SI)
            h_rot, s_rot = kT, KB_EV * (np.log(T / (sigma * theta)) + 1.0)
        else:
            theta = H_SI ** 2 / (8 * np.pi ** 2 * mom * KB_SI)
            h_rot = 1.5 * kT
            s_rot = KB_EV * (np.log(np.sqrt(np.pi) / sigma * np.sqrt(T ** 3 / np.prod(theta))) + 1.5)
    freqs = np.asarray(freqs_cm1, float)
    if freqs.size:
        th = thermo(freqs, list(T))
        zpe = th["zpe_eV"] + 0 * T
        u_vib = np.asarray(th["U_vib_eV"])
        s_vib = np.asarray(th["S_vib_meV_K"]) / 1000.0
    else:
        zpe, u_vib, s_vib = 0 * T, 0 * T, 0 * T
    s_el = KB_EV * np.log(mult) + 0 * T
    h = zpe + (u_vib - zpe) + h_trans + h_rot
    s = s_trans + s_rot + s_vib + s_el
    return {
        "T": T.tolist(), "geometry": geo, "sigma": int(sigma),
        "zpe": zpe.tolist(), "dU_vib": (u_vib - zpe).tolist(), "H_trans": h_trans.tolist(), "H_rot": h_rot.tolist(),
        "S_trans": s_trans.tolist(), "S_rot": s_rot.tolist(), "S_vib": s_vib.tolist(), "S_el": s_el.tolist(),
        "H_corr": h.tolist(), "S": s.tolist(), "G_corr": (h - T * s).tolist(),
    }
