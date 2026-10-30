"""
Nucleation and growth of a dot population after precursor injection: a continuum model coupled to the
precursor pool and to the nucleation grid's grand potential. The facet energetics and kinetics are placeholders
(PLACEHOLDER, `params`) until the calibration fits of the shape library replace them.

State: Se precursor P, monomer supply C_M (MA in solution: free and in the MA(MX_q)_n complexes), MX_q left in
solution C_X (free and complexed), and cohorts of dots, each born at t_b with a number concentration N_c (mol of
dots per L) and a size m (CdSe units); its shape is the Wulff construction at the current gamma_f(C) (h_f from it).

    precursor   Q = k_conv(T) P C_X  (MX_q + P=A -> MA + P X_q, grid_thermo)
    pool        mu_MA, mu_MX, free MX_q from (C_M, C_X)  (grid_thermo.pool: complexation equilibrium)
    nucleation  Omega(m,k) = G(m,k) - m mu_MA - k mu_MX  (thermo of the grid),  Omega(m) = -kT ln sum_k e^-Omega(m,k)/kT,
                barrier dOmega* = max_{m <= m0} Omega(m) - Omega(1),
                J = [sum_{m<m0} 1/(f_m c_eq(m))]^-1, c_eq(m) = c0 e^(-Omega(m)/kT), f_m Smoluchowski attachment
                (steady-state Becker-Doring flux, nucleation_flux), nuclei of m0
                CdSe and k*(m0) MX_q, Wulff-shaped at the current gamma_f(C)
    coverage    theta_f = K_f c / (1 + K_f c),  K_f = exp(-dG_bind,f / kT),  c = free MX_q (M)
    coverage    K_f c = theta_f/(1 - theta_f) e^(w_f theta_f/kT),  K_f = e^(-dG_bind,f/kT),
                dG_bind,f = -eps0_f + dG_th(T) - dG_solv(MX_q)  (Fowler-Guggenheim; fitted on the builder dots)
    surface     gamma_f(C) = gamma_f^0 + a_solv + Gamma_f^max [kT ln(1 - theta_f) - w_f theta_f^2 / 2]
    sites       a family with `sites` (site classes s of the half-swap dot, e.g. facet centre and the edge classes)
                instead: theta_s = K_s c / (1 + K_s c), K_s = e^(-dG_bind,s/kT) (release umbrellas), theta_f the
                class-weighted mean; the rate constant of the facet from its active classes (attachment and direct
                exchange umbrellas, 1D transition-state flux through the barrier of each profile),
                  k_f = sum_s rho_s [A_dir v_C e^(-dW_dir,s/kT) + K_s c A_disp v_CX e^(-dW_disp,s/kT)] / (1 + K_s c),
                v = (kT / 2 pi mass)^(1/2) (the complex; the complex against MX_q for the exchange), A the umbrellas'
                capture areas, rho_s per area (centre) or per length of the facet's perimeter (edges: rho_s = n_s P_f
                / A_f, so small facets have more), capped at the collision rate; MX_q blocks, the complex never stays.
                gamma_f = gamma_eff + a_solv - Gamma_eff kT ln c when gamma_eff is given (from the layer additions on
                the dots: Delta mu_f = 2 Omega_m gamma_f / h_f, no slab), else the classes' Langmuir sum on gamma0
    solvation   implicit only (grid_solvation): every solution species, monomer or cluster, and every facet
                (a_solv per area). EXPLICIT SOLVENT AND LIGANDS (oleic acid, amines, phosphines) ARE NOT INCLUDED:
                they are the next step and change the monomers, every (m, k) and the facets.
    hand-over   m0 = "auto": the size where the grid's two-step chemical potential meets the continuum (_handover)
    dot mu      mu_dot = g_bulk(T) + (2/3) Omega_m sum_f gamma_f A_f / V   (Gibbs-Thomson of the Wulff dot: one
                chemical potential per dot, = g_bulk + 2 gamma_f Omega_m / h_f on every facet of the Wulff shape)
                g_bulk(T) = E_bulk + F_vib(T) of the bulk phonons (material.g_bulk_T, free energy like the grid's
                G), or the Einstein estimate E_bulk + 3 n_at [hw_E/2 + kT ln(1 - e^-hw_E/kT)] without that table
                C_eq(m) = C_M e^((mu_dot - mu_MA)/kT)  (the monomer level in equilibrium with the dot, at fixed mu_MX)
                with grid nucleation, a dot within the grid's sizes takes the grid's step instead,
                mu_dot(m) = [G~(m) - G~(m - 2)] / 2, G~(m) = -kT ln sum_k e^-(G(m,k) - k mu_MX)/kT  (grid_mu)
    bridge      gamma_f^0 x s so that g_bulk(T) + 2 gamma_bar Omega_m / r0 of the m0 Wulff dot equals the grid step
                mu_MA + Omega(m0) - Omega(m0 - 1) at the pool of half conversion  (match_m0)
    size law    dm/dt = A n_A (C_M - C_eq) / (A / (4 pi D r) + 1 / k_eff)   (diffusion to the equivalent sphere of
                radius r = (3 m Omega_m / 4 pi)^(1/3) in series with the surface reaction, area A; n_A converts M to
                molecules/A^3), k_eff = k_r (1 - theta_bar) for growth, k_r for dissolution; k_r = None: the
                diffusion limit dm/dt = 4 pi D r n_A (C_M - C_eq). Limits: LSW ripening (diffusion: <r>^3 ~ t;
                reaction: <r>^2 ~ t) and size focusing at high supersaturation (Talapin-Peng)
    dot         P(h) = {r: n_i.r <= h_f(i)},  h = Wulff(gamma(C), m Omega_m),  m = V(h)/Omega_m,
                k = sum_f [theta_f Gamma_f^max + x_f / a_f] A_f(h)   (Z-type from coverage + X-type by electron counting)
    balance     dP/dt = -Q,  dC_M/dt = Q - m0 J - sum_c N_c dm_c/dt,  dC_X/dt = -Q - k0 J - sum_c N_c dk_c/dt

A dot at time t goes to the builder as a recipe (builder_recipe): facet_distances = h_f, target m, and
post_treatment.z_type_target with target_total = k and per_family weights theta_f Gamma_f^max A_f.

    python -m qd_workflows.growth params --out params.yaml                 (the placeholder parameters, annotated)
    python -m qd_workflows.growth run --thermo thermo_merged.json --params params.yaml --out growth.json
"""
from __future__ import annotations

import argparse
import copy
import json
import math
import sys
from pathlib import Path

import numpy as np
import yaml
from scipy.optimize import brentq

from .grid_solvation import SOLVATION, per_area
from .grid_thermo import KB_EV, k_conv_at, pool, thermo_at

A_ZB = 6.202964                                 # A: MACE-MH-1 zinc-blende CdSe (shapes fit)
SQ3 = math.sqrt(3.0)

# Every number below is a PLACEHOLDER (with where it came from) until the shape-library calibration replaces it.
PARAMS = {
    "conditions": {"T": 500.0, "c_P0": 0.1, "c_MX0": 0.2, "k_conv_500K": 1.0, "Ea_conv": 0.0, "t_end": 600.0,
                   "t_start": 1e-3, "n_out": 60,
                   "gamma_MA": 1.0, "gamma_MX": 1.0},     # activity coefficients (calibration, 1 = ideal)
    "material": {
        "Omega_m": A_ZB ** 3 / 4,               # A^3 per CdSe (4 formula units per cubic cell)
        "g_bulk": -5.4903,                      # eV per CdSe: MACE bulk zinc-blende energy at 0 K (shapes fit)
        "omega_E_cm1": 150.0,                   # Einstein phonon of bulk CdSe (Debye ~210 K) for g_bulk(T) (PLACEHOLDER)
        "atoms_per_unit": 2,
        "m0": "auto",                           # nucleus size: "auto" = where grid and continuum agree (_handover)
        "m0_min": 10,                           # smallest nucleus the automatic hand-over may pick
        "m_dissolve": 4,                        # a dot below this many CdSe units dissolves back into the pool
        "dt_rel_max": 0.01,                     # a step is at most 1 % of t (CdSe batch: within 7 % of 0.3 %)
        "dmu_step_kT": 0.2,                     # ... and moves the pool's mu_MA and mu_MX by at most this many kT (a
                                                # collapse of the free MX_q, e.g. carboxylate complexes, in many steps)
        "birth_interval": 0.05,                 # nuclei become a new cohort every 5 % of t ...
        "birth_share": 0.05,                    # ... or once they reach 5 % of the dots present
        "merge_tol": 0.01,                      # cohorts past 4 m_top within 1 % in m merge (N, m, k kept) ...
        "max_cohorts": 64,                      # ... and the closest such ones while there are more than this
        "N_min": 1e-15,                         # M: fewer nuclei than this in an output interval make no cohort
        "match_m0": False,                      # scale gamma^0 so the m0 dot's continuum chemical potential equals the
                                                # grid's free-energy step at m0 (the bridge used with placeholder gamma)
    },
    "nucleation": {"omega": "grid",             # grid: the cluster grid's Omega(m); continuum: classical nucleation
                                                # on the continuum's own free energy (the check of the kinetics)
                   "flux": "becker_doring",     # steady-state flux over the grid's Omega(m) (nucleation_flux)
                   "attach_barrier_eV": 0.0,    # incorporation barrier of a monomer on a cluster: 0 = diffusion limit
                   "J0": 1e-6},                 # M/s, used only with flux: constant (J = J0 e^(-barrier/kT))
    "kinetics": {"mu_small": "grid",            # grid nucleation: dots within the grid's sizes take the grid's own
                                                # chemical potential (grid_mu); "continuum": Gibbs-Thomson at every size
                 "transport": "stokes",         # stokes: contact radii and Stokes-Einstein D from the grid's structures
                                                # (transport.add_geometry), the dot's polyhedron capacitance; fixed: D
                                                # below with bulk-volume spheres (also when the thermo has no geometry)
                 "viscosity_mPa_s": 0.4,        # solvent viscosity at T (PLACEHOLDER: an n-alkane like octadecene near
                                                # 500 K, about 0.4 mPa s; to be set from data for the solvent used)
                 "stokes_c": 6.0,               # 6: stick boundary, 4: slip
                 "D": 1e-5,                     # cm^2/s, monomer diffusion coefficient, transport "fixed" (PLACEHOLDER)
                 "k_r": 1e4},                   # A/s, surface reaction rate constant of the size law (PLACEHOLDER);
                                                # None: the diffusion limit (D/r ~ 1e10 A/s), which turns every run into
                                                # one extreme LaMer burst (a few dots of ~100 nm)
    # Zinc-blende CdSe facets from the fit of 108 relaxed QD_builder dots (MACE-MH-1, shape ladder, 2026-10-25):
    #   gamma0     bare facet energy, the empty-shell limit of the passivated-dot fit (+- ~2 meV/A^2)
    #   eps0       binding energy of one CdCl2 on the empty facet (stripping intercept); the binding free energy at T
    #              is dG_bind = -eps0 + dG_th(T), dG_th from the grid's Hessians (binding_thermal)
    #   Gamma_max  CdCl2 per nm^2 at saturation on the builder's Wulff dots ({100} tight, the others wide)
    #   w          lateral interaction (Fowler-Guggenheim), w = 2 delta Gamma_max with delta = 12.1 eV A^2 of the fit
    #   tau        step energy, the dot fit's single edge tension (14 +- 10 meV/A; kept, not used by the size law)
    # site_area: area per surface site from the CIF (shapes); k0 is left from the facet law (not used).
    "families": {
        "100":  {"gamma0": 11.3, "eps0": 1.38, "Gamma_max": 3.05, "w": 0.738, "tau": 14.0, "k0": 1e4, "x_X": 0.0,
                 "site_area": A_ZB ** 2 / 2},
        "111A": {"gamma0": 25.3, "eps0": 1.06, "Gamma_max": 1.26, "w": 0.305, "tau": 14.0, "k0": 1e4, "x_X": 0.0,
                 "site_area": SQ3 * A_ZB ** 2 / 4},
        "111B": {"gamma0": 21.9, "eps0": 1.08, "Gamma_max": 0.66, "w": 0.160, "tau": 14.0, "k0": 1e4, "x_X": 0.0,
                 "site_area": SQ3 * A_ZB ** 2 / 4},
        "110":  {"gamma0": 21.5, "eps0": 1.28, "Gamma_max": 1.21, "w": 0.293, "tau": 14.0, "k0": 1e4, "x_X": 0.0,
                 "site_area": A_ZB ** 2 / math.sqrt(2)},
    },
    "provenance": {"families": "fit of 108 QD_builder dots (MACE-MH-1)"},
    # implicit solvation of every solution species (grid_solvation): monomers from the campaign's COSMO library,
    # clusters from their geometric surface; eps = 1 is the gas phase. Explicit solvent and ligands come later.
    "solvation": copy.deepcopy(SOLVATION),
    "units": {"gamma0": "meV/A^2 (bare facet)", "tau": "meV/A", "eps0": "eV per Z-type MX_q on the empty facet (> 0: bound)",
              "dG_bind": "eV per Z-type MX_q at 1 M (if given instead of eps0)", "Gamma_max": "Z-type MX_q per nm^2",
              "w": "eV (lateral interaction at full coverage)",
              "k0": "A/s", "k_r": "A/s", "D": "cm^2/s", "x_X": "X-type MX_q per surface site (electron counting)", "site_area": "A^2",
              "J0": "M/s", "attach_barrier_eV": "eV", "T": "K", "c_P0": "M", "c_MX0": "M", "k_conv_500K": "1/(M s)",
              "Ea_conv": "eV", "t_end": "s", "gamma_MA": "1 (activity coefficient)", "gamma_MX": "1 (activity coefficient)",
              "c_MX_sat": "M (solubility of the free MX_q; optional)",
              "g_MX_solid_T": "eV per MX_q formula unit of the solid at the gas-phase level of the grid's MX_q reference "
                              "(optional; sets c_MX_sat when that is not given; clamped outside its T range)",
              "g_MX_solid_shift_eV": "eV (added to g_MX_solid_T)"},
}


def binding_thermal(thermo: dict, t: float, mx: str | None = None) -> float:
    """Thermal part of binding one MX_q on a cluster at T, from the grid's own Hessians: the mean change of G - E per
    added MX_q between consecutive k of the best states of the upper half of the rows (k >= 2), minus G - E of the
    gas-phase monomer at 1 M. Adds to the binding energy: dG_bind = -eps0 + dG_th (0 when the thermo has no E)."""
    from .grid_thermo import _at, unit_names
    mx = mx or unit_names(thermo["nodes"])[0]
    ref = thermo["references"].get(mx, {})
    if "E_eV" not in ref:
        return 0.0
    temps = thermo["temperatures"]
    best = {(n["m"], n["k"]): n["states"][0] for n in thermo["nodes"]}
    m_lo = (max(m for m, _ in best) + 1) // 2 + 1
    d = [(_at(t, temps, best[(m, k + 1)]["G_eV"]) - best[(m, k + 1)]["E_eV"]) - (_at(t, temps, best[(m, k)]["G_eV"]) - best[(m, k)]["E_eV"])
         for (m, k) in best if m >= m_lo and k >= 2 and (m, k + 1) in best]
    if not d:
        return 0.0
    return float(np.mean(d) - (_at(t, temps, ref["G_eV"]) - ref["E_eV"]))


def gamma_max(p: dict) -> float:
    """Z-type MX_q per A^2 at saturation: Gamma_max (per nm^2), or sites_per_unit per surface site."""
    return float(p["Gamma_max"]) * 1e-2 if "Gamma_max" in p else float(p["sites_per_unit"]) / float(p["site_area"])


def facet_state(p: dict, c: float, kt: float) -> dict:
    """Coverage and surface energy of a facet in contact with free MX_q at c (M), 1 M standard state:
    K = e^(-dG_bind/kT),  K c = theta/(1 - theta) e^(w theta/kT)  (Fowler-Guggenheim; w = 0: Langmuir),
    gamma(c) = gamma0 + Gamma_max [kT ln(1 - theta) - w theta^2 / 2]  (Gibbs adsorption integrated from c = 0)."""
    gmax, w = gamma_max(p), float(p.get("w", 0.0))
    if "sites" in p:
        return sites_state(p, c, kt, gmax)
    kc = math.exp(min(-float(p["dG_bind"]) / kt, 700.0)) * c
    theta = theta_fg(kc, w / kt)
    return {"theta": theta, "gamma": float(p["gamma0"]) * 1e-3 + gmax * (kt * math.log1p(-theta) - 0.5 * w * theta ** 2),
            "Gamma_max": gmax}


def class_weights(p: dict) -> dict:
    """Share of the family's MX_q sites in each site class (from the class counts n on the reference dot)."""
    n = {s: float(q["n"]) for s, q in p["sites"].items()}
    tot = sum(n.values())
    return {s: x / tot for s, x in n.items()}


def sites_state(p: dict, c: float, kt: float, gmax: float) -> dict:
    """Coverage and surface energy of a family with site classes (facet_state): each class its own Langmuir
    coverage theta_s = K_s c / (1 + K_s c), the family's theta their class-weighted mean; gamma = gamma_eff - Gamma_eff
    kT ln c when gamma_eff is given (meV/A^2 at c = 1 M, Gamma_eff MX_q per nm^2); else, while the family still
    carries the shape-library fit (dG_bind, w), that fit's Fowler-Guggenheim gamma; else gamma0 + Gamma_max kT
    sum_s w_s ln(1 - theta_s)."""
    wts = class_weights(p)
    th = {}
    for s, q in p["sites"].items():
        kc = math.exp(min(-float(q["dG_bind"]) / kt, 700.0)) * c
        th[s] = kc / (1.0 + kc)
    theta = sum(wts[s] * th[s] for s in th)
    if "gamma_eff" in p:
        gam = float(p["gamma_eff"]) * 1e-3 - float(p.get("Gamma_eff", 0.0)) * 1e-2 * kt * math.log(max(c, 1e-300))
    elif "dG_bind" in p:
        gam = facet_state({k: v for k, v in p.items() if k != "sites"}, c, kt)["gamma"]
    else:
        gam = float(p["gamma0"]) * 1e-3 + gmax * kt * sum(wts[s] * math.log(max(1.0 - th[s], 1e-300)) for s in th)
    return {"theta": theta, "gamma": gam, "Gamma_max": gmax, "theta_s": th}


def theta_fg(kc: float, a: float) -> float:
    """Fowler-Guggenheim coverage: K c = theta/(1 - theta) e^(a theta), a = w/kT (a = 0: Langmuir)."""
    if kc <= 0:
        return 0.0
    target = math.log(kc)
    lo, hi = 1e-15, 1.0 - 1e-15
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if math.log(mid) - math.log1p(-mid) + a * mid < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# capacitances of the dots' rounded shapes (Growth.dot_capacitance), shared by all runs of a process; cleared when full
SHAPE_CAPS: dict = {}
SHAPE_CAPS_MAX = 20000


FA_FAMILIES = ("100", "111A", "111B")           # the facet-area law's families (the computed sitekin dots')

ZB_NORMALS = {                                  # outward normals of the zinc-blende families (Td: {111} and {-1-1-1} apart)
    "100": [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)],
    "111A": [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)],
    "111B": [(-1, -1, -1), (-1, 1, 1), (1, -1, 1), (1, 1, -1)],
    "110": [(a, b, 0) for a in (1, -1) for b in (1, -1)] + [(a, 0, b) for a in (1, -1) for b in (1, -1)]
           + [(0, a, b) for a in (1, -1) for b in (1, -1)],
}


# ---------------- geometry of a dot ----------------
class Shape:
    """Polyhedron of facet distances h_f: areas per family, volume, and each family's support distance (where its
    plane would touch the dot when the family is not exposed). The combinatorics (each vertex's three planes, each
    face's vertex cycle) are cached: while every vertex of a cached topology still satisfies every half-space, the
    vertices follow from 3x3 solves and the geometry is exact; otherwise the half-space intersection is redone."""

    def __init__(self, normals: dict):
        self.fams = list(normals)
        self.n = {f: np.array([np.asarray(v, float) / np.linalg.norm(v) for v in vs]) for f, vs in normals.items()}
        self.N = np.vstack([self.n[f] for f in self.fams])
        self.fam_of = np.concatenate([[i] * len(self.n[f]) for i, f in enumerate(self.fams)])
        self.topos = []

    def _topology(self, hv):
        from scipy.spatial import HalfspaceIntersection
        hs = np.hstack([self.N, -hv[:, None]])
        v = np.unique(np.round(HalfspaceIntersection(hs, np.zeros(3)).intersections, 8), axis=0)
        tight = np.abs(v @ self.N.T - hv[None, :]) < 1e-6 * max(1.0, hv.max())
        tri = []
        for row in tight:
            pick = []
            for j in np.flatnonzero(row):
                if np.linalg.matrix_rank(self.N[pick + [j]], tol=1e-8) == len(pick) + 1:
                    pick.append(j)
                if len(pick) == 3:
                    break
            tri.append(pick)
        a, b, face = [], [], []
        for j in range(len(self.N)):
            on = np.flatnonzero(tight[:, j])
            if len(on) < 3:
                continue
            p = v[on]
            c = p.mean(0)
            u = (p[0] - c) / np.linalg.norm(p[0] - c)
            w = np.cross(self.N[j], u)
            on = on[np.argsort(np.arctan2((p - c) @ w, (p - c) @ u))]
            a += list(on)
            b += list(np.roll(on, -1))
            face += [j] * len(on)
        return {"tri": np.array(tri), "a": np.array(a, int), "b": np.array(b, int), "face": np.array(face, int)}

    def _eval(self, t, hv, tol=1e-9):
        v = np.linalg.solve(self.N[t["tri"]], hv[t["tri"]][..., None])[..., 0]
        if (v @ self.N.T > hv[None, :] + tol * max(1.0, hv.max())).any():
            return None
        p, q = v[t["a"]], v[t["b"]]
        cr = np.stack([p[:, 1] * q[:, 2] - p[:, 2] * q[:, 1], p[:, 2] * q[:, 0] - p[:, 0] * q[:, 2],
                       p[:, 0] * q[:, 1] - p[:, 1] * q[:, 0]], axis=1)
        proj = (cr * self.N[t["face"]]).sum(1)
        a_plane = 0.5 * np.abs(np.bincount(t["face"], weights=proj, minlength=len(self.N)))
        return v, a_plane

    def __call__(self, h: dict):
        hv = np.concatenate([np.full(len(self.n[f]), float(h[f])) for f in self.fams])
        for i, t in enumerate(self.topos):
            out = self._eval(t, hv)
            if out is not None:
                if i:
                    self.topos.insert(0, self.topos.pop(i))
                break
        else:
            t = self._topology(hv)
            # a fresh topology at a degenerate vertex (more than three planes through it) may pick three planes whose
            # solve sits a hair outside a fourth: accept it at the tightness tolerance of _topology
            out = self._eval(t, hv, tol=1e-5)
            self.topos = [t] + self.topos[:7]
        v, a_plane = out
        a_fam = np.bincount(self.fam_of, weights=a_plane, minlength=len(self.fams))
        proj = v @ self.N.T
        area = {f: float(a_fam[i]) for i, f in enumerate(self.fams)}
        support = {f: float(proj[:, self.fam_of == i].max()) for i, f in enumerate(self.fams)}
        return area, float((hv * a_plane).sum() / 3.0), support

    def perimeter(self, h: dict) -> dict:
        """Total edge length of each family's faces (A) at facet distances h."""
        self(h)                                           # puts the topology of h first
        hv = np.concatenate([np.full(len(self.n[f]), float(h[f])) for f in self.fams])
        t = self.topos[0]
        v = np.linalg.solve(self.N[t["tri"]], hv[t["tri"]][..., None])[..., 0]
        seg = np.linalg.norm(v[t["a"]] - v[t["b"]], axis=1)
        per = np.bincount(self.fam_of[t["face"]], weights=seg, minlength=len(self.fams))
        return {f: float(per[i]) for i, f in enumerate(self.fams)}

    def edge_lengths(self, h: dict) -> dict:
        """Total length (A) of the polyhedron's edges per edge type "f|g" (family labels sorted) at distances h."""
        self(h)
        hv = np.concatenate([np.full(len(self.n[f]), float(h[f])) for f in self.fams])
        t = self.topos[0]
        v = np.linalg.solve(self.N[t["tri"]], hv[t["tri"]][..., None])[..., 0]
        faces = {}
        for a, b, p in zip(t["a"], t["b"], t["face"]):
            faces.setdefault((min(a, b), max(a, b)), []).append(int(p))
        out = {}
        for (a, b), ps in faces.items():
            if len(ps) != 2:
                continue
            seg = float(np.linalg.norm(v[a] - v[b]))
            if seg < 1e-9:
                continue
            key = "|".join(sorted(self.fams[self.fam_of[p]] for p in ps))
            out[key] = out.get(key, 0.0) + seg
        return out

    def wulff(self, gamma: dict, volume: float) -> dict:
        """Facet distances of the Wulff shape (h_f proportional to gamma_f) with the given volume."""
        h = {f: max(g, 1e-6) for f, g in gamma.items()}
        _a, vol, _s = self(h)
        s = (volume / vol) ** (1.0 / 3.0)
        return {f: x * s for f, x in h.items()}


# ---------------- nucleation flux ----------------
def nucleation_flux(om, ms, m0: int, kt: float, omega_m: float, P: dict, barrier: float | None = None,
                    radii=None, diff=None, e_att=None, act: float = 1.0) -> float:
    """Steady-state Becker-Doring flux (M/s) of clusters through the grid up to the nucleus m0:

        J = [ sum_{m=1}^{m0-1} 1 / (f_m c_eq(m)) ]^-1,   c_eq(m) = c0 exp(-Omega(m)/kT)   (1 M standard state),
        f_m = 4 pi (D_1 + D_m)(r_1 + r_m) N_A C_MA exp(-E_att/kT),   D_m = D_1 m^(-1/3),   r_m = (3 m Omega_m / 4 pi)^(1/3)

    Smoluchowski attachment of monomers to an m-cluster (relative diffusion, contact radius), with an optional
    incorporation barrier E_att. The monomer supply is c_eq(1) = exp(-Omega(1)/kT): the grid's m = 1 row and the pool
    share mu_MA, so Omega(1) = -kT ln(C_MA/c0). The largest term (the rarest cluster) dominates, J ~ f* c_eq(m*); the
    effective prefactor J exp(barrier/kT) is reported as J0. P['nucleation']['flux'] = 'constant' keeps J0 exp(-barrier/kT).
    radii, diff, e_att (arrays over ms, A, A^2/s, eV): the atomistic contact radii, diffusion coefficients and attachment
    barriers of each size (transport, facet_kinetics), in place of the bulk-volume sphere, D_1 m^(-1/3) and one E_att."""
    nuc = P["nucleation"]
    if nuc.get("flux", "becker_doring") == "constant":
        return float(nuc["J0"]) * math.exp(-(barrier or 0.0) / kt)
    om, ms = np.asarray(om, float), np.asarray(ms)
    sel = ms < m0
    if not sel.any() or ms[0] != 1:
        return 0.0
    # Omega(1) = -kT ln(gamma_MA C_MA): the monomer's activity; attachment counts molecules, so the concentration
    # C_MA = activity / gamma_MA (act) sets the attachment frequency and the m = 1 term
    om = om.copy()
    om[0] += kt * math.log(act)
    c1 = math.exp(-om[0] / kt)                                                    # M, the monomer concentration
    m = ms[sel].astype(float)
    if radii is not None:
        r = np.asarray(radii, float)[sel]
    else:
        r = (3.0 * m * omega_m / (4.0 * math.pi)) ** (1.0 / 3.0)                  # A
    if diff is not None:
        d = np.asarray(diff, float)[sel]
        dsum = d[0] + d
    else:
        d1 = float(P["kinetics"].get("D", 1e-5)) * 1e16                           # A^2/s
        dsum = d1 * (1.0 + m ** (-1.0 / 3.0))
    ea = (np.asarray(e_att, float)[sel] if e_att is not None
          else float(nuc.get("attach_barrier_eV", 0.0)))
    f = 4.0 * math.pi * dsum * (r[0] + r) * 6.02214076e-4 * c1 * np.exp(-np.asarray(ea) / kt)   # 1/s
    x = om[sel] / kt - np.log(f)                                                  # ln of each resistance 1/(f c_eq)
    xm = x.max()
    return float(math.exp(-(xm + math.log(np.exp(x - xm).sum()))))


# ---------------- the model ----------------
class Growth:
    def __init__(self, thermo: dict, params: dict, normals: dict | None = None):
        self.P = params                                  # the inputs: never changed (derived values live in self.F)
        self.thermo = thermo
        self.solv = params.get("solvation") or {"eps": 1.0}
        # atomistic contact and hydrodynamic radii of each node's lowest state (transport.add_geometry)
        self.geo = {}
        for nd in thermo.get("nodes", []):
            st = min(nd["states"], key=lambda s: s["E_eV"]) if nd.get("states") else None
            if st and "R_se_A" in st and "R_h_A" in st:
                self.geo[(nd["m"], nd["k"])] = (float(st["R_se_A"]), float(st["R_h_A"]))
        self.Om = float(params["material"]["Omega_m"])
        fams = FA_FAMILIES if params["kinetics"].get("size_law") == "facet_area" else params["families"]
        self.shape = Shape({f: (normals or ZB_NORMALS)[f] for f in fams})
        self._gscale = 1.0                               # gamma^0 scale of match_m0
        self.act_ma = float(params["conditions"].get("gamma_MA", 1.0))
        self.set_temperature(float(params["conditions"]["T"]))
        self.match = self._match_m0() if params["material"].get("match_m0", False) else None

    def set_temperature(self, T: float):
        """Everything that depends on T: the grid's free energies, k_conv, g_bulk, the facets' binding free energies
        and solvation-shifted gamma^0 (a copy of the input families), and the nucleus hand-over; the run calls it again
        when the temperature program (conditions.T_growth, t_switch_s, ramp_s or T_profile) changes T."""
        c = self.P["conditions"]
        self.T = float(T)
        self.kt = KB_EV * self.T
        self.th = thermo_at(self.thermo, self.T, self.solv)
        m0 = self.P["material"]["m0"]
        self.m0 = max(m for m, _k in self.th["G"]) if m0 == "auto" else int(m0)      # "auto": set by handover()
        self.k_conv = k_conv_at(self.T, float(c["k_conv_500K"]), float(c["Ea_conv"]))
        self.gb = self.g_bulk()
        self.c_sat = self.mx_sat()
        self.dg_th = binding_thermal(self.thermo, self.T)               # thermal part of binding one MX_q at T
        self.F = self._fa_families() if self.facet_area else copy.deepcopy(self.P["families"])
        for p in self.F.values():
            if "eps0" in p:                                              # binding free energy at T from the energy
                p["dG_bind"] = -float(p["eps0"]) + self.dg_th
        # in solution: a unit binding to a facet gives up the solvation of the free MX_q (the facet stays solvated),
        # and every facet gains the solvation free energy of a neutral nanocrystal surface per area
        eps = float(self.solv.get("eps", 1.0))
        self.solv_area = per_area(eps, self.solv) if eps > 1.0 else 0.0
        for p in self.F.values():
            # site classes: binding free energies from the umbrellas at 500 K (no thermal part added at other T)
            binds = [p] if "dG_bind" in p else []
            binds += list(p.get("sites", {}).values())
            g_key = "gamma_eff" if ("sites" in p and "gamma_eff" in p) else "gamma0"
            if eps > 1.0:
                for q in binds:
                    q["dG_bind"] = float(q["dG_bind"]) - self.th["dG_solv_mx"]
                p[g_key] = float(p[g_key]) + 1000.0 * self.solv_area
                if g_key == "gamma_eff":                 # the solvated MX_q: c e^(dG_solv_mx/kT) as for dG_bind
                    p[g_key] -= 10.0 * float(p.get("Gamma_eff", 0.0)) * self.th["dG_solv_mx"]
            p[g_key] = float(p[g_key]) * self._gscale
        # the per-facet law's NEB free energies are against the bare monomer in the gas phase: its solvation at T
        self.dg_solv_ma = 0.0
        if self.facet_law and eps > 1.0 and 0 in self.th["g_n"]:
            self.dg_solv_ma = self.th["g_n"][0] - thermo_at(self.thermo, self.T)["g_n"][0]
        if self.facet_law and self.P["kinetics"].get("growth_law") == "facet":
            raise ValueError("size_law = 'facets' works with growth_law 'size' or 'shape' (the facet law has its own)")
        self._R1 = 0.0                                   # the monomer's transport is set again at this T
        classical = self.P["nucleation"].get("omega", "grid") == "continuum"
        # the facet-area law hands over at the grid's top (no continuum surface energy to match)
        self.handover = self._handover() if m0 == "auto" and not classical and not self.facet_area else None
        self.sphere_used = self.sphere_setup() if self.sphere else None

    def ref_pool(self):
        """Monomer supply and MX_q left of the reference pool (half of all the precursor converted): with an
        injection, the precursor is what the flask ends up with (c_P0 V0 plus the syringe's delivery, over V_end)."""
        c = self.P["conditions"]
        inj = self.injection()
        if inj is None:
            c_p, c_x = float(c["c_P0"]), float(c["c_MX0"])
        else:
            v_end = inj["V0"] + inj["rate"] * inj["duration"]
            c_p = (float(c["c_P0"]) * inj["V0"] + inj["c"] * inj["rate"] * inj["duration"]) / v_end
            c_x = float(c["c_MX0"]) * inj["V0"] / v_end
        return 0.5 * c_p, c_x - 0.5 * c_p

    def injection(self):
        """conditions.injection: {c_syringe_M, rate_L_h, V0_L (flask volume at the start), t_start_s (default 0),
        and duration_s or volume_L (what the syringe delivers)} -> {c (M), rate (L/s), V0 (L), t0, duration (s)};
        None without an injection."""
        inj = self.P["conditions"].get("injection")
        if not inj:
            return None
        rate = float(inj["rate_L_h"]) / 3600.0
        dur = float(inj["duration_s"]) if "duration_s" in inj else float(inj["volume_L"]) / rate
        return {"c": float(inj["c_syringe_M"]), "rate": rate, "V0": float(inj["V0_L"]),
                "t0": float(inj.get("t_start_s", 0.0)), "duration": dur}

    def temperature_at(self, t: float) -> float:
        """The temperature program: conditions.T_profile [[t, T], ...] (piecewise linear), or T until t_switch_s,
        then a linear ramp over ramp_s (default 0: a step) to T_growth; plain T without either."""
        c = self.P["conditions"]
        prof = c.get("T_profile")
        if prof:
            ts, Ts = zip(*sorted((float(a), float(b)) for a, b in prof))
            return float(np.interp(t, ts, Ts))
        if c.get("T_growth") is None:
            return float(c["T"])
        t1, ramp = float(c.get("t_switch_s", 0.0)), float(c.get("ramp_s", 0.0))
        if t <= t1:
            return float(c["T"])
        if ramp <= 0 or t >= t1 + ramp:
            return float(c["T_growth"])
        return float(c["T"]) + (float(c["T_growth"]) - float(c["T"])) * (t - t1) / ramp

    def temperature_breaks(self):
        """Times where the temperature program has a kink (a step must not jump over them)."""
        c = self.P["conditions"]
        if c.get("T_profile"):
            return sorted(float(a) for a, _b in c["T_profile"])
        if c.get("T_growth") is None:
            return []
        t1 = float(c.get("t_switch_s", 0.0))
        return [t1, t1 + float(c.get("ramp_s", 0.0))]

    def _handover(self):
        """Nucleus size m0 where the grid hands over to the continuum. The grid's step chemical potential
        mu_MA + Omega(m) - Omega(m - 1) oscillates between even and odd m (magic clusters), so it is averaged over two
        steps, mu2(m) = mu_MA + [Omega(m) - Omega(m - 2)] / 2, and compared with the continuum's
        g_bulk(T) + 2 gamma_f Omega_m / h_f of a Wulff dot of m units (exposed facet, gamma at the pool's coverage),
        both at the pool of half conversion. m0 is the size, from m0_min (default: half the largest row) up, where the
        two agree best (a tie goes to the larger m); the fitted gamma is not rescaled, the mismatch left is reported."""
        try:
            pl, fam = self.chem(*self.ref_pool())
        except ValueError:
            return None
        self._arrays()
        _J, _b, _k, per_m = self.nucleation(pl["mu_MA"], pl["mu_MX"])
        gam = {f: fam[f]["gamma"] for f in self.F}
        top = max(per_m)
        rows = []
        m_lo = int(self.P["material"].get("m0_min") or top // 2)
        for m in range(max(3, min(m_lo, top)), top + 1):
            if m - 2 not in per_m:
                continue
            h = self.shape.wulff(gam, m * self.Om)
            area, _v, _s = self.shape(h)
            f = max(self.F, key=lambda x: area[x])
            cont = 2.0 * gam[f] * self.Om / h[f]
            two = pl["mu_MA"] + 0.5 * (per_m[m] - per_m[m - 2]) - self.gb
            one = pl["mu_MA"] + per_m[m] - per_m[m - 1] - self.gb
            rows.append({"m": m, "grid_step_eV": round(one, 4), "grid_two_step_eV": round(two, 4),
                         "continuum_eV": round(cont, 4), "diff_eV": round(two - cont, 4)})
        if not rows:
            return None
        best = min(rows, key=lambda r: (round(abs(r["diff_eV"]), 2), -r["m"]))
        self.m0 = best["m"]
        self._arrays()
        return {"m0": self.m0, "mismatch_eV": best["diff_eV"], "rows": rows}

    def g_bulk(self):
        """Bulk free energy per MA unit at T: material.g_bulk_T = {"T": [...], "G_eV": [...]} (the bulk phonon free
        energy, interpolated in T), else the 0 K energy plus an Einstein estimate of the bulk phonons (3 x
        atoms-per-unit modes at omega_E). material.g_bulk_shift_eV (default 0) is added: a grid at another level of
        theory than the table (e.g. with D3) takes the bulk energy difference at that level."""
        mat = self.P["material"]
        shift = float(mat.get("g_bulk_shift_eV", 0.0))
        tab = mat.get("g_bulk_T")
        if tab:
            return float(np.interp(self.T, tab["T"], tab["G_eV"])) + shift
        hw = float(mat.get("omega_E_cm1", 0.0)) * 1.239841984e-4                       # eV
        if hw <= 0:
            return float(mat["g_bulk"]) + shift
        modes = 3 * int(mat.get("atoms_per_unit", 2))
        return float(mat["g_bulk"]) + shift + modes * (0.5 * hw + self.kt * math.log1p(-math.exp(-hw / self.kt)))

    def mx_sat(self):
        """Solubility of the free MX_q at T (mol/L), or None (no limit): conditions.c_MX_sat (a measured solubility,
        any T), else from the solid, material.g_MX_solid_T = {"T": [...], "G_eV": [...]} per MX_q formula unit at the
        level of the grid's MX_q reference (+ material.g_MX_solid_shift_eV): the free MX_q saturates where its
        chemical potential reaches the solid's, c_sat = c0 e^((g_solid - g_MX)/kT) / gamma_MX."""
        c, mat = self.P["conditions"], self.P["material"]
        if c.get("c_MX_sat") is not None:
            return float(c["c_MX_sat"])
        tab = mat.get("g_MX_solid_T")
        if not tab:
            return None
        g_s = float(np.interp(self.T, tab["T"], tab["G_eV"])) + float(mat.get("g_MX_solid_shift_eV", 0.0))
        return math.exp(min((g_s - self.th["g_mx"]) / self.kt, 700.0)) / float(c.get("gamma_MX", 1.0))

    def _match_m0(self):
        """Cluster-continuum matching at the nucleus: scale every gamma^0 by s so that the continuum chemical
        potential of the m0 dot, g_bulk(T) + 2 s gamma_bar Omega_m / r0 (Wulff shape of volume m0 Omega_m, gamma_bar
        the area-weighted gamma at the reference coverage, r0 the area-weighted facet distance), equals the grid's
        step mu_step = mu_MA + Omega(m0) - Omega(m0 - 1) (best shells included) at reference chemical potentials: the
        pool at half conversion. Growth then starts where the grid says adding MA lowers Omega. PLACEHOLDER bridge."""
        try:
            pl = self.pool(*self.ref_pool())
        except ValueError:
            return None
        self._arrays()
        _J, _b, _k, per_m = self.nucleation(pl["mu_MA"], pl["mu_MX"])
        if self.m0 - 1 not in per_m:
            return None
        mu_step = pl["mu_MA"] + per_m[self.m0] - per_m[self.m0 - 1]
        _pl, fam = self.chem(*self.ref_pool())
        g = {f: fam[f]["gamma"] for f in self.F}
        h0 = self.shape.wulff(g, self.m0 * self.Om)
        area, _v, _s = self.shape(h0)
        gbar = sum(g[f] * area[f] for f in self.F) / sum(area.values())
        r0 = sum(h0[f] * area[f] for f in self.F) / sum(area.values())
        s = (mu_step - self.gb) * r0 / (2.0 * gbar * self.Om)
        if s > 0:
            self._gscale *= s                                # kept when set_temperature derives the families again
            for p in self.F.values():
                g_key = "gamma_eff" if ("sites" in p and "gamma_eff" in p) else "gamma0"
                p[g_key] = p[g_key] * s
        return {"mu_step_eV": mu_step, "g_bulk_eV": self.gb, "mu_MA_ref_eV": pl["mu_MA"], "scale": s, "applied": s > 0}

    # -- pool, coverage, surface energies --
    def pool(self, C_M, C_X):
        """The precursor pool at monomer supply C_M and MX_q left C_X, with the activity coefficients of
        conditions.gamma_MA and gamma_MX (calibration factors, 1 = ideal); MX_q above its solubility (mx_sat) is solid,
        counted in C_X as a reservoir in equilibrium with the solution."""
        c = self.P["conditions"]
        args = (max(C_M, 1e-30), max(C_X, 1e-30), float(c.get("gamma_MA", 1.0)), float(c.get("gamma_MX", 1.0)), self.c_sat)
        last = getattr(self, "_pool_last", None)              # a step's end pool is the next step's start pool
        if last is not None and last[0] is self.th and last[1] == args:
            return last[2]
        pl = pool(self.th, args[0], args[1], gamma_ma=args[2], gamma_mx=args[3], c_sat=args[4])
        self._pool_last = (self.th, args, pl)
        return pl

    def chem(self, C_M, C_X):
        pl = self.pool(C_M, C_X)
        self._pl = pl                                   # the pool of the per-facet law's rates (size_state)
        c = pl["a_MX_free"]                             # the facets' coverage sees the free MX_q's activity
        fam = {}
        for f, p in self.F.items():
            fam[f] = facet_state(p, c, self.kt)
        # surface free energy and facet areas of the Wulff dot of one unit: classical nucleation, and the continuum
        # tail of the grid's Omega(m) past m0
        if self.facet_area:                             # no surface energies: no continuum tail, no Gibbs-Thomson
            self._unit_area = self.shape(self.fa_start_h(1.0))[0]
            self._sigma = None
            return pl, fam
        gam = {f: fam[f]["gamma"] for f in self.F}
        area, _v, _s = self.shape(self.shape.wulff(gam, self.Om))
        self._unit_area = area
        self._sigma = sum(gam[f] * area[f] for f in self.F)
        return pl, fam

    # -- families with site classes (`sites`): coverage -> pool level, and the incorporation rate --
    def sited(self, f: str) -> bool:
        return "sites" in self.F[f]

    def _c_of_theta(self, f: str, th: float) -> float:
        """The free MX_q activity at which family f has the class-weighted coverage th (the dots' coverage may lag
        the pool's): bisection in ln c, cached per T and th (to 1e-6, near full coverage to 1e-4 in ln(1 - th)). At
        th = 1 to double precision (K c > ~1e16) c is the smallest such level: the open sites' share is then an upper
        bound, and the direct exchange carries the rate."""
        th = float(th)
        key = (f, round(th, 6) if th < 0.99 else round(math.log(max(1.0 - th, 1e-300)), 4), self.T)
        cache = self.__dict__.setdefault("_cth_cache", {})
        if key not in cache:
            p, kt = self.F[f], self.kt
            wts = class_weights(p)
            lk = {s: -float(q["dG_bind"]) / kt for s, q in p["sites"].items()}

            def theta(lnc):
                return sum(wts[s] / (1.0 + math.exp(min(-(lk[s] + lnc), 700.0))) for s in lk)
            lo, hi = -120.0, 120.0
            for _ in range(100):
                mid = 0.5 * (lo + hi)
                lo, hi = (mid, hi) if theta(mid) < th else (lo, mid)
            if len(cache) > 200000:
                cache.clear()
            cache[key] = math.exp(0.5 * (lo + hi))
        return cache[key]

    def _speed(self, amu: float) -> float:
        """(kT / 2 pi mass)^(1/2) in A/s: the one-sided flux per concentration through a plane."""
        return math.sqrt(1.380649e-23 * self.T / (2.0 * math.pi * amu * 1.66053907e-27)) * 1e10

    def site_rate(self, f: str, th: float, p_over_a: float | None = None) -> float:
        """Incorporation rate constant k_f (A/s) of a family with site classes at the dots' coverage th, from its
        active classes s (one-dimensional transition-state flux through each umbrella barrier):
            k_f = sum_s rho_s [A_dir v_C e^(-dW_dir,s/kT) + K_s c A_disp v_CX e^(-dW_disp,s/kT)] / (1 + K_s c)
        c the free MX_q activity of that coverage (_c_of_theta); v = (kT/2 pi mass)^(1/2) of the complex (attachment)
        and of the complex against the MX_q it displaces (reduced mass, direct exchange); A_dir, A_disp the capture
        areas of the umbrellas' restraint cylinders; rho_s per area (per_area) or, for perimeter classes, per_length
        times the family's perimeter per area p_over_a (this dot's, or perimeter_per_area_ref of the reference dot).
        A class without dW_disp has no direct exchange. Capped at the collision rate v_C."""
        p = self.F[f]
        ratio = float(p_over_a) if p_over_a is not None else float(p.get("perimeter_per_area_ref", 0.0))
        k_area, k_len, v_c = self._site_parts(f, th)
        return min(k_area + k_len * ratio, v_c)

    def _site_parts(self, f: str, th: float) -> tuple:
        """site_rate split into the per-area classes' sum and the per-length classes' sum (to be multiplied by the
        dot's perimeter per area), and the collision cap: the same for every dot at one coverage, so cached per
        (f, th, T)."""
        key = (f, float(th), self.T)
        cache = self.__dict__.setdefault("_site_cache", {})
        if key in cache:
            return cache[key]
        p, kt = self.F[f], self.kt
        ks = self.P["kinetics"].get("sites") or {}
        m_c, m_x = float(ks.get("m_complex_amu", 374.7)), float(ks.get("m_mx_amu", 183.3))
        a_dir, a_disp = float(ks.get("A_dir_A2", 14.9)), float(ks.get("A_disp_A2", 42.6))
        v_c, v_cx = self._speed(m_c), self._speed(m_c * m_x / (m_c + m_x))
        c = self._c_of_theta(f, th)
        k_area = k_len = 0.0
        for q in p["sites"].values():
            if not q.get("active", True):
                continue
            x = math.exp(min(-float(q["dG_bind"]) / kt, 700.0)) * c
            r = a_dir * v_c * math.exp(-float(q["dW_dir"]) / kt)
            if q.get("dW_disp") is not None:
                r += x * a_disp * v_cx * math.exp(-float(q["dW_disp"]) / kt)
            if "per_length" in q:
                k_len += float(q["per_length"]) * r / (1.0 + x)
            else:
                k_area += float(q["per_area"]) * r / (1.0 + x)
        if len(cache) > 100000:
            cache.clear()
        cache[key] = (k_area, k_len, v_c)
        return cache[key]

    def cohort_poa(self, co) -> dict:
        """p_over_a of a cohort, kept on it while its shape (the h object) stays the same (a frozen dot)."""
        c = co.get("_poa")
        if c is None or c[0] is not co.get("h"):
            c = co["_poa"] = (co.get("h"), self.p_over_a(co.get("h"), co["area"]))
        return c[1]

    def p_over_a(self, h: dict, area: dict) -> dict:
        """Perimeter per area of each exposed family of a dot (1/A), for the perimeter site classes ({} without
        facet distances: the reference dot's ratio is used)."""
        if h is None or not any(self.sited(f) for f in self.F):
            return {}
        # scale-free: P/A of the shape h / s times 1/s, cached per shape to 0.001 in h_f / s (s from the exposed
        # families, so a far hidden plane does not wash out their resolution)
        s = max([h[f] for f in self.F if area.get(f, 0.0) > 1e-9] or list(h.values()))
        key = tuple(round(h[f] / s, 3) for f in self.F)
        cache = self.__dict__.setdefault("_poa_cache", {})
        if key not in cache:
            hu = dict(zip(self.F, key))
            per, (au, _v, _s) = self.shape.perimeter(hu), self.shape(hu)
            if len(cache) > 50000:
                cache.clear()
            cache[key] = {f: per[f] / au[f] for f in self.F if au[f] > 1e-9}
        return {f: x / s for f, x in cache[key].items() if area.get(f, 0.0) > 1e-9}

    def k_block(self, f: str, th: float, k_f: float, p_over_a: float | None = None) -> float:
        """A family's incorporation rate (A/s) at the dots' coverage th: site_rate for a family with site classes,
        else k_f (1 - th)."""
        if self.sited(f):
            return self.site_rate(f, th, p_over_a)
        return k_f * max(1.0 - th, 1e-6)

    def k_of(self, m: float, fam: dict) -> float:
        """MX_q of a continuum Wulff dot of m units at the coverage of fam (areas scale as m^(2/3))."""
        a = getattr(self, "_unit_area", None) or {}
        s = m ** (2.0 / 3.0)
        return sum((fam[f]["theta"] * fam[f]["Gamma_max"] + self.F[f]["x_X"] / self.F[f]["site_area"]) * a.get(f, 0.0) * s
                   for f in self.F)

    # -- per-facet size law from the computed roles (kinetics.size_law = "facets") --
    @property
    def facet_law(self) -> bool:
        """Channel laws: per role ("facets") or per facet family on its area ("facet_area")."""
        return self.P["kinetics"].get("size_law") in ("facets", "facet_area") and not self.sphere

    @property
    def facet_area(self) -> bool:
        return self.P["kinetics"].get("size_law") == "facet_area" and not self.sphere

    ROLE_DIM = {"facet": 2, "edge": 1, "vertex": 0}

    @staticmethod
    def role_facets(role: str) -> dict:
        """The facet families a builder role touches and its weight w_rf on each: 1 (facet), 1/2 (edge), 1/3 (vertex)."""
        kind, fams = role.split(":")
        fs = fams.replace("+", "|").split("|")
        return {f: 1.0 / len(fs) for f in fs}

    def role_envs(self, m: float) -> list:
        """The role environments of a dot of m units: (role, N, dmu, bind) from kinetics.roles.envs, the computed dots'
        environments with their Cd counts N, Delta mu (a_X = 1) and CdCl2 binding (1 M). Between two computed sizes the
        dot is a mixture of both sizes' environments, weighted linearly in m^(-1/3) (exact at the computed sizes; a role
        missing at one size has no channel there); beyond the largest and below the smallest, that size's environments
        with N scaled as m^(d/3) (d = 2 facet, 1 edge, 0 vertex)."""
        tab = self.P["kinetics"]["roles"]["envs"]
        sizes = sorted({int(s) for r in tab.values() for s in r})
        bad = [r for r in tab if not set(self.role_facets(r)) <= set(self.F)]
        if not sizes or bad:
            raise ValueError(f"kinetics.roles.envs: no computed sizes, or roles on families not in the template: {bad}")
        x = m ** (-1.0 / 3.0)
        if m <= sizes[0] or m >= sizes[-1]:
            s0 = sizes[0] if m <= sizes[0] else sizes[-1]
            parts = [(s0, None)]
        else:
            j = next(i for i in range(1, len(sizes)) if sizes[i] >= m)
            x1, x2 = sizes[j - 1] ** (-1.0 / 3.0), sizes[j] ** (-1.0 / 3.0)
            t = (x1 - x) / (x1 - x2)
            parts = [(sizes[j - 1], 1.0 - t), (sizes[j], t)]
        out = []
        for role, by in tab.items():
            dim = self.ROLE_DIM[role.split(":")[0]]
            for s, w in parts:
                envs = by.get(s, by.get(str(s)))
                if not envs:
                    continue
                scale = (m / s) ** (dim / 3.0) if w is None else w
                out += [(role, n * scale, float(dmu), float(b)) for n, dmu, b in envs]
        return out

    def role_kappa(self, pl, a0: float) -> dict:
        """Per-site rate constants (A^3/s per monomer per A^3) of each role's two routes from the bare-CdSe NEBs
        (kinetics.roles.neb). CdSe and CdSe(CdCl2)_n with n = 1 to 3 all react with the bare-CdSe NEB kinetics, and
        the flux is on the total C_M (role_channels). No barrier -> the collision a_0 v_1 (v_1 = (kT/2 pi m_CdSe)^(1/2), a_0 the dot's area per
        site; note: v_1 with the bare CdSe mass for every monomer package for now); otherwise kinetics.roles.rate_form
        (no default: a choice for Ivan)
            second_order:  (kT/h) e^(-dG_TS/kT) / n_1M             (TS against the bare monomer at 1 M)
            first_order:   (kT/h) e^(-dG_TS/kT) / (n_1M (1 + K_P a_1)),  K_P = e^(-dG_P/kT)   (the precursor saturates)
        each capped at the collision; free energies against the bare monomer shifted by its solvation. A role without a
        NEB has no rate (None)."""
        R = self.P["kinetics"]["roles"]
        form = R.get("rate_form")
        if form not in ("second_order", "first_order"):
            raise ValueError("kinetics.roles.rate_form must be set to 'second_order' or 'first_order' (no default)")
        kt, n1, ds1 = self.kt, 6.02214076e-4, self.dg_solv_ma
        a1 = math.exp(min((pl["mu_MA"] - self.th["g_n"][0]) / kt, 700.0))
        k_coll = a0 * self._speed(191.385)

        def kappa(q):
            if q is None:
                return None
            if q.get("no_barrier"):
                return k_coll
            k = (kt / 4.135667696e-15) * math.exp(-min((float(q["dG_TS"]) - ds1) / kt, 700.0)) / n1
            if form == "first_order":
                k /= 1.0 + math.exp(min(-(float(q["dG_P"]) - ds1) / kt, 700.0)) * a1
            return min(k, k_coll)
        return {role: (kappa(q.get("dir")), kappa(q.get("disp"))) for role, q in (R.get("neb") or {}).items()}

    def role_channels(self, m: float, area: float, pl) -> dict:
        """The growth channels of a dot of m units and area A at the pool pl: one per builder role r and facet f it
        touches, N_r w_rf sites (N_r the role's Cd count), Delta mu_r(c_X) = <dmu_e>_r + kT ln a_X (a CdCl2 released per
        exchange; <>_r the plain Cd-weighted mean over the role's environments, the cost per unit of completing the
        role's layer, Ivan 2026-10-30: no channel per environment, no Boltzmann mean) and the per-site constant
            kappa_rf = (kappa_dir,r + kappa_disp,r <K_X,e>_r a_X) / (1 + K_X,eff^(f) a_X),
        K_X,e = e^(-bind_e/kT), K_X,eff^(f) the site-weighted mean K_X over the environments on f (kappa_rf is then the
        Cd-weighted mean of the environments' constants). CdSe and CdSe(CdCl2)_n with n = 1 to 3 all react with the
        bare-CdSe NEB kinetics, and the flux is on the total C_M (free CdSe plus every complex; Ivan, 2026-10-30): the
        package enters only through mu_CdSe (Delta mu_r and the pool). In solution the released CdCl2 and the bound one
        carry its solvation. Roles without a NEB route are left out (listed in self.role_flags)."""
        kt = self.kt
        dsx = self.th["dG_solv_mx"]
        envs = self.role_envs(m)
        n_tot = sum(e[1] for e in envs)
        a0 = area / n_tot if n_tot > 0 else 0.0
        kap = self.role_kappa(pl, a0)
        ax = float(pl["a_MX_free"])
        lnax = math.log(max(ax, 1e-300))
        num, den = {}, {}
        for role, n, _dmu, b in envs:
            K = math.exp(min(-(b - dsx) / kt, 700.0))
            for f, w in self.role_facets(role).items():
                num[f] = num.get(f, 0.0) + w * n * K
                den[f] = den.get(f, 0.0) + w * n
        kbar = {f: num[f] / den[f] for f in num if den[f] > 0}
        sums = {}                                       # role -> [Cd count, sum n dmu, sum n K]
        for role, n, dmu, b in envs:
            s = sums.setdefault(role, [0.0, 0.0, 0.0])
            s[0] += n
            s[1] += n * dmu
            s[2] += n * math.exp(min(-(b - dsx) / kt, 700.0))
        ch = {"role": [], "f": [], "N": [], "dmu": [], "kappa": [], "a0": a0}
        flags = self.__dict__.setdefault("role_flags", set())
        for role, (n, ndmu, nK) in sums.items():
            kd, kx = kap.get(role, (None, None))
            if kd is None and kx is None:
                flags.add(f"{role}: no NEB, no channel")
                continue
            if n <= 0:
                continue
            for f, w in self.role_facets(role).items():
                ch["role"].append(role)
                ch["f"].append(f)
                ch["N"].append(w * n)
                ch["dmu"].append(ndmu / n + dsx + kt * lnax)
                ch["kappa"].append(((kd or 0.0) + (kx or 0.0) * nK / n * ax) / (1.0 + kbar[f] * ax))
        for k_ in ("N", "dmu", "kappa"):
            ch[k_] = np.array(ch[k_], float)
        return ch

    def channel_conductance(self, ch, res_d: float) -> np.ndarray:
        """Per-channel conductance g_r = N_r / (res_d / a_0 + 1 / kappa_r) (A^3/s), res_d = delta/D the dot's diffusion
        resistance per area: the channel's flux is n g_r (C - C_eq,inf e^(Delta mu_r/kT)) units/s. growth_control
        "diffusion" drops 1/kappa_r, "reaction" drops res_d (as resistances)."""
        if not len(ch["N"]):
            return np.zeros(0)
        ctl = self.P["kinetics"].get("growth_control", "mixed")
        rd = 0.0 if ctl == "reaction" else res_d / max(ch["a0"], 1e-300)
        rr = 0.0 if ctl == "diffusion" else 1.0 / np.maximum(ch["kappa"], 1e-300)
        return ch["N"] / np.maximum(rd + rr, 1e-300)

    def facet_speeds(self, co, C, C_M, mu0, res_d) -> dict:
        """Facet velocities (A/s) of a dot in the per-facet law: v_f = Omega_m J_f / A_f, J_f = n sum_(r on f) g_r (C -
        C_eq,inf e^(Delta mu_r/kT)), C_eq,inf = C_M e^((g_bulk - mu)/kT); a facet that would shrink stands still (as in
        shape_speeds)."""
        n = 6.02214076e-4
        ch = co["_ch"]
        g = self.channel_conductance(ch, res_d)
        ceq = C_M * np.exp(np.minimum((self.gb + ch["dmu"] - mu0) / self.kt, 700.0))
        flux = n * g * (C - ceq)
        v = {}
        for f in self.F:
            sel = np.array([x == f for x in ch["f"]], bool)
            a_f = co["area"].get(f, 0.0)
            vf = self.Om * float(flux[sel].sum()) / a_f if (a_f > 1e-6 and sel.any()) else 0.0
            v[f] = vf if self.facet_area else max(vf, 0.0)          # facet-area law: a facet may dissolve
        return v

    # -- facet-area size law (kinetics.size_law = "facet_area"): one channel per facet family on its area --
    def _fa(self) -> dict:
        return self.P["kinetics"]["facet_area"]

    def fa_site_area(self, f: str) -> float:
        """Area per terrace Cd (A^2), from the lattice: a^2/2 on {100}, sqrt(3) a^2/4 on {111} (a^3 = 4 Omega_m)."""
        a = (4.0 * self.Om) ** (1.0 / 3.0)
        return a * a / 2.0 if f == "100" else SQ3 * a * a / 4.0

    def _fa_families(self) -> dict:
        """The facet-area law's families: the terrace binding pooled over the computed sizes (gas phase at 1 M;
        set_temperature adds the CdCl2 solvation), Langmuir (w = 0), one Z-type site per terrace Cd. No gamma."""
        B = self._fa()["facets"]
        out = {}
        for f in FA_FAMILIES:
            sa = self.fa_site_area(f)
            out[f] = {"dG_bind": float(B[f]["dG_bind_eV"]), "w": 0.0, "gamma0": 0.0, "Gamma_max": 100.0 / sa,
                      "site_area": sa, "x_X": 0.0}
        return out

    def fa_start_h(self, m: float) -> dict:
        """Plane distances of the start shape scaled to V = m Omega_m: the lattice-site planes of the builder's
        hand-over dot (kinetics.facet_area.start_h_A)."""
        h0 = {f: float(self._fa()["start_h_A"][f]) for f in FA_FAMILIES}
        _a, vol, _s = self.shape(h0)
        s = (max(m, 1e-12) * self.Om / vol) ** (1.0 / 3.0)
        return {f: h0[f] * s for f in FA_FAMILIES}

    def geometry(self, h: dict, m: float):
        """Facet areas (A^2), edge lengths per edge type (A) and the built dot's record of a dot with lattice-site
        planes h and m units: QD_builder's lattice-site polyhedron of the dot built on the fly (facet_area.geometry =
        "builder", builder_geometry.BuilderGeometry) or the model's own polyhedron ("shape"). A failed build falls
        back to the model's polyhedron and is counted (self.geometry_fallbacks)."""
        if self._fa().get("geometry", "builder") == "builder":
            bg = self.__dict__.get("_bgeo")
            if bg is None:
                from .builder_geometry import BuilderGeometry
                bg = self._bgeo = BuilderGeometry(self._fa()["builder"], FA_FAMILIES,
                                                  step_A=float(self._fa().get("geometry_step_A", 0.5)))
            try:
                return bg.get(h, m)
            except Exception as exc:
                if "lattice_fitted" in str(exc) and "QD_builder without" in str(exc):
                    raise
                self.geometry_fallbacks = getattr(self, "geometry_fallbacks", 0) + 1
                self.geometry_error = str(exc)[-300:]
        return self.shape(h)[0], self.shape.edge_lengths(h), None

    def fa_dmu(self, f: str, m: float) -> float:
        """<Delta mu_e>_T of facet f at m (eV, gas phase, a_X = 1; terrace environments only): the computed sizes
        interpolated linearly in m^(-1/3), held at the nearest computed size outside them (decision 2: below a
        family's first size, that size's value)."""
        tab = sorted((int(k), float(v)) for k, v in self._fa()["facets"][f]["dmu_eV"].items())
        ms, vs = [x for x, _ in tab], [y for _, y in tab]
        if m <= ms[0]:
            return vs[0]
        if m >= ms[-1]:
            return vs[-1]
        j = next(i for i in range(1, len(ms)) if ms[i] >= m)
        x, x1, x2 = m ** (-1.0 / 3.0), ms[j - 1] ** (-1.0 / 3.0), ms[j] ** (-1.0 / 3.0)
        t = (x1 - x) / (x1 - x2)
        return vs[j - 1] + t * (vs[j] - vs[j - 1])

    def _fa_rate(self, q, area_per_site: float, a1: float):
        """Rate constant per area (A/s) of one NEB route: (kT/h) e^(-dG'_TS/kT) / (n_1M (1 + K_P a_1)), K_P =
        e^(-dG'_P/kT) (facet_area.rate_form first_order; second_order without the precursor term), dG' the NEB free
        energy minus the bare monomer's solvation, per site over area_per_site, capped at the collision rate v_1;
        no barrier: v_1; None without a NEB."""
        if q is None:
            return None
        v1 = self._speed(191.385)
        if q.get("no_barrier"):
            return v1
        kt, ds1 = self.kt, self.dg_solv_ma
        k = (kt / 4.135667696e-15) * math.exp(-min((float(q["dG_TS"]) - ds1) / kt, 700.0)) / 6.02214076e-4
        if self._fa().get("rate_form", "first_order") == "first_order":
            k /= 1.0 + math.exp(min(-(float(q["dG_P"]) - ds1) / kt, 700.0)) * a1
        return min(k / area_per_site, v1)

    def fa_keff(self, f: str, area: dict, edges: dict, a_x: float, a1: float) -> float:
        """Rate constant of facet f (A/s): terrace part and the step band along its edges, each with its own Langmuir
        blocking (the same K above and below), the band's share capped at the facet,
            k_terr = (k_dir + k_disp K_f a_X) / (1 + K_f a_X),   k_step = sum_e L_e k_step,e / sum_e L_e,
            k_eff = (1 - phi) k_terr + phi k_step,   phi = min(1, lambda_0 L_f / A_f),   lambda_0 = a / sqrt(2),
        over the edge types of f that have a NEB (one site per lambda_0^2 of band)."""
        B, kt = self._fa(), self.kt
        nf = B["facets"][f].get("neb") or {}
        Kf = math.exp(min(-float(self.F[f]["dG_bind"]) / kt, 700.0))
        sa = self.fa_site_area(f)
        kd, kx = self._fa_rate(nf.get("dir"), sa, a1), self._fa_rate(nf.get("disp"), sa, a1)
        kterr = ((kd or 0.0) + (kx or 0.0) * Kf * a_x) / (1.0 + Kf * a_x)
        A_f = float(area.get(f, 0.0))
        if A_f <= 0.0:
            return kterr
        lam = (4.0 * self.Om) ** (1.0 / 3.0) / math.sqrt(2.0)
        dsx = self.th["dG_solv_mx"]
        num = den = 0.0
        for e, q in (B.get("edges") or {}).items():
            L = float(edges.get(e, 0.0))
            if f not in e.split("|") or not q.get("neb") or L <= 0.0:
                continue
            Ke = math.exp(min(-(float(q["dG_bind_eV"]) - dsx) / kt, 700.0))
            ed = self._fa_rate(q["neb"].get("dir"), lam * lam, a1)
            ex = self._fa_rate(q["neb"].get("disp"), lam * lam, a1)
            num += L * ((ed or 0.0) + (ex or 0.0) * Ke * a_x) / (1.0 + Ke * a_x)
            den += L
        if den <= 0.0:
            return kterr
        phi = min(1.0, lam * den / A_f)
        return (1.0 - phi) * kterr + phi * num / den

    def _fa_channels(self, co) -> dict:
        """One channel per facet family of a dot: N = A_f (A^2), a_0 = 1, kappa = k_eff^(f) (A/s), and
        Delta mu_f(m; a_X) = <Delta mu_e>_T(m) + dG_solv(CdCl2) + kT ln a_X (terrace environments only)."""
        pl, kt = self._pl, self.kt
        ax = float(pl["a_MX_free"])
        a1 = math.exp(min((pl["mu_MA"] - self.th["g_n"][0]) / kt, 700.0))
        dsx, lnax = self.th["dG_solv_mx"], math.log(max(ax, 1e-300))
        area, edges = co["area"], co.get("edges") or {}
        ch = {"role": [], "f": [], "N": [], "dmu": [], "kappa": [], "a0": 1.0}
        for f in FA_FAMILIES:
            ch["role"].append(f"facet:{f}")
            ch["f"].append(f)
            ch["N"].append(float(area.get(f, 0.0)))
            ch["dmu"].append(self.fa_dmu(f, co["m"]) + dsx + kt * lnax)
            ch["kappa"].append(self.fa_keff(f, area, edges, ax, a1))
        for k_ in ("N", "dmu", "kappa"):
            ch[k_] = np.array(ch[k_], float)
        return ch

    def _fa_step(self, co, m_new: float, v: dict) -> dict:
        """Plane distances of a dot taken to m_new units: h_f + s v_f, v_f = Omega_m j_f signed on exposed facets (a
        facet with j_f < 0 dissolves) and 0 on hidden planes (held), s from V(h) = m_new Omega_m (Newton). When the
        signed speeds cannot give the new volume the shape is kept and scaled."""
        h, target = co["h"], m_new * self.Om

        def keep():
            _a, vol, _s = self.shape(h)
            s = (target / vol) ** (1.0 / 3.0)
            return {f: h[f] * s for f in h}

        if abs(m_new - co["m"]) <= 1e-12 * max(co["m"], 1.0):
            return dict(h)
        if not v or not any(v.values()):
            return keep()
        area, vol, _s = self.shape(h)
        s = 0.0
        for _ in range(40):
            slope = sum(area[f] * v.get(f, 0.0) for f in h)
            if slope == 0.0:
                return keep()
            s_new = s + (target - vol) / slope
            hn = {f: h[f] + s_new * v.get(f, 0.0) for f in h}
            if s_new < 0.0 or min(hn.values()) <= 0.5:
                return keep()
            s = s_new
            area, vol, _s = self.shape(hn)
            if abs(vol - target) <= 1e-10 * target:
                return hn
        return keep()

    def _fa_dmu_dm(self, cohorts, A, r, th_bar) -> np.ndarray:
        """d mu_dot / dm of each dot under the facet law (eV per unit), by central difference at +-1 % in m with the
        areas scaled as m^(2/3) and the edges as m^(1/3): the stiffness of a dot near its own equilibrium (a dot whose
        mu rises with m sits in a stable well; run() then updates its size backward-Euler)."""
        res_d, _rr = self.resistances(A, r, th_bar)
        out = np.zeros(len(cohorts))
        for i, co in enumerate(cohorts):
            mu = []
            for fac in (0.99, 1.01):
                s = fac ** (1.0 / 3.0)
                c2 = {"m": co["m"] * fac, "area": {f: a * s * s for f, a in co["area"].items()},
                      "edges": {e: L * s for e, L in (co.get("edges") or {}).items()}}
                ch = self._fa_channels(c2)
                g = self.channel_conductance(ch, float(res_d[i]))
                if g.sum() <= 0:
                    break
                top = float(ch["dmu"].max())
                mu.append(top + self.kt * math.log(float((g * np.exp((ch["dmu"] - top) / self.kt)).sum() / g.sum())))
            if len(mu) == 2:
                out[i] = (mu[1] - mu[0]) / (0.02 * max(co["m"], 1e-12))
        return out

    FA_TAIL = 500                                       # the nucleation chain past the grid's top, in units

    def _fa_tail(self, m_e: int):
        """Delta mu_dot(m) - g_bulk-referenced chemical potential of a start-shape dot of m units under the facet law
        at the step's pool (m = m_e+1 .. FA_TAIL): kT ln [sum_f g_f e^(Delta mu_f/kT) / sum_f g_f], g_f = A_f / (r/D_1
        + 1/k_eff^(f)) (diffusion to the equivalent sphere in series with the facet), the model's own polyhedron.
        Cached on the pool's free MX_q and mu_MA (the step's speciation)."""
        pl = self._pl
        key = (m_e, round(math.log(max(float(pl["a_MX_free"]), 1e-300)), 6), round(float(pl["mu_MA"]), 4))
        c = self.__dict__.setdefault("_fa_tail_cache", {})
        if key in c:
            return c[key]
        if len(c) > 200:
            c.clear()
        kt, d1 = self.kt, (self._D1 or self.D)
        mt = np.arange(m_e + 1, max(self.FA_TAIL, m_e + 2) + 1)
        h1 = self.fa_start_h(1.0)
        area1, _v, _s = self.shape(h1)
        edges1 = self.shape.edge_lengths(h1)
        out = np.empty(len(mt))
        for i, m in enumerate(mt):
            s = m ** (1.0 / 3.0)
            co = {"m": float(m), "area": {f: a * s * s for f, a in area1.items()},
                  "edges": {e: L * s for e, L in edges1.items()}}
            ch = self._fa_channels(co)
            r = (3.0 * m * self.Om / (4.0 * math.pi)) ** (1.0 / 3.0)
            g = ch["N"] / (r / d1 + 1.0 / np.maximum(ch["kappa"], 1e-300))
            top = float(ch["dmu"].max())
            out[i] = top + kt * math.log(float((g * np.exp((ch["dmu"] - top) / kt)).sum() / max(g.sum(), 1e-300)))
        c[key] = (mt, out)
        return mt, out

    def _fa_capf(self) -> float:
        """Capacitance of the start shape over the radius of the sphere of the same volume (scale-free)."""
        if "fa" not in self._capf_cache:
            from .transport import capacitance_polyhedron
            h1 = self.fa_start_h(1.0)
            normals, dists = [], []
            for f in FA_FAMILIES:
                for nrm in ZB_NORMALS[f]:
                    normals.append(np.asarray(nrm, float) / np.linalg.norm(nrm))
                    dists.append(h1[f])
            c, _se = capacitance_polyhedron(normals, dists, n=20000)
            self._capf_cache["fa"] = c / (3.0 * self.Om / (4.0 * math.pi)) ** (1.0 / 3.0)
        return self._capf_cache["fa"]

    def _arrays(self):
        """Per-family constants as arrays (order self.fams) and the grid's G(m,k) per m."""
        self.fams = list(self.F)
        F = self.F
        self.D = float(self.P["kinetics"].get("D", 1e-5)) * 1e16                                          # A^2/s
        rows = {}
        for (m, k), g in self.th["G"].items():
            rows.setdefault(m, []).append((k, g))
        self.Gm = {m: (np.array([k for k, _ in sorted(v)]), np.array([g for _, g in sorted(v)])) for m, v in sorted(rows.items())}
        self.ms = np.array(sorted(rows))
        self.ks = np.array(sorted({k for (_m, k) in self.th["G"]}))
        self.Gp = np.full((len(self.ms), len(self.ks)), np.inf)                         # G(m, k), inf where no node
        for (m, k), g in self.th["G"].items():
            self.Gp[np.searchsorted(self.ms, m), np.searchsorted(self.ks, k)] = g
        # transport from the structures: Stokes-Einstein on the hydrodynamic radii, contact radii of the clusters,
        # the capacitance of the dot's polyhedron (kinetics.transport = "stokes"; without geometry: the fixed D)
        kin = self.P["kinetics"]
        self.stokes = kin.get("transport", "stokes") == "stokes" and bool(self.geo)
        self.eta = float(kin.get("viscosity_mPa_s", 0.4))
        self.c_st = float(kin.get("stokes_c", 6.0))
        self.Rse = np.full(self.Gp.shape, np.nan)
        self.Rh = np.full(self.Gp.shape, np.nan)
        for (m, k), (rs, rh) in self.geo.items():
            if m in self.Gm:
                self.Rse[np.searchsorted(self.ms, m), np.searchsorted(self.ks, k)] = rs
                self.Rh[np.searchsorted(self.ms, m), np.searchsorted(self.ks, k)] = rh
        self._R1, self._D1 = 0.0, self.D
        self._capf_cache = {}
        self._capf = 1.0
        self._cdot = self._ccon = None
        self.m_top = int(self.ms.max())
        # grid nucleation: dots within the grid's sizes take the grid's chemical potential (grid_mu), and the
        # nucleation chain runs on the grid up to its top before the continuum takes over
        self.grid_mu_on = (self.P["nucleation"].get("omega", "grid") == "grid"
                           and self.P["kinetics"].get("mu_small", "grid") == "grid")
        self.below = self.ms <= (self.m_top if self.grid_mu_on else self.m0)

    def grid_mu(self, mu_mx):
        """Chemical potential of a dot of m units within the grid: the two-step free-energy step of the grid at the
        free MX_q's mu, mu(m) = [G~(m) - G~(m - 2)] / 2 with G~(m) = -kT ln sum_k e^-(G(m,k) - k mu_MX)/kT (the even-odd
        oscillation of single steps averaged, as in _handover). Returns the sizes and mu (eV), m >= 3."""
        v = self.Gp - self.ks[None, :] * mu_mx
        vmin = v.min(axis=1)
        gt = vmin - self.kt * np.log(np.exp(-(v - vmin[:, None]) / self.kt).sum(axis=1))
        pos = {m: i for i, m in enumerate(self.ms.tolist())}
        ms = [m for m in self.ms.tolist() if m - 2 in pos]
        return np.array(ms, float), np.array([0.5 * (gt[pos[m]] - gt[pos[m - 2]]) for m in ms])

    def _transport(self, v, sizes):
        """Contact radii (A) and Stokes-Einstein diffusion coefficients (A^2/s) of the clusters of `sizes`: within the
        grid the shells weighted by e^(-Omega(m,k)/kT) at the pool, beyond it scaled as m^(1/3) from the grid's top.
        Also sets the monomer's R_1, D_1 (the complexes of the pool) for the growth of dots."""
        from .transport import stokes_d
        vmin = np.nanmin(np.where(np.isfinite(v), v, np.nan), axis=1)
        w = np.exp(-(v - vmin[:, None]) / self.kt)
        ok = np.isfinite(self.Rse) & np.isfinite(w)
        wsum = np.where(ok, w, 0.0).sum(1)
        rse = np.where(wsum > 0, np.where(ok, w * np.nan_to_num(self.Rse), 0.0).sum(1) / np.maximum(wsum, 1e-300), np.nan)
        rh = np.where(wsum > 0, np.where(ok, w * np.nan_to_num(self.Rh), 0.0).sum(1) / np.maximum(wsum, 1e-300), np.nan)
        have = np.flatnonzero(np.isfinite(rse))
        mg = self.ms[have]
        m = np.asarray(sizes, float)
        j = have[np.where(m >= mg[0], np.searchsorted(mg, m, side="right") - 1, 0)]
        s = np.where(m != self.ms[j], (m / self.ms[j]) ** (1.0 / 3.0), 1.0)
        if self.ms[have[0]] == 1:
            self._R1, self._D1 = float(rse[have[0]]), stokes_d(float(rh[have[0]]), self.T, self.eta, self.c_st)
        return rse[j] * s, stokes_d(rh[j] * s, self.T, self.eta, self.c_st)

    def _cap_factor(self, fam):
        """Capacitance of the dot's Wulff polyhedron over the radius of the sphere of the same volume (1 for a sphere,
        1.065 for a cube): scale-free, so one walk-on-spheres run per surface-energy ratio (cached to 1 %)."""
        gam = {f: fam[f]["gamma"] for f in self.F}
        g0 = max(min(gam.values()), 1e-9)
        key = tuple(round(gam[f] / g0, 2) for f in self.F)
        if key not in self._capf_cache:
            from .transport import capacitance_polyhedron
            h1 = self.shape.wulff(gam, self.Om)
            normals, dists = [], []
            for f in self.F:
                for nrm in (ZB_NORMALS[f] if f in ZB_NORMALS else []):
                    normals.append(np.asarray(nrm, float) / np.linalg.norm(nrm))
                    dists.append(h1[f])
            c, _se = capacitance_polyhedron(normals, dists, n=20000)
            self._capf_cache[key] = c / (3.0 * self.Om / (4.0 * math.pi)) ** (1.0 / 3.0)
        return self._capf_cache[key]

    CAP_RHO = 0.1                                       # grid in R_1 / r_eq of the contact capacitance (linear between)
    CAP_KEY = 0.05                                      # shape key of the capacitance cache: h_f / r_eq to 0.05

    def dot_capacitance(self, co):
        """The dot's hydrodynamic radius R_h = C(P) + dR_h and its contact capacitance C(P + B(R_1 + dR_se)) (A): P the
        polyhedron of its facet distances h, grown by the monomer's contact radius R_1 and the dot's own shell dR_se
        (where a monomer centre touches the dot). Scale-free: per shape (h_f / r_eq of the exposed families, to CAP_KEY)
        one walk-on-spheres run for C(P) and one per point of the (R_1 + dR_se) / r_eq grid. These depend on the
        rounded shape alone (fixed seed), so they are cached for every run of the process (SHAPE_CAPS)."""
        from .transport import capacitance_polyhedron
        r = (3.0 * max(co["m"], 1e-12) * self.Om / (4.0 * math.pi)) ** (1.0 / 3.0)
        # the dot's CdCl2 shell and van der Waals skin beyond the CdSe polyhedron (kinetics.shell_offset, a + b sigma
        # in A with sigma = k / A the MX_q per A^2, fitted on the shape library's dots): R_h = C(P) + dR_h, and the
        # monomer's contact surface P + B(R_1 + dR_se)
        so = self.P["kinetics"].get("shell_offset") or {}
        area = sum(co["area"].values())
        sig = co.get("k", 0.0) / area if area > 0 else 0.0
        d_h = (so["R_h_A"][0] + so["R_h_A"][1] * sig) if "R_h_A" in so else 0.0
        d_se = (so["R_se_A"][0] + so["R_se_A"][1] * sig) if "R_se_A" in so else 0.0
        # 0.05 in h_f / r_eq: C(P) moves by far less than the walk-on-spheres noise itself (a 10000-walker run's seed
        # shifts the diffusion-limited N_dots by up to 1 % at the end, this rounding by 0.1 %), and the default page
        # run meets 55 shapes (40 s instead of 190 s at 0.01)
        ck = co.get("_capkey")                             # kept while the shape (h object) and m stay the same
        if ck is not None and ck[0] is co["h"] and ck[1] == co["m"]:
            key = ck[2]
        else:
            key = tuple((f, round(round(co["h"][f] / r / self.CAP_KEY) * self.CAP_KEY, 6))
                        for f in self.F if co["area"].get(f, 0.0) > 0)
            co["_capkey"] = (co["h"], co["m"], key)
        e = SHAPE_CAPS.get(key)
        if e is None:
            normals, dists = [], []
            for f, hf in key:
                for nrm in ZB_NORMALS.get(f, []):
                    normals.append(np.asarray(nrm, float) / np.linalg.norm(nrm))
                    dists.append(hf)
            e = {"planes": (np.array(normals), np.array(dists)), "dil": {},
                 "c0": capacitance_polyhedron(normals, dists, n=10000)[0]}
            if len(SHAPE_CAPS) >= SHAPE_CAPS_MAX:
                SHAPE_CAPS.clear()
            SHAPE_CAPS[key] = e
        rho = (self._R1 + d_se) / r
        if rho <= 0:
            return r * e["c0"] + d_h, r * e["c0"]

        def at(i):
            if i == 0:
                return e["c0"]
            if i not in e["dil"]:
                e["dil"][i] = capacitance_polyhedron(*e["planes"], n=10000, dilate=i * self.CAP_RHO)[0]
            return e["dil"][i]
        j = int(rho / self.CAP_RHO)
        w = rho / self.CAP_RHO - j
        return r * e["c0"] + d_h, r * ((1.0 - w) * at(j) + w * at(j + 1))

    # -- facet-resolved growth (kinetics.growth_law = "facet") --
    def facet_mu(self, h, area, fam):
        """mu_f - g_bulk (eV) of each facet: the variational facet potential Omega_m (dE_s/dh_f) / (dV/dh_f),
        E_s = sum_g gamma_g A_g(h) (central differences on the polyhedron). On the Wulff shape it is 2 gamma_f
        Omega_m / h_f; off it a facet that lags behind (large area) has a low potential, one that runs ahead a high one.
        A hidden facet is probed just inside its support plane (where it would reappear)."""
        gam = {f: fam[f]["gamma"] for f in self.F}

        def es_vol(hh):
            a, v, s = self.shape(hh)
            return sum(gam[g] * a[g] for g in self.F), v, s

        _e0, _v0, sup = es_vol(h)
        # an exposed facet is probed with the hidden planes lifted clear of the dot: a hidden plane sitting right at
        # its support would otherwise cut the corners as the facet moves out and inflate its potential
        lifted = {g: (h[g] if area[g] > 1e-6 else max(h[g], sup[g] + 0.05 * max(h.values()))) for g in self.F}
        out = {}
        for f in self.F:
            if area[f] > 1e-6:
                d = 0.005 * h[f]
                hp, hm = dict(lifted), dict(lifted)
                hp[f] += d
                hm[f] -= d
            else:
                d = 0.25
                hp, hm = dict(h), dict(h)
                hp[f] = sup[f]
                hm[f] = sup[f] - 2 * d
            ep, vp, _ = es_vol(hp)
            em, vm, _ = es_vol(hm)
            out[f] = self.Om * (ep - em) / (vp - vm) if abs(vp - vm) > 1e-12 else 2.0 * gam[f] * self.Om / max(h[f], 1e-6)
        jn = self.junction(_v0 / self.Om)                  # the same junction term as the dot's mu (size_state)
        return {f: x + jn for f, x in out.items()}
    def _facet_rate(self, co, C, C_M, mu0, fam, th):
        """dm/dt and v_f of one dot at monomer level C (see _facet_setup)."""
        return self._facet_setup(co, C_M, mu0, fam, th)(C)

    def _facet_setup(self, co, C_M, mu0, fam, th):
        """A function C -> (dm/dt (units/s), facet velocities v_f (A/s)) for one dot, with its geometry, facet
        potentials and Gibbs-Thomson levels set up once (they do not depend on C). Every facet grows from
        the monomer level C_s at the dot's surface, all fed by one diffusion field, sum_f A_f v_f / Omega_m =
        4 pi (D_1 + D_dot)(C_dot + R_1) n (C - C_s), solved for C_s:
          fast facets (a lone unit sticks):  v_f = Omega_m n k_f (1 - theta_f)(C_s - C_eq,f)
          island facets (a new layer needs a 2D nucleus):  v_f = d_f min(J_f A_f, (J_f v_step^2)^(1/3)),
             J_f = (kT/h) rho_f (1 - theta_f) e^(-dG*_f/kT),  dG*_f = beta_f^2 / (4 kT ln S_f),  S_f = C_s / C_eq,f
             (beta_f: the step term of the island ladder dG(N) = g_inf N + beta_f sqrt(N)); dissolving: the fast law at
             the collision rate
        with C_eq,f = C_M e^((g_bulk + 2 gamma_f Omega_m / h_f - mu)/kT) the facet's own Gibbs-Thomson level."""
        n = 6.02214076e-4
        FK = self.P["kinetics"].get("facets") or {}
        h, area = co["h"], co["area"]
        r = (3.0 * co["m"] * self.Om / (4.0 * math.pi)) ** (1.0 / 3.0)
        if self.stokes:
            from .transport import KB
            cap, con = self.dot_capacitance(co)
            d_dot = KB * self.T / (self.c_st * math.pi * self.eta * 1e-3 * cap * 1e-10) * 1e20
            g_diff = 4.0 * math.pi * (self._D1 + d_dot) * con * n
        else:
            g_diff = 4.0 * math.pi * self.D * r * n
        kt = self.kt
        dmu = self.facet_mu(h, area, fam)
        ceq = {f: C_M * math.exp(min((self.gb + dmu[f] - mu0) / kt, 700.0)) for f in self.F}
        m_cx = (112.414 * 2 + 78.971 + 2 * 35.453) * 1.66053907e-27
        k_coll = math.sqrt(1.380649e-23 * self.T / (2 * math.pi * m_cx)) * 1e10           # A/s
        # families with site classes: the share of the collision flux that is taken up, from their active classes
        # (open sites and direct exchange); it plays the part of 1 - theta below (attachment, islands, dissolution)
        poa = self.p_over_a(h, area)
        free_s = {f: min(max(self.site_rate(f, th[f], poa.get(f)) / k_coll, 1e-300), 1.0) for f in self.F if self.sited(f)}

        def v_of(cs):
            v = {}
            for f in self.F:
                p = FK.get(f, {})
                free = free_s[f] if f in free_s else max(1.0 - th[f], 1e-6)
                drive = cs - ceq[f]
                beta = p.get("island_beta_eV")
                if beta is None or beta <= 0 or drive <= 0:
                    k = float(p.get("k_f", k_coll)) if drive > 0 and f not in free_s else k_coll
                    v[f] = self.Om * n * min(k, k_coll) * free * drive
                else:
                    d_f = self.Om / self.F[f]["site_area"]
                    lnS = math.log(cs / ceq[f])
                    dG = beta ** 2 / (4.0 * kt * lnS)
                    J = (kt / 4.135667696e-15) / self.F[f]["site_area"] * free * math.exp(-min(dG / kt, 700.0))
                    v_step = self.Om * n * k_coll * free * drive
                    v[f] = d_f * min(J * area[f], (J * v_step ** 2) ** (1.0 / 3.0)) if area[f] > 0 else 0.0
                if area[f] <= 1e-6:
                    v[f] = max(v[f], 0.0)                 # a hidden facet does not dissolve: it can only reappear
            return v

        ce = [ceq[f] for f in self.F if area[f] > 1e-6] or list(ceq.values())

        def rate(C):
            def mismatch(lncs):
                cs = math.exp(lncs)
                v = v_of(cs)
                return sum(area[f] * v[f] for f in self.F) / self.Om - g_diff * (C - cs)

            lo, hi = math.log(max(min(ce), 1e-30) * 1e-6), math.log(max(C, max(ce)) * 1.0001)
            flo, fhi = mismatch(lo), mismatch(hi)
            lncs = brentq(mismatch, lo, hi, xtol=1e-10) if flo < 0 < fhi else (lo if flo >= 0 else hi)
            v = v_of(math.exp(lncs))
            return sum(area[f] * v[f] for f in self.F) / self.Om, v
        return rate

    def nucleation(self, mu_ma, mu_mx):
        """Barrier, flux J (nucleation_flux) and the nuclei: on the grid's Omega(m) continued past m0 by the continuum
        (nucleation.omega = "grid", nuclei of m0 units with the shell k*(m0), or larger when the critical nucleus lies
        beyond the grid), or on the classical free energy of the continuum (omega = "continuum"):

            dG(m) = -(m - 1) dmu + sigma (m^(2/3) - 1),   dmu = mu_MA - g_bulk(T),
            sigma = sum_f gamma_f(C) A_f of the Wulff dot of one unit (the facet energies the growth uses)

        with the grid's own monomer level Omega(1) = -kT ln(C_MA/c0); nuclei are born at the first size past the
        critical one where dG has dropped by kT (self._m_birth), their MX_q from the continuum's coverage."""
        v = self.Gp - self.ms[:, None] * mu_ma - self.ks[None, :] * mu_mx
        vmin = v.min(axis=1)
        om = vmin - self.kt * np.log(np.exp(-(v - vmin[:, None]) / self.kt).sum(axis=1))
        if self.P["nucleation"].get("omega", "grid") == "continuum":
            return self._nucleation_continuum(mu_ma, float(om[0]), v)
        per_m = dict(zip(self.ms.tolist(), om.tolist()))
        ob, mb = om[self.below], self.ms[self.below]
        # past the grid part (m0, or the grid's top when dots there take the grid's mu) the continuum carries on,
        # Omega(m) = Omega(m_e) - (m - m_e) dmu + sigma (m^(2/3) - m_e^(2/3)), so a critical nucleus larger than the
        # grid is still crossed: nuclei are born at m0, or later where the joined Omega(m) has dropped by kT past its
        # maximum (as in the classical mode)
        sigma, dmu = getattr(self, "_sigma", None), mu_ma - self.gb
        m0, m_e = int(self.m0), int(mb[-1])
        # junction: where the dots take the grid's mu up to its top, the continuum's (2/3) sigma m^(-1/3) is joined to
        # the grid's two-step mu at m_top by a term delta m_top / m that dies away with size (growth uses the same,
        # self.junction), so no step in mu at the switch traps dots there
        self._delta, self._m_join = 0.0, m_e
        if sigma is not None and self.grid_mu_on:
            mg, mug = self.grid_mu(mu_mx)
            if len(mg) and mg[-1] == m_e:
                self._delta = float(mug[-1] - self.gb - (2.0 / 3.0) * sigma * m_e ** (-1.0 / 3.0))
                if self.sphere and hasattr(self, "sph_gam"):
                    # the growth's own mu at the junction: the sphere's 2 gamma Omega / r (the nucleation tail keeps
                    # sigma, the nucleation is the full model's)
                    r_e = (3.0 * m_e * self.Om / (4.0 * math.pi)) ** (1.0 / 3.0)
                    self._delta_sph = float(mug[-1] - self.gb - 2.0 * self.sph_gam * self.Om / r_e)
        if sigma is not None:
            if dmu <= 0 or (2.0 * sigma / (3.0 * dmu)) ** 3 > self.M_CLASSICAL:
                self._m_birth = m0
                return 0.0, float("inf"), None, per_m
            top = int(max(m_e, 2.0 * (2.0 * sigma / (3.0 * dmu)) ** 3, 3 * m_e, m_e + 10)) + 2
            mt = np.arange(m_e + 1, top + 1)
            ob = np.concatenate([ob, ob[-1] - (mt - m_e) * dmu + sigma * (mt ** (2.0 / 3.0) - m_e ** (2.0 / 3.0))
                                 + self._delta * m_e * np.log(mt / m_e)])
            mb = np.concatenate([mb, mt])
        elif self.facet_area:
            # facet-area law: past the grid's top Omega(m) continues with the dots' own chemical potential from the
            # facet law, Omega(m) = Omega(m_e) - sum_(j = m_e+1..m) [dmu - Delta mu_dot(j)] (no surface energy):
            # nuclei count only where the law lets them grow
            mt, dmu_t = self._fa_tail(m_e)
            ob = np.concatenate([ob, ob[-1] - np.cumsum(dmu - dmu_t)])
            mb = np.concatenate([mb, mt])
        # the barrier is the largest climb from a preceding minimum: a well of stable clusters in the grid (below the
        # monomers) has to be climbed out of, toward the continuum's maximum
        climb = ob - np.minimum.accumulate(ob)
        i_max = int(np.argmax(climb))
        barrier = max(0.0, float(climb[i_max]))
        # nuclei are born past the last hump: the first size beyond the barrier, kT below the highest Omega before
        # it, from which Omega never climbs again (within kT/4: the grid's even-odd ripple), and past the grid's top
        # when the dots there take the grid's mu (no dot is left in a shallow well of the grid; the Becker-Doring
        # chain runs through all of it up to the birth size)
        m_min = (m_e + 1) if (self.grid_mu_on and sigma is not None) else m0
        run_max = np.maximum.accumulate(ob)
        later_max = np.concatenate([np.maximum.accumulate(ob[::-1])[::-1][1:], [-np.inf]])
        ok = np.nonzero((np.arange(len(ob)) > i_max) & (mb >= m_min) & (ob <= run_max - self.kt)
                        & (later_max <= ob + 0.25 * self.kt))[0]
        m_b = int(mb[ok[0]]) if len(ok) else int(max(m_min, mb[-1]))
        self._m_birth = m_b
        if m_b in self.Gm and (m_b == m0 or self.grid_mu_on):                # a grid size: the grid's best shell
            k_star = int(self.ks[int(np.argmin(v[np.searchsorted(self.ms, m_b)]))])
        else:
            k_star = None                                                   # continuum coverage (k_of), set by run
        if self.stokes:
            radii, diff = self._transport(v, mb)
            J = nucleation_flux(ob, mb, m_b, self.kt, self.Om, self.P, barrier, radii=radii, diff=diff,
                                act=self.act_ma)
        else:
            J = nucleation_flux(ob, mb, m_b, self.kt, self.Om, self.P, barrier, act=self.act_ma)
        return J, barrier, k_star, per_m

    M_CLASSICAL = 3000                                  # a critical size beyond this nucleates nothing (dG* > 30 eV)

    def _nucleation_continuum(self, mu_ma, om1, v=None):
        sigma = getattr(self, "_sigma", None)
        dmu = mu_ma - self.gb
        if sigma is None or dmu <= 0:
            self._m_birth = 2
            return 0.0, float("inf"), None, {}
        m_star = (2.0 * sigma / (3.0 * dmu)) ** 3
        if m_star > self.M_CLASSICAL:
            self._m_birth = int(m_star) + 1
            return 0.0, 4.0 * sigma ** 3 / (27.0 * dmu ** 2), None, {}
        top = int(max(2.0 * m_star, m_star + 10)) + 2
        ms = np.arange(1, top + 1)
        dg = -(ms - 1) * dmu + sigma * (ms ** (2.0 / 3.0) - 1.0)
        i_star = int(np.argmax(dg))
        barrier = max(0.0, float(dg[i_star]))
        after = np.nonzero((ms > ms[i_star]) & (dg <= dg[i_star] - self.kt))[0]
        m_b = int(ms[after[0]]) if len(after) else int(ms[i_star]) + 1
        self._m_birth = max(2, m_b)
        if self.stokes and v is not None:
            radii, diff = self._transport(v, ms)
            J = nucleation_flux(om1 + dg, ms, self._m_birth, self.kt, self.Om, self.P, barrier, radii=radii, diff=diff,
                                act=self.act_ma)
        else:
            J = nucleation_flux(om1 + dg, ms, self._m_birth, self.kt, self.Om, self.P, barrier, act=self.act_ma)
        return J, barrier, None, {}

    def dot(self, h, fam):
        """m, k and the areas of a dot with facet distances h (ideal spheres carry no MX_q: k = 0)."""
        area, vol, support = self.shape(h)
        if self.facet_area:
            area = self.geometry(h, vol / self.Om)[0]
        k = 0.0 if self.sphere else sum((fam[f]["theta"] * fam[f]["Gamma_max"] + self.F[f]["x_X"] / self.F[f]["site_area"])
                                        * area[f] for f in self.F)
        return vol / self.Om, k, area, support

    # -- ideal growth (kinetics.growth_law = "sphere") --
    @property
    def sphere(self) -> bool:
        return self.P["kinetics"].get("growth_law") == "sphere"

    def sphere_setup(self):
        """The ideal-growth constants: the same nucleation, hand-over and pool as the full model, but every dot a
        sphere of one surface energy gamma (kinetics.sphere.gamma_J_m2; default: the facets' Wulff surface energy of
        the reference pool, sigma_W / ((36 pi)^(1/3) Omega^(2/3))), no MX_q shell on the dots, and one surface
        reaction rate k_r (kinetics.sphere.k_r_A_s; default: the facet rates over that Wulff dot)."""
        sp = self.P["kinetics"].get("sphere") or {}
        pl, fam = self.chem(*self.ref_pool())
        g_ref = self._sigma / ((36.0 * math.pi) ** (1.0 / 3.0) * self.Om ** (2.0 / 3.0))     # eV/A^2
        self.sph_gam = float(sp["gamma_J_m2"]) * 0.0624150907 if sp.get("gamma_J_m2") is not None else g_ref
        if sp.get("k_r_A_s") is not None:
            self.sph_kr = float(sp["k_r_A_s"])
        else:
            h = self.shape.wulff({f: fam[f]["gamma"] for f in self.F}, self.Om)
            k = self.mean_keff([{"N": 1.0, "area": self.shape(h)[0], "h": h}], {f: fam[f]["theta"] for f in self.F})
            self.sph_kr = float(k) if k is not None and math.isfinite(k) else 1e4
        return {"gamma_J_m2": self.sph_gam / 0.0624150907, "k_r_A_s": self.sph_kr}

    def size_state(self, cohorts, fam, th):
        """Per cohort: total area A (A^2), equivalent radius r (A), the dot's chemical potential mu_dot - g_bulk (eV)
        of its Wulff shape at gamma(C), mu_dot = g_bulk + (2/3) Omega_m sum_f gamma_f A_f / V, and the area-weighted
        Z-type coverage of its facets (th: the dots' coverage per family). With kinetics.k_r = "facets" also each
        dot's surface reaction from the facet rates over the facets that survive in the kinetic shape (facet_keff),
        k_eff = sum_f A_f k_f (1 - theta_f) / sum_f A_f (self._keff, A/s)."""
        n = len(cohorts)
        A, r, dmu_dot, th_bar = (np.zeros(n) for _ in range(4))
        facets_kr = self.P["kinetics"].get("k_r") == "facets" and not self.facet_law
        if facets_kr:
            FK = self.P["kinetics"].get("facets") or {}
            missing = [f for f in self.F if "k_f" not in FK.get(f, {}) and not self.sited(f)]
            if missing:
                raise ValueError(f"k_r = 'facets' needs kinetics.facets with k_f for every family (missing {missing})")
        self._keff = np.zeros(n) if facets_kr else None
        surv = self.kinetic_survivors(th) if facets_kr else None
        if self.sphere:                                   # ideal growth: spheres of one gamma, no shell
            m = np.array([max(co["m"], 1e-12) for co in cohorts])
            r = (3.0 * m * self.Om / (4.0 * math.pi)) ** (1.0 / 3.0)
            dmu = 2.0 * self.sph_gam * self.Om / r + np.array([self.junction(x) for x in m])
            self._keff = None
            self._cdot, self._ccon = r, r + self._R1                  # a sphere's capacitance is its radius
            return 4.0 * math.pi * r * r, r, dmu, np.zeros(n)
        if self.stokes:                                   # each dot's own polyhedron (dot_capacitance)
            caps = [self.dot_capacitance(co) for co in cohorts]
            self._cdot = np.array([c[0] for c in caps])
            self._ccon = np.array([c[1] for c in caps])
        for i, co in enumerate(cohorts):
            area = co["area"]
            a = sum(area.values())
            A[i] = a
            r[i] = (3.0 * co["m"] * self.Om / (4.0 * math.pi)) ** (1.0 / 3.0)
            dmu_dot[i] = (2.0 / 3.0) * sum(fam[f]["gamma"] * area[f] for f in self.F) / max(co["m"], 1e-12) \
                + self.junction(co["m"])
            th_bar[i] = sum(th[f] * area[f] for f in self.F) / a if a > 0 else 0.0
            if facets_kr:
                self._keff[i] = self.facet_keff(area, th, self.cohort_poa(co), surv)
        if self.facet_law and n:
            # Ivan's per-facet law from the computed roles: each channel r (builder role on facet f) has its own
            # Delta mu_r (Cd-weighted mean) and kappa_r, in series with the dot's diffusion resistance (res_d per area, res_d / a_0 per
            # site); the dot's rate n sum_r g_r (C - C_eq,r) is G (C - C_eq) with C_eq at the g_r-weighted mean of
            # e^(Delta mu_r/kT) (resistances returns the matching reaction resistance)
            self._chs = [self._fa_channels(co) if self.facet_area else self.role_channels(co["m"], A[i], self._pl)
                         for i, co in enumerate(cohorts)]
            res_d, _res_r = self.resistances(A, r, th_bar)
            for i, co in enumerate(cohorts):
                ch = co["_ch"] = self._chs[i]
                g = self.channel_conductance(ch, res_d[i])
                if g.sum() > 0:
                    top = float(ch["dmu"].max())
                    dmu_dot[i] = top + self.kt * math.log(float((g * np.exp((ch["dmu"] - top) / self.kt)).sum() / g.sum()))
        return A, r, dmu_dot, th_bar

    def mean_keff(self, cohorts, th):
        """Number-weighted surface reaction rate of the dots (A/s): k_r, or with k_r = "facets" each dot's
        sum_f A_f k_f (1 - theta_f) / A; None without dots."""
        N = sum(co["N"] for co in cohorts)
        if N <= 0:
            return None
        if self.facet_area:                              # facet-area law: sum_f A_f k_eff^(f) / A of each dot
            tot = 0.0
            for co in cohorts:
                a = sum(co["area"].values())
                if a > 0:
                    ch = self._fa_channels(co)
                    tot += co["N"] * float((ch["N"] * ch["kappa"]).sum()) / a
            return tot / N
        if self.facet_law:                               # the per-facet law: sum_e N_e kappa_e / A of each dot
            tot = 0.0
            for co in cohorts:
                a = sum(co["area"].values())
                if a > 0:
                    ch = self.role_channels(co["m"], a, self._pl)
                    tot += co["N"] * float((ch["N"] * ch["kappa"]).sum()) / a
            return tot / N
        k_r = self.P["kinetics"].get("k_r")
        if k_r != "facets":
            return float(k_r) if k_r is not None else None
        tot = 0.0
        for co in cohorts:
            if sum(co["area"].values()) > 0:
                tot += co["N"] * self.facet_keff(co["area"], th, self.p_over_a(co.get("h"), co["area"]))
        return tot / N

    KIN_SPREAD = 1e3                                     # facets more than this much faster than the slowest: far out

    def kinetic_survivors(self, th):
        """The facet families of the kinetic shape: the polyhedron with each family at a distance proportional to its
        incorporation rate v_f = k_f (1 - theta_f). A family with no area there grows out of a dot (it advances faster
        than its neighbours allow) and takes no part in the dot's surface reaction. Cached per coverage (to 0.001)."""
        # a family with site classes sits near full coverage, where 1 - theta sets the rate: keyed on its log
        key = tuple(round(math.log(max(1.0 - float(th[f]), 1e-300)), 1) if self.sited(f) else round(float(th[f]), 3)
                    for f in self.F) + (self.T,)                 # site rates depend on T as well
        cache = self.__dict__.setdefault("_kin_cache", {})
        if key not in cache:
            FK = self.P["kinetics"].get("facets") or {}
            v = {f: max(self.site_rate(f, th[f]), 1e-300) if self.sited(f) else float(FK[f]["k_f"]) * max(1.0 - th[f], 1e-6)
                 for f in self.F}
            vmin = min(v.values())
            area, _v, _s = self.shape({f: min(v[f] / vmin, self.KIN_SPREAD) for f in self.F})
            tot = sum(area.values())
            cache[key] = {f for f in self.F if area[f] > 1e-9 * tot}
        return cache[key]

    def facet_keff(self, area, th, poa=None, surv=None):
        """A dot's surface reaction rate (A/s) from the facet incorporation rates, over the facets it exposes that
        survive in the kinetic shape (kinetic_survivors): sum A_f k_f (1 - theta_f) / sum A_f (site_rate for a family
        with site classes, at this dot's perimeter per area poa). A facet much faster than its neighbours grows out
        and does not add its rate; when none of the exposed facets survives (a dot held in a shape the kinetics would
        change), the slowest exposed facet sets the rate."""
        FK = self.P["kinetics"].get("facets") or {}
        poa = poa or {}
        v = {f: self.site_rate(f, th[f], poa.get(f)) if self.sited(f) else float(FK[f]["k_f"]) * (1.0 - th[f])
             for f in self.F}
        exposed = [f for f in self.F if area.get(f, 0.0) > 0]
        if not exposed:
            return 0.0
        surv = self.kinetic_survivors(th) if surv is None else surv     # the same for every dot at th
        keep = [f for f in exposed if f in surv]
        if not keep:
            return min(v[f] for f in exposed)
        return sum(v[f] * area[f] for f in keep) / sum(area[f] for f in keep)

    def _equilibrate_coverage(self, cohorts: list, C_M: float, C_X: float, th: dict, fam: dict):
        """The dots' Z-type coverage in equilibrium with the free MX_q it leaves in solution: when moving the
        coverage to the Langmuir value of the step's pool would take more than 2 % of the solution's MX_q (the dots'
        surfaces hold more MX_q than the solution), MX_q is exchanged between the solution and the dots until
        theta_f = theta_f(a_MX_free(C_M, C_X')), C_X' = C_X - sum_f (theta_f - th_f) Gamma_max,f A_f (A_f the dots'
        total area of f). The step's explicit move to the start-of-step Langmuir value overshoots there, and no dt
        keeps the pool within dmu_step_kT. Updates the cohorts' k in place (amounts conserved exactly); returns
        (C_X', theta) or None when not needed."""
        A = {f: sum(co["N"] * co["area"].get(f, 0.0) for co in cohorts) for f in self.F}
        gm = {f: fam[f]["Gamma_max"] for f in self.F}
        up = sum((fam[f]["theta"] - th[f]) * gm[f] * A[f] for f in self.F)
        if abs(up) <= 0.02 * C_X:
            return None

        def cx_of(lna):
            c = math.exp(lna)
            t = {f: facet_state(self.F[f], c, self.kt)["theta"] for f in self.F}
            return C_X - sum((t[f] - th[f]) * gm[f] * A[f] for f in self.F), t

        def res(lna):
            cx, _t = cx_of(lna)
            if cx <= 0:
                return -1e3
            return math.log(max(self.pool(C_M, cx)["a_MX_free"], 1e-300)) - lna
        a0 = math.log(max(self.pool(C_M, C_X)["a_MX_free"], 1e-300))
        lo, hi = a0 - 5.0, a0 + 5.0
        for _ in range(40):
            if res(lo) > 0:
                break
            lo -= 5.0
        for _ in range(40):
            if res(hi) < 0:
                break
            hi += 5.0
        lna = brentq(res, lo, hi, xtol=1e-6)
        cx, t = cx_of(lna)
        for co in cohorts:
            co["k"] += sum((t[f] - th[f]) * gm[f] * co["area"].get(f, 0.0) for f in self.F)
        return cx, t

    def _merge(self, cohorts: list):
        """Merge neighbouring cohorts (in ln m) that are within material.merge_tol (1 %) of each other, and the
        closest ones while there are more than material.max_cohorts (64): N adds up, m and k are the number-weighted
        means (amounts conserved exactly), h the number-weighted mean scaled to the merged volume."""
        tol = float(self.P["material"].get("merge_tol", 0.01))
        cap = int(self.P["material"].get("max_cohorts", 64))
        # only dots well past the nucleus (m > 4 m_top): near the birth size, whether a cohort grows or dissolves
        # hangs on small differences in m, and merging there changes the outcome
        m_min = 4.0 * max(getattr(self, "m_top", 0), self._m_birth if hasattr(self, "_m_birth") else 0)
        while True:
            big = sorted((co for co in cohorts if co["m"] > m_min), key=lambda co: co["m"])
            if len(big) < 2:
                break
            gaps = [math.log(big[i + 1]["m"] / big[i]["m"]) for i in range(len(big) - 1)]
            i = int(np.argmin(gaps))
            if gaps[i] > tol and len(cohorts) <= cap:
                break
            a, b = big[i], big[i + 1]
            cohorts.remove(b)
            n = a["N"] + b["N"]
            m = (a["N"] * a["m"] + b["N"] * b["m"]) / n
            h = {f: (a["N"] * a["h"][f] + b["N"] * b["h"][f]) / n for f in self.F}
            _area, vol, _sup = self.shape(h)
            sc = (m * self.Om / vol) ** (1.0 / 3.0)
            h = {f: x * sc for f, x in h.items()}
            area = self.shape(h)[0]
            if self.facet_area:
                area, edges, built = self.geometry(h, m)
                a.update(edges=edges, build=built and {"dc": built["dc"], "m": built["m"], "h_err_A": built["h_err_A"]})
            a.update(N=n, m=m, k=(a["N"] * a["k"] + b["N"] * b["k"]) / n, h=h, area=area,
                     t_b=merged_birth(a, b))
            a.pop("_speeds", None)                         # a new shape: its facet speeds are computed again

    def junction(self, m: float) -> float:
        """The continuum's mu past the grid's top joined to the grid's (nucleation): delta m_top / m (eV), 0 within
        the grid or without grid mu."""
        d, mj = getattr(self, "_delta", 0.0), getattr(self, "_m_join", None)
        if self.sphere:
            d = getattr(self, "_delta_sph", 0.0)
        return d * mj / m if (d and mj and m > mj) else 0.0

    def size_rates(self, A, r, c_eq, th_bar, C):
        """dm/dt (CdSe units per second per dot) at monomer level C (M): diffusion to the dot (equivalent sphere) in
        series with the surface reaction, for growth and dissolution alike,

            dm/dt = A n (C - C_eq) / (A / (4 pi D r) + 1 / k_eff),   n = 6.022e-4 C molecules/A^3 per M,

        k_eff = k_r (1 - theta_bar): the Z-type shell blocks the free sites, for attachment and detachment alike
        (C_eq carries the thermodynamics, the prefactor is the same both ways); k_r None: the diffusion limit."""
        return self.conductance(A, r, th_bar) * (C - c_eq)

    def conductance(self, A, r, th_bar):
        """dm/dt per M of driving concentration, A n / (A / (4 pi D r) + 1 / k_eff) (size_rates). With the atomistic
        transport (self.stokes) the sphere becomes the dot's polyhedron: 4 pi (D_1 + D_dot) C(P + B(R_1)), the
        capacitance of each dot's own polyhedron grown by the monomer's contact radius (dot_capacitance), D_dot by
        Stokes-Einstein on the dot's capacitance C_dot, R_1 and D_1 the monomer complexes'."""
        res_d, res_r = self.resistances(A, r, th_bar)
        return A * 6.02214076e-4 / (res_d + res_r)

    def resistances(self, A, r, th_bar):
        """The two resistances in series of the size law (s/A): diffusion to the dot, A / (4 pi D r) (or the
        polyhedron's 4 pi (D_1 + D_dot)(C_dot + R_1)), and the surface reaction 1 / (k_r (1 - theta_bar))."""
        r = np.maximum(r, 1e-6)
        if self.stokes:
            from .transport import KB
            if self._cdot is not None and np.size(self._cdot) == np.size(r):
                cap, con = self._cdot, self._ccon                  # the dots' own polyhedra (size_state)
            else:
                cap = self._capf * r                              # the Wulff shape's factor on the equal-volume sphere
                con = cap + self._R1
            d_dot = KB * self.T / (self.c_st * math.pi * self.eta * 1e-3 * cap * 1e-10) * 1e20     # A^2/s
            res_d = A / (4.0 * math.pi * (self._D1 + d_dot) * con)
        else:
            res_d = A / (4.0 * math.pi * self.D * r)                                     # s/A
        k_r = self.P["kinetics"].get("k_r")                                             # A/s
        ctl = self.P["kinetics"].get("growth_control", "mixed")
        chs = getattr(self, "_chs", None)
        if self.facet_law and np.size(res_d) == 0:
            res_r = 0.0 * res_d                                                           # no dots yet
        elif self.facet_law and chs is not None and len(chs) == np.size(res_d):
            # per-facet law: sum_r g_r = A / (res_d + res_r), res_r the dot's reaction resistance
            rd = np.broadcast_to(0.0 * res_d if ctl == "reaction" else res_d, (len(chs),))
            g = np.array([self.channel_conductance(ch, res_d[i]).sum() for i, ch in enumerate(chs)])
            res_r = np.maximum(np.asarray(A) / np.maximum(g, 1e-300) - rd, 0.0)
        elif self.facet_law:
            raise ValueError("size_law = 'facets': size_state must run first for these dots")
        elif self.sphere:
            res_r = 0.0 * res_d + 1.0 / self.sph_kr                                       # ideal: one k_r, no shell
        elif k_r == "facets":
            # each dot's own surface reaction from the facet incorporation rates (size_state): computed, not set
            keff = getattr(self, "_keff", None)
            if keff is None or len(keff) != np.size(res_d):
                raise ValueError("k_r = 'facets': size_state must run first for these dots")
            res_r = 1.0 / np.maximum(keff, 1e-30)
        else:
            res_r = 1.0 / (float(k_r) * np.maximum(1.0 - th_bar, 1e-6)) if k_r is not None else 0.0 * res_d
        # kinetics.growth_control: "mixed" (both in series, the default), "diffusion" (the surface reaction is fast:
        # diffusion-limited growth) or "reaction" (transport is fast: reaction-limited growth)
        if ctl == "diffusion":
            res_r = 0.0 * res_d
        elif ctl == "reaction":
            res_d = 0.0 * res_d
        return res_d, res_r

    # -- facet shape within the size law (kinetics.growth_law = "shape") --
    def shape_speeds(self, co, dmu_s, fam, th):
        """Relative speeds of the facets of one dot (the size law sets how much it grows, these only how the growth
        is shared): v_f = rho_f k_f (1 - theta_f) [1 - e^(-(mu_s - mu_f)/kT)] P_f, with mu_s - g_bulk = dmu_s the
        monomer's potential at the dot's surface, mu_f the variational facet potential (facet_mu), k_f the facet's
        incorporation rate (kinetics.facets, default the same for every facet) and P_f the layer-nucleation factor
        e^(-dG*_f/kT), dG*_f = dG(min(n_face, N*_f)), dG(n) = beta_f sqrt(n) - n kT ln S_f, N*_f = (beta_f /
        (2 kT ln S_f))^2: the critical island, or the whole face when it is smaller (continuous in the face's size)."""
        FK = self.P["kinetics"].get("facets") or {}
        # the speeds change slowly: reused while the dot's m moves < 5 % and mu_s < 2 meV (the facet potentials
        # take nine polyhedron evaluations a dot)
        # ... and, for families with site classes (whose rate is steep near full coverage), while ln(1 - theta) moves
        # less than 0.1
        th_key = tuple(round(math.log(max(1.0 - float(th[f]), 1e-300)), 1) for f in self.F if self.sited(f))
        c = co.get("_speeds")
        if c and abs(math.log(co["m"] / c[0])) < 0.05 and abs(dmu_s - c[1]) < 2e-3 and c[2] == self.T and c[4] == th_key:
            return c[3]
        dmu = self.facet_mu(co["h"], co["area"], fam)
        poa = self.p_over_a(co["h"], co["area"])
        v = {}
        for f in self.F:
            p = FK.get(f, {})
            drive = (dmu_s - dmu[f]) / self.kt
            if drive <= 0:
                v[f] = 0.0
                continue
            v[f] = self.k_block(f, th[f], float(p.get("k_f", 1.0)), poa.get(f)) * -math.expm1(-drive) / self.F[f]["site_area"]
            beta = p.get("island_beta_eV")
            if beta:
                beta *= float(self.P["kinetics"].get("beta_scale", 1.0))   # calibration against experiments
                # a new layer climbs dG(n) = beta sqrt(n) - n kT ln S_f up to the critical island N* or, on a smaller
                # face, up to the face's own n units (a hidden facet: one unit); barrier dG(min(n_face, N*)) >= 0
                n_face = max(co["area"][f] / (len(ZB_NORMALS.get(f, [1])) * self.F[f]["site_area"]), 1.0)
                n = min(n_face, (beta / (2.0 * self.kt * drive)) ** 2)
                barrier = max(beta * math.sqrt(n) - n * self.kt * drive, 0.0)
                v[f] *= math.exp(-min(barrier / self.kt, 700.0))
        co["_speeds"] = (co["m"], dmu_s, self.T, v, th_key)
        return v

    def shape_step(self, co, m_new, v):
        """Facet distances of a dot grown to m_new units, each facet advanced in proportion to v_f (shape_speeds):
        h_f + s v_f / max(v) with s set by the volume; a dissolving dot, or one whose facets all stand still, keeps
        its shape. A hidden facet stays at most 1 A outside its support plane."""
        h = co["h"]
        if abs(m_new - co["m"]) <= 1e-10 * co["m"]:
            return h, co["area"]                         # a frozen dot: the shape moves by < 1e-10 (3 polyhedra saved)
        exposed = [f for f in self.F if co["area"][f] > 1e-6]
        vm = max((v[f] for f in exposed), default=0.0)        # only exposed facets add volume

        def keep_shape():
            s = (m_new / co["m"]) ** (1.0 / 3.0)
            hn = {f: h[f] * s for f in self.F}
            return hn, self.shape(hn)[0]

        if m_new <= co["m"] or vm <= 0:
            return keep_shape()
        # relative to the fastest exposed facet; a hidden one faster than 10x that stays hidden anyway (a plane outruns
        # its support at most sqrt(3) times faster), and the cap keeps the polyhedron finite
        u = {f: min(v[f] / vm, 10.0) for f in self.F}
        target = m_new * self.Om

        # Newton on the volume: dV/ds = sum_f A_f u_f, with the areas the polyhedron returns anyway (2-3 evaluations)
        area, vol, support = self.shape(h)
        s = 0.0
        for _ in range(30):
            slope = sum(area[f] * u[f] for f in self.F)
            if slope <= 0:
                return keep_shape()
            s_new = s + (target - vol) / slope
            if s_new < 0 or s_new > 10.0 * max(h.values()):
                return keep_shape()
            s = s_new
            hn = {f: h[f] + s * u[f] for f in self.F}
            area, vol, support = self.shape(hn)
            if abs(vol - target) <= 1e-10 * target:
                break
        else:
            return keep_shape()
        return {f: min(hn[f], support[f] + 1.0) for f in self.F}, area

    def _stop_at_equilibrium(self, m, m_new, dmu_pool, mg, mug, sigma):
        """A dot moving from m to m_new in one step stops where its own mu_dot meets the pool's mu (dmu_pool =
        mu - g_bulk): within the grid mu_dot(m) has kinks (grid_mu, linear between the grid sizes) and a jump to the
        continuum's (2/3) sigma m^(-1/3) past the grid's top, and an explicit step would jump across a stable size."""
        def f(x):
            x = np.asarray(x, float)
            cont = (2.0 / 3.0) * sigma * np.maximum(x, 1e-6) ** (-1.0 / 3.0)
            grid = np.interp(x, mg, mug - self.gb)
            return dmu_pool - np.where((x >= mg[0]) & (x <= self.m_top), grid, cont)
        m, m_new = float(m), float(m_new)
        lo, hi = sorted((m, m_new))
        inner = [x for x in range(int(math.floor(lo)) + 1, int(math.ceil(hi))) if lo < x < hi]
        if lo < self.m_top < hi:
            inner += [self.m_top + 1e-9]
        pts = np.array(sorted(inner, reverse=m_new < m))
        pts = np.concatenate([[m], pts, [m_new]])
        fv = f(pts)
        s0 = np.sign(fv[0])
        for j in range(1, len(pts)):
            if np.sign(fv[j]) != s0:
                a, b = pts[j - 1], pts[j]
                fa, fb = fv[j - 1], fv[j]
                return float(a + (b - a) * fa / (fa - fb)) if fa != fb else float(a)
        return m_new

    # -- integration --
    def run(self, log=None, seed=None) -> dict:
        """Explicit in the precursor and the dots' sizes, implicit in the monomer pool: each step solves
        C_M = C_M0 + dt [Q - m_birth J(mu(C_M)) - sum_c N_c dm_c/dt(C_M)] with the free MX_q held at its value at the
        start of the step (mu_MA(C) = mu_MA0 + kT ln(C / C_M0)); every dot keeps the Wulff shape of the current
        gamma(C); the pools are then updated from the actual change of m and k of every cohort, so the totals are
        conserved exactly. seed: dots present at t = 0, [{"N": M, "m": units}, ...] (Wulff at the initial gamma and
        coverage, their MX_q on top of c_MX0), for ripening and focusing checks without nucleation."""
        self._arrays()
        self.steps = []                                                                  # (t, dt, limit) per step
        self.pools = []                                     # (t, T, C_M, C_X) at the end of every step (a record only)
        c = self.P["conditions"]
        P_, C_M, C_X = float(c["c_P0"]), float(c.get("c_M0", 1e-12)), float(c["c_MX0"])     # c_M0: monomer at t = 0
        t, t_end = 0.0, float(c["t_end"])
        outs = list(np.geomspace(float(c["t_start"]), t_end, int(c["n_out"])))
        cohorts, series, snaps = [], [], []
        pending = pending_m = pending_k = 0.0           # this interval's nuclei (M), the MA and MX_q they took (M)
        th_dots = None                                  # Z-type coverage of the dots' facets per family
        m_dis = self.P["material"]["m_dissolve"]
        if self.grid_mu_on and self.P["nucleation"].get("omega", "grid") == "grid":
            # nuclei are born past the grid's top (nucleation); a dot that shrinks back into the grid's sizes goes back
            # to the pool, where the Becker-Doring chain accounts for clusters, instead of sitting in a shallow well
            m_dis = max(m_dis, self.m_top + 1)
        if self.facet_area:
            # facet-area law: nuclei are born at the grid's top and live from there (their mu is the facet law's at
            # every size; the grid only nucleates); a dot that shrinks below the top goes back to the pool
            m_dis = max(float(self.P["material"]["m_dissolve"]), float(getattr(self, "m_top", 0)) - 0.5)
        grid_dots = m_dis <= getattr(self, "m_top", 0) and not self.facet_area  # dots may live at grid sizes
        if seed:
            _pl, fam = self.chem(C_M, C_X)
            gam = {f: fam[f]["gamma"] for f in self.F}
            for sd in seed:
                h = (self.fa_start_h(float(sd["m"])) if self.facet_area
                     else self.shape.wulff(gam, float(sd["m"]) * self.Om))
                m_, k_, area, _support = self.dot(h, fam)
                co = {"t_b": 0.0, "N": float(sd["N"]), "h": h, "m": m_, "k": k_, "area": area}
                if self.facet_area:
                    co["edges"] = self.geometry(h, m_)[1]
                cohorts.append(co)
        seeded = {"m": sum(co["N"] * co["m"] for co in cohorts), "k": sum(co["N"] * co["k"] for co in cohorts)}
        inj = self.injection()
        V0 = inj["V0"] if inj else 1.0                  # L (concentrations are per litre of the flask at time t)
        V, n_inj, n_conv = V0, 0.0, 0.0                 # flask volume, moles injected, moles of precursor converted
        n_nuc, track = 0.0, [(0.0, 0.0, 0.0)]           # moles of nuclei so far; (t, moles of nuclei, J V) per step
        t_last_birth = 0.0
        birth_frac = float(self.P["material"].get("birth_interval", 0.05))   # new cohort every 5 % of t
        # the step solves the pool implicitly in C_M at the step's start speciation (mu_MX, the MX_q complexes): a step
        # whose end pool moves mu_MA or mu_MX by more than dmu_step_kT kT is halved, and the next one starts from it
        dmu_max = float(self.P["material"].get("dmu_step_kT", 0.2))
        dt_pool = None
        birth_share = float(self.P["material"].get("birth_share", 0.05))     # or at 5 % of the dots present
        t_breaks = self.temperature_breaks() + ([inj["t0"], inj["t0"] + inj["duration"]] if inj else [])
        while outs:
            T_now = self.temperature_at(t)
            if abs(T_now - self.T) > 1.0 or (T_now != self.T and any(abs(t - b) < 1e-12 for b in t_breaks)):
                self.set_temperature(T_now)                 # the temperature program: free energies, rates at T
                self._arrays()
            pl, fam = self.chem(C_M, C_X)
            mu0, mu_mx = pl["mu_MA"], pl["mu_MX"]
            if self.stokes:
                self._capf = 1.0 if self.sphere else (self._fa_capf() if self.facet_area else self._cap_factor(fam))
                if cohorts and not self._R1:
                    self.nucleation(mu0, mu_mx)                                       # sets the monomer's R_1, D_1
            if th_dots is None:
                th_dots = {f: fam[f]["theta"] for f in self.F}
            elif cohorts and not self.sphere:
                eq = self._equilibrate_coverage(cohorts, C_M, C_X, th_dots, fam)
                if eq is not None:
                    C_X, th_dots = eq
                    pl, fam = self.chem(C_M, C_X)
                    mu0, mu_mx = pl["mu_MA"], pl["mu_MX"]
            Q = self.k_conv * P_ * C_X
            Nc = np.array([co["N"] for co in cohorts])
            # size law: each dot's m changes by dm/dt at one chemical potential; at fixed mu_MX the monomer level in
            # equilibrium with a dot, c_eq = C_M e^((mu_dot - mu)/kT), is the same for every trial C of the solve
            Ac, rc, dmu_dot, th_bar = self.size_state(cohorts, fam, th_dots)
            mc = np.array([co["m"] for co in cohorts])
            dmu_dm = -dmu_dot / (3.0 * np.maximum(mc, 1.0))                  # continuum: mu_dot - g_bulk ~ m^(-1/3)
            if self.facet_area and len(cohorts):
                dmu_dm = self._fa_dmu_dm(cohorts, Ac, rc, th_bar)          # the facet law's own d mu_dot / dm
            mg = mug = None
            if self.grid_mu_on and len(cohorts) and not self.facet_area:
                # within the grid's sizes a dot's chemical potential is the grid's own (grid_mu, linear in m between
                # the grid sizes); beyond the grid's top, the continuum's Gibbs-Thomson
                mg, mug = self.grid_mu(mu_mx)
                inside = (mc >= mg[0]) & (mc <= self.m_top)
                if inside.any():
                    i = np.clip(np.searchsorted(mg, mc[inside]), 1, len(mg) - 1)
                    slope = (mug[i] - mug[i - 1]) / (mg[i] - mg[i - 1])
                    dmu_dot[inside] = mug[i - 1] + slope * (mc[inside] - mg[i - 1]) - self.gb
                    dmu_dm[inside] = slope
            c_eq = np.exp(np.clip(math.log(C_M) + (self.gb + dmu_dot - mu0) / self.kt, -700.0, 5.0))
            rate0 = self.size_rates(Ac, rc, c_eq, th_bar, C_M) if len(cohorts) else np.zeros(0)
            # facet-resolved dots (beyond the grid's sizes): dm/dt from the facet law, linear in C for the pool solve
            facet_on = self.P["kinetics"].get("growth_law", "size") == "facet" and bool(self.P["kinetics"].get("facets"))
            fmask = (mc > self.m_top + 1) if (facet_on and len(cohorts)) else np.zeros(len(cohorts), bool)
            shape_on = self.P["kinetics"].get("growth_law", "size") == "shape"
            smask = (mc > self.m_top + 1) if (shape_on and len(cohorts)) else np.zeros(len(cohorts), bool)
            if (shape_on or self.facet_area) and len(cohorts):
                res_d, res_r = self.resistances(Ac, rc, th_bar)
            fr, fv = {}, {}                                    # per facet dot: C -> (dm/dt, v_f), set up once a step
            for i in np.flatnonzero(fmask):
                fr[i] = self._facet_setup(cohorts[i], C_M, mu0, fam, th_dots)
                rate0[i], fv[i] = fr[i](C_M)
            vmax = float(np.abs(rate0 / np.maximum(mc, 1e-12)).max(initial=0.0))           # fastest relative change
            dt = outs[0] - t
            h1 = area1 = sigma_now = None
            if not self.facet_area:
                gam_now = {f: fam[f]["gamma"] for f in self.F}
                # the Wulff shape at fixed gamma is scale-free: the dot of one unit, scaled by m^(1/3) (h) and m^(2/3) (A)
                h1 = self.shape.wulff(gam_now, self.Om)
                area1, _v1, _s1 = self.shape(h1)
                sigma_now = sum(gam_now[f] * area1[f] for f in self.F)                # continuum: (2/3) sigma m^(-1/3)
            with np.errstate(over="ignore", divide="ignore"):
                g_, d_ = rate0 > 0, rate0 < 0
                if mg is not None and grid_dots:
                    # a grid-size dot that meets its own equilibrium within the allowed change stops there exactly
                    # (_stop_at_equilibrium): it does not limit the step
                    for i in np.nonzero(mc <= self.m_top + 1)[0]:
                        frac = 0.02 if g_[i] else 0.2
                        m_try = mc[i] * (1 + frac) if g_[i] else mc[i] * (1 - frac)
                        if rate0[i] != 0 and abs(self._stop_at_equilibrium(mc[i], m_try, mu0 - self.gb, mg, mug, sigma_now) - m_try) > 1e-12:
                            g_[i] = d_[i] = False
                if self.facet_area and len(cohorts):
                    # facet-area law (collision-limited facets: dots exchange with the pool within microseconds): a dot
                    # in a stable well (mu_dot rising with m, s < 0) whose own equilibrium lies within the allowed change
                    # of its size, |rate / s| <= 2 % (growing) or 20 % (dissolving) of m, cannot overshoot it under the
                    # backward-Euler size update below: it does not limit the step
                    s_w = -self.conductance(Ac, rc, th_bar) * c_eq * dmu_dm / self.kt
                    near = (s_w < 0) & (np.abs(rate0) <= np.where(rate0 > 0, 0.02, 0.2) * np.maximum(mc, 1.0) * -s_w)
                    g_ &= ~near
                    d_ &= ~near
                if g_.any():
                    # a growing dot changes m by < 2 % a step, measured against the larger of its own size and 5 % of
                    # the largest dot: newborn nuclei that grow fast (diffusion-limited) take larger relative steps,
                    # amounts stay exact, and the dots that hold the material keep the 2 %
                    m_ref = np.maximum(mc[g_], 0.05 * float(mc.max()))
                    dt = min(dt, float((0.02 * np.maximum(m_ref, 1.0) / rate0[g_]).min()))
                # a dissolving one by < 20 %, except the smallest (below 3 m_dissolve): subcritical nuclei that may
                # vanish within a step go back to the pool whole (exact), instead of forcing microsecond steps when
                # growth is fast (diffusion-limited) and the newborn nuclei sit at the critical size
                d_ &= mc > 3.0 * m_dis
                if d_.any():
                    dt = min(dt, float((0.2 * np.maximum(mc[d_], 1.0) / -rate0[d_]).min()))
                # each dot's own feedback at fixed C: s = d(dm/dt)/dm = -g c_eq (d mu_dot/dm) / kT, in the continuum
                # g c_eq (mu_dot - g_bulk) / (3 m kT). s < 0 (gamma < 0, or a grid size where mu rises with m: the dot's
                # size is a stable equilibrium) is stiff: there
                # the size update is linearised backward Euler, dm = dt rate / (1 - s dt), stable for any dt. s > 0 (the
                # Ostwald instability) stays explicit: with the pool implicit, one dot in a closed pool relaxes as
                # x' = x (1 + s dt) / (1 + N g dt), stable whenever the physics is (N g > s), and the 2 % / 20 % limits
                # bound how far ripening moves in one step (checked: LSW and Wagner rate constants and distributions)
                s_fb = np.minimum(-self.conductance(Ac, rc, th_bar) * c_eq * dmu_dm / self.kt, 0.0)
                s_fb[fmask] = 0.0
                for i, v0 in fv.items():                                    # a facet moves by < 2 % of its h a step
                    hs = [cohorts[i]["h"][f] / abs(v0[f]) for f in self.F if cohorts[i]["area"][f] > 0 and v0[f]]
                    if hs:
                        dt = min(dt, 0.02 * min(hs))
            if Q > 0:
                # MX_q < 2 %; the precursor itself decays exactly over the step (below), so a fast conversion of a
                # slowly injected precursor does not limit the step
                dt = min(dt, 0.02 * C_X / Q)
            dt_rel = self.P["material"].get("dt_rel_max", 0.01)
            if dt_rel and t > 0:
                dt = min(dt, float(dt_rel) * t)             # resolution in log time (material.dt_rel_max)
            if dt_pool:
                dt = min(dt, dt_pool)                       # the last pool-limited step, doubling back each step
            nxt = [b for b in t_breaks if b > t + 1e-12]
            if nxt:
                dt = min(dt, nxt[0] - t)                    # land on the temperature and injection breakpoints
            dTdt = abs(self.temperature_at(t + 1e-6) - self.temperature_at(t)) / 1e-6
            if dTdt > 0:
                dt = min(dt, 1.0 / dTdt)                    # a ramp moves T by < 1 K a step
            injecting = bool(inj) and inj["t0"] - 1e-12 <= t < inj["t0"] + inj["duration"] - 1e-12
            if injecting:
                dt = min(dt, 0.02 * V / inj["rate"])        # the volume grows by < 2 % a step
            q_inst = Q
            pool_limited = False
            while True:
                # a dot gives back at most what it holds within the step (the update below clips m at 1e-6): the pool
                # solve takes the same clipped change, or a dot dissolving faster than 1/dt drives the solve far above
                # what is actually released
                r_min = (1e-6 - mc) / dt if dt > 0 else np.full(len(mc), -np.inf)
                # conversion over the step at the step's MX_q: P -> P e^(-k C_X dt), Q the mean rate
                kx = self.k_conv * C_X
                Q = P_ * -math.expm1(-kx * dt) / dt if (q_inst > 0 and dt > 0) else 0.0

                def balance(lnC):
                    mu = mu0 + self.kt * (lnC - math.log(C_M))
                    J_, *_ = self.nucleation(mu, mu_mx)
                    if len(cohorts):
                        rr = np.maximum(self.size_rates(Ac, rc, c_eq, th_bar, math.exp(lnC)) / (1.0 - s_fb * dt), r_min)
                        for i, f_ in fr.items():
                            rr[i] = f_(math.exp(lnC))[0]
                        U = float((Nc * rr).sum())
                    else:
                        U = 0.0
                    return math.exp(lnC) - C_M - dt * (Q - self._m_birth * J_ - U)
                lo, hi = math.log(1e-30), math.log(C_M + Q * dt + 1e-12)
                b_lo, b_hi = balance(lo), balance(hi)
                # dissolving dots can push C above C_M + Q dt: widen the bracket (not for facet-law dots, whose capped
                # inward moves do not follow the solve's rates: their steps rely on the old bracket's limit)
                for _ in range(0 if fr else 200):
                    if b_hi >= 0:
                        break
                    hi += math.log(2.0)
                    b_hi = balance(hi)
                lnC = brentq(balance, lo, hi, xtol=1e-12, rtol=1e-12) if b_lo < 0 < b_hi else (hi if b_hi < 0 else lo)
                mu = mu0 + self.kt * (lnC - math.log(C_M))
                J, barrier, k_star, _per_m = self.nucleation(mu, mu_mx)
                m_b = self._m_birth                                                    # units of each new nucleus
                if k_star is None:                                                     # classical: continuum coverage
                    k_star = self.k_of(m_b, fam)
                rate = (np.maximum(self.size_rates(Ac, rc, c_eq, th_bar, math.exp(lnC)) / (1.0 - s_fb * dt), r_min)
                        if len(cohorts) else np.zeros(0))
                # Z-type coverage on the dots moves toward the Langmuir value of the free MX_q, but only as far as the
                # MX_q left in solution allows (the dots' surfaces are a reservoir of MX_q too)
                lam = 1.0
                vs = {}                                                    # facet velocities at the solved C
                while True:
                    th_eff = {f: th_dots[f] + lam * (fam[f]["theta"] - th_dots[f]) for f in self.F}
                    k_area = {f: 0.0 if self.sphere else th_eff[f] * fam[f]["Gamma_max"] + self.F[f]["x_X"] / self.F[f]["site_area"]
                              for f in self.F}
                    new, up_m, up_k = [], 0.0, 0.0
                    for i, co in enumerate(cohorts):
                        if self.facet_area:
                            # facet-area law: m from the channels' rate; each exposed plane moves at Omega_m j_f (a
                            # facet with j_f < 0 dissolves), hidden planes are held; areas and edges from the dot
                            # built on the fly (or the model's polyhedron), k = sum_f A_f rho_f theta_f
                            m_new = max(co["m"] + rate[i] * dt, 1e-6)
                            if i not in vs:
                                vs[i] = self.facet_speeds(co, math.exp(lnC), C_M, mu0, res_d[i])
                            h_new = self._fa_step(co, m_new, vs[i])
                            area, edges, built = self.geometry(h_new, m_new)
                            k_new = sum(k_area[f] * area[f] for f in self.F)
                            up_m += co["N"] * (m_new - co["m"])
                            up_k += co["N"] * (k_new - co["k"])
                            new.append({"h": h_new, "m": m_new, "k": k_new, "area": area, "edges": edges,
                                        "build": built and {"dc": built["dc"], "m": built["m"],
                                                            "h_err_A": built["h_err_A"]}})
                            continue
                        if fmask[i]:                                       # facet law: each facet moves by v_f dt
                            if i not in vs:
                                vs[i] = fr[i](math.exp(lnC))[1]
                            # hidden facets move too: a fast one stays hidden, an island one (no area, no nuclei)
                            # stands still and comes back as the dot grows past it
                            # (the step limit used the velocities at the old C: at the solved C a facet may dissolve
                            # faster, so the inward move is capped at 10 % of h; m and the uptake follow from V(h))
                            h_new = {f: co["h"][f] + max(vs[i][f] * dt, -0.1 * co["h"][f]) for f in self.F}
                            area, vol, support = self.shape(h_new)
                            h_new = {f: min(h_new[f], support[f] + 1.0) for f in self.F}
                            m_new = vol / self.Om
                            k_new = sum(k_area[f] * area[f] for f in self.F)
                            up_m += co["N"] * (m_new - co["m"])
                            up_k += co["N"] * (k_new - co["k"])
                            new.append({"h": h_new, "m": m_new, "k": k_new, "area": area})
                            continue
                        m_new = max(co["m"] + rate[i] * dt, 1e-6)
                        if mg is not None and grid_dots and min(co["m"], m_new) <= self.m_top + 1:
                            m_new = self._stop_at_equilibrium(co["m"], m_new, mu - self.gb, mg, mug, sigma_now)
                        if smask[i] and m_new > self.m_top + 1:
                            # shape rule: the size law's growth shared among the facets by their relative speeds,
                            # at the monomer level of the dot's surface (diffusion and reaction resistances in series)
                            if i not in vs and abs(m_new - co["m"]) > 1e-10 * co["m"]:    # frozen: shape_step keeps h
                                C_now = math.exp(lnC)
                                c_s = C_now - (C_now - c_eq[i]) * res_d[i] / (res_d[i] + res_r[i])
                                dmu_s = mu0 - self.gb + self.kt * math.log(max(c_s, 1e-300) / C_M)
                                vs[i] = (self.facet_speeds(co, C_now, C_M, mu0, res_d[i]) if self.facet_law
                                         else self.shape_speeds(co, dmu_s, fam, th_dots))
                            h_new, area = self.shape_step(co, m_new, vs.get(i))
                        else:
                            s = m_new ** (1.0 / 3.0)                                   # shape: Wulff at gamma(C)
                            h_new = {f: h1[f] * s for f in self.F}
                            area = {f: area1[f] * s * s for f in self.F}
                        k_new = sum(k_area[f] * area[f] for f in self.F)               # as self.dot at theta_eff
                        up_m += co["N"] * (m_new - co["m"])
                        up_k += co["N"] * (k_new - co["k"])
                        new.append({"h": h_new, "m": m_new, "k": k_new, "area": area})
                    C_X_new = C_X + (-Q - k_star * J) * dt - up_k
                    if C_X_new > 0 or lam < 1e-6:
                        break
                    lam *= 0.25
                C_M_new = C_M + (Q - m_b * J) * dt - up_m
                ok = C_M_new > 0 and C_X_new > 0
                if ok and dmu_max > 0 and dt > 1e-9:
                    pn = self.pool(C_M_new, C_X_new)
                    if max(abs(pn["mu_MA"] - mu0), abs(pn["mu_MX"] - mu_mx)) > dmu_max * self.kt:
                        ok, pool_limited = False, True
                if ok or dt < 1e-9:
                    break
                dt *= 0.5
            dt_pool = dt if pool_limited else (2.0 * dt_pool if dt_pool and 2.0 * dt_pool < t else None)
            th_dots = th_eff
            for co, n in zip(cohorts, new):
                co.update(n)
            P_ -= Q * dt
            C_M, C_X = max(C_M_new, 1e-30), max(C_X_new, 1e-30)
            pending += J * dt
            pending_m += m_b * J * dt
            pending_k += k_star * J * dt
            for co in [co for co in cohorts if co["m"] < m_dis]:                         # dissolved dots give back
                C_M += co["N"] * co["m"]
                C_X += co["N"] * co["k"]
                cohorts.remove(co)
            n_conv += Q * dt * V
            n_nuc += J * dt * V
            track.append((t + dt, n_nuc, J * V))
            if injecting:
                # syringe: the flask's volume grows by rate dt; everything already there is diluted by V / V_new and
                # the precursor gains c_syringe rate dt moles (amounts conserved exactly)
                dV = inj["rate"] * dt
                f_d = V / (V + dV)
                P_, C_M, C_X = P_ * f_d + inj["c"] * dV / (V + dV), C_M * f_d, C_X * f_d
                pending, pending_m, pending_k = pending * f_d, pending_m * f_d, pending_k * f_d
                for co in cohorts:
                    co["N"] *= f_d
                n_inj += inj["c"] * dV
                V += dV
            self.steps.append((t, dt, vmax, Q, C_M, C_X_new, J))
            t += dt
            self.pools.append((t, self.T, C_M, C_X))
            at_out = t >= outs[0] - 1e-12
            n_dots = sum(co["N"] for co in cohorts)
            if at_out or t - t_last_birth >= birth_frac * t or pending >= birth_share * (n_dots + pending):
                # this interval's nuclei become a cohort at every output, every birth_interval of t, and whenever they
                # reach birth_share of all dots (a burst is resolved into many cohorts; results do not depend on n_out)
                if pending > self.P["material"].get("N_min", 1e-15):
                    t_last_birth = t
                    m_new = pending_m / pending                          # the interval's mean birth size
                    h0 = (self.fa_start_h(m_new) if self.facet_area
                          else self.shape.wulff({f: fam[f]["gamma"] for f in self.F}, m_new * self.Om))
                    # the nuclei took m0 MA and k*(t) MX_q each from the pool; a Wulff dot of that volume may hold a
                    # slightly different m and k: the difference goes back to (or comes from) the pool. When the MX_q
                    # left cannot cover the new dots at the dots' coverage, the shared coverage of all dots is lowered
                    # (their MX_q returns) until the balance keeps 10 % of the MX_q left
                    def at_scale(s):
                        f_s = {f: {**fam[f], "theta": s * th_dots[f]} for f in self.F}
                        ks = [self.dot(co["h"], f_s)[1] for co in cohorts]
                        m_s, k_s, area_s, sup_s = self.dot(h0, f_s)
                        dk = sum(co["N"] * (kk - co["k"]) for co, kk in zip(cohorts, ks)) + pending * k_s - pending_k
                        return dk, ks, (m_s, k_s, area_s, sup_s)
                    s = 1.0
                    dk, ks, new_dot = at_scale(s)
                    if dk > 0.9 * C_X:
                        lo, hi = 0.0, 1.0
                        for _ in range(40):
                            mid = 0.5 * (lo + hi)
                            lo, hi = (mid, hi) if at_scale(mid)[0] <= 0.9 * C_X else (lo, mid)
                        s = lo
                        dk, ks, new_dot = at_scale(s)
                        th_dots = {f: s * th_dots[f] for f in self.F}
                        for co, kk in zip(cohorts, ks):
                            co["k"] = kk
                        C_X -= dk - (pending * new_dot[1] - pending_k)               # the old dots' share
                    m_, k_, area, support = new_dot
                    if dk > C_X:                    # even bare, the new dots' X-type MX_q exceeds what is left
                        self.mx_shortfall = getattr(self, "mx_shortfall", 0.0) + dk - C_X
                    C_M += pending_m - pending * m_
                    C_X = max(C_X + pending_k - pending * k_, 1e-30)   # a shortfall is counted in mx_shortfall
                    born = {"t_b": t, "N": pending, "h": {f: min(h0[f], support[f]) for f in self.F},
                            "m": m_, "k": k_, "area": area}
                    if self.facet_area:                          # hidden planes are held where the start shape has them
                        _a, edges, built = self.geometry(h0, m_)
                        born.update(h=dict(h0), edges=edges,
                                    build=built and {"dc": built["dc"], "m": built["m"], "h_err_A": built["h_err_A"]})
                    cohorts.append(born)
                    pending = pending_m = pending_k = 0.0
                    self._merge(cohorts)
            if at_out:
                N = sum(co["N"] for co in cohorts)
                mm = np.array([co["m"] for co in cohorts]) if cohorts else np.zeros(0)
                w = np.array([co["N"] for co in cohorts]) if cohorts else np.zeros(0)
                mean = float((w * mm).sum() / N) if N > 0 else 0.0
                std = float(np.sqrt((w * (mm - mean) ** 2).sum() / N)) if N > 0 else 0.0
                pl_out = self.pool(C_M, C_X) if self.c_sat is not None else pl  # the solid at the output time
                row = {"t": t, "P": P_, "C_M": C_M, "C_X": C_X, "c_MX_free": pl_out["c_MX_free"],
                       "c_MX_solid": pl_out.get("c_MX_solid", 0.0), "mu_MA": mu,
                       "mu_MX": mu_mx, "S_bulk": math.exp(min((mu - self.gb) / self.kt, 700)),
                       "T": self.T, "V_L": V, "n_inj_mol": n_inj, "n_conv_mol": n_conv, "nucleated_M": n_nuc / V,
                       "k_eff": self.mean_keff(cohorts, th_dots),
                       "barrier": barrier, "J": J, "J0": J * math.exp(min(barrier / self.kt, 700)), "k_star": k_star,
                       "m_birth": m_b,
                       "N_dots": N, "m_mean": mean, "m_std": std,
                       "theta": dict(th_dots), "theta_eq": {f: fam[f]["theta"] for f in self.F},
                       "gamma": {f: fam[f]["gamma"] for f in self.F}}
                series.append(row)
                snaps.append({"t": t, "cohorts": [{"t_b": co["t_b"], "N": co["N"], "m": co["m"], "k": co["k"],
                                                   "h": co["h"], "area": co["area"],
                                                   **({"edges": co.get("edges"), "build": co.get("build")}
                                                      if self.facet_area else {})} for co in cohorts]})
                if log:
                    log(f"[growth] t = {t:9.3g} s  P = {P_:.3g} M  C_M = {C_M:.3g} M  barrier {barrier:.3f} eV  "
                        f"J = {J:.3g} M/s  dots {N:.3g} M  <m> = {mean:.1f}")
                outs.pop(0)
        bg = self.__dict__.get("_bgeo")
        geo = ({"source": self._fa().get("geometry", "builder"), "builds": bg.builds if bg else 0,
                "build_seconds": round(bg.seconds, 2) if bg else 0.0, "dots": len(bg.cache) if bg else 0,
                "fallbacks": getattr(self, "geometry_fallbacks", 0), "last_error": getattr(self, "geometry_error", None)}
               if self.facet_area else None)
        return {"params": self.P, "match_m0": self.match, "handover": self.handover, "m0": self.m0, "series": series,
                "geometry": geo,
                "snapshots": snaps, "seeded": seeded, "V0_L": V0, "nucleation_window": nucleation_window(track),
                "sphere": getattr(self, "sphere_used", None),
                "survivor_window": survivor_window(cohorts)}


def merged_birth(a: dict, b: dict) -> float:
    """Birth time of two merged cohorts: their number-weighted mean in ln t (arithmetic when one was seeded at 0)."""
    N = a["N"] + b["N"]
    if a["t_b"] > 0 and b["t_b"] > 0:
        return math.exp((a["N"] * math.log(a["t_b"]) + b["N"] * math.log(b["t_b"])) / N)
    return (a["N"] * a["t_b"] + b["N"] * b["t_b"]) / N


def survivor_window(cohorts) -> dict | None:
    """When the dots left at the end were born: the birth times (cohorts' t_b, number-weighted) by which 10 %, 50 %
    and 90 % of them had formed, and the fraction of the end's dots this covers. None without dots."""
    co = sorted((c for c in cohorts if c["N"] > 0), key=lambda c: c["t_b"])
    if not co:
        return None
    t = np.array([c["t_b"] for c in co])
    w = np.cumsum([c["N"] for c in co])
    frac = w / w[-1]
    out = {}
    for q in (0.1, 0.5, 0.9):
        out[f"t{int(100 * q)}"] = float(t[min(int(np.searchsorted(frac, q)), len(t) - 1)])
    return {**out, "dots_M": float(w[-1])}


def nucleation_window(track) -> dict | None:
    """When the nuclei formed, from the per-step record (t, moles of nuclei so far, J V): the times by which 10 %,
    50 % and 90 % of all nuclei had formed and the time of the largest nucleation rate. None without nuclei."""
    t = np.array([x[0] for x in track])
    n = np.array([x[1] for x in track])
    rate = np.array([x[2] for x in track])
    if len(t) < 2 or n[-1] <= 0:
        return None
    frac = n / n[-1]
    at = {f"t{int(100 * q)}": float(np.interp(q, frac, t)) for q in (0.1, 0.5, 0.9)}
    return {**at, "t_peak": float(t[int(np.argmax(rate))]), "nuclei_mol": float(n[-1])}


def conservation(res) -> dict:
    """Total MA (precursor + pool + dots, less what the syringe added) and total MX_q (left + converted + dots) at
    each output, in moles per litre of the starting flask: both constant (0 here)."""
    c = res["params"]["conditions"]
    s0 = res.get("seeded") or {"m": 0.0, "k": 0.0}                                       # dots present at t = 0
    V0 = res.get("V0_L", 1.0)
    out = []
    for row, snap in zip(res["series"], res["snapshots"]):
        dots_m = sum(co["N"] * co["m"] for co in snap["cohorts"])
        dots_k = sum(co["N"] * co["k"] for co in snap["cohorts"])
        V = row.get("V_L", V0)
        n_conv = row.get("n_conv_mol", (c["c_P0"] - row["P"]) * V0)
        out.append({"t": row["t"],
                    "MA": (V * (row["P"] + row["C_M"] + dots_m) - row.get("n_inj_mol", 0.0)) / V0
                    - c["c_P0"] - c.get("c_M0", 1e-12) - s0["m"],
                    "MX": (V * (row["C_X"] + dots_k) + n_conv) / V0 - c["c_MX0"] - s0["k"]})
    return out


def facet_area_from_roles(roles: dict, start_h: dict, builder: dict, rate_form: str = "first_order") -> dict:
    """kinetics.facet_area block from the per-role block (kinetics.roles: envs {role: {m: [[n, dmu, bind], ...]}},
    neb {role: {dir, disp}}), computed values only:
      facets   per family: dG_bind_eV = Cd-weighted plain mean of the terrace (facet:f) binding over every computed
               size (gas phase, 1 M); dmu_eV = per size, the Cd-weighted plain mean of the terrace Delta mu_e (a_X = 1);
               neb = the facet NEBs
      edges    per edge type: dG_bind_eV pooled the same way over the edge:e environments; neb = the edge NEBs (null
               until computed: the edge is then left out of the step band)
    Vertices and the edges' Delta mu are not used (Delta mu from the terraces only)."""
    def pooled(by):
        n = b = 0.0
        for envs in by.values():
            for x in envs:
                n += x[0]
                b += x[0] * x[2]
        return round(b / n, 4)

    facets, edges = {}, {}
    for role, by in roles["envs"].items():
        kind, fams = role.split(":")
        if kind == "facet":
            dmu = {int(m): round(sum(x[0] * x[1] for x in envs) / sum(x[0] for x in envs), 4) for m, envs in by.items()}
            facets[fams] = {"dG_bind_eV": pooled(by), "dmu_eV": dict(sorted(dmu.items())),
                            "neb": (roles.get("neb") or {}).get(role)}
        elif kind == "edge":
            edges[fams] = {"dG_bind_eV": pooled(by), "neb": (roles.get("neb") or {}).get(role)}
    return {"T_K": roles.get("T_K"), "rate_form": rate_form, "geometry": "builder", "geometry_step_A": 0.5,
            "start_h_A": dict(start_h), "builder": dict(builder), "facets": facets, "edges": edges}


def builder_recipe(res, t: float, preset: dict, which: str = "mode", cohort: int | None = None) -> dict:
    """QD_builder recipe of the dot at time t: the cohort with the most dots ('mode'), the one closest to the
    number-weighted mean size ('mean'), the one holding the most material ('mass'), or cohort number `cohort` of
    the snapshot. preset: charges, passivation, symmetry of the material config."""
    snap = min(res["snapshots"], key=lambda s: abs(s["t"] - t))
    if not snap["cohorts"]:
        raise ValueError(f"no dots at t = {snap['t']:.3g} s")
    cos = snap["cohorts"]
    if cohort is not None:
        co = cos[cohort]
    elif which == "mass":
        co = max(cos, key=lambda c: c["N"] * c["m"])
    elif which == "mean":
        N = sum(c["N"] for c in cos)
        mean = sum(c["N"] * c["m"] for c in cos) / N
        co = min(cos, key=lambda c: abs(c["m"] - mean))
    else:
        co = max(cos, key=lambda c: c["N"])
    row = min(res["series"], key=lambda r: abs(r["t"] - snap["t"]))
    if res["params"]["kinetics"].get("size_law") == "facet_area":
        # the dot built on the fly for this cohort (its construction distances), then the model's Z-type CdCl2:
        # target k, split over the families as sum_f A_f rho_f theta_f (the terms of k)
        b = co.get("build")
        if not b:
            raise ValueError("no built dot for this cohort: kinetics.facet_area.geometry must be 'builder'")
        spec = res["params"]["kinetics"]["facet_area"]["builder"]
        a = (4.0 * float(res["params"]["material"]["Omega_m"])) ** (1.0 / 3.0)
        z = {f: row["theta"][f] / (a * a / 2.0 if f == "100" else SQ3 * a * a / 4.0) * co["area"].get(f, 0.0)
             for f in FA_FAMILIES}
        recipe = {"charges": copy.deepcopy(spec["charges"]), "passivation": copy.deepcopy(spec["passivation"]),
                  "symmetry": copy.deepcopy(spec.get("symmetry", {})),
                  "construction_origin": {"center_on_species": spec["centre"]},
                  "facet_distances": {f: round(float(v), 4) for f, v in b["dc"].items()},
                  "post_treatment": {"z_type_target": {"target_total": int(round(co["k"])),
                                                       "per_family": {f: round(float(v), 4) for f, v in z.items()
                                                                      if v > 0}}}}
        import os
        return {"t": snap["t"], "cohort_born": co["t_b"], "m": co["m"], "k": co["k"], "recipe": recipe,
                "cif": os.path.expandvars(spec["cif"]), "positive_q_mode": spec.get("positive_q_mode", "add")}
    F = res["params"]["families"]
    z = {f: row["theta"][f] * gamma_max(F[f]) * co["area"].get(f, 0.0) for f in F}
    recipe = {k: copy.deepcopy(v) for k, v in preset.items() if k not in ("facets", "size_unit_cells")}
    recipe["facet_distances"] = {f: round(float(co["h"][f]), 3) for f in F}
    recipe["target"] = {"m": int(round(co["m"])), "w": 1.0}
    recipe.setdefault("post_treatment", {})["z_type_target"] = {
        "target_total": int(round(co["k"])), "per_family": {f: round(float(v), 4) for f, v in z.items() if v > 0}}
    return {"t": snap["t"], "cohort_born": co["t_b"], "m": co["m"], "k": co["k"], "recipe": recipe}


def build_frames(res, out_dir, which: str = "mode", even: bool = True, log=print) -> list:
    """The growth film of a facet-area run: at every output time with dots, the chosen cohort ('mode', 'mean', 'mass')
    built by QD_builder: the charge-balanced dot at the cohort's construction distances (the dot built on the fly for
    its geometry), then z_type_target adds or displaces neutral CdCl2 units to the model's k, split over the families
    as A_f rho_f theta_f and placed evenly over the facets, strongest-bound family first (ztype_even). Writes
    frame_NNN.xyz and frames.json (model against built: m, k, lattice-site areas, units per family); returns the
    records."""
    import tempfile
    from builder.scripts.generate_library import _build
    from .ztype_even import even_placement
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    B = res["params"]["kinetics"]["facet_area"]["facets"]
    anion = res["params"]["kinetics"]["facet_area"]["builder"].get("anion", "Se")
    frames = []
    for i, snap in enumerate(res["snapshots"]):
        if not snap["cohorts"]:
            continue
        r = builder_recipe(res, snap["t"], {}, which)
        zt = r["recipe"]["post_treatment"]["z_type_target"]
        zt["per_family"] = dict(sorted(zt["per_family"].items(), key=lambda kv: float(B[kv[0]]["dG_bind_eV"])))
        rec = {"frame": i, "t": r["t"], "cohort_born": r["cohort_born"],
               "model": {"m": r["m"], "k": r["k"]}, "recipe": r["recipe"]}
        with tempfile.TemporaryDirectory() as tmp, even_placement(even):
            o = _build(r["cif"], r["recipe"], r["positive_q_mode"], Path(tmp) / "b")
        if not o or "error" in o:
            rec["error"] = (o or {}).get("error", "no output")
            frames.append(rec)
            log(f"[film] t = {r['t']:.4g} s: build failed: {rec['error']}")
            continue
        syms = list(o["symbols"])
        pts = np.asarray(o["pts"], float)
        man = o.get("manifest") or {}
        lf = (man.get("facet_geometry") or {}).get("lattice_fitted") or {}
        led = man.get("z_type_target_ledger") or {}
        cation = next(e for e, q in r["recipe"]["charges"].items() if q > 0)
        rec["built"] = {"m": syms.count(anion), "k": syms.count(cation) - syms.count(anion), "n_atoms": len(syms),
                        "composition": {e: syms.count(e) for e in sorted(set(syms))},
                        "areas_A2": lf.get("areas_A2"), "edges": lf.get("edges")}
        rec["z_type"] = {k: led.get(k) for k in ("target_total", "k_X", "direction", "units", "achieved_total",
                                                  "shortfall", "error")}
        rec["z_type"]["per_family"] = {f: {k: v.get(k) for k in ("requested", "added", "removed", "even_count",
                                                                  "facets", "per_facet") if k in v}
                                       for f, v in (led.get("per_family") or {}).items()}
        xyz = out / f"frame_{i:03d}.xyz"
        with open(xyz, "w") as fh:
            fh.write(f"{len(syms)}\nt = {r['t']:.4g} s, model m = {r['m']:.1f}, k = {r['k']:.1f}\n")
            for e, p in zip(syms, pts):
                fh.write(f"{e:2s} {p[0]:12.6f} {p[1]:12.6f} {p[2]:12.6f}\n")
        rec["xyz"] = xyz.name
        frames.append(rec)
        log(f"[film] t = {r['t']:.4g} s: model m {r['m']:.1f} k {r['k']:.1f} -> built m {rec['built']['m']} "
            f"k {rec['built']['k']} ({rec['z_type'].get('direction')} {rec['z_type'].get('units')})")
    (out / "frames.json").write_text(json.dumps(frames, indent=1, default=float))
    return frames


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python -m qd_workflows.growth")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("params", help="write the placeholder parameters")
    p.add_argument("--out", required=True)
    p = sub.add_parser("run", help="integrate nucleation and growth")
    p.add_argument("--thermo", required=True, help="thermo_merged.json of the nucleation workflow")
    p.add_argument("--params", help="parameters (default: the placeholders)")
    p.add_argument("--out", required=True)
    p.add_argument("--recipe-at", type=float, help="also write the builder recipe of the dot at this time (s)")
    p.add_argument("--preset", help="material config whose preset (charges, passivation) goes into the recipe")
    p.add_argument("--films", help="facet-area law: directory for the growth film (the most numerous cohort built by "
                                   "QD_builder at every output time, with the model's CdCl2 placed evenly)")
    a = ap.parse_args(argv)
    if a.cmd == "params":
        Path(a.out).write_text("# Growth-model parameters: every value is a PLACEHOLDER until the shape-library fits replace it.\n"
                               + yaml.safe_dump(PARAMS, sort_keys=False))
        return 0
    params = yaml.safe_load(Path(a.params).read_text()) if a.params else PARAMS
    res = Growth(json.loads(Path(a.thermo).read_text()), params).run(log=lambda x: print(x, flush=True))
    res["conservation"] = conservation(res)
    if a.recipe_at is not None:
        preset = yaml.safe_load(Path(a.preset).read_text())["preset"] if a.preset else {}
        res["recipe"] = builder_recipe(res, a.recipe_at, preset)
    if a.films:
        res["films"] = str(Path(a.films).resolve())
        build_frames(res, a.films)
    Path(a.out).write_text(json.dumps(res, default=float))
    return 0


if __name__ == "__main__":
    sys.exit(main())
