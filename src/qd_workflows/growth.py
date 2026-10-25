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
    solvation   implicit only (grid_solvation): every solution species, monomer or cluster, and every facet
                (a_solv per area). EXPLICIT SOLVENT AND LIGANDS (oleic acid, amines, phosphines) ARE NOT INCLUDED:
                they are the next step and change the monomers, every (m, k) and the facets.
    hand-over   m0 = "auto": the size where the grid's two-step chemical potential meets the continuum (_handover)
    dot mu      mu_dot = g_bulk(T) + (2/3) Omega_m sum_f gamma_f A_f / V   (Gibbs-Thomson of the Wulff dot: one
                chemical potential per dot, = g_bulk + 2 gamma_f Omega_m / h_f on every facet of the Wulff shape)
                g_bulk(T) = E_bulk + 3 n_at [hw_E/2 + kT ln(1 - e^-hw_E/kT)]  (free energy, like the grid's G)
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
                   "t_start": 1e-3, "n_out": 60},
    "material": {
        "Omega_m": A_ZB ** 3 / 4,               # A^3 per CdSe (4 formula units per cubic cell)
        "g_bulk": -5.4903,                      # eV per CdSe: MACE bulk zinc-blende energy at 0 K (shapes fit)
        "omega_E_cm1": 150.0,                   # Einstein phonon of bulk CdSe (Debye ~210 K) for g_bulk(T) (PLACEHOLDER)
        "atoms_per_unit": 2,
        "m0": "auto",                           # nucleus size: "auto" = where grid and continuum agree (_handover)
        "m0_min": 10,                           # smallest nucleus the automatic hand-over may pick
        "m_dissolve": 4,                        # a dot below this many CdSe units dissolves back into the pool
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
                 "D": 1e-5,                     # cm^2/s, monomer diffusion coefficient in the solvent (PLACEHOLDER)
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
              "Ea_conv": "eV", "t_end": "s"},
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
    kc = math.exp(min(-float(p["dG_bind"]) / kt, 700.0)) * c
    theta = theta_fg(kc, w / kt)
    return {"theta": theta, "gamma": float(p["gamma0"]) * 1e-3 + gmax * (kt * math.log1p(-theta) - 0.5 * w * theta ** 2),
            "Gamma_max": gmax}


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

    def wulff(self, gamma: dict, volume: float) -> dict:
        """Facet distances of the Wulff shape (h_f proportional to gamma_f) with the given volume."""
        h = {f: max(g, 1e-6) for f, g in gamma.items()}
        _a, vol, _s = self(h)
        s = (volume / vol) ** (1.0 / 3.0)
        return {f: x * s for f, x in h.items()}


# ---------------- nucleation flux ----------------
def nucleation_flux(om, ms, m0: int, kt: float, omega_m: float, P: dict, barrier: float | None = None) -> float:
    """Steady-state Becker-Doring flux (M/s) of clusters through the grid up to the nucleus m0:

        J = [ sum_{m=1}^{m0-1} 1 / (f_m c_eq(m)) ]^-1,   c_eq(m) = c0 exp(-Omega(m)/kT)   (1 M standard state),
        f_m = 4 pi (D_1 + D_m)(r_1 + r_m) N_A C_MA exp(-E_att/kT),   D_m = D_1 m^(-1/3),   r_m = (3 m Omega_m / 4 pi)^(1/3)

    Smoluchowski attachment of monomers to an m-cluster (relative diffusion, contact radius), with an optional
    incorporation barrier E_att. The monomer supply is c_eq(1) = exp(-Omega(1)/kT): the grid's m = 1 row and the pool
    share mu_MA, so Omega(1) = -kT ln(C_MA/c0). The largest term (the rarest cluster) dominates, J ~ f* c_eq(m*); the
    effective prefactor J exp(barrier/kT) is reported as J0. P['nucleation']['flux'] = 'constant' keeps J0 exp(-barrier/kT)."""
    nuc = P["nucleation"]
    if nuc.get("flux", "becker_doring") == "constant":
        return float(nuc["J0"]) * math.exp(-(barrier or 0.0) / kt)
    om, ms = np.asarray(om, float), np.asarray(ms)
    sel = ms < m0
    if not sel.any() or ms[0] != 1:
        return 0.0
    c1 = math.exp(-om[0] / kt)                                                    # M, the monomer supply
    d1 = float(P["kinetics"].get("D", 1e-5)) * 1e16                               # A^2/s
    m = ms[sel].astype(float)
    r = (3.0 * m * omega_m / (4.0 * math.pi)) ** (1.0 / 3.0)                      # A
    f = (4.0 * math.pi * d1 * (1.0 + m ** (-1.0 / 3.0)) * (r[0] + r) * 6.02214076e-4 * c1
         * math.exp(-float(nuc.get("attach_barrier_eV", 0.0)) / kt))              # 1/s; 6.022e-4 molecules/A^3 per M
    x = om[sel] / kt - np.log(f)                                                  # ln of each resistance 1/(f c_eq)
    xm = x.max()
    return float(math.exp(-(xm + math.log(np.exp(x - xm).sum()))))


# ---------------- the model ----------------
class Growth:
    def __init__(self, thermo: dict, params: dict, normals: dict | None = None):
        self.P = params
        c = params["conditions"]
        self.T = float(c["T"])
        self.kt = KB_EV * self.T
        self.solv = params.get("solvation") or {"eps": 1.0}
        self.th = thermo_at(thermo, self.T, self.solv)
        m0 = params["material"]["m0"]
        self.m0 = max(m for m, _k in self.th["G"]) if m0 == "auto" else int(m0)      # "auto": set by handover()
        self.Om = float(params["material"]["Omega_m"])
        self.F = params["families"]
        self.shape = Shape({f: (normals or ZB_NORMALS)[f] for f in self.F})
        self.k_conv = k_conv_at(self.T, float(c["k_conv_500K"]), float(c["Ea_conv"]))
        self.gb = self.g_bulk()
        self.dg_th = binding_thermal(thermo, self.T)                    # thermal part of binding one MX_q at T
        for p in self.F.values():
            if "eps0" in p:                                              # binding free energy at T from the energy
                p["dG_bind"] = -float(p["eps0"]) + self.dg_th
        # in solution: a unit binding to a facet gives up the solvation of the free MX_q (the facet stays solvated),
        # and every facet gains the solvation free energy of a neutral nanocrystal surface per area
        eps = float(self.solv.get("eps", 1.0))
        self.solv_area = per_area(eps, self.solv) if eps > 1.0 else 0.0
        if eps > 1.0:
            for p in self.F.values():
                p["dG_bind"] = float(p["dG_bind"]) - self.th["dG_solv_mx"]
                p["gamma0"] = float(p["gamma0"]) + 1000.0 * self.solv_area
        classical = params["nucleation"].get("omega", "grid") == "continuum"
        self.handover = self._handover() if m0 == "auto" and not classical else None
        self.match = self._match_m0() if params["material"].get("match_m0", False) else None

    def _handover(self):
        """Nucleus size m0 where the grid hands over to the continuum. The grid's step chemical potential
        mu_MA + Omega(m) - Omega(m - 1) oscillates between even and odd m (magic clusters), so it is averaged over two
        steps, mu2(m) = mu_MA + [Omega(m) - Omega(m - 2)] / 2, and compared with the continuum's
        g_bulk(T) + 2 gamma_f Omega_m / h_f of a Wulff dot of m units (exposed facet, gamma at the pool's coverage),
        both at the pool of half conversion. m0 is the size, from m0_min (default: half the largest row) up, where the
        two agree best (a tie goes to the larger m); the fitted gamma is not rescaled, the mismatch left is reported."""
        c = self.P["conditions"]
        try:
            pl, fam = self.chem(0.5 * c["c_P0"], c["c_MX0"] - 0.5 * c["c_P0"])
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
        """Bulk free energy per MA unit at T: the 0 K energy plus an Einstein estimate of the bulk phonons
        (3 x atoms-per-unit modes at omega_E), PLACEHOLDER until a bulk phonon calculation replaces it."""
        mat = self.P["material"]
        hw = float(mat.get("omega_E_cm1", 0.0)) * 1.239841984e-4                       # eV
        if hw <= 0:
            return float(mat["g_bulk"])
        modes = 3 * int(mat.get("atoms_per_unit", 2))
        return float(mat["g_bulk"]) + modes * (0.5 * hw + self.kt * math.log1p(-math.exp(-hw / self.kt)))

    def _match_m0(self):
        """Cluster-continuum matching at the nucleus: scale every gamma^0 by s so that the continuum chemical
        potential of the m0 dot, g_bulk(T) + 2 s gamma_bar Omega_m / r0 (Wulff shape of volume m0 Omega_m, gamma_bar
        the area-weighted gamma at the reference coverage, r0 the area-weighted facet distance), equals the grid's
        step mu_step = mu_MA + Omega(m0) - Omega(m0 - 1) (best shells included) at reference chemical potentials: the
        pool at half conversion. Growth then starts where the grid says adding MA lowers Omega. PLACEHOLDER bridge."""
        c = self.P["conditions"]
        try:
            pl = pool(self.th, 0.5 * c["c_P0"], c["c_MX0"] - 0.5 * c["c_P0"])
        except ValueError:
            return None
        self._arrays()
        _J, _b, _k, per_m = self.nucleation(pl["mu_MA"], pl["mu_MX"])
        if self.m0 - 1 not in per_m:
            return None
        mu_step = pl["mu_MA"] + per_m[self.m0] - per_m[self.m0 - 1]
        _pl, fam = self.chem(0.5 * c["c_P0"], c["c_MX0"] - 0.5 * c["c_P0"])
        g = {f: fam[f]["gamma"] for f in self.F}
        h0 = self.shape.wulff(g, self.m0 * self.Om)
        area, _v, _s = self.shape(h0)
        gbar = sum(g[f] * area[f] for f in self.F) / sum(area.values())
        r0 = sum(h0[f] * area[f] for f in self.F) / sum(area.values())
        s = (mu_step - self.gb) * r0 / (2.0 * gbar * self.Om)
        if s > 0:
            for p in self.F.values():
                p["gamma0"] = p["gamma0"] * s
        return {"mu_step_eV": mu_step, "g_bulk_eV": self.gb, "mu_MA_ref_eV": pl["mu_MA"], "scale": s, "applied": s > 0}

    # -- pool, coverage, surface energies --
    def chem(self, C_M, C_X):
        pl = pool(self.th, max(C_M, 1e-30), max(C_X, 1e-30))
        c = pl["c_MX_free"]
        fam = {}
        for f, p in self.F.items():
            fam[f] = facet_state(p, c, self.kt)
        # surface free energy and facet areas of the Wulff dot of one unit: classical nucleation, and the continuum
        # tail of the grid's Omega(m) past m0
        gam = {f: fam[f]["gamma"] for f in self.F}
        area, _v, _s = self.shape(self.shape.wulff(gam, self.Om))
        self._unit_area = area
        self._sigma = sum(gam[f] * area[f] for f in self.F)
        return pl, fam

    def k_of(self, m: float, fam: dict) -> float:
        """MX_q of a continuum Wulff dot of m units at the coverage of fam (areas scale as m^(2/3))."""
        a = getattr(self, "_unit_area", None) or {}
        s = m ** (2.0 / 3.0)
        return sum((fam[f]["theta"] * fam[f]["Gamma_max"] + self.F[f]["x_X"] / self.F[f]["site_area"]) * a.get(f, 0.0) * s
                   for f in self.F)

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
            return self._nucleation_continuum(mu_ma, float(om[0]))
        per_m = dict(zip(self.ms.tolist(), om.tolist()))
        ob, mb = om[self.below], self.ms[self.below]
        # past the grid part (m0, or the grid's top when dots there take the grid's mu) the continuum carries on,
        # Omega(m) = Omega(m_e) - (m - m_e) dmu + sigma (m^(2/3) - m_e^(2/3)), so a critical nucleus larger than the
        # grid is still crossed: nuclei are born at m0, or later where the joined Omega(m) has dropped by kT past its
        # maximum (as in the classical mode)
        sigma, dmu = getattr(self, "_sigma", None), mu_ma - self.gb
        m0, m_e = int(self.m0), int(mb[-1])
        if sigma is not None:
            if dmu <= 0 or (2.0 * sigma / (3.0 * dmu)) ** 3 > self.M_CLASSICAL:
                self._m_birth = m0
                return 0.0, float("inf"), None, per_m
            top = int(max(m_e, 2.0 * (2.0 * sigma / (3.0 * dmu)) ** 3, m_e + 10)) + 2
            mt = np.arange(m_e + 1, top + 1)
            ob = np.concatenate([ob, ob[-1] - (mt - m_e) * dmu + sigma * (mt ** (2.0 / 3.0) - m_e ** (2.0 / 3.0))])
            mb = np.concatenate([mb, mt])
        # the barrier is the largest climb from a preceding minimum: a well of stable clusters in the grid (below the
        # monomers) has to be climbed out of, toward the continuum's maximum
        climb = ob - np.minimum.accumulate(ob)
        i_max = int(np.argmax(climb))
        barrier = max(0.0, float(climb[i_max]))
        after = np.nonzero((mb > mb[i_max]) & (ob <= ob[i_max] - self.kt))[0]
        m_b = max(m0, int(mb[after[0]]) if len(after) else m0)
        self._m_birth = m_b
        if m_b in self.Gm and (m_b == m0 or self.grid_mu_on):                # a grid size: the grid's best shell
            k_star = int(self.ks[int(np.argmin(v[np.searchsorted(self.ms, m_b)]))])
        else:
            k_star = None                                                   # continuum coverage (k_of), set by run
        J = nucleation_flux(ob, mb, m_b, self.kt, self.Om, self.P, barrier)
        return J, barrier, k_star, per_m

    M_CLASSICAL = 3000                                  # a critical size beyond this nucleates nothing (dG* > 30 eV)

    def _nucleation_continuum(self, mu_ma, om1):
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
        J = nucleation_flux(om1 + dg, ms, self._m_birth, self.kt, self.Om, self.P, barrier)
        return J, barrier, None, {}

    def dot(self, h, fam):
        """m, k and the areas of a dot with facet distances h."""
        area, vol, support = self.shape(h)
        k = sum((fam[f]["theta"] * fam[f]["Gamma_max"] + self.F[f]["x_X"] / self.F[f]["site_area"]) * area[f] for f in self.F)
        return vol / self.Om, k, area, support

    def size_state(self, cohorts, fam, th):
        """Per cohort: total area A (A^2), equivalent radius r (A), the dot's chemical potential mu_dot - g_bulk (eV)
        of its Wulff shape at gamma(C), mu_dot = g_bulk + (2/3) Omega_m sum_f gamma_f A_f / V, and the area-weighted
        Z-type coverage of its facets (th: the dots' coverage per family)."""
        n = len(cohorts)
        A, r, dmu_dot, th_bar = (np.zeros(n) for _ in range(4))
        for i, co in enumerate(cohorts):
            area = co["area"]
            a = sum(area.values())
            A[i] = a
            r[i] = (3.0 * co["m"] * self.Om / (4.0 * math.pi)) ** (1.0 / 3.0)
            dmu_dot[i] = (2.0 / 3.0) * sum(fam[f]["gamma"] * area[f] for f in self.F) / max(co["m"], 1e-12)
            th_bar[i] = sum(th[f] * area[f] for f in self.F) / a if a > 0 else 0.0
        return A, r, dmu_dot, th_bar

    def size_rates(self, A, r, c_eq, th_bar, C):
        """dm/dt (CdSe units per second per dot) at monomer level C (M): diffusion to the dot (equivalent sphere) in
        series with the surface reaction, for growth and dissolution alike,

            dm/dt = A n (C - C_eq) / (A / (4 pi D r) + 1 / k_eff),   n = 6.022e-4 C molecules/A^3 per M,

        k_eff = k_r (1 - theta_bar): the Z-type shell blocks the free sites, for attachment and detachment alike
        (C_eq carries the thermodynamics, the prefactor is the same both ways); k_r None: the diffusion limit."""
        return self.conductance(A, r, th_bar) * (C - c_eq)

    def conductance(self, A, r, th_bar):
        """dm/dt per M of driving concentration, A n / (A / (4 pi D r) + 1 / k_eff) (size_rates)."""
        res = A / (4.0 * math.pi * self.D * np.maximum(r, 1e-6))                       # s/A
        k_r = self.P["kinetics"].get("k_r")                                             # A/s
        if k_r is not None:
            res = res + 1.0 / (float(k_r) * np.maximum(1.0 - th_bar, 1e-6))
        return A * 6.02214076e-4 / res

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
        c = self.P["conditions"]
        P_, C_M, C_X = float(c["c_P0"]), float(c.get("c_M0", 1e-12)), float(c["c_MX0"])     # c_M0: monomer at t = 0
        t, t_end = 0.0, float(c["t_end"])
        outs = list(np.geomspace(float(c["t_start"]), t_end, int(c["n_out"])))
        cohorts, series, snaps = [], [], []
        pending = pending_m = pending_k = 0.0           # this interval's nuclei (M), the MA and MX_q they took (M)
        th_dots = None                                  # Z-type coverage of the dots' facets per family
        m_dis = self.P["material"]["m_dissolve"]
        if seed:
            _pl, fam = self.chem(C_M, C_X)
            gam = {f: fam[f]["gamma"] for f in self.F}
            for sd in seed:
                h = self.shape.wulff(gam, float(sd["m"]) * self.Om)
                m_, k_, area, _support = self.dot(h, fam)
                cohorts.append({"t_b": 0.0, "N": float(sd["N"]), "h": h, "m": m_, "k": k_, "area": area})
        seeded = {"m": sum(co["N"] * co["m"] for co in cohorts), "k": sum(co["N"] * co["k"] for co in cohorts)}
        while outs:
            pl, fam = self.chem(C_M, C_X)
            mu0, mu_mx = pl["mu_MA"], pl["mu_MX"]
            if th_dots is None:
                th_dots = {f: fam[f]["theta"] for f in self.F}
            Q = self.k_conv * P_ * C_X
            Nc = np.array([co["N"] for co in cohorts])
            # size law: each dot's m changes by dm/dt at one chemical potential; at fixed mu_MX the monomer level in
            # equilibrium with a dot, c_eq = C_M e^((mu_dot - mu)/kT), is the same for every trial C of the solve
            Ac, rc, dmu_dot, th_bar = self.size_state(cohorts, fam, th_dots)
            mc = np.array([co["m"] for co in cohorts])
            dmu_dm = -dmu_dot / (3.0 * np.maximum(mc, 1.0))                  # continuum: mu_dot - g_bulk ~ m^(-1/3)
            mg = mug = None
            if self.grid_mu_on and len(cohorts):
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
            vmax = float(np.abs(rate0 / np.maximum(mc, 1e-12)).max(initial=0.0))           # fastest relative change
            dt = outs[0] - t
            gam_now = {f: fam[f]["gamma"] for f in self.F}
            # the Wulff shape at fixed gamma is scale-free: the dot of one unit, scaled by m^(1/3) (h) and m^(2/3) (A)
            h1 = self.shape.wulff(gam_now, self.Om)
            area1, _v1, _s1 = self.shape(h1)
            sigma_now = sum(gam_now[f] * area1[f] for f in self.F)                    # continuum: (2/3) sigma m^(-1/3)
            with np.errstate(over="ignore", divide="ignore"):
                g_, d_ = rate0 > 0, rate0 < 0
                if mg is not None:
                    # a grid-size dot that meets its own equilibrium within the allowed change stops there exactly
                    # (_stop_at_equilibrium): it does not limit the step
                    for i in np.nonzero(mc <= self.m_top + 1)[0]:
                        frac = 0.02 if g_[i] else 0.2
                        m_try = mc[i] * (1 + frac) if g_[i] else mc[i] * (1 - frac)
                        if rate0[i] != 0 and abs(self._stop_at_equilibrium(mc[i], m_try, mu0 - self.gb, mg, mug, sigma_now) - m_try) > 1e-12:
                            g_[i] = d_[i] = False
                if g_.any():                                                # a growing dot changes m by < 2 % a step
                    dt = min(dt, float((0.02 * np.maximum(mc[g_], 1.0) / rate0[g_]).min()))
                if d_.any():                                                # a dissolving one by < 20 %
                    dt = min(dt, float((0.2 * np.maximum(mc[d_], 1.0) / -rate0[d_]).min()))
                # each dot's own feedback at fixed C: s = d(dm/dt)/dm = -g c_eq (d mu_dot/dm) / kT, in the continuum
                # g c_eq (mu_dot - g_bulk) / (3 m kT). s < 0 (gamma < 0, or a grid size where mu rises with m: the dot's
                # size is a stable equilibrium) is stiff: there
                # the size update is linearised backward Euler, dm = dt rate / (1 - s dt), stable for any dt. s > 0 (the
                # Ostwald instability) stays explicit: with the pool implicit, one dot in a closed pool relaxes as
                # x' = x (1 + s dt) / (1 + N g dt), stable whenever the physics is (N g > s), and the 2 % / 20 % limits
                # bound how far ripening moves in one step (checked: LSW and Wagner rate constants and distributions)
                s_fb = np.minimum(-self.conductance(Ac, rc, th_bar) * c_eq * dmu_dm / self.kt, 0.0)
            if Q > 0:
                dt = min(dt, 0.02 * P_ / Q, 0.02 * C_X / Q)                               # precursor and MX_q < 2 %
            while True:
                def balance(lnC):
                    mu = mu0 + self.kt * (lnC - math.log(C_M))
                    J_, *_ = self.nucleation(mu, mu_mx)
                    U = float((Nc * self.size_rates(Ac, rc, c_eq, th_bar, math.exp(lnC)) / (1.0 - s_fb * dt)).sum()) if len(cohorts) else 0.0
                    return math.exp(lnC) - C_M - dt * (Q - self._m_birth * J_ - U)
                lo, hi = math.log(1e-30), math.log(C_M + Q * dt + 1e-12)
                b_lo, b_hi = balance(lo), balance(hi)
                lnC = brentq(balance, lo, hi, xtol=1e-12, rtol=1e-12) if b_lo < 0 < b_hi else (hi if b_hi < 0 else lo)
                mu = mu0 + self.kt * (lnC - math.log(C_M))
                J, barrier, k_star, _per_m = self.nucleation(mu, mu_mx)
                m_b = self._m_birth                                                    # units of each new nucleus
                if k_star is None:                                                     # classical: continuum coverage
                    k_star = self.k_of(m_b, fam)
                rate = self.size_rates(Ac, rc, c_eq, th_bar, math.exp(lnC)) / (1.0 - s_fb * dt) if len(cohorts) else np.zeros(0)
                # Z-type coverage on the dots moves toward the Langmuir value of the free MX_q, but only as far as the
                # MX_q left in solution allows (the dots' surfaces are a reservoir of MX_q too)
                lam = 1.0
                while True:
                    th_eff = {f: th_dots[f] + lam * (fam[f]["theta"] - th_dots[f]) for f in self.F}
                    k_area = {f: th_eff[f] * fam[f]["Gamma_max"] + self.F[f]["x_X"] / self.F[f]["site_area"] for f in self.F}
                    new, up_m, up_k = [], 0.0, 0.0
                    for i, co in enumerate(cohorts):
                        m_new = max(co["m"] + rate[i] * dt, 1e-6)
                        if mg is not None and min(co["m"], m_new) <= self.m_top + 1:
                            m_new = self._stop_at_equilibrium(co["m"], m_new, mu - self.gb, mg, mug, sigma_now)
                        s = m_new ** (1.0 / 3.0)                                       # shape: Wulff at gamma(C)
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
                if (C_M_new > 0 and C_X_new > 0) or dt < 1e-9:
                    break
                dt *= 0.5
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
            self.steps.append((t, dt, vmax, Q, C_M, C_X_new, J))
            t += dt
            if t >= outs[0] - 1e-12:
                if pending > self.P["material"].get("N_min", 1e-15):                   # this interval's nuclei
                    m_new = pending_m / pending                          # the interval's mean birth size
                    h0 = self.shape.wulff({f: fam[f]["gamma"] for f in self.F}, m_new * self.Om)
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
                    cohorts.append({"t_b": t, "N": pending, "h": {f: min(h0[f], support[f]) for f in self.F},
                                    "m": m_, "k": k_, "area": area})
                    pending = pending_m = pending_k = 0.0
                N = sum(co["N"] for co in cohorts)
                mm = np.array([co["m"] for co in cohorts]) if cohorts else np.zeros(0)
                w = np.array([co["N"] for co in cohorts]) if cohorts else np.zeros(0)
                mean = float((w * mm).sum() / N) if N > 0 else 0.0
                std = float(np.sqrt((w * (mm - mean) ** 2).sum() / N)) if N > 0 else 0.0
                row = {"t": t, "P": P_, "C_M": C_M, "C_X": C_X, "c_MX_free": pl["c_MX_free"], "mu_MA": mu,
                       "mu_MX": mu_mx, "S_bulk": math.exp(min((mu - self.gb) / self.kt, 700)),
                       "barrier": barrier, "J": J, "J0": J * math.exp(min(barrier / self.kt, 700)), "k_star": k_star,
                       "m_birth": m_b,
                       "N_dots": N, "m_mean": mean, "m_std": std,
                       "theta": dict(th_dots), "theta_eq": {f: fam[f]["theta"] for f in self.F},
                       "gamma": {f: fam[f]["gamma"] for f in self.F}}
                series.append(row)
                snaps.append({"t": t, "cohorts": [{"t_b": co["t_b"], "N": co["N"], "m": co["m"], "k": co["k"],
                                                   "h": co["h"], "area": co["area"]} for co in cohorts]})
                if log:
                    log(f"[growth] t = {t:9.3g} s  P = {P_:.3g} M  C_M = {C_M:.3g} M  barrier {barrier:.3f} eV  "
                        f"J = {J:.3g} M/s  dots {N:.3g} M  <m> = {mean:.1f}")
                outs.pop(0)
        return {"params": self.P, "match_m0": self.match, "handover": self.handover, "m0": self.m0, "series": series,
                "snapshots": snaps, "seeded": seeded}


def conservation(res) -> dict:
    """Total MA (precursor + pool + dots) and total MX_q (left + converted + dots) at each output: both constant."""
    c = res["params"]["conditions"]
    s0 = res.get("seeded") or {"m": 0.0, "k": 0.0}                                       # dots present at t = 0
    out = []
    for row, snap in zip(res["series"], res["snapshots"]):
        dots_m = sum(co["N"] * co["m"] for co in snap["cohorts"])
        dots_k = sum(co["N"] * co["k"] for co in snap["cohorts"])
        out.append({"t": row["t"], "MA": row["P"] + row["C_M"] + dots_m - c["c_P0"] - c.get("c_M0", 1e-12) - s0["m"],
                    "MX": row["C_X"] + (c["c_P0"] - row["P"]) + dots_k - c["c_MX0"] - s0["k"]})
    return out


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
    F = res["params"]["families"]
    z = {f: row["theta"][f] * gamma_max(F[f]) * co["area"].get(f, 0.0) for f in F}
    recipe = {k: copy.deepcopy(v) for k, v in preset.items() if k not in ("facets", "size_unit_cells")}
    recipe["facet_distances"] = {f: round(float(co["h"][f]), 3) for f in F}
    recipe["target"] = {"m": int(round(co["m"])), "w": 1.0}
    recipe.setdefault("post_treatment", {})["z_type_target"] = {
        "target_total": int(round(co["k"])), "per_family": {f: round(float(v), 4) for f, v in z.items() if v > 0}}
    return {"t": snap["t"], "cohort_born": co["t_b"], "m": co["m"], "k": co["k"], "recipe": recipe}


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
    Path(a.out).write_text(json.dumps(res, default=float))
    return 0


if __name__ == "__main__":
    sys.exit(main())
