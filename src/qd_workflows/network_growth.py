"""
Network growth model: dots grow over the relaxed structures of the layer-path network, and nothing else.

A dot is one of the network's relaxed states (network_export.py), or below the seed one (m, k) node of the nucleation
grid; the population N_i (mol/L) moves over the steps between states with rates fixed by the grand potential

    Omega_i(t) = G_i(T, eps) - m_i mu_MA(t) - k_i mu_MX(t)          (m = CdSe units, k = CdCl2 units)

where mu_MA, mu_MX come from the solution (grid_thermo.pool: the CdSe(CdCl2)_n complexes share the CdCl2 left with
the free CdCl2). Each step a -> b (b the larger) has one top, and both directions climb to it:

    Omega_top = max(Omega_a + w_s, Omega_b, Omega_a + w_s + B_ab),   w_s = -kT ln(c_s / c0)
    k_ab = nu_ab exp(-(Omega_top - Omega_a) / kT),    k_ba = nu_ab exp(-(Omega_top - Omega_b) / kT)

w_s puts the arriving unit (the CdSe carrier, a given complex, or free CdCl2) at the dot: the forward rate is then
proportional to its concentration and the backward one is not. B_ab is the 2D-island barrier of a layer step, from the
relaxed islands of the step's parent and facet (or the median curve of its facet and size bin), recomputed at every
pool with the shell free at each added CdSe; it is zero for every other step. nu_ab is the collision rate at
c0 = 1 M in series with diffusion (Stokes-Einstein), times the one kinetic scale s_nu.
Detailed balance holds by construction: at a fixed pool the network relaxes to N_i = c0 exp(-Omega_i / kT).

The flask: x(t) CdSe converted from the precursor, dx/dt = Q = k_conv (c_P0 - x) C_X (k_conv at T) with C_X the
salt left in solution, and the pool takes what the dots do not hold,

    C_M = x - sum_i m_i N_i,    C_X = c_MX0 - x - sum_i k_i N_i,

solved together with the population at every implicit (backward Euler) step: sparse LU for N, Newton on (C_M, C_X).
Over a step the salt held in the states is taken at the step's end, so x(t + dt) is grid_thermo.converted over dt
from (c_P0 - x, C_X(t + dt) + x(t + dt) - x): the closed form while the states hold no salt, and no conversion
once the shells have taken it all.
Nothing is taken from growth.py (no k_f, island beta, facet gamma, coverage isotherms or NEB rates).

    python -m qd_workflows.network_growth run --network DIR --thermo thermo_merged.json --params cdse_zb.yaml [--T ..]
"""
from __future__ import annotations

import argparse
import copy
import json
import math
import os
import sys
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp
from scipy.optimize import brentq
from scipy.sparse.linalg import splu

from . import grid_solvation as gs
from .grid_thermo import KB_EV, converted, k_conv_at, pool, thermo_at

NA = 6.02214076e23
KB_J = 1.380649e-23
H_EV = 4.135667696e-15
AMU = 1.66053907e-27
MASS = {"Cd": 112.414, "Se": 78.971, "Cl": 35.453, "C": 12.011, "O": 15.999, "H": 1.008}
TOL = 1e-8                                   # pool residual (relative to the totals)
NEWTON_MAX = 0.2                             # largest Newton step in ln C_M, ln C_X
CH_MAX = 0.2                                 # largest change per time step of ln C and of ln(number of dots)
OUT_CH = (0.3, 0.5, 0.5)                     # dense outputs: change of ln C_M, ln(number of dots), <m> since the last
M_STEP = 1.0                                 # largest change of the dots' <m> per step (CdSe units), once they hold
DOT_SEEN = 1e-12                             # this share of the CdSe (the burst: <m> 13 -> 58 in ~0.3 s)
STALL_TRIES = 400                            # step attempts between two outputs before a run is stopped as stalled
M_DOT = 13                                   # a dot: m >= the seed's 13 CdSe (the grid -> network hand-over)
FAMILIES = ("100", "111A", "111B", "110")
M_BINS = (0, 30, 60, 100, 150, 10 ** 6)      # parent size bins of the median island curves (CdSe units)
SP_NONE, SP_MA, SP_MX = 0, 1, 2              # arriving species of a step; 3 + n: the complex CdSe(CdCl2)_n
KIND = {"iso": 0, "shell": 1, "attach": 2, "layer": 3, "grid": 4, "nucleus": 5}

DEFAULTS = {
    "T": 500.0, "eps": 2.0, "c_P0": 0.025, "c_MX0": 0.03, "k_conv_500K": 1.0, "Ea_conv": 0.0,
    "t_start": 1e-3, "t_end": 600.0, "n_out": 60,
    "s_nu": 1.0,                             # the one kinetic scale: collision prefactor x s_nu
    "eta_mPas": 0.5,                         # solvent viscosity (1-octadecene near 500 K), Stokes-Einstein
    "a_site_A2": 15.0,                       # area of one surface site (the collision cross-section of a site)
    "islands": True,                         # False: no island barriers (attachment-limited layers)
    "prune_eV": 2.0,                         # leave out dots this far above the best of their size (0: exact);
                                             # 2 eV: 5-11x faster, CdSe beyond m = 13 within 1 % of the exact run
    "max_wall_s": 1800.0,                    # a run stops here (and reports why) with the outputs it has
    "dense_out": True,                       # also an output at a step that moves ln C_M, ln(dots) or <m> (OUT_CH)
    "frontier_stop": 0.5,                    # stop once this fraction of the dots sits on states without a growth
                                             # step (the network's edge), the dots holding >= 1% of the CdSe; 1: never
}
LIMITS = {"T": (400.0, 650.0), "eps": (1.0, 10.0), "c_P0": (1e-4, 0.5), "c_MX0": (1e-4, 1.0),
          "k_conv_500K": (1e-4, 1e3), "Ea_conv": (0.0, 2.0), "t_end": (1.0, 7200.0), "s_nu": (1e-6, 1e3),
          "eta_mPas": (0.05, 20.0), "a_site_A2": (1.0, 100.0), "n_out": (5, 200), "prune_eV": (0.0, 20.0),
          "max_wall_s": (1.0, 86400.0), "frontier_stop": (0.01, 1.0)}


def _complex(n: int, solv: dict, m_x: float = MASS["Cl"]) -> tuple[float, float]:
    """(mass amu, radius A) of the CdSe(CdCl2)_n carrier, n = -1: free CdCl2 (MX = Cd X2, X of mass m_x: Cl, or the
    carboxylate of an acetate network); radii from the COSMO surfaces of the run's solvation block (the module's
    CdCl2 defaults where it has none)."""
    mono = solv.get("monomers") or gs.SOLVATION["monomers"]
    cplx = solv.get("complexes", gs.SOLVATION["complexes"]) if "monomers" in solv else gs.SOLVATION["complexes"]
    mx = MASS["Cd"] + 2 * m_x
    if n < 0:
        return mx, math.sqrt(mono["CdCl2"]["S_A2"] / (4 * math.pi))
    area = cplx[str(n)]["S_A2"] if str(n) in cplx else mono["CdSe"]["S_A2"]
    return MASS["Cd"] + MASS["Se"] + n * mx, math.sqrt(area / (4 * math.pi))


def _x_mass(folder: Path, meta: dict) -> float:
    """Mass of one X of the MX = Cd X2 unit, from the seed's relaxed structure: (M - m (Cd + Se) - k Cd) / 2k. Cl when
    the export has no geometry or the seed no MX."""
    st = meta["states"][meta["seed"]]
    if not (folder / "geometry.npz").exists() or st["k"] <= 0:
        return MASS["Cl"]
    g = np.load(folder / "geometry.npz")
    els = meta.get("elements", ("Cd", "Se", "Cl"))
    off = g["offsets"]
    mass = sum(MASS[els[int(a)]] for a in g["symbols"][off[meta["seed"]]:off[meta["seed"] + 1]])
    return (mass - st["m"] * (MASS["Cd"] + MASS["Se"]) - st["k"] * MASS["Cd"]) / (2 * st["k"])


class Network:
    """The states and steps: the grid's nodes with m >= 2 (m = 1 is the pool), then the network's relaxed states
    (those with m <= 13 whose (m, k) the grid has are that grid node), and the steps between them."""

    def __init__(self, folder, thermo: dict, solvation: dict, thermo_dir=None):
        folder = Path(folder)
        self.folder, self.thermo, self.solv, self.thermo_dir = folder, thermo, dict(solvation), thermo_dir
        meta = json.loads((folder / "network.json").read_text())
        arr = np.load(folder / "arrays.npz")
        self.meta, self.G_gas, self.S_row = meta, arr["G_gas"], arr["S"]
        self.m_x = _x_mass(folder, meta)                        # X of the MX unit (Cl, or the carboxylate)
        self.temps = np.asarray(meta["temperatures"], float)
        if not np.allclose(self.temps, thermo["temperatures"]):
            raise ValueError("the network and the grid thermo have different temperature grids")
        grid = sorted((n["m"], n["k"]) for n in thermo["nodes"] if n["m"] >= 2)
        self.n_grid = len(grid)
        gidx = {mk: i for i, mk in enumerate(grid)}
        best = {(n["m"], n["k"]): min(n["states"], key=lambda s: s["E_eV"]) for n in thermo["nodes"]}
        st = meta["states"]
        state_of_row, extra = {}, []
        for r, s in enumerate(st):
            if (s["m"], s["k"]) in gidx and s["m"] <= M_DOT:
                state_of_row[r] = gidx[(s["m"], s["k"])]
            else:
                state_of_row[r] = self.n_grid + len(extra)
                extra.append(r)
        self.state_of_row = state_of_row
        self.n = self.n_grid + len(extra)
        self.row = np.full(self.n, -1)                          # the state's row in arrays.npz (-1: grid node)
        self.row[self.n_grid:] = extra
        self.m = np.array([mk[0] for mk in grid] + [st[r]["m"] for r in extra])
        self.k = np.array([mk[1] for mk in grid] + [st[r]["k"] for r in extra])
        self.grid_keys = grid
        self.S = np.array([float(best[mk].get("S_vdW_A2", 0.0)) for mk in grid] + [self.S_row[r] for r in extra])
        self.label = [f"grid ({m}, {k})" for m, k in grid] + [st[r]["formula"] for r in extra]
        self.reg_id = np.array([-1] * self.n_grid + [st[r]["id"] for r in extra])
        self.seed = state_of_row[meta["seed"]]
        mono_k = sorted(k for (m, k) in best if m == 1)
        self.mono_k = mono_k
        # ---- steps: the network's (re-pointed onto grid nodes), then the grid's Becker-Doering steps, then the
        # nuclei formed from two pool complexes
        A, B, DM, DK, KD, FAM, SITES, SPC, PAR = [], [], [], [], [], [], [], [], []

        def add(a, b, dm, dk, kind, fam=-1, sites=0, spc=SP_NONE, par=-1):
            A.append(a); B.append(b); DM.append(dm); DK.append(dk); KD.append(KIND[kind]); FAM.append(fam)
            SITES.append(sites); SPC.append(spc); PAR.append(par)

        seen = set()
        for e in meta["steps"]:
            a, b = state_of_row[e["a"]], state_of_row[e["b"]]
            if a == b or (a, b) in seen:
                continue
            seen.add((a, b))
            spc = SP_MA if e["dm"] > 0 else (SP_MX if e["dk"] > 0 else (SP_NONE if e["dk"] == 0 else SP_MX))
            if e["dm"] == 0 and e["dk"] < 0:                     # cannot happen (b is the larger), kept explicit
                a, b = b, a
            fam = FAMILIES.index(e["facet"]) if e["facet"] in FAMILIES else -1
            add(a, b, e["dm"], e["dk"], e["kind"], fam, e["sheet"], spc, e["a"] if e["kind"] == "layer" else -1)
        for (m, k), i in gidx.items():
            for n in mono_k:                                     # + CdSe(CdCl2)_n
                j = gidx.get((m + 1, k + n))
                if j is not None:
                    add(i, j, 1, n, "grid", spc=3 + n)
            j = gidx.get((m, k + 1))                             # + CdCl2
            if j is not None:
                add(i, j, 0, 1, "grid", spc=SP_MX)
        self.n_net_steps = len(A)
        src = []                                                 # (n1, n2, b): CdSe(CdCl2)_n1 + CdSe(CdCl2)_n2
        for x, n1 in enumerate(mono_k):
            for n2 in mono_k[x:]:
                j = gidx.get((2, n1 + n2))
                if j is not None:
                    src.append((n1, n2, j))
        self.src = np.array(src, int).reshape(-1, 3)
        self.a, self.b = np.array(A), np.array(B)
        self.dm, self.dk, self.kind = np.array(DM), np.array(DK), np.array(KD)
        self.fam, self.sites, self.spc, self.par = np.array(FAM), np.array(SITES), np.array(SPC), np.array(PAR)
        self.grows = np.zeros(self.n, bool)
        self.grows[self.a[self.dm > 0]] = True                   # states with a step that adds CdSe
        # ---- islands: per (parent row, family) the rows of the island structures and their dm, dk
        self.isl = {}
        for r in meta["islands"]:
            self.isl.setdefault((r["parent"], r["family"]), []).append((r["dm"], r["dk"], r["row"]))
        self.isl_m = {key: st[key[0]]["m"] for key in self.isl}
        self._g_cache = {}

    def subset(self, keep):
        """The network on the states in `keep` (bool per state; every grid node must be kept): states and steps
        re-indexed, the islands and the thermochemistry shared. .full maps the subset's states to this network's."""
        keep = np.asarray(keep, bool)
        if not keep[:self.n_grid].all() or not keep[self.seed]:
            raise ValueError("a subset keeps every grid node and the seed")
        sub = copy.copy(self)
        new = np.full(self.n, -1)
        new[keep] = np.arange(int(keep.sum()))
        sub.full = np.where(keep)[0]
        sub.n = int(keep.sum())
        for name in ("m", "k", "S", "row", "reg_id", "grows"):
            setattr(sub, name, getattr(self, name)[keep])
        sub.label = [lab for lab, q in zip(self.label, keep) if q]
        sub.seed = int(new[self.seed])
        st = keep[self.a] & keep[self.b]
        sub.a, sub.b = new[self.a[st]], new[self.b[st]]
        for name in ("dm", "dk", "kind", "fam", "sites", "spc", "par"):
            setattr(sub, name, getattr(self, name)[st])
        sub.n_net_steps = int(st[:self.n_net_steps].sum())
        sub.grows = np.zeros(sub.n, bool)
        sub.grows[sub.a[sub.dm > 0]] = True
        for name in ("lay", "groups", "group_idx", "curve_keys", "_adj", "_indices", "_indptr", "_slot", "_diag"):
            sub.__dict__.pop(name, None)
        sub._g_cache = {}
        return sub

    def prune(self, T, eps, c_mx_hi, c_mx_lo, window):
        """States to keep for a run whose CdCl2 left goes from c_mx_hi down to c_mx_lo: every grid node, the seed, and
        each dot within `window` eV of the lowest state of its size at some point of that range. Within one size
        mu_CdSe cancels, so the order depends on mu_CdCl2 only, which is checked at both ends and the middle."""
        g, th, _pa = self.G(T, eps)
        keep = np.zeros(self.n, bool)
        keep[:self.n_grid] = True
        keep[self.seed] = True
        for c in (c_mx_hi, math.sqrt(c_mx_hi * c_mx_lo), c_mx_lo):
            om = g - self.k * _pool(th, 1e-4, c)["mu_MX"]
            low = np.full(int(self.m.max()) + 1, np.inf)
            np.minimum.at(low, self.m, om)
            keep |= om - low[self.m] <= window
        return keep

    # ---------------- thermochemistry
    def G(self, T: float, eps: float):
        """(G of every state, the thermo_at() of the grid) at T and eps."""
        key = (round(T, 6), round(eps, 6))
        if key not in self._g_cache:
            th = thermo_at(self.thermo, T, solvation={**self.solv, "eps": eps})
            pa = gs.per_area(eps, self.solv)
            g = np.empty(self.n)
            g[:self.n_grid] = [th["G"][mk] for mk in self.grid_keys]
            rows = self.row[self.n_grid:]
            g[self.n_grid:] = self._interp(T, rows) + pa * self.S_row[rows]
            if len(self._g_cache) > 8:
                self._g_cache.clear()
            self._g_cache[key] = (g, th, pa)
        return self._g_cache[key]

    def _interp(self, T, rows):
        j = int(np.clip(np.searchsorted(self.temps, T) - 1, 0, len(self.temps) - 2))
        w = (T - self.temps[j]) / (self.temps[j + 1] - self.temps[j])
        return (1 - w) * self.G_gas[rows, j] + w * self.G_gas[rows, j + 1]

    def row_G(self, T, eps, rows):
        return self._interp(T, np.asarray(rows)) + gs.per_area(eps, self.solv) * self.S_row[np.asarray(rows)]

    # ---------------- island barriers
    def island_curves(self, T, eps, mu_ma, mu_mx):
        """Per (parent row, family): Omega*(dm) = lowest island or rung with dm CdSe added (any CdCl2), from the parent
        (eV), for dm >= 1."""
        key = ("isl", round(T, 6), round(eps, 6))
        if key not in self._g_cache:                             # G of the islands from the parent: T and eps only
            tab = {}
            for (p, f), L in self.isl.items():
                sel = [x for x in L if x[0] >= 1]
                if sel:
                    dm = np.array([x[0] for x in sel])
                    g = self.row_G(T, eps, [x[2] for x in sel]) - self.row_G(T, eps, [p])[0]
                    tab[(p, f)] = (dm, np.array([x[1] for x in sel]), g, np.unique(dm))
            self._g_cache[key] = tab
        out = {}
        for pf, (dm, dk, g, u) in self._g_cache[key].items():
            om = g - dm * mu_ma - dk * mu_mx
            out[pf] = {int(d): float(om[dm == d].min()) for d in u}
        return out

    @staticmethod
    def barrier(curve: dict) -> float:
        """Highest Omega*(dm) before the curve first drops below the parent; the highest point when it never does."""
        top = 0.0
        for dm in sorted(curve):
            if curve[dm] < 0:
                return top
            top = max(top, curve[dm])
        return top

    def _groups(self):
        """Which curve each layer step takes: its parent's own, the median of its facet and size bin, of its facet, or
        of all islands ({110}, which has none)."""
        self.lay = np.where(self.kind == KIND["layer"])[0]
        have = {key for key, L in self.isl.items() if any(dm >= 1 for dm, _dk, _r in L)}
        self.curve_keys = sorted(have)
        groups, gid, idx = [], {}, np.zeros(len(self.lay), int)
        bin_of = lambda m: max(x for x in M_BINS if x <= m)                        # noqa: E731
        fam_bins = {(f, bin_of(self.isl_m[key])) for key in have for f in [key[1]]}
        fams = {key[1] for key in have}
        for j, i in enumerate(self.lay):
            f = FAMILIES[self.fam[i]]
            key = (int(self.par[i]), f)
            if key in have:
                g = ("own", key)
            elif (f, bin_of(self.m[self.a[i]])) in fam_bins:
                g = ("bin", f, bin_of(self.m[self.a[i]]))
            elif f in fams:
                g = ("family", f)
            else:
                g = ("all",)
            if g not in gid:
                gid[g] = len(groups); groups.append(g)
            idx[j] = gid[g]
        self.groups, self.group_idx = groups, idx

    def barriers(self, T, eps, mu_ma, mu_mx):
        """B of every layer step (eV), and the table of the curves' barriers (own per parent, median per facet and
        size bin, per facet, over all)."""
        if not hasattr(self, "groups"):
            self._groups()
        curves = self.island_curves(T, eps, mu_ma, mu_mx)

        def median_curve(keys):
            pts = {}
            for key in keys:
                for dm, v in curves[key].items():
                    pts.setdefault(dm, []).append(v)
            return {dm: _median(v) for dm, v in pts.items()}

        bin_of = lambda m: max(x for x in M_BINS if x <= m)                        # noqa: E731
        val = {}
        for g in set(self.groups) | {("family", f) for f in ("100", "111A", "111B")} | {("all",)} | \
                {("bin", key[1], bin_of(self.isl_m[key])) for key in curves}:
            if g[0] == "own":
                keys = [g[1]] if g[1] in curves else []
            elif g[0] == "bin":
                keys = [key for key in curves if key[1] == g[1] and bin_of(self.isl_m[key]) == g[2]]
            elif g[0] == "family":
                keys = [key for key in curves if key[1] == g[1]]
            else:
                keys = list(curves)
            val[g] = self.barrier(median_curve(keys)) if keys else 0.0
        B = np.zeros(len(self.a))
        B[self.lay] = np.array([val[g] for g in self.groups])[self.group_idx] if len(self.lay) else 0.0
        table = {"own": {f"{g[1][0]}|{g[1][1]}": v for g, v in val.items() if g[0] == "own"},
                 "bin": {f"{g[1]}|{g[2]}": v for g, v in val.items() if g[0] == "bin"},
                 "family": {g[1]: v for g, v in val.items() if g[0] == "family"}, "all": val[("all",)]}
        return B, table

    # ---------------- prefactors
    def prefactors(self, T, eta_mPas, s_nu, a_site_A2):
        """nu of every step (1/s at c0 = 1 M) and of the nucleus sources: surface collisions in series with
        diffusion to the dot; isomerisations at kT/h."""
        R = np.sqrt(np.maximum(self.S, 1.0) / (4 * math.pi)) * 1e-10          # m
        eta = eta_mPas * 1e-3
        D = KB_J * T / (6 * math.pi * eta * R)

        def nu(i_state, n_sites, mass, r_s):
            r_s = r_s * 1e-10
            d_s = KB_J * T / (6 * math.pi * eta * r_s)
            k_d = 4 * math.pi * (D[i_state] + d_s) * (R[i_state] + r_s) * NA * 1e3
            v = math.sqrt(KB_J * T / (2 * math.pi * mass * AMU))
            k_r = n_sites * a_site_A2 * 1e-20 * v * NA * 1e3
            return s_nu / (1.0 / k_r + 1.0 / k_d)

        sites_all = np.maximum(self.S / a_site_A2, 1.0)
        out = np.empty(len(self.a))
        carriers = {n: _complex(n, self.solv, self.m_x) for n in self.mono_k}
        ma = carriers.get(1, carriers[self.mono_k[0]])
        mx = _complex(-1, self.solv, self.m_x)
        for spc in np.unique(self.spc):
            sel = self.spc == spc
            if spc == SP_NONE:
                out[sel] = KB_EV * T / H_EV
                continue
            mass, r_s = ma if spc == SP_MA else (mx if spc == SP_MX else carriers[spc - 3])
            n_sites = np.where(self.kind[sel] == KIND["layer"], np.maximum(self.sites[sel], 1), sites_all[self.a[sel]])
            out[sel] = nu(self.a[sel], n_sites, mass, r_s)
        src = np.array([nu(np.array([j]), np.array([1.0]), *carriers[n2])[0] for _n1, n2, j in self.src])
        return out, src

    # ---------------- rates
    def rates(self, T, eps, pl, nu, nu_src, islands=True):
        """Forward and backward rate constants of every step, the source terms of the nuclei, Omega and B."""
        g, th, _pa = self.G(T, eps)
        kt = KB_EV * T
        mu_ma, mu_mx = pl["mu_MA"], pl["mu_MX"]
        om = g - self.m * mu_ma - self.k * mu_mx
        c_m = pl["c_MA"]
        x_n = pl["speciation"]
        w = np.zeros(len(self.a))
        w[self.spc == SP_MA] = -kt * math.log(c_m)
        w[self.spc == SP_MX] = -kt * math.log(pl["c_MX_free"])
        for n in self.mono_k:
            sel = self.spc == 3 + n
            if sel.any():
                w[sel] = -kt * math.log(max(c_m * x_n.get(n, 0.0), 1e-300))
        if islands:
            B, table = self.barriers(T, eps, mu_ma, mu_mx)
        else:
            B, table = np.zeros(len(self.a)), None
        oa, ob = om[self.a], om[self.b]
        top = np.maximum(np.maximum(oa + w, ob), oa + w + B)
        kf = nu * np.exp(-(top - oa) / kt)
        kr = nu * np.exp(-(top - ob) / kt)
        # nuclei: CdSe(CdCl2)_n1 (Omega = -kT ln c_n1) + CdSe(CdCl2)_n2 (w) -> (2, n1 + n2)
        s_f, s_r = np.zeros(len(self.src)), np.zeros(len(self.src))
        for i, (n1, n2, j) in enumerate(self.src):
            c1, c2 = c_m * x_n.get(n1, 0.0), c_m * x_n.get(n2, 0.0)
            if c1 <= 0 or c2 <= 0:
                continue
            o1, w2 = -kt * math.log(c1), -kt * math.log(c2)
            t_ = max(o1 + w2, om[j])
            sym = 0.5 if n1 == n2 else 1.0
            s_f[i] = sym * nu_src[i] * c1 * math.exp(-(t_ - o1) / kt)   # mol/L/s of nuclei (c1 = N of the source)
            s_r[i] = sym * nu_src[i] * math.exp(-(t_ - om[j]) / kt)     # 1/s back to the pool
        return {"kf": kf, "kr": kr, "s_f": s_f, "s_r": s_r, "Omega": om, "B": B, "table": table, "w": w}

    def _pattern(self):
        """CSC pattern of M: every step's four entries, the sources' and the whole diagonal, zeros kept, so every M
        (and I - dt M) of a run has one structure; and the slot of each of matrix()'s entries in it."""
        n, a, b, j, d = self.n, self.a, self.b, self.src[:, 2], np.arange(self.n)
        key, slot = np.unique(np.concatenate([a, b, a, b, j, d]) * n + np.concatenate([b, a, a, b, j, d]),
                              return_inverse=True)
        self._indices = (key % n).astype(np.int32)
        self._indptr = np.searchsorted(key // n, np.arange(n + 1)).astype(np.int32)
        self._slot, self._diag = slot[:-n], slot[-n:]

    def matrix(self, r):
        """dN/dt = M N + s."""
        n = self.n
        if not hasattr(self, "_slot"):
            self._pattern()
        data = np.concatenate([r["kf"], r["kr"], -r["kf"], -r["kr"], -r["s_r"]])
        M = sp.csc_matrix((np.bincount(self._slot, weights=data, minlength=len(self._indices)), self._indices,
                           self._indptr), shape=(n, n))
        s = np.bincount(self.src[:, 2], weights=r["s_f"], minlength=n) if len(self.src) else np.zeros(n)
        return M, s

    def system(self, M, dt):
        """I - dt M on M's pattern."""
        data = -dt * M.data
        data[self._diag] += 1.0
        return sp.csc_matrix((data, self._indices, self._indptr), shape=M.shape)

    def fluxes(self, r, N):
        """Net flux of every step (mol/L/s, a -> b positive) and of the nucleus sources."""
        return r["kf"] * N[self.a] - r["kr"] * N[self.b], r["s_f"] - r["s_r"] * N[self.src[:, 2]]


def _median(v):
    """np.median of a short list (the same value), without numpy's per-call overhead."""
    v = sorted(v)
    h = len(v) // 2
    return float(v[h]) if len(v) % 2 else float((v[h - 1] + v[h]) / 2)


def _threads() -> int:
    """Threads of the sparse LU: QDW_LU_THREADS, else the cores this process may use, at most 4 (a 7.6k-state LU
    is no faster on more: M5 Pro 4 threads 0.019 s, 8 threads 0.026 s)."""
    if os.environ.get("QDW_LU_THREADS"):
        return max(1, int(os.environ["QDW_LU_THREADS"]))
    n = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else (os.cpu_count() or 1)
    return min(4, n)


class _Solver:
    """Solves (I - dt M) x = b by sparse LU, one factorization per matrix, kept for further solves with the same
    matrix (the Jacobian's). Every matrix of a run has M's pattern (Network.system), so the fill-reducing ordering
    and the symbolic analysis are done once. MUMPS (python-mumps, conda-forge: METIS ordering, OpenMP threads) where
    it is installed, else SciPy's SuperLU (serial, ordering and analysis at every factorization; 5x slower on one core).
    QDW_LU=superlu forces SuperLU. GMRES on the previous LU is not tried: it almost never reached 1e-11 on a new
    matrix (619 LU for 610 matrices on the acetate network) and cost ~10 % of the run."""

    def __init__(self, n, threads=None):
        self.n, self.A, self.lu, self.n_lu, self.n_it = n, None, None, 0, 0
        self.threads = threads or _threads()
        self.mumps = None
        if os.environ.get("QDW_LU", "mumps").lower() != "superlu":
            try:
                import mumps
                self.mumps = mumps.Context()
            except ImportError:
                pass
        self.backend = f"MUMPS ({self.threads} threads)" if self.mumps is not None else "SuperLU"
        self._analysed, self.n_fallback, self.fallback = False, 0, None

    def factor(self, A):
        self.A, self.lu = A, None
        self.n_lu += 1
        if self.mumps is not None:
            inst = self.mumps
            try:
                inst.set_matrix(A)
                if not self._analysed:
                    inst.mumps_instance.icntl[6] = 0             # no value-based column permutation: the analysis
                    inst.mumps_instance.icntl[16] = self.threads  # is structural, valid for every matrix of the run
                    inst.analyze(ordering="metis")               # (threads set before it: MUMPS error -58 else)
                    self._analysed = True
                inst.factor(reuse_analysis=True)
                return
            except Exception as e:                               # noqa: BLE001 (mumps.MUMPSError): this matrix by
                self.n_fallback += 1                             # SuperLU
                self.fallback = f"{self.n_fallback} matrices by SuperLU after MUMPS errors, the last: {e}"
        self.lu = splu(A, permc_spec="MMD_AT_PLUS_A", options={"SymmetricMode": True})

    def solve(self, A, b):
        if A is not self.A:
            self.factor(A)
        return self.lu.solve(b) if self.lu is not None else self.mumps.solve(b)


def _pool(th, c_m, c_x):
    pl = pool(th, c_m, c_x)
    pl["c_MA"] = c_m
    pl["c_X_left"] = c_x
    return pl


def simulate(net: Network, cond: dict | None = None, log=None, stop=None) -> dict:
    """Backward Euler over [t_start, t_end] with the pool solved at every step; outputs on n_out log-spaced times,
    and with dense_out also at every step that moves ln C_M, ln(number of dots) or <m> by OUT_CH since the last
    output (the dots' once they hold DOT_SEEN of the CdSe). A step moves the dots' <m> by at most M_STEP, so the
    nucleation burst is resolved.
    The run ends early, with the outputs it has, its last state and the reason in out["stopped"], when stop() is
    true, after max_wall_s, when STALL_TRIES step attempts do not reach the next output, when the pool does not
    converge, or when frontier_stop of the dots sit on states without a growth step: past that the dots pile up at
    the network's edge, C_M collapses and the run says nothing about growth.
    prune_eV > 0: on the subset of Network.prune (the CdCl2 left from c_MX0 down to c_MX0 - c_P0, at least 1e-6 M),
    populations returned on the whole network."""
    c = {**DEFAULTS, **(cond or {})}
    for key, (lo, hi) in LIMITS.items():
        if key in c and not lo <= float(c[key]) <= hi:
            raise ValueError(f"{key} = {c[key]} is outside [{lo}, {hi}]")
    if float(c["prune_eV"]) > 0:
        keep = net.prune(float(c["T"]), float(c["eps"]), float(c["c_MX0"]),
                         max(float(c["c_MX0"]) - float(c["c_P0"]), 1e-6), float(c["prune_eV"]))
        sub = net.subset(keep)
        r = simulate(sub, {**c, "prune_eV": 0.0}, log, stop)
        r["cond"] = c
        full = []
        for Ns in r["N"]:
            Nf = np.zeros(net.n, np.float32)
            Nf[sub.full] = Ns
            full.append(Nf)
        r["N"], r["states"] = full, sub.n
        return r
    T, eps = float(c["T"]), float(c["eps"])
    c_p0, c_mx0 = float(c["c_P0"]), float(c["c_MX0"])
    k_conv = k_conv_at(T, float(c["k_conv_500K"]), float(c["Ea_conv"]))
    t0, t_end = float(c["t_start"]), float(c["t_end"])
    t_out = list(np.geomspace(t0, t_end, int(c["n_out"])))
    nu, nu_src = net.prefactors(T, float(c["eta_mPas"]), float(c["s_nu"]), float(c["a_site_A2"]))
    _g, th, _pa = net.G(T, eps)
    isl = bool(c["islands"])

    def conv(x_, c_x1, dt):
        """x at the step's end: the conversion over dt with the salt the states hold fixed at its end value,
        x1 - x = converted(dt, P = c_P0 - x, MX = C_X1 + (x1 - x)); monotone in x1 - x on [0, c_P0 - x]."""
        p = c_p0 - x_
        if p <= 0.0:
            return x_
        g = lambda d: d - converted(dt, k_conv, p, c_x1 + d)            # noqa: E731
        return x_ + (p if g(p) <= 0.0 else brentq(g, 0.0, p, xtol=1e-18, rtol=1e-15))
    mvec, kvec = net.m.astype(float), net.k.astype(float)
    # the CdSe and X each step and nucleus source moves from the pool (b - a; the fast isomerisations move none)
    d_mk = np.array([mvec[net.b] - mvec[net.a], kvec[net.b] - kvec[net.a]])
    js = net.src[:, 2]
    mk_src = np.array([mvec[js], kvec[js]])
    N = np.zeros(net.n)
    t = t0
    x = converted(t, k_conv, c_p0, c_mx0)                 # nothing held before t_start
    u = np.array([math.log(x), math.log(c_mx0 - x)])
    out = {"cond": c, "t": [], "alpha": [], "c_MA": [], "c_X_left": [], "c_MX_free": [], "mu_MA": [], "mu_MX": [],
           "n_dots": [], "n_clusters": [], "m_mean": [], "m_std": [], "nucleation": [], "frontier": [],
           "layer_flux": {f: [] for f in FAMILIES}, "cdse_flux": {f: [] for f in FAMILIES}, "attach_flux": [],
           "barriers": [], "N": [], "steps": 0, "lu": 0, "wall_s": 0.0, "stopped": None}
    w0 = time.time()
    solver = _Solver(net.n)
    dudt = None
    diag = {"newton": 0, "backtrack": 0, "failed": 0, "rejected": 0}
    dot = net.m >= M_DOT

    def assemble(u_):
        pl = _pool(th, math.exp(u_[0]), math.exp(u_[1]))
        r = net.rates(T, eps, pl, nu, nu_src, isl)
        return pl, r, *net.matrix(r)

    def moved(r, N1):
        """CdSe and X taken from the pool over the step, dt x (sum of the steps' and sources' net fluxes x what
        they move): equal to (m, k) . (N1 - N) for the exact N1, without the round-off of the isomerisations'
        ~1e13 1/s terms, which leaves (m, k) . N1 noisy at ~1e-7 of the totals, above TOL."""
        J, j_s = net.fluxes(r, N1)
        return d_mk @ J + mk_src @ j_s

    def solve(u_, dt):
        """The step's N and conversion x1 for a pool guess, and the pool residual F (relative to the totals)."""
        pl, r, M, s = assemble(u_)
        A = net.system(M, dt)
        N1 = solver.solve(A, N + dt * s)
        x1 = conv(x, math.exp(u_[1]), dt)
        sc = np.array([x1, c_mx0 - x1])
        F = (np.exp(u_) + np.array([mvec @ N, kvec @ N]) + dt * moved(r, N1) - sc) / sc
        return {"F": F, "N": N1, "pl": pl, "r": r, "A": A, "M": M, "s": s, "sc": sc, "dt": dt, "x1": x1}

    def jacobian(st, u_, dt):
        """Exact Jacobian of F in ln C_M, ln C_X: the rates' derivatives by differences of the assembly, the
        population's through the same matrix, and the conversion's through C_X."""
        J = np.diag(np.exp(u_) / st["sc"])
        dx1 = (conv(x, math.exp(u_[1] + 1e-6), dt) - st["x1"]) / 1e-6
        J[:, 1] -= (st["F"] + 1.0) / st["sc"] * np.array([1.0, -1.0]) * dx1
        r, N1 = st["r"], st["N"]
        for q in range(2):
            h = np.zeros(2); h[q] = 1e-6
            _p, rh, Mh, sh = assemble(u_ + h)
            dN = solver.solve(st["A"], dt * (((Mh - st["M"]) @ N1) + (sh - st["s"])) / 1e-6)
            dr = {key: (rh[key] - r[key]) / 1e-6 for key in ("kf", "kr", "s_f", "s_r")}
            d_flux = moved(dr, N1) + moved({**r, "s_f": np.zeros_like(r["s_f"])}, dN)  # rates' and N1's change
            J[:, q] += dt * d_flux / st["sc"]
        return J

    def record(t_, N_, pl, r):
        nd = float(N_[dot].sum())
        out["t"].append(t_); out["alpha"].append(x / c_p0)
        out["c_MA"].append(pl["c_MA"]); out["c_X_left"].append(pl["c_X_left"]); out["c_MX_free"].append(pl["c_MX_free"])
        out["mu_MA"].append(pl["mu_MA"]); out["mu_MX"].append(pl["mu_MX"])
        out["n_dots"].append(nd); out["n_clusters"].append(float(N_[~dot].sum()))
        mm = float((N_[dot] * net.m[dot]).sum() / nd) if nd > 0 else 0.0
        out["m_mean"].append(mm)
        out["m_std"].append(float(math.sqrt(max((N_[dot] * (net.m[dot] - mm) ** 2).sum() / nd, 0.0))) if nd > 0 else 0.0)
        J, _js = net.fluxes(r, N_)
        cross = (net.m[net.a] < M_DOT) & (net.m[net.b] >= M_DOT)
        out["nucleation"].append(float(J[cross].sum()))
        out["frontier"].append(float(N_[dot & ~net.grows].sum() / nd) if nd > 0 else 0.0)
        lay = net.kind == KIND["layer"]
        for i, f in enumerate(FAMILIES):
            sel = lay & (net.fam == i)
            out["layer_flux"][f].append(float(J[sel].sum()))
            out["cdse_flux"][f].append(float((J[sel] * net.dm[sel]).sum()))
        out["attach_flux"].append(float((J[(net.kind == KIND["attach"])] * net.dm[net.kind == KIND["attach"]]).sum()))
        out["barriers"].append(r["table"])
        out["N"].append(N_.astype(np.float32))

    pl0 = _pool(th, math.exp(u[0]), math.exp(u[1]))
    acc = (pl0, net.rates(T, eps, pl0, nu, nu_src, isl))   # the last accepted step's pool and rates
    record(t, N, *acc)
    t_out.pop(0)
    dt = t0 * 1e-3
    tries, wall_max, f_stop = 0, float(c["max_wall_s"]), float(c["frontier_stop"])
    dense = bool(c["dense_out"])

    def _moved_since(N_, pl):
        """The step moved ln C_M, ln(number of dots) or <m> by OUT_CH since the last output."""
        if abs(math.log(pl["c_MA"] / out["c_MA"][-1])) >= OUT_CH[0]:
            return True
        nd = N_[dot].sum()
        if nd <= 0 or mvec[dot] @ N_[dot] < DOT_SEEN * x:
            return False
        nd0 = out["n_dots"][-1]
        return (nd0 <= 0 or abs(math.log(nd / nd0)) >= OUT_CH[1]
                or abs((N_[dot] @ mvec[dot]) / nd - out["m_mean"][-1]) >= OUT_CH[2])

    def end(why):
        out["stopped"] = why
        if out["t"][-1] < t:
            record(t, N, *acc)
        if log:
            log(why)

    while t_out:
        tries += 1
        why = ("stopped from the page" if stop is not None and stop() else
               f"wall-time limit {wall_max:.0f} s" if time.time() - w0 > wall_max else
               f"stalled: {tries - 1} step attempts without reaching t = {t_out[0]:.4g} s" if tries > STALL_TRIES else None)
        if why:
            end(f"{why} (at t = {t:.4g} s, dt = {dt:.3g} s, {diag})")
            break
        dt = min(dt, t_out[0] - t)
        t1 = t + dt
        if dudt is None:                                 # predictor: the pool's last rate of change, else the supply
            x1 = conv(x, math.exp(u[1]), dt)
            v = np.array([math.log(max(math.exp(u[0]) + (x1 - x), x1 * 1e-12)), u[1] - (x1 - x) / max(c_mx0 - x, 1e-30)])
        else:
            v = u + dudt * dt
        ok = False
        st = solve(v, dt)
        for _it in range(25):
            if np.abs(st["F"]).max() < TOL:
                ok = True
                break
            diag["newton"] += 1
            try:
                d = -np.linalg.solve(jacobian(st, v, dt), st["F"])
            except np.linalg.LinAlgError:
                break
            d *= min(1.0, NEWTON_MAX / max(np.abs(d).max(), 1e-30))   # the clusters hold ~exp(m ln C_M)
            lam, better = 1.0, False
            for _ls in range(8):
                st2 = solve(v + lam * d, dt)
                if np.abs(st2["F"]).max() < np.abs(st["F"]).max() * (1 - 0.2 * lam) or np.abs(st2["F"]).max() < TOL:
                    better = True
                    break
                lam *= 0.3
                diag["backtrack"] += 1
            if not better:
                break
            v, st = v + lam * d, st2
        F, N1, pl, r, x1 = st["F"], np.maximum(st["N"], 0.0), st["pl"], st["r"], st["x1"]
        if not ok:
            diag["failed"] += 1
            dt *= 0.25
            if dt < 1e-14 * max(t, 1.0):
                end(f"the pool does not converge (at t = {t:.4g} s, {diag})")
                break
            continue
        ch = max(abs(v[0] - u[0]), abs(v[1] - u[1]))
        nd0, nd1 = N[dot].sum(), N1[dot].sum()
        if nd0 > 0 and (mvec[dot] @ N1[dot]) > 1e-3 * x1:      # the dots' number, once they hold CdSe
            ch = max(ch, abs(math.log(max(nd1, 1e-300)) - math.log(nd0)))
        if nd0 > 0 and nd1 > 0 and (mvec[dot] @ N1[dot]) > DOT_SEEN * x1:   # their size, from the first ones
            ch = max(ch, CH_MAX * abs(mvec[dot] @ N1[dot] / nd1 - mvec[dot] @ N[dot] / nd0) / M_STEP)
        if ch > CH_MAX and dt > t0 * 1e-6:
            diag["rejected"] += 1
            dt *= 0.3
            continue
        dudt = (v - u) / dt
        t, x, u, N = t1, x1, v, N1
        acc = (pl, r)
        out["steps"] += 1
        if t_out and t >= t_out[0] * (1 - 1e-12):
            record(t, N, pl, r)
            t_out.pop(0)
            tries = 0
        elif dense and _moved_since(N, pl):
            record(t, N, pl, r)
            if log:
                log(f"t = {t:.4g} s: alpha {out['alpha'][-1]:.3f}, C_M {pl['c_MA']:.3g} M, dots {out['n_dots'][-1]:.3g} M, "
                    f"<m> {out['m_mean'][-1]:.1f}, {out['steps']} steps, {solver.n_lu} LU, {solver.n_it} GMRES it, {diag}, {time.time() - w0:.1f} s")
        nd = N[dot].sum()
        if f_stop < 1.0 and nd > 0 and mvec[dot] @ N[dot] >= 1e-2 * x:
            fr = N[dot & ~net.grows].sum() / nd
            if fr > f_stop:
                end(f"network edge: {fr:.0%} of the dots are on states without a growth step (the network ends at "
                    f"m = {int(net.m.max())}; at t = {t:.4g} s, <m> {(N[dot] @ net.m[dot]) / nd:.1f})")
                break
        dt *= min(3.0, max(1.1, 0.5 * CH_MAX / max(ch, 1e-6)))
    out["lu"] = solver.n_lu
    out["gmres_it"] = solver.n_it
    out["lu_backend"] = solver.backend + (f" ({solver.fallback})" if solver.fallback else "")
    out["diag"] = diag
    out["wall_s"] = time.time() - w0
    return out


def stationary_residual(net: Network, T, eps, pl, nu=None, nu_src=None) -> float:
    """Largest |(M N_eq + s)_i| / (outflow_i N_eq_i) with N_eq = c0 exp(-Omega / kT), on the states whose N_eq is
    representable: 0 up to round-off when the rates obey detailed balance."""
    if nu is None:
        nu, nu_src = net.prefactors(T, DEFAULTS["eta_mPas"], 1.0, DEFAULTS["a_site_A2"])
    r = net.rates(T, eps, pl, nu, nu_src)
    kt = KB_EV * T
    z = -r["Omega"] / kt
    ok = np.abs(z) < 600
    sub = np.where(ok)[0]
    idx = -np.ones(net.n, int); idx[sub] = np.arange(len(sub))
    neq = np.zeros(net.n); neq[ok] = np.exp(z[ok])
    M, s = net.matrix(r)
    res = M @ neq + s
    out = -M.diagonal() * neq
    keep = ok.copy()
    for a_, b_ in ((net.a, net.b), (net.b, net.a)):           # a state next to one outside the range is not tested
        keep[a_[~ok[b_]]] = False
    keep &= out > 0
    return float(np.max(np.abs(res[keep]) / out[keep])) if keep.any() else 0.0


def load(network, thermo_path, params_path) -> Network:
    import yaml
    thermo = json.loads(Path(thermo_path).read_text())
    solv = yaml.safe_load(Path(params_path).read_text())["solvation"]
    return Network(network, thermo, solv, Path(thermo_path).resolve().parent)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python -m qd_workflows.network_growth")
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("run", "check"):
        p = sub.add_parser(name)
        p.add_argument("--network", required=True, help="network_export folder")
        p.add_argument("--thermo", required=True, help="the grid's thermo_merged.json")
        p.add_argument("--params", required=True, help="growth template with the solvation block")
        for key, v in DEFAULTS.items():
            if not isinstance(v, bool):
                p.add_argument(f"--{key.replace('_', '-')}", dest=key, type=type(v), default=v)
        p.add_argument("--no-islands", action="store_true")
        p.add_argument("--json", help="write the run here (without the populations)")
    a = ap.parse_args(argv)
    net = load(a.network, a.thermo, a.params)
    print(f"{net.n} states ({net.n_grid} grid nodes), {len(net.a)} steps, {len(net.src)} nucleus sources, "
          f"{len(net.isl)} island curves", flush=True)
    cond = {key: getattr(a, key) for key in DEFAULTS if hasattr(a, key)}
    cond["islands"] = not a.no_islands
    if a.cmd == "check":
        _g, th, _pa = net.G(a.T, a.eps)
        for alpha in (0.05, 0.2, 0.5):
            x = alpha * a.c_P0
            pl = _pool(th, x, a.c_MX0 - x)
            print(f"alpha {alpha}: mu_MA {pl['mu_MA']:.3f}, mu_MX {pl['mu_MX']:.3f} eV; stationary residual "
                  f"{stationary_residual(net, a.T, a.eps, pl):.2e}")
        return 0
    res = simulate(net, cond, log=lambda s: print(s, flush=True))
    print(f"{res['steps']} steps, {res['lu']} LU solves, {res['wall_s']:.1f} s"
          + (f"; stopped early: {res['stopped']}" if res["stopped"] else ""))
    if a.json:
        Path(a.json).write_text(json.dumps({k: v for k, v in res.items() if k != "N"}, default=float))
    return 0


if __name__ == "__main__":
    sys.exit(main())
