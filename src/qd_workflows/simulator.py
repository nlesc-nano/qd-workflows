"""
Nucleation and growth simulator: a local web page on top of the growth model (growth.py) and QD_builder.

The page sets the conditions (T, precursor and MX_q concentrations, conversion rate, time), runs the growth model on
the grid thermochemistry, plots the pools, the nucleation and the size distribution, and builds any dot of any time
with QD_builder (target mode on the dot's facet distances, then the Z-type target), shown atom by atom.

    python -m qd_workflows.simulator --thermo thermo_merged.json --params growth_params.yaml \
        --preset templates/nucleation/cdse_zb.yaml [--port 8772]

Run it in an environment with QD_builder, FastAPI and uvicorn (the QDSpace backend's).
"""
from __future__ import annotations

import argparse
import copy
import io
import json
import math
import tempfile
import threading
import time
import uuid
from collections import OrderedDict
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import yaml

from . import growth as gr
from .grid_thermo import KB_EV, converted, k_conv_at, pool, thermo_at

PAGE = Path(__file__).parent / "templates" / "simulator" / "index.html"
CLASSICAL_PAGE = Path(__file__).parent / "templates" / "simulator" / "classical.html"
IDEAL_PAGE = Path(__file__).parent / "templates" / "simulator" / "ideal.html"
CLASSICAL_LIMITS = {                        # material sliders of the classical twin
    "gamma_J_m2": [0.01, 2.0], "Vm_A3": [20.0, 200.0], "C_inf_M": [1e-15, 1.0], "dH_sol_eV": [0.0, 3.0],
    "D_cm2_s": [1e-8, 1e-2], "k_r_A_s": [1e-2, 1e12], "J0_M_s": [1e-10, 1e40],
}
LIMITS = {                                  # slider ranges of the page (the grid thermo spans its own T range)
    "T": [300.0, 700.0], "c_P0": [0.005, 0.5], "c_MX0": [0.01, 1.0], "k_conv_500K": [1e-8, 1e6],
    "t_end": [1.0, 36000.0], "Ea_conv": [0.0, 2.0], "t": [1e-3, 1e5], "eps": [1.0, 40.0],
    "k_r": [1.0, 1e10],                     # A/s, surface reaction of the size law (1e10: the diffusion limit)
    # calibration factors (tuned to target experiments): activity coefficients and the island step scale
    "gamma_MA": [0.3, 3.0], "gamma_MX": [0.3, 3.0], "beta_scale": [0.5, 1.5],
    # syringe injection of the precursor (volumes in mL on the page, L in the model)
    "c_syringe_M": [1e-3, 1.0], "rate_mL_h": [0.01, 1000.0], "V0_mL": [1.0, 1000.0], "inj_duration_s": [1.0, 36000.0],
    # temperature program: growth temperature, switch time, ramp
    "T_growth": [300.0, 700.0], "t_switch_s": [0.0, 36000.0], "ramp_s": [0.0, 3600.0],
}
PAGE_DEFAULTS = {"gamma_MA": 1.0, "gamma_MX": 1.0, "c_syringe_M": 0.03, "rate_mL_h": 0.6, "V0_mL": 10.0,
                 "inj_duration_s": 3600.0, "t_switch_s": 120.0, "ramp_s": 60.0}


def solvent_eps(x) -> float:
    """Solvent dielectric constant as the model takes it: 1 (gas phase) or from 1.8 (the cavity term's range) up."""
    x = min(max(float(x), 1.0), 40.0)
    return 1.0 if x <= 1.0 else max(x, 1.8)


def _finite(x):
    """JSON has no inf/nan: an infinite barrier (no supercritical nucleus) goes to the page as null."""
    return None if isinstance(x, float) and not math.isfinite(x) else x


SETTINGS = ("T", "c_P0", "c_MX0", "k_conv_500K", "Ea_conv", "t_end")   # the conditions block of growth.PARAMS
EDITABLE = {                                # facet parameters the page may change (advanced, checked on the server)
    "gamma0": (1.0, 200.0, "meV/A^2"),      # bare facet energy
    "eps0": (0.0, 3.0, "eV"),               # binding energy of one MX_q on the empty facet (> 0: bound)
}
MAX_ATOMS = 8000                            # the builder button refuses dots above this many atoms (23k atoms: 4 min)
# X-type unit MX_q of the grid: the page's label and the X in the grid's formulas, which keep the Cl placeholders
# (the nucleation workflow dresses them with the carboxylate; its "CdCl2" reference is then CdX2)
LIGANDS = {"CdCl2": {"label": "CdCl₂", "x": "Cl"}, "acetate": {"label": "Cd(OAc)₂", "x": "(OAc)"},
           "oleate": {"label": "Cd(OA)₂", "x": "(OA)"}}


def grid_ligand(thermo: dict) -> str:
    """The grid's X-type ligand: the nucleation workflow names the MX_q reference monomer CdCl2-<ligand>."""
    src = str((thermo.get("references") or {}).get("CdCl2", {}).get("source", ""))
    name = Path(src).parent.name
    return name.split("-", 1)[1] if name.startswith("CdCl2-") and name.split("-", 1)[1] in LIGANDS else "CdCl2"


class Simulator:
    def __init__(self, thermo: dict, params: dict, preset_cfg: dict, max_atoms: int = MAX_ATOMS,
                 thermo_dir: str | Path | None = None, ligand: str | None = None, nucleation_only: bool = False,
                 m0: int | None = None):
        """nucleation_only: the page's nucleation panel alone (a grid without growth parameters of its own, e.g. a
        carboxylate grid while the facet parameters are CdCl2's); m0 then fixes the nucleus size (the automatic
        hand-over needs the facets)."""
        self.thermo, self.params, self.cfg = thermo, params, preset_cfg
        self.ligand = ligand or grid_ligand(thermo)
        self.xl = LIGANDS[self.ligand]
        self.nucleation_only = bool(nucleation_only)
        if m0 is not None:
            self.params = params = copy.deepcopy(params)
            params["material"]["m0"] = int(m0)
        self.thermo_dir = Path(thermo_dir) if thermo_dir else None      # the grid's structure files are relative to it
        self.max_atoms = max_atoms
        # solvent-excluded radius of each grid node's lowest state (transport.add_geometry), A
        self._rse = {}
        for nd in thermo.get("nodes", []):
            st = min(nd["states"], key=lambda s: s["E_eV"]) if nd.get("states") else None
            if st and "R_se_A" in st:
                self._rse[(nd["m"], nd["k"])] = float(st["R_se_A"])
        self.runs: OrderedDict[str, dict] = OrderedDict()
        self._grids: dict = {}
        self._handovers: dict = {}
        self.build_lock = threading.Lock()  # the builder redirects stdout: one build at a time
        t = thermo["temperatures"]
        trange = [max(LIMITS["T"][0], min(t)), min(LIMITS["T"][1], max(t))]
        self.limits = dict(LIMITS, T=trange, T_growth=trange)

    def defaults(self) -> dict:
        c = self.params["conditions"]
        F = self.params["families"]
        eps = float(self.params.get("solvation", {}).get("eps", 1.0))
        k_r = self.params["kinetics"].get("k_r")
        kin = self.params["kinetics"]
        return {"conditions": {**{k: c[k] for k in SETTINGS}, "t": 10.0, "eps": eps,
                               "k_r": (None if k_r == "facets" else
                                       float(k_r) if k_r is not None else LIMITS["k_r"][1]),
                               **PAGE_DEFAULTS, "T_growth": c["T"],
                               "beta_scale": float(kin.get("beta_scale", 1.0))},
                "growth_law": kin.get("growth_law", "size"), "facet_kinetics": kin.get("facets") or {},
                "growth_control": kin.get("growth_control", "mixed"),
                "k_r_mode": "facets" if k_r == "facets" else "set",
                "m0": self._m0({**c, "eps": eps})[0],
                "solvation_note": self.params.get("solvation", {}).get("note", ""),
                "limits": self.limits, "families": list(F), "max_atoms": self.max_atoms,
                "material": self.cfg.get("material", ""), "phase": self.cfg.get("phase", ""),
                "ligand": {"name": self.ligand, "label": self.xl["label"]}, "nucleation_only": self.nucleation_only,
                "surface": {"source": self.params.get("provenance", {}).get("families", "placeholder"),
                            "match_m0": bool(self.params["material"].get("match_m0", False)),
                            "editable": {k: {"min": lo, "max": hi, "unit": u} for k, (lo, hi, u) in EDITABLE.items()},
                            "families": {f: {k: F[f][k] for k in ("gamma0", "eps0", "Gamma_max", "w", "tau", "k0",
                                                                  "site_area") if k in F[f]} for f in F}}}

    def _overrides(self, P: dict, fams: dict | None) -> bool:
        """Facet parameters typed on the page (gamma0, eps0 per family), checked against EDITABLE. Typed values
        are used as they are: the cluster-continuum rescaling of gamma0 (match_m0) is switched off for that run."""
        if not fams:
            return False
        changed = False
        for f, vals in fams.items():
            if f not in P["families"]:
                raise ValueError(f"unknown facet family {f!r}")
            for k, v in (vals or {}).items():
                if k not in EDITABLE:
                    raise ValueError(f"{k} cannot be changed from the page")
                lo, hi, unit = EDITABLE[k]
                v = float(v)
                if not (lo <= v <= hi) or v != v:
                    raise ValueError(f"{k} of {{{f}}} must be between {lo} and {hi} {unit} (got {v})")
                if v != P["families"][f][k]:
                    P["families"][f][k] = v
                    changed = True
        if changed:
            P["material"]["match_m0"] = False
        return changed

    def _clip(self, cond: dict, k: str, default=None):
        x = cond.get(k, default)
        if x is None:
            return None
        lo, hi = self.limits[k]
        return min(max(float(x), lo), hi)

    def params_for(self, cond: dict) -> dict:
        """The growth parameters for the page's conditions: T, concentrations, conversion, solvent, k_r; the
        calibration factors (gamma_MA, gamma_MX, beta_scale); the syringe (inject: c_syringe_M, rate_mL_h, V0_mL,
        inj_duration_s; the flask's own precursor is then c_P0, 0 for a pure injection) and the temperature program
        (two_T: T_growth, t_switch_s, ramp_s); growth_law "shape", "size" or "sphere" (ideal growth, with
        gamma_J_m2 and k_r_A_s; left out, both come from the facets)."""
        P = copy.deepcopy(self.params)
        c = P["conditions"]
        for k in SETTINGS:
            if cond.get(k) is not None:
                c[k] = self._clip(cond, k)
        if cond.get("eps") is not None and "solvation" in P:
            P["solvation"]["eps"] = solvent_eps(cond["eps"])
        if cond.get("k_r") is not None and P["kinetics"].get("k_r") != "facets":
            P["kinetics"]["k_r"] = self._clip(cond, "k_r")      # only where k_r is a set number, not computed
        for k in ("gamma_MA", "gamma_MX"):
            c[k] = self._clip(cond, k, 1.0)
        P["kinetics"]["beta_scale"] = self._clip(cond, "beta_scale", P["kinetics"].get("beta_scale", 1.0))
        if cond.get("growth_control") in ("mixed", "diffusion", "reaction"):
            P["kinetics"]["growth_control"] = cond["growth_control"]   # transport and surface reaction, or one of them
        if cond.get("nucleation_omega") in ("grid", "continuum"):
            P["nucleation"]["omega"] = cond["nucleation_omega"]   # continuum: CNT with the Wulff sigma (comparison)
        if cond.get("growth_law") in ("shape", "size", "sphere"):
            P["kinetics"]["growth_law"] = cond["growth_law"]
        if P["kinetics"].get("growth_law") == "sphere":     # ideal growth: one gamma, one k_r (None: from the facets)
            P["kinetics"]["sphere"] = {k: (min(max(float(cond[k]), CLASSICAL_LIMITS[k][0]), CLASSICAL_LIMITS[k][1])
                                           if cond.get(k) is not None else None) for k in ("gamma_J_m2", "k_r_A_s")}
        if cond.get("inject"):
            c["c_P0"] = float(cond.get("c_P0_flask", 0.0) or 0.0)
            c["injection"] = {"c_syringe_M": self._clip(cond, "c_syringe_M", PAGE_DEFAULTS["c_syringe_M"]),
                              "rate_L_h": self._clip(cond, "rate_mL_h", PAGE_DEFAULTS["rate_mL_h"]) * 1e-3,
                              "V0_L": self._clip(cond, "V0_mL", PAGE_DEFAULTS["V0_mL"]) * 1e-3,
                              "duration_s": self._clip(cond, "inj_duration_s", PAGE_DEFAULTS["inj_duration_s"])}
        if cond.get("two_T"):
            c["T_growth"] = self._clip(cond, "T_growth", c["T"])
            c["t_switch_s"] = self._clip(cond, "t_switch_s", PAGE_DEFAULTS["t_switch_s"])
            c["ramp_s"] = self._clip(cond, "ramp_s", PAGE_DEFAULTS["ramp_s"])
        inj = c.get("injection")
        n_p = c["c_P0"] * (inj["V0_L"] if inj else 1.0) + (inj["c_syringe_M"] * inj["rate_L_h"] / 3600.0
                                                           * inj["duration_s"] if inj else 0.0)
        n_x = c["c_MX0"] * (inj["V0_L"] if inj else 1.0)
        if n_x <= n_p:
            raise ValueError(f"{self.xl['label']} must exceed all the Se precursor (in the flask and from the syringe): "
                             "every conversion takes one MX_q")
        return P

    def _growth_allowed(self):
        if self.nucleation_only:
            raise ValueError(f"nucleation only: there are no growth (facet) parameters for the {self.ligand} grid yet")

    def _formula(self, f: str) -> str:
        return f if self.xl["x"] == "Cl" else f.replace("Cl", self.xl["x"])

    def run(self, cond: dict) -> dict:
        self._growth_allowed()
        P = self.params_for(cond)
        modified = self._overrides(P, cond.get("families"))
        t0 = time.time()
        g = gr.Growth(self.thermo, P)
        res = g.run()
        rid = uuid.uuid4().hex[:10]
        self.runs[rid] = {"res": res, "growth": g}
        while len(self.runs) > 8:
            self.runs.popitem(last=False)
        Om = P["material"]["Omega_m"]
        snaps = []
        for s in res["snapshots"]:
            cos = []
            for i, co in enumerate(s["cohorts"]):
                A = sum(co["area"].values()) or 1.0
                cos.append({"i": i, "t_b": co["t_b"], "N": co["N"], "m": co["m"], "k": co["k"],
                            "d_nm": (6 * co["m"] * Om / math.pi) ** (1 / 3) / 10,
                            "n_atoms": int(round(2 * co["m"] + 3 * co["k"])),
                            "facets": {f: a / A for f, a in co["area"].items() if a > 0}})
            snaps.append({"t": s["t"], "cohorts": cos})
        keep = ("t", "T", "V_L", "P", "C_M", "C_X", "c_MX_free", "S_bulk", "barrier", "J", "J0", "k_star", "N_dots",
                "m_mean", "m_std", "theta", "theta_eq", "gamma", "k_eff", "nucleated_M")
        used = g.F                                           # derived: solvation, binding at T, match_m0 rescaling
        return {"run_id": rid, "seconds": round(time.time() - t0, 2), "conditions": P["conditions"],
                "growth_law": P["kinetics"].get("growth_law", "size"),
                "growth_control": P["kinetics"].get("growth_control", "mixed"), "sphere": res.get("sphere"),
                "nucleation_window": res.get("nucleation_window"), "survivor_window": res.get("survivor_window"),
                "match_m0": res["match_m0"], "modified": modified,
                "surface_used": {f: {"gamma0": used[f]["gamma0"], "eps0": used[f].get("eps0"),
                                     "dG_bind": used[f]["dG_bind"]} for f in used},
                "dG_thermal": g.dg_th, "m0": res["m0"], "handover": res["handover"],
                "series": [{k: _finite(r[k]) for k in keep} for r in res["series"]],
                "snapshots": snaps, "conservation": gr.conservation(res)[-1],
                "omega_run": self._omega_along(P, res)}

    def _omega_along(self, P: dict, res: dict) -> dict:
        """The grid's Omega(m) along a growth run, at the pool the run leaves (precursor converted minus what the
        nuclei and dots took), with the same definition as the live panel (grid up to m0, barrier = largest
        Omega(m) - Omega(1)): the barrier and m* per output time and Delta Omega(m, t) for the map."""
        c = P["conditions"]
        eps = float(P.get("solvation", {}).get("eps", 1.0))
        m0 = self._m0({**c, "eps": eps})[0]
        out = {"t": [], "barrier": [], "m_star": [], "C_M": [], "S": [], "m": None, "dOmega": [], "R_se_A": []}
        cols, rcols = [], []
        for row in res["series"]:
            if row["C_M"] <= 1e-20:
                continue
            try:
                o = self._omega(row.get("T", c["T"]), row["C_M"], row["C_X"], m0, eps,
                                float(c.get("gamma_MA", 1.0)), float(c.get("gamma_MX", 1.0)))
            except ValueError:
                continue
            rows = [r for r in o["per_m"] if r["m"] <= o["m0"]]
            out["t"].append(row["t"])
            out["barrier"].append(o["barrier"])
            out["m_star"].append(o["m_star"])
            out["C_M"].append(row["C_M"])
            out["S"].append(o["S"])
            out["m"] = [r["m"] for r in rows]
            cols.append([r["Omega_m"] - rows[0]["Omega_m"] for r in rows])
            rcols.append([r["R_se_A"] for r in rows])
        out["dOmega"] = [list(col) for col in zip(*cols)]              # rows m = 1..m0, columns the run's times
        out["R_se_A"] = [list(col) for col in zip(*rcols)]             # the clusters' radius, same layout
        return out

    def build(self, rid: str, t: float, cohort: int | None, which: str = "mode", full: bool = False) -> dict:
        """Build the dot with QD_builder's target mode: the construction centres only (about a second), or with
        full=True also the layer roundings h_f + r delta_f, r in {-1, 0, 1} (closer to m, ten times slower)."""
        from builder.scripts.generate_library import _build
        self._growth_allowed()
        if rid not in self.runs:
            raise KeyError("run not found (the server keeps the last 8 runs): run the model again")
        res = self.runs[rid]["res"]
        r = gr.builder_recipe(res, t, self.cfg["preset"], which, cohort)
        est = int(round(2 * r["m"] + 3 * r["k"]))
        if est > self.max_atoms:
            raise ValueError(f"this dot has about {est} atoms, above the {self.max_atoms}-atom limit of the page")
        recipe = r["recipe"]
        if not full:
            recipe["target"]["roundings"] = [0]
        t0 = time.time()
        with self.build_lock, tempfile.TemporaryDirectory() as tmp:
            out = _build(self.cfg["cif"], recipe, self.cfg.get("positive_q_mode", "add"), Path(tmp) / "b")
        secs = round(time.time() - t0, 2)
        if "error" in out:
            raise RuntimeError(out["error"] + "\n" + "\n".join(out.get("log", "").splitlines()[-15:]))
        syms = list(out["symbols"])
        pts = np.asarray(out["pts"], float)
        cation, anion = self.cfg["native_elements"]
        xyz = io.StringIO()
        xyz.write(f"{len(syms)}\nt = {r['t']:.4g} s, cohort born at {r['cohort_born']:.4g} s, growth model m = "
                  f"{r['m']:.1f}, k = {r['k']:.1f}\n")
        for s, p in zip(syms, pts):
            xyz.write(f"{s:2s} {p[0]:12.6f} {p[1]:12.6f} {p[2]:12.6f}\n")
        comp = {e: syms.count(e) for e in sorted(set(syms))}
        man = out.get("manifest") or {}
        best = man.get("target_search", {}).get("best")
        led = man.get("z_type_target_ledger") or {}
        ztype = {k: led.get(k) for k in ("target_total", "k_X", "direction", "units", "shortfall", "error")}
        ztype["per_family"] = {f: {k: v.get(k) for k in ("requested", "added", "removed", "candidates") if k in v}
                               for f, v in (led.get("per_family") or {}).items()}
        return {"t": r["t"], "cohort_born": r["cohort_born"], "model": {"m": r["m"], "k": r["k"]},
                "built": {"m": syms.count(anion), "k": syms.count(cation) - syms.count(anion), "n_atoms": len(syms),
                          "composition": comp}, "recipe": recipe, "seconds": secs,
                "search": best and {"centre": best.get("centre"), "roundings": best.get("roundings")}, "z_type": ztype,
                "xyz": xyz.getvalue()}

    def _g_bulk(self, T: float) -> float:
        """Bulk free energy per MA unit at T, as the growth model has it (growth.Growth.g_bulk)."""
        return gr.Growth.g_bulk(SimpleNamespace(P=self.params, T=float(T), kt=KB_EV * float(T)))

    def _solv(self, eps: float) -> dict:
        return {**self.params.get("solvation", {"eps": 1.0}), "eps": solvent_eps(eps)}

    def _grid(self, T: float, eps: float = 1.0) -> dict:
        """Free energies of the grid at T (implicit solvent at eps) as arrays: G(m, k) (inf where no node), m and k."""
        key = (round(float(T), 3), solvent_eps(eps))
        if key not in self._grids:
            th = thermo_at(self.thermo, key[0], self._solv(eps))
            ms = np.array(sorted({m for m, _ in th["G"]}))
            ks = np.array(sorted({k for _, k in th["G"]}))
            Gp = np.full((len(ms), len(ks)), np.inf)
            for (m, k), g in th["G"].items():
                Gp[np.searchsorted(ms, m), np.searchsorted(ks, k)] = g
            self._grids[key] = {"th": th, "ms": ms, "ks": ks, "Gp": Gp}
            while len(self._grids) > 32:
                self._grids.pop(next(iter(self._grids)))
        return self._grids[key]

    def _m0(self, cond: dict):
        """Nucleus size for these conditions: the fixed m0 of the parameters, or growth's automatic hand-over (the
        size where the grid's two-step chemical potential meets the continuum), cached per T and concentrations."""
        if self.params["material"]["m0"] != "auto":
            return int(self.params["material"]["m0"]), None
        key = tuple(round(float(cond.get(k, self.params["conditions"][k])), 9) for k in ("T", "c_P0", "c_MX0"))
        eps = solvent_eps(cond.get("eps", self.params.get("solvation", {}).get("eps", 1.0)))
        key = key + (eps,)
        if key not in self._handovers:
            P = copy.deepcopy(self.params)
            P["conditions"].update(T=key[0], c_P0=key[1], c_MX0=key[2])
            if "solvation" in P:
                P["solvation"]["eps"] = eps
            g = gr.Growth(self.thermo, P)
            self._handovers[key] = (g.m0, g.handover)
            while len(self._handovers) > 64:
                self._handovers.pop(next(iter(self._handovers)))
        return self._handovers[key]

    def _omega(self, T: float, C_M: float, C_X: float, m0: int | None = None, eps: float = 1.0,
               gamma_ma: float = 1.0, gamma_mx: float = 1.0) -> dict:
        """Omega(m, k) = G(m, k) - m mu_MA - k mu_MX of the grid at the pool (C_M monomer supply, C_X MX_q left),
        Omega(m) = -kT ln sum_k e^(-Omega(m,k)/kT), the barrier max_{m <= m0} Omega(m) - Omega(1), m*, k*(m),
        P(k|m) and J (the nucleation step of the growth model); every species in implicit solvent at eps. The
        supersaturation S = e^((mu_MA - g_bulk(T))/kT), as the growth run's S_bulk."""
        gr_ = self._grid(T, eps)
        th, ms, ks, Gp = gr_["th"], gr_["ms"], gr_["ks"], gr_["Gp"]
        kt = KB_EV * T
        pl = pool(th, max(C_M, 1e-30), max(C_X, 1e-30), gamma_ma=gamma_ma, gamma_mx=gamma_mx)
        mu_ma, mu_mx = pl["mu_MA"], pl["mu_MX"]
        gb = self._g_bulk(T)
        v = Gp - ms[:, None] * mu_ma - ks[None, :] * mu_mx
        vmin = v.min(axis=1)
        P = np.exp(-(v - vmin[:, None]) / kt)
        om_m = vmin - kt * np.log(P.sum(axis=1))
        P /= P.sum(axis=1, keepdims=True)
        m0 = int(m0 if m0 is not None else self._m0({"T": T, "eps": eps})[0])
        rows = []
        for i, m in enumerate(ms.tolist()):
            j = int(np.argmin(v[i]))
            node = th["nodes"][(m, int(ks[j]))]
            # the cluster's size: solvent-excluded radius (CdCl2 shell included) weighted over its shells k at this pool
            rk = [(P[i][q], self._rse.get((m, int(k)))) for q, k in enumerate(ks.tolist()) if P[i][q] > 0]
            wsum = sum(w for w, r in rk if r is not None)
            r_se = sum(w * r for w, r in rk if r is not None) / wsum if wsum > 0 else None
            rows.append({"m": m, "Omega_m": float(om_m[i]), "k_star": int(ks[j]), "Omega_min": float(v[i, j]),
                         "mean_k": float((P[i] * ks).sum()), "formula": self._formula(node["formula"]),
                         "source": node["states"][0].get("source", "grid"), "R_se_A": r_se})
        below = [r for r in rows if r["m"] <= m0]
        top = max(below, key=lambda r: r["Omega_m"])
        barrier = max(0.0, top["Omega_m"] - below[0]["Omega_m"])
        k0 = next((r["k_star"] for r in rows if r["m"] == m0), None)
        return {"T": T, "C_M": C_M, "C_X": C_X, "mu_MA": mu_ma, "mu_MX": mu_mx, "c_MX_free": pl["c_MX_free"],
                "g_bulk": gb, "S": math.exp(min((mu_ma - gb) / kt, 700.0)),
                "mean_n": pl["mean_n"], "speciation": {str(n): float(x) for n, x in pl["speciation"].items()},
                "m0": m0, "barrier": barrier, "m_star": top["m"], "k_star_m0": k0,
                "J": (J := gr.nucleation_flux(om_m, ms, m0, kt, float(self.params["material"]["Omega_m"]), self.params, barrier,
                                              act=gamma_ma)),
                "J0": J * math.exp(min(barrier / kt, 700.0)),
                "ms": ms.tolist(), "ks": ks.tolist(), "Omega": [[None if not np.isfinite(x) else float(x) for x in r] for r in v],
                "P_k": P.round(6).tolist(), "per_m": rows}

    def pre_pool(self, P: dict, times) -> list:
        """The pool before any nucleus at each time (s): precursor conversion P -> MA at k_conv(T(t)) P C_X, the
        syringe adding precursor and diluting the flask, the temperature program. Without injection and at one
        temperature the closed form (grid_thermo.converted). Returns (T, C_M, C_X) per time."""
        from scipy.integrate import solve_ivp
        c = P["conditions"]
        prog = SimpleNamespace(P=P)
        inj = c.get("injection")
        times = [float(x) for x in times]
        if not inj and c.get("T_growth") is None:
            k = k_conv_at(c["T"], c["k_conv_500K"], c["Ea_conv"])
            out = []
            for tt in times:
                x = converted(tt, k, c["c_P0"], c["c_MX0"])
                out.append((c["T"], x, c["c_MX0"] - x))
            return out
        q = inj["rate_L_h"] / 3600.0 if inj else 0.0
        t_stop = inj["duration_s"] if inj else 0.0

        def rhs(tt, y):
            p_, cm, cx, v = y
            k = k_conv_at(gr.Growth.temperature_at(prog, tt), c["k_conv_500K"], c["Ea_conv"])
            qq = q if tt < t_stop else 0.0
            conv = k * max(p_, 0.0) * max(cx, 0.0)
            feed = qq * (inj["c_syringe_M"] - p_) / v if inj else 0.0
            return [feed - conv, conv - qq * cm / v, -conv - qq * cx / v, qq]

        ts = sorted(set(t for t in times if t > 0))
        breaks = sorted({b for b in [t_stop] + gr.Growth.temperature_breaks(prog) if 0 < b < ts[-1]})
        y = np.array([c["c_P0"], 0.0, c["c_MX0"], inj["V0_L"] if inj else 1.0])
        got, t_now = {0.0: y}, 0.0
        for seg_end in breaks + [ts[-1]]:
            pts = sorted({t for t in ts if t_now < t <= seg_end} | {seg_end})
            sol = solve_ivp(rhs, (t_now, seg_end), y, method="LSODA", t_eval=pts, rtol=1e-8, atol=1e-14)
            for tt, yy in zip(sol.t, sol.y.T):
                got[float(tt)] = yy
            y, t_now = sol.y[:, -1], seg_end
        out = []
        for tt in times:
            yy = got[max(tt, 0.0)] if max(tt, 0.0) in got else got[min(got, key=lambda x: abs(x - tt))]
            out.append((gr.Growth.temperature_at(prog, tt), max(float(yy[1]), 0.0), max(float(yy[2]), 1e-30)))
        return out

    def omega_live(self, cond: dict, n_curve: int = 41) -> dict:
        """Omega at time t from the start, from the precursor alone (no nuclei yet: the pre-nucleation pool, the
        qdw nucleation Omega page), with the syringe, the temperature program and the activity coefficients; and
        the barrier along t on a log grid for the same conditions."""
        P = self.params_for(cond)
        c = P["conditions"]
        t = min(max(float(cond.get("t", 10.0)), self.limits["t"][0]), self.limits["t"][1])
        eps = solvent_eps(cond.get("eps", self.params.get("solvation", {}).get("eps", 1.0)))
        m0, handover = self._m0({**c, "eps": eps})
        g_ma, g_mx = c["gamma_MA"], c["gamma_MX"]
        grid_t = list(np.geomspace(1e-3, 1e4, n_curve))
        pools = self.pre_pool(P, [t] + grid_t)
        T, cm, cx = pools[0]
        out = self._omega(T, cm, cx, m0, eps, g_ma, g_mx)
        inj = c.get("injection")
        if inj:
            t_in = min(t, inj["duration_s"])
            n_p = c["c_P0"] * inj["V0_L"] + inj["c_syringe_M"] * inj["rate_L_h"] / 3600.0 * t_in
            v_t = inj["V0_L"] + inj["rate_L_h"] / 3600.0 * t_in
        else:
            n_p, v_t = c["c_P0"], 1.0
        out.update(t=t, T_now=T, alpha=out["C_M"] * v_t / n_p if n_p > 0 else 0.0, conditions={**c, "eps": eps},
                   handover=handover, injection=bool(inj))
        curve, heat = [], []
        ms = [r["m"] for r in out["per_m"] if r["m"] <= out["m0"]]
        for tt, (Tt, cm_t, cx_t) in zip(grid_t, pools[1:]):
            if cm_t <= 1e-20:
                continue
            try:
                o = self._omega(Tt, cm_t, cx_t, m0, eps, g_ma, g_mx)
            except ValueError:
                continue
            curve.append({"t": float(tt), "T": Tt, "barrier": o["barrier"], "m_star": o["m_star"], "J": o["J"],
                          "C_M": o["C_M"], "S": o["S"]})
            om1 = o["per_m"][0]["Omega_m"]
            heat.append([r["Omega_m"] - om1 for r in o["per_m"] if r["m"] <= out["m0"]])
        out["curve"] = curve
        # Delta Omega(m, t) = Omega(m) - Omega(1) along the time: rows m = 1..m0, columns the curve's t
        out["heat"] = {"t": [c["t"] for c in curve], "m": ms, "dOmega": [list(col) for col in zip(*heat)]}
        return out

    def omega(self, rid: str, t: float) -> dict:
        """Omega at time t of a growth run: the run's pool (precursor converted, minus what nuclei and dots took),
        at the run's temperature then and with its activity coefficients."""
        if rid not in self.runs:
            raise KeyError("run not found (the server keeps the last 8 runs): run the model again")
        res, g = self.runs[rid]["res"], self.runs[rid]["growth"]
        row = min(res["series"], key=lambda r: abs(r["t"] - t))
        c = g.P["conditions"]
        om = self._omega(row.get("T", g.T), row["C_M"], row["C_X"], g.m0, float(g.solv.get("eps", 1.0)),
                         float(c.get("gamma_MA", 1.0)), float(c.get("gamma_MX", 1.0)))
        return {**om, "t": row["t"]}

    # -- the classical twin (classical.py) --
    def from_atomistic(self, cond: dict) -> dict:
        """The classical parameters the atomistic model implies at these conditions (reference pool: half of the
        precursor converted): gamma_eff = sigma_W / ((36 pi)^(1/3) V_m^(2/3)) of the Wulff dot at the pool's
        coverage; C_inf where the pool's mu_MA meets g_bulk(T), and its enthalpy from T +- 10 K; D of the monomer
        complexes (Stokes-Einstein); k_r of the Wulff dot from the facet rates; J0 = J e^(dOmega*/kT) of the grid."""
        from .classical import J_M2
        self._growth_allowed()
        P = self.params_for(cond)
        g = gr.Growth(self.thermo, P)
        g._arrays()
        cm, cx = g.ref_pool()

        def c_inf():
            lo, hi = math.log(1e-30), math.log(10.0)
            for _ in range(200):
                mid = 0.5 * (lo + hi)
                lo, hi = (mid, hi) if g.pool(math.exp(mid), cx)["mu_MA"] < g.gb else (lo, mid)
            return math.exp(0.5 * (lo + hi))

        T = g.T
        cinf = {}
        for dT in (-10.0, 10.0, 0.0):                              # ends at T
            g.set_temperature(T + dT)
            g._arrays()
            cinf[dT] = c_inf()
        dH = -KB_EV * math.log(cinf[10.0] / cinf[-10.0]) / (1.0 / (T + 10.0) - 1.0 / (T - 10.0))
        pl, fam = g.chem(cm, cx)
        J, barrier, _k, _pm = g.nucleation(pl["mu_MA"], pl["mu_MX"])
        gam_eff = g._sigma / ((36.0 * math.pi) ** (1.0 / 3.0) * g.Om ** (2.0 / 3.0))          # eV/A^2
        h = g.shape.wulff({f: fam[f]["gamma"] for f in g.F}, g.Om)
        area = g.shape(h)[0]
        th = {f: fam[f]["theta"] for f in g.F}
        k_r = g.mean_keff([{"N": 1.0, "area": area}], th)
        if k_r is None or not math.isfinite(k_r):
            k_r = 1e4
        D = (g._D1 if (g.stokes and g._D1) else g.D) / 1e16                                    # cm^2/s
        out = {"gamma_J_m2": gam_eff / J_M2, "Vm_A3": g.Om, "C_inf_M": cinf[0.0], "T_ref": T, "dH_sol_eV": max(dH, 0.0),
               "D_cm2_s": D, "k_r_A_s": k_r, "J0_mode": "zeldovich",
               "J0_M_s": (J * math.exp(min(barrier / g.kt, 700.0))) if math.isfinite(barrier) else None}
        raw = dict(out)
        notes = []
        for k, (lo, hi) in CLASSICAL_LIMITS.items():          # a classical sphere needs values in the sliders' range
            if out.get(k) is not None and not lo <= out[k] <= hi:
                notes.append(f"{k} = {out[k]:.3g} is outside the classical range [{lo:g}, {hi:g}]: clipped")
                out[k] = min(max(out[k], lo), hi)
        return {"material": out, "raw": raw, "notes": notes,
                "reference_pool": {"C_M": cm, "C_X": cx, "S_at_ref": cm / cinf[0.0],
                                                    "barrier_eV": barrier if math.isfinite(barrier) else None},
                "sources": {"gamma_J_m2": "Wulff surface energy of the facets at the pool's CdCl2 coverage",
                            "Vm_A3": "volume per CdSe (MACE zinc blende)",
                            "C_inf_M": "monomer level where the pool's mu_MA = g_bulk(T)",
                            "dH_sol_eV": "from C_inf at T - 10 K and T + 10 K",
                            "D_cm2_s": "monomer complexes, Stokes-Einstein on their hydrodynamic radius",
                            "k_r_A_s": "facet incorporation rates over the Wulff dot's facets, k_f (1 - theta_f)",
                            "J0_M_s": "Becker-Doring prefactor of the grid at the reference pool (for J0_mode set)"}}

    def classical_run(self, cond: dict, material: dict | None) -> dict:
        """The classical twin at the page's conditions (syringe and temperature program included) and the
        material parameters of the sliders (clipped to CLASSICAL_LIMITS)."""
        from . import classical as cl
        self._growth_allowed()
        P = copy.deepcopy(cl.PARAMS)
        P["conditions"].update({k: v for k, v in self.params_for(cond)["conditions"].items()
                                if k not in ("gamma_MA", "gamma_MX")})
        for k, v in (material or {}).items():
            if k in CLASSICAL_LIMITS and v is not None:
                lo, hi = CLASSICAL_LIMITS[k]
                P["material"][k] = min(max(float(v), lo), hi)
            elif k in ("J0_mode",) and v in ("zeldovich", "set"):
                P["material"][k] = v
            elif k == "control" and v in ("mixed", "diffusion", "reaction"):
                P["material"][k] = v
            elif k == "T_ref" and v is not None:
                P["material"][k] = float(v)
        t0 = time.time()
        m = cl.Classical(P)
        res = m.run()
        return {"seconds": round(time.time() - t0, 2), "material": P["material"], "conditions": P["conditions"],
                "nucleation_window": res.get("nucleation_window"), "survivor_window": res.get("survivor_window"),
                "series": [{k: _finite(v) for k, v in r.items()} for r in res["series"]],
                "snapshots": res["snapshots"], "conservation": cl.conservation(res)[-1], "steps": m.steps}

    def structure(self, m: int, k: int) -> dict:
        """xyz of the grid node's lowest state (the structure file named in the thermo, relative to its folder)."""
        node = next((n for n in self.thermo["nodes"] if n["m"] == m and n["k"] == k), None)
        if node is None:
            raise KeyError(f"no grid node m = {m}, k = {k}")
        st = node["states"][0]
        path = Path(st["file"]) if self.thermo_dir is None else self.thermo_dir / st["file"]
        if not path.exists():
            raise FileNotFoundError(f"structure file {st['file']} is not in {self.thermo_dir}")
        return {"m": m, "k": k, "formula": self._formula(node["formula"]), "E_eV": st["E_eV"], "source": st.get("source", "grid"),
                "n_states": len(node["states"]), "file": st["file"], "xyz": path.read_text()}


def make_app(sim: Simulator):
    from fastapi import Body, FastAPI, HTTPException
    from fastapi.responses import HTMLResponse

    app = FastAPI(title="QD nucleation and growth simulator")

    @app.get("/", response_class=HTMLResponse)
    def page():
        return PAGE.read_text()

    @app.get("/classical", response_class=HTMLResponse)
    def ideal_page():
        return IDEAL_PAGE.read_text()                # the same nucleation, ideal growth, against the full model

    @app.get("/cnt", response_class=HTMLResponse)
    def classical_page():
        return CLASSICAL_PAGE.read_text()            # textbook CNT and LaMer with their own parameters

    @app.get("/api/classical/defaults")
    def classical_defaults():
        d = sim.defaults()
        try:
            mapped = sim.from_atomistic(d["conditions"])
        except ValueError as exc:
            raise HTTPException(400, str(exc))
        return {**d, "classical_limits": CLASSICAL_LIMITS, "mapped": mapped}

    @app.get("/api/classical/presets")
    def classical_presets():
        """Literature values of the classical parameters per material (survey of II-VI, III-V, IV-VI nanocrystal
        models; null where not reported), with their sources and flags."""
        return json.loads((CLASSICAL_PAGE.parent / "literature_presets.json").read_text())

    @app.post("/api/classical/from_atomistic")
    def classical_from_atomistic(body: dict = Body(...)):
        try:
            return sim.from_atomistic(body)
        except (KeyError, ValueError, TypeError) as exc:
            raise HTTPException(400, str(exc).strip("'\""))

    @app.post("/api/classical/run")
    def classical_run(body: dict = Body(...)):
        try:
            return sim.classical_run(body.get("conditions") or {}, body.get("material"))
        except (KeyError, ValueError, TypeError) as exc:
            raise HTTPException(400, str(exc).strip("'\""))

    @app.get("/api/defaults")
    def defaults():
        return sim.defaults()

    @app.post("/api/run")
    def run(body: dict = Body(...)):
        try:
            return sim.run(body)
        except ValueError as exc:
            raise HTTPException(400, str(exc))

    @app.post("/api/omega")
    def omega(body: dict = Body(...)):
        try:
            return sim.omega(str(body["run_id"]), float(body["t"]))
        except (KeyError, ValueError, TypeError) as exc:
            raise HTTPException(400, str(exc).strip("'\""))

    @app.post("/api/omega_live")
    def omega_live(body: dict = Body(...)):
        try:
            return sim.omega_live(body)
        except (KeyError, ValueError, TypeError) as exc:
            raise HTTPException(400, str(exc).strip("'\""))

    @app.get("/api/structure")
    def structure(m: int, k: int):
        try:
            return sim.structure(m, k)
        except (KeyError, FileNotFoundError) as exc:
            raise HTTPException(404, str(exc).strip("'\""))

    @app.post("/api/build")
    def build(body: dict = Body(...)):
        try:
            co = body.get("cohort")
            return sim.build(str(body["run_id"]), float(body["t"]), None if co is None else int(co),
                             str(body.get("which", "mode")), bool(body.get("full", False)))
        except (KeyError, ValueError, IndexError, TypeError) as exc:
            raise HTTPException(400, str(exc).strip("'\""))
        except RuntimeError as exc:
            raise HTTPException(500, str(exc))

    return app


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python -m qd_workflows.simulator")
    ap.add_argument("--thermo", required=True, help="thermo_merged.json of the nucleation workflow")
    ap.add_argument("--params", help="growth parameters (yaml); default: the placeholders of growth.PARAMS")
    ap.add_argument("--preset", required=True, help="material config: cif, native_elements, preset (charges, passivation)")
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--port", type=int, default=8772)
    ap.add_argument("--max-atoms", type=int, default=MAX_ATOMS)
    ap.add_argument("--ligand", choices=sorted(LIGANDS), help="X-type ligand of the grid (default: from its references)")
    ap.add_argument("--nucleation-only", action="store_true",
                    help="the nucleation panel alone: no growth run, build or classical twin (a grid whose own growth "
                         "parameters do not exist yet)")
    ap.add_argument("--m0", type=int, help="fixed nucleus size (default: the parameters' m0; needed with "
                                           "--nucleation-only when m0 is auto)")
    a = ap.parse_args(argv)
    import uvicorn
    params = yaml.safe_load(Path(a.params).read_text()) if a.params else copy.deepcopy(gr.PARAMS)
    cfg = yaml.safe_load(Path(a.preset).read_text())
    if a.nucleation_only and a.m0 is None and params["material"]["m0"] == "auto":
        ap.error("--nucleation-only needs --m0: the automatic hand-over uses the growth (facet) parameters")
    sim = Simulator(json.loads(Path(a.thermo).read_text()), params, cfg, a.max_atoms, Path(a.thermo).resolve().parent,
                    a.ligand, a.nucleation_only, a.m0)
    uvicorn.run(make_app(sim), host=a.host, port=a.port, log_level="warning")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
