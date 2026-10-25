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

import numpy as np
import yaml

from . import growth as gr
from .grid_thermo import KB_EV, converted, k_conv_at, pool, thermo_at

PAGE = Path(__file__).parent / "templates" / "simulator" / "index.html"
LIMITS = {                                  # slider ranges of the page (the grid thermo spans its own T range)
    "T": [300.0, 700.0], "c_P0": [0.005, 0.5], "c_MX0": [0.01, 1.0], "k_conv_500K": [0.01, 100.0],
    "t_end": [1.0, 3600.0], "Ea_conv": [0.0, 2.0], "t": [1e-3, 1e4], "eps": [1.0, 40.0],
    "k_r": [1.0, 1e10],                     # A/s, surface reaction of the size law (1e10: the diffusion limit)
}


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


class Simulator:
    def __init__(self, thermo: dict, params: dict, preset_cfg: dict, max_atoms: int = MAX_ATOMS,
                 thermo_dir: str | Path | None = None):
        self.thermo, self.params, self.cfg = thermo, params, preset_cfg
        self.thermo_dir = Path(thermo_dir) if thermo_dir else None      # the grid's structure files are relative to it
        self.max_atoms = max_atoms
        self.runs: OrderedDict[str, dict] = OrderedDict()
        self._grids: dict = {}
        self._handovers: dict = {}
        self.build_lock = threading.Lock()  # the builder redirects stdout: one build at a time
        t = thermo["temperatures"]
        self.limits = dict(LIMITS, T=[max(LIMITS["T"][0], min(t)), min(LIMITS["T"][1], max(t))])

    def defaults(self) -> dict:
        c = self.params["conditions"]
        F = self.params["families"]
        eps = float(self.params.get("solvation", {}).get("eps", 1.0))
        k_r = self.params["kinetics"].get("k_r")
        return {"conditions": {**{k: c[k] for k in SETTINGS}, "t": 10.0, "eps": eps,
                               "k_r": float(k_r) if k_r is not None else LIMITS["k_r"][1]},
                "m0": self._m0({**c, "eps": eps})[0],
                "solvation_note": self.params.get("solvation", {}).get("note", ""),
                "limits": self.limits, "families": list(F), "max_atoms": self.max_atoms,
                "material": self.cfg.get("material", ""), "phase": self.cfg.get("phase", ""),
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

    def run(self, cond: dict) -> dict:
        P = copy.deepcopy(self.params)
        for k in SETTINGS:
            if k in cond and cond[k] is not None:
                lo, hi = self.limits[k]
                P["conditions"][k] = min(max(float(cond[k]), lo), hi)
        if cond.get("eps") is not None and "solvation" in P:
            P["solvation"]["eps"] = solvent_eps(cond["eps"])
        if cond.get("k_r") is not None:
            lo, hi = self.limits["k_r"]
            P["kinetics"]["k_r"] = min(max(float(cond["k_r"]), lo), hi)
        if P["conditions"]["c_MX0"] <= P["conditions"]["c_P0"]:
            raise ValueError("MX_q must exceed the precursor (every conversion takes one MX_q)")
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
        keep = ("t", "P", "C_M", "C_X", "c_MX_free", "S_bulk", "barrier", "J", "J0", "k_star", "N_dots", "m_mean",
                "m_std", "theta", "theta_eq", "gamma")
        used = res["params"]["families"]                     # gamma0 after the match_m0 rescaling, when it applied
        return {"run_id": rid, "seconds": round(time.time() - t0, 2), "conditions": P["conditions"],
                "match_m0": res["match_m0"], "modified": modified,
                "surface_used": {f: {"gamma0": used[f]["gamma0"], "eps0": used[f].get("eps0"),
                                     "dG_bind": used[f]["dG_bind"]} for f in used},
                "dG_thermal": g.dg_th, "m0": res["m0"], "handover": res["handover"],
                "series": [{k: _finite(r[k]) for k in keep} for r in res["series"]],
                "snapshots": snaps, "conservation": gr.conservation(res)[-1]}

    def build(self, rid: str, t: float, cohort: int | None, which: str = "mode", full: bool = False) -> dict:
        """Build the dot with QD_builder's target mode: the construction centres only (about a second), or with
        full=True also the layer roundings h_f + r delta_f, r in {-1, 0, 1} (closer to m, ten times slower)."""
        from builder.scripts.generate_library import _build
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

    def _omega(self, T: float, C_M: float, C_X: float, m0: int | None = None, eps: float = 1.0) -> dict:
        """Omega(m, k) = G(m, k) - m mu_MA - k mu_MX of the grid at the pool (C_M monomer supply, C_X MX_q left),
        Omega(m) = -kT ln sum_k e^(-Omega(m,k)/kT), the barrier max_{m <= m0} Omega(m) - Omega(1), m*, k*(m),
        P(k|m) and J (the nucleation step of the growth model); every species in implicit solvent at eps."""
        gr_ = self._grid(T, eps)
        th, ms, ks, Gp = gr_["th"], gr_["ms"], gr_["ks"], gr_["Gp"]
        kt = KB_EV * T
        pl = pool(th, max(C_M, 1e-30), max(C_X, 1e-30))
        mu_ma, mu_mx = pl["mu_MA"], pl["mu_MX"]
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
            rows.append({"m": m, "Omega_m": float(om_m[i]), "k_star": int(ks[j]), "Omega_min": float(v[i, j]),
                         "mean_k": float((P[i] * ks).sum()), "formula": node["formula"],
                         "source": node["states"][0].get("source", "grid")})
        below = [r for r in rows if r["m"] <= m0]
        top = max(below, key=lambda r: r["Omega_m"])
        barrier = max(0.0, top["Omega_m"] - below[0]["Omega_m"])
        k0 = next((r["k_star"] for r in rows if r["m"] == m0), None)
        return {"T": T, "C_M": C_M, "C_X": C_X, "mu_MA": mu_ma, "mu_MX": mu_mx, "c_MX_free": pl["c_MX_free"],
                "mean_n": pl["mean_n"], "speciation": {str(n): float(x) for n, x in pl["speciation"].items()},
                "m0": m0, "barrier": barrier, "m_star": top["m"], "k_star_m0": k0,
                "J": (J := gr.nucleation_flux(om_m, ms, m0, kt, float(self.params["material"]["Omega_m"]), self.params, barrier)),
                "J0": J * math.exp(min(barrier / kt, 700.0)),
                "ms": ms.tolist(), "ks": ks.tolist(), "Omega": [[None if not np.isfinite(x) else float(x) for x in r] for r in v],
                "P_k": P.round(6).tolist(), "per_m": rows}

    def omega_live(self, cond: dict, n_curve: int = 41) -> dict:
        """Omega at time t after injection from the precursor conversion alone (no nuclei yet: the pre-nucleation
        pool, the qdw nucleation Omega page), and the barrier along t on a log grid for the same conditions."""
        c = {k: float(cond.get(k, self.params["conditions"][k])) for k in SETTINGS if k != "t_end"}
        t = float(cond.get("t", 10.0))
        if c["c_MX0"] <= c["c_P0"]:
            raise ValueError("CdCl2 must exceed the precursor (every conversion takes one MX_q)")
        k = k_conv_at(c["T"], c["k_conv_500K"], c["Ea_conv"])
        c["eps"] = solvent_eps(cond.get("eps", self.params.get("solvation", {}).get("eps", 1.0)))
        m0, handover = self._m0(c)

        def at(tt):
            x = converted(tt, k, c["c_P0"], c["c_MX0"])
            return self._omega(c["T"], x, c["c_MX0"] - x, m0, c["eps"])

        out = at(t)
        out.update(t=t, alpha=out["C_M"] / c["c_P0"], conditions=c, handover=handover)
        curve, heat = [], []
        ms = [r["m"] for r in out["per_m"] if r["m"] <= out["m0"]]
        for tt in np.geomspace(1e-3, 1e4, n_curve):
            try:
                o = at(tt)
            except ValueError:
                continue
            curve.append({"t": float(tt), "barrier": o["barrier"], "m_star": o["m_star"], "J": o["J"], "C_M": o["C_M"]})
            om1 = o["per_m"][0]["Omega_m"]
            heat.append([r["Omega_m"] - om1 for r in o["per_m"] if r["m"] <= out["m0"]])
        out["curve"] = curve
        # Delta Omega(m, t) = Omega(m) - Omega(1) along the time after injection: rows m = 1..m0, columns the curve's t
        out["heat"] = {"t": [c["t"] for c in curve], "m": ms, "dOmega": [list(col) for col in zip(*heat)]}
        return out

    def omega(self, rid: str, t: float) -> dict:
        """Omega at time t of a growth run: the run's pool (precursor converted, minus what nuclei and dots took)."""
        if rid not in self.runs:
            raise KeyError("run not found (the server keeps the last 8 runs): run the model again")
        res, g = self.runs[rid]["res"], self.runs[rid]["growth"]
        row = min(res["series"], key=lambda r: abs(r["t"] - t))
        return {**self._omega(g.T, row["C_M"], row["C_X"], g.m0, float(g.solv.get("eps", 1.0))), "t": row["t"]}

    def structure(self, m: int, k: int) -> dict:
        """xyz of the grid node's lowest state (the structure file named in the thermo, relative to its folder)."""
        node = next((n for n in self.thermo["nodes"] if n["m"] == m and n["k"] == k), None)
        if node is None:
            raise KeyError(f"no grid node m = {m}, k = {k}")
        st = node["states"][0]
        path = Path(st["file"]) if self.thermo_dir is None else self.thermo_dir / st["file"]
        if not path.exists():
            raise FileNotFoundError(f"structure file {st['file']} is not in {self.thermo_dir}")
        return {"m": m, "k": k, "formula": node["formula"], "E_eV": st["E_eV"], "source": st.get("source", "grid"),
                "n_states": len(node["states"]), "file": st["file"], "xyz": path.read_text()}


def make_app(sim: Simulator):
    from fastapi import Body, FastAPI, HTTPException
    from fastapi.responses import HTMLResponse

    app = FastAPI(title="QD nucleation and growth simulator")

    @app.get("/", response_class=HTMLResponse)
    def page():
        return PAGE.read_text()

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
    a = ap.parse_args(argv)
    import uvicorn
    params = yaml.safe_load(Path(a.params).read_text()) if a.params else copy.deepcopy(gr.PARAMS)
    cfg = yaml.safe_load(Path(a.preset).read_text())
    sim = Simulator(json.loads(Path(a.thermo).read_text()), params, cfg, a.max_atoms, Path(a.thermo).resolve().parent)
    uvicorn.run(make_app(sim), host=a.host, port=a.port, log_level="warning")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
