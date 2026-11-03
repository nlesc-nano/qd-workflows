"""
Network growth simulator: a local web page on top of the network growth model (network_growth.py).

The page sets the conditions (T, solvent, precursor and CdCl2 concentrations, conversion rate, the kinetic scale s_nu,
viscosity), shows step 1 as the nucleation simulator has it (the grid's Omega(m) up to the seed's m = 13 from the
precursor alone, simulator.Simulator.omega_live) and the grand-potential map of the network live at a time, runs the
model, and plots the pool, the dots, the size distribution, the growth of each facet family, the island barriers and
the most populated structures at any time.

    python -m qd_workflows.network_sim --network DIR --thermo thermo_merged.json --params cdse_zb.yaml [--port 8773]

Run it in an environment with FastAPI, uvicorn and scipy.
"""
from __future__ import annotations

import argparse
import heapq
import math
import threading
import uuid
from collections import OrderedDict
from pathlib import Path

import numpy as np

from . import network_growth as ng
from .grid_thermo import converted, k_conv_at

PAGE = Path(__file__).parent / "templates" / "network" / "index.html"
VOLUME_A3 = 59.667                       # zinc-blende CdSe volume per CdSe unit (a = 6.203 A): d = (6 m V / pi)^(1/3)
MAX_RUNS = 6
PATH_MAX = 500                           # states on a run's path (the growth movie)
KIND_NAME = {v: k for k, v in ng.KIND.items()}
DISSOLVE_MIN_N = 1e-4                    # dissolving dots: N >= this share of the most populated dot's


def diameter_nm(m):
    return (6.0 * np.asarray(m, float) * VOLUME_A3 / math.pi) ** (1 / 3) / 10.0


class NetSim:
    def __init__(self, net: ng.Network, params: dict | None = None):
        """params: the full growth template; with it, step 1 (the nucleation simulator's grid Omega, nucleus at the
        seed's m) on the same thermo, solvation and pool as the network."""
        self.net = net
        self.nuc = None
        if params is not None:
            from .simulator import Simulator
            self.nuc = Simulator(net.thermo, params, {}, thermo_dir=net.thermo_dir, nucleation_only=True, m0=ng.M_DOT)
            for k in ("T", "c_P0", "c_MX0", "k_conv_500K", "Ea_conv", "eps"):
                if k in ng.LIMITS:
                    self.nuc.limits[k] = list(ng.LIMITS[k])        # the network's ranges, not the old page's
        self.runs: OrderedDict = OrderedDict()
        self.jobs: dict = {}
        g = np.load(net.folder / "geometry.npz")
        self.geo = (g["symbols"], g["positions"], g["offsets"])
        self.elements = net.meta.get("elements", ("Cd", "Se", "Cl"))   # exports before the carboxylate networks: no list
        self.comp = {}                                           # (m, k) -> states
        for i, (m, k) in enumerate(zip(net.m.tolist(), net.k.tolist())):
            self.comp.setdefault((m, k), []).append(i)
        self.comp_keys = sorted(self.comp)
        self.comp_of = np.empty(net.n, int)
        for j, key in enumerate(self.comp_keys):
            self.comp_of[self.comp[key]] = j

    # ---------------- page data
    def defaults(self):
        net = self.net
        kinds = {name: int((net.kind == v).sum()) for name, v in ng.KIND.items()}
        return {"conditions": {k: v for k, v in ng.DEFAULTS.items()}, "limits": ng.LIMITS,
                "network": {"states": int(net.n), "grid_nodes": int(net.n_grid), "steps": int(len(net.a)),
                            "kinds": kinds, "island_curves": len(net.isl), "m_max": int(net.m.max()),
                            "d_max_nm": float(diameter_nm(net.m.max())), "seed": net.label[net.seed],
                            "compositions": len(self.comp_keys)},
                "nucleation": self.nuc is not None,
                "no_thermo": sorted([m, k] for m, k in self.nuc.no_thermo) if self.nuc else []}

    def _cond(self, body):
        c = {**ng.DEFAULTS}
        for k, v in (body or {}).items():
            if k in ng.DEFAULTS:
                c[k] = type(ng.DEFAULTS[k])(v) if not isinstance(ng.DEFAULTS[k], bool) else bool(v)
        return c

    def start(self, body):
        """Run in a thread (a run takes minutes on the full network); poll with job()."""
        c = self._cond(body.get("conditions"))
        jid = uuid.uuid4().hex[:10]
        job = {"state": "running", "log": "", "result": None, "error": None, "stop": False}
        self.jobs[jid] = job

        def work():
            try:
                job["result"] = self.run(c, log=lambda s: job.__setitem__("log", s), stop=lambda: job["stop"])
                job["state"] = "done"
            except Exception as exc:                            # reported to the page
                job["state"], job["error"] = "error", str(exc)
        threading.Thread(target=work, daemon=True).start()
        return {"job": jid}

    def stop(self, jid):
        """Ask a running job to end at its next step; it returns the outputs it has."""
        self.jobs[jid]["stop"] = True
        return {"job": jid}

    def job(self, jid):
        j = self.jobs[jid]
        out = {"state": j["state"], "log": j["log"], "error": j["error"]}
        if j["state"] == "done":
            out["result"] = j["result"]
            self.jobs.pop(jid, None)
        return out

    def run(self, c, log=None, stop=None):
        r = ng.simulate(self.net, c, log=log, stop=stop)
        rid = uuid.uuid4().hex[:10]
        self.runs[rid] = r
        while len(self.runs) > MAX_RUNS:
            self.runs.popitem(last=False)
        net = self.net
        dot = net.m >= ng.M_DOT
        ms = np.arange(ng.M_DOT, int(net.m.max()) + 1)
        psd = []
        top = []
        for N in r["N"]:
            h = np.bincount(net.m[dot] - ng.M_DOT, weights=N[dot], minlength=len(ms))
            psd.append(h.tolist())
            order = np.argsort(-N)[:8]
            top.append([[int(i), net.label[i], int(net.m[i]), int(net.k[i]), float(N[i])] for i in order if N[i] > 0])
        series = {k: r[k] for k in ("t", "alpha", "c_MA", "c_X_left", "c_MX_free", "mu_MA", "mu_MX", "n_dots",
                                    "n_clusters", "m_mean", "m_std", "nucleation", "frontier", "layer_flux",
                                    "cdse_flux", "attach_flux", "barriers")}
        series["d_mean_nm"] = diameter_nm(r["m_mean"]).tolist()
        return {"run_id": rid, "cond": c, "series": series, "psd": {"m": ms.tolist(), "d_nm": diameter_nm(ms).tolist(),
                                                                    "N": psd}, "top": top,
                "path_i": len(r["t"]) - 1, "omega_run": self._omega_run(r),
                "stats": {"steps": r["steps"], "lu": r["lu"], "wall_s": r["wall_s"], "stopped": r["stopped"]}}

    # ---------------- step 1: the nucleation simulator's grid Omega
    @staticmethod
    def _nuc_cond(c):
        return {k: c[k] for k in ("T", "c_P0", "c_MX0", "k_conv_500K", "Ea_conv", "eps")}

    def nucleation(self, body):
        """Omega(m) of the grid up to m = 13 at time t from the precursor alone, and the barrier along the time
        (simulator.Simulator.omega_live)."""
        if self.nuc is None:
            raise ValueError("started without the growth template: no nucleation panel")
        c = self._cond(body.get("conditions"))
        return self.nuc.omega_live({**self._nuc_cond(c), "t": float(body.get("t", 10.0))})

    def _omega_run(self, r):
        """The grid's Omega(m) at the pool a network run leaves at each output (CdSe in solution and CdCl2 left),
        as the nucleation simulator's growth-run curves: barrier, m*, S and Delta Omega(m, t)."""
        if self.nuc is None:
            return None
        c = r["cond"]
        out = {"t": [], "barrier": [], "m_star": [], "C_M": [], "C_X": [], "c_X_free": [], "S": [], "m": None,
               "dOmega": []}
        cols = []
        for t, cm, cx in zip(r["t"], r["c_MA"], r["c_X_left"]):
            o = self.nuc._omega_at(float(c["T"]), float(cm), float(cx), ng.M_DOT, float(c["eps"]), {})
            if o is None:
                continue
            rows = [q for q in o["per_m"] if q["m"] <= o["m0"]]
            out["t"].append(float(t))
            out["barrier"].append(o["barrier"])
            out["m_star"].append(o["m_star"])
            out["C_M"].append(o["C_M"])
            out["C_X"].append(o["C_X"])
            out["c_X_free"].append(o["c_MX_free"])
            out["S"].append(o["S"])
            out["m"] = [q["m"] for q in rows]
            cols.append([q["Omega_m"] - rows[0]["Omega_m"] for q in rows])
        out["dOmega"] = [list(col) for col in zip(*cols)]
        return out

    def profile(self, body):
        """Omega(m) = -kT ln sum_states e^(-Omega/kT) of every size, from the pool's complexes (m = 1) through the grid
        (m <= 13) to the network's dots, at the precursor-only pool at time t, or at the pool a run leaves at its
        output i; with a run, the population of each size as -kT ln(N(m)/c0), equal to Omega(m) where the size is in
        equilibrium with the pool and above it where the run has not filled it. Relative to Omega(1)."""
        net = self.net
        if body.get("run_id"):
            r = self.runs[str(body["run_id"])]
            i = int(body["i"])
            c = r["cond"]
            _g, th, _pa = net.G(c["T"], c["eps"])
            pl, N, t = ng._pool(th, r["c_MA"][i], r["c_X_left"][i]), r["N"][i].astype(float), r["t"][i]
        else:
            c = self._cond(body.get("conditions"))
            t = float(body.get("t", 10.0))
            x = converted(t, k_conv_at(c["T"], c["k_conv_500K"], c["Ea_conv"]), c["c_P0"], c["c_MX0"])
            _g, th, _pa = net.G(c["T"], c["eps"])
            pl, N = ng._pool(th, x, c["c_MX0"] - x), None
        kt = ng.KB_EV * c["T"]
        g = _g
        om = g - net.m * pl["mu_MA"] - net.k * pl["mu_MX"]
        om1 = np.array([th["G"][(1, n)] - pl["mu_MA"] - n * pl["mu_MX"] for n in net.mono_k])
        ref = float(om1.min() - kt * np.log(np.exp(-(om1 - om1.min()) / kt).sum()))
        ms = np.arange(2, int(net.m.max()) + 1)
        omm, low, pop = [], [], []
        for m in ms:
            sel = np.where(net.m == m)[0]
            if not len(sel):
                omm.append(None); low.append(None); pop.append(None)
                continue
            o = om[sel]
            omm.append(float(o.min() - kt * np.log(np.exp(-(o - o.min()) / kt).sum()) - ref))
            j = sel[np.argmin(o)]
            low.append([float(om[j] - ref), int(net.k[j]), net.label[j], "grid" if j < net.n_grid else "network"])
            if N is not None:
                n_m = float(N[sel].sum())
                pop.append(float(-kt * math.log(n_m) - ref) if n_m > 1e-30 else None)
        out = {"t": t, "m": [1] + ms.tolist(), "dOmega": [0.0] + omm, "lowest": [[0.0, None, "pool", "pool"]] + low,
               "mu_MA": pl["mu_MA"], "mu_MX": pl["mu_MX"], "c_MA": pl["c_MA"], "m_grid": ng.M_DOT, "run": N is not None}
        if N is not None:
            out["pop"] = [float(-kt * math.log(pl["c_MA"]) - ref)] + pop
            out["N"] = [pl["c_MA"]] + [None if q is None else math.exp(-(q + ref) / kt) for q in pop]
        return out

    def grid_structure(self, m, k):
        if self.nuc is None:
            raise ValueError("started without the growth template")
        return self.nuc.structure(m, k)

    def _pool_at(self, body):
        """The pool of a run's output i, or of the conditions with no dots (live) at time t from the start (the
        precursor converted by then), or at conversion alpha."""
        if body.get("run_id"):
            r = self.runs[str(body["run_id"])]
            i = int(body["i"])
            c = r["cond"]
            _g, th, _pa = self.net.G(c["T"], c["eps"])
            return c, ng._pool(th, r["c_MA"][i], r["c_X_left"][i]), r["N"][i], r["t"][i]
        c = self._cond(body.get("conditions"))
        if body.get("t") is not None:
            t = float(body["t"])
            x = converted(t, k_conv_at(c["T"], c["k_conv_500K"], c["Ea_conv"]), c["c_P0"], c["c_MX0"])
        else:
            t = None
            x = float(body.get("alpha", 0.2)) * c["c_P0"]
        if not 0 < x < c["c_P0"]:
            raise ValueError("the conversion must be in (0, 1)")
        _g, th, _pa = self.net.G(c["T"], c["eps"])
        return c, ng._pool(th, x, c["c_MX0"] - x), None, t

    def map(self, body):
        """Omega per composition (the lowest state; eV, from the seed's), the population per composition, and the
        steps carrying the most net flux (between compositions)."""
        net = self.net
        c, pl, N, t = self._pool_at(body)
        nu, nu_src = net.prefactors(c["T"], c["eta_mPas"], c["s_nu"], c["a_site_A2"])
        r = net.rates(c["T"], c["eps"], pl, nu, nu_src, c["islands"])
        om = r["Omega"]
        best = np.full(len(self.comp_keys), np.inf)
        np.minimum.at(best, self.comp_of, om)
        ref = float(om[net.seed])
        pop = np.zeros(len(self.comp_keys))
        routes = []
        if N is not None:
            np.add.at(pop, self.comp_of, N)
            J, _js = net.fluxes(r, N)
            ca, cb = self.comp_of[net.a], self.comp_of[net.b]
            sel = ca != cb
            agg = {}
            for x, y, j in zip(ca[sel], cb[sel], J[sel]):
                agg[(x, y)] = agg.get((x, y), 0.0) + j
            tot = sum(abs(v) for v in agg.values()) or 1.0
            for (x, y), j in sorted(agg.items(), key=lambda kv: -abs(kv[1]))[:150]:
                if abs(j) < 1e-6 * tot:
                    break
                (m1, k1), (m2, k2) = self.comp_keys[x], self.comp_keys[y]
                routes.append([m1, k1, m2, k2, float(j)])
        cells = [[m, k, float(best[j] - ref), float(pop[j]), len(self.comp[(m, k)])]
                 for j, (m, k) in enumerate(self.comp_keys)]
        # the m = 1 row: the pool's complexes CdSe(CdCl2)_n (not states: the pool), Omega at the same mu, their
        # concentration c_MA x_n during a run
        _g, th, _pa = net.G(c["T"], c["eps"])
        for n in net.mono_k:
            om1 = th["G"][(1, n)] - pl["mu_MA"] - n * pl["mu_MX"]
            cells.append([1, n, float(om1 - ref), float(pl["c_MA"] * pl["speciation"].get(n, 0.0)) if N is not None
                          else 0.0, 1])
        return {"t": t, "alpha": float(pl["c_MA"] / c["c_P0"]) if N is None else None, "mu_MA": pl["mu_MA"], "mu_MX": pl["mu_MX"], "c_MA": pl["c_MA"], "c_MX_free": pl["c_MX_free"],
                "seed": net.label[net.seed], "cells": cells, "routes": routes, "barriers": r["table"]}

    def _transfer(self, run_id):
        """Net transfer of every step and nucleus source up to each output (M; the fluxes at the outputs integrated
        by the trapezoid rule)."""
        r = self.runs[run_id]
        if "_transfer" not in r:
            J = [self.net.fluxes(self._rates_at(run_id, i), r["N"][i]) for i in range(len(r["t"]))]
            t = np.asarray(r["t"])
            acc, out = [np.zeros_like(J[0][0]), np.zeros_like(J[0][1])], []
            for i in range(len(t)):
                if i:
                    for q in range(2):
                        acc[q] = acc[q] + 0.5 * (J[i][q] + J[i - 1][q]) * (t[i] - t[i - 1])
                out.append((acc[0].copy(), acc[1].copy()))
            r["_transfer"] = out
        return r["_transfer"]

    def _rates_at(self, run_id, i):
        c, pl, _N, _t = self._pool_at({"run_id": run_id, "i": i})
        nu, nu_src = self.net.prefactors(c["T"], c["eta_mPas"], c["s_nu"], c["a_site_A2"])
        return self.net.rates(c["T"], c["eps"], pl, nu, nu_src, c["islands"])

    def _edges(self):
        """Per state, its steps as (step, other state, sign): up (to the larger state, sign +1) and down (sign -1);
        isomerisations both ways, in both lists."""
        if not hasattr(self, "_up"):
            net = self.net
            up, down = [[] for _ in range(net.n)], [[] for _ in range(net.n)]
            for e, (a, b) in enumerate(zip(net.a.tolist(), net.b.tolist())):
                up[a].append((e, b, 1.0)); down[b].append((e, a, -1.0))
                if net.dm[e] == 0 and net.dk[e] == 0:
                    up[b].append((e, a, -1.0)); down[a].append((e, b, 1.0))
            self._up, self._down = up, down
        return self._up, self._down

    @staticmethod
    def _widest(start, edges, F):
        """Max-bottleneck paths from start {state: width} over edges (a step's width: sign * F[step], when > 0):
        the width and the parent (previous state, step) of every state reached."""
        width = {}
        par = {}
        heap = [(-w, x) for x, w in start.items()]
        heapq.heapify(heap)
        best = dict(start)
        while heap:
            w, x = heapq.heappop(heap)
            if x in width:
                continue
            width[x] = -w
            for e, y, sign in edges[x]:
                f = min(-w, sign * F[e])
                if f > 0 and y not in width and f > best.get(y, 0.0):
                    best[y], par[y] = f, (x, e, sign * F[e])
                    heapq.heappush(heap, (-f, y))
        return width, par

    def path(self, run_id, i, state=None, mode="growth"):
        """A path of the run at output i, frames: state, label, m, k, the step into it (its net transfer in M, or its
        net flux in M/s), the step's kind, N (M) at output i.
        growth: from the pool to state (default the most populated dot, or state before there are dots), the chain
        of steps along which the most material has moved so far, i.e. whose narrowest step has the largest net
        transfer (the flux integrated over the run; steps to the larger state, isomerisations either way, the
        nucleus sources first).
        dissolution: from state (default the dot losing the most to smaller states) down the steps carrying net flux
        to the smaller state at output i, the chain with the largest narrowest flux, to the pool when the flux reaches
        it (a nucleus going back to the pool's complexes), else to the smallest state it reaches."""
        net = self.net
        i = int(i)
        r = self.runs[run_id]
        N = np.asarray(r["N"][i], float)
        t = float(r["t"][i])
        up, down = self._edges()
        frame = lambda x, f, how: [x, net.label[x], int(net.m[x]), int(net.k[x]), float(f), how, float(N[x])]  # noqa: E731
        if mode == "dissolution":
            J, js = net.fluxes(self._rates_at(run_id, i), N)
            if state is None:
                lst = self.dissolving(run_id, i)["states"]
                if not lst:
                    return {"t": t, "i": i, "mode": mode, "frames": []}
                state = lst[0][0]
            state = int(state)
            width, par = self._widest({state: np.inf}, down, J)
            out = {}                                                     # nuclei with net flux back to the pool
            for q, (_n1, _n2, j) in enumerate(net.src.tolist()):
                if j in width and js[q] < 0:
                    out[j] = out.get(j, 0.0) - js[q]
            if out:
                end = max(out, key=lambda j: min(width[j], out[j]))
            else:
                end = min(width, key=lambda x: (net.m[x], -width[x] if x != state else 0.0))
            frames = []
            cur = end
            while cur != state and len(frames) < PATH_MAX:
                prev, e, f = par[cur]
                frames.append(frame(cur, f, KIND_NAME[int(net.kind[e])]))
                cur = prev
            frames.append(frame(state, 0.0, "start"))
            frames = frames[::-1]
            if out:
                frames.append([-1, "the pool", 1, int(net.src[net.src[:, 2] == end][0][1]), float(out[end]), "to the pool",
                               float(r["c_MA"][i])])
            return {"t": t, "i": i, "mode": mode, "frames": frames}
        if N.max() <= 0:
            return {"t": t, "i": i, "mode": mode, "frames": []}
        X, xs = self._transfer(run_id)[i]
        start = {}
        for q, (_n1, _n2, j) in enumerate(net.src.tolist()):                # the pool's nuclei
            if xs[q] > start.get(j, 0.0):
                start[j] = float(xs[q])
        width, par = self._widest(start, up, X)
        if state is None:
            reach = np.zeros(net.n, bool); reach[list(width)] = True
            dot = net.m >= ng.M_DOT
            cand = dot & reach if (dot & reach & (N > 0)).any() else reach
            if not cand.any():
                return {"t": t, "i": i, "mode": mode, "frames": []}
            state = int(np.argmax(np.where(cand, N, -1.0)))
        state = int(state)
        if state not in width:
            raise ValueError(f"nothing has moved from the pool to {net.label[state]} by t = {t:.3g} s")
        frames, cur = [], state
        while cur in par and len(frames) < PATH_MAX:
            prev, e, f = par[cur]
            frames.append(frame(cur, f, KIND_NAME[int(net.kind[e])]))
            cur = prev
        frames.append(frame(cur, start[cur], "nucleus"))
        return {"t": t, "i": i, "mode": mode, "frames": frames[::-1]}

    def dissolving(self, run_id, i, n=8):
        """The dots shrinking most at output i: the CdSe they lose per second to smaller states (steps that add CdSe,
        run backwards), among the dots holding >= DISSOLVE_MIN_N of the most populated dot's N whose loss exceeds
        their growth; with the CdSe they gain from steps up for comparison (M/s of CdSe)."""
        net = self.net
        i = int(i)
        N = np.asarray(self.runs[run_id]["N"][i], float)
        J, _js = net.fluxes(self._rates_at(run_id, i), N)
        g = net.dm > 0
        w = J[g] * net.dm[g]
        dn = np.bincount(net.b[g], weights=np.maximum(-w, 0.0), minlength=net.n)
        upf = np.bincount(net.a[g], weights=np.maximum(w, 0.0), minlength=net.n)
        dot = net.m >= ng.M_DOT
        big = N >= DISSOLVE_MIN_N * (N[dot].max() if dot.any() else 0.0)
        score = np.where(dot & big & (dn > upf), dn, 0.0)
        order = [int(x) for x in np.argsort(-score)[:n] if score[x] > 0]
        return {"t": float(self.runs[run_id]["t"][i]), "i": i,
                "states": [[x, net.label[x], int(net.m[x]), int(net.k[x]), float(N[x]), float(dn[x]), float(upf[x])]
                           for x in order]}

    def composition(self, m, k, run_id=None, i=None):
        """The states of one composition (lowest G first) with their populations at a run's output."""
        net = self.net
        ids = self.comp.get((m, k))
        if m == 1 and not ids:
            return {"m": m, "k": k, "states": [], "pool": True}
        if not ids:
            raise KeyError(f"no state with m = {m}, k = {k}")
        N = self.runs[run_id]["N"][int(i)] if run_id else None
        rows = [[int(s), net.label[s], int(net.reg_id[s]), None if N is None else float(N[s])] for s in ids]
        return {"m": m, "k": k, "states": rows}

    def structure(self, state: int) -> dict:
        net = self.net
        if not 0 <= state < net.n:
            raise KeyError(f"no state {state}")
        if state < net.n_grid:
            m, k = net.grid_keys[state]
            node = next(n for n in net.thermo["nodes"] if n["m"] == m and n["k"] == k)
            st = min(node["states"], key=lambda s: s["E_eV"])
            path = Path(st["file"]) if net.thermo_dir is None else net.thermo_dir / st["file"]
            if not path.exists():
                raise FileNotFoundError(f"structure file {st['file']} is not in {net.thermo_dir}")
            return {"state": state, "label": net.label[state], "m": m, "k": k, "source": "grid", "xyz": path.read_text()}
        row = int(net.row[state])
        sym, pos, off = self.geo
        z, x = sym[off[row]:off[row + 1]], pos[off[row]:off[row + 1]]
        lines = [str(len(z)), f"{net.label[state]} registry id {int(net.reg_id[state])}"]
        lines += [f"{self.elements[int(a)]} {p[0]:.4f} {p[1]:.4f} {p[2]:.4f}" for a, p in zip(z, x)]
        return {"state": state, "label": net.label[state], "m": int(net.m[state]), "k": int(net.k[state]),
                "source": "network", "registry_id": int(net.reg_id[state]), "xyz": "\n".join(lines) + "\n"}


def make_app(sim: NetSim):
    from fastapi import Body, FastAPI, HTTPException
    from fastapi.responses import HTMLResponse

    app = FastAPI(title="QD network growth simulator")

    def guard(fn, *args):
        try:
            return fn(*args)
        except (KeyError, ValueError, TypeError, IndexError) as exc:
            raise HTTPException(400, str(exc).strip("'\""))
        except FileNotFoundError as exc:
            raise HTTPException(404, str(exc))
        except RuntimeError as exc:
            raise HTTPException(500, str(exc))

    @app.get("/", response_class=HTMLResponse)
    def page():
        return PAGE.read_text()

    @app.get("/api/defaults")
    def defaults():
        return sim.defaults()

    @app.post("/api/run")
    def run(body: dict = Body(...)):
        return guard(sim.start, body)

    @app.post("/api/stop")
    def stop(id: str):
        return guard(sim.stop, id)

    @app.get("/api/job")
    def job(id: str):
        return guard(sim.job, id)

    @app.post("/api/map")
    def map_(body: dict = Body(...)):
        return guard(sim.map, body)

    @app.post("/api/nucleation")
    def nucleation(body: dict = Body(...)):
        return guard(sim.nucleation, body)

    @app.post("/api/profile")
    def profile(body: dict = Body(...)):
        return guard(sim.profile, body)

    @app.get("/api/grid_structure")
    def grid_structure(m: int, k: int):
        return guard(sim.grid_structure, m, k)

    @app.get("/api/path")
    def path(run_id: str, i: int, state: int | None = None, mode: str = "growth"):
        return guard(sim.path, run_id, i, state, mode)

    @app.get("/api/dissolving")
    def dissolving(run_id: str, i: int):
        return guard(sim.dissolving, run_id, i)

    @app.get("/api/composition")
    def composition(m: int, k: int, run_id: str | None = None, i: int | None = None):
        return guard(sim.composition, m, k, run_id, i)

    @app.get("/api/structure")
    def structure(state: int):
        return guard(sim.structure, state)

    return app


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python -m qd_workflows.network_sim")
    ap.add_argument("--network", required=True, help="network_export folder")
    ap.add_argument("--thermo", required=True, help="the grid's thermo_merged.json")
    ap.add_argument("--params", required=True, help="growth template with the solvation block")
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--port", type=int, default=8773)
    a = ap.parse_args(argv)
    import uvicorn
    import yaml
    sim = NetSim(ng.load(a.network, a.thermo, a.params), yaml.safe_load(Path(a.params).read_text()))
    uvicorn.run(make_app(sim), host=a.host, port=a.port, log_level="warning")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
