"""
Classical twin of the growth model: what the textbook picture predicts for the same synthesis.

Spherical particles of one isotropic surface energy gamma; LaMer-type burst nucleation from classical nucleation
theory; growth and ripening by the Talapin-Sugimoto size law with Gibbs-Thomson solubility. The synthesis is the
growth model's (precursor conversion k_conv P C_X with its exact decay over a step, the syringe, the temperature
program); the material parameters are set by hand or mapped from the atomistic model (simulator.from_atomistic).

    S = C / C_inf(T),   C_inf(T) = C_inf(T_ref) e^(-dH_sol/k (1/T - 1/T_ref))
    nucleus of n units:  dG(n) = -n kT ln S + gamma a n^(2/3),   a = (36 pi)^(1/3) V_m^(2/3)
    n* = (2 gamma a / (3 kT ln S))^3,   dG* = 4 gamma^3 a^3 / (27 (kT ln S)^2),   r* = 2 gamma V_m / (kT ln S)
    J = Z f* C e^(-dG*/kT)   (Zeldovich Z = sqrt(dG* / (3 pi kT n*^2)), f* = monomer attachment to the nucleus,
                              diffusion in series with the surface reaction), or J0 e^(-dG*/kT) with J0 set
    dr/dt = V_m n_A (C - C_inf e^(2 gamma V_m / (r kT))) / (r / D + 1 / k_r)
    nuclei are born where dG has dropped by kT past n*; dots below n_min units dissolve into the pool.

Units: A, s, M (n_A = 6.022e-4 molecules per A^3 per M), eV; gamma in J/m^2 and D in cm^2/s at the interface.

    python -m qd_workflows.classical run [--params p.yaml] [--json out.json]
"""
from __future__ import annotations

import argparse
import copy
import json
import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import yaml
from scipy.optimize import brentq

from .grid_thermo import KB_EV, k_conv_at

NA = 6.02214076e-4                          # molecules per A^3 per M
J_M2 = 0.0624150907                         # eV/A^2 per J/m^2

PARAMS = {
    "conditions": {"T": 500.0, "c_P0": 0.1, "c_MX0": 0.2, "k_conv_500K": 1.0, "Ea_conv": 0.0, "t_end": 600.0,
                   "t_start": 1e-3, "n_out": 60, "c_M0": 1e-12},
    "material": {
        "gamma_J_m2": 0.2,                  # isotropic surface energy
        "Vm_A3": 6.202964 ** 3 / 4,         # A^3 per CdSe (zinc blende, MACE lattice)
        "C_inf_M": 1e-6,                    # bulk solubility of the monomer at T_ref
        "T_ref": 500.0,
        "dH_sol_eV": 0.0,                   # dissolution enthalpy: C_inf(T)
        "D_cm2_s": 1e-5,                    # monomer diffusion
        "k_r_A_s": 1e4,                     # surface reaction
        "J0_mode": "zeldovich",             # "zeldovich": J = Z f* C e^(-dG*/kT); "set": J = J0 e^(-dG*/kT)
        "J0_M_s": 1e10,
    },
    "numerics": {"dt_rel_max": 0.01, "birth_interval": 0.05, "birth_share": 0.05, "merge_tol": 0.01,
                 "max_cohorts": 64, "n_min": 2.0},
}


class Classical:
    def __init__(self, params: dict):
        self.P = copy.deepcopy(params)
        self.T = None
        self.set_temperature(float(self.P["conditions"]["T"]))

    # -- the synthesis program (shared with the growth model) --
    def _prog(self):
        return SimpleNamespace(P=self.P)

    def temperature_at(self, t):
        from .growth import Growth
        return Growth.temperature_at(self._prog(), t)

    def temperature_breaks(self):
        from .growth import Growth
        return Growth.temperature_breaks(self._prog())

    def injection(self):
        from .growth import Growth
        return Growth.injection(self._prog())

    def set_temperature(self, T: float):
        c, m = self.P["conditions"], self.P["material"]
        self.T = float(T)
        self.kt = KB_EV * self.T
        self.k_conv = k_conv_at(self.T, float(c["k_conv_500K"]), float(c["Ea_conv"]))
        self.c_inf = float(m["C_inf_M"]) * math.exp(-float(m.get("dH_sol_eV", 0.0)) / KB_EV
                                                    * (1.0 / self.T - 1.0 / float(m.get("T_ref", self.T))))
        self.gam = float(m["gamma_J_m2"]) * J_M2                                    # eV/A^2
        self.Vm = float(m["Vm_A3"])
        self.D = float(m["D_cm2_s"]) * 1e16                                         # A^2/s
        self.k_r = float(m["k_r_A_s"])
        self.a = (36.0 * math.pi) ** (1.0 / 3.0) * self.Vm ** (2.0 / 3.0)           # A^2 per unit^(2/3)

    def radius(self, n):
        return (3.0 * np.asarray(n, float) * self.Vm / (4.0 * math.pi)) ** (1.0 / 3.0)

    # -- nucleation --
    def cnt(self, C: float) -> dict:
        """Classical nucleation at monomer level C (M): S, n*, r* (A), dG* (eV), J (M/s) and the birth size."""
        S = C / self.c_inf
        n_min = float(self.P["numerics"].get("n_min", 2.0))
        if S <= 1.0:
            return {"S": S, "n_star": None, "r_star": None, "dG_star": math.inf, "J": 0.0, "n_birth": n_min}
        L = self.kt * math.log(S)
        ga = self.gam * self.a
        n_star = (2.0 * ga / (3.0 * L)) ** 3
        dG_star = 4.0 * ga ** 3 / (27.0 * L * L)
        r_star = 2.0 * self.gam * self.Vm / L

        def dG(n):
            return -n * L + ga * n ** (2.0 / 3.0)

        lo, hi = n_star, max(2.0 * n_star, n_star + 2.0)                             # born kT below the top
        while dG(hi) > dG_star - self.kt:
            hi *= 2.0
        n_b = brentq(lambda n: dG(n) - (dG_star - self.kt), lo, hi) if dG(lo) > dG_star - self.kt else lo
        m = self.P["material"]
        x = min(dG_star / self.kt, 700.0)
        if m.get("J0_mode", "zeldovich") == "set":
            J = float(m["J0_M_s"]) * math.exp(-x)
        else:
            Z = math.sqrt(dG_star / (3.0 * math.pi * self.kt * n_star * n_star)) if n_star > 0 else 0.0
            f_star = 4.0 * math.pi * r_star ** 2 * NA * C / (r_star / self.D + 1.0 / self.k_r)     # 1/s
            J = Z * f_star * C * math.exp(-x)
        return {"S": S, "n_star": n_star, "r_star": r_star, "dG_star": dG_star, "J": J, "n_birth": max(n_b, n_min)}

    # -- growth --
    def rate(self, n, C: float):
        """dn/dt (units/s) of dots of n units at monomer level C: diffusion in series with the surface reaction,
        Gibbs-Thomson solubility C_inf e^(2 gamma V_m / (r kT))."""
        r = self.radius(n)
        gt = np.exp(np.minimum(2.0 * self.gam * self.Vm / (np.maximum(r, 1e-6) * self.kt), 700.0))
        flux = NA * (C - self.c_inf * gt) / (r / self.D + 1.0 / self.k_r)          # molecules / A^2 / s
        return 4.0 * math.pi * r * r * flux

    def _merge(self, co: list, n_b: float):
        num = self.P["numerics"]
        tol, cap = float(num.get("merge_tol", 0.01)), int(num.get("max_cohorts", 64))
        m_min = 4.0 * n_b                                    # near the birth size a cohort's fate hangs on small dn
        while True:
            big = sorted((c for c in co if c["n"] > m_min), key=lambda c: c["n"])
            if len(big) < 2:
                return
            gaps = [math.log(big[i + 1]["n"] / big[i]["n"]) for i in range(len(big) - 1)]
            i = int(np.argmin(gaps))
            if gaps[i] > tol and len(co) <= cap:
                return
            a, b = big[i], big[i + 1]
            co.remove(b)
            N = a["N"] + b["N"]
            a.update(n=(a["N"] * a["n"] + b["N"] * b["n"]) / N, N=N, t_b=min(a["t_b"], b["t_b"]))

    # -- integration --
    def run(self, seed=None) -> dict:
        """Cohorts of dots in n, the monomer pool implicit each step (nucleation and growth at the solved C), the
        precursor decaying exactly over the step, the syringe diluting the flask (amounts conserved exactly)."""
        c, num = self.P["conditions"], self.P["numerics"]
        P_, C, C_X = float(c["c_P0"]), float(c.get("c_M0", 1e-12)), float(c["c_MX0"])
        t, t_end = 0.0, float(c["t_end"])
        outs = list(np.geomspace(float(c.get("t_start", 1e-3)), t_end, int(c.get("n_out", 60))))
        inj = self.injection()
        V0 = inj["V0"] if inj else 1.0
        V, n_inj, n_conv = V0, 0.0, 0.0
        co = [{"N": float(s["N"]), "n": float(s["n"]), "t_b": 0.0} for s in (seed or [])]
        seeded = sum(x["N"] * x["n"] for x in co)
        pend = pend_n = 0.0
        t_last = 0.0
        breaks = self.temperature_breaks() + ([inj["t0"], inj["t0"] + inj["duration"]] if inj else [])
        dt_rel = float(num.get("dt_rel_max", 0.01))
        n_min = float(num.get("n_min", 2.0))
        series, snaps, self.steps = [], [], 0
        nuc = self.cnt(C)
        while outs:
            T_now = self.temperature_at(t)
            if abs(T_now - self.T) > 1.0 or (T_now != self.T and any(abs(t - b) < 1e-12 for b in breaks)):
                self.set_temperature(T_now)
            Ns = np.array([x["N"] for x in co])
            ns = np.array([x["n"] for x in co])
            q_inst = self.k_conv * P_ * C_X
            dt = outs[0] - t
            if t > 0:
                dt = min(dt, dt_rel * t)
            if len(co):
                r0 = self.rate(ns, C)
                with np.errstate(divide="ignore"):
                    g, d = r0 > 0, r0 < 0
                    if g.any():
                        dt = min(dt, float((0.02 * np.maximum(ns[g], 1.0) / r0[g]).min()))
                    if d.any():
                        dt = min(dt, float((0.2 * np.maximum(ns[d], 1.0) / -r0[d]).min()))
            if q_inst > 0:
                dt = min(dt, 0.02 * C_X / q_inst)
            nxt = [b for b in breaks if b > t + 1e-12]
            if nxt:
                dt = min(dt, nxt[0] - t)
            dTdt = abs(self.temperature_at(t + 1e-6) - self.temperature_at(t)) / 1e-6
            if dTdt > 0:
                dt = min(dt, 1.0 / dTdt)
            injecting = bool(inj) and inj["t0"] - 1e-12 <= t < inj["t0"] + inj["duration"] - 1e-12
            if injecting:
                dt = min(dt, 0.02 * V / inj["rate"])
            while True:
                Q = P_ * -math.expm1(-self.k_conv * C_X * dt) / dt if dt > 0 else 0.0

                def balance(lnC):
                    Cc = math.exp(lnC)
                    nc = self.cnt(Cc)
                    U = float((Ns * self.rate(ns, Cc)).sum()) if len(co) else 0.0
                    return Cc - C - dt * (Q - nc["n_birth"] * nc["J"] - U)

                lo, hi = math.log(1e-30), math.log(C + Q * dt + 1e-12)
                b_lo, b_hi = balance(lo), balance(hi)
                for _ in range(200):
                    if b_hi >= 0:
                        break
                    hi += math.log(2.0)
                    b_hi = balance(hi)
                lnC = brentq(balance, lo, hi, xtol=1e-12, rtol=1e-12) if b_lo < 0 < b_hi else (hi if b_hi < 0 else lo)
                Cs = math.exp(lnC)
                nuc = self.cnt(Cs)
                n_new = np.maximum(ns + self.rate(ns, Cs) * dt, 0.0) if len(co) else ns
                up = float((Ns * (n_new - ns)).sum()) if len(co) else 0.0
                C_new = C + (Q - nuc["n_birth"] * nuc["J"]) * dt - up
                if C_new > 0 or dt < 1e-12:
                    break
                dt *= 0.5
            for x, nn in zip(co, n_new):
                x["n"] = float(nn)
            P_ -= Q * dt
            C_X = max(C_X - Q * dt, 1e-30)
            C = max(C_new, 1e-30)
            n_conv += Q * dt * V
            pend += nuc["J"] * dt
            pend_n += nuc["n_birth"] * nuc["J"] * dt
            for x in [x for x in co if x["n"] < n_min]:                       # dissolved dots give back
                C += x["N"] * x["n"]
                co.remove(x)
            if injecting:
                dV = inj["rate"] * dt
                f = V / (V + dV)
                P_, C, C_X = P_ * f + inj["c"] * dV / (V + dV), C * f, C_X * f
                pend, pend_n = pend * f, pend_n * f
                for x in co:
                    x["N"] *= f
                n_inj += inj["c"] * dV
                V += dV
            t += dt
            self.steps += 1
            at_out = t >= outs[0] - 1e-12
            n_dots = sum(x["N"] for x in co)
            if at_out or t - t_last >= float(num.get("birth_interval", 0.05)) * t \
                    or pend >= float(num.get("birth_share", 0.05)) * (n_dots + pend):
                if pend > 1e-15:
                    co.append({"N": pend, "n": pend_n / pend, "t_b": t})
                    pend = pend_n = 0.0
                    t_last = t
                    self._merge(co, nuc["n_birth"])
            if at_out:
                N = sum(x["N"] for x in co)
                d = np.array([2.0 * float(self.radius(x["n"])) / 10.0 for x in co])
                w = np.array([x["N"] for x in co])
                md = float((w * d).sum() / N) if N > 0 else 0.0
                sd = float(np.sqrt((w * (d - md) ** 2).sum() / N) / md) if N > 0 and md > 0 else 0.0
                mm = float(sum(x["N"] * x["n"] for x in co) / N) if N > 0 else 0.0
                series.append({"t": t, "T": self.T, "V_L": V, "P": P_, "C_M": C, "C_X": C_X, "C_inf": self.c_inf,
                               "S": C / self.c_inf, "J": nuc["J"],
                               "barrier": nuc["dG_star"] if math.isfinite(nuc["dG_star"]) else None,
                               "n_star": nuc["n_star"], "r_star_nm": nuc["r_star"] / 10.0 if nuc["r_star"] else None,
                               "n_birth": nuc["n_birth"], "N_dots": N, "m_mean": mm, "d_mean_nm": md, "d_sd_rel": sd,
                               "n_inj_mol": n_inj, "n_conv_mol": n_conv})
                snaps.append({"t": t, "cohorts": [{"N": x["N"], "n": x["n"], "t_b": x["t_b"],
                                                   "d_nm": 2.0 * float(self.radius(x["n"])) / 10.0} for x in co]})
                outs.pop(0)
        return {"params": self.P, "series": series, "snapshots": snaps, "seeded": seeded, "V0_L": V0}


def conservation(res) -> list:
    """Total CdSe (precursor + monomer + dots, less what the syringe added) per litre of the starting flask: 0."""
    c = res["params"]["conditions"]
    V0 = res.get("V0_L", 1.0)
    out = []
    for row, snap in zip(res["series"], res["snapshots"]):
        dots = sum(x["N"] * x["n"] for x in snap["cohorts"])
        out.append({"t": row["t"], "MA": (row["V_L"] * (row["P"] + row["C_M"] + dots) - row["n_inj_mol"]) / V0
                    - c["c_P0"] - c.get("c_M0", 1e-12) - res.get("seeded", 0.0)})
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python -m qd_workflows.classical")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("run", help="run the classical model and print the final dots")
    p.add_argument("--params", help="yaml with conditions / material / numerics overrides")
    p.add_argument("--json", help="write the result here")
    a = ap.parse_args(argv)
    P = copy.deepcopy(PARAMS)
    if a.params:
        for k, v in (yaml.safe_load(Path(a.params).read_text()) or {}).items():
            P.setdefault(k, {}).update(v)
    res = Classical(P).run()
    last = res["series"][-1]
    print(f"t {last['t']:.4g} s: {last['N_dots']:.3g} M dots, <d> {last['d_mean_nm']:.2f} nm, "
          f"sd/<d> {last['d_sd_rel']:.2f}, S {last['S']:.3g}")
    if a.json:
        Path(a.json).write_text(json.dumps(res))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
