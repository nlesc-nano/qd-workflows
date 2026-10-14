"""Summary, static plots (matplotlib) and the interactive page (Plotly + plain JS) of a reaction run."""
from __future__ import annotations

import html
import json
import math
from pathlib import Path

import numpy as np

from . import model, solvation

KCAL = 23.060548
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
INK, INK2, INK3, GRID, SURF = "#0b0b0b", "#52514e", "#8a8984", "#e6e5e0", "#fcfcfb"
PARTS = [("dE_DFT", "ΔE DFT"), ("dZPE", "ΔZPE"), ("dH_thermal", "ΔH thermal"), ("-TdS_trans", "−TΔS trans"),
         ("-TdS_rot", "−TΔS rot"), ("-TdS_vib", "−TΔS vib"), ("-TdS_el", "−TΔS el"),
         ("ddG_solv_el", "ΔΔG solv, electrostatic"), ("ddG_solv_nonel", "ΔΔG solv, cavity+dispersion"),
         ("ligand", "ΔG lig (chain)")]


def _ti(res, t=298.15):
    return int(np.argmin(np.abs(np.asarray(res["T"]) - t)))


def _crossing(T, y):
    if not np.all(np.isfinite(y)):
        return None
    j = np.where(np.sign(y[:-1]) != np.sign(y[1:]))[0]
    if not j.size:
        return None
    i = int(j[0])
    return float(T[i] + (0 - y[i]) * (T[i + 1] - T[i]) / (y[i + 1] - y[i]))


def models(res):
    c = res["conditions"]
    cav = c.get("cavity", "beta0")
    return {
        "solution, COSMO": {"eps": c["eps"], "electrostatics": "cosmo", "cavity": cav},
        "solution, SCCS": {"eps": c["eps"], "electrostatics": "sccs", "cavity": cav},
        "gas phase": {"eps": 1.0, "gas": True},
    }


def summary(res: dict) -> dict:
    i = _ti(res)
    out = {"equation": res["equation"], "nu": res["nu"], "T_K": res["T"][i], "per": res.get("per"),
           "conditions": res["conditions"], "concentrations_M": {n: s["c"] for n, s in res["species"].items()},
           "models": {}, "flags": {n: s["flags"] for n, s in res["species"].items() if s["flags"]}}
    for name, cond in models(res).items():
        r = model.reaction(res, cond)
        T = np.asarray(r["T"])
        out["models"][name] = {
            "dG0_eV": float(r["dG0"][i]), "dG_eV": float(r["dG"][i]), "dH0_eV": float(r["dH0"][i]),
            "TdS0_eV": float(T[i] * r["dS0"][i]), "kTlnQ_eV": float(r["kTlnQ"][i]),
            "T_dG_zero_K": _crossing(T, r["dG"]),
            "parts_eV": {k: float(v[i]) for k, v in r["parts"].items()},
        }
    k = model.per_factor(res)
    out["dE_alternatives_eV_per"] = {lab: (None if v is None else v / k) for lab, v in (res.get("dE_alternatives_eV") or {}).items()}
    out["solvation_eV"] = {}
    for n, s in res["species"].items():
        eps = res["conditions"]["eps"]
        out["solvation_eV"][n] = {
            "eps": eps, "cosmo_el": solvation.electrostatic(s, eps, "cosmo"), "sccs_el": solvation.electrostatic(s, eps, "sccs"),
            "gb_el": solvation.electrostatic(s, eps, "gb"),
            "cavity_beta0": solvation.dg_solv(s, eps, "cosmo", "beta0") - solvation.electrostatic(s, eps, "cosmo"),
            "dipole_D": s.get("dipole_D"), "alpb": (s.get("alpb") or {}).get("dG_eV"),
        }
    out["energy_sets"] = res.get("energy_sets")
    return out


# --------------------------------------------------------------------------
# Static plots
# --------------------------------------------------------------------------

def _style(ax, xlabel, ylabel, title):
    ax.set_xlabel(xlabel, color=INK2)
    ax.set_ylabel(ylabel, color=INK2)
    ax.set_title(title, color=INK, loc="left", fontsize=11)
    ax.grid(color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK3)
    ax.tick_params(colors=INK2)


def plots(res: dict, out: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

    per = f" per {res['per']}" if res.get("per") else ""
    unit = f"ΔG (eV{per})"
    T = np.asarray(res["T"])
    ms = models(res)

    fig, ax = plt.subplots(figsize=(7, 4.2))
    for k, (name, cond) in enumerate(ms.items()):
        r = model.reaction(res, cond)
        ax.plot(T, r["dG"], color=SERIES[k], lw=2, label=f"{name}, at the given c")
        ax.plot(T, r["dG0"], color=SERIES[k], lw=1.2, ls="--", label=f"{name}, 1 M (ΔG°)")
    ax.axhline(0, color=INK3, lw=1)
    _style(ax, "T (K)", unit, f"Reaction free energy, ε = {res['conditions']['eps']:g}, given energies")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "dG_T.png", dpi=160)
    plt.close(fig)

    # the single-model plots: the first model with numbers (solution models need the CP2K results)
    main = next((n for n, c in ms.items() if np.isfinite(model.reaction(res, c)["dG"]).all()), "gas phase")
    i = _ti(res)
    r = model.reaction(res, ms[main])
    vals = [float(r["parts"][k][i]) for k, _ in PARTS] + [float(r["kTlnQ"][i]), float(r["dG"][i])]
    labels = [l for _, l in PARTS] + ["kT ln Q", "ΔG"]
    fig, ax = plt.subplots(figsize=(7, 4.4))
    y = np.arange(len(vals))[::-1]
    ax.barh(y, vals, color=[SERIES[0]] * (len(vals) - 1) + [INK2], height=0.6)
    for yy, v in zip(y, vals):
        ax.text(v, yy, f" {v:+.3f} " if np.isfinite(v) else " n/a ", va="center",
                ha="left" if v >= 0 else "right", fontsize=8, color=INK)
    ax.set_yticks(y, labels)
    ax.axvline(0, color=INK3, lw=1)
    fin = [v for v in vals if np.isfinite(v)]
    lo, hi = min(fin + [0.0]), max(fin + [0.0])
    pad = 0.22 * (hi - lo or 1.0)
    ax.set_xlim(lo - pad, hi + pad)                 # room for the value labels
    _style(ax, f"eV{per}", "", f"Decomposition at {res['T'][i]:g} K ({main})")
    fig.tight_layout()
    fig.savefig(out / "decomposition.png", dpi=160)
    plt.close(fig)

    x = res["scan"]["species"]
    lc = np.linspace(*res["scan"]["log10_c"], 181)
    mp = model.scan_map(res, ms[main], x, lc)
    cmap = LinearSegmentedColormap.from_list("div", [SERIES[0], "#e8e7e3", SERIES[1]])
    z = mp["dG"]
    lim = np.nanmax(np.abs(z)) if np.isfinite(z).any() else 1.0
    fig, ax = plt.subplots(figsize=(7, 4.4))
    if np.isfinite(z).any():
        pc = ax.pcolormesh(lc, T, z, cmap=cmap, norm=TwoSlopeNorm(0.0, -lim, lim), shading="auto")
        ax.contour(lc, T, z, levels=[0.0], colors=INK, linewidths=1.5)
        fig.colorbar(pc, ax=ax, label=unit)
    _style(ax, f"log10 c({x} / M)", "T (K)", f"ΔG over T and μ({x}), {main} (blue: products favoured)")
    fig.tight_layout()
    fig.savefig(out / "map_T_c.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4.0))
    for k, (name, cond) in enumerate(ms.items()):
        eq = model.equilibrium(res, cond)
        ax.plot(T, 100 * eq["conversion"], color=SERIES[k], lw=2, label=name)
    _style(ax, "T (K)", "conversion (%)", "Equilibrium conversion from the given reactant concentrations")
    ax.set_ylim(-2, 102)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "conversion.png", dpi=160)
    plt.close(fig)

    if res.get("recipe"):
        rc = res["recipe"]
        eqs = np.geomspace(1, 1e4, 41)
        sets = [k for k in ("given", "pbe0_mace", "pbe0d3_mace") if k in (res.get("energy_sets") or {"given": "given"})]
        fig, axes = plt.subplots(1, len(sets), figsize=(4.6 * len(sets), 3.8), squeeze=False, sharey=True)
        adds = [0.0, 1e-3, 1e-1, 1.0] if rc.get("additive") else [0.0]
        for ax, es in zip(axes.flat, sets):
            for k, ca in enumerate(adds):
                th = [model.speciate(res, {**ms[main], "energy": es}, equivalents=e, c_additive=ca, ti=i)["theta"][i] for e in eqs]
                ax.plot(eqs, 100 * np.array(th), color=SERIES[k], lw=2,
                        label=f"{rc['additive']} {ca:g} M" if ca else f"no {rc.get('additive') or 'additive'} added")
            ax.set_xscale("log")
            ax.set_ylim(-2, 102)
            _style(ax, f"{rc['reagent']} equivalents per dot", "dots converted (%)",
                   (res.get("energy_sets") or {}).get(es, es))
        axes.flat[0].legend(frameon=False, fontsize=8)
        fig.suptitle(f"Recipe at {res['T'][i]:g} K, dots {rc['c_dot']:g} M, {main}", fontsize=10, x=0.01, ha="left")
        fig.tight_layout()
        fig.savefig(out / "recipe.png", dpi=160)
        plt.close(fig)

    eps = np.geomspace(solvation.NONEL_MIN_EPS, 80, 120)
    names = list(res["species"])
    ncol = min(3, len(names))
    nrow = math.ceil(len(names) / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 3.2 * nrow), squeeze=False)
    for ax, n in zip(axes.flat, names):
        s = res["species"][n]
        for k, (lab, el, cv) in enumerate([("COSMO el. (RESP)", "cosmo", "none"), ("COSMO + cav/disp (β=0)", "cosmo", "beta0"),
                                           ("SCCS el.", "sccs", "none"), ("GB el. (RESP)", "gb", "none")]):
            ax.plot(eps, [solvation.dg_solv(s, e, el, cv) for e in eps], color=SERIES[k], lw=2, label=lab)
        if s.get("sccs"):
            pts = s["sccs"]["points"]
            ax.plot([float(e) for e in pts], [p["dG_el_eV"] for p in pts.values()], "o", ms=5, color=SERIES[0])
        alpb = (s.get("alpb") or {}).get("dG_eV") or {}
        xs = [solvation.ALPB_EPS[k] for k, v in alpb.items() if v is not None]
        ys = [v for v in alpb.values() if v is not None]
        if xs:
            ax.plot(xs, ys, "s", ms=5, color=SERIES[3], label="ALPB (GFN2)")
        ax.set_xscale("log")
        ax.set_xlim(1.7, 85)             # also when a species has no data (all NaN)
        _style(ax, "ε", "ΔG solv (eV)", n)
    for ax in axes.flat[len(names):]:
        ax.axis("off")
    axes.flat[0].legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(out / "solvation_eps.png", dpi=160)
    plt.close(fig)

    ti = _ti(res)
    fig, ax = plt.subplots(figsize=(7, 4.0))
    for k, (lab, el, cv) in enumerate([("COSMO", "cosmo", "none"), ("COSMO + cav/disp (β=0)", "cosmo", "beta0"),
                                       ("SCCS", "sccs", "none"), ("SCCS + cav/disp (β=0)", "sccs", "beta0")]):
        ys = [model.reaction(res, {"eps": e, "electrostatics": el, "cavity": cv})["dG"][ti] for e in eps]
        ax.plot(eps, ys, color=SERIES[k], lw=2, label=lab)
    ax.axhline(model.reaction(res, {"eps": 1.0, "gas": True})["dG"][ti], color=INK3, lw=1, ls="--", label="gas phase")
    ax.set_xscale("log")
    ax.set_xlim(1.7, 85)
    _style(ax, "ε", unit, f"ΔG at {res['T'][ti]:g} K against the dielectric constant")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "dG_eps.png", dpi=160)
    plt.close(fig)


# --------------------------------------------------------------------------
# Interactive page
# --------------------------------------------------------------------------

JS_MODEL = r"""
const KB = 8.617333262e-5, RGAS = 8.314462618, C0 = 1000.0;
const MN = 1e-3 / 1.602176634e-19 * 1e-20, GPA = 1e9 / 1.602176634e-19 * 1e-30;
const fOf = e => 1 - 1 / e;
function interpF(knots, eps) {
  const f = fOf(eps), n = knots.length;
  if (n === 1) return knots[0][1];
  let j = 0; while (j < n - 2 && knots[j + 1][0] < f) j++;
  return knots[j][1] + (knots[j + 1][1] - knots[j][1]) * (f - knots[j][0]) / (knots[j + 1][0] - knots[j][0]);
}
function cavTerm(eps, S, V, model) {   // Hille et al. 2019: beta0 = their beta = 0 fit (surface only), hille = Eq. 19
  if (eps <= 1 || model === "none") return 0;
  if (eps < NONEL.min_eps) return NaN;
  if (model === "hille") return (NONEL.alpha_gamma_mN_m[0] * Math.log(eps) + NONEL.alpha_gamma_mN_m[1]) * MN * S + NONEL.beta_GPa * GPA * V;
  return (NONEL.beta0_alpha_gamma_mN_m[0] * Math.log(eps) + NONEL.beta0_alpha_gamma_mN_m[1]) * MN * S;
}
function elec(sp, eps, model) {     // cosmo | sccs | gb
  if (eps <= 1) return 0;
  if (model === "sccs") return sp.sccs ? interpF(sp.sccs.el, eps) : NaN;
  const v = model === "cosmo" ? sp.cosmo_inf_eV : sp.gb_inf_eV;
  return (v === null || v === undefined) ? NaN : fOf(eps) * v;
}
function dgSolv(sp, c) {          // a gas species is in equilibrium with its headspace: no solvation term
  if (sp.phase === "gas" || c.gas || c.eps <= 1) return 0;
  const el = elec(sp, c.eps, c.el);
  if (c.cav === "none") return el;
  if (!sp.cavity) return NaN;
  return el + cavTerm(c.eps, sp.cavity.S_A2, sp.cavity.V_A3, c.cav);
}
function E(sp) { return st.eset === "given" ? sp.energy_eV : ((sp.energies || {})[st.eset] ?? NaN); }
function lnA(sp, t) {             // ln a: c / 1 M (dissolved), c_liquid / 1 M (pure liquid), p / (c0 R T) (gas)
  if (sp.phase === "solvent") return Math.log(sp.c_liquid);
  if (sp.phase === "gas") return Math.log(sp.p_bar * 1e5 / (C0 * RGAS * t));
  return Math.log(sp.c);
}
const PK = ["dE_DFT", "dZPE", "dH_thermal", "-TdS_trans", "-TdS_rot", "-TdS_vib", "-TdS_el", "ddG_solv_el", "ddG_solv_nonel", "ligand", "kTlnQ"];
function reaction(c) {            // species as set on the page: SP(n)
  const T = R.T, out = {dG0: [], dG: [], parts: {}};
  PK.forEach(k => out.parts[k] = []);
  const sp = {}, solvAll = {}, solvEl = {};
  for (const n in R.nu) { sp[n] = SP(n); solvAll[n] = dgSolv(sp[n], c); solvEl[n] = dgSolv(sp[n], {...c, cav: "none"}); }
  const lig = R.ligand ? R.nu[R.ligand.species] * st.dglig : 0;     // nu_L dG_lig: released ligands (nu_L > 0)
  T.forEach((t, i) => {
    let g0 = lig, lq = 0; const p = Object.fromEntries(PK.map(k => [k, 0]));
    for (const n in R.nu) {
      const v = R.nu[n], s = sp[n], th = s.thermo;
      g0 += v * (E(s) + th.G_corr[i] + solvAll[n]);
      lq += v * lnA(s, t);
      p["dE_DFT"] += v * E(s); p["dZPE"] += v * th.zpe[i];
      p["dH_thermal"] += v * (th.dU_vib[i] + th.H_trans[i] + th.H_rot[i]);
      p["-TdS_trans"] -= v * t * th.S_trans[i]; p["-TdS_rot"] -= v * t * th.S_rot[i];
      p["-TdS_vib"] -= v * t * th.S_vib[i]; p["-TdS_el"] -= v * t * th.S_el[i];
      p["ddG_solv_el"] += v * solvEl[n]; p["ddG_solv_nonel"] += v * (solvAll[n] - solvEl[n]);
    }
    p["ligand"] = lig; p["kTlnQ"] = KB * t * lq;
    out.dG0.push(g0 / PER); out.dG.push((g0 + KB * t * lq) / PER);
    PK.forEach(k => out.parts[k].push(p[k] / PER));
  });
  return out;
}
// Recipe: the reaction and its side equilibria (e.g. amine + HCl -> adduct) from the amounts put in.
//   c_i = c_i0 + sum_r nu_ri xi_r (dissolved), a_i fixed (gas, pure liquid); g_r = dG_r0/kT + sum_i nu_ri ln a_i = 0.
// Main extent by bisection; for each trial the side extents: closed form for one A + B -> C, else bisection sweeps.
function recipeSetup(c) {
  const rx = [R.nu].concat((R.equilibria || []).map(e => e.nu));
  const names = [...new Set(rx.flatMap(nu => Object.keys(nu)))];
  const sp = Object.fromEntries(names.map(n => [n, SP(n)]));
  const bal = names.filter(n => sp[n].phase === "solute");
  const G = Object.fromEntries(names.map(n => [n, dgSolv(sp[n], c)]));
  const lig = R.ligand ? R.nu[R.ligand.species] * st.dglig : 0;
  const dg0 = rx.map((nu, r) => R.T.map((t, i) => { let g = r === 0 ? lig : st.dgbind;
    for (const n in nu) g += nu[n] * (E(sp[n]) + sp[n].thermo.G_corr[i] + G[n]); return g; }));
  const fixed = rx.map(nu => R.T.map(t => { let f = 0; for (const n in nu) if (!bal.includes(n)) f += nu[n] * lnA(sp[n], t); return f; }));
  const s1 = rx[1] || {}, k1 = Object.keys(s1);
  const quad = rx.length === 2 && k1.length === 3 && k1.every(n => bal.includes(n) && Math.abs(s1[n]) === 1) && k1.filter(n => s1[n] < 0).length === 2;
  return {rx, bal, dg0, fixed, quad};
}
function speciate(S, i, cdot, eqv, cadd) {
  const kT = KB * R.T[i], rx = S.rx, rc = R.recipe;
  const n0 = Object.fromEntries(S.bal.map(n => [n, 0]));
  n0[rc.dot] = cdot; n0[rc.reagent] += eqv * cdot; if (rc.additive in n0) n0[rc.additive] += cadd;
  const shift = (a0, nu, x) => { const a = {...a0}; for (const n in nu) if (n in a) a[n] += nu[n] * x; return a; };
  const amounts = xi => rx.reduce((a, nu, r) => xi[r] ? shift(a, nu, xi[r]) : a, n0);
  const g = (r, a) => { let f = S.dg0[r][i] / kT + S.fixed[r][i]; for (const n in rx[r]) if (n in a) f += rx[r][n] * Math.log(Math.max(a[n], 1e-300)); return f; };
  const bounds = (nu, a) => { let lo = -Infinity, hi = Infinity;
    for (const n in nu) { if (!(n in a)) continue; const v = nu[n]; if (v > 0) lo = Math.max(lo, -a[n] / v); else hi = Math.min(hi, a[n] / -v); } return [lo, hi]; };
  const bis = (f, lo, hi) => { for (let q = 0; q < 64; q++) { const m = 0.5 * (lo + hi); if (f(m) < 0) lo = m; else hi = m; } return 0.5 * (lo + hi); };
  const sides = x0 => {
    const xi = [x0].concat(rx.slice(1).map(() => 0));
    if (rx.length === 1) return xi;
    if (S.quad) {             // (c + y) = K (a - y)(b - y): the root in [-c, min(a, b)]
      const nu = rx[1], base = amounts([x0]), K = Math.exp(-(S.dg0[1][i] / kT + S.fixed[1][i]));
      const [A, B] = Object.keys(nu).filter(n => nu[n] < 0), P = Object.keys(nu).find(n => nu[n] > 0);
      const a = Math.max(base[A], 0), b = Math.max(base[B], 0), cc = Math.max(base[P], 0);
      xi[1] = 2 * (K * a * b - cc) / (K * (a + b) + 1 + Math.sqrt(K * K * (a - b) ** 2 + 2 * K * (a + b) + 1 + 4 * K * cc));
      return xi;
    }
    for (let sweep = 0; sweep < (rx.length > 2 ? 30 : 1); sweep++)
      for (let r = 1; r < rx.length; r++) {
        const rest = xi.slice(); rest[r] = 0; const a0 = amounts(rest), [lo, hi] = bounds(rx[r], a0);
        if (isFinite(lo) && isFinite(hi) && hi > lo) xi[r] = bis(x => g(r, shift(a0, rx[r], x)), lo, hi);
      }
    return xi;
  };
  const [lo, hi] = bounds(rx[0], n0);
  if (!isFinite(lo) || !isFinite(hi) || hi <= lo || !isFinite(S.dg0[0][i])) return {theta: NaN, a: {}};
  const x0 = bis(x => g(0, amounts(sides(x))), lo, hi);
  return {theta: x0 / cdot, a: amounts(sides(x0))};
}
function equilibrium(c) {         // mass balance over the dissolved species; gas and pure liquid at fixed activity
  const r = reaction(c), names = Object.keys(R.nu), sp = Object.fromEntries(names.map(n => [n, SP(n)]));
  const sol = names.filter(n => sp[n].phase === "solute"), oth = names.filter(n => !sol.includes(n));
  const init = {}; sol.forEach(n => init[n] = R.nu[n] < 0 ? sp[n].c : (R.species[n].c0 ?? 0));
  let lo = 0, hi = Infinity;
  sol.forEach(n => { const v = R.nu[n]; if (v > 0) lo = Math.max(lo, -init[n] / v); else hi = Math.min(hi, init[n] / -v); });
  return R.T.map((t, i) => {
    if (!isFinite(hi) || hi <= lo || !isFinite(r.dG0[i])) return NaN;
    const kT = KB * t;
    let target = -r.dG0[i] * PER; oth.forEach(n => target -= R.nu[n] * kT * lnA(sp[n], t));
    target /= kT;
    let a = lo, b = hi;
    for (let k = 0; k < 200; k++) {
      const m = 0.5 * (a + b); let f = -target;
      sol.forEach(n => f += R.nu[n] * Math.log(Math.max(init[n] + R.nu[n] * m, 1e-300)));
      if (f < 0) a = m; else b = m;
    }
    return 100 * (0.5 * (a + b) - lo) / (hi - lo);
  });
}
"""

PAGE = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__</title>
<script src="https://cdn.jsdelivr.net/npm/plotly.js-dist-min@2.35.2/plotly.min.js"></script>
<style>
:root { --bg: #fcfcfb; --ink: #0b0b0b; --ink2: #52514e; --ink3: #8a8984; --rule: #e6e5e0; --panel: #f4f4f1; --flag: #8a5a00; }
body { margin: 0; background: var(--bg); color: var(--ink); font: 14px/1.5 system-ui, -apple-system, "Segoe UI", sans-serif; }
.wrap { max-width: 1180px; margin: 0 auto; padding-inline: 16px; padding-block: 20px 40px; }
h1 { font-size: 1.35rem; margin: 0 0 .2rem; } h2 { font-size: 1.1rem; margin: 28px 0 .2rem; } .eq { color: var(--ink2); font-family: ui-monospace, Menlo, monospace; font-size: .9rem; overflow-wrap: anywhere; }
.controls { display: grid; grid-template-columns: repeat(auto-fit, minmax(230px, 1fr)); gap: 10px 22px; background: var(--panel); padding: 12px 14px; border-radius: 8px; margin: 14px 0; }
.controls label { display: grid; gap: 2px; font-size: .82rem; color: var(--ink2); }
.controls .v { color: var(--ink); font-variant-numeric: tabular-nums; font-weight: 600; }
input[type=range] { width: 100%; }
.grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(min(100%, 520px), 1fr)); gap: 16px; }
.card { min-width: 0; } .plot { width: 100%; height: 340px; }
.kpi { display: flex; flex-wrap: wrap; gap: 8px 26px; margin: 8px 0 4px; font-variant-numeric: tabular-nums; }
.kpi div span { display: block; font-size: .75rem; color: var(--ink2); text-transform: uppercase; letter-spacing: .05em; }
.kpi div b { font-size: 1.15rem; }
.flags { color: var(--flag); font-size: .82rem; } table { border-collapse: collapse; font-size: .82rem; font-variant-numeric: tabular-nums; }
.tbl { overflow-x: auto; } td, th { padding: 4px 8px; border-bottom: 1px solid var(--rule); text-align: right; } th:first-child, td:first-child { text-align: left; }
.note { color: var(--ink2); font-size: .8rem; max-width: 80ch; }
</style></head><body><div class="wrap">
<h1>__TITLE__</h1><div class="eq">__EQ__</div>
<div class="kpi" id="kpi"></div>
<div class="note" id="alt"></div>
<div class="controls" id="ctl"></div>
<div class="grid">
 <div class="card"><div class="plot" id="pT"></div></div>
 <div class="card"><div class="plot" id="pDec"></div></div>
 <div class="card"><div class="plot" id="pMap"></div></div>
 <div class="card"><div class="plot" id="pEps"></div></div>
 <div class="card"><div class="plot" id="pConv"></div></div>
 <div class="card"><div class="tbl" id="tab"></div></div>
</div>
<div id="recipe" hidden>
<h2>Recipe: what is put in the flask</h2>
<div class="note" id="rcNote"></div>
<div class="controls" id="rctl"></div>
<div class="kpi" id="rkpi"></div>
<div class="grid">
 <div class="card"><div class="plot" id="pRecEq"></div></div>
 <div class="card"><div class="plot" id="pRecMap"></div></div>
 <div class="card"><div class="plot" id="pRecSpec"></div></div>
 <div class="card"><div class="plot" id="pRecT"></div></div>
</div>
</div>
<div class="flags" id="flags"></div>
<p class="note">ΔG = Σ ν<sub>i</sub> μ<sub>i</sub> + ν<sub>L</sub> ΔG<sub>lig</sub>. Dissolved: μ = E<sup>DFT</sup> + G<sup>corr</sup>(T) + ΔG<sub>solv</sub>(ε) + kT ln(c/1 M);
gas in equilibrium with the headspace: μ = E<sup>DFT</sup> + G<sup>corr</sup>(T) + kT ln(p/c°RT), no solvation term; pure liquid: as dissolved, at c = c<sub>liquid</sub>.
ΔG<sub>lig</sub>: free energy per released ligand that the model leaves out (chain packing in the shell, conformational entropy).
G<sup>corr</sup>: ideal gas at 1 M, MACE-MH-1 Hessian at the MACE minimum, quasi-RRHO entropy.
Electrostatics: COSMO (conductor, scaled by 1 − 1/ε) on the DFT electrostatic potential of the gas-phase single point, with the cavity on the electron-density isosurface ρ₀ = e<sup>−7.6</sup> = 5×10<sup>−4</sup> e/bohr³ (the centre of the SCCS switching); or CP2K SCCS interpolated in 1 − 1/ε (molecules); or Generalized Born on RESP charges.
Cavity + dispersion (Hille et al., J. Chem. Phys. 150, 041710, 2019): β = 0 fit, (α+γ) = 2.16 ln ε − 14.79 mN/m fitted here on their Table S2, times the area S of the same isosurface (within 1–4 % of the SCCS quantum surface); or their main fit (3.0 ln ε + 36.0) S − 0.5 GPa V, molecules only. Defined for ε ≥ 1.8.</p>
</div>
<script>
const R = __DATA__;
const NONEL = R.nonel_model;
const PER = __PER__;
__JS__
const C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"];
const UNIT = "eV" + (R.per ? " per " + R.per : "");
const PHASES = {solute: ["dissolved", "c", "M", -12, 1.5], gas: ["gas (headspace)", "p", "bar", -10, 1.5], solvent: ["liquid (neat or diluted)", "c", "M", -3, 1.5]};
const st = {eps: R.conditions.eps, gas: false, el: R.conditions.electrostatics || "cosmo",
            cav: R.conditions.cavity || "beta0", eset: "given",
            ti: R.T.reduce((b, t, i) => Math.abs(t - 298.15) < Math.abs(R.T[b] - 298.15) ? i : b, 0),
            x: R.scan.species, dglig: R.ligand ? (R.ligand.dG_eV || 0) : 0, dgbind: 0,
            rc: R.recipe ? {lgdot: Math.log10(R.recipe.c_dot), lgeq: Math.log10(R.recipe.equivalents),
                            add: R.recipe.c_additive > 0, lgadd: R.recipe.c_additive > 0 ? Math.log10(R.recipe.c_additive) : -2} : null,
            phase: Object.fromEntries(Object.entries(R.species).map(([n, s]) => [n, s.phase])),
            // each phase keeps its own value (log10 of c in M, p in bar, or c of the pure liquid in M)
            lg: Object.fromEntries(Object.entries(R.species).map(([n, s]) => [n, {solute: Math.log10(s.c), gas: Math.log10(s.p_bar || 1),
                                                                                  solvent: Math.log10(s.c_liquid || 1)}]))};
function SP(n) {                  // species n with the phase and amount chosen on the page
  const s = R.species[n], ph = st.phase[n], v = Math.pow(10, st.lg[n][ph]);
  return Object.assign({}, s, {phase: ph, c: ph === "solute" ? v : s.c, p_bar: ph === "gas" ? v : s.p_bar,
                               c_liquid: ph === "solvent" ? v : s.c_liquid});
}
const LAY = (xt, yt, extra = {}) => Object.assign({template: "plotly_white", margin: {l: 64, r: 16, t: 36, b: 48},
  xaxis: {title: {text: xt}, gridcolor: "#e6e5e0", zeroline: false}, yaxis: {title: {text: yt}, gridcolor: "#e6e5e0", zeroline: false},
  font: {family: "system-ui, sans-serif", size: 12, color: "#0b0b0b"}, legend: {orientation: "h", y: -0.22}, hovermode: "x unified"}, extra);
const fmt = v => isFinite(v) ? (v >= 0 ? "+" : "") + v.toFixed(3) : "n/a";
function ctl() {
  const c = document.getElementById("ctl");
  let h = `<label><span>Dielectric constant ε <span class="v" id="vEps"></span></span><input id="eps" type="range" min="${Math.log10(1.8)}" max="${Math.log10(80)}" step="0.005" value="${Math.log10(st.eps)}"></label>`;
  h += `<label><span>Medium</span><select id="gas"><option value="0">solution</option><option value="1">gas phase</option></select></label>`;
  h += `<label><span>Electronic energy</span><select id="eset">${Object.entries(R.energy_sets || {given: "given"}).map(([k, l]) => `<option value="${k}">${l}</option>`).join("")}</select></label>`;
  h += `<label><span>Electrostatics</span><select id="el"><option value="cosmo">COSMO on the DFT potential, density isosurface (any size)</option><option value="sccs">CP2K SCCS (electronic polarisation; molecules)</option><option value="gb">Generalized Born on RESP charges</option></select></label>`;
  h += `<label><span>Cavity + dispersion (Hille 2019)</span><select id="cav"><option value="beta0">β = 0 fit: (α+γ)(ε) S</option><option value="hille">β = −0.5 GPa fit (molecules only)</option><option value="none">left out</option></select></label>`;
  h += `<label><span>Temperature for the decomposition <span class="v" id="vT"></span></span><input id="ti" type="range" min="0" max="${R.T.length - 1}" step="1" value="${st.ti}"></label>`;
  h += `<label><span>Chemical-potential axis</span><select id="x"></select></label>`;
  if (R.ligand) h += `<label title="Free energy per released ${R.ligand.species} beyond the model: chain packing in the ligand shell, conformational entropy (positive: binding stronger than the model)"><span>ΔG<sub>lig</sub> per released ${R.ligand.species} (eV) <span class="v" id="vLig"></span></span><input id="dglig" type="range" min="-0.5" max="0.5" step="0.01" value="${st.dglig}"></label>`;
  for (const n in R.nu)
    h += `<label><span>${n} <select id="ph_${n}">${Object.entries(PHASES).map(([k, d]) => `<option value="${k}">${d[0]}</option>`).join("")}</select> log<sub>10</sub> <span id="q_${n}"></span> <span class="v" id="v_${n}"></span></span><input id="lg_${n}" type="range" step="0.05"></label>`;
  c.innerHTML = h;
  document.getElementById("el").value = st.el; document.getElementById("cav").value = st.cav; document.getElementById("eset").value = st.eset;
  for (const n in R.nu) { document.getElementById("ph_" + n).value = st.phase[n]; slider(n); }
  axisOptions();
  c.addEventListener("input", e => {
    const id = e.target.id, v = e.target.value;
    if (id === "eps") st.eps = Math.pow(10, +v); else if (id === "gas") st.gas = v === "1"; else if (id === "el") st.el = v;
    else if (id === "cav") st.cav = v; else if (id === "eset") st.eset = v; else if (id === "ti") st.ti = +v; else if (id === "x") st.x = v;
    else if (id === "dglig") st.dglig = +v;
    else if (id.startsWith("ph_")) { const n = id.slice(3); st.phase[n] = v; slider(n); axisOptions(); }
    else if (id.startsWith("lg_")) { const n = id.slice(3); st.lg[n][st.phase[n]] = +v; }
    draw();
  });
}
function slider(n) {              // the slider of species n for its current phase
  const d = PHASES[st.phase[n]], el = document.getElementById("lg_" + n);
  el.min = d[3]; el.max = d[4]; el.value = st.lg[n][st.phase[n]];
  document.getElementById("q_" + n).textContent = `${d[1]} / ${d[2]}`;
}
function axisOptions() {          // every species can be the map's axis (a liquid such as the amine: diluted in the solvent)
  const ok = Object.keys(R.nu);
  if (!ok.includes(st.x)) st.x = ok[0];
  document.getElementById("x").innerHTML = ok.map(n => `<option ${n === st.x ? "selected" : ""}>${n}</option>`).join("");
}
function cond() { return {eps: st.gas ? 1 : st.eps, gas: st.gas, el: st.el, cav: st.cav}; }
function draw() {
  document.getElementById("vEps").textContent = st.gas ? "(gas)" : st.eps.toFixed(2);
  document.getElementById("vT").textContent = R.T[st.ti] + " K";
  for (const n in R.nu) document.getElementById("v_" + n).textContent = st.lg[n][st.phase[n]].toFixed(2);
  if (R.ligand) document.getElementById("vLig").textContent = (st.dglig >= 0 ? "+" : "") + st.dglig.toFixed(2);
  const c = cond(), r = reaction(c), i = st.ti, t = R.T[i];
  let cross = null; for (let k = 0; k < R.T.length - 1; k++) if (Math.sign(r.dG[k]) !== Math.sign(r.dG[k + 1])) { cross = R.T[k] + (0 - r.dG[k]) * (R.T[k + 1] - R.T[k]) / (r.dG[k + 1] - r.dG[k]); break; }
  const alt = Object.entries(R.dE_alternatives_eV || {}).map(([l, v]) => `${l}: ${v === null ? "n/a" : fmt(v / PER)}`).join(" · ");
  document.getElementById("alt").textContent = alt ? "ΔE (" + UNIT + ") " + alt + "  (the reaction uses the given energies)" : "";
  document.getElementById("kpi").innerHTML = `<div><span>ΔG at ${t} K</span><b>${fmt(r.dG[i])} ${UNIT}</b></div><div><span>ΔG° (1 M)</span><b>${fmt(r.dG0[i])}</b></div>` +
    `<div><span>kcal/mol</span><b>${fmt(r.dG[i] * 23.060548)}</b></div><div><span>ΔG = 0 at</span><b>${cross ? cross.toFixed(0) + " K" : "none in range"}</b></div>`;
  Plotly.react("pT", [
    {x: R.T, y: r.dG, name: "ΔG at the chosen c", line: {color: C[0], width: 2}, hovertemplate: "%{y:.3f} " + UNIT},
    {x: R.T, y: r.dG0, name: "ΔG° (1 M)", line: {color: C[0], width: 1.4, dash: "dash"}, hovertemplate: "%{y:.3f}"},
    {x: R.T, y: reaction({eps: 1, gas: true}).dG, name: "gas phase", line: {color: "#8a8984", width: 1.4, dash: "dot"}, hovertemplate: "%{y:.3f}"}],
    LAY("T (K)", "ΔG (" + UNIT + ")", {title: {text: "Reaction free energy", x: 0, font: {size: 13}}, shapes: [{type: "line", xref: "paper", x0: 0, x1: 1, y0: 0, y1: 0, line: {color: "#8a8984", width: 1}}]}), {displaylogo: false, responsive: true});
  const PL = [["dE_DFT", "ΔE DFT"], ["dZPE", "ΔZPE"], ["dH_thermal", "ΔH thermal"], ["-TdS_trans", "−TΔS trans"], ["-TdS_rot", "−TΔS rot"], ["-TdS_vib", "−TΔS vib"], ["-TdS_el", "−TΔS el"], ["ddG_solv_el", "ΔΔG solv el."], ["ddG_solv_nonel", "ΔΔG cav+disp"], ["ligand", "ΔG lig (chain)"], ["kTlnQ", "kT ln Q"]];
  const dv = PL.map(([k]) => r.parts[k][i]).concat([r.dG[i]]), dl = PL.map(p => p[1]).concat(["ΔG"]);
  Plotly.react("pDec", [{type: "bar", orientation: "h", y: dl.slice().reverse(), x: dv.slice().reverse(), marker: {color: dl.map((_, k) => k === 0 ? "#52514e" : C[0])},
    text: dv.slice().reverse().map(fmt), textposition: "outside", cliponaxis: false, hovertemplate: "%{y}: %{x:.4f} " + UNIT + "<extra></extra>"}],
    LAY(UNIT, "", {title: {text: `Decomposition at ${t} K`, x: 0, font: {size: 13}}, margin: {l: 120, r: 60, t: 36, b: 48}, hovermode: "closest", showlegend: false}), {displaylogo: false, responsive: true});
  const dx = PHASES[st.phase[st.x]], nx = R.nu[st.x], sx = SP(st.x);
  const lc = Array.from({length: 121}, (_, k) => dx[3] + (dx[4] - dx[3]) * k / 120);
  const z = R.T.map((tt, k) => { const base = r.dG[k] - nx * KB * tt * lnA(sx, tt) / PER;
    return lc.map(l => base + nx * KB * tt * lnA({...sx, c: Math.pow(10, l), p_bar: Math.pow(10, l), c_liquid: Math.pow(10, l)}, tt) / PER); });
  const zmax = Math.max(...z.flat().filter(isFinite).map(Math.abs), 1e-6);
  const qx = `log10 ${dx[1]}(${st.x} / ${dx[2]})`;
  Plotly.react("pMap", [{type: "heatmap", x: lc, y: R.T, z: z, zmin: -zmax, zmax: zmax, colorscale: [[0, C[0]], [0.5, "#e8e7e3"], [1, C[1]]],
    colorbar: {title: {text: UNIT, side: "right"}}, hovertemplate: qx + " %{x:.2f}<br>T %{y} K<br>ΔG %{z:.3f}<extra></extra>"},
    {type: "contour", x: lc, y: R.T, z: z, contours: {start: 0, end: 0, size: 1, coloring: "none", showlabels: false}, line: {color: "#0b0b0b", width: 2}, showscale: false, hoverinfo: "skip"}],
    LAY(qx, "T (K)", {title: {text: `ΔG over T and μ(${st.x}), ${dx[0]}; line: ΔG = 0, blue: products favoured`, x: 0, font: {size: 13}}, hovermode: "closest"}), {displaylogo: false, responsive: true});
  const eg = Array.from({length: 100}, (_, k) => Math.pow(10, Math.log10(1.8) + (Math.log10(80) - Math.log10(1.8)) * k / 99));
  const ser = [["cosmo", "beta0", "COSMO + cav/disp (β=0)"], ["cosmo", "none", "COSMO"], ["sccs", "beta0", "SCCS + cav/disp (β=0)"], ["sccs", "none", "SCCS"]];
  const tr = ser.map(([el, cv, nm], k) => ({x: eg, y: eg.map(e => reaction({eps: e, gas: false, el: el, cav: cv}).dG[i]), name: nm, line: {color: C[k], width: 2}, hovertemplate: "%{y:.3f}"}));
  tr.push({x: [st.eps], y: [r.dG[i]], mode: "markers", marker: {size: 10, color: "#0b0b0b"}, name: "current", hoverinfo: "skip"});
  Plotly.react("pEps", tr, LAY("ε", "ΔG (" + UNIT + ")", {title: {text: `ΔG at ${t} K against ε`, x: 0, font: {size: 13}}, xaxis: {type: "log", title: {text: "ε"}, gridcolor: "#e6e5e0"}}), {displaylogo: false, responsive: true});
  Plotly.react("pConv", [{x: R.T, y: equilibrium(c), name: "conversion", line: {color: C[2], width: 2}, hovertemplate: "%{y:.1f} %"}],
    LAY("T (K)", "conversion (%)", {title: {text: "Equilibrium conversion (dissolved reactants at their c)", x: 0, font: {size: 13}}, yaxis: {range: [-2, 102], title: {text: "conversion (%)"}, gridcolor: "#e6e5e0"}, showlegend: false}), {displaylogo: false, responsive: true});
  let tb = `<table><tr><th>species</th><th>ν</th><th>atoms</th><th>E (${st.eset}) eV</th><th>G corr ${t} K</th><th>ΔG solv</th><th>COSMO el.</th><th>SCCS el.</th><th>GB el.</th><th>cav+disp</th><th>dipole DFT / RESP / Hirshfeld (D)</th><th>RMSD Å</th></tr>`;
  for (const n in R.species) { const s = SP(n), m = s.mace || {}, d = s.dipole_D || {};
    const co = elec(s, c.eps, "cosmo"), sc = elec(s, c.eps, "sccs"), gb = elec(s, c.eps, "gb"), all = dgSolv(s, c);
    const cv = s.cavity ? cavTerm(c.gas ? 1 : c.eps, s.cavity.S_A2, s.cavity.V_A3, st.cav) : NaN;
    const dd = [d.dft, d.resp, d.hirshfeld].map(v => v === null || v === undefined ? "–" : v.toFixed(2)).join(" / ");
    tb += `<tr><td>${n}</td><td>${R.nu[n] ?? "equil."}</td><td>${s.n_atoms}</td><td>${E(s).toFixed(4)}</td><td>${fmt(s.thermo.G_corr[i])}</td><td>${fmt(all)}</td><td>${fmt(co)}</td><td>${fmt(sc)}</td><td>${fmt(gb)}</td><td>${fmt(cv)}</td><td>${dd}</td><td>${m.rmsd_heavy_A !== undefined ? m.rmsd_heavy_A.toFixed(3) : "n/a"}</td></tr>`; }
  document.getElementById("tab").innerHTML = tb + "</table>";
  if (R.recipe) drawRecipe(c);
  const fl = Object.entries(R.species).flatMap(([n, s]) => s.flags.map(f => `${n}: ${f}`));
  document.getElementById("flags").innerHTML = fl.length ? "<b>Flags</b><br>" + fl.map(f => f.replace(/</g, "&lt;")).join("<br>") : "";
}
function rctl() {                 // the recipe controls (dots, reagent equivalents, additive, side-equilibrium correction)
  const rc = R.recipe, side = (R.equilibria || []).map(e => e.equation).join("; ");
  document.getElementById("recipe").hidden = false;
  let h = `<label><span>Dots (${rc.dot}): log<sub>10</sub> c / M <span class="v" id="vDot"></span></span><input id="r_lgdot" type="range" min="-8" max="-3" step="0.05" value="${st.rc.lgdot}"></label>`;
  h += `<label><span>${rc.reagent} per dot: log<sub>10</sub> equivalents <span class="v" id="vEq"></span></span><input id="r_lgeq" type="range" min="0" max="4" step="0.02" value="${st.rc.lgeq}"></label>`;
  if (rc.additive) {
    h += `<label><span>${rc.additive} added</span><select id="r_add"><option value="0">none</option><option value="1">added</option></select></label>`;
    h += `<label><span>${rc.additive} added: log<sub>10</sub> c / M <span class="v" id="vAdd"></span></span><input id="r_lgadd" type="range" min="-6" max="0.5" step="0.05" value="${st.rc.lgadd}"></label>`;
  }
  if (side) h += `<label title="Added to ΔG° of ${side}: e.g. ion pairing (ammonium chloride) beyond the H-bonded complex computed here"><span>ΔG° correction of ${side} (eV) <span class="v" id="vBind"></span></span><input id="r_dgbind" type="range" min="-1" max="0.5" step="0.01" value="${st.dgbind}"></label>`;
  const c = document.getElementById("rctl"); c.innerHTML = h;
  if (rc.additive) document.getElementById("r_add").value = st.rc.add ? "1" : "0";
  c.addEventListener("input", e => {
    const id = e.target.id, v = e.target.value;
    if (id === "r_lgdot") st.rc.lgdot = +v; else if (id === "r_lgeq") st.rc.lgeq = +v; else if (id === "r_add") st.rc.add = v === "1";
    else if (id === "r_lgadd") st.rc.lgadd = +v; else if (id === "r_dgbind") st.dgbind = +v;
    draw();
  });
  document.getElementById("rcNote").innerHTML = `Equilibrium in solution of ${R.equation}` + (side ? ` together with ${side}` : "") +
    `, from c(${rc.dot}), ${rc.reagent} = equivalents × c(${rc.dot})` + (rc.additive ? ` and the ${rc.additive} added (0 allowed: the dots release their own)` : "") +
    `. Species set to <i>gas</i> above stay at their partial pressure (open flask), <i>pure liquid</i> at fixed activity; every dissolved species follows the mass balance (the concentration sliders above do not apply here). ` +
    `θ: fraction of the dots converted (all ${-R.nu[rc.reagent]} ${rc.reagent} per dot in one step). Temperature: the decomposition slider; energy set, ε and solvation as chosen above.`;
}
const sci = v => isFinite(v) ? v.toExponential(2) : "n/a";
function drawRecipe(c) {
  const rc = R.recipe, S = recipeSetup(c), i = st.ti, t = R.T[i];
  const cdot = Math.pow(10, st.rc.lgdot), eqv = Math.pow(10, st.rc.lgeq), cadd = rc.additive && st.rc.add ? Math.pow(10, st.rc.lgadd) : 0;
  document.getElementById("vDot").textContent = sci(cdot) + " M";
  document.getElementById("vEq").textContent = eqv.toFixed(eqv < 10 ? 1 : 0) + ` (${sci(eqv * cdot)} M)`;
  if (rc.additive) { document.getElementById("vAdd").textContent = st.rc.add ? sci(cadd) + " M" : "(none)"; document.getElementById("r_lgadd").disabled = !st.rc.add; }
  if (document.getElementById("vBind")) document.getElementById("vBind").textContent = (st.dgbind >= 0 ? "+" : "") + st.dgbind.toFixed(2);
  const now = speciate(S, i, cdot, eqv, cadd), a = now.a;
  const sideP = [...new Set((R.equilibria || []).flatMap(e => Object.keys(e.nu).filter(n => e.nu[n] > 0)))];
  const shown = [rc.reagent, ...sideP, ...(rc.additive ? [rc.additive] : [])].filter(n => n in a);
  document.getElementById("rkpi").innerHTML = `<div><span>dots converted at ${t} K</span><b>${isFinite(now.theta) ? (100 * now.theta).toFixed(1) + " %" : "n/a"}</b></div>` +
    shown.map(n => `<div><span>${n === rc.reagent ? "free " : ""}${n}</span><b>${sci(a[n])} M</b></div>`).join("");
  const cfg = {displaylogo: false, responsive: true};
  const lgE = Array.from({length: 81}, (_, k) => 4 * k / 80);
  const curves = [[0, `no ${rc.additive || "additive"} added`, "dash"]].concat(cadd > 0 ? [[cadd, `${rc.additive} ${sci(cadd)} M`, "solid"]] : []);
  Plotly.react("pRecEq", curves.map(([ca, nm, ds], k) => ({x: lgE.map(l => Math.pow(10, l)), y: lgE.map(l => 100 * speciate(S, i, cdot, Math.pow(10, l), ca).theta),
      name: nm, line: {color: C[k], width: 2, dash: ds}, hovertemplate: "%{y:.1f} %"}))
    .concat([{x: [eqv], y: [100 * now.theta], mode: "markers", marker: {size: 10, color: "#0b0b0b"}, name: "recipe", hoverinfo: "skip"}]),
    LAY(`${rc.reagent} equivalents per dot`, "dots converted (%)", {title: {text: `Conversion against ${rc.reagent} at ${t} K, dots ${sci(cdot)} M`, x: 0, font: {size: 13}},
      xaxis: {type: "log", title: {text: `${rc.reagent} equivalents per dot`}, gridcolor: "#e6e5e0"}, yaxis: {range: [-2, 102], title: {text: "dots converted (%)"}, gridcolor: "#e6e5e0"}}), cfg);
  if (rc.additive) {
    const mx = Array.from({length: 41}, (_, k) => 4 * k / 40), my = Array.from({length: 27}, (_, k) => -6 + 6.5 * k / 26);
    const z = my.map(ly => mx.map(lx => 100 * speciate(S, i, cdot, Math.pow(10, lx), Math.pow(10, ly)).theta));
    Plotly.react("pRecMap", [{type: "heatmap", x: mx, y: my, z: z, zmin: 0, zmax: 100, colorscale: [[0, "#f1efe8"], [1, C[0]]],
        colorbar: {title: {text: "%", side: "right"}}, hovertemplate: `log10 eq %{x:.2f}<br>log10 ${rc.additive} %{y:.2f}<br>%{z:.1f} %<extra></extra>`},
      {type: "contour", x: mx, y: my, z: z, contours: {start: 50, end: 50, size: 1, coloring: "none"}, line: {color: "#0b0b0b", width: 2}, showscale: false, hoverinfo: "skip"}]
      .concat(cadd > 0 ? [{x: [st.rc.lgeq], y: [st.rc.lgadd], mode: "markers", marker: {size: 10, color: "#0b0b0b", line: {color: "#fff", width: 2}}, hoverinfo: "skip", showlegend: false}] : []),
      LAY(`log10 ${rc.reagent} equivalents per dot`, `log10 ${rc.additive} added (M)`, {title: {text: `Dots converted at ${t} K (line: 50 %); no ${rc.additive}: left plot`, x: 0, font: {size: 13}}, hovermode: "closest"}), cfg);
  }
  const sp = lgE.map(l => speciate(S, i, cdot, Math.pow(10, l), cadd).a);
  Plotly.react("pRecSpec", shown.map((n, k) => ({x: lgE.map(l => Math.pow(10, l)), y: sp.map(q => q[n]), name: n === rc.reagent ? `free ${n}` : n,
      line: {color: C[k + 2], width: 2}, hovertemplate: "%{y:.2e} M"})),
    LAY(`${rc.reagent} equivalents per dot`, "c (M)", {title: {text: `In solution at equilibrium, ${t} K` + (cadd > 0 ? `, ${rc.additive} ${sci(cadd)} M added` : `, no ${rc.additive || "additive"} added`), x: 0, font: {size: 13}},
      xaxis: {type: "log", title: {text: `${rc.reagent} equivalents per dot`}, gridcolor: "#e6e5e0"}, yaxis: {type: "log", exponentformat: "power", title: {text: "c (M)"}, gridcolor: "#e6e5e0"}}), cfg);
  Plotly.react("pRecT", curves.map(([ca, nm, ds], k) => ({x: R.T, y: R.T.map((_, j) => 100 * speciate(S, j, cdot, eqv, ca).theta), name: nm,
      line: {color: C[k], width: 2, dash: ds}, hovertemplate: "%{y:.1f} %"})),
    LAY("T (K)", "dots converted (%)", {title: {text: `Conversion against T at ${eqv.toFixed(eqv < 10 ? 1 : 0)} ${rc.reagent} per dot`, x: 0, font: {size: 13}},
      yaxis: {range: [-2, 102], title: {text: "dots converted (%)"}, gridcolor: "#e6e5e0"}}), cfg);
}
// no solvation data yet (CP2K still running): open on the gas phase
if (!reaction({eps: st.eps, gas: false, el: st.el, cav: st.cav}).dG.some(isFinite)) st.gas = true;
ctl(); if (R.recipe) rctl(); document.getElementById("gas").value = st.gas ? "1" : "0"; draw();
</script></body></html>
"""


def page(res: dict) -> str:
    slim = json.loads(json.dumps(res, default=float))
    for s in slim["species"].values():
        s.pop("frequencies_cm1", None)
    # NaN stays NaN: the data is a JS literal, where NaN is valid (null would count as 0 in the sums)
    data = json.dumps(slim, separators=(",", ":"), allow_nan=True).replace("</", "<\\/")
    title = html.escape(f"Reaction free energy: {res['name']}")
    return (PAGE.replace("__TITLE__", title).replace("__EQ__", html.escape(" + ".join(
        f"{-v} {n}" for n, v in res["nu"].items() if v < 0) + " → " + " + ".join(f"{v} {n}" for n, v in res["nu"].items() if v > 0)))
        .replace("__DATA__", data).replace("__PER__", repr(model.per_factor(res))).replace("__JS__", JS_MODEL))
