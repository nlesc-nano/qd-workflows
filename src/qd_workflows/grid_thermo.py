"""
Thermochemistry of an [MA]_m(MX_q)_k grid (lattice_grid) and its grand potential (cd16_lattice_grid/grand_potential.md).

    hessians  (GPU)  MACE analytic Hessian of every kept structure of every node (nodes/m<m>_k<k>[_i].xyz), its
                     rotational symmetry number, and G°(T) at 1 M as an ideal gas (Orchestr.AI ideal_gas_g:
                     translation, rigid rotor, quasi-RRHO vibrations; the treatment of the monomer library), plus
                     the library monomers MX_q and MA: thermo.json, checkpointed after every node.
    omega     (CPU)  at a temperature, initial precursor and MX_q concentrations and a precursor-to-monomer
                     conversion: G(m,k) (sum over the node's isomers), the chemical potentials of MA (supply alpha c_P0,
                     from the complexes MA(MX_q)_n = the grid's m = 1 row, n = 0 the MA monomer when the row lacks it)
                     and of MX_q (what monomer formation left, shared with the complexes; optionally ligand-bound), the grand
                     potential Omega(m,k), the best shell k*(m), P(k|m), <k>_m, Omega(m), the steps along both axes
                     and the formation free energy of the dot the grid started from.

    python -m qd_workflows.grid_thermo hessians <grid dir> --model <mace model> [--refs <refs dir>]
    python -m qd_workflows.grid_thermo omega <grid dir>/thermo.json --T 300 --c-p0 1e-2 --c-mx0 5e-2 --alpha 0.1 [--k-l .. --c-l ..]

Isomers without a Hessian (the node's other evaluated configurations, energies in grid.json) take the thermal part
of the node's best structure: G_i = E_i + (G_best - E_best).
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np

KB_EV = 8.617333262e-5
TEMPERATURES = [float(t) for t in range(50, 801, 25)]    # the monomer library's grid


def _read_xyz(path):
    lines = Path(path).read_text().splitlines()
    n = int(lines[0])
    rows = [l.split() for l in lines[2:2 + n]]
    return [r[0] for r in rows], np.array([[float(x) for x in r[1:4]] for r in rows])


def _gas(symbols, pts, energy, freqs, temperatures):
    from ase import Atoms
    from orchestr_ai.qd.references import ideal_gas_g, rotational_symmetry_number
    sigma = rotational_symmetry_number(symbols, pts)
    atoms = Atoms(list(symbols), positions=np.asarray(pts, float))
    freqs = np.asarray(freqs, float)
    if len(freqs) and len(freqs) < 3 * len(symbols) - 6:
        # partial Hessian (carboxylate inner region, ligands.frequencies): translation, rotation and pV of the whole
        # species, plus the quasi-RRHO F_vib of the inner modes (the full path's sum, with fewer modes)
        from orchestr_ai.qd.steps.hessian import thermo, thermo_modes
        g0 = ideal_gas_g(atoms, energy, np.zeros(0), temperatures, sigma)
        fv = thermo(thermo_modes(freqs), temperatures)["F_vib_eV"]
        g = [a + b for a, b in zip(g0, fv)]
    else:
        g = ideal_gas_g(atoms, energy, freqs, temperatures, sigma)
    return sigma, [round(float(x), 6) for x in g]


def _monomers(refs, head, names, temperatures, vibrations=True):
    """The library monomers (MACE energies and frequencies), G° recomputed on `temperatures` (without vibrations
    when the clusters have none either)."""
    out = {}
    root = Path(refs) / head / "monomers"
    for name in names:
        hits = sorted(root.glob(f"{name}-*/monomer.json"))
        if not hits:
            continue
        m = json.loads(hits[0].read_text())
        sigma, g = _gas(m["symbols"], m["positions"], m["energy_eV"],
                        m["frequencies_cm1"] if vibrations else [], temperatures)
        out[name] = {"E_eV": m["energy_eV"], "sigma": sigma, "G_eV": g, "frequencies_cm1": m["frequencies_cm1"],
                     "source": str(hits[0])}
    return out


def unit_names(nodes):
    """[MX_q, MA] monomer names from the node formulas (M first, A second, X third, as the builder orders them;
    q = n_X / (n_M - n_A))."""
    import re
    for n in nodes:
        parts = re.findall(r"([A-Z][a-z]?)(\d*)", n["formula"])
        if len(parts) == 3:
            (m, nm), (a, na), (x, nx) = [(e, int(c or 1)) for e, c in parts]
            q = nx // (nm - na)
            return [f"{m}{x}{q if q > 1 else ''}", f"{m}{a}"]
    raise ValueError("no node with M, A and X to name the units")


# ---------------- hessians ----------------
def hessians(grid_dir, settings, refs=None, temperatures=TEMPERATURES, log=print, vibrations=True, out_path=None,
              only=None):
    """vibrations=False: no Hessian, G° = E + translation + rotation (+ pV) for clusters and monomers alike
    (a first look before the Hessians; CPU only). only: the (m, k) nodes to treat (nucleation: those the row anneal
    selected); the others are left out of thermo.json and listed under `skipped`."""
    grid_dir = Path(grid_dir)
    g = json.loads((grid_dir / "grid.json").read_text())
    ref_names = unit_names(g["nodes"])
    if only is not None:
        only = {tuple(x) for x in only}
        g = g | {"nodes": [n for n in g["nodes"] if (n["m"], n["k"]) in only],
                 "skipped": [[n["m"], n["k"]] for n in g["nodes"] if (n["m"], n["k"]) not in only]}
    if vibrations:
        from orchestr_ai.qd.steps import detachment as dt
        work = grid_dir / "thermo"
        work.mkdir(exist_ok=True)
        cache = dt.RelaxCache(SimpleNamespace(props=work, settings=settings, record_dir=grid_dir))
    out_path = Path(out_path) if out_path else grid_dir / ("thermo.json" if vibrations else "thermo_novib.json")
    old = json.loads(out_path.read_text()) if out_path.is_file() else {}
    done = {(n["m"], n["k"]): n for n in old.get("nodes", [])}
    nodes = []
    t0 = time.time()
    for n in g["nodes"]:
        key = (n["m"], n["k"])
        if key in done and len(done[key]["states"]) == len(n["states"]):
            nodes.append(done[key])
            continue
        states = []
        for i, st in enumerate(n["states"]):
            tag = f"m{n['m']:02d}_k{n['k']:02d}{'' if i == 0 else f'_{i}'}"
            syms, pts = _read_xyz(grid_dir / "nodes" / f"{tag}.xyz")
            if vibrations and getattr(settings, "ligand", None):      # carboxylate X: partial Hessian (ligands)
                freqs = np.asarray(_partial_frequencies(cache, syms, pts, settings))
            else:
                freqs = np.asarray(cache.frequencies(syms, pts, settings)) if vibrations else np.zeros(0)
            sigma, gt = _gas(syms, pts, st["E_eV"], freqs, temperatures)
            states.append({"file": f"nodes/{tag}.xyz", "E_eV": st["E_eV"], "fingerprint": st["fingerprint"],
                           "sigma": sigma, "n_imaginary": int((freqs < -10.0).sum()),
                           "lowest_cm1": round(float(freqs.min()), 2) if freqs.size else None,
                           "frequencies_cm1": freqs.round(3).tolist(), "G_eV": gt,
                           **({"source": st["source"]} if "source" in st else {})})
        nodes.append({"m": n["m"], "k": n["k"], "formula": n["formula"], "ensemble_E_eV": n["ensemble_E_eV"],
                      "states": states})
        log(f"[grid-thermo] m={n['m']} k={n['k']}: {len(states)} Hessian(s), "
            f"{sum(s['n_imaginary'] for s in states)} imaginary modes, {time.time() - t0:.0f} s")
        _write(out_path, grid_dir, temperatures, nodes, None, False, g, vibrations)
    references = _monomers(refs, settings.head, ref_names, temperatures, vibrations) if refs else {}
    _write(out_path, grid_dir, temperatures, nodes, references, len(nodes) == len(g["nodes"]), g, vibrations)
    return out_path


def _partial_frequencies(cache, syms, pts, settings):
    """Frequencies of the inner region (qd_workflows.ligands.frequencies), cached in the relaxation cache."""
    from orchestr_ai.qd.steps import detachment as dt
    from .ligands import frequencies
    key = f"pfreq|{cache.tag}|{settings.ligand.get('inner', 'cbeta')}|{dt.geometry_key(syms, pts)}"
    if key not in cache.data:
        cache.data[key] = np.asarray(frequencies(syms, pts, settings)).round(4).tolist()
        cache._save()
    return cache.data[key]


def _write(path, grid_dir, temperatures, nodes, references, complete, g, vibrations=True):
    out = {"grid": str(grid_dir), "grid_complete": g.get("complete"), "complete": complete,
           "vibrations": vibrations, "temperatures": temperatures, "references": references or {}, "nodes": nodes,
           **({"skipped": g["skipped"]} if g.get("skipped") else {})}
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(out))
    os.replace(tmp, path)


# ---------------- omega ----------------
def _at(t, temperatures, values):
    return float(np.interp(t, temperatures, values))


def node_g(node, t, temperatures, solv=None):
    """G(m,k;T) = -kT ln sum_i exp(-G_i/kT) over the node's isomers: those with a Hessian, and the other
    evaluated configurations with the thermal part of the best one. solv(node, state) -> dG_solv adds the
    implicit solvation of each state (the other configurations take the best state's)."""
    kt = KB_EV * t
    dsol = (lambda st: solv(node, st)) if solv else (lambda st: 0.0)
    # one term per state: the states are distinct minima already (merged by energy), and two of them can share a
    # lattice fingerprint (the lattice and relaxed starts, or two grids), so they are not keyed by it
    gs = [_at(t, temperatures, s["G_eV"]) + dsol(s) for s in node["states"]]
    best = min(node["states"], key=lambda s: s["E_eV"])
    shift = _at(t, temperatures, best["G_eV"]) - best["E_eV"] + dsol(best)
    known = sorted(s["E_eV"] for s in node["states"])
    extra = list(node["ensemble_E_eV"])
    for e in known:                                    # the states are in the ensemble too: drop one copy each
        j = int(np.argmin([abs(x - e) for x in extra])) if extra else None
        if j is not None and abs(extra[j] - e) < 1e-4:
            extra.pop(j)
    allg = np.array(gs + [e + shift for e in extra])
    g0 = allg.min()
    return float(g0 - kt * np.log(np.exp(-(allg - g0) / kt).sum()))


T_REF_CONV = 500.0      # K: temperature at which k_conv is given when an activation energy scales it


def k_conv_at(t, k_ref, ea_ev=0.0, t_ref=T_REF_CONV):
    """Arrhenius scaling of the conversion rate constant from its value at t_ref."""
    return k_ref * np.exp(-ea_ev / KB_EV * (1.0 / t - 1.0 / t_ref))


def converted(time_s, k_conv, c_p0, c_mx0):
    """x(t) = integral of Q = k_conv c_MX c_P dt for MX_q + P=A -> MA + ... (1:1, second order), mol/L:
    c0^2 k t / (1 + c0 k t) for equal initial concentrations, else the closed form; x -> min(c_p0, c_mx0)."""
    d = c_mx0 - c_p0
    kt_ = k_conv * time_s
    if abs(d) * kt_ < 1e-8:
        c0 = 0.5 * (c_mx0 + c_p0)
        return c0 * c0 * kt_ / (1.0 + c0 * kt_)
    e = np.exp(-abs(d) * kt_)                               # stable for either sign of d
    hi, lo = max(c_mx0, c_p0), min(c_mx0, c_p0)
    return hi * lo * (1.0 - e) / (hi - lo * e)


def _speciation(g_n, mu_mx, kt):
    """(ln sum_n exp[-(G_n - n mu_MX)/kT], {n: x_n}, <n>) of the MA(MX_q)_n complexes."""
    a = {n: -(g - n * mu_mx) / kt for n, g in g_n.items()}
    amax = max(a.values())
    lse = amax + np.log(sum(np.exp(v - amax) for v in a.values()))
    x = {n: float(np.exp(a[n] - lse)) for n in sorted(a)}
    return lse, x, sum(n * v for n, v in x.items())


def thermo_at(thermo, t, solvation: dict | None = None):
    """Free energies at one temperature: G(m,k) of every node, G of the MX_q monomer, G of the MA(MX_q)_n precursor
    complexes (the grid's m = 1 row; n = 0 from the library monomer when the row lacks it). solvation: the
    grid_solvation block (eps, ...): every species gets its implicit solvation free energy (none at eps = 1)."""
    temps = thermo["temperatures"]
    nodes = {(n["m"], n["k"]): n for n in thermo["nodes"]}
    eps = float((solvation or {}).get("eps", 1.0))
    if solvation is not None and "complexes" in thermo.get("solvation", {}):
        # a thermo file that carries its own m = 1 complexes' COSMO data (nucleation, seeded) uses those, not the
        # parameters' (which belong to another grid's structures); {} means none: the cluster model
        solvation = {**solvation, "complexes": thermo["solvation"]["complexes"]}
    mx_name, ma_name = unit_names(thermo["nodes"])
    solv, d_mx, d_ma = None, 0.0, 0.0
    if eps > 1.0:
        from . import grid_solvation as gs_
        names = (ma_name, mx_name)
        solv = lambda node, st: gs_.state_shift(node, st, eps, solvation, names)  # noqa: E731
        d_mx = gs_.monomer(mx_name, eps, solvation) if mx_name in solvation["monomers"] else 0.0
        d_ma = gs_.monomer(ma_name, eps, solvation) if ma_name in solvation["monomers"] else 0.0
    G = {key: node_g(n, t, temps, solv) for key, n in nodes.items()}
    refs = thermo["references"]
    ma_name = ma_name if ma_name in refs else None
    if mx_name not in refs:
        raise ValueError("thermo.json has no MX_q monomer reference (run hessians with --refs)")
    g_n = {k: G[(1, k)] for (m, k) in G if m == 1}
    if 0 not in g_n and ma_name:
        g_n[0] = _at(t, temps, refs[ma_name]["G_eV"]) + d_ma
    if not g_n:
        raise ValueError("no MA precursor: the grid has no m = 1 row and thermo.json no MA monomer")
    return {"t": t, "nodes": nodes, "G": G, "g_mx": _at(t, temps, refs[mx_name]["G_eV"]) + d_mx, "g_n": g_n,
            "eps": eps, "dG_solv_mx": d_mx}


def pool(th, c_ma, avail, k_l=0.0, c_l=0.0, c0=1.0, bare=False, gamma_ma=1.0, gamma_mx=1.0, c_sat=None):
    """Chemical potentials of the precursor pool (grand_potential.md section 4b) for any MA supply c_ma and MX_q left
    in solution `avail` (mol/L): the MX_q left is shared between the free (or ligand-bound) MX_q and the MA(MX_q)_n
    complexes, and the free MX_q solves that balance. bare=True: the supply is the uncomplexed MA monomer and all
    the MX_q left is free. th: thermo_at(). gamma_ma, gamma_mx: activity coefficients of the MA species (every
    MA(MX_q)_n alike) and of the free MX_q, mu = mu0 + kT ln(gamma c / c0), 1 = ideal; the complexes bind MX_q at its
    activity, the balance stays in concentrations.
    c_sat: solubility of the free MX_q (mol/L, None = no limit): the free MX_q never exceeds it, and whatever the
    balance leaves above it is solid MX_q (c_MX_solid, mol/L of solution), in equilibrium with the solution (it
    dissolves again as the solution takes MX_q up).
    Returns mu_MA, mu_MX, c_MX_free, a_MX_free (= gamma_mx c_MX_free, what the facets' coverage sees), c_MX_solid,
    speciation, mean_n."""
    kt = KB_EV * th["t"]
    g_n, g_mx = th["g_n"], th["g_mx"]
    if c_ma <= 0 or avail <= 0:
        raise ValueError(f"no monomer or no MX_q left: c_MA = {c_ma:g} M, MX_q left = {avail:g} M")
    lg_ma, lg_mx = float(np.log(gamma_ma)), float(np.log(gamma_mx))
    bind = 1.0 + k_l * (c_l / c0) ** 2
    if bare:
        if 0 not in g_n:
            raise ValueError("bare monomer: no MA monomer (m = 1, k = 0 node or library reference)")
        c_mx_free = avail / bind if c_sat is None else min(avail / bind, c_sat)
        mu_mx = g_mx + kt * (np.log(c_mx_free / c0) + lg_mx)
        return {"mu_MA": float(g_n[0] + kt * (np.log(c_ma / c0) + lg_ma)), "mu_MX": float(mu_mx),
                "c_MX_free": float(c_mx_free), "a_MX_free": float(gamma_mx * c_mx_free),
                "c_MX_solid": max(float(avail - c_mx_free * bind), 0.0), "speciation": {0: 1.0}, "mean_n": 0.0}

    def excess(ln_x):                                   # MX_q balance (in concentrations), increasing in ln_x
        mu = g_mx + kt * (ln_x + lg_mx - np.log(c0))    # the complexes see the free MX_q's activity
        return np.exp(ln_x) * bind + c_ma * _speciation(g_n, mu, kt)[2] - avail

    lo, hi = np.log(avail) - 80.0, np.log(avail / bind)
    if excess(lo) > 0:
        raise ValueError("the MA complexes need more MX_q than is left (no free MX_q solves the balance)")
    solid = 0.0
    if c_sat is not None and np.log(c_sat) < hi and excess(np.log(c_sat)) < 0:
        # saturated: the free MX_q sits at its solubility, the rest of the balance is solid
        solid = float(-excess(np.log(c_sat)))
        c_mx_free = float(c_sat)
    else:
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if not lo < mid < hi:                       # converged to double precision (~60 halvings)
                break
            lo, hi = (mid, hi) if excess(mid) < 0 else (lo, mid)
        c_mx_free = float(np.exp(0.5 * (lo + hi)))
    mu_mx = g_mx + kt * (np.log(c_mx_free / c0) + lg_mx)
    lse, speciation, mean_n = _speciation(g_n, mu_mx, kt)
    return {"mu_MA": float(kt * (np.log(c_ma / c0) + lg_ma) - kt * lse), "mu_MX": float(mu_mx), "c_MX_free": c_mx_free,
            "a_MX_free": float(gamma_mx * c_mx_free), "c_MX_solid": solid, "speciation": speciation, "mean_n": mean_n}


def omega(thermo, t, c_p0, c_mx0, alpha, k_l=0.0, c_l=0.0, c0=1.0, bare=False, c_sat=None):
    """Sections 4-6 of grand_potential.md. Concentrations in mol/L (c0 = 1 M standard state).
    c_p0: initial MA precursor (e.g. the Se precursor); c_mx0: initial MX_q; alpha: precursor-to-monomer conversion
    (from the conversion kinetics: alpha = converted(time, k_conv, c_p0, c_mx0) / c_p0).
    Monomer formation MX_q + P=A -> MA + P X_q supplies c_MA = alpha c_p0 and consumes as much MX_q; the MX_q left is
    shared between the free (or ligand-bound) MX_q and the MA(MX_q)_n complexes: the free MX_q solves that balance.
    bare=True: the whole supply is the uncomplexed MA monomer and all the MX_q left is free (no complexes), the
    limit of the pre-nucleation formula Omega = G_f - m kT ln(c_MA/c0) - k kT ln(c_MX,left/c0)."""
    kt = KB_EV * t
    th = thermo_at(thermo, t)
    nodes, G = th["nodes"], th["G"]
    c_ma = alpha * c_p0
    avail = c_mx0 - c_ma
    pl = pool(th, c_ma, avail, k_l, c_l, c0, bare, c_sat=c_sat)              # c_sat: MX_q solubility, as pool()
    mu_ma, mu_mx, c_mx_free, speciation, mean_n = pl["mu_MA"], pl["mu_MX"], pl["c_MX_free"], pl["speciation"], pl["mean_n"]
    Om = {key: G[key] - key[0] * mu_ma - key[1] * mu_mx for key in G}
    ms = sorted({m for m, _ in Om}, reverse=True)
    per_m = {}
    for m in ms:
        ks = sorted(k for mm, k in Om if mm == m)
        v = np.array([Om[(m, k)] for k in ks])
        w = np.exp(-(v - v.min()) / kt)
        p = w / w.sum()
        per_m[m] = {"k_star": ks[int(v.argmin())], "Omega_min_eV": float(v.min()),
                    "Omega_m_eV": float(v.min() - kt * np.log(w.sum())), "mean_k": float((p * np.array(ks)).sum()),
                    "P_k": {k: float(x) for k, x in zip(ks, p)},
                    "best_structure": nodes[(m, ks[int(v.argmin())])]["states"][0]["file"]}
    steps = {f"{m},{k}": {"dOmega_MA_eV": (Om[(m, k)] - Om[(m - 1, k)]) if (m - 1, k) in Om else None,
                          "dOmega_MX_eV": (Om[(m, k)] - Om[(m, k - 1)]) if (m, k - 1) in Om else None}
             for (m, k) in sorted(Om)}
    return {"T": t, "c_P0_M": c_p0, "c_MX0_M": c_mx0, "alpha": alpha, "c_MA_M": c_ma, "c_MX_left_M": avail,
            "c_MX_in_complexes_M": c_ma * mean_n, "c_MX_free_M": c_mx_free, "K_L": k_l, "c_L_M": c_l,
            "mu_MX_eV": float(mu_mx), "mu_MA_eV": float(mu_ma), "precursor_speciation": speciation, "mean_n": mean_n,
            "G_eV": {f"{m},{k}": v for (m, k), v in sorted(G.items())},
            "Omega_eV": {f"{m},{k}": v for (m, k), v in sorted(Om.items())}, "per_m": per_m, "steps": steps}


def _report(res, start=None):
    print(f"T = {res['T']:.0f} K, precursor {res['c_P0_M']:g} M at conversion {res['alpha']:g}: MA supply "
          f"{res['c_MA_M']:.3g} M; MX {res['c_MX0_M']:g} M, left {res['c_MX_left_M']:.3g} M, in complexes "
          f"{res['c_MX_in_complexes_M']:.3g} M, free {res['c_MX_free_M']:.3g} M: "
          f"mu_MX = {res['mu_MX_eV']:.4f} eV, mu_MA = {res['mu_MA_eV']:.4f} eV")
    print("precursor speciation MA(MX)_n: " + ", ".join(f"n={n}: {x:.3g}" for n, x in res["precursor_speciation"].items()))
    print("   m   k*   Omega_min (eV)   Omega(m) (eV)   <k>    best structure")
    for m, r in sorted(res["per_m"].items(), key=lambda x: -x[0]):
        print(f"  {m:2d}  {r['k_star']:3d}   {r['Omega_min_eV']:12.3f}   {r['Omega_m_eV']:12.3f}   {r['mean_k']:5.2f}   "
              f"{r['best_structure']}")
    if start:
        print(f"Omega{start} = {res['Omega_eV'].get(f'{start[0]},{start[1]}', float('nan')):.3f} eV "
              f"(formation from the precursors at these conditions)")


def page(thermo_path, out, title="Grid Grand Potential", start=(13, 3)) -> str:
    """The interactive grand-potential page (templates/grid_omega.html) with the thermochemistry embedded: the same
    equations as omega(), recomputed in the browser as the sliders move."""
    from . import TEMPLATES
    th = json.loads(Path(thermo_path).read_text())
    mx, ma = unit_names(th["nodes"])
    data = {"title": title, "vibrations": th.get("vibrations", True), "grid_complete": th.get("grid_complete"),
            "temperatures": th["temperatures"], "units": {"MX": mx, "MA": ma}, "start": list(start),
            "references": {k: {"G": v["G_eV"]} for k, v in th["references"].items()},
            "nodes": [{"m": n["m"], "k": n["k"], "formula": n["formula"], "ens": n["ensemble_E_eV"],
                       "states": [{"E": s["E_eV"], "G": s["G_eV"], "file": s["file"], "src": s.get("source", "grid")}
                                  for s in n["states"]]}
                      for n in th["nodes"]]}
    html = (TEMPLATES / "grid" / "omega.html").read_text()
    html = html.replace("{{TITLE}}", title).replace("{{DATA}}", json.dumps(data).replace("</", "<\\/"))
    Path(out).write_text(html)
    return str(out)


def prepare(grid: str, cluster: dict, code: str, *, partition: str | None = None, time_: str = "01:00:00",
            pieces: bool = True) -> dict:
    """GPU job script in <grid>/thermo/ running `hessians` on the grid (the cluster profile's props job)."""
    import copy
    from .props import RUN_PIECES, _job, _seconds
    cl = copy.deepcopy(cluster)
    pp = cl["props"]
    if partition:
        pp["partition"] = partition
    out = Path(grid) / "thermo"
    out.mkdir(parents=True, exist_ok=True)
    run = (f"env PYTHONPATH={code}:$PYTHONPATH python -u -m qd_workflows.grid_thermo hessians {grid} "
           f"--model {pp['model']} --head {pp.get('head', 'omat_pbe')} --refs {pp['refs']}")
    body = RUN_PIECES.format(seconds=max(_seconds(time_) - 60, 60), run=run, job="job.sh", then="") if pieces else run
    (out / "job.sh").write_text(_job(cl, "thermo-" + Path(grid).name, time_, None, body))
    return {"job": str(out / "job.sh"), "partition": pp["partition"], "run": run}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python -m qd_workflows.grid_thermo")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("prepare", help="GPU job for the Hessians of a grid (in <grid>/thermo/)")
    p.add_argument("grid")
    p.add_argument("--cluster", default="default")
    p.add_argument("--code", required=True, help="PYTHONPATH entries on the cluster (qd-workflows src)")
    p.add_argument("--partition")
    p.add_argument("--time", default="01:00:00")
    p = sub.add_parser("hessians", help="Hessians and G°(T) of every kept grid structure (GPU)")
    p.add_argument("grid", help="lattice_grid output directory (grid.json, nodes/)")
    p.add_argument("--model", default=os.environ.get("QDPROPS_MACE_MODEL"))
    p.add_argument("--head", default="omat_pbe")
    p.add_argument("--refs", default=os.environ.get("QDPROPS_REFS"), help="references dir with the monomer library")
    p.add_argument("--no-vibrations", action="store_true",
                   help="no Hessian: translation and rotation only (thermo_novib.json; CPU)")
    p = sub.add_parser("page", help="interactive page: sliders for T, precursor concentrations, conversion")
    p.add_argument("thermo", help="thermo.json (or thermo_novib.json)")
    p.add_argument("--out", required=True, help="HTML file to write")
    p.add_argument("--title", default="Grid Grand Potential")
    p.add_argument("--start", default="13,3", help="m,k of the dot the grid started from")
    p = sub.add_parser("omega", help="grand potential at T and precursor concentrations")
    p.add_argument("thermo", help="thermo.json")
    p.add_argument("--T", type=float, nargs="+", default=[300.0])
    p.add_argument("--c-p0", type=float, default=1e-2, help="initial MA precursor (e.g. Se precursor), M")
    p.add_argument("--c-mx0", type=float, default=5e-2, help="initial MX_q, M")
    p.add_argument("--alpha", type=float, default=0.1, help="precursor-to-monomer conversion (0-1)")
    p.add_argument("--time", type=float, help="reaction time (s): alpha from the conversion kinetics instead")
    p.add_argument("--k-conv", type=float, default=1.0, help=f"conversion rate constant at {T_REF_CONV:g} K, 1/(M s)")
    p.add_argument("--ea", type=float, default=0.0, help="activation energy of the conversion (eV), Arrhenius in T")
    p.add_argument("--bare", action="store_true", help="bare monomer only: no MA(MX_q)_n complexes")
    p.add_argument("--k-l", type=float, default=0.0, help="binding constant of MX_q to two ligands L (dimensionless)")
    p.add_argument("--c-l", type=float, default=0.0, help="ligand concentration, M")
    p.add_argument("--start", default="13,3", help="m,k of the dot the grid started from (formation free energy)")
    p.add_argument("--json", help="write the full result here")
    a = ap.parse_args(argv)
    if a.cmd == "prepare":
        from .cluster import load_cluster
        print(json.dumps(prepare(a.grid, load_cluster(a.cluster), a.code, partition=a.partition, time_=a.time),
                         indent=1))
        return 0
    if a.cmd == "hessians":
        from orchestr_ai.qd.run import Settings
        s = Settings(model=a.model or "none", head=a.head, hessian="analytic")
        out = hessians(a.grid, s, a.refs, log=lambda x: print(x, flush=True), vibrations=not a.no_vibrations)
        print(f"[grid-thermo] {out}")
        return 0
    if a.cmd == "page":
        print(page(a.thermo, a.out, a.title, tuple(int(x) for x in a.start.split(","))))
        return 0
    th = json.loads(Path(a.thermo).read_text())
    results = []
    for t in a.T:
        alpha = a.alpha
        if a.time is not None:
            k = k_conv_at(t, a.k_conv, a.ea)
            alpha = converted(a.time, k, a.c_p0, a.c_mx0) / a.c_p0
            print(f"t = {a.time:g} s, k_conv({t:.0f} K) = {k:.3g} 1/(M s): conversion alpha = {alpha:.4g}")
        r = omega(th, t, a.c_p0, a.c_mx0, alpha, a.k_l, a.c_l, bare=a.bare)
        _report(r, tuple(int(x) for x in a.start.split(",")))
        results.append(r)
    if a.json:
        Path(a.json).write_text(json.dumps(results, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
