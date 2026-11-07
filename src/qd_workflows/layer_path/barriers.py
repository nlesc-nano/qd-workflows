"""
2D-island barriers from the relaxed builder islands and the registry's ladder rungs (the scripts island_barriers.py,
island_types.py and island_cdse.py).

island_barriers: per parent and facet, Delta Omega(n) = G_island - G_parent - dSe mu_CdSe - (dCd - dSe) mu_CdCl2 with
n = added Cd; the lowest state at each n (islands, and the registry's one-sheet rungs of the same parent and family)
gives the curve. Barrier = the highest point of the curve before it first drops below the parent (the classical
nucleation barrier on the lowest-state curve; it assumes a step from each size to the next). Compared with ladder_fit.

island_types: the island states split by what they add, dm = dSe CdSe units, dk = dCd - dSe CdCl2 units and lattice
Cl -> Se swaps (dCl < 2 dk); Omega at alpha = 0.2. Writes island_states.json.

island_cdse: the CdSe-growth barrier per parent and facet with the shell let to equilibrate at each crystal size:
Omega*(dm) = lowest state with dm added CdSe (any dk) from the parent, its islands and the rungs, measured from
Omega*(0) (the parent or a better shell of it); barrier = highest Omega* before it first drops below Omega*(0).
"""
from __future__ import annotations

import collections
import json
from pathlib import Path

import numpy as np

from .relax import g_model
from .steps import POOLS


def island_states(net, islands_dir, results, g_cache, thermo, params, states_files):
    """Relaxed pure island states with 1 <= dCd <= 10 and dG = G_island - G_parent."""
    reg = json.loads((Path(net) / "registry.json").read_text())
    byid = {v["id"]: k for k, v in reg["dots"].items()}
    res = json.loads(Path(results).read_text())
    gc = json.loads(Path(g_cache).read_text())
    T, g_of = g_model(thermo, params)

    def g_dot(dk):
        s = res[dk]
        ck = f"{T}|{dk}|{s['E_eV']}"
        return gc[ck] if ck in gc else g_of(s)
    st = {}
    for f in states_files:
        st.update(json.loads(Path(f).read_text())["states"])
    isl = []
    for r in json.loads((Path(islands_dir) / "states_index.json").read_text()):
        if r.get("status") not in ("todo", "relaxed") or not r.get("pure") or not 1 <= r.get("dCd", 0) <= 10:
            continue
        s = st.get(r.get("name")) if r["status"] == "todo" else res.get(r["key"])
        if s is None:
            continue
        pk = byid[r["parent"]]
        if pk not in res:
            continue
        isl.append({**{k: r[k] for k in ("tag", "parent", "family", "dSe", "dCd", "dCl")}, "dG": g_of(s) - g_dot(pk),
                    "n_imag": s.get("n_imaginary", 0), "converged": s.get("converged", True)})
    print(f"{len(isl)} relaxed island states ({sum(x['n_imag'] > 0 for x in isl)} with imaginary modes, "
          f"{sum(not x['converged'] for x in isl)} not converged)")
    return isl


def island_barriers(isl, ladders_file, ladder_fit_file, out, pools=POOLS):
    lad = json.loads(Path(ladders_file).read_text())["ladders"]
    lf = json.loads(Path(ladder_fit_file).read_text())
    res = {}
    for al, (mm, mk) in pools.items():
        fitd = {(int(x["a"]), f): x for f, Lf in lf[str(al)].items() for x in Lf}
        by = {}
        for x in isl:
            om = x["dG"] - x["dSe"] * mm - (x["dCd"] - x["dSe"]) * mk
            by.setdefault((x["parent"], x["family"]), {}).setdefault(x["dCd"], []).append(om)
        rows = []
        for (p, f), d in by.items():
            curve = {n: min(v) for n, v in d.items()}
            for v in lad.get(f"{p}|{f}", {}).values():
                n = v["dm"] + v["dk"]
                curve[n] = min(curve.get(n, np.inf), v["dG"] - v["dm"] * mm - v["dk"] * mk)
            ns = sorted(curve)
            top, first_down = 0.0, None
            for n in ns:
                if curve[n] < 0:
                    first_down = n
                    break
                top = max(top, curve[n])
            F = fitd.get((p, f))
            rows.append({"parent": p, "family": f, "curve": {str(n): curve[n] for n in ns}, "n_island_sizes": len(d),
                         "barrier": top if first_down is not None else None, "n_crit": first_down,
                         "barrier_fit": F["barrier_fit"] if F else None, "s_star_fit": F["s_star"] if F else None,
                         "m_a": F["m_a"] if F else None})
        res[al] = rows
        print(f"alpha = {al}:")
        for f in ("111B", "100", "111A"):
            R = [r for r in rows if r["family"] == f]
            if not R:
                continue
            B = [r["barrier"] for r in R if r["barrier"] is not None]
            pr = [(r["barrier"], r["barrier_fit"]) for r in R if r["barrier"] is not None
                  and r["barrier_fit"] is not None and np.isfinite(r["barrier_fit"])]
            print(f"  {f}: {len(R)} parents, {len(B)} reach a downhill state; barrier median "
                  f"{np.median(B) if B else np.nan:.2f} eV (min {min(B) if B else np.nan:.2f}, "
                  f"max {max(B) if B else np.nan:.2f})"
                  + (f"; vs ladder fit {np.median([q for _p, q in pr]):.2f} eV, median diff "
                     f"{np.median([a - q for a, q in pr]):+.2f}" if pr else ""))
            for r in sorted(R, key=lambda r: r["m_a"] or 0):
                c = " ".join(f"{n}:{v:+.2f}" for n, v in list(r["curve"].items())[:8])
                print(f"     p{r['parent']} m_a {r['m_a']}: barrier "
                      f"{r['barrier'] if r['barrier'] is None else round(r['barrier'], 2)} (fit "
                      f"{r['barrier_fit'] if r['barrier_fit'] is None else round(r['barrier_fit'], 2)}), "
                      f"n_crit {r['n_crit']} | {c}")
    out_ = {str(k): v for k, v in res.items()}
    if out:
        Path(out).write_text(json.dumps(out_))
    return out_


def island_types(isl, out="island_states.json", alpha=0.2, pools=POOLS):
    mm, mk = pools[alpha]
    for x in isl:
        x["dm"], x["dk"] = x["dSe"], x["dCd"] - x["dSe"]
        x["swap"] = 2 * x["dk"] - x["dCl"]          # lattice Cl turned into Se
        x["om"] = x["dG"] - x["dm"] * mm - x["dk"] * mk
    for f in ("111B", "100", "111A"):
        X = [x for x in isl if x["family"] == f]
        print(f"{f}: {len(X)} states")
        c = collections.defaultdict(list)
        for x in X:
            c[(x["dm"], x["dk"], x["dCl"])].append(x["om"])
        for k, v in sorted(c.items(), key=lambda kv: (kv[0][0] + kv[0][1], kv[0])):
            if len(v) >= 3:
                print(f"   +{k[0]} CdSe +{k[1]} CdCl2 dCl {k[2]:+d}: n {len(v):3d}, Omega median {np.median(v):+.2f}, "
                      f"min {min(v):+.2f} eV")
    Path(out).write_text(json.dumps(isl))
    return isl


def island_cdse(island_states_file, ladders_file, ladder_fit_file, out, pools=POOLS, show_alpha=0.2):
    isl = json.loads(Path(island_states_file).read_text())
    lad = json.loads(Path(ladders_file).read_text())["ladders"]
    lf = json.loads(Path(ladder_fit_file).read_text())
    res_all = {}
    for al, (mm, mk) in pools.items():
        fitd = {(int(x["a"]), f): x for f, Lf in lf[str(al)].items() for x in Lf}
        by = collections.defaultdict(lambda: collections.defaultdict(list))
        for x in isl:
            by[(x["parent"], x["family"])][x["dSe"]].append(
                ("isl", x["dCd"] - x["dSe"], x["dG"] - x["dSe"] * mm - (x["dCd"] - x["dSe"]) * mk))
        for key in list(by):
            by[key][0].append(("parent", 0, 0.0))
            for v in lad.get(f"{key[0]}|{key[1]}", {}).values():
                by[key][v["dm"]].append(("rung", v["dk"], v["dG"] - v["dm"] * mm - v["dk"] * mk))
        res = []
        for (p, f), d in by.items():
            ref = min(o for _t, _k, o in d[0])
            cur = {dm: min(o for _t, _k, o in v) - ref for dm, v in d.items() if dm >= 0}
            top, down = 0.0, None
            for dm in sorted(cur):
                if dm == 0:
                    continue
                if cur[dm] < 0:
                    down = dm
                    break
                top = max(top, cur[dm])
            F = fitd.get((p, f))
            res.append({"parent": p, "family": f, "m_a": F["m_a"] if F else None, "shell_gain": ref,
                        "curve": {str(k): cur[k] for k in sorted(cur)}, "barrier": top if down else None,
                        "dm_crit": down, "fit": F["barrier_fit"] if F else None,
                        "sources": {str(k): sorted({t for t, _k, _o in d[k]}) for k in sorted(cur)}})
        res_all[str(al)] = res
        print(f"alpha = {al}")
        for f in ("111B", "100", "111A"):
            R = sorted([r for r in res if r["family"] == f], key=lambda r: r["m_a"] or 0)
            B = [r["barrier"] for r in R if r["barrier"] is not None]
            print(f"  {f}: {len(R)} parents, {len(B)} cross; CdSe barrier median "
                  f"{np.median(B) if B else float('nan'):.2f} eV, shell gain median "
                  f"{np.median([r['shell_gain'] for r in R]) if R else float('nan'):+.2f} eV")
            if al == show_alpha:
                for r in R:
                    c = " ".join(f"{k}:{v:+.2f}" for k, v in list(r["curve"].items())[1:8])
                    b = "none" if r["barrier"] is None else f"{r['barrier']:.2f}"
                    fit = "-" if r["fit"] is None else f"{r['fit']:.2f}"
                    print(f"     p{r['parent']} m {r['m_a']}: shell {r['shell_gain']:+.2f}, barrier {b} (fit {fit}), "
                          f"crit dm {r['dm_crit']} | {c}")
    Path(out).write_text(json.dumps(res_all))
    return res_all
