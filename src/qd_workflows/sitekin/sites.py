"""Step 3, pick: one site per symmetry-distinct environment of each builder role, its CdCl2 unit and three states.

Environment: two surface Cd are equivalent when their surroundings in the as-built dot (exact lattice positions plus
the placed Cl) agree, i.e. the same sorted list of (element, distance in 0.05 A bins) of every atom within r_env.
Each environment is computed once; its Cd count n is its weight. Per environment the first Cd (in dot order) that
carries a whole unit is the site; an environment with no whole unit is listed, not dropped.
Unit: the Cd and two of its Cl. Cl bound to that Cd alone first, then the outermost shared one; where it has only one,
a bridge or atop Cl within borrow_A is borrowed (not on facet:111B). Rejected when it would leave a Cl with no Cd.
The Se goes on the more strongly bound Cl (anion > bridge > atop, then the nearer; on facet:111B the one nearest the
role's middle). States (Delta k = -1, no CdCl2 on the new Se):
  ref    the relaxed dot, unit on the site            atoms: context + [Cd, Cl, Cl]
  empty  the unit removed (dot*)                       atoms: context
  M      the unit's Cd kept, one Cl -> Se, other gone  atoms: context + [Cd, Se]
local: context atoms within r_local of the unit, the partial-Hessian set of every state.
Writes sites/<tag>.json (the job of the states step) and sites/<tag>_summary.json (roles, environments, sites)."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from . import recipe as rc

RANK = {"anion": 3, "bridge": 2, "atop": 1, "loose": 0}


def site_name(role: str, env: str, cd: int) -> str:
    return f"{role.replace(':', '_').replace('|', '_').replace('+', '_')}_{env.replace('+', 'p')}_cd{cd}"


def pick(dot: dict, ref: dict, r_env: float = 4.5, r_local: float = 6.0, borrow: float = 5.5, cl: str = "Cl",
         cation: str = "Cd", anion: str = "Se"):
    """Sites of one dot. dot: the build (symbols, as-built positions, roles); ref: the relaxed reference (positions,
    E). Returns (job, summary rows)."""
    s = np.array(dot["symbols"])
    x = np.array(ref["positions"], float)
    roles = dot["roles"]
    t = cKDTree(x)
    nb = [[j for j in t.query_ball_point(x[i], 3.2) if j != i] for i in range(len(s))]
    cdcl = {int(k): [j for j in nb[k] if s[j] == cation] for k in np.flatnonzero(s == cl)}
    kind = {k: ("anion" if len(c) >= 3 else "bridge" if len(c) == 2 else "atop" if c else "loose") for k, c in cdcl.items()}
    ctr = x.mean(0)
    cds = [int(i) for i in np.flatnonzero(s == cation) if roles[i] != "interior"]
    edgey = np.array([i for i in cds if not roles[i].startswith("facet")])
    te = cKDTree(x[edgey]) if len(edgey) else None

    def env(i):
        c = sorted((k for k, v in cdcl.items() if i in v), key=lambda k: -RANK[kind[k]])
        return f"{sum(s[j] == anion for j in nb[i])}{anion}" + "".join(f"+{kind[k]}" for k in c)

    def depth_bin(i):
        return 0 if not roles[i].startswith("facet") or te is None else int(float(te.query(x[i])[0]) // 2)

    def unit_of(cd):
        own = [k for k, c in cdcl.items() if cd in c]
        only = [k for k in own if cdcl[k] == [cd]]
        rest = sorted((k for k in own if k not in only), key=lambda k: -((x[k] - ctr) @ (x[cd] - ctr)))
        got = (only + rest)[:2]
        if len(got) < 2 and not roles[cd].startswith("facet:111B"):
            extra = sorted((k for k, c in cdcl.items() if k not in got and kind[k] in ("bridge", "atop") and cd not in c
                            and np.linalg.norm(x[k] - x[cd]) < borrow), key=lambda k: np.linalg.norm(x[k] - x[cd]))
            got = got + extra[:1]
        gone = {cd, *got}
        ok = len(got) == 2 and not [k for k, c in cdcl.items() if k not in gone and set(c) <= gone]
        return got if ok else None

    x_b = np.array(dot["positions"], float)                      # environments on the as-built dot
    tb = cKDTree(x_b)

    def fingerprint(i):
        return tuple(sorted((str(s[j]), round(float(np.linalg.norm(x_b[j] - x_b[i])) / 0.05))
                            for j in tb.query_ball_point(x_b[i], r_env) if j != i))

    job = {"dot": dot["tag"], "size": dot["size"], "surface": dot["surface"], "reference": "relaxed",
           "E_ref": ref["E"], "sites": {}}
    table = []
    for role in sorted({roles[i] for i in cds}):
        members = [i for i in cds if roles[i] == role]
        whole = {i: u for i in members if (u := unit_of(i))}
        envs = {}
        for i in members:
            envs.setdefault(fingerprint(i), []).append(i)
        chosen, rows = [], []
        for k, (_f, grp) in enumerate(sorted(envs.items(), key=lambda kv: -len(kv[1]))):
            ok = [i for i in grp if i in whole]
            site = ok[0] if ok else None
            rows.append({"env_id": k, "n": len(grp), "cd": grp, "site": site, "no_whole_unit": site is None,
                         "site_name": site_name(role, env(site), site) if site is not None else None})
            if site is not None:
                chosen.append(site)
        table.append({"role": role, "n": len(members), "whole_units": len(whole), "environments": rows,
                      "sites": [{"cd": i, "env": env(i), "depth_bin": depth_bin(i)} for i in chosen]})
        mid = x[members].mean(0)
        for i in chosen:
            unit_cl = whole[i]
            if role.startswith("facet:111B"):
                se_at = min(unit_cl, key=lambda k: np.linalg.norm(x[k] - mid))
            else:
                se_at = max(unit_cl, key=lambda k: (RANK[kind[k]], -np.linalg.norm(x[k] - x[i])))
            removed = [k for k in unit_cl if k != se_at][0]
            unit = [i, *unit_cl]
            ctx = [j for j in range(len(s)) if j not in set(unit)]
            sc, xc = s[ctx].tolist(), x[ctx]
            tc = cKDTree(xc)
            local = sorted({j for u in unit for j in tc.query_ball_point(x[u], r_local)})
            states = {
                "ref": {"symbols": sc + [cation, cl, cl], "positions": np.vstack([xc, x[unit]]).round(5).tolist()},
                "empty": {"symbols": sc, "positions": xc.round(5).tolist()},
                "M": {"symbols": sc + [cation, anion], "positions": np.vstack([xc, x[[i, se_at]]]).round(5).tolist()},
            }
            job["sites"][site_name(role, env(i), i)] = {
                "role": role, "env": env(i), "cd": i, "cl": unit_cl, "se_at": se_at, "cl_removed": removed,
                "cl_kinds": [kind[k] for k in unit_cl], "n_context": len(ctx), "local": local, "states": states}
    return job, table


def run(r: dict, out: Path, sizes=None, log=print) -> list[Path]:
    out = Path(out)
    (out / "sites").mkdir(parents=True, exist_ok=True)
    paths = []
    for e in r["sizes"]:
        t = rc.tag(e["size"])
        if sizes and t not in sizes:
            continue
        dot = json.loads((out / "dots" / f"{t}.json").read_text())
        ref = json.loads((out / "reference" / f"{t}.json").read_text())
        sp = r["sites"]
        job, table = pick(dot, ref, float(sp["r_env"]), float(sp["r_local"]), float(e.get("borrow_A", sp["borrow_A"])),
                          r["ligand"]["anion"], r["material"]["cation"], r["material"]["anion"])
        p = out / "sites" / f"{t}.json"
        p.write_text(json.dumps(job))
        (out / "sites" / f"{t}_summary.json").write_text(json.dumps(table, indent=1))
        log(f"== {t}: {len(dot['symbols'])} atoms, {len(job['sites'])} sites")
        for row in table:
            log(f"   {row['role']:26s} Cd {row['n']:3d}  whole units {row['whole_units']:3d}  environments "
                + "  ".join(f"[{v['n']}: {'NO UNIT' if v['no_whole_unit'] else v['site_name']}]" for v in row["environments"]))
        paths.append(p)
    return paths
