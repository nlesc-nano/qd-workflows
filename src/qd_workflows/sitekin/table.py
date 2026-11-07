"""Step 4b, table: per size, builder role and environment, the CdCl2 binding and the facet chemical potential.

    dG_bind,X = G(ref) - G(empty) - G(CdCl2, 1 M),    K_X = e^(-dG_bind,X/kT)  (M^-1)
    Delta mu_f = G(M) - G(ref) + G(CdCl2, 1 M) - g_bulk          (a_X = 1, Delta k = -1; + kT ln a_X otherwise)
Each environment carries its Cd count n (its weight); per role the n-weighted K_X and Delta mu_f, and the Cd of
environments without a whole unit (no channel). Same fields as the 2026-10 role tables (n, dmu, bind, sites).
Writes role_table.json."""
from __future__ import annotations

import json
import math
from pathlib import Path

from . import recipe as rc


def size_table(summary: list, states: dict, g_bulk: float, kt: float) -> dict:
    gx = states["CdCl2"]["G"]
    out = {}
    for row in summary:
        envs, missing = [], 0
        for v in row["environments"]:
            if v["no_whole_unit"]:
                missing += v["n"]
                continue
            st = states["sites"].get(v["site_name"])
            if st is None or "G" not in st:
                raise KeyError(f"no states for {v['site_name']}")
            bind = st["dG_bind_X"]
            dmu = st["G"]["M"] - st["G"]["ref"] + gx - g_bulk
            envs.append({"n": v["n"], "sites": [v["site_name"]], "bind": bind, "K_X": math.exp(-bind / kt), "dmu": dmu})
        if not envs:
            out[row["role"]] = {"envs": [], "n_cd": row["n"], "cd_no_unit": missing}
            continue
        n = sum(e["n"] for e in envs)
        out[row["role"]] = {
            "envs": envs, "n_cd": n, "cd_no_unit": missing, "n_env": len(envs),
            "K_X_mean": sum(e["n"] * e["K_X"] for e in envs) / n,
            "bind_mean": sum(e["n"] * e["bind"] for e in envs) / n,
            "dmu_mean": sum(e["n"] * e["dmu"] for e in envs) / n,
            "dmu_min": min(e["dmu"] for e in envs), "dmu_max": max(e["dmu"] for e in envs)}
    return out


def run(r: dict, out: Path, log=print) -> Path:
    out = Path(out)
    kt, gb = rc.kt(r), float(r["material"]["g_bulk_eV"])
    table = {"T_K": float(r["T_K"]), "g_bulk_eV": gb, "standard_state_M": 1.0, "sizes": {}}
    for e in r["sizes"]:
        t = rc.tag(e["size"])
        sp, stp = out / "sites" / f"{t}_summary.json", out / "states" / f"{t}.json"
        if not (sp.exists() and stp.exists()):
            log(f"{t}: no sites or states yet, left out")
            continue
        st = json.loads(stp.read_text())
        if not st.get("done"):
            log(f"{t}: states not finished (rerun the states step, it resumes), left out")
            continue
        table["sizes"][f"{float(e['size']):.1f}"] = size_table(json.loads(sp.read_text()), st, gb, kt)
    p = out / "role_table.json"
    p.write_text(json.dumps(table, indent=1))
    for size, roles in table["sizes"].items():
        for role, v in roles.items():
            if v["envs"]:
                log(f"{size} {role:26s} Cd {v['n_cd']:3d} in {v['n_env']} env.  dG_bind {v['bind_mean']:+.2f}  "
                    f"K_X {v['K_X_mean']:.2e} /M  dmu {v['dmu_mean']:+.2f} [{v['dmu_min']:+.2f}, {v['dmu_max']:+.2f}] eV")
            else:
                log(f"{size} {role:26s} no whole CdCl2 unit ({v['cd_no_unit']} Cd)")
    return p
