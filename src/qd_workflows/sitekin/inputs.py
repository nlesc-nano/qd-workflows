"""Step 6, inputs: the role table and the NEB results in one file for the size law, computed values only.

kinetics_inputs.json = {about, T_K, role_table (role_table.json), neb: {size: {role: {k_disp, k_dir}}}}. A role
without a NEB result is absent, never filled; the velocity form of k_eff (k1 x site density with precursor saturation,
or the second-order rate from dG_TS vs 1 M, capped at the collision rate) is left to the growth model."""
from __future__ import annotations

import json
from pathlib import Path

ABOUT = ("Site-kinetics protocol (qdw sites). role_table: per size, builder role and symmetry-distinct environment, "
         "the Cd count n, bind = dG_bind,X (eV, 1 M), K_X = exp(-bind/kT) (1/M) and dmu = Delta mu_f (eV, a_X = 1, "
         "Delta k = -1). neb: per size and role, k_disp and k_dir from CI-NEB on the role's largest environment, bare "
         "CdSe: dG_P_1M and dG_TS_1M against G(site) + G(CdSe, 1 M), dG_act = G(TS) - G(P), k1 = (kT/h) e^(-dG_act/kT); "
         "no_barrier = collision-limited; converged_CI false = TS from the refined loose window (about 0.1 eV).")


def run(out: Path, log=print) -> Path:
    out = Path(out)
    table = json.loads((out / "role_table.json").read_text())
    neb = {}
    for res in sorted((out / "neb").glob("*/*/*/result.json")):
        R = json.loads(res.read_text())
        if "result" not in R:
            continue
        ep = json.loads((res.parent / "endpoints.json").read_text())["meta"]
        size = ep["size_tag"].replace("p", ".")
        entry = {"site": ep["site"], **R["result"]}
        neb.setdefault(size, {}).setdefault(ep["role"], {})[f"k_{ep['kind']}"] = entry
    data = {"about": ABOUT, "T_K": table["T_K"], "role_table": table, "neb": neb}
    p = out / "kinetics_inputs.json"
    p.write_text(json.dumps(data, indent=1))
    for size, roles in neb.items():
        for role, d in roles.items():
            for k, v in sorted(d.items()):
                txt = "no barrier" if v["no_barrier"] else (f"dG_act {v['dG_act']:+.2f} eV, k1 {v['k1_s']:.1e} /s, "
                                                             f"TS {v['dG_TS_1M']:+.2f} eV vs 1 M")
                log(f"{size} {role:26s} {k:6s} P {v['dG_P_1M']:+.2f} eV vs 1 M; {txt}")
    return p
