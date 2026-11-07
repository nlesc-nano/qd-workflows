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


def roles_block(kinetics_inputs: dict, m_of_size: dict) -> dict:
    """The growth template's per-role block (kinetics.roles) from kinetics_inputs.json: envs {role: {m: [[n, dmu,
    bind], ...]}} (environments with a whole unit) and neb {role: {dir, disp}} with dG_TS/dG_P against 1 M or
    no_barrier, from the NEB size (one size per role). m_of_size: recipe size ("2.5") -> CdSe units of that dot."""
    rt = kinetics_inputs["role_table"]
    envs = {}
    for size, roles in rt["sizes"].items():
        m = int(m_of_size[size])
        for role, v in roles.items():
            if v.get("envs"):
                envs.setdefault(role, {})[m] = [[e["n"], round(float(e["dmu"]), 4), round(float(e["bind"]), 4)]
                                                for e in v["envs"]]
    neb = {}
    for _size, roles in (kinetics_inputs.get("neb") or {}).items():
        for role, d in roles.items():
            q = {}
            for kind in ("dir", "disp"):
                r = d.get(f"k_{kind}")
                if r is None:
                    continue
                q[kind] = ({"no_barrier": True} if r.get("no_barrier")
                           else {"dG_P": round(float(r["dG_P_1M"]), 4), "dG_TS": round(float(r["dG_TS_1M"]), 4)})
            if q:
                neb[role] = q
    return {"T_K": kinetics_inputs.get("T_K"), "rate_form": None, "envs": envs, "neb": neb}


def facet_area_block(out: Path, start_h: dict, builder: dict, rate_form: str = "first_order") -> dict:
    """kinetics.facet_area for the growth template from a qdw sites run: kinetics_inputs.json, with m of each size
    from dots/<tag>.json (the anion count of the dot the sites were computed on), via growth.facet_area_from_roles
    (terrace binding pooled over sizes, terrace Delta mu per size, facet and edge NEBs)."""
    from ..growth import facet_area_from_roles
    out = Path(out)
    ki = json.loads((out / "kinetics_inputs.json").read_text())
    anion = builder.get("anion", "Se")
    m_of = {}
    for size in ki["role_table"]["sizes"]:
        d = json.loads((out / "dots" / f"{size.replace('.', 'p')}.json").read_text())
        m_of[size] = int(d["composition"].get(anion, 0))
    return facet_area_from_roles(roles_block(ki, m_of), start_h, builder, rate_form)
