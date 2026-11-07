"""
Per-plane QD_builder builds of a layer state (the scripts' build2.build).

Each state goes to the builder as a Wulff cut: one facet per plane (scope: facet) at its distance in A from the origin
(the state's central Se), the planes that carve nothing far out (FAR, faceless families), radius = the closest plane.
The builder does its own role-aware prepass and charge balance (Cl ligand, positive Q made up by adding Cl). With
k_total, the builder's z_type_target sets the number of CdCl2 units (the shells k0 +- 1, +- 2 of the fill).

Needs QD_builder with fix/followups (6561a07: faceless Wulff planes no longer decide {111} site roles).
"""
from __future__ import annotations

import contextlib
import io
from pathlib import Path

import yaml

PRESET = {"charges": {"Cd": 2, "Se": -2, "Cl": -1},
          "passivation": {"ligand": "Cl", "surf_tol": 2.0, "prepass_mode": "role-aware",
                          "prepass_min_cn_terrace": 3, "prepass_min_cn_edge": 2, "prepass_min_cn_vertex": 1},
          "symmetry": {"proper_rotations_only": True}}


def recipe(planes_A, k_total=None) -> dict:
    facets = [{"hkl": list(n), "gamma": round(float(h), 6), "scope": "facet"} for n, h in planes_A]
    rec = dict(PRESET, facets=facets, construction_origin={"center_on_species": "Se"})
    if k_total is not None:
        rec["post_treatment"] = {"z_type_target": {"target_total": int(k_total)}}
    return rec


def build(planes_A, wd, cif, k_total=None) -> bool:
    """Build one state in wd (recipe.yaml, out.xyz, out_cut.xyz, builder.log); True when out.xyz was written."""
    from builder.main import main as builder_main
    wd = Path(wd)
    wd.mkdir(parents=True, exist_ok=True)
    (wd / "recipe.yaml").write_text(yaml.safe_dump(recipe(planes_A, k_total), sort_keys=False))
    log = io.StringIO()
    argv = [str(cif), str(wd / "recipe.yaml"), "-o", str(wd / "out.xyz"), "--positive-q-mode", "add",
            "-r", f"{min(h for _n, h in planes_A):.6f}", "--write-all", "--verbose"]
    try:
        with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
            rc = builder_main(argv)
    except SystemExit as exc:
        rc = exc.code
    (wd / "builder.log").write_text(log.getvalue())
    return rc in (0, None) and (wd / "out.xyz").exists()
