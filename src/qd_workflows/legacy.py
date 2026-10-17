"""Records for the earlier DFT dots (source: legacy), so that they go through the same stages as every other dot.

A legacy dot has no record.json / start.xyz (only the QD_builder writes those). `make_records` takes its DFT
geometry from the ledger (`legacy.geo`) or from the first search template that exists, checks that the atoms match
the formula in the dot's id, writes `<out>/<id>/{record.json,start.xyz}` and lists the ledger's `record`. The
structure is then re-relaxed by MACE (stage 1) like every other dot; the DFT geometry is only the starting point.
Dots without a usable geometry are returned in `missing`: they are rebuilt from the QD_builder library instead.
"""
import json
import re
from collections import Counter
from pathlib import Path

import numpy as np
import yaml

from . import campaign
from .xyz import read_frames, write_xyz

CHARGES = {"Cd": 2, "Zn": 2, "Hg": 2, "Pb": 2, "Al": 3, "Ga": 3, "In": 3,
           "S": -2, "Se": -2, "Te": -2, "P": -3, "As": -3, "Sb": -3, "F": -1, "Cl": -1, "Br": -1, "I": -1}
LIGANDS = ("F", "Cl", "Br", "I")
FAMILY = {"CdS": "II-VI", "CdSe": "II-VI", "CdTe": "II-VI", "HgS": "II-VI", "HgSe": "II-VI", "HgTe": "II-VI",
          "ZnS": "II-VI", "ZnSe": "II-VI", "ZnTe": "II-VI",
          "AlAs": "III-V", "AlP": "III-V", "AlSb": "III-V", "GaAs": "III-V", "GaP": "III-V", "GaSb": "III-V",
          "InAs": "III-V", "InP": "III-V", "InSb": "III-V", "PbS": "IV-VI", "PbSe": "IV-VI", "PbTe": "IV-VI"}
PHASE = {"II-VI": "zinc-blende", "III-V": "zinc-blende", "IV-VI": "rock-salt"}


def formula_of(dot_id: str) -> Counter:
    """Composition from the formula part of a dot id: 'CdSe-Cd-Cd139Se92Cl94-reconstructed' -> Cd 139, Se 92, Cl 94."""
    parts = dot_id.split("-")
    if len(parts) < 3:
        raise ValueError(f"cannot read a formula from {dot_id!r}")
    return Counter({el: int(n) for el, n in re.findall(r"([A-Z][a-z]?)(\d+)", parts[2])})


def find_geometry(d: dict, searches: list[str]) -> Path | None:
    """The DFT geometry of a legacy dot: ledger `legacy.geo`, else the first template (`{id}`, `{stem}`) that exists."""
    geo = (d.get("legacy") or {}).get("geo")
    if geo and Path(geo).is_file():
        return Path(geo)
    stem = re.sub(r"-(reconstructed|clean)$", "", d["id"])        # `{stem}`: the id without the surface suffix
    for t in searches:
        p = Path(t.format(id=d["id"], stem=stem))
        if p.is_file():
            return p
    return None


def _size(symbols, coords, core: set) -> dict:
    xyz = coords[[s in core for s in symbols]]
    axes = np.sort(np.ptp(xyz, axis=0)) / 10.0                      # nm (bounding box of the core)
    return {"n_core": int(len(xyz)), "d_nm": round(float(axes.mean()), 3),
            "axes_nm": [round(float(a), 3) for a in axes], "aspect_ratio": round(float(axes[2] / max(axes[0], 1e-6)), 3)}


def record_of(d: dict, symbols: list, coords, geo: Path, stripped: dict | None = None) -> dict:
    """record.json of a legacy dot from its geometry (charges by the builder recipe)."""
    material = d.get("material") or d["id"].split("-")[0]
    comp = Counter(symbols)
    unknown = sorted(set(comp) - set(CHARGES))
    if unknown:
        raise ValueError(f"{d['id']}: no charge for {', '.join(unknown)}")
    core = {el: n for el, n in comp.items() if el not in LIGANDS}
    ligands = {el: n for el, n in comp.items() if el in LIGANDS}
    parts = d["id"].split("-")
    surface = d.get("surface") or "clean"
    rec = {"schema_version": 1, "id": d["id"], "material": material, "family": FAMILY.get(material),
           "phase": PHASE.get(FAMILY.get(material), ""), "surface": surface, "source": "legacy",
           "formula": "".join(f"{el}{n}" for el, n in comp.items()), "composition": dict(comp),
           "core": core, "ligands": ligands, "n_atoms": len(symbols),
           "total_charge": int(sum(CHARGES[s] for s in symbols)),
           "centre": parts[1] if len(parts) > 1 else None, "centre_offset_A": 0.0,
           "size": _size(symbols, coords, set(core)),
           "stages": [{"stage": "start", "file": "start.xyz"}],
           "origin": {"generator": "qdw campaign legacy-records", "geometry": str(geo),
                      "recipe": {"charges": {el: CHARGES[el] for el in comp},
                                 "passivation": {"ligand": "/".join(ligands) or None}}}}
    if d.get("reconstruction"):
        rec["reconstruction"] = d["reconstruction"]
    if stripped:
        rec["origin"]["stripped_atoms"] = stripped                      # solvent/ligand atoms of the DFT run, not in the id
    return rec


def make_records(campaign_dir, out_dir, searches=(), *, max_atoms: int = 2000, update_ledger: bool = False,
                 strip_extra: bool = False) -> dict:
    """
    Write records for the legacy dots of the ledger that have none. Returns
    {"made": [ids], "missing": {id: reason}, "skipped": {id: reason}}; a dot is missing when no geometry exists
    or its atoms do not match the formula in its id (the geometry may belong to another dot), skipped when it has
    more than `max_atoms` atoms or already has a record. With `strip_extra`, atoms of elements that are not in the
    id's formula (the organic ligands or solvent of the earlier DFT runs) are removed first; what is left must match.
    """
    led = campaign.load(campaign_dir)
    out = Path(out_dir)
    res: dict = {"made": [], "missing": {}, "skipped": {}}
    done = {}
    for d in led["dots"]:
        if d.get("source") != "legacy":
            continue
        if d.get("record") and Path(d["record"]).is_dir():
            res["skipped"][d["id"]] = "has a record"
            continue
        geo = find_geometry(d, list(searches))
        if geo is None:
            res["missing"][d["id"]] = "no geometry found"
            continue
        symbols, coords, _ = read_frames(geo)[-1]                    # the last frame is the final geometry
        want = formula_of(d["id"])
        stripped = None
        if strip_extra:
            keep = [i for i, el in enumerate(symbols) if el in want]
            if len(keep) < len(symbols):
                stripped = dict(Counter(el for el in symbols if el not in want))
                symbols, coords = [symbols[i] for i in keep], coords[keep]
        got = Counter(symbols)
        if got != want:
            res["missing"][d["id"]] = f"{geo.name} has {sum(got.values())} atoms, id says {sum(want.values())}"
            continue
        if len(symbols) > max_atoms:
            res["skipped"][d["id"]] = f"{len(symbols)} atoms above {max_atoms}"
            continue
        rd = out / d["id"]
        rd.mkdir(parents=True, exist_ok=True)
        write_xyz(rd / "start.xyz", symbols, coords, d["id"])
        (rd / "record.json").write_text(json.dumps(record_of(d, symbols, coords, geo, stripped), indent=1))
        done[d["id"]] = str(rd)
        res["made"].append(d["id"])
    if update_ledger and done:
        for d in led["dots"]:
            if d["id"] in done:
                d["record"] = done[d["id"]]
        campaign.ledger_path(campaign_dir).write_text(yaml.safe_dump(led, sort_keys=False, allow_unicode=True))
    return res
