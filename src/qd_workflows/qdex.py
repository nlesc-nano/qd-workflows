"""QDEX on the TREXIO orbitals of a CP2K chain: config from the record and the MO window, and a summary."""
from __future__ import annotations

import csv
import json
import re
from pathlib import Path

import yaml

MAX_STATES = 1000       # occupied and virtual MOs in the diagonal sBSE, at most
BRIGHT_F = 0.05         # oscillator strength of the first "bright" exciton
ORGANIC = {"H", "C"}


def _roles(record: dict) -> tuple[list, list, list]:
    """Core cations, core anions and ligand binding elements of a library record."""
    core = list(record.get("core") or {})
    charges = ((record.get("origin") or {}).get("recipe") or {}).get("charges") or {}
    cations = [e for e in core if charges.get(e, 0) > 0]
    anions = [e for e in core if charges.get(e, 0) < 0]
    if not cations or not anions:          # no recipe charges: first element(s) cations, last anion
        cations, anions = core[:-1], core[-1:]
    ligands = [e for e in (record.get("composition") or {}) if e not in core and e not in ORGANIC]
    return cations, anions, ligands


def coop_pairs(record: dict) -> list:
    """Cation-anion, like-element and cation-ligand pairs (the ligand bond is part of the COOP)."""
    cations, anions, ligands = _roles(record)
    pairs = [f"{c}-{a}" for c in cations for a in anions]
    pairs += [f"{e}-{e}" for e in cations + anions]
    pairs += [f"{c}-{l}" for c in cations for l in ligands]
    return pairs


def find_cif(name: str, dirs) -> str:
    for d in dirs:
        if (Path(d) / name).is_file():
            return str(Path(d) / name)
    raise FileNotFoundError(f"{name} not in {', '.join(map(str, dirs))}")


def prepare(record_path: str, window_path: str, out_dir: str, cluster: dict, *, orbitals: str = "../orbitals.h5",
            geom: str = "../relaxed.xyz", threads: int = 8, max_states: int | None = None) -> dict:
    """
    Write `<out_dir>/config.yaml` with Ivan's QDEX settings (bulk QP with the scaled
    vertex, diagonal sBSE with the Resta kernel, SOC, toluene eps_out, fuzzy bands,
    PDOS and COOP). The sBSE space is balanced: min(occupied, virtual) MOs of the
    TREXIO window, at most `max_states`.
    """
    q, cp = cluster.get("qdex", {}), cluster["cp2k"]
    record = json.loads(Path(record_path).read_text())
    window = json.loads(Path(window_path).read_text())
    n = min(window["occupied_in_window"], window["virtual_in_window"], int(max_states or q.get("max_states", MAX_STATES)))
    cif_name = (record.get("origin") or {}).get("cif")
    cif_dirs = q.get("cif_dirs") or cluster.get("props", {}).get("cif_dirs", [])
    elements = list(record.get("composition") or {})
    config = {
        "system": {"mo_file": orbitals, "xyz": geom, "basis_txt": f"{cp['data_dir']}/{cp['basis_file']}",
                   "basis_name": cp["basis"], "material": str(record["material"]).upper(), "nthreads": int(threads)},
        "environment": {"eps_out": float(q.get("eps_out", 2.24))},
        "quasiparticles": {"model": "bulk", "bulk_vertex": "scaled"},
        "excitations": {"mode": "diagonal_sbse", "kernel": "resta", "nhomos": n, "nlumos": n},
        "fuzzy": {"run_fuzzy": True, "cif": find_cif(cif_name, cif_dirs) if cif_name else None, "soc_window": 10.0,
                  "pdos_atoms": elements, "coop_pairs": coop_pairs(record), "fuzzy_sigma": 0.01,
                  "pdos_sigma": 0.08, "ewin": [-5.0, 5.0]},
        "soc": {"soc_flag": True, "gth_file": q["gth_soc_file"]},
        "output": {"cube": False, "sigma": 0.03, "plot": True, "write_csv": True, "csv_roots": 200},
    }
    if not config["fuzzy"]["cif"]:
        config["fuzzy"]["run_fuzzy"] = False
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    (out / "config.yaml").write_text(yaml.safe_dump(config, sort_keys=False))
    return {"config": str(out / "config.yaml"), "material": config["system"]["material"], "states": n,
            "coop_pairs": config["fuzzy"]["coop_pairs"]}


def _excitons(path: Path, bright_f: float) -> dict:
    rows = list(csv.DictReader(path.open())) if path.is_file() else []
    if not rows:
        return {}
    first = rows[0]
    bright = next((r for r in rows if float(r["f_osc"]) >= bright_f), None)
    out = {"lowest_ev": float(first["Energy_eV"]), "lowest_f": float(first["f_osc"])}
    if bright:
        out.update(first_bright_ev=float(bright["Energy_eV"]), first_bright_f=float(bright["f_osc"]),
                   first_bright_state=int(bright["State"]))
    return out


def _number(pattern: str, text: str):
    m = re.search(pattern, text)
    return float(m.group(1)) if m else None


def summary(qdex_dir: str, window_path: str, bright_f: float = BRIGHT_F) -> dict:
    """
    Frontier levels, gaps and excitons for the webapp, spin-free and with SOC.
    DFT levels are CP2K's (window.json); QDEX prints the SOC spinor edges relative
    to the spin-free mid-gap, and its SOC QP gap is the spin-free QP gap minus the
    SOC shrinking of the DFT gap, so the SOC QP edges are the spin-free QP edges
    moved by the SOC shifts of the DFT edges.
    """
    d = Path(qdex_dir)
    text = (d / "qdex.out").read_text()
    w = json.loads(Path(window_path).read_text())
    mid = 0.5 * (w["e_homo_ev"] + w["e_lumo_ev"])
    sf = {"dft_homo_ev": w["e_homo_ev"], "dft_lumo_ev": w["e_lumo_ev"], "dft_gap_ev": w["gap_ev"],
          "qp_homo_ev": _number(r"QP HOMO \(IP\)\s*:\s*(-?[\d.]+)", text),
          "qp_lumo_ev": _number(r"QP LUMO \(EA\)\s*:\s*(-?[\d.]+)", text),
          "qp_gap_ev": _number(r"\[QP\]\s+Target Gap\s*:\s*(-?[\d.]+)", text)}
    sf.update(_excitons(d / "exciton_results_sf.csv", bright_f))
    soc = {}
    h_rel = _number(r"Spinor HOMO \(Idx \d+\):\s*(-?[\d.]+)", text)
    l_rel = _number(r"Spinor LUMO \(Idx \d+\):\s*(-?[\d.]+)", text)
    if h_rel is not None and l_rel is not None:
        soc = {"dft_homo_ev": mid + h_rel, "dft_lumo_ev": mid + l_rel, "dft_gap_ev": l_rel - h_rel}
        if sf["qp_homo_ev"] is not None:
            soc["qp_homo_ev"] = sf["qp_homo_ev"] + soc["dft_homo_ev"] - sf["dft_homo_ev"]
            soc["qp_lumo_ev"] = sf["qp_lumo_ev"] + soc["dft_lumo_ev"] - sf["dft_lumo_ev"]
        soc["qp_gap_ev"] = _number(r"Final SOC QP Gap\s*:\s*(-?[\d.]+)", text)
    soc.update(_excitons(d / "exciton_results_soc.csv", bright_f))
    sf, soc = ({k: (round(v, 6) if isinstance(v, float) else v) for k, v in x.items()} for x in (sf, soc))
    return {"bright_f": bright_f, "spin_free": sf, "soc": soc,
            "files": sorted(p.name for p in d.iterdir() if p.suffix in (".csv", ".html"))}
