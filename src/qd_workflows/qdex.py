"""QDEX on the TREXIO orbitals of a CP2K chain: config from the record and the MO window, and a summary."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import yaml

MAX_STATES = 1000       # occupied and virtual MOs in the diagonal sBSE, at most
BRIGHT_F = 0.05         # oscillator strength of the first "bright" exciton
EXCITATIONS_EMAX = 4.5  # eV, excitons stored (about 275 nm; ligands absorb above ~4 eV)
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
        # database output: two HDF5 files and four coarse MO cubes, no dashboards / CSV / plots
        "output": {"sigma": 0.03, "h5": True, "html": False, "plot": False, "write_csv": False,
                   "cube": False, "mo_cubes": True, "mo_cube_spacing": float(q.get("mo_cube_spacing", 0.8)),
                   "cube_nhomos": 2, "cube_nlumos": 2,
                   "excitations_emax": float(q.get("excitations_emax", EXCITATIONS_EMAX)),
                   "verbosity": "quiet"},
    }
    if not config["fuzzy"]["cif"]:
        config["fuzzy"]["run_fuzzy"] = False
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    (out / "config.yaml").write_text(yaml.safe_dump(config, sort_keys=False))
    return {"config": str(out / "config.yaml"), "material": config["system"]["material"], "states": n,
            "coop_pairs": config["fuzzy"]["coop_pairs"]}


def _excitons(h5, group: str, bright_f: float) -> dict:
    if group not in h5:
        return {}
    E, f = h5[f"{group}/energy_ev"][:], h5[f"{group}/f_osc"][:]
    if not len(E):
        return {}
    out = {"lowest_ev": round(float(E[0]), 6), "lowest_f": float(f[0]),
           "n_states_stored": int(len(E)), "emax_ev": float(h5[group].attrs["excitations_emax_ev"])}
    bright = np.flatnonzero(f >= bright_f)
    if bright.size:
        k = int(bright[0])
        out.update(first_bright_ev=round(float(E[k]), 6), first_bright_f=float(f[k]), first_bright_state=k + 1)
    return out


def summary(qdex_dir: str, bright_f: float = BRIGHT_F) -> dict:
    """
    What the webapp shows, from QDEX's qdex_electronic.h5 and qdex_excitations.h5:
    DFT and QP HOMO, LUMO and gap (vs vacuum), lowest and first bright exciton,
    spin-free and with SOC, and the files to publish.
    """
    import h5py
    d = Path(qdex_dir)
    with h5py.File(d / "qdex_electronic.h5", "r") as el, h5py.File(d / "qdex_excitations.h5", "r") as ex:
        qp = {k: (float(v) if isinstance(v, (int, float, np.number)) else v) for k, v in el["qp"].attrs.items()}
        h = int(el["sf/mo"].attrs["homo_index"])
        e_dft = el["sf/mo/energy_dft_abs_ev"][:]
        sf = {"dft_homo_ev": float(e_dft[h]), "dft_lumo_ev": float(e_dft[h + 1]), "dft_gap_ev": qp["dft_gap_ev"],
              "qp_homo_ev": qp["qp_homo_ev"], "qp_lumo_ev": qp["qp_lumo_ev"], "qp_gap_ev": qp["qp_gap_ev"]}
        sf.update(_excitons(ex, "sf", bright_f))
        soc = {}
        if "soc_dft_homo_ev" in qp:
            soc = {"dft_homo_ev": qp["soc_dft_homo_ev"], "dft_lumo_ev": qp["soc_dft_lumo_ev"],
                   "dft_gap_ev": qp["soc_dft_gap_ev"], "qp_homo_ev": qp["soc_qp_homo_ev"],
                   "qp_lumo_ev": qp["soc_qp_lumo_ev"], "qp_gap_ev": qp["soc_qp_gap_ev"]}
            soc.update(_excitons(ex, "soc", bright_f))
        commit = str(el.attrs.get("commit", ""))
    rnd = lambda x: {k: (round(v, 6) if k.endswith("_ev") and isinstance(v, float) else v) for k, v in x.items()}
    return {"qdex_commit": commit, "bright_f": bright_f, "spin_free": rnd(sf), "soc": rnd(soc),
            "files": sorted(p.name for p in d.iterdir() if p.suffix in (".h5", ".cube") and p.name.startswith(("qdex_", "spatial_")))}
