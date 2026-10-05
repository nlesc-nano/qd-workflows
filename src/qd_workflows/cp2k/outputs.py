"""Reading CP2K outputs: MO window from .pdos files, energies, forces, relaxed geometry."""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np

from ..xyz import read_frames

HARTREE_EV = 27.211386245988
HARTREE_BOHR_TO_EV_A = 51.422086190832
WINDOW_EV = 5.0


def read_pdos(path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """MO index (1-based), eigenvalue (eV) and occupation from a CP2K .pdos file."""
    rows = [line.split()[:3] for line in Path(path).read_text().splitlines() if line.strip() and not line.startswith("#")]
    data = np.array(rows, dtype=float)
    return data[:, 0].astype(int), data[:, 1] * HARTREE_EV, data[:, 2]


def mo_window(pdos_files, below_ev: float = WINDOW_EV, above_ev: float = WINDOW_EV) -> dict:
    """
    MOs from `below_ev` under the valence band edge (HOMO) to `above_ev` over the
    conduction band edge (LUMO). Every .pdos file of one run lists the same MOs;
    the first is read. `added_mos` is the number of virtual MOs in the window,
    which the TREXIO step asks CP2K for.
    """
    files = sorted(str(f) for f in pdos_files)
    if not files:
        raise ValueError("no .pdos file given")
    idx, eps, occ = read_pdos(files[0])
    occupied = occ > 1e-3
    if not occupied.any() or occupied.all():
        raise ValueError(f"{files[0]}: needs occupied and unoccupied MOs (print the DOS with NLUMO -1)")
    homo = int(idx[occupied].max())
    lumo = int(idx[~occupied].min())
    e_homo, e_lumo = float(eps[idx == homo][0]), float(eps[idx == lumo][0])
    lo_e, hi_e = e_homo - below_ev, e_lumo + above_ev
    inside = (eps >= lo_e) & (eps <= hi_e)
    first, last = int(idx[inside].min()), int(idx[inside].max())
    if eps.max() < hi_e:
        raise ValueError(f"{files[0]}: the highest MO ({eps.max():.2f} eV) is below the window top "
                         f"({hi_e:.2f} eV); the DOS print lists too few unoccupied MOs")
    return {
        "source": files[0], "below_ev": below_ev, "above_ev": above_ev,
        "homo_index": homo, "lumo_index": lumo,
        "e_homo_ev": round(e_homo, 6), "e_lumo_ev": round(e_lumo, 6), "gap_ev": round(e_lumo - e_homo, 6),
        "e_min_ev": round(lo_e, 6), "e_max_ev": round(hi_e, 6),
        "first_index": first, "last_index": last,
        "occupied_in_window": homo - first + 1, "virtual_in_window": last - homo,
        "added_mos": last - homo, "n_mo_listed": int(len(idx)),
    }


_ENERGY = re.compile(r"ENERGY\|\s*Total FORCE_EVAL \( QS \) energy\s*\[(?:a\.u\.|hartree)\]:\s*(-?\d+\.\d+(?:[Ee][-+]?\d+)?)")


def read_energy(cp2k_out) -> float:
    """Last total energy in a CP2K output, in eV."""
    hits = _ENERGY.findall(Path(cp2k_out).read_text())
    if not hits:
        raise ValueError(f"no total energy in {cp2k_out}")
    return float(hits[-1]) * HARTREE_EV


def read_forces(forces_file, n_atoms: int) -> np.ndarray:
    """
    Forces (eV/Å) from a CP2K FORCES print (Hartree/bohr). Refuses a file whose
    atom count differs from the structure's: a sample without every force is not
    a training label.
    """
    rows = []
    for line in Path(forces_file).read_text().splitlines():
        tok = line.split()
        if len(tok) >= 6 and tok[0].isdigit() and tok[1].isdigit():
            rows.append([float(x) for x in tok[3:6]])
    if len(rows) != n_atoms:
        raise ValueError(f"{forces_file}: {len(rows)} forces for {n_atoms} atoms")
    return np.array(rows) * HARTREE_BOHR_TO_EV_A


def last_frame(pos_file):
    """(symbols, coordinates) of the last frame of a geo_opt trajectory (<PROJECT>-pos-1.xyz)."""
    symbols, coords, _ = read_frames(pos_file)[-1]
    return symbols, coords
