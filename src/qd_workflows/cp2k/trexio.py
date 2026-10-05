"""Trim a TREXIO HDF5 file to the MO window, so QDEX reads (and the archive keeps) only those MOs."""
from __future__ import annotations

from pathlib import Path

import numpy as np

HARTREE_EV = 27.211386245988
TOL_EV = 1e-3


def trim(src, dst, window: dict) -> dict:
    """
    Copy `src` to `dst` keeping the MOs whose energy lies in the window
    [e_min_ev, e_max_ev]. Every other group (basis, AOs, nuclei, cell, metadata)
    is copied unchanged; QDEX reads any MO subset with its energies and occupations.
    """
    import h5py

    src, dst = Path(src), Path(dst)
    if src.resolve() == dst.resolve():
        raise ValueError("trim writes a new file; give a different destination")
    with h5py.File(src, "r") as f, h5py.File(dst, "w") as g:
        for key, value in f.attrs.items():
            g.attrs[key] = value
        for name in f:
            if name != "mo":
                f.copy(f[name], g, name=name)
        mo = f["mo"]
        eps_ev = np.asarray(mo["mo_energy"][:]) * HARTREE_EV
        n_mo = len(eps_ev)
        keep = np.flatnonzero((eps_ev >= window["e_min_ev"] - TOL_EV) & (eps_ev <= window["e_max_ev"] + TOL_EV))
        if keep.size == 0:
            raise ValueError(f"{src}: no MO inside the window")
        out = g.create_group("mo")
        for key, value in mo.attrs.items():
            out.attrs[key] = value
        out.attrs["mo_num"] = np.int64(keep.size)
        for name, ds in mo.items():
            if ds.shape and ds.shape[0] == n_mo:
                out.create_dataset(name, data=ds[keep, ...] if ds.ndim > 1 else ds[:][keep], compression="gzip")
            else:
                out.create_dataset(name, data=ds[()])
        occ = np.asarray(mo["mo_occupation"][:])[keep]
    return {"kept": int(keep.size), "of": int(n_mo), "occupied": int((occ > 1e-3).sum()),
            "virtual": int((occ <= 1e-3).sum()), "bytes": dst.stat().st_size}
