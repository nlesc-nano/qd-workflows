"""
Solvation of a library dot from its DFT chain, with the model of qdw reaction (reaction/solvation.py):

    dG_solv(eps) = (1 - 1/eps) E_inf + (alpha + gamma)(eps) S

E_inf: conductor (COSMO) energy of the DFT electrostatic potential on the electron-density
isosurface rho_0 = e^-7.6 = 5e-4 e/bohr^3; S: the area of that isosurface (Hille et al. 2019,
beta = 0 fit for the cavity and dispersion term, applied by Orchestr.AI's solvation step).

Runs, all with the chain's own single-point settings (the field step's base input: GAPW,
MULTIPOLE, OT/DIIS), printing the Hartree-potential and density cubes:
    dot_0      the intact dot, restarted from the geo_opt .wfn (converges in a few SCF steps)
    dot_<K>    the fully stripped dot, the last state of the MACE desorption ladder
               (props/detach_<K>.xyz), from scratch: the stripped states in between are
               interpolated in k by Orchestr.AI
    MA, MX     the monomers of props/monomers.xyz, from scratch in a small cell

    qdw cp2k esp prepare <field/base/cp2k.inp> <esp dir> [--stripped detach_K.xyz] [--monomers monomers.xyz]
    qdw cp2k esp status <esp dir>          # "all done" or the runs still to do
    qdw cp2k esp collect <esp dir> -o solvation.json
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np

from ..xyz import read_xyz
from .inputs import cell_abc

CUBES = """      &V_HARTREE_CUBE
        STRIDE 2 2 2
        FILENAME =v_hartree
      &END V_HARTREE_CUBE
      &E_DENSITY_CUBE
        STRIDE 2 2 2
        FILENAME =density
      &END E_DENSITY_CUBE
"""
MONOMER_MIN_CELL_A = 15.0
MAX_POINTS = 20000          # surface elements: the COSMO matrix is dense (20,000: 3.2 GB, a minute)


def _sub(text: str, pattern: str, repl: str) -> str:
    new, n = re.subn(pattern, repl, text, count=1, flags=re.M)
    if n == 0:
        raise ValueError(f"base input has no match for {pattern!r}")
    return new


def esp_input(base: str, project: str, *, scratch: bool = False, cell: str | None = None) -> str:
    """The chain's single point as an energy run printing the potential and density cubes;
    scratch: atomic guess instead of the geo_opt .wfn (another structure); cell: its own box."""
    t = _sub(base, r"^(\s*PROJECT\s+)\S+", rf"\g<1>{project}")
    t = _sub(t, r"^(\s*RUN_TYPE\s+)\S+", r"\g<1>ENERGY")
    if re.search(r"^    &PRINT\s*$", t, flags=re.M):          # a &DFT &PRINT already there
        t = _sub(t, r"^(    &PRINT\s*\n)", r"\g<1>" + CUBES.replace("\\", "\\\\"))
    else:
        t = _sub(t, r"^(    &SCF\s*)$", "    &PRINT\n" + CUBES + "    &END PRINT\n" + r"\g<1>")
    if scratch:
        t = _sub(t, r"^(\s*SCF_GUESS\s+)\S+", r"\g<1>atomic")
        t = re.sub(r"^\s*WFN_RESTART_FILE_NAME.*\n", "", t, flags=re.M)
    if cell:
        t = _sub(t, r"^(\s*ABC\s+).*$", rf"\g<1>{cell}")
    return t


def _frames(path: Path) -> list:
    lines, out, i = Path(path).read_text().splitlines(), [], 0
    while i < len(lines) and lines[i].strip():
        n = int(lines[i])
        rows = [l.split() for l in lines[i + 2:i + 2 + n]]
        out.append((lines[i + 1].split()[0], [r[0] for r in rows], np.array([[float(x) for x in r[1:4]] for r in rows])))
        i += n + 2
    return out


def _write(d: Path, inp: str, sym, pts, comment: str) -> None:
    d.mkdir(parents=True, exist_ok=True)
    (d / "cp2k.inp").write_text(inp)
    (d / "geom.xyz").write_text(f"{len(sym)}\n{comment}\n" + "".join(
        f"{s} {x:.8f} {y:.8f} {z:.8f}\n" for s, (x, y, z) in zip(sym, np.asarray(pts, float))))


def prepare(base_inp: str, out_dir: str, stripped: str | None = None, monomers: str | None = None) -> list[str]:
    """<out_dir>/<run>/{cp2k.inp, geom.xyz} and runs.txt; returns the run names."""
    base_path = Path(base_inp)
    base = base_path.read_text()
    out = Path(out_dir)
    names = []
    sym, pts = read_xyz(base_path.parent / "geom.xyz")
    _write(out / "dot_0", esp_input(base, "dot_0"), sym, pts, "intact dot (geo_opt)")
    names.append("dot_0")
    if stripped:
        m = re.search(r"detach_(\d+)\.xyz$", str(stripped))
        k = int(m.group(1)) if m else 0
        sym, pts = read_xyz(stripped)
        name = f"dot_{k}"
        _write(out / name, esp_input(base, name, scratch=True, cell=cell_abc(pts)), sym, pts,
               f"fully stripped dot, {k} units removed (MACE ladder)")
        names.append(name)
    if monomers:
        for label, sym, pts in _frames(Path(monomers)):
            ext = np.ptp(pts, axis=0) + 12.0
            cell = " ".join(f"{max(x, MONOMER_MIN_CELL_A):.3f}" for x in ext)
            _write(out / label, esp_input(base, label, scratch=True, cell=cell), sym, pts, f"{label} monomer (MACE)")
            names.append(label)
    (out / "runs.txt").write_text("\n".join(names) + "\n")
    return names


def _done(d: Path) -> bool:
    o = d / "cp2k.out"
    if not o.is_file():
        return False
    t = o.read_text(errors="ignore")
    return ("SCF run converged" in t and "SCF run NOT converged" not in t
            and any(d.glob("v_hartree*")) and any(d.glob("density*")))


def status(out_dir: str) -> str:
    out = Path(out_dir)
    names = (out / "runs.txt").read_text().split()
    todo = [n for n in names if not _done(out / n)]
    return " ".join(todo) if todo else "all done"


def collect(out_dir: str, max_points: int = MAX_POINTS) -> dict:
    """COSMO energy (1/eps -> 0), isosurface area and volume of every run."""
    from ..reaction import solvation as sv
    out = Path(out_dir)
    res = {}
    for name in (out / "runs.txt").read_text().split():
        d = out / name
        if not _done(d):
            raise RuntimeError(f"{d}: run not finished (cubes and a converged SCF needed)")
        pot, den = sv.read_cube(sorted(d.glob("v_hartree*"))[0]), sv.read_cube(sorted(d.glob("density*"))[0])
        spacing = sv.ISO_SPACING_A
        n = len(sv.isosurface(den, sv.RHO_ISO, spacing)[0])
        if n > max_points:                                   # coarser surface elements for large dots
            spacing *= float(np.sqrt(n / max_points))
        r = sv.cosmo_isodensity(pot, den, sv.RHO_ISO, spacing)
        sym, _ = read_xyz(d / "geom.xyz")
        res[name] = {"E_inf_eV": r["E_inf_eV"], "S_A2": r["S_A2"], "V_A3": r["V_A3"], "n_points": r["n_points"],
                     "spacing_A": spacing, "n_atoms": len(sym), "screening_charge": r["screening_charge"]}
        m = re.match(r"dot_(\d+)$", name)
        if m:
            res[name]["k"] = int(m.group(1))
    return {"model": "COSMO on the DFT electrostatic potential, cavity on the density isosurface "
                     f"{sv.RHO_ISO:.2e} e/bohr^3 (qd_workflows.reaction.solvation)",
            "rho0_e_bohr3": sv.RHO_ISO, "beta0_alpha_gamma_mN_m": [sv.BETA0_SLOPE, sv.BETA0_OFFSET],
            "nonel_min_eps": sv.NONEL_MIN_EPS, "species": res}


def write(out_dir: str, path: str) -> dict:
    data = collect(out_dir)
    Path(path).write_text(json.dumps(data, indent=1))
    return data
