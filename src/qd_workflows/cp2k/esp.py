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
    dot_add<K> the fully added dot, the last state of the MACE adsorption ladder
               (props/adsorb_<K>.xyz), from scratch: the added states (ladder up) are
               interpolated in k between dot_0 and it (the webapp's solution-phase ladder)
    MA, MX     the monomers of props/monomers.xyz, from scratch in a small cell

    qdw cp2k esp prepare <field/base/cp2k.inp> <esp dir> [--stripped detach_K.xyz] [--added adsorb_K.xyz]
                         [--monomers monomers.xyz]
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


RESTART_EACH = 50      # SCF cycles between .wfn writes of the field and dot ESP runs


def restart_print(t: str) -> str:
    """Write the run's own .wfn every RESTART_EACH SCF cycles (no backups), so a run stopped at MAX_SCF
    can go on from where it stopped (scf_restart) instead of from the start."""
    block = (r"\g<1>&RESTART\n\g<1>  BACKUP_COPIES 0\n\g<1>  &EACH\n"
             rf"\g<1>    QS_SCF {RESTART_EACH}\n\g<1>  &END EACH")
    t, n = re.subn(r"^(\s*)&RESTART OFF\s*$", block, t, count=1, flags=re.M)
    if n:
        return t
    return _sub(t, r"^(\s*)&SCF\s*$", r"\g<0>\n\g<1>  &PRINT\n" + block.replace(r"\g<1>", r"\g<1>    ")
                + r"\n\g<1>    &END RESTART\n\g<1>  &END PRINT")


def scf_restart(run_dir) -> str:
    """A run that stopped with "SCF run NOT converged": restart it once from its own .wfn (SCF_GUESS restart);
    its output is kept as cp2k.out.unconverged. Returns "restarted" or why not."""
    d = Path(run_dir)
    out = d / "cp2k.out"
    if not out.is_file() or "SCF run NOT converged" not in out.read_text(errors="ignore"):
        return "not an SCF failure"
    if (d / "cp2k.out.unconverged").exists():
        return "already restarted once"
    t = (d / "cp2k.inp").read_text()
    project = re.search(r"^\s*PROJECT\s+(\S+)", t, flags=re.M).group(1)
    wfn = f"{project}-RESTART.wfn"
    if not (d / wfn).is_file():
        return f"no {wfn}"
    t = _sub(t, r"^(\s*SCF_GUESS\s+)\S+", r"\g<1>restart")
    t, n = re.subn(r"^(\s*WFN_RESTART_FILE_NAME\s+).*$", rf"\g<1>{wfn}", t, count=1, flags=re.M)
    if not n:
        t = _sub(t, r"^(\s*)&DFT\s*$", rf"\g<0>\n\g<1>  WFN_RESTART_FILE_NAME {wfn}")
    (d / "cp2k.inp").write_text(t)
    out.rename(d / "cp2k.out.unconverged")
    if (d / "srun.log").is_file():
        (d / "srun.log").rename(d / "srun.log.unconverged")
    return "restarted"


def esp_input(base: str, project: str, *, scratch: bool = False, cell: str | None = None,
              restartable: bool = False) -> str:
    """The chain's single point as an energy run printing the potential and density cubes;
    scratch: atomic guess instead of the geo_opt .wfn (another structure); cell: its own box;
    restartable: write its own .wfn (restart_print; the dots, not the monomers of the library)."""
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
    return restart_print(t) if restartable else t


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


def prepare(base_inp: str, out_dir: str, stripped: str | None = None, monomers: str | None = None,
            added: str | None = None) -> list[str]:
    """<out_dir>/<run>/{cp2k.inp, geom.xyz} and runs.txt; returns the run names."""
    base_path = Path(base_inp)
    base = base_path.read_text()
    out = Path(out_dir)
    names = []
    sym, pts = read_xyz(base_path.parent / "geom.xyz")
    _write(out / "dot_0", esp_input(base, "dot_0", restartable=True), sym, pts, "intact dot (geo_opt)")
    names.append("dot_0")
    if stripped:
        m = re.search(r"detach_(\d+)\.xyz$", str(stripped))
        k = int(m.group(1)) if m else 0
        sym, pts = read_xyz(stripped)
        name = f"dot_{k}"
        _write(out / name, esp_input(base, name, scratch=True, cell=cell_abc(pts), restartable=True), sym, pts,
               f"fully stripped dot, {k} units removed (MACE ladder)")
        names.append(name)
    if added:
        m = re.search(r"adsorb_(\d+)\.xyz$", str(added))
        k = int(m.group(1)) if m else 0
        sym, pts = read_xyz(added)
        name = f"dot_add{k}"
        _write(out / name, esp_input(base, name, scratch=True, cell=cell_abc(pts), restartable=True), sym, pts,
               f"fully added dot, {k} units added (MACE adsorption ladder)")
        names.append(name)
    if monomers:
        for label, sym, pts in _frames(Path(monomers)):
            cell = cell_abc(pts, min_length=MONOMER_MIN_CELL_A)        # 6 Å each side, common grid
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
        m = re.match(r"dot_add(\d+)$", name)
        if m:                                                # not "k": Orchestr.AI reads the removals only
            res[name]["k_added"] = int(m.group(1))
    return {"model": "COSMO on the DFT electrostatic potential, cavity on the density isosurface "
                     f"{sv.RHO_ISO:.2e} e/bohr^3 (qd_workflows.reaction.solvation)",
            "rho0_e_bohr3": sv.RHO_ISO, "beta0_alpha_gamma_mN_m": [sv.BETA0_SLOPE, sv.BETA0_OFFSET],
            "nonel_min_eps": sv.NONEL_MIN_EPS, "species": res}


def write(out_dir: str, path: str) -> dict:
    data = collect(out_dir)
    Path(path).write_text(json.dumps(data, indent=1))
    return data
