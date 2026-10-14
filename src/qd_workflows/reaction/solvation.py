"""
Solvation free energy of one species as a continuous function of the dielectric constant,
so the reaction page can move eps with a slider:

    dG_solv(eps) = dG_el(eps) + dG_cav(eps)

Electrostatics (cond["electrostatics"]):
  cosmo  conductor-like screening of the DFT potential of the gas-phase single point on the
         density isosurface rho_0 = 5e-4 e/bohr^3 (cosmo_isodensity; RESP charges in 1.2 x vdW
         spheres only as the fallback): (1 - 1/eps) E_inf. Post-processing, any size; no
         electronic polarisation of the solute. Default.
  sccs   CP2K SCCS single points (GPW, periodic cell, single-threaded: see cp2k.py and
         job_cp2k.sh) at a few eps, dG_el = E(eps) - E(vacuum), interpolated piecewise-linearly in
         f = 1 - 1/eps through (0, 0); includes the electronic polarisation. Molecules; for the
         InAs dots it converged at eps 2 only.
  gb     Generalized Born, (1 - 1/eps) G_GB^inf, on the RESP charges (the CP2K Hirshfeld
         charges, with Gaussian shape functions, overestimate dipoles 2-30x).

Cavity + repulsion + dispersion (cond["cavity"]), Hille et al., J. Chem. Phys. 150, 041710 (2019):
  beta0  their beta = 0 fit (SI Table S2): (alpha + gamma) S, with (alpha + gamma) fitted over the 67
         non-aqueous solvents as 2.16 ln(eps) - 14.79 dyn/cm (this work, R^2 0.71; data/). A surface
         term only: fine for nanocrystals. Default.
  hille  their main fit: (3.0 ln(eps) + 36.0) S - 0.5 GPa V (Eq. 19). The volume term was fitted on
         molecules; for a 5,700 A^3 dot it gives -18 eV: molecules only.
  none
S and V are those of the same isosurface (within 1-4 % of the SCCS quantum surface, the InAs
dot included); 1.3 x vdW spheres as the fallback. Defined for eps >= 1.8 (pentane 1.84); eps = 1 is the gas phase.

    python -m qd_workflows.reaction.solvation alpb <run dir>     GFN2-xTB ALPB check
"""
from __future__ import annotations

import json
import math
import re
import sys
from pathlib import Path

import numpy as np

HARTREE_EV = 27.211386245988
MN_M_EV_A2 = 1e-3 / 1.602176634e-19 * 1e-20        # 1 mN/m in eV/Å^2
GPA_EV_A3 = 1e9 / 1.602176634e-19 * 1e-30          # 1 GPa in eV/Å^3
HILLE_SLOPE, HILLE_OFFSET, HILLE_BETA_GPA = 3.0, 36.0, -0.5
BETA0 = json.loads((Path(__file__).parent / "data" / "hille2019_beta0.json").read_text())
BETA0_SLOPE, BETA0_OFFSET = BETA0["fit_nonaqueous"]["slope"], BETA0["fit_nonaqueous"]["intercept"]
NONEL_MIN_EPS = 1.8
CAVITY_SCALE = 1.3                                 # vdW scale of the cavity for S and V
RHO_MAX = math.exp(-7.6 + 1.0)                     # 1.36e-3 e/bohr^3
RHO_MIN = math.exp(-7.6 - 1.0)                     # 1.84e-4 e/bohr^3
ALPB_EPS = {"hexane": 1.88, "toluene": 2.38, "chcl3": 4.71, "thf": 7.43, "acetone": 20.7, "dmso": 46.7,
            "water": 78.4}


def alpha_gamma(eps: float, model: str = "beta0") -> float:
    """(alpha + gamma) in mN/m: the beta = 0 fit (default) or Hille's Eq. 19 (beta = -0.5 GPa)."""
    if model == "hille":
        return HILLE_SLOPE * math.log(eps) + HILLE_OFFSET
    return BETA0_SLOPE * math.log(eps) + BETA0_OFFSET


def cavity_term(eps: float, surface_A2: float, volume_A3: float, model: str = "beta0") -> float:
    """Cavity + repulsion + dispersion (eV); 0 in the gas phase or for model none, nan for 1 < eps < 1.8."""
    if eps <= 1.0 or model == "none":
        return 0.0
    if eps < NONEL_MIN_EPS:
        return float("nan")
    e = alpha_gamma(eps, model) * MN_M_EV_A2 * surface_A2
    return e + (HILLE_BETA_GPA * GPA_EV_A3 * volume_A3 if model == "hille" else 0.0)


def nonel(eps: float, surface_A2: float, volume_A3: float) -> float:
    """Hille's main fit (beta = -0.5 GPa); kept for the tests and the comparison."""
    return cavity_term(eps, surface_A2, volume_A3, "hille")


def f_of(eps: float) -> float:
    return 1.0 - 1.0 / eps


def interp_f(knots: list, eps: float) -> float:
    """Piecewise-linear in f = 1 - 1/eps through knots [(f, value), ...] (sorted); linear beyond the ends."""
    f = f_of(eps)
    xs = [k[0] for k in knots]
    ys = [k[1] for k in knots]
    if len(xs) == 1:
        return ys[0]
    j = int(np.clip(np.searchsorted(xs, f) - 1, 0, len(xs) - 2))
    return ys[j] + (ys[j + 1] - ys[j]) * (f - xs[j]) / (xs[j + 1] - xs[j])


def sccs_knots(points: dict) -> dict:
    """From {eps: {dG_el, S, V}}: knots of dG_el (with the vacuum point) and of S and V, in f."""
    eps = sorted(float(e) for e in points)
    el = [(0.0, 0.0)] + [(f_of(e), points[_key(points, e)]["dG_el_eV"]) for e in eps]
    s = [(f_of(e), points[_key(points, e)]["S_A2"]) for e in eps]
    v = [(f_of(e), points[_key(points, e)]["V_A3"]) for e in eps]
    return {"el": el, "S": s, "V": v}


def _key(points, e):
    return next(k for k in points if abs(float(k) - e) < 1e-9)


def electrostatic(sp: dict, eps: float, model: str = "cosmo") -> float:
    """dG_el (eV) of a results.json species at eps by `model` (cosmo | sccs | gb); nan if not available."""
    if eps <= 1.0:
        return 0.0
    if model == "sccs":
        return interp_f(sp["sccs"]["el"], eps) if sp.get("sccs") else float("nan")
    key = "cosmo_inf_eV" if model == "cosmo" else "gb_inf_eV"
    return f_of(eps) * sp[key] if sp.get(key) is not None else float("nan")


def dg_solv(sp: dict, eps: float, electrostatics: str = "cosmo", cavity: str | bool = "beta0") -> float:
    """dG_solv (eV) of one species entry of results.json at eps (1 = gas)."""
    if cavity is True:
        cavity = "beta0"
    elif cavity is False:
        cavity = "none"
    if eps <= 1.0:
        return 0.0
    el = electrostatic(sp, eps, electrostatics)
    if cavity == "none":
        return el
    cav = sp.get("cavity")
    if not cav:
        return float("nan")
    return el + cavity_term(eps, cav["S_A2"], cav["V_A3"], cavity)


# --------------------------------------------------------------------------
# Conductor-like screening on point charges (post-processing of one gas-phase single point)
# --------------------------------------------------------------------------
#
#   cavity     union of spheres R_a = scale * r_vdW(a) (Bondi), points t_k on its surface with areas a_k
#   A_kl       1 / |t_k - t_l|,   A_kk = 1.07 sqrt(4 pi / a_k)               (Klamt, COSMO)
#   Phi_k      sum_a q_a / |t_k - R_a|                                          potential of the solute
#   sigma      = -A^-1 Phi                                                     conductor surface charges
#   dG_el(eps) = (1 - 1/eps) * (1/2) sum_k sigma_k Phi_k                       (scaled conductor energy)
#
# The charges are RESP (fitted to the DFT potential around the molecule): CP2K's default Hirshfeld
# charges (Gaussian shape functions) overestimate dipoles 2-30x and make Generalized Born far too
# large. Fixed charges: no electronic polarisation of the solute (SCCS has it).

VDW = {"H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47, "P": 1.80, "S": 1.80, "Cl": 1.75, "Br": 1.85,
       "I": 1.98, "Se": 1.90, "Te": 2.06, "As": 1.85, "Sb": 2.06, "Zn": 1.39, "Cd": 1.58, "Hg": 1.55,
       "Ga": 1.87, "In": 1.93, "Al": 1.84, "Pb": 2.02, "Cs": 3.43}
COSMO_SCALE = 1.2
COULOMB_EV_A = 14.3996454784


def _sphere_points(n: int) -> np.ndarray:
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    th = np.pi * (1 + 5 ** 0.5) * i
    return np.stack([np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)], axis=1)


def cavity(symbols, positions, scale: float = COSMO_SCALE, density: float = 4.0):
    """Surface points (A), their areas (A^2) and the enclosed volume (A^3) of the scaled-vdW cavity;
    `density` points per A^2 on each sphere."""
    x = np.asarray(positions, float)
    r = np.array([VDW.get(s, 2.0) for s in symbols]) * scale
    pts, area = [], []
    for a in range(len(x)):
        n = max(12, int(4 * np.pi * r[a] ** 2 * density))
        p = x[a] + r[a] * _sphere_points(n)
        d = np.linalg.norm(p[:, None, :] - x[None, :, :], axis=2)
        d[:, a] = np.inf
        keep = (d > r[None, :]).all(axis=1)
        pts.append(p[keep])
        area.append(np.full(keep.sum(), 4 * np.pi * r[a] ** 2 / n))
    pts, area = np.vstack(pts), np.concatenate(area)
    # volume: Monte Carlo on a box around the spheres
    rng = np.random.default_rng(0)
    lo, hi = (x - r[:, None]).min(0), (x + r[:, None]).max(0)
    m = 200000
    s = lo + (hi - lo) * rng.random((m, 3))
    inside = np.zeros(m, bool)
    for a in range(len(x)):
        inside |= ((s - x[a]) ** 2).sum(1) < r[a] ** 2
    vol = float(np.prod(hi - lo) * inside.mean())
    return pts, area, vol


def cosmo_inf(symbols, positions, charges, scale: float = COSMO_SCALE) -> dict:
    """Conductor energy (eV) of the point charges in the scaled-vdW cavity: dG_el = (1 - 1/eps) * E_inf;
    plus the cavity surface S (A^2) and volume V (A^3)."""
    t, a, vol = cavity(symbols, positions, scale)
    q = np.asarray(charges, float)
    x = np.asarray(positions, float)
    phi = (q[None, :] / np.linalg.norm(t[:, None, :] - x[None, :, :], axis=2)).sum(1)
    d = np.linalg.norm(t[:, None, :] - t[None, :, :], axis=2)
    np.fill_diagonal(d, 1.0)
    A = 1.0 / d
    np.fill_diagonal(A, 1.07 * np.sqrt(4 * np.pi / a))
    sigma = -np.linalg.solve(A, phi)
    return {"E_inf_eV": float(0.5 * COULOMB_EV_A * sigma @ phi), "S_A2": float(a.sum()), "V_A3": vol,
            "n_points": int(len(t)), "scale": scale}


# --------------------------------------------------------------------------
# COSMO on the DFT potential, cavity = electron-density isosurface (post-processing of the ESP run)
# --------------------------------------------------------------------------
#
# Spheres around atoms are the wrong boundary for an ionic nanocrystal surface: they leave
# crevices between surface atoms where the potential is strong, and the conductor energy of the
# InAs dot changed 4x between 1.1 and 1.4 x vdW (and RESP point charges, singular at the atoms,
# made it worse). Here the cavity is the density isosurface rho = rho_0 that SCCS uses as the
# centre of its switching (Hille's n_c = 7.6: rho_0 = e^-7.6 = 5.0e-4 e/bohr^3), the potential on
# it is the DFT one (Hartree-potential cube; CP2K's sign is the opposite of the usual), and its
# area is the cavity surface of the cavity/dispersion term.

BOHR_A = 0.529177210903
RHO_ISO = math.exp(-7.6)
ISO_SPACING_A = 0.55        # marching-cubes step: about one surface element per 0.3 A^2


def read_cube(path) -> dict:
    """Gaussian cube (CP2K): atomic numbers, positions (A), grid axes (A), data (as written)."""
    lines = Path(path).read_text().split("\n")
    n = abs(int(lines[2].split()[0]))
    origin = np.array(lines[2].split()[1:4], float) * BOHR_A
    ax = [lines[3 + i].split() for i in range(3)]
    npts = [int(a[0]) for a in ax]
    step = np.array([[float(v) for v in a[1:4]] for a in ax]) * BOHR_A
    if np.abs(step - np.diag(np.diag(step))).max() > 1e-8:
        raise ValueError(f"{path}: only orthogonal cube grids are read")
    atoms = [lines[6 + k].split() for k in range(n)]
    data = np.array(re.findall(r"-?\d\.\d+E[-+]\d+", "\n".join(lines[6 + n:])), float).reshape(npts)
    axes = [origin[i] + step[i, i] * np.arange(npts[i]) for i in range(3)]
    return {"Z": [int(a[0]) for a in atoms], "positions": np.array([[float(v) for v in a[2:5]] for a in atoms]) * BOHR_A,
            "axes": axes, "data": data}


def isosurface(density: dict, rho0: float = RHO_ISO, spacing: float = ISO_SPACING_A):
    """Surface elements of the isosurface rho = rho0: centroids (A), areas (A^2); enclosed volume (A^3)."""
    from skimage.measure import marching_cubes
    h = [ax[1] - ax[0] for ax in density["axes"]]
    stride = max(1, int(round(spacing / min(h))))
    verts, faces, _, _ = marching_cubes(density["data"], level=rho0, spacing=tuple(h), step_size=stride)
    verts = verts + np.array([ax[0] for ax in density["axes"]])
    tri = verts[faces]
    tri_area = 0.5 * np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1)
    # surface elements on the mesh vertices, each with a third of its triangles' area: more uniform
    # than triangle centroids, whose slivers (tiny areas, self-terms ~ 1/sqrt(a)) made the conductor
    # energy of the InAs dot jump by 15 % between grid spacings
    area = np.zeros(len(verts))
    for k in range(3):
        np.add.at(area, faces[:, k], tri_area / 3)
    keep = area > 1e-8
    vol = float((density["data"] > rho0).sum() * np.prod(h))
    return verts[keep], area[keep], vol


def cosmo_isodensity(potential: dict, density: dict, rho0: float = RHO_ISO, spacing: float | None = None) -> dict:
    """Conductor energy (eV) of the DFT potential on the density isosurface: dG_el = (1 - 1/eps) E_inf;
    plus the cavity surface S (A^2) and volume V (A^3). spacing: marching-cubes step (A)."""
    from scipy.interpolate import RegularGridInterpolator
    t, a, vol = isosurface(density, rho0, spacing or ISO_SPACING_A)
    v = RegularGridInterpolator(potential["axes"], potential["data"], bounds_error=False, fill_value=np.nan)(t)
    if np.isnan(v).any():
        raise ValueError("isosurface outside the potential cube")
    phi = -v / BOHR_A                                   # Hartree (CP2K sign) -> e/A, as sum q/r
    d = np.linalg.norm(t[:, None, :] - t[None, :, :], axis=2)
    np.fill_diagonal(d, 1.0)
    A = 1.0 / d
    np.fill_diagonal(A, 1.07 * np.sqrt(4 * np.pi / a))
    sigma = -np.linalg.solve(A, phi)
    return {"E_inf_eV": float(0.5 * COULOMB_EV_A * sigma @ phi), "S_A2": float(a.sum()), "V_A3": vol,
            "n_points": int(len(t)), "rho0": rho0, "screening_charge": float(sigma.sum())}


def esp_cubes(esp_dir) -> tuple[Path, Path] | None:
    """(potential, density) cube files written by the ESP run (CP2K names them without .cube)."""
    d = Path(esp_dir)
    pot = sorted(d.glob("v_hartree*"))
    den = sorted(d.glob("density*"))
    return (pot[0], den[0]) if pot and den else None


def read_resp(path) -> np.ndarray | None:
    """RESP charges written by CP2K (RESP_CHARGES_TO_FILE)."""
    p = Path(path)
    if not p.is_file():
        return None
    q = re.findall(r"^\s*RESP\s+\d+\s+[A-Z][a-z]?\s+(-?\d+\.\d+)", p.read_text(), re.M)
    return np.array(q, float) if q else None


# --------------------------------------------------------------------------
# CP2K outputs
# --------------------------------------------------------------------------

def converged(out_text: str) -> bool:
    """The last SCF of the run converged."""
    m = re.findall(r"SCF run (NOT )?converged", out_text)
    return bool(m) and m[-1] == ""


def read_sccs(out_path) -> dict | None:
    """Total energy (eV) and, for a solvated run, cavity S (Å^2), V (Å^3) and the polarisation energy."""
    p = Path(out_path)
    if not p.is_file():
        return None
    t = p.read_text(errors="ignore")
    # a converged SCF and the total energy printed after it are enough: a job killed while CP2K
    # writes its closing references (no PROGRAM ENDED) still has the final energy
    if not converged(t):
        return None
    e = re.findall(r"ENERGY\| Total FORCE_EVAL \( QS \) energy[^\n-]*?(-?\d+\.\d+)\s*$", t, re.M)
    if not e:
        return None
    res = {"E_eV": float(e[-1]) * HARTREE_EV}
    vol = re.findall(r"SCCS\| Volume of the solute cavity\s+\[bohr\^3\]\s+\S+\s*\n\s*SCCS\|\s+\[angstrom\^3\]\s+(\S+)", t)
    sur = re.findall(r"SCCS\| Surface of the solute cavity \[bohr\^2\]\s+\S+\s*\n\s*SCCS\|\s+\[angstrom\^2\]\s+(\S+)", t)
    pol = re.findall(r"SCCS\| Polarisation energy\s+\[Hartree\]\s+(-?\S+)", t)
    if vol and sur:
        res.update({"V_A3": float(vol[-1]), "S_A2": float(sur[-1])})
    if pol:
        res["E_pol_eV"] = float(pol[-1]) * HARTREE_EV
    return res


def gb_inf(symbols, positions, charges) -> float:
    """-(1/2) sum q_i q_j / f_GB in eV (dG_el = (1 - 1/eps) times this)."""
    from orchestr_ai.qd.steps.solvation import gb_conductor_energy
    return gb_conductor_energy(list(symbols), np.asarray(positions, float), np.asarray(charges, float))


# --------------------------------------------------------------------------
# GFN2-xTB ALPB check
# --------------------------------------------------------------------------

def run_alpb(run_dir: str) -> None:
    """species/<name>/alpb.json: GFN2-xTB ALPB (gsolv: 1 M gas -> 1 M solution) at the DFT geometry."""
    from orchestr_ai.qd.engines import XtbRunner, available_cpus
    from ..xyz import read_xyz
    run = Path(run_dir)
    rx = json.loads((run / "reaction.json").read_text())
    xtb = XtbRunner("gfn2", threads=max(1, available_cpus()))
    for name, sp in rx["species"].items():
        out = run / "species" / name / "alpb.json"
        if out.is_file():
            continue
        sym, pos = read_xyz(run / "species" / name / "dft.xyz")
        kw = dict(charge=sp["charge"], uhf=sp["mult"] - 1)
        res = {"method": "GFN2-xTB ALPB (gsolv)", "eps": ALPB_EPS, "dG_eV": {}}
        try:
            gas = xtb.run(sym, pos, **kw).energy_eV
            for s in ALPB_EPS:
                try:
                    res["dG_eV"][s] = xtb.run(sym, pos, solvation=["--alpb", s, "gsolv"], **kw).energy_eV - gas
                except RuntimeError as exc:
                    res["dG_eV"][s] = None
                    print(f"[alpb] {name} {s}: {exc}".splitlines()[0])
        except RuntimeError as exc:
            res["error"] = str(exc).splitlines()[0]
        out.write_text(json.dumps(res, indent=1))
        print(f"[alpb] {name}: " + ", ".join(f"{s} {v:+.3f}" for s, v in res["dG_eV"].items() if v is not None), flush=True)


if __name__ == "__main__":
    if len(sys.argv) != 3 or sys.argv[1] != "alpb":
        raise SystemExit("usage: python -m qd_workflows.reaction.solvation alpb <run dir>")
    run_alpb(sys.argv[2])
