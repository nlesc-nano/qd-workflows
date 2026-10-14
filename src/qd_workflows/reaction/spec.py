"""The reaction file (YAML): species, equation, units; balance checks; the normalised reaction.json."""
from __future__ import annotations

import re
from collections import Counter
from fractions import Fraction
from pathlib import Path

import numpy as np
import yaml

from ..xyz import read_frames

HARTREE_EV = 27.211386245988
PHASES = ("solute", "solvent", "gas")
DEFAULT_SETTINGS = {
    "head": "omat_pbe",            # MACE-MH-1 head
    "fmax": 0.005,                 # eV/Å, MACE relaxation from the DFT structure
    "max_steps": 3000,
    "hessian": "auto",             # auto | analytic | fd
    "analytic_max_atoms": 2000,
    "fd_step": 0.01,               # Å
    "rmsd_warn": 0.1,              # Å, DFT -> MACE minimum: above this the species is flagged
    "sccs": True,                  # CP2K SCCS (GPW) electrostatics + cavity/dispersion descriptors
    "sccs_eps": [2.0, 4.0, 10.0, 80.0],
    "sp_functional": "PBE",        # GAPW single point at the given geometry: Hirshfeld charges, .wfn (PBE | HLE17)
    "sccs_functional": "PBE",      # GPW vacuum + SCCS runs, restarted from that .wfn (PBE | HLE17 | PBE0)
    "energy_check_functional": None,  # an isolated GPW single point at this level too (e.g. PBE0): dE check
    "hfx_memory_mb": 6000,         # PBE0: MB of stored integrals per rank
    "omp_threads": 1,              # NHR: with OpenMP threads every rank ended up pinned to one core
    "sccs_max_atoms": 2000,
    "alpb": ["hexane", "toluene", "chcl3", "thf", "acetone", "dmso", "water"],   # GFN2-xTB check
}
DEFAULT_T = [250.0, 800.0, 10.0]
TERM = re.compile(r"^\s*(\d+(?:\.\d+)?)?\s*\*?\s*([A-Za-z_][\w\-.]*)\s*$")
XYZ_ENERGY = re.compile(r"\bE\s*=\s*(-?\d+\.\d+(?:[eE][-+]?\d+)?)")


def parse_equation(text: str) -> dict:
    """'a + 12 b -> c + 4 d' as {name: nu}, nu < 0 for reactants."""
    parts = re.split(r"->|=>|=", text)
    if len(parts) != 2:
        raise ValueError(f"equation needs one '->': {text!r}")
    nu: dict = {}
    for side, sign in zip(parts, (-1, 1)):
        for term in side.split("+"):
            m = TERM.match(term)
            if not m:
                raise ValueError(f"cannot read the term {term.strip()!r} of {text!r}")
            coef = float(m.group(1)) if m.group(1) else 1.0
            nu[m.group(2)] = nu.get(m.group(2), 0.0) + sign * coef
    return {k: (int(v) if float(v).is_integer() else v) for k, v in nu.items() if v != 0}


def xyz_energy(path) -> float | None:
    """Energy in the comment line of the first XYZ frame ('E = ...', as CP2K's pos files write it)."""
    m = XYZ_ENERGY.search(read_frames(path)[0][2])
    return float(m.group(1)) if m else None


def balance(species: dict, nu: dict) -> dict:
    """Element and charge residuals sum_i nu_i n_i,Z (all zero for a balanced reaction)."""
    res: Counter = Counter()
    for name, v in nu.items():
        for el, n in species[name]["composition"].items():
            res[el] += v * n
        res["charge"] += v * species[name]["charge"]
    return {k: v for k, v in res.items() if abs(v) > 1e-9}


def suggest_coefficients(species: dict, reactants: list, products: list) -> dict | None:
    """Smallest integer coefficients that balance elements and charge, if there is exactly one way."""
    names = list(reactants) + list(products)
    rows = sorted({el for n in names for el in species[n]["composition"]}) + ["charge"]
    a = np.array([[(-1 if n in reactants else 1) * (species[n]["charge"] if r == "charge"
                   else species[n]["composition"].get(r, 0)) for n in names] for r in rows], float)
    _, sv, vt = np.linalg.svd(a, full_matrices=True)
    rank = int((sv > 1e-9 * sv.max()).sum())
    null = vt[rank:]
    if len(null) != 1:
        return None
    v = null[0] * np.sign(null[0][np.argmax(np.abs(null[0]))])
    if np.any(v <= 1e-9):                    # every species must take part with a positive coefficient
        return None
    fr = [Fraction(x).limit_denominator(1000) for x in v / v.min()]
    lcm = np.lcm.reduce([f.denominator for f in fr])
    return {n: int(f * lcm) for n, f in zip(names, fr)}


def load(path: str) -> dict:
    """Read and check a reaction file; returns the normalised reaction (energies in eV, paths absolute)."""
    src = Path(path).resolve()
    raw = yaml.safe_load(src.read_text())
    unit = str(raw.get("energy_unit", "")).lower()
    if unit not in ("hartree", "ev"):
        raise ValueError("energy_unit must be given: hartree or eV")
    scale = HARTREE_EV if unit == "hartree" else 1.0
    nu = parse_equation(raw["equation"])
    species = {}
    for name, sp in (raw.get("species") or {}).items():
        xyz = (src.parent / sp["xyz"]).resolve()
        symbols, _, _ = read_frames(xyz)[0]
        e = sp.get("energy", "xyz")
        if e == "xyz":
            e = xyz_energy(xyz)
            if e is None:
                raise ValueError(f"{name}: no energy given and none in the comment line of {xyz}")
        phase = sp.get("phase", "solute")
        if phase not in PHASES:
            raise ValueError(f"{name}: phase {phase!r} not supported (one of {', '.join(PHASES)})")
        if phase == "solvent" and "c_liquid" not in sp:
            raise ValueError(f"{name}: a solvent reactant needs c_liquid (mol/L of the pure liquid, water 55.3)")
        entry = {
            "name": name, "xyz": str(xyz), "symbols": symbols, "n_atoms": len(symbols),
            "composition": dict(sorted(Counter(symbols).items())),
            "energy_eV": float(e) * scale, "energy_input": float(e), "charge": int(sp.get("charge", 0)),
            "mult": int(sp.get("mult", 1)), "phase": phase,
            "c": float(sp.get("c", 1.0)), "c_liquid": sp.get("c_liquid"), "p_bar": sp.get("p_bar", 1.0),
            "symmetry_number": sp.get("symmetry_number"), "label": sp.get("label"),
            "charges_file": str((src.parent / sp["charges"]).resolve()) if sp.get("charges") else None,
        }
        if (entry["mult"] - 1) % 2 != (sum(Z(s) for s in symbols) - entry["charge"]) % 2:
            raise ValueError(f"{name}: multiplicity {entry['mult']} does not fit {sum(Z(s) for s in symbols)} "
                             f"electrons at charge {entry['charge']}")
        species[name] = entry
    missing = [n for n in nu if n not in species]
    if missing:
        raise ValueError(f"species in the equation without an entry: {', '.join(missing)}")
    sides = []
    for text in raw.get("equilibria") or []:
        snu = parse_equation(text)
        miss = [n for n in snu if n not in species]
        if miss:
            raise ValueError(f"species in the equilibrium {text!r} without an entry: {', '.join(miss)}")
        if balance(species, snu):
            raise ValueError(f"equilibrium {text!r} not balanced: residual {balance(species, snu)}")
        sides.append({"equation": text, "nu": snu})
    used = list(nu) + [n for s in sides for n in s["nu"] if n not in nu]
    used = list(dict.fromkeys(used))
    unused = [n for n in species if n not in used]
    res = balance(species, nu)
    if res:
        reac = [n for n, v in nu.items() if v < 0]
        prod = [n for n, v in nu.items() if v > 0]
        hint = suggest_coefficients(species, reac, prod)
        msg = f"reaction not balanced: residual {res}"
        if hint:
            msg += ("; balanced: " + " + ".join(_term(hint[n], n) for n in reac) + " -> "
                    + " + ".join(_term(hint[n], n) for n in prod))
        raise ValueError(msg)
    settings = {**DEFAULT_SETTINGS, **(raw.get("settings") or {})}
    t0, t1, dt = raw.get("T", DEFAULT_T)
    temps = [round(float(t), 6) for t in np.arange(float(t0), float(t1) + 0.5 * float(dt), float(dt))]
    if 298.15 not in temps:
        temps = sorted(temps + [298.15])
    cond = raw.get("conditions") or {}
    scan = raw.get("scan") or {}
    scan_sp = scan.get("species") or next(n for n, v in nu.items() if v < 0)
    if scan_sp not in nu:
        raise ValueError(f"scan species {scan_sp!r} is not in the equation")
    lig = raw.get("ligand")
    if isinstance(lig, str):
        lig = {"species": lig}
    if lig:
        if lig["species"] not in nu:
            raise ValueError(f"ligand {lig['species']!r} is not in the equation")
        lig = {"species": lig["species"], "dG_eV": float(lig.get("dG_eV", 0.0))}
    recipe = _recipe(raw.get("recipe"), nu, used)
    return {
        "name": raw.get("name") or src.stem, "equation": raw["equation"], "nu": nu, "ligand": lig,
        "energy_method": str(raw.get("energy_method", "PBE")), "equilibria": sides, "recipe": recipe,
        "species": {n: species[n] for n in used}, "unused_species": unused, "T": temps,
        "conditions": {"eps": cond.get("eps", 2.38), "cavity": cond.get("cavity", "beta0"),
                       "electrostatics": cond.get("electrostatics", "cosmo")},
        "scan": {"species": scan_sp, "log10_c": list(scan.get("log10_c", [-8.0, 1.0]))},
        "per": raw.get("per"), "settings": settings, "source": str(src),
    }


def _recipe(rc: dict | None, nu: dict, used: list) -> dict | None:
    """The experiment as prepared: dots at c_dot, the reagent in equivalents per dot, an optional additive."""
    if not rc:
        return None
    out = {"dot": rc.get("dot") or next(n for n, v in nu.items() if v < 0),
           "reagent": rc.get("reagent"), "additive": rc.get("additive"),
           "c_dot": float(rc.get("c_dot", 1.0e-5)), "equivalents": float(rc.get("equivalents", 100.0)),
           "c_additive": float(rc.get("c_additive", 0.0))}
    if nu.get(out["dot"], 0) >= 0:
        raise ValueError(f"recipe: the dot {out['dot']!r} must be a reactant of the equation")
    if out["reagent"] is None or nu.get(out["reagent"], 0) >= 0:
        raise ValueError(f"recipe: the reagent {out['reagent']!r} must be a reactant of the equation")
    if out["additive"] is not None and out["additive"] not in used:
        raise ValueError(f"recipe: the additive {out['additive']!r} is in neither the equation nor an equilibrium")
    return out


def _term(c, name):
    return f"{c} {name}" if c != 1 else name


ELEMENTS = ("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr "
            "Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm "
            "Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi").split()


def Z(symbol: str) -> int:
    return ELEMENTS.index(symbol) + 1
