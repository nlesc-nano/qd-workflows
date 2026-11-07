"""MACE evaluation for the GPU steps (reference, states, neb) and the free-energy pieces they share.

Engine wraps two batched MACE evaluators (float32 for relaxations and NEB forces, float64 for energies and Hessians)
from Orchestr.AI; every call takes [(symbols, positions)] and returns [(energy, forces)]. Tests pass any object with
ev32, ev64 and relax_many methods instead. QDW_OEQ=1 swaps in OpenEquivariance kernels (no cuEquivariance on ROCm)."""
from __future__ import annotations

import math
import os

import numpy as np

from .. import sxl_engine

KB, NA = 8.617333262e-5, 6.02214076e23
CM_PER_SQRT_EV_A2_AMU = 521.47090            # sqrt(eV / (A^2 amu)) in cm^-1


def install_oeq():
    """Every MACE calculator built from now on uses OpenEquivariance kernels when QDW_OEQ=1 (LUMI's MI250X)."""
    if os.environ.get("QDW_OEQ", "0") != "1":
        return
    import mace.calculators as mc
    if getattr(mc.MACECalculator, "_oeq", False):
        return
    base = mc.MACECalculator

    class MACECalculatorOEQ(base):
        _oeq = True

        def __init__(self, *a, **kw):
            if not kw.get("enable_cueq"):
                kw["enable_oeq"] = True
            super().__init__(*a, **kw)

    mc.MACECalculator = MACECalculatorOEQ


class Engine:
    """Batched MACE-MH-1 on one GPU (head, model path from the recipe)."""

    def __init__(self, model: str, head: str = "omat_pbe", device: str = "cuda", max_atoms: int = 8000):
        install_oeq()
        sxl_engine.configure()
        import orchestr_ai.qd.batch_relax as br
        from orchestr_ai.qd.engines import mace_calculator
        import torch
        self._torch, self._br = torch, br
        self._b32 = br.BatchMACE(mace_calculator(head, model, device, "float32"), max_atoms=max_atoms)
        self._b64 = br.BatchMACE(mace_calculator(head, model, device, "float64"), max_atoms=max_atoms)

    def ev32(self, items):
        self._torch.set_default_dtype(self._torch.float32)
        return self._b32(items)

    def ev64(self, items):
        self._torch.set_default_dtype(self._torch.float64)
        return self._b64(items)

    def relax_many(self, items, fmax, max_steps, precision="float32"):
        """L-BFGS of every (symbols, positions) together; results with .positions, .energy, .steps, .reason."""
        ev = self.ev32 if precision == "float32" else self.ev64
        return self._br.relax_many(items, ev, fmax=fmax, max_steps=max_steps)


def masses(symbols):
    from ase import Atoms
    return Atoms(symbols).get_masses()


def phva(engine, symbols, x, local, h=0.01, vectors=False):
    """Frequencies (cm^-1, negative = imaginary) of the partial Hessian over the atoms `local`: central differences of
    +-h A on each of their coordinates, float64, mass-weighted, symmetrised. vectors=True: also the mass-weighted
    normal modes (columns, over the 3 len(local) coordinates)."""
    x = np.asarray(x, float)
    items = []
    for i in local:
        for k in range(3):
            for sgn in (1, -1):
                y = x.copy()
                y[i, k] += sgn * h
                items.append((symbols, y))
    F = [np.asarray(f, float) for _e, f in engine.ev64(items)]
    nl = len(local)
    hm = np.zeros((3 * nl, 3 * nl))
    c = 0
    for a in range(nl):
        for k in range(3):
            hm[3 * a + k] = -(F[c][local] - F[c + 1][local]).reshape(-1) / (2 * h)
            c += 2
    m = np.repeat(masses(symbols)[list(local)], 3)
    w, v = np.linalg.eigh(0.5 * (hm + hm.T) / np.sqrt(np.outer(m, m)))
    fr = (np.sign(w) * np.sqrt(np.abs(w)) * CM_PER_SQRT_EV_A2_AMU).tolist()
    return (fr, v) if vectors else fr


def vib_free(freqs, T, n_free=0, floor=50.0):
    """Harmonic vibrational free energy (eV) over the real modes, the n_free smallest |freq| dropped (translation,
    rotation), modes below `floor` cm^-1 held at the floor."""
    fr = sorted(freqs, key=abs)[int(n_free):]
    f = np.maximum(np.asarray([v for v in fr if v > 0], float), floor) * 1.239841984e-4
    kt = KB * T
    return float(np.sum(0.5 * f + kt * np.log1p(-np.exp(-f / kt))))


def gas_trans_rot(symbols, x, T, c0=1.0):
    """Translational (at c0 M) plus rotational free energy (eV) of a free molecule, as the 2026-10 runs: the
    rotational partition function is halved for every molecule (symmetry number 2: right for CdCl2, but CdSe has
    sigma = 1, so G(CdSe) is kT ln 2 = 0.030 eV too high at 500 K). Kept so the results match; to be decided."""
    from ase import Atoms
    at = Atoms(symbols, positions=x)
    m = at.get_masses().sum() * 1.66053907e-27
    kt_j = 1.380649e-23 * T
    lam = 6.62607015e-34 / math.sqrt(2 * math.pi * m * kt_j)
    g_trans = -KB * T * math.log(1.0 / (c0 * NA * 1e3) / lam ** 3)
    inertia = np.asarray(at.get_moments_of_inertia()) * 1.66053907e-47
    inertia = inertia[inertia > 1e-50]
    hb = 1.054571817e-34
    q = (2 * inertia.mean() * kt_j / hb ** 2 / 2) if len(inertia) == 2 else \
        math.sqrt(math.pi * np.prod(2 * inertia * kt_j / hb ** 2)) / 2
    return g_trans - KB * T * math.log(q)


def free_molecule(engine, symbols, x0, T, h=0.01, floor=50.0, fmax=0.001):
    """Relaxed free molecule (float64): E, F_vib over 3N - 5 modes, translation at 1 M and rotation, G."""
    x = np.asarray(engine.relax_many([(symbols, np.asarray(x0, float) + 1e-3)], fmax=fmax, max_steps=500,
                                     precision="float64")[0].positions, float)
    fr = phva(engine, symbols, x, list(range(len(symbols))), h)
    e = float(engine.ev64([(symbols, x)])[0][0])
    out = {"E": e, "F_vib": vib_free(fr, T, 5, floor), "G_trans_rot_1M": gas_trans_rot(symbols, x, T), "freqs": fr,
           "positions": np.round(x, 5).tolist()}
    out["G"] = out["E"] + out["F_vib"] + out["G_trans_rot_1M"]
    return out
