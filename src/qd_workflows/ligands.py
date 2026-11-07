"""
Carboxylate X-type ligands on the nucleation grid (CdSe pilot, ligand block of nucleation.json).

The lattice configurations keep a Cl placeholder on every X site: the moves (lattice_grid) are decided on them, as for
CdCl2. Every structure that is relaxed, annealed or given a Hessian carries the carboxylate instead, placed by
QD_builder's X-type exchange (donor O on the X site, tail along the facet normal, clashes avoided). The dressed
structure has a fixed atom order: the lattice configuration's atoms in their order, each placeholder replaced by the
molecule's atoms in the template order (donor O first as the builder writes it), so a lattice index maps to a fixed
slice of the dressed structure (offsets) and a removal deletes whole molecules.

    ligand block   {"name": "oleate", "smiles": <optional, else LIGANDS[name]>, "seed": 7, "ff": "mmff", "d3": true,
                    "md_timestep_fs": 1.0, "inner": "cbeta", "relax_batch_atoms": 12000,
                    optional: "h_mass": 4.0 (amu, anneal MD only, with md_timestep_fs 2.0), "quench_screen_steps": 50 and
                    "quench_keep": 8 (two-stage anneal quench), "checkpoint_calls": 1}
    D3             D3(BJ), PBE parameters (torch-dftd), added to every batched MACE evaluation (install_d3 patches
                   orchestr_ai.qd.batch_relax.BatchMACE, which relaxations, MD and FD Hessians all build at call time)
    Hessians       partial (PHVA): the block of the inner region, core Cd/Se, every COO, C_alpha and C_beta with their H;
                   the rest of each chain fixed but in the energy (mechanical embedding). When the inner region is the
                   whole structure (acetate) the full Hessian with translations and rotations projected out.
    references     the Cd(L)2 monomer, dressed from linear CdCl2, relaxed from several conformers, the lowest kept,
                   its frequencies with the same inner rule; written as the grid's MX_q reference (bookkeeping name
                   CdCl2, the placeholder's) next to the campaign CdSe monomer (+ its D3 energy)

    python -m qd_workflows.ligands refs --wd WD --model M --refs CAMPAIGN_REFS --out REFS   (GPU)
    python -m qd_workflows.ligands dress --wd WD IN.xyz OUT.xyz                            (laptop, check)
"""
from __future__ import annotations

import argparse
import contextlib
import functools
import io
import json
import os
import sys
import tempfile
from pathlib import Path

import numpy as np

LIGANDS = {
    "acetate": "CC(=O)O",
    "butyrate": "CCCC(=O)O",
    "hexanoate": "CCCCCC(=O)O",
    "oleate": "CCCCCCCC/C=C\\CCCCCCCC(=O)O",
}
CLASH_A = 1.0          # Å: a dressed addition with a new atom this close to an old one falls back to the lattice start
D_MATCH = 1e-3         # Å: the donor O sits exactly on the placeholder (builder)
FD_MAX_ATOMS = 12000   # atoms per MACE call in the FD Hessian (float64, ~2.7 MB per atom measured)


class DressError(RuntimeError):
    pass


@contextlib.contextmanager
def _quiet():
    with contextlib.redirect_stdout(io.StringIO()):
        yield


def _cached_prepare():
    """The builder re-embeds the ligand for every exchange; one conformer per (SMILES, seed, ff) is enough."""
    from builder import ligand_exchange_posttreat as lep
    if not getattr(lep._prepare_charged_ligand, "_cached", False):
        inner = lep._prepare_charged_ligand

        @functools.lru_cache(maxsize=16)
        def cached(smiles, charge, seed, ff):
            return inner(smiles, charge, seed=seed, ff=ff)

        def wrapper(smiles, charge, seed, ff):
            return cached(smiles, charge, seed, ff)

        wrapper._cached = True
        lep._prepare_charged_ligand = wrapper


class Dressing:
    """Cl placeholders -> carboxylate molecules with QD_builder's exchange, in the fixed atom order (module doc)."""

    def __init__(self, cif, charges, name="acetate", smiles=None, seed=7, ff="mmff", placeholder="Cl",
                 surf_tol=2.0, facets=None, **_ignored):
        import yaml
        from pymatgen.core import Structure
        from builder.config import parse_yaml_config
        from builder.facets import expand_facets
        self.name, self.smiles = name, smiles or LIGANDS[name]
        self.placeholder, self.cif = placeholder, str(cif)
        self.charges = {k: int(v) for k, v in charges.items()}
        self.native = [e for e in self.charges if e != placeholder]
        facets = facets or [{"hkl": "100", "scope": "family", "termination": "cation_rich", "gamma": 1.0},
                            {"hkl": "111", "scope": "family", "termination": "cation_rich", "gamma": 1.0},
                            {"hkl": "-1-1-1", "scope": "family", "termination": "anion_rich", "gamma": 1.0}]
        recipe = {"size_unit_cells": [2, 2, 2], "facets": facets,
                  "passivation": {"ligand": placeholder, "surf_tol": float(surf_tol)},
                  "charges": self.charges, "symmetry": {"proper_rotations_only": True},
                  "post_treatment": {"ligand_exchange": {"enabled": True, "seed": int(seed), "ff": ff, "passes": [
                      {"replace": placeholder, "charge": -1, "smiles": self.smiles, "ratio": 1.0,
                       "distribution": "uniform"}]}}}
        with tempfile.NamedTemporaryFile("w", suffix=".yaml", delete=False) as fh:
            yaml.safe_dump(recipe, fh)
            path = fh.name
        try:
            self.cfg = parse_yaml_config(path)
        finally:
            os.unlink(path)
        self.struct = Structure.from_file(self.cif)
        self.facet_seeds = expand_facets(self.struct, self.cfg.seeds, proper_only=self.cfg.proper_only)
        _cached_prepare()
        self.template = None                     # (symbols of one molecule, donor index)
        s, p = self.dress([self.native[0], placeholder, placeholder],
                          np.array([[0.0, 0.0, 0.0], [2.3, 0.0, 0.0], [-2.3, 0.0, 0.0]]), planes=[])
        self.monomer = s
        self.template_pts = np.asarray(p[1:1 + self.n_lig], float)      # one molecule's coordinates (graft)

    # -- the fixed atom order --
    @property
    def n_lig(self):
        return len(self.template[0])

    def expand(self, syms):
        """Dressed symbols of a lattice configuration (placeholders -> the template's symbols)."""
        if self.template is None:
            raise DressError("no template yet")
        out = []
        for s in syms:
            out += list(self.template[0]) if s == self.placeholder else [s]
        return out

    def offsets(self, syms):
        """[(start, end)] in the dressed structure of every lattice atom."""
        n = [self.n_lig if s == self.placeholder else 1 for s in syms]
        ends = np.cumsum(n)
        return [(int(e - k), int(e)) for e, k in zip(ends, n)]

    def slices(self, syms, idx):
        """Dressed atom indices of the lattice atoms `idx`."""
        off = self.offsets(syms)
        return [j for i in sorted(idx) for j in range(*off[i])]

    # -- the exchange --
    def planes(self, syms, pts):
        from builder.facets import detect_facets_from_nc
        core = [i for i, s in enumerate(syms) if s in self.native]
        if len(core) < 4:
            return []
        return detect_facets_from_nc([syms[i] for i in core], np.asarray(pts, float)[core], self.struct.lattice,
                                     self.charges, self.facet_seeds, self.cfg.passivation.surf_tol)[1]

    def dress(self, syms, pts, planes=None):
        """(dressed symbols, positions) of a structure whose X are placeholders (all of them exchanged; atoms that are
        already ligand atoms stay where they are and count in the sterics)."""
        from builder.ligand_exchange_posttreat import run_ligand_exchange_posttreatment
        syms, pts = list(syms), np.asarray(pts, float)
        ph = [i for i, s in enumerate(syms) if s == self.placeholder]
        if not ph:
            return syms, pts
        planes = self.planes(syms, pts) if planes is None else planes
        with _quiet():
            out_s, out_p, ledger = run_ligand_exchange_posttreatment(syms, pts, self.cfg, self.struct, planes, self.cif)
        out_s, out_p = list(out_s), np.asarray(out_p, float)
        base = [i for i in range(len(syms)) if i not in set(ph)]
        nb = len(base)
        if out_s[:nb] != [syms[i] for i in base] or len(ledger) != len(ph):
            raise DressError(f"{len(ledger)} of {len(ph)} placeholders exchanged")
        rest_s, rest_p = out_s[nb:], out_p[nb:]
        if self.template is None:
            n_lig = len(rest_s) // len(ph)
            mol = rest_s[:n_lig]
            d1 = next(j for j in range(n_lig) if min(np.linalg.norm(rest_p[j] - pts[i]) for i in ph) < D_MATCH)
            self.template = (mol, d1)
        n_lig, d1 = self.n_lig, self.template[1]
        if len(rest_s) != n_lig * len(ph):
            raise DressError(f"{len(rest_s)} ligand atoms for {len(ph)} placeholders of {n_lig}")
        blocks = {}
        for b in range(len(ph)):
            blk = slice(b * n_lig, (b + 1) * n_lig)
            if rest_s[blk] != list(self.template[0]):
                raise DressError("a molecule's atoms are not in the template order")
            dists = [float(np.linalg.norm(rest_p[b * n_lig + d1] - pts[i])) for i in ph]
            i = ph[int(np.argmin(dists))]
            if min(dists) > D_MATCH or i in blocks:
                raise DressError(f"molecule {b} does not sit on one placeholder (closest {min(dists):.3f} Å)")
            blocks[i] = rest_p[blk]
        out_syms, out_pts, k = [], [], 0
        for i, s in enumerate(syms):
            if s == self.placeholder:
                out_syms += list(self.template[0])
                out_pts.append(blocks[i])
            else:
                out_syms.append(s)
                out_pts.append(out_p[k][None])
                k += 1
        return out_syms, np.vstack(out_pts)

    def grow(self, rel_syms, rel_pts, new_syms, new_pts):
        """A dressed parent (no placeholders left) plus new lattice atoms (cation, placeholders) at the end: the new
        placeholders exchanged against the relaxed shell. None when a new atom lands too close to an old one."""
        s, p = self.dress(list(rel_syms) + list(new_syms), np.vstack([rel_pts, new_pts]))
        n = len(rel_syms)
        if np.linalg.norm(p[n:, None] - np.asarray(rel_pts)[None], axis=-1).min() < CLASH_A:
            return None
        return s, p

    # -- a short carboxylate grown into this one (acetate grid -> oleate) --
    def graft(self, syms, pts, n_rot=12, passes=3):
        """A structure carrying a short carboxylate (acetate: each molecule's atoms contiguous) -> the same structure
        with this Dressing's carboxylate on every one of them, in this Dressing's atom order (each old molecule's
        block replaced by a template block). The binding is kept: per molecule the template's O, O, C and C_alpha
        fitted onto the old ones (both O assignments), the chain turned about the C-C_alpha and the C_alpha-C_beta
        axes in n_rot steps each, the pose with the least overlap (core, the heads of all molecules, the chains placed so far) kept; `passes`
        sweeps, the later ones against every other chain. Returns (symbols, positions, worst overlap score)."""
        from scipy.spatial import cKDTree
        syms, pts = list(syms), np.asarray(pts, float)
        mols = carboxylates(syms, pts)
        if not mols:
            return syms, pts, 0.0
        T = self.template_pts
        tm = carboxylates(list(self.template[0]), T)
        if len(tm) != 1:
            raise DressError("template is not one carboxylate")
        t = tm[0]
        moving = [j for j in range(len(T)) if j not in (t["O"][0], t["O"][1], t["C"])]
        beyond = _side(list(self.template[0]), T, t["Ca"], t["Cb"]) if t.get("Cb") is not None else []
        ts = np.array(self.template[0])
        heavy_t = ts != "H"
        core = [i for i, s in enumerate(syms) if s not in ORGANIC]
        heads = np.vstack([pts[[m["O"][0], m["O"][1], m["C"], m["Ca"]]] for m in mols])
        placed = [None] * len(mols)

        hm = heavy_t[moving]

        def overlap(X, tree, env_heavy):
            """sum of (r0 - d)^2 over close pairs, r0 3.0 (C/O-C/O, metal), 2.4 (with one H), 2.0 A (H-H)"""
            d = cKDTree(X).sparse_distance_matrix(tree, 3.0, output_type="ndarray")
            if not len(d):
                return 0.0
            j, i, r = d["i"], d["j"], d["v"]
            r0 = np.where(hm[j] & env_heavy[i], 3.0, np.where(hm[j] | env_heavy[i], 2.4, 2.0))
            return float((np.clip(r0 - r, 0, None) ** 2).sum())

        def poses(m):
            ref = pts[[m["O"][0], m["O"][1], m["C"], m["Ca"]]]
            for oa in ((t["O"][0], t["O"][1]), (t["O"][1], t["O"][0])):
                src = T[[oa[0], oa[1], t["C"], t["Ca"]]]
                R, sh = _kabsch(src, ref)
                X0 = T @ R.T + sh
                c, axis = X0[t["C"]], X0[t["Ca"]] - X0[t["C"]]
                axis = axis / np.linalg.norm(axis)
                for a in np.linspace(0, 2 * np.pi, n_rot, endpoint=False):
                    X1 = X0.copy()
                    X1[moving] = _rotate(X0[moving] - c, axis, a) + c
                    if not beyond:
                        yield X1
                        continue
                    c2, ax2 = X1[t["Ca"]], X1[t["Cb"]] - X1[t["Ca"]]
                    ax2 = ax2 / np.linalg.norm(ax2)
                    for b2 in np.linspace(0, 2 * np.pi, n_rot, endpoint=False):
                        X = X1.copy()
                        X[beyond] = _rotate(X1[beyond] - c2, ax2, b2) + c2
                        yield X

        score = [0.0] * len(mols)
        for sweep in range(passes):
            for b, m in enumerate(mols):
                env = [pts[core]] + [heads]
                env_s = [np.array(syms)[core], np.array(["O", "O", "C", "C"] * len(mols))]
                for c_, X in enumerate(placed):
                    if X is not None and c_ != b:
                        env.append(X[moving])
                        env_s.append(ts[moving])
                envp, envh = np.vstack(env), np.concatenate(env_s) != "H"
                own = set(range(4 * b, 4 * b + 4))              # this molecule's own head in `heads`
                keep = np.array([i not in {len(core) + o for o in own} for i in range(len(envp))])
                tree, eh = cKDTree(envp[keep]), envh[keep]
                best = min(((overlap(X[moving], tree, eh), k_, X) for k_, X in enumerate(poses(m))), key=lambda x: x[0])
                placed[b], score[b] = best[2], best[0]
        out_s, out_p, first = [], [], {m["atoms"][0]: b for b, m in enumerate(mols)}
        in_mol = {i for m in mols for i in m["atoms"]}
        for i, s in enumerate(syms):
            if i in first:
                out_s += list(self.template[0])
                out_p.append(placed[first[i]])
            elif i not in in_mol:
                out_s.append(s)
                out_p.append(pts[i][None])
        return out_s, np.vstack(out_p), float(max(score))


ORGANIC = ("C", "H", "O")


def carboxylates(syms, pts):
    """Every carboxylate molecule (connected C/H/O piece holding a C with two O): its atoms (sorted, contiguous
    expected), the two O (donor first: the one closest to a metal, else the lower index), the carboxyl C and C_alpha."""
    from ase import Atoms
    from ase.neighborlist import NeighborList, natural_cutoffs
    from scipy.sparse.csgraph import connected_components
    syms, pts = list(syms), np.asarray(pts, float)
    org = [i for i, s in enumerate(syms) if s in ORGANIC]
    if not org:
        return []
    a = Atoms([syms[i] for i in org], positions=pts[org])
    nl = NeighborList(natural_cutoffs(a), self_interaction=False, bothways=True)
    nl.update(a)
    n, lab = connected_components(nl.get_connectivity_matrix())
    metal = np.array([p for s, p in zip(syms, pts) if s not in ORGANIC]).reshape(-1, 3)
    out = []
    for c in range(n):
        loc = [j for j in range(len(org)) if lab[j] == c]
        cc = None
        for j in loc:
            nb = [int(x) for x in nl.get_neighbors(j)[0]]
            os_ = [x for x in nb if a[x].symbol == "O"]
            cs = [x for x in nb if a[x].symbol == "C"]
            if a[j].symbol == "C" and len(os_) == 2 and len(cs) == 1:
                cc = (j, os_, cs[0])
                break
        if cc is None:
            raise DressError(f"organic piece of {len(loc)} atoms without a carboxylate")
        j, os_, ca = cc
        if len(metal):
            os_ = sorted(os_, key=lambda x: float(np.linalg.norm(metal - a.positions[x], axis=1).min()))
        cb = [int(x) for x in nl.get_neighbors(ca)[0] if a[x].symbol == "C" and x != j]
        out.append({"atoms": sorted(org[x] for x in loc), "O": [org[os_[0]], org[os_[1]]], "C": org[j], "Ca": org[ca],
                    "Cb": org[cb[0]] if cb else None})
    return sorted(out, key=lambda m: m["atoms"][0])


def _side(syms, pts, a, b):
    """Atoms on b's side of the bond a-b (b included) of one molecule."""
    from ase import Atoms
    from ase.neighborlist import NeighborList, natural_cutoffs
    at = Atoms(syms, positions=pts)
    nl = NeighborList(natural_cutoffs(at), self_interaction=False, bothways=True)
    nl.update(at)
    seen, todo = {b}, [b]
    while todo:
        i = todo.pop()
        for j in nl.get_neighbors(i)[0]:
            j = int(j)
            if j not in seen and not (i == b and j == a):
                seen.add(j)
                todo.append(j)
    return sorted(seen)


def _kabsch(P, Q):
    cp, cq = P.mean(0), Q.mean(0)
    U, _s, Vt = np.linalg.svd((P - cp).T @ (Q - cq))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return R, cq - R @ cp


def _rotate(v, axis, a):
    """Rodrigues rotation of the rows of v about the unit axis by a."""
    k = np.asarray(axis)
    return v * np.cos(a) + np.cross(k, v) * np.sin(a) + np.outer(v @ k, k) * (1 - np.cos(a))


# ---------------- D3 ----------------
class _D3:
    """D3(BJ) energy and forces of single structures on the GPU (torch-dftd, PBE parameters)."""

    def __init__(self, device="cuda", dtype="float32"):
        import torch
        from torch_dftd.torch_dftd3_calculator import TorchDFTD3Calculator
        self.calc = TorchDFTD3Calculator(atoms=None, device=device, damping="bj", xc="pbe",
                                         dtype=torch.float64 if dtype == "float64" else torch.float32)

    def __call__(self, items):
        from ase import Atoms
        out = []
        for syms, pts in items:
            a = Atoms(list(syms), positions=np.asarray(pts, float))
            a.calc = self.calc
            out.append((float(a.get_potential_energy()), np.asarray(a.get_forces(), float)))
        return out


def install_d3():
    """Make every BatchMACE built from now on add D3(BJ) to its energies and forces (relaxations, MD, FD Hessians)."""
    from orchestr_ai.qd import batch_relax as br
    if getattr(br.BatchMACE, "_d3", False):
        return
    base = br.BatchMACE

    class BatchMACED3(base):
        _d3 = True

        def __init__(self, calc, max_atoms: int = 4000, descriptors: bool = False):
            super().__init__(calc, max_atoms, descriptors)
            self.d3 = _D3(str(calc.device), calc.default_dtype)

        def __call__(self, structures):
            res = super().__call__(structures)
            tup = [(j, x) for j, x in enumerate(structures) if isinstance(x, tuple)]
            if not tup:
                return res                          # graphs: displaced() adds D3 itself
            d3 = self.d3([x for _j, x in tup])
            res = list(res)
            for (j, _x), (e, f) in zip(tup, d3):
                res[j] = (res[j][0] + e, res[j][1] + f, *res[j][2:])
            return res

        def displaced(self, symbols, pos0, displacements):
            res = super().displaced(symbols, pos0, displacements)
            d3 = self.d3([(symbols, np.asarray(pos0, float) + d) for d in displacements])
            return [(r[0] + e, r[1] + f, *r[2:]) for r, (e, f) in zip(res, d3)]

    br.BatchMACE = BatchMACED3


def install_oeq():
    """QDW_OEQ=1: every MACE calculator built from now on uses OpenEquivariance kernels (the CUDA-only cuEquivariance
    has no ROCm build; OpenEquivariance runs on AMD GPUs, e.g. LUMI's MI250X). orchestr_ai's mace_calculator imports
    MACECalculator at call time, so replacing the class covers relaxations, MD and Hessians."""
    import mace.calculators as mc
    if os.environ.get("QDW_OEQ", "0") != "1" or getattr(mc.MACECalculator, "_oeq", False):
        return
    base = mc.MACECalculator

    class MACECalculatorOEQ(base):
        _oeq = True

        def __init__(self, *a, **kw):
            if not kw.get("enable_cueq"):
                kw["enable_oeq"] = True
            super().__init__(*a, **kw)

    mc.MACECalculator = MACECalculatorOEQ


def _gpu_worker(conn, spec):
    """One GPU of install_multigpu: a process that sees only its device (ROCR/HIP/CUDA_VISIBLE_DEVICES set before it
    started), with its own MACE (+ OpenEquivariance, + D3 as in the parent), evaluating what the parent sends."""
    install_oeq()
    if spec["d3"]:
        install_d3()
    from orchestr_ai.qd import batch_relax as br
    from orchestr_ai.qd.engines import mace_calculator
    calc = mace_calculator(spec["head"], spec["model"], "cuda", spec["dtype"])
    evs = {}
    while True:
        msg = conn.recv()
        if msg is None:
            return
        kind, key, args = msg
        if key not in evs:
            evs[key] = br.BatchMACE(calc, *key)
        try:
            conn.send(("ok", evs[key](*args) if kind == "call" else evs[key].displaced(*args)))
        except Exception as exc:                 # noqa: BLE001  (reported to the parent, which raises)
            conn.send(("err", f"{type(exc).__name__}: {exc}"))


_POOLS = {}


def _pool(spec, n_dev):
    """Worker processes for devices 1 .. n_dev - 1 (the parent keeps device 0), one pool per model/head/dtype/D3,
    started once and kept for the run."""
    import atexit
    import multiprocessing as mp
    key = (spec["model"], spec["head"], spec["dtype"], spec["d3"])
    if key in _POOLS:
        return _POOLS[key]
    ctx = mp.get_context("spawn")
    vis = os.environ.get("ROCR_VISIBLE_DEVICES") or os.environ.get("HIP_VISIBLE_DEVICES")
    ids = vis.split(",") if vis else [str(i) for i in range(n_dev)]
    saved = {v: os.environ.get(v) for v in ("ROCR_VISIBLE_DEVICES", "HIP_VISIBLE_DEVICES", "CUDA_VISIBLE_DEVICES")}
    # spawn re-imports the parent's main module in every worker, which reruns a script without a __main__ guard
    # (redress.py, the benchmarks): hide it while the workers start
    main = sys.modules["__main__"]
    main_file, main_spec = getattr(main, "__file__", None), getattr(main, "__spec__", None)
    if main_file is not None:
        del main.__file__
    main.__spec__ = None
    conns = []
    for i in range(1, n_dev):
        a, b = ctx.Pipe()
        if os.environ.get("ROCR_VISIBLE_DEVICES") is not None or "rocm" in __import__("torch").__version__:
            os.environ["ROCR_VISIBLE_DEVICES"] = ids[i]
            os.environ.pop("HIP_VISIBLE_DEVICES", None)
        else:
            os.environ["CUDA_VISIBLE_DEVICES"] = ids[i]
        pr = ctx.Process(target=_gpu_worker, args=(b, spec), daemon=True)
        pr.start()
        conns.append((a, pr))
    for v, x in saved.items():
        if x is None:
            os.environ.pop(v, None)
        else:
            os.environ[v] = x
    if main_file is not None:
        main.__file__ = main_file
    main.__spec__ = main_spec

    def stop():
        for c, pr in conns:
            try:
                c.send(None)
            except Exception:                    # noqa: BLE001
                pass
            pr.join(5)

    atexit.register(stop)
    _POOLS[key] = conns
    return conns


def install_multigpu(n_dev=None):
    """QDW_NGPU=n (> 1): every BatchMACE built from now on spreads each call over n GPUs of the node (LUMI: the 8
    MI250X GCDs): the structures split into n groups of about equal atom count, device 0 in this process and one worker
    process per other device (each sees only its GPU and builds its own MACE, OpenEquivariance and D3 as here), the
    groups evaluated at the same time, results in input order. The grid's relaxations, the anneal MD and the FD
    Hessians all go through BatchMACE, so one grid runs on the whole node. Processes, not threads: OpenEquivariance's
    kernels fault when one process drives several GCDs."""
    import torch
    from orchestr_ai.qd import batch_relax as br
    n_dev = int(n_dev or os.environ.get("QDW_NGPU", "1"))
    n_dev = min(n_dev, torch.cuda.device_count()) if torch.cuda.is_available() else 1
    if n_dev <= 1 or getattr(br.BatchMACE, "_multi", False):
        return
    base = br.BatchMACE
    import mace.calculators as mc
    if not getattr(mc.MACECalculator, "_records", False):
        calc_cls = mc.MACECalculator

        class MACECalculatorRec(calc_cls):       # remembers how it was built, for the worker processes
            _records = True

            def __init__(self, *a, **kw):
                self._qdw_init = (a, dict(kw))
                super().__init__(*a, **kw)

        mc.MACECalculator = MACECalculatorRec

    class MultiBatchMACE:
        _multi = True
        _d3 = getattr(base, "_d3", False)        # install_d3 must not wrap it again
        _size = staticmethod(base._size)

        def __init__(self, calc, max_atoms: int = 4000, descriptors: bool = False):
            self.calc, self.model = calc, calc.models[0]
            self.max_atoms, self.descriptors = int(max_atoms), descriptors
            self.local = base(calc, max_atoms, descriptors)
            self.workers = []
            if hasattr(calc, "_qdw_init"):
                kw = calc._qdw_init[1]
                spec = {"model": kw.get("model_paths"), "head": kw.get("head"), "dtype": kw.get("default_dtype"),
                        "d3": bool(self._d3)}
                if spec["model"] and spec["head"]:
                    self.workers = _pool(spec, n_dev)

        def __getattr__(self, name):             # anything else: the local evaluator
            if name == "local":
                raise AttributeError(name)
            return getattr(self.local, name)

        def _split(self, sizes):
            """Contiguous groups of about equal atom count (at most 1 + workers, none empty)."""
            n = min(1 + len(self.workers), len(sizes))
            per, groups, cur, acc = sum(sizes) / n, [], [], 0
            for i, s in enumerate(sizes):
                cur.append(i)
                acc += s
                if acc >= per * (len(groups) + 1) and len(groups) < n - 1:
                    groups.append(cur)
                    cur = []
            if cur:
                groups.append(cur)
            return groups

        def _spread(self, kind, groups, args_of, local_fn):
            key = (self.max_atoms, self.descriptors)
            for (c, _pr), idx in zip(self.workers, groups[1:]):
                c.send((kind, key, args_of(idx)))
            res = [local_fn(*args_of(groups[0]))]
            for (c, _pr), _idx in zip(self.workers, groups[1:]):
                st, r = c.recv()
                if st != "ok":
                    raise RuntimeError(f"GPU worker: {r}")
                res.append(r)
            return res

        def __call__(self, structures):
            structures = list(structures)
            if not structures or not self.workers or not all(isinstance(x, tuple) for x in structures):
                return self.local(structures)
            groups = self._split([base._size(x) for x in structures])
            if len(groups) == 1:
                return self.local(structures)
            res = self._spread("call", groups, lambda idx: ([(list(structures[i][0]), np.asarray(structures[i][1]))
                                                              for i in idx],), self.local)
            out = [None] * len(structures)
            for idx, r in zip(groups, res):
                for i, x in zip(idx, r):
                    out[i] = x
            return out

        def displaced(self, symbols, pos0, displacements):
            displacements = [np.asarray(d) for d in displacements]
            if not self.workers or len(displacements) < 2:
                return self.local.displaced(symbols, pos0, displacements)
            groups = self._split([len(symbols)] * len(displacements))
            res = self._spread("displaced", groups,
                               lambda idx: (list(symbols), np.asarray(pos0), [displacements[i] for i in idx]),
                               self.local.displaced)
            return [x for r in res for x in r]

    br.BatchMACE = MultiBatchMACE


def configure(settings, lig: dict):
    """Settings for a carboxylate grid: the ligand block on the settings (lattice_grid and grid_thermo read it),
    D3 installed when asked."""
    settings.ligand = dict(lig)
    if lig.get("relax_batch_atoms"):            # atoms per batched MACE call (float64 ~2.7 MB per atom)
        settings.relax_batch_atoms = int(lig["relax_batch_atoms"])
    # oleate relaxations take thousands of LBFGS steps: checkpoint every batched call, or a group of quenches
    # (4 calls' worth by default) outlasts a 1 h piece and starts over in every piece
    from orchestr_ai.qd.steps import detachment as dt
    dt.CHECKPOINT_CALLS = int(lig.get("checkpoint_calls", 1))
    from . import sxl_engine
    if sxl_engine.enabled():                    # QDW_SXL=1: MACE from Symmetrix-XL, D3 from the chunked torch-dftd copy
        sxl_engine.install(d3=lig.get("d3", True))
        return settings
    install_oeq()
    if lig.get("d3", True):
        install_d3()
    install_multigpu()                          # after D3: each device's evaluator carries its own D3
    return settings


# ---------------- partial Hessian ----------------
def inner_atoms(syms, pts, depth="cbeta"):
    """Core Cd/Se (every non C/H/O atom), every carboxylate COO, C_alpha (and C_beta for depth 'cbeta') and the H on
    them. All atoms for acetate."""
    from ase import Atoms
    from ase.neighborlist import NeighborList, natural_cutoffs
    a = Atoms(list(syms), positions=np.asarray(pts, float))
    org = {"C", "H", "O"}
    nl = NeighborList(natural_cutoffs(a), self_interaction=False, bothways=True)
    nl.update(a)
    nb = [[int(j) for j in nl.get_neighbors(i)[0] if syms[j] in org] for i in range(len(syms))]
    keep = {i for i, s in enumerate(syms) if s not in org}
    levels = {"calpha": 1, "cbeta": 2}[depth]
    for c in (i for i, s in enumerate(syms) if s == "C"):
        os_ = [j for j in nb[c] if syms[j] == "O"]
        if len(os_) < 2:
            continue
        keep |= {c, *os_}
        front, seen = [c], {c}
        for _ in range(levels):
            nxt = [j for i in front for j in nb[i] if syms[j] == "C" and j not in seen]
            seen |= set(nxt)
            keep |= set(nxt)
            front = nxt
        keep |= {j for i in seen for j in nb[i] if syms[j] == "H"}
    return sorted(keep)


def frequencies(syms, pts, settings, depth=None):
    """Frequencies (cm^-1) of the inner region (partial FD Hessian, MACE + D3 if installed); the whole structure with
    translations and rotations projected out when the inner region is all of it."""
    from ase import Atoms
    from orchestr_ai.qd import batch_relax as br
    from orchestr_ai.qd.engines import mace_calculator
    from orchestr_ai.qd.steps.hessian import vibrations
    depth = depth or (getattr(settings, "ligand", None) or {}).get("inner", "cbeta")
    syms, pts = list(syms), np.asarray(pts, float)
    inner = inner_atoms(syms, pts, depth)
    ev = br.BatchMACE(mace_calculator(settings.head, settings.model, settings.device, "float64"), FD_MAX_ATOMS)
    cols = [3 * i + a for i in inner for a in range(3)]
    h = br.fd_hessian(syms, pts, ev, delta=settings.fd_step, columns=None if len(inner) == len(syms) else cols)
    masses = Atoms(syms, positions=pts).get_masses()
    if len(inner) == len(syms):
        return vibrations(h, pts, masses)[0]
    hb = h[np.ix_(cols, cols)]
    hb = 0.5 * (hb + hb.T)
    w = 1.0 / np.sqrt(np.repeat(masses[inner], 3))
    lam = np.linalg.eigvalsh(hb * w[:, None] * w[None, :])
    from orchestr_ai.qd.steps.hessian import EIG_TO_CM1
    return np.sign(lam) * np.sqrt(np.abs(lam)) * EIG_TO_CM1


# ---------------- references ----------------
def _wd_dressing(wd, lig):
    rec = Path(wd) / "start"
    r = json.loads((rec / "record.json").read_text())
    recipe = r["origin"]["recipe"]
    return Dressing(rec / r["origin"]["cif"], recipe["charges"], surf_tol=recipe["passivation"].get("surf_tol", 2.0),
                    facets=recipe.get("facets"), **lig), recipe


def references(wd, model, head, campaign_refs, out, n_conf=4, log=print):
    """The pilot's reference folder: CdCl2-<ligand> = the Cd(L)2 monomer (lowest of n_conf conformers, relaxed, rattled,
    re-relaxed and polished like the grid, frequencies by the inner rule) and the campaign CdSe monomer (D3 added)."""
    from orchestr_ai.qd.steps import detachment as dt
    from types import SimpleNamespace
    from .nucleation import _block, _settings
    wd, out = Path(wd), Path(out)
    P = _block(json.loads((wd / "nucleation.json").read_text()))
    lig = P["ligand"]
    s = configure(_settings(model, head), lig)
    root = out / head / "monomers"
    root.mkdir(parents=True, exist_ok=True)
    # CdSe: the campaign monomer, D3 added at its geometry
    src = sorted((Path(campaign_refs) / head / "monomers").glob("CdSe-*/monomer.json"))[0]
    m = json.loads(src.read_text())
    if lig.get("d3", True):
        from orchestr_ai.qd import batch_relax as br
        from orchestr_ai.qd.engines import mace_calculator
        ev = br.BatchMACE(mace_calculator(head, model, s.device, s.dtype), 4000)
        m["energy_eV_mace"] = m["energy_eV"]
        m["energy_eV"] = float(ev([(m["symbols"], np.asarray(m["positions"], float))])[0][0])
        m["d3"] = "bj-pbe"
    (root / "CdSe-pilot").mkdir(exist_ok=True)
    (root / "CdSe-pilot" / "monomer.json").write_text(json.dumps(m))
    # Cd(L)2 from linear CdCl2, n_conf builder conformers (seeds)
    cache = dt.RelaxCache(SimpleNamespace(props=out, settings=s, record_dir=out))
    best = None
    for c in range(n_conf):
        dz, _recipe = _wd_dressing(wd, {**lig, "seed": int(lig.get("seed", 7)) + 100 * c})
        syms, pts = dz.dress(["Cd", "Cl", "Cl"], np.array([[0.0, 0, 0], [2.3, 0, 0], [-2.3, 0, 0]]), planes=[])
        (a, e, _ok), = cache.relax_many([(syms, pts)], s)
        p = a.get_positions() + np.random.default_rng(c).normal(0.0, 0.05, (len(syms), 3))
        (a, e, _ok), = cache.relax_many([(syms, p)], s)
        (a, e, _ok), = cache.relax_many([(syms, a.get_positions())], s, polish=True)
        log(f"[ligands] Cd({lig['name']})2 conformer {c}: E = {e:.4f} eV")
        if best is None or e < best[2]:
            best = (syms, a.get_positions(), e, c)
    syms, pts, e, c = best
    freqs = frequencies(syms, pts, s)
    mono = {"symbols": syms, "positions": np.asarray(pts).round(6).tolist(), "energy_eV": float(e),
            "frequencies_cm1": np.asarray(freqs).round(4).tolist(), "conformer": c, "ligand": lig,
            "note": "the grid's MX_q reference: bookkeeping name CdCl2 (placeholder), structure Cd(L)2",
            "inner_atoms": len(inner_atoms(syms, pts, lig.get("inner", "cbeta")))}
    (root / f"CdCl2-{lig['name']}").mkdir(exist_ok=True)
    (root / f"CdCl2-{lig['name']}" / "monomer.json").write_text(json.dumps(mono))
    log(f"[ligands] references in {root}: Cd({lig['name']})2 E = {e:.4f} eV (conformer {c}), "
        f"{sum(np.asarray(freqs) < -10)} imaginary modes")
    return root


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python -m qd_workflows.ligands")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("refs", help="the Cd(L)2 and CdSe references of a carboxylate grid (GPU)")
    p.add_argument("--wd", required=True)
    p.add_argument("--model", required=True)
    p.add_argument("--head", default="omat_pbe")
    p.add_argument("--refs", required=True, help="the campaign references (CdSe monomer)")
    p.add_argument("--out", required=True)
    p.add_argument("--conformers", type=int, default=4)
    p = sub.add_parser("dress", help="dress one xyz with the wd's ligand (laptop check)")
    p.add_argument("--wd", required=True)
    p.add_argument("inp")
    p.add_argument("out")
    a = ap.parse_args(argv)
    if a.cmd == "refs":
        references(a.wd, a.model, a.head, a.refs, a.out, a.conformers, log=lambda x: print(x, flush=True))
    else:
        from .nucleation import _block
        from .xyz import read_xyz, write_xyz
        P = _block(json.loads((Path(a.wd) / "nucleation.json").read_text()))
        dz, _r = _wd_dressing(a.wd, P["ligand"])
        s, p = read_xyz(a.inp)
        ds, dp = dz.dress(list(s), np.asarray(p, float))
        write_xyz(a.out, ds, dp, f"dressed with {dz.name}")
        print(f"{len(s)} -> {len(ds)} atoms")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
