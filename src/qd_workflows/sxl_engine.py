"""
Symmetrix-XL as the MACE engine of the batched relaxations (orchestr_ai.qd.batch_relax.BatchMACE), with D3(BJ) from a
chunked, bounded-memory copy of torch-dftd. Opt-in: QDW_SXL=1 (qd_workflows.ligands.configure calls install()).

Same model, same numbers: MACE-MH-1 omat_pbe in float64 through Symmetrix-XL's ASE calculator (same checkpoint
exported to its compact JSON; against mace-torch e3nn -1.5e-7 eV/atom, max |dF| 2e-4 eV/A, from its radial splines),
plus D3 equal to torch-dftd's TorchDFTD3Calculator(damping="bj", xc="pbe", float64) to 1e-13 eV and 2e-15 eV/A.
Phase 0 (4 Nov 2026, 4 dressed acetate dots, production relax + polish): LUMI GCD 776 -> 242 s, NHR A100 228 -> 123 s
(one dot at a time; install() makes relax_many do that), final energies within -11..+46 meV of mace-torch.

Two things Symmetrix-XL needs:
  - a cell: non-periodic atoms have a zero cell, which it inverts. Each dot goes into a periodic box with VACUUM Å
    around it (MACE cutoff 6 Å: no image within reach; identical to the open cluster to 5e-13 eV/A). A fully periodic
    cell also lets it keep the neighbour graph on the device between steps.
  - no built-in D3: Symmetrix(dispersion=True) adds the D3 energy but no D3 forces (v0.1.1) and is slow on HIP.

D3 in bounded memory: torch-dftd keeps every pair within 95 bohr (all pairs of a 4-5 nm dot) with 25 reference C6
values in one autograd graph (23 GB at 4,632 atoms; out of memory at 6,957 on an 80 GB A100). ChunkedD3 makes the
pairs in row blocks and differentiates each block before the next:
  1. CN_i = sum_j damp(r_ij) over pairs within cnthr                                      (no graph)
  2. E = sum_pairs e(r_ij, CN_i, CN_j); per block, backward gives dE/dx (direct) and dE/dCN
  3. per block of CN pairs, backward of sum damp(r_ij) (dE/dCN_i + dE/dCN_j) adds the CN chain rule to dE/dx
2.5 GB at 4,632 atoms, 6 GB at 6,957, 3.3 GB at 13,914.

Environment: QDW_SXL=1 switches the engine; QDW_SXL_JSON = the exported model (symmetrix_extract_mace --model M
--head omat_pbe --output JSON); QDW_SXL_VACUUM (default 40 Å).
"""
from __future__ import annotations

import os

import numpy as np

VACUUM = 40.0                     # Å around each dot in its periodic box


class ChunkedD3:
    """D3(BJ) energy (eV) and forces (eV/Å) of non-periodic structures, equal to torch-dftd's defaults: two-body only,
    pairs with d < 95 bohr (ase Bohr, in Å), coordination numbers over pairs with r <= 40 bohr, no smoothing. The
    tables, constants and C6 interpolation are torch-dftd's own."""

    def __init__(self, device="cuda", dtype="float64", block_pairs: int = 4_000_000):
        import torch
        from ase.units import Bohr
        from torch_dftd.dftd3_xc_params import get_dftd3_default_params
        from torch_dftd.functions.dftd3 import d3_autoang, d3_autoev, d3_k1, d3_k3
        from torch_dftd.nn.dftd3_module import DFTD3Module
        self.torch = torch
        self.device = torch.device(device)
        self.dtype = torch.float64 if dtype in ("float64", torch.float64) else torch.float32
        self.params = get_dftd3_default_params(damping="bj", xc="pbe")
        m = DFTD3Module(self.params, dtype=self.dtype)
        self.c6ab = m.c6ab.to(self.device)
        self.rcov = m.rcov.to(self.device)
        self.r2r4 = m.r2r4.to(self.device)
        self.cut2_A = float(m.cutoff) ** 2
        self.cnthr_b = float(m.cnthr) / Bohr
        self.autoang, self.autoev, self.k1, self.k3 = d3_autoang, d3_autoev, d3_k1, d3_k3
        self.block_pairs = int(block_pairs)

    def _blocks(self, pos_A):
        """Pairs i < j with d < cutoff, a block of rows at a time."""
        torch = self.torch
        n = pos_A.shape[0]
        rows = max(1, self.block_pairs // max(n, 1))
        for a in range(0, n, rows):
            b = min(n, a + rows)
            d2 = ((pos_A[a:b, None, :] - pos_A[None, :, :]) ** 2).sum(-1)
            ii = torch.arange(a, b, device=pos_A.device)[:, None]
            jj = torch.arange(n, device=pos_A.device)[None, :]
            i, j = torch.nonzero((d2 < self.cut2_A) & (jj > ii), as_tuple=True)
            if i.numel():
                yield i + a, j

    def _damp_cn(self, Z, pos_b, i, j):
        r = self.torch.linalg.norm(pos_b[i] - pos_b[j], dim=1)
        keep = r <= self.cnthr_b
        i, j, r = i[keep], j[keep], r[keep]
        rr = (self.rcov[Z[i]] + self.rcov[Z[j]]).type(r.dtype) / r
        return i, j, 1.0 / (1.0 + self.torch.exp(-self.k1 * (rr - 1.0)))

    def _c6(self, Zi, Zj, nci, ncj):
        """torch_dftd.functions.dftd3._getc6_impl."""
        torch = self.torch
        cn0, cn1, cn2 = self.c6ab.reshape(-1, 25, 3).split(1, dim=2)
        index = Zi * self.c6ab.size(1) + Zj
        cn0 = cn0.squeeze(2)[index].type(nci.dtype)
        cn1 = cn1.squeeze(2)[index].type(nci.dtype)
        cn2 = cn2.squeeze(2)[index].type(nci.dtype)
        r = (cn1 - nci[:, None]) ** 2 + (cn2 - ncj[:, None]) ** 2
        w = torch.softmax(torch.where(cn0 > 0.0, self.k3 * r, -1.0e20), dim=1)
        return (w * cn0).sum(dim=1)

    def _pair_energy(self, Z, pos_b, cn, i, j):
        """Hartree, summed over the block (torch-dftd's -0.5 per direction, both directions)."""
        r = self.torch.linalg.norm(pos_b[i] - pos_b[j], dim=1)
        Zi, Zj = Z[i], Z[j]
        c6 = self._c6(Zi, Zj, cn[i], cn[j])
        c8 = 3 * c6 * self.r2r4[Zi].type(c6.dtype) * self.r2r4[Zj].type(c6.dtype)
        r2 = r ** 2
        r6 = r2 ** 3
        r8 = r6 * r2
        p = self.params
        tmp = p["rs6"] * self.torch.sqrt(c8 / c6) + p["rs18"]
        tmp2 = tmp ** 2
        tmp6 = tmp2 ** 3
        tmp8 = tmp6 * tmp2
        e = -(p["s6"] * c6 / (r6 + tmp6) + p["s18"] * c8 / (r8 + tmp8))
        return e.to(self.torch.float64).sum()

    def energy_forces(self, Z, pos):
        torch = self.torch
        Z = torch.as_tensor(np.asarray(Z), device=self.device, dtype=torch.long)
        pos_A = torch.as_tensor(np.asarray(pos, float), device=self.device, dtype=self.dtype)
        pos_b = (pos_A / self.autoang).detach().requires_grad_(True)
        with torch.no_grad():                                  # 1. coordination numbers
            cn = pos_b.new_zeros(Z.shape[0])
            for i, j in self._blocks(pos_A):
                i, j, d = self._damp_cn(Z, pos_b, i, j)
                cn.index_add_(0, i, d)
                cn.index_add_(0, j, d)
        cn_leaf = cn.clone().requires_grad_(True)
        e_h = 0.0
        for i, j in self._blocks(pos_A):                       # 2. energy, direct gradient, dE/dCN
            e = self._pair_energy(Z, pos_b, cn_leaf, i, j)
            e.backward()
            e_h += float(e.detach())
        if pos_b.grad is None:                                 # a single atom
            return 0.0, np.zeros((Z.shape[0], 3))
        g_cn = cn_leaf.grad.detach()
        for i, j in self._blocks(pos_A):                       # 3. CN chain rule
            i, j, d = self._damp_cn(Z, pos_b, i, j)
            if i.numel():
                (d * (g_cn[i] + g_cn[j])).sum().backward()
        f = -(self.autoev / self.autoang) * pos_b.grad.detach()
        return self.autoev * e_h, f.cpu().numpy().astype(float)

    def __call__(self, items):
        """[(symbols, positions)] -> [(energy, forces)], as qd_workflows.ligands._D3."""
        from ase.data import atomic_numbers
        return [self.energy_forces([atomic_numbers[s] for s in syms], pts) for syms, pts in items]


class SxlBatch:
    """Drop-in for orchestr_ai.qd.batch_relax.BatchMACE on Symmetrix-XL: (symbols, positions) pairs ->
    [(energy, forces)], one dot per call (Symmetrix-XL has no batch API), each in its periodic box; D3 from ChunkedD3
    when d3. displaced() serves the FD Hessian. The mace-torch calculator it is built with is not used for energies
    (orchestr_ai still loads it to make one)."""

    _d3 = True                                       # install_d3 must not wrap it with torch-dftd
    _sxl = True

    def __init__(self, calc=None, max_atoms: int = 4000, descriptors: bool = False):
        if descriptors:
            raise NotImplementedError("Symmetrix-XL engine: no MACE descriptors")
        self.calc_torch = calc
        self.max_atoms = int(max_atoms)
        self.descriptors = descriptors
        self.sxl = _calculator()
        self.d3 = _CONFIG["d3_engine"]() if _CONFIG["d3"] else None
        self.vacuum = float(os.environ.get("QDW_SXL_VACUUM", VACUUM))
        self.n_calls = 0
        self.n_oom = 0

    @staticmethod
    def _size(item) -> int:
        return len(item[0])

    def _boxed(self, symbols, x):
        from ase import Atoms
        lo = x.min(axis=0) - self.vacuum / 2
        return Atoms(list(symbols), positions=x - lo, cell=np.diag(x.max(axis=0) - lo + self.vacuum / 2), pbc=True)

    def _one(self, symbols, x):
        x = np.asarray(x, float)
        a = self._boxed(symbols, x)
        a.calc = self.sxl
        e, f = float(a.get_potential_energy()), np.asarray(a.get_forces(), float)
        if self.d3 is not None:
            (ed, fd), = self.d3([(symbols, x)])
            e, f = e + ed, f + fd
        self.n_calls += 1
        return e, f

    def __call__(self, structures):
        out = []
        for item in structures:
            if not isinstance(item, tuple):
                raise TypeError("Symmetrix-XL engine takes (symbols, positions) pairs, not prebuilt graphs")
            out.append(self._one(*item))
        return out

    def displaced(self, symbols, pos0, displacements):
        pos0 = np.asarray(pos0, float)
        return [self._one(symbols, pos0 + np.asarray(d)) for d in displacements]


_CONFIG = {"d3": True, "d3_engine": None, "json": None, "calc": None}


def _calculator():
    """One Symmetrix-XL calculator per process (float64, omat_pbe), shared by every SxlBatch."""
    if _CONFIG["calc"] is None:
        from symmetrix import Symmetrix
        path = _CONFIG["json"] or os.environ.get("QDW_SXL_JSON")
        if not path:
            raise RuntimeError("QDW_SXL_JSON: path of the Symmetrix-XL export of the MACE model")
        _CONFIG["calc"] = Symmetrix(path, dtype="float64", head=os.environ.get("QDW_SXL_HEAD", "omat_pbe"))
    return _CONFIG["calc"]


def enabled() -> bool:
    return os.environ.get("QDW_SXL", "0") == "1"


def _one_at_a_time(relax_many):
    """relax_many over an SxlBatch, one structure after the other. relax_many treats each structure on its own (its
    own L-BFGS, its own stopping test), so this is the same algorithm; on Symmetrix-XL alternating dots every step costs
    ~35 % more (NHR A100, 4 dressed dots: 178 s batched, 131 s one at a time; mace-torch 228 s)."""
    def wrapped(structures, evaluator, **kw):
        if not getattr(evaluator, "_sxl", False) or len(structures) < 2:
            return relax_many(structures, evaluator, **kw)
        return [relax_many([s], evaluator, **kw)[0] for s in structures]
    wrapped._sxl = True
    return wrapped


def install(d3: bool = True, json_path: str | None = None, device: str = "cuda"):
    """Every BatchMACE built from now on is an SxlBatch (MACE from Symmetrix-XL, D3 from ChunkedD3 when d3), and
    relax_many relaxes its structures one at a time when given one."""
    from orchestr_ai.qd import batch_relax as br
    if getattr(br.BatchMACE, "_sxl", False):
        return
    _CONFIG.update(d3=bool(d3), json=json_path, d3_engine=lambda: ChunkedD3(device, "float64"))
    br.BatchMACE = SxlBatch
    if not getattr(br.relax_many, "_sxl", False):
        br.relax_many = _one_at_a_time(br.relax_many)


def provenance() -> dict:
    """What produced the energies, for states.json settings."""
    import importlib.metadata as md
    vers = {}
    for p in ("symmetrix-xl", "symmetrix-xl-cuda12-sm80", "symmetrix-xl-rocm6-gfx90a", "torch-dftd"):
        try:
            vers[p] = md.version(p)
        except md.PackageNotFoundError:
            pass
    return {"engine": "symmetrix-xl", "versions": vers, "json": _CONFIG["json"] or os.environ.get("QDW_SXL_JSON"),
            "vacuum_A": float(os.environ.get("QDW_SXL_VACUUM", VACUUM)),
            "d3": "chunked torch-dftd bj/pbe" if _CONFIG["d3"] else None}
