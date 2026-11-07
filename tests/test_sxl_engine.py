"""ChunkedD3 reproduces torch-dftd's TorchDFTD3Calculator(damping="bj", xc="pbe", float64) on a small cluster, with
the pairs in one block and in many small ones."""
import numpy as np
import pytest

pytest.importorskip("torch_dftd")


def _cluster(n=40, seed=3):
    rng = np.random.default_rng(seed)
    symbols = list(rng.choice(["Cd", "Se", "O", "C", "H"], size=n))
    pos = []
    while len(pos) < n:                               # random points at least 1.0 Å apart in a 9 Å box
        p = rng.uniform(0.0, 9.0, 3)
        if all(np.linalg.norm(p - q) > 1.0 for q in pos):
            pos.append(p)
    return symbols, np.array(pos)


@pytest.mark.parametrize("block_pairs", [4_000_000, 37])
def test_chunked_d3_matches_torch_dftd(block_pairs):
    import torch
    from ase import Atoms
    from torch_dftd.torch_dftd3_calculator import TorchDFTD3Calculator

    from qd_workflows.sxl_engine import ChunkedD3

    symbols, pos = _cluster()
    ref = TorchDFTD3Calculator(atoms=None, device="cpu", damping="bj", xc="pbe", dtype=torch.float64)
    a = Atoms(symbols, positions=pos)
    a.calc = ref
    e_ref, f_ref = a.get_potential_energy(), a.get_forces()
    (e, f), = ChunkedD3("cpu", "float64", block_pairs=block_pairs)([(symbols, pos)])
    assert abs(e - e_ref) < 1e-10
    assert np.abs(f - f_ref).max() < 1e-10
