"""CP2K cells on one common grid: every edge n x h, with CP2K choosing exactly n points."""
import re

import numpy as np
import pytest

from qd_workflows import TEMPLATES
from qd_workflows.cp2k.inputs import (CUTOFF_RY, aligned_length, cell_abc, fft_lengths, grid_points,
                                      grid_spacing)


def test_every_template_has_the_cutoff_of_the_grid():
    for f in (TEMPLATES / "cp2k").glob("*.inp"):
        cut = re.findall(r"^\s*CUTOFF\s+(\d+)", f.read_text(), flags=re.M)
        assert cut == [str(CUTOFF_RY)], f


def test_fft_lengths_are_cp2k_fftw3():
    n = fft_lengths(1000)
    assert {2, 3, 5, 7, 11, 27, 25, 77, 990} <= set(n)
    assert not {13, 49, 81, 125, 121} & set(n)          # 13 not allowed; at most 7^1, 3^3, 5^2, 11^1


def test_grid_points_follow_cp2k():
    # n_min = 2 floor(L sqrt(Ecut) / 2 pi) + 1 (bohr, Ry), then the next FFT length
    L = 36.557
    n_min = 2 * int(L / 0.529177210903 * 20 / (2 * np.pi)) + 1
    assert grid_points(L) == min(n for n in fft_lengths() if n >= n_min)


@pytest.mark.parametrize("want", [15.0, 16.564, 23.3, 36.557, 36.502, 58.0, 97.13])
def test_aligned_edges_share_one_spacing(want):
    h = grid_spacing()
    L = aligned_length(want)
    n = grid_points(L)
    assert want <= L < want * 1.06
    assert abs(L / n - h) < 1e-12                       # CP2K's own spacing is the common one
    assert grid_points(float(f"{L:.6f}")) == n          # still so after printing ABC with 6 decimals


def test_cell_abc_of_a_dot_and_of_a_monomer():
    pts = np.array([[0, 0, 0], [24.5, 0, 0], [0, 24.4, 0], [0, 0, 24.4]])
    abc = [float(x) for x in cell_abc(pts).split()]
    assert all(grid_points(x) * grid_spacing() == pytest.approx(x, abs=1e-6) for x in abc)
    raw = [float(x) for x in cell_abc(pts, cutoff_ry=None).split()]
    assert raw == [36.5, 36.4, 36.4] and all(a >= r for a, r in zip(abc, raw))


def test_monomer_cells_equal_orchestr_ai():
    field = pytest.importorskip("orchestr_ai.qd.field")
    base = "      CUTOFF  400\n"
    for pts in (np.array([[0, 0, 0], [0, 0, 2.41]]), np.array([[0, 0, 0], [0, 0, 2.28], [0, 0, -2.28]]),
                np.array([[0, 0, 0], [2.3, 0, 0], [-1.15, 1.99, 0], [-1.15, -1.99, 0]])):
        assert field.monomer_cell(base, pts) == cell_abc(pts, min_length=15.0)
