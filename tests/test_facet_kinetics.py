"""Facet kinetics geometry: next-layer sites and the monomer unit on lattice sites, on an ideal zinc-blende cube."""
import math

import numpy as np
import pytest

pytest.importorskip("scipy")
from qd_workflows import facet_kinetics as fk  # noqa: E402

A = 6.203


def _cube(n=2):
    """Cd on the fcc sites of a (2n+1)^3 cell block, Se at Cd + a/4 (1,1,-1) style bonds; {100} planes."""
    fcc = np.array([[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0]])
    cd = np.array([(np.array([i, j, k]) + f) * A for i in range(-n, n) for j in range(-n, n) for k in range(-n, n)
                   for f in fcc])
    se = cd + A / 4 * np.array([1, 1, -1])
    half = n * A - 0.5
    cd = cd[np.all(np.abs(cd) <= half, axis=1)]
    se = se[np.all(np.abs(se) <= half, axis=1)]
    sym = ["Cd"] * len(cd) + ["Se"] * len(se)
    pos = np.vstack([cd, se])
    planes = [([s * (i == 0), s * (i == 1), s * (i == 2)], "100") for i in range(3) for s in (1, -1)]
    return sym, pos, planes


def test_bond_vectors_are_tetrahedral():
    sym, pos, _ = _cube()
    b = fk.bond_vectors(sym, pos)
    assert np.allclose(np.linalg.norm(b, axis=1), fk.BOND)
    assert np.allclose(b.sum(0), 0.0, atol=1e-9)                       # tetrahedron
    assert np.allclose([b[i] @ b[j] for i in range(4) for j in range(i + 1, 4)], -fk.BOND ** 2 / 3, atol=1e-6)


def test_unit_on_lattice_sites_keeps_the_composition():
    sym, pos, planes = _cube()
    surf = fk.Surface(sym, pos, pos, planes, "100")
    picked = surf.pick(2, 0)
    assert picked
    site, (pre, inc) = picked[0]
    assert site["class"] == "terrace" and site["anchor"] >= 1
    assert pre.shape == inc.shape == (len(sym) + 5, 3)
    assert np.allclose(pre[:len(sym)], inc[:len(sym)])                  # no Cl to displace on the bare cube
    new = inc[len(sym):]
    # Se and the two Cd sit on lattice sites: bonded to the cube or to each other at the bulk bond
    d = np.linalg.norm(new[0] - new[1:3], axis=1)
    assert np.allclose(d, fk.BOND, atol=1e-6)
    # nothing closer than the clash distance, and the precursor sits outward of the incorporated unit
    from scipy.spatial import cKDTree
    assert cKDTree(inc).query(inc, k=2)[0][:, 1].min() >= fk.CLASH - 1e-9
    n = np.asarray(surf.n)
    assert (pre[len(sym):] @ n).mean() > (inc[len(sym):] @ n).mean() + 2.0
    assert math.isclose(np.linalg.norm(n), 1.0)


def test_far_shell_moves_are_detected_only_away_from_the_unit():
    a = np.zeros((6, 3))
    a[:, 0] = [0, 1, 2, 20, 21, 30]                                   # atoms 0-2: the unit's neighbourhood
    b = a.copy()
    b[1, 1] += 2.0                                                    # near the unit: part of the event
    b[4, 1] += 3.0                                                    # far from it: a shell that rearranged
    b[5, 1] += 0.5                                                    # far, but small
    assert fk.Kinetics._far_moved(a, b, [0]) == [4]
