"""Transport geometry: walk-on-spheres capacitance against exact values, solvent-excluded radius, Stokes-Einstein."""
import math

import numpy as np
import pytest

pytest.importorskip("scipy")
from qd_workflows import transport as tr  # noqa: E402

CUBE = [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]]


def test_capacitance_of_a_sphere_is_its_radius():
    c, se = tr.capacitance_spheres(["Se"], [[0.0, 0.0, 0.0]], n=40000)
    assert abs(c - 1.90) < 3 * se + 0.01


def test_capacitance_of_a_cube():
    c, se = tr.capacitance_polyhedron(CUBE, [1.0] * 6, n=60000)      # edge a = 2
    assert abs(c / 2.0 - 0.6607) < 3 * se / 2.0 + 0.005                # Hwang & Mascagni: 0.66068 a


def test_capacitance_of_a_grown_polyhedron():
    # a point-like cube grown by R is the sphere of radius R
    c, se = tr.capacitance_polyhedron(CUBE, [1e-3] * 6, n=20000, dilate=1.0)
    assert abs(c - 1.0) < 3 * se + 0.01
    # the cube grown by R lies between C + R (exact for a sphere) and the cube of planes moved out by R
    c0, _ = tr.capacitance_polyhedron(CUBE, [1.0] * 6, n=20000)
    cd, se = tr.capacitance_polyhedron(CUBE, [1.0] * 6, n=20000, dilate=0.5)
    cs, _ = tr.capacitance_polyhedron(CUBE, [1.5] * 6, n=20000)
    assert c0 + 0.5 < cd - 2 * se and cd + 2 * se < cs


def test_grown_polyhedron_is_reproducible():
    # the zinc-blende families at the distances 1, 1.1, 1.05, 1.2 ({100}, {111}A, {111}B, {110}): seed 0 gives
    # these values bit for bit (the growth runs' capacitances must not change with the implementation)
    from qd_workflows.growth import ZB_NORMALS
    normals, dists = [], []
    for f, hf in zip(ZB_NORMALS, (1.0, 1.1, 1.05, 1.2)):
        normals += [np.asarray(v, float) / np.linalg.norm(v) for v in ZB_NORMALS[f]]
        dists += [hf] * len(ZB_NORMALS[f])
    assert tr.capacitance_polyhedron(normals, dists, n=10000)[0] == 1.1431583179794313
    assert tr.capacitance_polyhedron(normals, dists, n=10000, dilate=0.3)[0] == 1.4638862127473176
    assert tr.capacitance_polyhedron(normals, dists, n=10000, dilate=1.0)[0] == 2.160811915462148


def test_two_touching_spheres():
    # two equal touching spheres: C = 2 a ln 2 = 1.3863 a (exact)
    a = 1.75
    c, se = tr.capacitance_spheres(["Cl", "Cl"], [[-a, 0, 0], [a, 0, 0]], n=60000)
    assert abs(c - 2 * a * math.log(2)) < 3 * se + 0.01


def test_solvent_excluded_radius():
    # a single atom: the probe closing changes nothing, R_se = r_vdW
    assert abs(tr.ses_radius(["Se"], [[0, 0, 0]], probe=2.0, spacing=0.1) - 1.90) < 0.05
    # two atoms 3 A apart: the probe fills the neck, so the excluded volume exceeds the plain union
    two = tr.ses_volume(["Cd", "Se"], [[0, 0, 0], [3.0, 0, 0]], probe=2.0, spacing=0.15)
    bare = tr.ses_volume(["Cd", "Se"], [[0, 0, 0], [3.0, 0, 0]], probe=0.0, spacing=0.15)
    assert two > bare


def test_stokes_einstein():
    d = tr.stokes_d(3.5, 500.0, 0.4)                                   # A^2/s
    assert math.isclose(d * 1e-16, 1.380649e-23 * 500 / (6 * math.pi * 0.4e-3 * 3.5e-10) * 1e4, rel_tol=1e-9)
    assert math.isclose(tr.stokes_d(3.5, 500.0, 0.4, c=4.0) / d, 1.5)
