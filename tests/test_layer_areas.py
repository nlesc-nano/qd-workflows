"""Collision areas of the layer steps on the built models (network_export.layer_areas)."""
import itertools

import numpy as np

from qd_workflows import network_export as nx

U = 16                                    # lattice units per cubic cell (the registry keys': Cd at Se + (4, 4, 4))


def _cube(n):
    """Zinc-blende sites of n x n x n cells (Se on the fcc lattice, Cd a quarter diagonal away)."""
    se = [(x, y, z) for x, y, z in itertools.product(range(0, n * U, U // 2), repeat=3) if (x + y + z) % U == 0]
    return {("Se", *c) for c in se} | {("Cd", c[0] + 4, c[1] + 4, c[2] + 4) for c in se}


def _sheet(z, xs, ys):
    """One Se + Cd sheet above the plane z (inside the x, y ranges, in lattice units)."""
    out = set()
    for x, y in itertools.product(range(xs[0], xs[1], U // 2), range(ys[0], ys[1], U // 2)):
        if (x + y + z) % U == 0:
            out |= {("Se", x, y, z), ("Cd", x + 4, y + 4, z + 4)}
    return out


def _key(sites):
    return "(" + ", ".join(f"('{e}', {x}, {y}, {z})" for e, x, y, z in sorted(sites)) + ", ('Cl', 0, 0, -40))"


def test_lattice_scale_from_the_relaxed_structures():
    a = _cube(2)
    sites = sorted(nx.sites(_key(a)))
    assert ("Cl", 0, 0, -40) not in nx.sites(_key(a)) and len(sites) == len(a)
    rng = np.random.default_rng(1)
    pos = np.array([c[1:] for c in sites], float) * 0.38 + rng.normal(0, 0.02, (len(sites), 3))
    assert abs(nx.lattice_scale([_key(a)], [([c[0] for c in sites], pos)]) - 0.38) < 0.01


def test_sheet_areas_and_their_overlap():
    """A sheet covers part of the parent's surface: two alternative sheets on one facet overlap (their summed areas
    exceed what they cover together), sheets on opposite facets do not."""
    a = _cube(3)
    top = 3 * U
    s1, s2 = _sheet(top, (0, 24), (0, 48)), _sheet(top, (16, 48), (0, 48))
    s3 = _sheet(-U // 2, (0, 48), (0, 48))
    keys = {0: _key(a), 1: _key(a | s1), 2: _key(a | s2), 3: _key(a | s3)}
    pairs = [(0, 1, [0, 0, 1], "100"), (0, 2, [0, 0, 1], "100"), (0, 3, [0, 0, -1], "100")]
    areas, overlap = nx.layer_areas(pairs, keys, 0.38, 1.4, 1)
    assert all(v > 0 for v in areas.values()) and areas[(0, 3)] > areas[(0, 2)] > 0
    up = next(r for r in overlap if r["group"] == "normal" and r["normal"] == [0, 0, 1])
    fam = next(r for r in overlap if r["group"] == "family")
    assert up["n_steps"] == 2 and up["union_A2"] < 0.9 * up["sum_A2"]
    assert fam["n_steps"] == 3 and abs(fam["sum_A2"] - up["sum_A2"] - areas[(0, 3)]) < 1e-9
    assert up["union_A2"] + areas[(0, 3)] - fam["union_A2"] < 1e-6 * fam["union_A2"]   # opposite facets: no overlap


def test_child_in_its_own_frame():
    """The registry keeps each dot in its own frame: a child rotated and shifted gives the same area."""
    a = _cube(3)
    b = a | _sheet(3 * U, (0, 24), (0, 48))
    g = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])                        # 90 degrees about z, then a shift
    moved = {(e, *map(int, g @ np.array(c) + np.array([16, -32, 8]))) for e, *c in b}
    keys = {0: _key(a), 1: _key(b), 2: _key(moved)}
    areas, _ov = nx.layer_areas([(0, 1, [0, 0, 1], "100"), (0, 2, [0, 0, 1], "100")], keys, 0.38, 1.4, 1)
    assert areas[(0, 1)] > 0 and abs(areas[(0, 2)] / areas[(0, 1)] - 1) < 0.02     # a symmetric equivalent placement:
    back = nx.align(nx.sites(_key(a)), nx.sites(_key(moved)), [0, 0, 1])            # the same up to point sampling
    new = np.array([c[1:] for c in back if c not in nx.sites(_key(a))])
    assert len(back) == len(b) and len(new) == len(b - a) and new[:, 2].min() >= 3 * U
