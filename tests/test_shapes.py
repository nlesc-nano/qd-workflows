"""Shape library: polyhedron geometry and the gamma-ratio deduplication (no builder needed)."""
import numpy as np

from qd_workflows.shapes import gamma_sets, polyhedron


def test_cube_geometry():
    planes = [(np.array(v, float), 1.0, "100") for v in
              ([1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1])]
    g = polyhedron(planes)
    assert np.isclose(g["areas_A2"]["100"], 24.0)
    assert g["edges"]["100|100"]["n"] == 12 and np.isclose(g["edges"]["100|100"]["length_A"], 24.0)
    assert g["n_vertices"] == 8 and np.isclose(g["volume_A3"], 8.0)


def test_truncated_cube_mixes_facet_types():
    planes = [(np.array(v, float), 1.0, "100") for v in
              ([1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1])]
    planes += [(np.array(v, float) / np.sqrt(3), 2.5 / np.sqrt(3), "111A") for v in
               ([1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1])]
    g = polyhedron(planes)
    assert set(g["areas_A2"]) == {"100", "111A"} and "100|111A" in g["edges"]
    assert g["edges"]["100|111A"]["n"] == 12                       # four corner triangles, three edges each


def test_gamma_sets_drop_equal_ratios():
    fam = [{"hkl": "100", "gamma": [None, 1.0, 1.2]}, {"hkl": "111", "termination": "cation_rich", "gamma": [None, 1.0, 1.2]}]
    sets = gamma_sets(fam)
    keys = {tuple((f["hkl"], f["gamma"]) for f in s) for s in sets}
    assert len(sets) == len(keys) == 5          # {100}, {111}, 1:1, 1:1.2, 1.2:1 (1.2:1.2 is 1:1)
