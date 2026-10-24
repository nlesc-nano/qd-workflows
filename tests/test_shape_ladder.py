"""Sampled CdCl2 ladders: removable units, site roles and stripping to the bare core (geometry only, no MACE).
Fixture: builder tetrahedra of zinc-blende CdSe, Cd-terminated (Cd35Se20Cl30) and Se-terminated (Cd35Se18Cl34)."""
import json
from collections import Counter
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("scipy")
from qd_workflows import shape_ladder as sl  # noqa: E402

DOTS = json.loads((Path(__file__).parent / "data" / "shapes" / "tetrahedra.json").read_text())


def bonded_cl_ok(symbols, positions):
    nb = sl.neighbours(np.asarray(positions, float), sl.BOND)
    return all(any(symbols[j] == "Cd" for j in nb[i]) for i, s in enumerate(symbols) if s == "Cl")


@pytest.mark.parametrize("name", ["tetra_A", "tetra_B"])
def test_units_are_neutral_cdcl2(name):
    d = DOTS[name]
    us = sl.units(d["symbols"], d["positions"], d["planes"])
    assert us
    for u in us:
        assert [d["symbols"][a] for a in u["atoms"]] == ["Cd", "Cl", "Cl"]
        s, p = sl.remove(d["symbols"], d["positions"], u["atoms"])
        assert bonded_cl_ok(s, p)                           # no Cl left without a Cd
        assert sl.connected(sl.neighbours(p, sl.LINK), set())


def test_roles_follow_the_facets():
    """Cd-terminated tetrahedron: its facet Cd carry no Cl, so CdCl2 comes off edges and vertices only;
    Se-terminated tetrahedron: units on {-1-1-1}B facets, edges and the four vertices."""
    a = Counter(u["role"] for u in sl.units(DOTS["tetra_A"]["symbols"], DOTS["tetra_A"]["positions"], DOTS["tetra_A"]["planes"]))
    b = Counter(u["role"] for u in sl.units(DOTS["tetra_B"]["symbols"], DOTS["tetra_B"]["positions"], DOTS["tetra_B"]["planes"]))
    assert set(a) == {"edge:111A|111A", "vertex"}
    assert set(b) == {"facet:111B", "edge:111B|111B", "vertex"} and b["vertex"] == 4


def test_strip_reaches_the_neutral_bare_core():
    d = DOTS["tetra_A"]
    s, p, n = d["symbols"], d["positions"], 0
    while True:
        us = sl.units(s, p, d["planes"])
        if not us:
            break
        s, p = sl.remove(s, p, us[0]["atoms"])
        n += 1
        assert 2 * s.count("Cd") - 2 * s.count("Se") - s.count("Cl") == 0
    assert n == d["k"] and s.count("Cl") == 0 and s.count("Cd") == s.count("Se")


def test_sampling_spreads_over_roles():
    d = DOTS["tetra_B"]
    us = sl.units(d["symbols"], d["positions"], d["planes"])
    picked = sl.sample(us, 3, np.random.default_rng(0), d["positions"])
    roles = Counter(u["role"] for u in picked)
    assert set(roles) == {"facet:111B", "edge:111B|111B", "vertex"}
    assert roles["vertex"] == 1                          # the four vertices of a Td tetrahedron are equivalent
    p = np.asarray(d["positions"], float)
    keys = [sl.site_key(p[u["atoms"][0]], p) for u in picked]
    assert len(set(keys)) == len(keys)                   # no symmetry copies
    rr = sl.round_robin(us, 4, np.random.default_rng(0))
    assert len({u["role"] for u in rr}) == 3
