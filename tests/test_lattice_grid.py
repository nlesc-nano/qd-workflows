"""Lattice moves of the MA detachment grid on the builder's Cd16Se13Cl6 (QD_builder rules, no MACE)."""
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("builder.z_type_ligand_posttreat")
pytest.importorskip("orchestr_ai.qd.records")
from qd_workflows.lattice_grid import Lattice  # noqa: E402
from qd_workflows.xyz import read_xyz  # noqa: E402

DATA = Path(__file__).parent / "data" / "cd16"


@pytest.fixture(scope="module")
def cd16():
    syms, pts = read_xyz(DATA / "start.xyz")
    return Lattice(DATA, DATA / "CdSe_zb.cif"), list(syms), np.asarray(pts)


@pytest.mark.parametrize("pick", [0, -1])
def test_ladder_up_fills_every_se_deficit(cd16, capsys, pick):
    """Whichever allowed addition is taken, the ladder fills the 12 Se deficits: one host site fewer per unit."""
    lat, s, p = cd16
    hosts = []
    while True:
        hosts.append(len(lat.z_sites(s, p)[0]))
        adds = lat.add_mx(s, p)
        if not adds:
            break
        _, s, p = adds[pick]
    capsys.readouterr()
    assert hosts == list(range(12, -1, -1))            # 12 Se sites at k = 3, none at k = 15
    assert (s.count("Cd"), s.count("Se"), s.count("Cl")) == (28, 13, 30)


def test_an_addition_that_blocks_a_host_is_dropped(cd16, capsys):
    lat, s, p = cd16
    sites, free = lat.z_sites(s, p)
    capsys.readouterr()
    assert lat.fillable(free)
    starved = [list(f) for f in free]
    starved[0] = starved[0][:1]                        # a host left with one free X site cannot take a CdCl2
    assert not lat.fillable(starved)


def test_cdse_removals_keep_the_inner_cd_the_graph_and_every_cd_on_se(cd16, capsys):
    """On the builder dot every surface Se can go with one of its two outer Cd (24 pairs, never a Cd of the central
    Se); each Z-type CdCl2 added protects its Se (two pairs fewer), so the fully passivated k = 15 has none."""
    lat, s, p = cd16
    counts = []
    while True:
        ma = lat.remove_ma(s, p)
        counts.append(len(ma))
        if len(counts) == 1:
            nb = lat.bonds(s, p)
            assert len({info["A"] for info, *_ in ma}) == 12
            assert all(sum(1 for j in nb[info["M"]] if s[j] == "Se") < 4 for info, *_ in ma)
            for _info, s2, p2 in ma:
                assert (s2.count("Cd"), s2.count("Se"), s2.count("Cl")) == (15, 12, 6)
                nb2 = lat.bonds(s2, p2)
                assert lat.connected(nb2, set())
                assert all(any(s2[j] == "Se" for j in nb2[i]) for i, e in enumerate(s2) if e == "Cd")
        adds = lat.add_mx(s, p)
        if not adds:
            break
        _, s, p = adds[0]
    capsys.readouterr()
    assert counts == list(range(24, -1, -2))


def test_relaxed_start_follows_the_parent_geometry(cd16, capsys):
    """start='relaxed': a removal deletes the atoms from the parent's relaxed coordinates; an addition puts the new
    atoms at their lattice sites moved by the host anion's displacement. start='lattice': the lattice positions."""
    from types import SimpleNamespace
    from qd_workflows.lattice_grid import LatticeGrid
    lat, s, p = cd16
    shift = np.array([0.1, -0.2, 0.05])
    parent = SimpleNamespace(syms=s, lat=p, relaxed=p + shift)
    g = LatticeGrid.__new__(LatticeGrid)
    info, s2, p2 = lat.remove_ma(s, p)[0]
    add_info, s3, p3 = lat.add_mx(s, p)[0]
    capsys.readouterr()
    g.start = "lattice"
    assert g._start(parent, info, s2, p2)[0] is p2
    g.start = "relaxed"
    start, how = g._start(parent, info, s2, p2)
    assert how == "relaxed" and np.allclose(start, p2 + shift)
    start, how = g._start(parent, add_info, s3, p3)
    assert how == "relaxed" and np.allclose(start, p3 + shift)


def test_ladder_down_from_the_builder_dot(cd16, capsys):
    lat, s, p = cd16
    counts = []
    while s.count("Cl"):
        rem = lat.remove_mx(s, p)
        counts.append(len(rem))
        _, s, p = rem[0]
    capsys.readouterr()
    assert counts == [24, 12, 4] and (s.count("Cd"), s.count("Se")) == (13, 13)
