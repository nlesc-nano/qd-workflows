"""Monomer library (qdw monomers): names, inputs as a chain writes them, linking a dot's runs, jobs."""
import json
from pathlib import Path

import pytest

from qd_workflows import monomers
from qd_workflows.cluster import load_cluster

BASIS = "".join(f"{e} DZVP-MOLOPT-GGA-GTH-q{q}\n1\n2 0 1 1 1 1\n 1.0 1.0 1.0\n"
                for e, q in (("Cd", 12), ("Se", 6), ("Cl", 7)))
POTENTIAL = "".join(f"{e} GTH-GGA-q{q}\n" for e, q in (("Cd", 12), ("Se", 6), ("Cl", 7)))
MONOMERS_XYZ = ("2\nMA CdSe\nCd 0.0 0.0 0.0\nSe 0.0 0.0 2.41\n"
                "3\nMX CdCl2\nCd 0.0 0.0 0.0\nCl 0.0 0.0 2.28\nCl 0.0 0.0 -2.28\n")
DOT_XYZ = "5\ndot\nCd 0 0 0\nSe 2.6 0 0\nCd 5.2 0 0\nCl -2.4 0 0\nCl 7.6 0 0\n"


def test_species_once_per_formula():
    sp = {s["name"]: s for s in monomers.species(["CdSe", "CdTe", "InP", "PbS"])}
    assert set(sp) == {"CdSe", "CdTe", "InP", "PbS", "CdCl2", "InCl3", "PbCl2"}
    assert sp["CdCl2"]["materials"] == ["CdSe", "CdTe"] and sp["CdCl2"]["label"] == "MX"
    assert sp["InP"]["label"] == "MA"
    with pytest.raises(ValueError):
        monomers.split_material("Cd2Se")


def test_normalized_keeps_the_monomer_kinds_only():
    inp = ("  PROJECT x  # note\n    &KIND Cd\n      BASIS_SET  B-q12\n    &END\n"
           "    &KIND Se\n      BASIS_SET  B-q6\n    &END\n\n   CUTOFF   400\n")
    assert monomers.normalized(inp, ["Cd"]) == "PROJECT x\n&KIND Cd\nBASIS_SET B-q12\n&END\nCUTOFF 400"


@pytest.fixture
def chain_and_library(tmp_path, monkeypatch):
    """A dot's field and esp monomer runs prepared as the chain does, and the library's runs."""
    field_mod = pytest.importorskip("orchestr_ai.qd.field")
    from qd_workflows.cp2k import esp
    from qd_workflows.cp2k.inputs import prepare
    data = tmp_path / "data"
    data.mkdir()
    (data / "BASIS_MOLOPT_UZH").write_text(BASIS)
    (data / "POTENTIAL_UZH").write_text(POTENTIAL)
    monkeypatch.setenv("QDW_CP2K_DATA", str(data))
    cluster = load_cluster("nhr")
    mono = tmp_path / "monomers.xyz"
    mono.write_text(MONOMERS_XYZ)
    # the chain of one dot
    run = tmp_path / "cp2k"
    (run).mkdir()
    (tmp_path / "dot.xyz").write_text(DOT_XYZ)
    prepare(str(tmp_path / "dot.xyz"), "sample", "Dot-field", str(run / "field" / "base"), cluster,
            wfn="../../geo_opt/Dot-RESTART.wfn")
    field_mod.prepare(str(run / "field" / "base" / "cp2k.inp"), str(run / "field"), str(mono))
    esp.prepare(str(run / "field" / "base" / "cp2k.inp"), str(run / "esp"), monomers=str(mono))
    # the library (stage 1)
    lib = tmp_path / "lib"
    monomers.write_cp2k(str(lib), str(mono), cluster)
    return run, lib, mono, cluster


def _finish(d: Path, cubes: bool):
    (d / "cp2k.out").write_text(" Hirshfeld Charges\n *** SCF run converged in 9 steps ***\n PROGRAM ENDED AT now\n")
    if cubes:
        (d / "v_hartree-1.cube").write_text("x")
        (d / "density-1.cube").write_text("x")


def test_library_inputs_equal_the_chains(chain_and_library):
    run, lib, mono, _ = chain_and_library
    for name, label, sym in (("CdSe", "MA", ["Cd", "Se"]), ("CdCl2", "MX", ["Cd", "Cl", "Cl"])):
        for mine, entry in ((run / "field" / f"mono_{label}", lib / name / "hirshfeld"),
                            (run / "esp" / label, lib / name / "esp")):
            assert monomers.normalized((mine / "cp2k.inp").read_text(), sym) == \
                   monomers.normalized((entry / "cp2k.inp").read_text(), sym)


def test_link_replaces_finished_entries_only(chain_and_library):
    run, lib, mono, _ = chain_and_library
    rows = monomers.link(str(run / "field"), str(run / "esp"), str(mono), str(lib))
    assert {r["status"] for r in rows} == {"library run not finished"}
    for name in ("CdSe", "CdCl2"):
        _finish(lib / name / "hirshfeld", cubes=False)
        _finish(lib / name / "esp", cubes=True)
    assert monomers.status(str(lib))["done"] == ["CdCl2", "CdSe"]
    rows = monomers.link(str(run / "field"), str(run / "esp"), str(mono), str(lib))
    assert [r["status"] for r in rows] == ["linked"] * 4
    assert (run / "field" / "mono_MX").is_symlink() and (run / "esp" / "MA").resolve() == (lib / "CdSe" / "esp").resolve()
    # the steps' own checks see them finished
    from qd_workflows.cp2k import esp
    assert "MA" not in esp.status(str(run / "esp")).split()


def test_link_refuses_another_input_or_geometry(chain_and_library):
    run, lib, mono, _ = chain_and_library
    for name in ("CdSe", "CdCl2"):
        _finish(lib / name / "hirshfeld", cubes=False)
        _finish(lib / name / "esp", cubes=True)
    inp = lib / "CdCl2" / "hirshfeld" / "cp2k.inp"
    inp.write_text(inp.read_text().replace("CUTOFF  400", "CUTOFF  600"))
    g = lib / "CdSe" / "esp" / "geom.xyz"
    g.write_text(g.read_text().replace("2.41000000", "2.45000000"))
    rows = {r["run"]: r["status"] for r in monomers.link(str(run / "field"), str(run / "esp"), str(mono), str(lib))}
    assert rows["field/mono_MX"].startswith("differs: input") and rows["esp/MA"] == "differs: geometry"
    assert rows["field/mono_MA"] == "linked" and not (run / "field" / "mono_MX").is_symlink()


def test_finished_entries_are_never_rewritten(chain_and_library):
    _run, lib, mono, cluster = chain_and_library
    for k in ("hirshfeld", "esp"):
        _finish(lib / "CdSe" / k, cubes=k == "esp")
    (lib / "CdSe" / "esp" / "cp2k.inp").write_text("kept")
    monomers.write_cp2k(str(lib), str(mono), cluster)
    assert (lib / "CdSe" / "esp" / "cp2k.inp").read_text() == "kept"
    assert monomers.todo(str(lib)) == [str(lib / "CdCl2" / "hirshfeld"), str(lib / "CdCl2" / "esp")]


def test_prepare_jobs(tmp_path):
    cluster = load_cluster("nhr")
    res = monomers.prepare_jobs(str(tmp_path / "job"), cluster, "nhr", ["CdSe", "CdTe", "InP"])
    assert [s["name"] for s in res["species"]] == ["CdSe", "CdCl2", "CdTe", "InP", "InCl3"]
    gpu = (tmp_path / "job" / "gpu.sh").read_text()
    assert "orchestr_ai.qd.references monomers CdSe CdCl2 CdTe InP InCl3 --ensemble" in gpu
    cpu = (tmp_path / "job" / "cpu.sh").read_text()
    assert "#SBATCH -N 1" in cpu and f"--ntasks-per-node={cluster['cores_per_node']}" in cpu
    assert "cp2k_n 8" in cpu and "--dependency=afterok" in (tmp_path / "job" / "submit.sh").read_text()
    assert json.loads((tmp_path / "job" / "species.json").read_text())[1]["materials"] == ["CdSe", "CdTe"]


def test_test_queue_and_campaign_profile(tmp_path):
    cluster = load_cluster("nhr_campaign")
    monomers.prepare_jobs(str(tmp_path / "t"), cluster, "nhr_campaign", ["CdSe"], test_queue=True)
    assert "#SBATCH -p gpu-a100:test" in (tmp_path / "t" / "gpu.sh").read_text()
    cpu = (tmp_path / "t" / "cpu.sh").read_text()
    assert "#SBATCH -p cpu-clx:test" in cpu and "--ntasks-per-node=96" in cpu and "--mem=330G" in cpu
    assert "PYTHONPATH=/scratch/usr/becivinf/qdw_campaign/code/Orchestr.AI/src" in (tmp_path / "t" / "gpu.sh").read_text()


def test_job_paths_are_absolute(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monomers.prepare_jobs("rel/job", load_cluster("nhr_campaign"), "nhr_campaign", ["InP"])
    for f in ("gpu.sh", "cpu.sh"):
        assert f"{tmp_path.resolve()}/rel/job/monomers.xyz" in (tmp_path / "rel" / "job" / f).read_text()
