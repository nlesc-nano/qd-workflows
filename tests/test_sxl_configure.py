"""sxl_engine.configure(): nothing without QDW_SXL=1; with it, every BatchMACE built afterwards is an SxlBatch, MACE
alone unless QDW_SXL_D3=1, and relax_many goes one structure at a time."""
import sys
import types

from qd_workflows import sxl_engine


def _stub_batch_relax(monkeypatch):
    br = types.ModuleType("orchestr_ai.qd.batch_relax")
    br.BatchMACE = type("BatchMACE", (), {})
    br.relax_many = lambda structures, evaluator, **kw: list(structures)
    pkg, qd = types.ModuleType("orchestr_ai"), types.ModuleType("orchestr_ai.qd")
    pkg.qd, qd.batch_relax = qd, br
    for name, mod in (("orchestr_ai", pkg), ("orchestr_ai.qd", qd), ("orchestr_ai.qd.batch_relax", br)):
        monkeypatch.setitem(sys.modules, name, mod)
    for key in ("d3", "d3_engine", "json", "calc"):
        monkeypatch.setitem(sxl_engine._CONFIG, key, sxl_engine._CONFIG[key])
    return br


def test_configure_is_off_by_default(monkeypatch):
    br = _stub_batch_relax(monkeypatch)
    monkeypatch.delenv("QDW_SXL", raising=False)
    assert sxl_engine.configure() is None
    assert br.BatchMACE is not sxl_engine.SxlBatch


def test_configure_installs_mace_alone(monkeypatch):
    br = _stub_batch_relax(monkeypatch)
    monkeypatch.setenv("QDW_SXL", "1")
    monkeypatch.setenv("QDW_SXL_JSON", "/models/mace_omat.json")
    monkeypatch.delenv("QDW_SXL_D3", raising=False)
    prov = sxl_engine.configure()
    assert br.BatchMACE is sxl_engine.SxlBatch and br.relax_many._sxl
    assert prov["engine"] == "symmetrix-xl" and prov["d3"] is None and prov["json"] == "/models/mace_omat.json"
    monkeypatch.setenv("QDW_SXL_D3", "1")
    br.BatchMACE = type("BatchMACE", (), {})                       # a fresh process, D3 asked for
    assert sxl_engine.configure()["d3"] == "chunked torch-dftd bj/pbe"
