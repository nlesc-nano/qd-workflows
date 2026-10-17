"""qdw qdex prepare: the QDEX config of the chain's two runs (all transitions, traps filtered)."""
import json

import yaml

from qd_workflows import qdex

CLUSTER = {"cp2k": {"data_dir": "/data", "basis_file": "BASIS_MOLOPT_UZH", "basis": "DZVP-MOLOPT-GGA-GTH"},
           "qdex": {"gth_soc_file": "/data/GTH_SOC_POTENTIALS"}}


def _prepare(tmp_path, out, **kw):
    record, window = tmp_path / "record.json", tmp_path / "window.json"
    record.write_text(json.dumps({"material": "CdSe", "composition": {"Cd": 16, "Se": 13, "Cl": 6}}))
    window.write_text(json.dumps({"occupied_in_window": 120, "virtual_in_window": 90}))
    qdex.prepare(str(record), str(window), str(tmp_path / out), CLUSTER, **kw)
    return yaml.safe_load((tmp_path / out / "config.yaml").read_text())


def test_trap_filter_only_in_the_trap_run(tmp_path):
    full, trap = _prepare(tmp_path, "qdex"), _prepare(tmp_path, "qdex_trap", trap_filter=True)
    assert "trap_filter" not in full["excitations"]
    assert trap["excitations"].pop("trap_filter") is True
    assert trap == full                                  # otherwise the same run
    assert full["excitations"]["nhomos"] == full["excitations"]["nlumos"] == 90
    assert full["output"]["h5"] and not full["output"]["html"]     # dashboards come from the HDF5 files
    assert full["system"]["mo_file"] == "../orbitals.h5"           # both runs read the chain's orbitals
