"""qdw cp2k esp: the solvation runs' inputs from the chain's single point."""
import re

from qd_workflows.cp2k import esp

BASE = """&FORCE_EVAL
  &DFT
    CHARGE 0
    &SCF
      EPS_SCF  1e-6
      SCF_GUESS  restart
    &END
    WFN_RESTART_FILE_NAME ../../geo_opt/X-RESTART.wfn
  &END
  &SUBSYS
    &CELL
      ABC 20.000 20.000 20.000
      PERIODIC None
    &END
  &END
&END
&GLOBAL
  PROJECT  X-field
  RUN_TYPE ENERGY_FORCE
&END
"""


def test_esp_input_restart_and_scratch():
    t = esp.esp_input(BASE, "dot_0")
    assert re.search(r"^\s*RUN_TYPE ENERGY$", t, re.M) and "PROJECT  dot_0" in t
    assert t.index("&V_HARTREE_CUBE") < t.index("    &SCF") and "&E_DENSITY_CUBE" in t
    assert "WFN_RESTART_FILE_NAME ../../geo_opt/X-RESTART.wfn" in t and "SCF_GUESS  restart" in t
    s = esp.esp_input(BASE, "MX", scratch=True, cell="15.000 15.000 15.000")
    assert "WFN_RESTART_FILE_NAME" not in s and "SCF_GUESS  atomic" in s and "ABC 15.000 15.000 15.000" in s


def test_esp_input_into_existing_print():
    base = BASE.replace("    &SCF\n", "    &PRINT\n      &MOMENTS\n      &END MOMENTS\n    &END PRINT\n    &SCF\n", 1)
    t = esp.esp_input(base, "dot_0")
    assert t.count("    &PRINT\n") == 1 and t.index("&V_HARTREE_CUBE") < t.index("&MOMENTS")


def test_prepare_stripped_and_added_dots(tmp_path):
    base = tmp_path / "field" / "base"
    base.mkdir(parents=True)
    (base / "cp2k.inp").write_text(BASE)
    (base / "geom.xyz").write_text("2\nintact\nCd 0 0 0\nCl 2.4 0 0\n")
    (tmp_path / "detach_12.xyz").write_text("1\nstripped\nCd 0 0 0\n")
    (tmp_path / "adsorb_12.xyz").write_text("3\nadded\nCd 0 0 0\nCl 2.4 0 0\nCl -2.4 0 0\n")
    names = esp.prepare(str(base / "cp2k.inp"), str(tmp_path / "esp"), stripped=str(tmp_path / "detach_12.xyz"),
                        added=str(tmp_path / "adsorb_12.xyz"))
    assert names == ["dot_0", "dot_12", "dot_add12"]
    t = (tmp_path / "esp" / "dot_add12" / "cp2k.inp").read_text()
    assert "PROJECT  dot_add12" in t and "SCF_GUESS  atomic" in t and "WFN_RESTART_FILE_NAME" not in t
    assert (tmp_path / "esp" / "runs.txt").read_text().split() == names
    # the stripped/added dots (atomic guess) get an OT outer loop with the same EPS_SCF; the intact dot does not
    assert "&OUTER_SCF" in t and "&OUTER_SCF" in (tmp_path / "esp" / "dot_12" / "cp2k.inp").read_text()
    assert "&OUTER_SCF" not in (tmp_path / "esp" / "dot_0" / "cp2k.inp").read_text()
    o = t.split("&OUTER_SCF")[1]
    assert "MAX_SCF 10" in o.split("&END OUTER_SCF")[0] and "EPS_SCF 1e-6" in o.split("&END OUTER_SCF")[0]
    assert t.index("    &SCF") < t.index("&OUTER_SCF") and "MAX_SCF 50" in t
