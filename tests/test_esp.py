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
