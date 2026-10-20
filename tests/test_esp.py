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
    # the dots write their own .wfn every 50 cycles (to restart at MAX_SCF); no outer loop
    for n in names:
        tn = (tmp_path / "esp" / n / "cp2k.inp").read_text()
        assert "QS_SCF 50" in tn and "BACKUP_COPIES 0" in tn and "&OUTER_SCF" not in tn



def test_scf_restart_goes_on_from_the_runs_own_wfn(tmp_path):
    d = tmp_path / "dot_12"
    d.mkdir()
    inp = esp.esp_input(BASE, "dot_12", scratch=True, restartable=True)
    (d / "cp2k.inp").write_text(inp)
    (d / "cp2k.out").write_text(" *** SCF run converged in 9 steps ***\n")
    assert esp.scf_restart(d) == "not an SCF failure"
    (d / "cp2k.out").write_text(" SCF run NOT converged. To continue the calculation regardless\n")
    assert esp.scf_restart(d) == "no dot_12-RESTART.wfn"
    (d / "dot_12-RESTART.wfn").write_text("x")
    assert esp.scf_restart(d) == "restarted"
    t = (d / "cp2k.inp").read_text()
    assert "SCF_GUESS  restart" in t and "WFN_RESTART_FILE_NAME dot_12-RESTART.wfn" in t
    # an ESP run accepts the SCF of its restart: EPS_SCF 1e-5 and IGNORE_CONVERGENCE_FAILURE, so the cubes come out
    assert re.search(r"^\s*EPS_SCF\s+1e-5\s*$", t, flags=re.M) and re.search(r"^\s*IGNORE_CONVERGENCE_FAILURE$", t, flags=re.M)
    assert t.count("IGNORE_CONVERGENCE_FAILURE\n") == 1
    assert (d / "cp2k.out.unconverged").is_file() and not (d / "cp2k.out").exists()
    (d / "cp2k.out").write_text(" SCF run NOT converged.\n PROGRAM ENDED AT\n")
    assert not esp._done(d)                                    # no cubes yet
    (d / "v_hartree-1.cube").write_text("x"); (d / "density-1.cube").write_text("x")
    assert esp._done(d) and not esp._scf_converged(d)          # accepted, flagged as not converged
    assert esp.scf_restart(d) == "not an SCF failure"          # no further restart
    for f in ("v_hartree-1.cube", "density-1.cube"):
        (d / f).unlink()
    for k in (2, 3):                                        # up to three restarts, every output kept
        (d / "cp2k.out").write_text(" SCF run NOT converged.\n")
        assert esp.scf_restart(d) == "restarted" and (d / f"cp2k.out.unconverged{k}").is_file()
    assert (d / "cp2k.inp").read_text().count("IGNORE_CONVERGENCE_FAILURE\n") == 1
    (d / "cp2k.out").write_text(" SCF run NOT converged.\n")
    assert esp.scf_restart(d) == "already restarted 3 times"


def test_scf_restart_of_a_field_run_keeps_its_scf_settings(tmp_path):
    d = tmp_path / "fxp"
    d.mkdir()
    inp = esp.restart_print(BASE.replace("X-field", "fxp"))       # a field run: no cubes
    (d / "cp2k.inp").write_text(inp)
    (d / "cp2k.out").write_text(" SCF run NOT converged.\n")
    (d / "fxp-RESTART.wfn").write_text("x")
    assert esp.scf_restart(d) == "restarted"
    t = (d / "cp2k.inp").read_text()
    assert "IGNORE_CONVERGENCE_FAILURE" not in t
    assert re.findall(r"^\s*EPS_SCF\s+(\S+)", t, flags=re.M) == re.findall(r"^\s*EPS_SCF\s+(\S+)", inp, flags=re.M)
