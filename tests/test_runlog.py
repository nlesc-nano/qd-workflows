"""Progress log of a run (qdw log): lines, CP2K and PROPS step results, the PROPS follower, status."""
import json
import threading
import time
from pathlib import Path

from qd_workflows import runlog
from qd_workflows.cli import main

GEO_OUT = """\
 OPTIMIZATION STEP:      1
 *** SCF run converged in    14 steps ***
 ENERGY| Total FORCE_EVAL ( QS ) energy [hartree]           -1234.567890123456
  Max. gradient              =         0.0012345600
 OPTIMIZATION STEP:      2
 *** SCF run converged in     9 steps ***
 ENERGY| Total FORCE_EVAL ( QS ) energy [hartree]           -1234.600000000000
  Max. gradient              =         0.0003000000
 *******************************************************************************
 ***                    GEOMETRY OPTIMIZATION COMPLETED                      ***
 *******************************************************************************
 CP2K                                 1  1.0    0.123    0.456 7384.512 7384.512
  **** **** ******  **  PROGRAM ENDED AT                 2026-10-18 10:00:00.000
"""


def test_line_has_time_stage_job_and_message(tmp_path, monkeypatch):
    log = tmp_path / "pilot.log"
    monkeypatch.setenv("QDW_LOG", str(log))
    monkeypatch.setenv("SLURM_JOB_ID", "123")
    monkeypatch.setenv("SLURM_JOB_NODELIST", "c[01-02]")
    runlog.write("chain", "geo_opt started")
    parts = [p.strip() for p in log.read_text().strip().split("|")]
    assert parts[1:] == ["chain", "123@c[01-02]", "geo_opt started"]
    assert len(parts[0]) == 19                                   # YYYY-mm-dd HH:MM:SS


def test_cp2k_result(tmp_path):
    (tmp_path / "cp2k.out").write_text(GEO_OUT)
    r = runlog.cp2k_result(str(tmp_path))
    assert r == "2 steps converged, E -1234.600000 Ha, max gradient 3.0e-04, 2h03m"
    (tmp_path / "cp2k.out").write_text(GEO_OUT.split(" ****")[0].replace("GEOMETRY OPTIMIZATION COMPLETED", ""))
    assert runlog.cp2k_result(str(tmp_path)).endswith("not converged, E -1234.600000 Ha, max gradient 3.0e-04, not finished")


def test_props_result_and_follow(tmp_path):
    rec = tmp_path / "rec" / "Cd1"
    (rec / "props").mkdir(parents=True)
    (rec / "props" / "relax.json").write_text(json.dumps({"summary": {
        "steps": 111, "fmax_eV_A": 0.0087, "relaxation_energy_meV_atom": -112.7, "rmsd_A": 0.53}}))
    assert runlog.props_result(str(rec), "relax") == "steps 111, fmax 0.009 eV/Å, ΔE -113 meV/atom, RMSD 0.53 Å"
    out, log, stop = tmp_path / "slurm-1.out", tmp_path / "pilot.log", tmp_path / "stop"
    out.write_text("2026 - INFO - [qdprops] Cd1: 381 atoms; steps relax, hessian\n")
    t = threading.Thread(target=runlog.follow, args=(str(out), "gpu", [str(rec)], str(stop)), kwargs={"poll": 0.05})
    import os
    os.environ["QDW_LOG"] = str(log)
    try:
        t.start()
        time.sleep(0.2)
        with open(out, "a") as f:
            f.write("2026 - INFO - [qdprops]   relax: done in 754.2 s\n2026 - INFO - [qdprops]   hessian: cached\n"
                    "[qdprops] Cd1: FAILED ValueError: boom\n")
        stop.touch()
        t.join(5)
    finally:
        del os.environ["QDW_LOG"]
    msgs = [ln.split(" | ")[-1] for ln in log.read_text().strip().split("\n")]
    assert msgs[0] == "Cd1: 381 atoms, steps relax, hessian"
    assert msgs[1] == "relax done in 12m34s: steps 111, fmax 0.009 eV/Å, ΔE -113 meV/atom, RMSD 0.53 Å -> hessian"
    assert msgs[2].startswith("ERROR Cd1: FAILED ValueError: boom (read ")
    assert len(msgs) == 3                                        # cached steps are not logged


def test_cli_and_status(tmp_path, monkeypatch, capsys):
    run = tmp_path / "run"
    c = run / "rec" / "Cd1" / "cp2k"
    (c / "geo_opt").mkdir(parents=True)
    (run / "rec" / "Cd1" / "record.json").write_text("{}")
    (c / "geo_opt" / "cp2k.out").write_text(GEO_OUT)
    (c / "window.json").write_text(json.dumps({"gap_ev": 1.39, "occupied_in_window": 400, "virtual_in_window": 300,
                                               "e_min_ev": -11.0, "e_max_ev": -1.0}))
    log = run / "pilot.log"
    monkeypatch.setenv("QDW_LOG", str(log))
    main(["log", "cp2k", "chain", "geo_opt", str(c / "geo_opt")])
    main(["log", "window", "chain", str(c / "window.json")])
    main(["log", "event", "chain", "field", "done", "--since", str(time.time() - 75)])
    msgs = [ln.split(" | ")[-1] for ln in log.read_text().strip().split("\n")]
    assert msgs == ["geo_opt done: 2 steps converged, E -1234.600000 Ha, max gradient 3.0e-04, 2h03m",
                    "window: KS gap 1.390 eV, 400 occupied + 300 virtual MOs (-11.0 to -1.0 eV)",
                    "field done (1m15s)"]
    main(["status", str(run)])
    text = capsys.readouterr().out
    assert "geo_opt   done: 2 steps" in text and "window    done" in text and "field done (1m15s)" in text


def test_props_jobs_log_and_follow(tmp_path, monkeypatch):
    from qd_workflows.cluster import load_cluster
    from qd_workflows.props import prepare
    rec = tmp_path / "Cd1"
    rec.mkdir()
    (rec / "record.json").write_text("{}")
    (rec / "start.xyz").write_text("1\n\nCd 0 0 0\n")
    monkeypatch.setenv("QDW_LOG", "/scratch/run/pilot.log")
    cl = load_cluster("nhr")
    cl["props"]["chain_time"] = "08:00:00"
    info = prepare([str(rec)], str(tmp_path / "batch"), cl)
    gpu, cpu = Path(info["job"]).read_text(), Path(info["job_cpu"]).read_text()
    for job, stage in ((gpu, "gpu"), (cpu, "cpu")):
        assert "export QDW_LOG=${QDW_LOG:-/scratch/run/pilot.log}" in job
        assert f'qdwlog follow "$OUT" --stage {stage} --records {rec.resolve()}' in job
        assert job.index("trap on_exit EXIT") < job.index("python -m orchestr_ai.postprocessing")
    assert 'qdwlog event gpu "submitted the CP2K chain(s)' in gpu and "--submit --time 08:00:00 |" in gpu
