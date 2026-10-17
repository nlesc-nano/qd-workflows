"""Campaign ledger and status (qdw campaign): dot phases from files, problems, changes, budget, poll."""
import json
import os
import subprocess
from datetime import datetime, timedelta

import pytest
import yaml

from qd_workflows import budget, campaign
from qd_workflows.cli import main

NOW = datetime(2026, 10, 18, 12, 0, 0)

GEO_RUNNING = """\
 OPTIMIZATION STEP:      1
 ENERGY| Total FORCE_EVAL ( QS ) energy [hartree]           -1234.5
  Max. gradient              =         0.0050000000
 OPTIMIZATION STEP:      2
  Max. gradient              =         0.0020000000
"""
GEO_DONE = GEO_RUNNING + " ***                    GEOMETRY OPTIMIZATION COMPLETED                      ***\n"


def _props(rec, steps, **summaries):
    (rec / "props").mkdir(parents=True, exist_ok=True)
    for i, s in enumerate(steps):
        f = rec / "props" / f"{s}.json"
        f.write_text(json.dumps({"summary": summaries.get(s, {})}))
        os.utime(f, (1_000_000 + i, 1_000_000 + i))


def _ledger(camp, dots, **extra):
    camp.mkdir(exist_ok=True)
    (camp / "campaign.yaml").write_text(yaml.safe_dump({"name": "t", "job_prefixes": ["cd176-"], "dots": dots, **extra}))


def _age(path, hours):
    t = (NOW - timedelta(hours=hours)).timestamp()
    os.utime(path, (t, t))


@pytest.fixture
def camp(tmp_path):
    gpu = tmp_path / "rec" / "A"                      # GPU stage half way
    _props(gpu, ["relax", "structure", "hessian"], hessian={"n_imaginary": 2, "lowest_cm1": -30.0})
    chain = tmp_path / "rec" / "B"                    # chain in geo_opt
    _props(chain, ["relax", "structure", "hessian", "stability", "wigner"])
    (chain / "cp2k" / "geo_opt").mkdir(parents=True)
    (chain / "cp2k" / "geo_opt" / "cp2k.out").write_text(GEO_RUNNING)
    cpu = tmp_path / "rec" / "C"                      # chain finished, CPU stage to do
    _props(cpu, ["relax", "hessian", "field"])
    (cpu / "cp2k" / "geo_opt").mkdir(parents=True)
    (cpu / "cp2k" / "geo_opt" / "cp2k.out").write_text(GEO_DONE)
    (cpu / "cp2k" / "qdex").mkdir()
    (cpu / "cp2k" / ".done").write_text("")
    done = tmp_path / "rec" / "D"
    _props(done, ["relax", "hessian", "report"])
    log = tmp_path / "pilot.log"
    log.write_text("2026-10-18 10:00:00 | gpu   | 1@x | start\n")
    _ledger(tmp_path / "camp", [
        {"id": "A", "record": str(gpu), "log": str(log)},
        {"id": "B", "record": str(chain), "log": str(log)},
        {"id": "C", "record": str(cpu)},
        {"id": "D", "record": str(done)},
        {"id": "L", "source": "legacy", "legacy": {"geo": str(tmp_path / "old.xyz")}},
        {"id": "N", "source": "qd_builder"},
    ])
    return tmp_path


def _st(camp_root, jobs=(), acct=(), update=True, cluster=None):
    return campaign.status(camp_root / "camp", update=update, now=NOW, jobs=list(jobs), acct=list(acct),
                           use_slurm=False, cluster=cluster or {})


def test_phases_and_steps(camp):
    st = _st(camp, update=False)
    d = {x["id"]: x for x in st["dots"]}
    assert (d["A"]["phase"], d["A"]["step"]) == ("gpu", "hessian")
    assert (d["B"]["phase"], d["B"]["step"]) == ("chain", "geo_opt")
    assert d["B"]["chain"]["geo_opt"] == {"steps": 2, "max_gradient": 0.002, "converged": False,
                                          "last_write_h": d["B"]["chain"]["geo_opt"]["last_write_h"]}
    assert (d["C"]["phase"], d["C"]["step"]) == ("cpu", "vibspec")
    assert d["D"]["phase"] == "done"
    assert d["L"]["phase"] == "legacy" and d["L"]["legacy"] == {"geo": False}
    assert d["N"]["phase"] == "new"
    assert st["totals"]["dots"] == 6 and st["totals"]["chain"] == 1


def test_problems(camp):
    st = _st(camp, update=False)
    keys = {(a["dot"], a["key"]) for a in st["attention"]}
    assert ("A", "imaginary") in keys
    assert ("C", "no_webapp") in keys                 # chain .done but no qdex/webapp.json
    assert ("L", "legacy_missing") in keys


def test_stalled_chain_only_without_a_job(camp):
    _age(camp / "rec" / "B" / "cp2k" / "geo_opt" / "cp2k.out", 5)
    keys = {(a["dot"], a["key"]) for a in _st(camp, update=False)["attention"]}
    assert ("B", "stalled") in keys
    job = {"id": "9", "name": "cp2k-B", "partition": "p", "state": "RUNNING", "reason": "", "waiting_h": None}
    keys = {(a["dot"], a["key"]) for a in _st(camp, jobs=[job], update=False)["attention"]}
    assert ("B", "stalled") not in keys


def test_changes_between_snapshots(camp):
    first = _st(camp)
    assert "changes" not in first                     # nothing to compare with yet
    assert (camp / "camp" / "status" / "last.json").is_file()
    # the chain moves on, a job ends badly, the log gets an error
    (camp / "rec" / "B" / "cp2k" / "geo_opt" / "cp2k.out").write_text(GEO_DONE)
    (camp / "rec" / "B" / "cp2k" / "field" / "x").mkdir(parents=True)
    with open(camp / "pilot.log", "a") as f:
        f.write("2026-10-18 11:00:00 | chain | 2@y | ERROR CP2K failed (read z)\n")
        f.write("2026-10-18 11:00:01 | chain | 2@y | geo_opt started\n")
    ended = {"id": "77", "user": os.environ.get("USER", "u"), "name": "cp2k-B", "partition": "cpu", "state": "TIMEOUT",
             "start": NOW - timedelta(hours=2), "end": NOW, "elapsed_s": 7200, "cpus": 96, "gpus": 0, "nodes": 1}
    import getpass
    ended["user"] = getpass.getuser()
    second = _st(camp, acct=[ended])
    kinds = {(c["kind"], c.get("dot") or c.get("job")) for c in second["changes"]}
    assert ("progress", "B") in kinds
    assert ("job_ended", "77") in kinds
    assert [c["line"] for c in second["changes"] if c["kind"] == "log"] == [
        "2026-10-18 11:00:00 | chain | 2@y | ERROR CP2K failed (read z)"]
    third = _st(camp, acct=[ended])                   # the same job is not reported twice
    assert not [c for c in third["changes"] if c["kind"] == "job_ended"]
    events = (camp / "camp" / "status" / "events.jsonl").read_text().splitlines()
    assert len(events) == len(second["changes"])


def test_pending_job_of_the_batch(camp):
    job = {"id": "5", "name": "cd176-gpu-cpu", "partition": "cpu-genoa", "state": "PENDING", "reason": "Priority",
           "waiting_h": 30.0}
    other = dict(job, id="6", name="usman_md")      # not the campaign's
    st = _st(camp, jobs=[job, other], update=False)
    assert [j["id"] for j in st["jobs"]] == ["5"]
    assert any(a["dot"] is None and a["key"] == "pending:5" for a in st["attention"])


def test_add_refuses_duplicates(tmp_path):
    campaign.add(tmp_path, {"id": "X", "material": "InP", "record": None})
    assert campaign.load(tmp_path)["dots"] == [{"id": "X", "material": "InP"}]
    with pytest.raises(ValueError):
        campaign.add(tmp_path, {"id": "X"})


# ---------------- budget ----------------
NHR = {"name": "nhr", "threads_per_core": 2, "budget": {
    "account": "bep00157", "period_start": "2026-10-01", "period_end": "2026-12-31", "allocation_core_h": 10000,
    "charge": {"cpu-genoa": {"per_node_h": 192}, "gpu-a100:shared": {"per_gpu_h": 150}, "gpu-a100": {"per_node_h": 600}}}}


def _job(i, part, hours, nodes=1, cpus=384, gpus=0, user="me", name="cp2k-B", ago_h=None):
    end = NOW - timedelta(hours=ago_h) if ago_h is not None else NOW
    return {"id": str(i), "user": user, "name": name, "partition": part, "state": "COMPLETED",
            "start": end - timedelta(hours=hours), "end": end, "elapsed_s": int(hours * 3600),
            "cpus": cpus, "gpus": gpus, "nodes": nodes}


def test_charge_rates():
    ch, tpc = NHR["budget"]["charge"], 2
    assert budget.job_core_h(_job(1, "cpu-genoa", 2, nodes=2), ch, tpc, 2) == 768
    assert budget.job_core_h(_job(2, "gpu-a100:shared", 1, cpus=36, gpus=1), ch, tpc, 1) == 150   # longest prefix
    assert budget.job_core_h(_job(3, "gpu-a100", 1, gpus=4), ch, tpc, 1) == 600
    assert budget.job_core_h(_job(4, "cpu-clx", 1, cpus=192), ch, tpc, 1) == 96    # no rate: cores x hours


def test_budget_totals_and_warnings():
    jobs = [_job(1, "cpu-genoa", 40, ago_h=48), _job(2, "cpu-genoa", 2, user="other", name="x"),
            _job(3, "gpu-a100:shared", 2, cpus=36, gpus=1)]
    b = budget.budget(NHR, now=NOW, jobs=jobs, belongs=lambda n: n.startswith("cp2k-"), user="me")
    assert b["used_core_h"] == 192 * 40 + 192 * 2 + 300
    assert b["mine_core_h"] == 192 * 40 + 300 and b["campaign_core_h"] == 192 * 40 + 300
    assert b["last_24h_core_h"] == 192 * 2 + 300       # the 40 h job ended two days ago
    assert b["fraction_used"] == pytest.approx(0.836, abs=1e-3)
    assert {p["key"] for p in b["problems"]} == {"budget_75", "budget_rate"}


def test_no_allocation_no_fraction():
    prof = {"budget": {"account": "bcmaterials", "period_start": "2026-10-01", "allocation_core_h": None, "charge": {}}}
    b = budget.budget(prof, now=NOW, jobs=[_job(1, "general", 10, cpus=48)], user="me")
    assert b["used_core_h"] == 480 and "fraction_used" not in b and not b["problems"]


def test_status_carries_the_budget(camp):
    import getpass
    jobs = [_job(1, "cpu-genoa", 1, user=getpass.getuser(), name="cp2k-B")]
    st = _st(camp, acct=jobs, cluster=NHR, update=False)
    assert st["budget"]["used_core_h"] == 192
    assert st["totals"]["campaign_core_h"] == 192
    assert "budget bep00157" in campaign.text(st)


# ---------------- poll ----------------
def test_poll_collects_every_cluster(tmp_path, monkeypatch):
    d = tmp_path / "clusters"
    d.mkdir()
    (d / "a.yaml").write_text(yaml.safe_dump({"name": "a", "campaign": {"host": "ha", "path": "/c", "python": "py"}}))
    (d / "b.yaml").write_text(yaml.safe_dump({"name": "b", "campaign": {"label": "bee", "host": "hb", "path": "/c",
                                                                         "python": "py", "pythonpath": "/src"}}))
    (d / "c.yaml").write_text(yaml.safe_dump({"name": "c"}))            # no campaign: left out
    monkeypatch.setenv("QDW_CLUSTERS", str(d))
    calls = []

    def fake_run(cmd, **kw):
        calls.append(cmd)
        if cmd[-2] == "hb":
            return subprocess.CompletedProcess(cmd, 255, "", "ssh: connect to host hb: timed out\n")
        return subprocess.CompletedProcess(cmd, 0, json.dumps({"attention": [{"dot": "X", "level": "warn", "key": "k",
                                                                               "what": "w"}], "changes": []}), "")
    monkeypatch.setattr(campaign.subprocess, "run", fake_run)
    res = campaign.poll(["a", "b"])
    assert set(res["clusters"]) == {"a", "bee"}
    assert "error" in res["clusters"]["bee"]
    assert {(a["cluster"], a["key"]) for a in res["attention"]} == {("a", "k"), ("bee", "unreachable")}
    assert calls[1][-1].startswith("PYTHONPATH=/src py -m qd_workflows.cli campaign status /c --json --update")


def test_cli_status_text(camp, capsys, monkeypatch):
    monkeypatch.setattr(campaign, "queue", lambda now: [])
    monkeypatch.setattr(budget, "sacct_jobs", lambda a, s: [])
    assert main(["campaign", "status", str(camp / "camp")]) == 0
    out = capsys.readouterr().out
    assert "needs attention" in out and "B" in out


def _xyz(path, atoms):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f"{len(atoms)}\nx\n" + "\n".join(f"{s} {i * 2.5} {i * 0.5} {-i * 0.3}" for i, s in enumerate(atoms)) + "\n")


def test_legacy_records(tmp_path, capsys):
    from qd_workflows import legacy
    camp = tmp_path / "camp"
    atoms_ok = ["In"] * 3 + ["P"] * 2 + ["Cl"] * 3           # In3P2Cl3: charge 9 - 6 - 3 = 0
    _xyz(tmp_path / "geo" / "a.xyz", ["In"] * 2 + atoms_ok)    # an earlier frame with another stoichiometry...
    with (tmp_path / "geo" / "a.xyz").open("a") as f:           # ...then the final frame
        f.write(f"{len(atoms_ok)}\nfinal\n" + "\n".join(f"{s} {i} 0 0" for i, s in enumerate(atoms_ok)) + "\n")
    _xyz(tmp_path / "trim" / "InP-P-In3P2Cl3-reconstructed" / "InP-P-In3P2Cl3-opt.xyz", atoms_ok)
    _xyz(tmp_path / "trim" / "InP-In-In3P2Cl4-reconstructed" / "InP-In-In3P2Cl4-opt.xyz", atoms_ok)
    led = {"name": "t", "dots": [
        {"id": "InP-In-In3P2Cl3", "material": "InP", "source": "legacy", "surface": "reconstructed",
         "reconstruction": "vacancy", "legacy": {"geo": str(tmp_path / "geo" / "a.xyz")}},
        {"id": "InP-P-In3P2Cl3-reconstructed", "material": "InP", "source": "legacy", "legacy": {"qdex": "/none"}},
        {"id": "InP-In-In3P2Cl4-reconstructed", "material": "InP", "source": "legacy"},        # atoms do not match
        {"id": "InP-In-In9P9Cl9-reconstructed", "material": "InP", "source": "legacy"},        # no geometry
        {"id": "InP-In-In3P2Cl3-big", "material": "InP", "source": "legacy", "legacy": {"geo": str(tmp_path / "geo" / "a.xyz")}},
        {"id": "CdSe-pilot", "material": "CdSe", "source": "pilot"}]}
    camp.mkdir()
    (camp / "campaign.yaml").write_text(yaml.safe_dump(led))
    res = legacy.make_records(camp, tmp_path / "rec", [str(tmp_path / "trim" / "{id}" / "{stem}-opt.xyz")],
                              max_atoms=8, update_ledger=True)
    assert res["made"] == ["InP-In-In3P2Cl3", "InP-P-In3P2Cl3-reconstructed", "InP-In-In3P2Cl3-big"]
    assert "id says 9" in res["missing"]["InP-In-In3P2Cl4-reconstructed"]
    assert res["missing"]["InP-In-In9P9Cl9-reconstructed"] == "no geometry found"
    rec = json.loads((tmp_path / "rec" / "InP-In-In3P2Cl3" / "record.json").read_text())
    assert rec["n_atoms"] == 8 and rec["total_charge"] == 0 and rec["family"] == "III-V"
    assert rec["composition"] == {"In": 3, "P": 2, "Cl": 3} and rec["ligands"] == {"Cl": 3} and rec["reconstruction"] == "vacancy"
    assert (tmp_path / "rec" / "InP-In-In3P2Cl3" / "start.xyz").read_text().startswith("8\n")      # last frame, not the first
    dots = {d["id"]: d for d in campaign.load(camp)["dots"]}
    assert dots["InP-In-In3P2Cl3"]["record"] == str(tmp_path / "rec" / "InP-In-In3P2Cl3")
    assert "record" not in dots["CdSe-pilot"]
    again = legacy.make_records(camp, tmp_path / "rec", [], max_atoms=8)           # a second call leaves them alone
    assert again["made"] == [] and again["skipped"]["InP-In-In3P2Cl3"] == "has a record"
    assert main(["campaign", "legacy-records", str(camp), "--out", str(tmp_path / "rec2"), "--max-atoms", "8"]) == 0
    assert json.loads(capsys.readouterr().out)["made"] == []


def test_legacy_records_strip_extra(tmp_path):
    from qd_workflows import legacy
    camp = tmp_path / "camp"
    camp.mkdir()
    atoms = ["In"] * 3 + ["P"] * 2 + ["Cl"] * 3 + ["N", "C", "H", "H", "H"]
    _xyz(tmp_path / "g.xyz", atoms)
    (camp / "campaign.yaml").write_text(yaml.safe_dump({"name": "t", "dots": [
        {"id": "InP-In-In3P2Cl3", "material": "InP", "source": "legacy", "legacy": {"geo": str(tmp_path / "g.xyz")}}]}))
    res = legacy.make_records(camp, tmp_path / "rec")
    assert "id says 8" in res["missing"]["InP-In-In3P2Cl3"]                       # without the option: not made
    res = legacy.make_records(camp, tmp_path / "rec", strip_extra=True)
    assert res["made"] == ["InP-In-In3P2Cl3"]
    rec = json.loads((tmp_path / "rec" / "InP-In-In3P2Cl3" / "record.json").read_text())
    assert rec["n_atoms"] == 8 and rec["total_charge"] == 0
    assert rec["origin"]["stripped_atoms"] == {"N": 1, "C": 1, "H": 3}
    assert (tmp_path / "rec" / "InP-In-In3P2Cl3" / "start.xyz").read_text().startswith("8\n")
