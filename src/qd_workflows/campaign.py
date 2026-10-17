"""
A campaign: the dots of the structure database and where each one's files are.

The ledger, <campaign>/campaign.yaml, lists the dots; nobody writes their state by hand.
`qdw campaign status` reads it from each record's files (props/*.json, the CP2K chain), the queue,
the accounting and the run logs, and with --update keeps a snapshot (<campaign>/status/last.json)
so that each call also says what changed since the previous one (appended to status/events.jsonl).
The watchdog and the agents read `qdw campaign status --json --update`; people read the text form.

    name: qd-database
    cluster: nhr                   # profile (clusters/nhr.yaml): account, charge rates, allocation
    started: 2026-10-17            # jobs ending from this day on are reported
    job_prefixes: [cd176-, cp2k-]  # Slurm job names that belong to the campaign
    dots:
      - id: CdSe-Se-Cd176Se147Cl58-clean
        material: CdSe
        source: pilot              # qd_builder | legacy | pilot
        surface: clean             # clean | reconstructed
        reconstruction: null       # reconstructed only: vacancy | swap (QD_builder reconstruction_pattern)
        record: /scratch/.../rec/CdSe-Se-Cd176Se147Cl58-clean
        log: /scratch/.../pilot.log
      - id: InP-In-In149P108Cl123
        material: InP
        source: legacy             # DFT done outside the workflow: to be imported
        legacy: {geo: /scratch/.../pbe_final.xyz, qdex: /scratch/.../qdex/In-In149P108Cl123}
"""
from __future__ import annotations

import json
import re
import subprocess
from datetime import datetime
from pathlib import Path

import yaml

import getpass

from . import budget as budget_mod
from . import runlog
from .cluster import load_cluster

PHASES = ["legacy", "new", "gpu", "chain", "cpu", "done"]
CHAIN_STEPS = ["geo_opt", "field", "esp", "pdos", "window", "trexio", "trim", "qdex"]
STALL_HOURS = 3.0          # a running step whose output has not changed for this long, with no job, is stalled
PENDING_HOURS = 24.0       # a job pending longer than this is worth a look
GEO_MAX_STEPS = 150        # geo_opt steps without convergence before it is worth a look


# ---------------- ledger ----------------
def ledger_path(campaign: str | Path) -> Path:
    return Path(campaign) / "campaign.yaml"


def load(campaign: str | Path) -> dict:
    led = yaml.safe_load(ledger_path(campaign).read_text()) or {}
    led.setdefault("dots", [])
    ids = [d["id"] for d in led["dots"]]
    dup = {i for i in ids if ids.count(i) > 1}
    if dup:
        raise ValueError(f"dots listed twice in {ledger_path(campaign)}: {', '.join(sorted(dup))}")
    return led


def add(campaign: str | Path, entry: dict) -> dict:
    """Append one dot to the ledger (refuses a duplicate id)."""
    p = ledger_path(campaign)
    led = load(campaign) if p.is_file() else {"name": Path(campaign).name, "dots": []}
    if any(d["id"] == entry["id"] for d in led["dots"]):
        raise ValueError(f"{entry['id']} is already in {p}")
    led["dots"].append({k: v for k, v in entry.items() if v is not None})
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(yaml.safe_dump(led, sort_keys=False, allow_unicode=True))
    return led


# ---------------- Slurm ----------------
def _run(cmd: list[str]) -> str | None:
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.TimeoutExpired):
        return None
    return r.stdout if r.returncode == 0 else None


def _time(s: str) -> datetime | None:
    try:
        return datetime.fromisoformat(s)
    except (TypeError, ValueError):
        return None


def queue(now: datetime) -> list[dict] | None:
    """The user's jobs in the queue (None when squeue cannot be run)."""
    out = _run(["squeue", "--me", "-h", "-o", "%i|%j|%P|%T|%M|%l|%D|%R|%V"])
    if out is None:
        return None
    jobs = []
    for ln in out.splitlines():
        f = ln.split("|")
        if len(f) < 9:
            continue
        submit = _time(f[8])
        jobs.append({"id": f[0], "name": f[1], "partition": f[2], "state": f[3], "elapsed": f[4], "limit": f[5],
                     "nodes": int(f[6]) if f[6].isdigit() else f[6], "reason": f[7].strip("()"),
                     "waiting_h": round((now - submit).total_seconds() / 3600, 1) if submit and f[3] == "PENDING" else None})
    return jobs


def _belongs(name: str, led: dict) -> bool:
    if any(name.startswith(p) for p in led.get("job_prefixes", [])):
        return True
    return any(d["id"] in name for d in led["dots"])


def _jobs_of(dot_id: str, jobs: list[dict]) -> list[dict]:
    return [j for j in jobs if dot_id in j["name"]]


# ---------------- one dot ----------------
def _summary(path: Path) -> dict:
    try:
        return json.loads(path.read_text()).get("summary") or {}
    except (OSError, ValueError, AttributeError):
        return {}


def _age_h(path: Path, now: datetime) -> float:
    return round((now - datetime.fromtimestamp(path.stat().st_mtime)).total_seconds() / 3600, 1)


def chain_state(cp2k: Path, now: datetime) -> dict:
    """Each CP2K chain step done or not, the first one not done, and what is wrong with them."""
    st: dict = {"steps": {}, "problems": []}
    g = cp2k / "geo_opt" / "cp2k.out"
    if g.is_file():
        t = g.read_text(errors="replace")
        steps, grad = runlog._OPT.findall(t), runlog._GRAD.findall(t)
        done = "GEOMETRY OPTIMIZATION COMPLETED" in t
        st["geo_opt"] = {"steps": int(steps[-1]) if steps else 0, "max_gradient": float(grad[-1]) if grad else None,
                         "converged": done, "last_write_h": _age_h(g, now)}
        st["steps"]["geo_opt"] = done
        if not done and steps and int(steps[-1]) > GEO_MAX_STEPS:
            st["problems"].append({"level": "warn", "key": "geo_opt_long",
                                   "what": f"geo_opt at step {steps[-1]} and not converged (max gradient {grad[-1] if grad else '?'})"})
    else:
        st["steps"]["geo_opt"] = False
    for name in ("field", "esp"):
        n_done, n = runlog._count_done(cp2k / name)
        st["steps"][name] = bool(n) and n_done == n
        if n:
            st[name] = {"finished": n_done, "runs": n}
    st["steps"]["pdos"] = any(cp2k.glob("pdos/*.pdos"))
    st["steps"]["window"] = (cp2k / "window.json").is_file()
    st["steps"]["trexio"] = any(cp2k.glob("trexio/orbitals*"))
    st["steps"]["trim"] = (cp2k / "orbitals.h5").is_file()
    st["steps"]["qdex"] = (cp2k / "qdex" / "webapp.json").is_file()
    if st["steps"]["qdex"]:
        st["qdex"] = runlog.qdex_result(str(cp2k / "qdex"))
    st["finished"] = (cp2k / ".done").is_file()
    st["next"] = next((s for s in CHAIN_STEPS if not st["steps"][s]), None)
    if st["finished"] and not st["steps"]["qdex"]:
        st["problems"].append({"level": "error", "key": "no_webapp", "what": "chain finished without qdex/webapp.json"})
    for out in sorted(cp2k.glob("*/cp2k.out")) + sorted(cp2k.glob("*/*/cp2k.out")):
        tail = out.read_text(errors="replace")[-20000:]
        if "ABORT" in tail:
            st["problems"].append({"level": "error", "key": f"abort:{out.parent.relative_to(cp2k)}",
                                   "what": f"CP2K aborted in {out.parent.relative_to(cp2k)} (read {out})"})
    return st


def props_state(props: Path) -> dict:
    """The GPU / CPU stage steps written so far and the warnings their summaries give."""
    files = sorted(props.glob("*.json"), key=lambda p: p.stat().st_mtime) if props.is_dir() else []
    steps = [p.stem for p in files]                     # in the order they were written
    st: dict = {"steps": steps, "problems": []}
    h = _summary(props / "hessian.json")
    if h.get("n_imaginary"):
        st["problems"].append({"level": "warn", "key": "imaginary",
                               "what": f"{h['n_imaginary']} imaginary modes (lowest {h.get('lowest_cm1')} cm-1)"})
    r = _summary(props / "relax.json")
    if (r.get("fmax_eV_A") or 0) > 0.05:
        st["problems"].append({"level": "warn", "key": "relax_fmax", "what": f"MACE relax ended at fmax {r['fmax_eV_A']:.3f} eV/Å"})
    s = _summary(props / "structure.json")
    if s.get("ligands_detached"):
        st["problems"].append({"level": "warn", "key": "ligands_detached",
                               "what": f"{s['ligands_detached']} ligands detached in the MACE relax"})
    return st


def dot_state(d: dict, jobs: list[dict], now: datetime) -> dict:
    out: dict = {"id": d["id"], "material": d.get("material"), "source": d.get("source"),
                 "surface": d.get("surface"), "reconstruction": d.get("reconstruction"), "problems": []}
    rec = Path(d["record"]) if d.get("record") else None
    if rec is None or not rec.is_dir():
        if d.get("source") == "legacy":
            leg = d.get("legacy") or {}
            out["phase"] = "legacy"
            out["legacy"] = {k: Path(v).exists() for k, v in leg.items()}
            missing = [k for k, ok in out["legacy"].items() if not ok]
            if missing:
                out["problems"].append({"level": "warn", "key": "legacy_missing", "what": f"legacy {', '.join(missing)} not found"})
        else:
            out["phase"] = "new"
            if rec is not None:
                out["problems"].append({"level": "error", "key": "no_record", "what": f"record directory {rec} does not exist"})
        out["jobs"] = [j["id"] for j in _jobs_of(d["id"], jobs)]
        return out

    props = props_state(rec / "props")
    cp2k = rec / "cp2k"
    chain = chain_state(cp2k, now) if cp2k.is_dir() else None
    out["props"] = props["steps"]
    out["problems"] += props["problems"]
    if chain:
        out["chain"] = {k: v for k, v in chain.items() if k != "problems"}
        out["problems"] += chain["problems"]
    if "report" in props["steps"]:
        out["phase"], out["step"] = "done", None
    elif chain and chain["finished"]:
        out["phase"], out["step"] = "cpu", next((s for s in ("field", "vibspec", "solvation", "sites", "report")
                                                 if s not in props["steps"]), "report")
    elif chain:
        out["phase"], out["step"] = "chain", chain["next"]
    elif props["steps"]:
        out["phase"], out["step"] = "gpu", props["steps"][-1]
    else:
        out["phase"], out["step"] = "new", None

    mine = _jobs_of(d["id"], jobs)
    out["jobs"] = [j["id"] for j in mine]
    running = any(j["state"] == "RUNNING" for j in mine)
    if chain and out["phase"] == "chain" and not running:
        g = chain.get("geo_opt") or {}
        newest = max((p.stat().st_mtime for p in cp2k.rglob("cp2k.out")), default=None)
        idle_h = round((now.timestamp() - newest) / 3600, 1) if newest else None
        if idle_h is not None and idle_h > STALL_HOURS and not any(j["state"] == "PENDING" for j in mine):
            out["problems"].append({"level": "error", "key": "stalled",
                                    "what": f"chain at {chain['next']}, no job in the queue, output unchanged for {idle_h} h"
                                            + (f" (geo_opt step {g.get('steps')})" if g and not g.get("converged") else "")})
    return out


# ---------------- logs ----------------
def log_lines(path: Path, offset: int) -> tuple[list[str], int]:
    """Lines added to a log since byte `offset` (the whole log when it shrank) and the new offset."""
    if not path.is_file():
        return [], 0
    size = path.stat().st_size
    if size < offset:
        offset = 0
    with open(path, "rb") as f:
        f.seek(offset)
        text = f.read().decode(errors="replace")
    return [ln for ln in text.splitlines() if ln.strip()], size


# ---------------- status ----------------
def status(campaign: str | Path, update: bool = False, now: datetime | None = None,
           jobs: list[dict] | None = None, acct: list[dict] | None = None, use_slurm: bool = True,
           cluster: dict | None = None) -> dict:
    """The state of every dot, the campaign's jobs, the budget, what needs attention and what changed."""
    camp = Path(campaign)
    led = load(camp)
    now = now or datetime.now()
    if cluster is None and led.get("cluster"):
        cluster = load_cluster(led["cluster"])
    cluster = cluster or {}
    account = (cluster.get("budget") or {}).get("account") or cluster.get("account")
    if use_slurm and jobs is None:
        jobs = queue(now)
    if use_slurm and acct is None:
        start = budget_mod._date(led.get("started")) or now.date()
        if cluster.get("budget"):
            start = min(start, budget_mod._date(cluster["budget"].get("period_start")) or start)
        acct = budget_mod.sacct_jobs(account, start)
    slurm_ok = jobs is not None
    belongs = lambda name: _belongs(name, led)  # noqa: E731
    me = getpass.getuser()
    all_acct = acct or []
    jobs = [j for j in (jobs or []) if belongs(j["name"])]
    started = budget_mod._date(led.get("started"))
    acct = [j for j in all_acct if belongs(j["name"]) and j.get("user", me) == me
            and (started is None or (j.get("start") and j["start"].date() >= started))]
    for j in acct:
        j["core_h"] = round(budget_mod.job_core_h(j, (cluster.get("budget") or {}).get("charge") or {},
                                                  int(cluster.get("threads_per_core", 1)), j["elapsed_s"] / 3600), 1)
    bud = budget_mod.budget(cluster, now=now, jobs=all_acct if use_slurm or all_acct else None,
                            belongs=belongs, user=me) if cluster.get("budget") and (use_slurm or all_acct) else None

    dots = [dot_state(d, jobs, now) for d in led["dots"]]
    general = []                       # problems of no single dot (batch jobs)
    for j in jobs:
        if j["state"] == "PENDING" and (j.get("waiting_h") or 0) > PENDING_HOURS:
            p = {"level": "warn", "key": f"pending:{j['id']}",
                 "what": f"job {j['id']} ({j['name']}) pending {j['waiting_h']} h ({j['reason']}) on {j['partition']}"}
            mine = [d for d in dots if j["id"] in d.get("jobs", [])]
            for d in mine:
                d["problems"].append(p)
            if not mine:
                general.append(p)
    general += (bud or {}).get("problems", [])

    snap_dir = camp / "status"
    last_file = snap_dir / "last.json"
    prev = json.loads(last_file.read_text()) if last_file.is_file() else None
    prev_offsets = (prev or {}).get("log_offsets", {})

    # new log lines (ERROR / WARN / done) since the last snapshot
    offsets, log_events = {}, []
    for lg in sorted({d["log"] for d in led["dots"] if d.get("log")}):
        new, offsets[lg] = log_lines(Path(lg), prev_offsets.get(lg, 0) if prev else 0)
        if prev:
            log_events += [{"log": lg, "line": ln} for ln in new
                           if re.search(r"\| (ERROR|WARN)\b|\bERROR\b|run ended|chain finished|end:", ln)]

    by_phase = {p: sum(d["phase"] == p for d in dots) for p in PHASES}
    attention = [{"dot": None, **p} for p in general] + [{"dot": d["id"], **p} for d in dots for p in d["problems"]]
    res = {
        "campaign": led.get("name", camp.name),
        "path": str(camp),
        "time": now.isoformat(timespec="seconds"),
        "slurm": slurm_ok,
        "totals": {"dots": len(dots), **by_phase,
                   "running_jobs": sum(j["state"] == "RUNNING" for j in jobs),
                   "pending_jobs": sum(j["state"] == "PENDING" for j in jobs),
                   "campaign_core_h": round(sum(j["core_h"] for j in acct))},
        "budget": bud,
        "attention": attention,
        "general": general,
        "jobs": jobs,
        "dots": dots,
        "log_tail": {lg: (log_lines(Path(lg), 0)[0][-3:]) for lg in offsets},
    }
    if prev is not None:
        res["since"] = prev["time"]
        res["changes"] = changes(prev, res, acct) + [{"kind": "log", **e} for e in log_events]
    if update:
        snap_dir.mkdir(parents=True, exist_ok=True)
        snap = {k: res[k] for k in ("time", "totals", "attention", "jobs", "dots")}
        snap["general"] = general
        snap["log_offsets"] = offsets
        snap["ended_jobs"] = sorted({j["id"] for j in acct if j["end"]} | set((prev or {}).get("ended_jobs", [])))
        snap["budget"] = bud
        last_file.write_text(json.dumps(snap, indent=1))
        if res.get("changes"):
            with open(snap_dir / "events.jsonl", "a") as f:
                for c in res["changes"]:
                    f.write(json.dumps({"time": res["time"], **c}) + "\n")
    return res


def changes(prev: dict, cur: dict, acct: list[dict]) -> list[dict]:
    """What moved between two status snapshots: dot phases and steps, problems, jobs started and ended."""
    out = []
    before = {d["id"]: d for d in prev.get("dots", [])}
    for d in cur["dots"]:
        b = before.get(d["id"])
        if b is None:
            out.append({"kind": "dot_added", "dot": d["id"], "phase": d["phase"]})
            continue
        if (b.get("phase"), b.get("step")) != (d.get("phase"), d.get("step")):
            out.append({"kind": "progress", "dot": d["id"], "from": f"{b.get('phase')} {b.get('step') or ''}".strip(),
                        "to": f"{d.get('phase')} {d.get('step') or ''}".strip()})
        keys_b = {p["key"] for p in b.get("problems", [])}
        keys_d = {p["key"] for p in d.get("problems", [])}
        out += [{"kind": "problem", "dot": d["id"], **p} for p in d["problems"] if p["key"] not in keys_b]
        out += [{"kind": "resolved", "dot": d["id"], "key": k} for k in sorted(keys_b - keys_d)]
    keys_b = {p["key"] for p in prev.get("general", [])}
    out += [{"kind": "problem", "dot": None, **p} for p in cur.get("general", []) if p["key"] not in keys_b]
    prev_state = {j["id"]: j["state"] for j in prev.get("jobs", [])}
    for j in cur["jobs"]:
        if j["state"] == "RUNNING" and prev_state.get(j["id"]) != "RUNNING":
            out.append({"kind": "job_started", "job": j["id"], "name": j["name"], "partition": j["partition"]})
    seen = set(prev.get("ended_jobs", []))
    for j in acct:
        if j["end"] and j["id"] not in seen:
            out.append({"kind": "job_ended", "job": j["id"], "name": j["name"], "state": j["state"],
                        "core_h": j["core_h"], "level": "info" if j["state"] in ("COMPLETED", "CANCELLED") else "error"})
    bp, bc = prev.get("budget") or {}, cur.get("budget") or {}
    if bc.get("used_core_h") is not None:
        out.append({"kind": "budget", "used_core_h": bc["used_core_h"],
                    "since_last": bc["used_core_h"] - (bp.get("used_core_h") or 0) if bp else None,
                    "fraction_used": bc.get("fraction_used")})
    return out


def text(res: dict) -> str:
    """The status for people: totals, what needs attention, one line per dot, what changed."""
    t = res["totals"]
    lines = [f"=== campaign {res['campaign']} · {res['time']}" + ("" if res["slurm"] else "  (Slurm not reachable)"),
             f"dots {t['dots']}: " + ", ".join(f"{p} {t[p]}" for p in PHASES if t[p])
             + f" · jobs running {t['running_jobs']}, pending {t['pending_jobs']} · campaign {t['campaign_core_h']:,} core-h"]
    b = res.get("budget")
    if b and b.get("used_core_h") is not None:
        lines.append(f"budget {b['account']} since {b['period_start']}: {b['used_core_h']:,} core-h"
                     + (f" of {b['allocation_core_h']:,} ({b['fraction_used']:.0%}), {b['remaining_core_h']:,} left"
                        if b.get("allocation_core_h") else " (no allocation in the profile)")
                     + f" · last 24 h {b['last_24h_core_h']:,}"
                     + (f" · {b['days_left_at_24h_rate']} days at that rate" if b.get("days_left_at_24h_rate") else "")
                     + (f" · fair share {b['fairshare']:.3f}" if b.get("fairshare") is not None else ""))
    lines.append("")
    if res["attention"]:
        lines.append("--- needs attention")
        lines += [f"  {a['level'].upper():5s} {a['dot'] or 'campaign'}: {a['what']}" for a in res["attention"]]
        lines.append("")
    lines.append("--- dots")
    for d in res["dots"]:
        extra = ""
        g = (d.get("chain") or {}).get("geo_opt")
        if d["phase"] == "chain" and d.get("step") == "geo_opt" and g:
            extra = f" (step {g['steps']}, max gradient {g['max_gradient']})"
        lines.append(f"  {d['id']:<44s} {d['phase']:<6s} {d.get('step') or '':<9s}{extra}"
                     + (f"  jobs {','.join(d['jobs'])}" if d.get("jobs") else ""))
    if "changes" in res:
        lines += ["", f"--- since {res['since']}"]
        for c in res["changes"] or [{"kind": "none"}]:
            k = c["kind"]
            lines.append("  " + (
                f"{c['dot']}: {c['from']} -> {c['to']}" if k == "progress" else
                f"{c['level'].upper()} {c['dot'] or 'campaign'}: {c['what']}" if k == "problem" else
                f"resolved {c['dot']}: {c['key']}" if k == "resolved" else
                f"job {c['job']} ({c['name']}) started on {c['partition']}" if k == "job_started" else
                f"job {c['job']} ({c['name']}) ended {c['state']}, {c['core_h']} core-h" if k == "job_ended" else
                f"log: {c['line']}" if k == "log" else
                f"{c['dot']} added ({c['phase']})" if k == "dot_added" else
                f"budget: +{c['since_last']:,} core-h" if k == "budget" and c.get("since_last") is not None else
                "budget: first reading" if k == "budget" else "nothing changed"))
    return "\n".join(lines)


# ---------------- every cluster, from the laptop ----------------
def profiles_with_campaign(names: list[str] | None = None) -> list[dict]:
    """Cluster profiles that have a `campaign` block (where the campaign lives there and how to reach it)."""
    from .cluster import cluster_dirs
    seen, out = set(), []
    for d in cluster_dirs():
        for f in sorted(d.glob("*.yaml")):
            if f.stem in seen:
                continue
            seen.add(f.stem)
            prof = yaml.safe_load(f.read_text()) or {}
            if prof.get("campaign") and (not names or f.stem in names or prof.get("name") in names
                                         or prof["campaign"].get("label") in names):
                out.append({"profile": f.stem, **prof})
    return out


def remote_command(prof: dict, update: bool = True) -> list[str]:
    c = prof["campaign"]
    env = f"PYTHONPATH={c['pythonpath']} " if c.get("pythonpath") else ""
    cmd = (f"{env}{c['python']} -m qd_workflows.cli campaign status {c['path']} --json"
           + (" --update" if update else "") + f" --cluster {prof['profile']}")
    return ["ssh", "-x", "-o", "BatchMode=yes", "-o", "ConnectTimeout=20", c["host"], cmd]


def poll(names: list[str] | None = None, update: bool = True, timeout: float = 300) -> dict:
    """`qdw campaign status --json` on every cluster (over ssh), side by side, with what needs attention."""
    res: dict = {"time": datetime.now().isoformat(timespec="seconds"), "clusters": {}, "attention": [], "changes": []}
    for prof in profiles_with_campaign(names):
        label = prof["campaign"].get("label") or prof["profile"]
        try:
            r = subprocess.run(remote_command(prof, update), capture_output=True, text=True, timeout=timeout)
            if r.returncode != 0:
                raise RuntimeError((r.stderr.strip().splitlines() or [f"exit {r.returncode}"])[-1])
            st = json.loads(r.stdout)
        except (OSError, subprocess.TimeoutExpired, RuntimeError, ValueError) as e:
            res["clusters"][label] = {"error": f"{type(e).__name__}: {e}"}
            res["attention"].append({"cluster": label, "level": "error", "key": "unreachable",
                                     "what": f"status on {prof['campaign']['host']} failed: {e}"})
            continue
        res["clusters"][label] = st
        res["attention"] += [{"cluster": label, **a} for a in st.get("attention", [])]
        res["changes"] += [{"cluster": label, **c} for c in st.get("changes", [])]
    return res


def poll_text(res: dict) -> str:
    parts = []
    for label, st in res["clusters"].items():
        parts.append(f"##### {label}\n" + (f"ERROR {st['error']}" if "error" in st else text(st)))
    return f"=== {res['time']}\n\n" + "\n\n".join(parts)
