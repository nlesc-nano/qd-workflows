"""
Compute budget of a cluster account: what its jobs have been charged since the allocation period
started, against the allocation, and how fast it is going.

From the profile's `budget` block (clusters/<name>.yaml):

    budget:
      account: bep00157
      period_start: 2026-10-01      # the allocation period (NHR: quarterly)
      period_end: 2026-12-31
      allocation_core_h: 500000     # from the cluster's portal; null when there is none (fair share only)
      charge:                       # core-h charged, by partition (longest prefix wins)
        cpu-genoa: {per_node_h: 192}
        gpu-a100:shared: {per_gpu_h: 150}

Slurm's accounting (sacct) of every user of the account, so the charge is the project's, not only
ours; a partition without a rate is charged its cores x hours (CPUs / threads_per_core).
"""
from __future__ import annotations

import getpass
import re
import subprocess
from datetime import date, datetime, timedelta

WARN_FRACTION = 0.75
ERROR_FRACTION = 0.90


def _run(cmd: list[str]) -> str | None:
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    except (OSError, subprocess.TimeoutExpired):
        return None
    return r.stdout if r.returncode == 0 else None


def _t(s: str) -> datetime | None:
    try:
        return datetime.fromisoformat(s)
    except (TypeError, ValueError):
        return None


def _date(v) -> date | None:
    if v is None or isinstance(v, date):
        return v
    return date.fromisoformat(str(v))


def rate(charge: dict, partition: str) -> dict | None:
    keys = [k for k in charge if partition == k or partition.startswith(k)]
    return charge[max(keys, key=len)] if keys else None


def job_core_h(job: dict, charge: dict, threads_per_core: int, hours: float) -> float:
    """Core-h charged for `hours` of one job."""
    r = rate(charge, job["partition"])
    if r and "per_gpu_h" in r and job["gpus"]:
        return r["per_gpu_h"] * job["gpus"] * hours
    if r and "per_node_h" in r:
        return r["per_node_h"] * job["nodes"] * hours
    if r and "per_core_h" in r:
        return r["per_core_h"] * job["cpus"] / threads_per_core * hours
    return job["cpus"] / threads_per_core * hours


def sacct_jobs(account: str | None, start: date) -> list[dict] | None:
    """Jobs started since `start`: every user's of `account`, or only ours without one."""
    who = ["-a", "-A", account] if account else []
    out = _run(["sacct", "-X", "-n", "-P", *who, "-S", start.isoformat(),
                "-o", "JobID,User,JobName%80,Partition,State,Start,End,ElapsedRaw,AllocTRES,NNodes"])
    if out is None:
        return None
    jobs = []
    for ln in out.splitlines():
        f = ln.split("|")
        if len(f) < 10 or not f[7].isdigit() or int(f[7]) == 0:
            continue
        tres = dict(kv.split("=", 1) for kv in f[8].split(",") if "=" in kv)
        gpus = int(tres.get("gres/gpu") or next((v for k, v in tres.items() if k.startswith("gres/gpu:")), 0) or 0)
        jobs.append({"id": f[0], "user": f[1], "name": f[2], "partition": f[3], "state": f[4].split()[0],
                     "start": _t(f[5]), "end": _t(f[6]), "elapsed_s": int(f[7]),
                     "cpus": int(re.sub(r"\D", "", tres.get("cpu", "0")) or 0), "gpus": gpus,
                     "nodes": int(f[9]) if f[9].isdigit() else 1})
    return jobs


def fairshare(account: str, user: str) -> float | None:
    out = _run(["sshare", "-n", "-P", "-A", account, "-u", user, "-o", "Account,User,FairShare"])
    for ln in (out or "").splitlines():
        f = ln.split("|")
        if len(f) == 3 and f[1].strip() == user:
            try:
                return float(f[2])
            except ValueError:
                return None
    return None


def budget(cluster: dict, now: datetime | None = None, jobs: list[dict] | None = None,
           belongs=None, user: str | None = None) -> dict | None:
    """Charge of the account since the period started, the campaign's share, the last 24 h and the outlook."""
    b = cluster.get("budget")
    if not b:
        return None
    now = now or datetime.now()
    user = user or getpass.getuser()
    start = _date(b.get("period_start")) or date(now.year, 3 * ((now.month - 1) // 3) + 1, 1)
    end = _date(b.get("period_end"))
    charge = b.get("charge") or {}
    tpc = int(cluster.get("threads_per_core", 1))
    if jobs is None:
        jobs = sacct_jobs(b["account"], start)
    res: dict = {"cluster": cluster.get("name"), "account": b["account"], "period_start": start.isoformat(),
                 "period_end": end.isoformat() if end else None, "allocation_core_h": b.get("allocation_core_h"),
                 "problems": []}
    if jobs is None:
        res["problems"].append({"level": "warn", "key": "sacct", "what": "sacct could not be run: no budget figures"})
        return res
    day_ago = now - timedelta(hours=24)
    used = mine = campaign = last24 = 0.0
    by_part: dict[str, float] = {}
    by_user: dict[str, float] = {}
    for j in jobs:
        h = j["elapsed_s"] / 3600
        c = job_core_h(j, charge, tpc, h)
        used += c
        by_part[j["partition"]] = by_part.get(j["partition"], 0.0) + c
        by_user[j["user"]] = by_user.get(j["user"], 0.0) + c
        if j["user"] == user:
            mine += c
        if belongs and belongs(j["name"]):
            campaign += c
        s = j["start"]
        e = j["end"] or (s + timedelta(seconds=j["elapsed_s"]) if s else None)
        if s and e and e > day_ago:
            last24 += job_core_h(j, charge, tpc, (e - max(s, day_ago)).total_seconds() / 3600)
    res.update({"used_core_h": round(used), "mine_core_h": round(mine), "campaign_core_h": round(campaign),
                "last_24h_core_h": round(last24),
                "by_partition": {k: round(v) for k, v in sorted(by_part.items(), key=lambda x: -x[1])},
                "by_user": {k: round(v) for k, v in sorted(by_user.items(), key=lambda x: -x[1])},
                "fairshare": fairshare(b["account"], user)})
    alloc = b.get("allocation_core_h")
    if alloc:
        frac = used / alloc
        res["fraction_used"] = round(frac, 3)
        res["remaining_core_h"] = round(alloc - used)
        if last24 > 0:
            res["days_left_at_24h_rate"] = round((alloc - used) / last24, 1)
        if frac >= ERROR_FRACTION:
            res["problems"].append({"level": "error", "key": "budget_90",
                                    "what": f"{b['account']}: {frac:.0%} of {alloc:,} core-h used"})
        elif frac >= WARN_FRACTION:
            res["problems"].append({"level": "warn", "key": "budget_75",
                                    "what": f"{b['account']}: {frac:.0%} of {alloc:,} core-h used"})
        if end and last24 > 0:
            empty = now + timedelta(days=(alloc - used) / last24)
            if empty.date() < end:
                res["problems"].append({"level": "warn", "key": "budget_rate",
                                        "what": f"{b['account']}: at the last 24 h rate ({last24:,.0f} core-h/day) "
                                                f"the allocation runs out on {empty:%Y-%m-%d}, before {end}"})
    return res
