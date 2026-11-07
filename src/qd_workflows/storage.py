"""
Storage of a cluster account: space and file count of each quota area against its limits, from the
cluster's own quota tool, named in the profile's `storage` block:

    storage:
      - {tool: nhr, command: show-quota --csv --home --work}         # NHR: CSV of every area
      - {tool: beegfs, command: getquota-scratch, name: scratch}       # Hyperion: one area per command
      - {tool: beegfs, command: getquota-home, name: home}

An area at WARN_FRACTION of its limit (the soft limit when there is one) is a warning, at
ERROR_FRACTION an error: the watchdog brief reports both. Areas the campaign writes to fill fastest
in files (a CP2K run directory holds ~10 files), so the file count is checked as well as the space.
"""
from __future__ import annotations

import re
import subprocess

WARN_FRACTION = 0.80
ERROR_FRACTION = 0.95
UNITS = {"B": 1, "KiB": 2 ** 10, "MiB": 2 ** 20, "GiB": 2 ** 30, "TiB": 2 ** 40, "PiB": 2 ** 50}


def _run(cmd: str) -> str | None:
    try:
        r = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=120)
    except (OSError, subprocess.TimeoutExpired):
        return None
    return r.stdout if r.returncode == 0 or r.stdout.strip() else None


def parse_nhr(text: str) -> list[dict]:
    """show-quota --csv: filesystem,label,allocation,blockUsage,blockSoft,blockHard,...,inodeUsage,inodeSoft,inodeHard."""
    rows = [ln.split(",") for ln in text.strip().splitlines()]
    if not rows or rows[0][0] != "filesystem":
        return []
    head = rows[0]
    out = []
    for r in rows[1:]:
        if len(r) != len(head):
            continue
        d = dict(zip(head, r))
        num = lambda k: int(d[k]) if d.get(k, "").isdigit() else 0  # noqa: E731
        out.append({"name": f"{d['filesystem'].lower()}:{d['allocation']}", "path": d["label"],
                    "bytes": num("blockUsage"), "bytes_soft": num("blockSoft"), "bytes_hard": num("blockHard"),
                    "files": num("inodeUsage"), "files_soft": num("inodeSoft"), "files_hard": num("inodeHard")})
    return out


def _size(s: str) -> int:
    m = re.match(r"([\d.]+)\s*([KMGTP]?i?B)", s.strip())
    return int(float(m.group(1)) * UNITS.get(m.group(2), 1)) if m else 0


def parse_beegfs(text: str, name: str) -> list[dict]:
    """getquota-*: '<user>|<id>||  1.38 TiB|  2.00 TiB||  739077| 10000000' (hard limits only)."""
    for ln in text.splitlines():
        parts = [p.strip() for p in ln.split("|")]
        if len(parts) >= 8 and parts[1].isdigit() and parts[3].endswith("B"):
            return [{"name": name, "path": name, "bytes": _size(parts[3]), "bytes_soft": 0, "bytes_hard": _size(parts[4]),
                     "files": int(parts[6]) if parts[6].isdigit() else 0, "files_soft": 0,
                     "files_hard": int(parts[7]) if parts[7].isdigit() else 0}]
    return []


def _frac(used: int, soft: int, hard: int) -> float | None:
    lim = soft or hard
    return round(used / lim, 3) if lim else None


def storage(cluster: dict, run=_run) -> dict | None:
    """Every quota area of the profile's `storage` block, with the fraction used and the problems."""
    spec = cluster.get("storage")
    if not spec:
        return None
    areas, problems = [], []
    for s in spec:
        text = run(s["command"])
        if text is None:
            problems.append({"level": "warn", "key": f"storage_tool:{s['command']}",
                             "what": f"storage: `{s['command']}` could not be run"})
            continue
        areas += parse_nhr(text) if s["tool"] == "nhr" else parse_beegfs(text, s.get("name", s["command"]))
    for a in areas:
        a["bytes_frac"] = _frac(a["bytes"], a["bytes_soft"], a["bytes_hard"])
        a["files_frac"] = _frac(a["files"], a["files_soft"], a["files_hard"])
        for kind in ("bytes", "files"):
            f = a[f"{kind}_frac"]
            if f is None or f < WARN_FRACTION:
                continue
            level = "error" if f >= ERROR_FRACTION else "warn"
            used = f"{a['bytes'] / 2 ** 40:.2f} TiB" if kind == "bytes" else f"{a['files']:,} files"
            problems.append({"level": level, "key": f"storage_{kind}_{level}:{a['path']}",
                             "what": f"storage {a['path']}: {used}, {f:.0%} of the "
                                     f"{'soft ' if a[f'{kind}_soft'] else ''}limit"})
    return {"areas": areas, "problems": problems}
