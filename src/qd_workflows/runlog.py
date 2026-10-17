"""
Progress log of a run (`$QDW_LOG`, e.g. <run>/pilot.log): one line per event, written by the job
scripts (GPU stage, CP2K chain, CPU stage) as they go, so a run can be followed step by step:

    2026-10-18 09:12:03 | chain | 9582101@c0123-c0124 | geo_opt done: 87 steps, E -10234.567891 Ha, fmax 4.1e-04, 3h12m

Each line: time, stage, Slurm job id and nodes, then what happened and the step's key result.
Warnings and errors start with WARN / ERROR and name the file to read. Without $QDW_LOG nothing
is written. `qdw status <run>` summarises the queue, the log and the progress of the running steps.
"""
from __future__ import annotations

import json
import os
import re
import subprocess
import time
from datetime import datetime
from pathlib import Path


def log_path() -> Path | None:
    p = os.environ.get("QDW_LOG")
    return Path(p) if p else None


def _job() -> str:
    jid, nodes = os.environ.get("SLURM_JOB_ID"), os.environ.get("SLURM_JOB_NODELIST")
    return f"{jid}@{nodes}" if jid else "login"


def write(stage: str, msg: str, path: Path | None = None) -> str:
    """Append one line to the log (and return it)."""
    line = f"{datetime.now():%Y-%m-%d %H:%M:%S} | {stage:<5s} | {_job():<22s} | {msg}"
    p = path or log_path()
    if p:
        p.parent.mkdir(parents=True, exist_ok=True)
        with open(p, "a") as f:
            f.write(line + "\n")
    return line


def hms(seconds: float) -> str:
    s = int(round(seconds))
    if s < 60:
        return f"{s}s"
    if s < 3600:
        return f"{s // 60}m{s % 60:02d}s"
    return f"{s // 3600}h{(s % 3600) // 60:02d}m"


# ---------------- CP2K ----------------
_E = re.compile(r"ENERGY\|\s*Total FORCE_EVAL \( \w+ \) energy\s*\[(?:a\.u\.|hartree)\]:?\s*(-?\d+\.\d+)")
_SCF = re.compile(r"SCF run converged in\s+(\d+) steps")
_OPT = re.compile(r"OPTIMIZATION STEP:\s+(\d+)")
_GRAD = re.compile(r"(?:Max\. gradient\s+=|OPT\| Maximum gradient)\s+([\d.Ee+-]+)")   # CP2K < 2025 / 2026
_TIME = re.compile(r"^\s*CP2K\s+1\s+1\.0\s+\S+\s+\S+\s+\S+\s+([\d.]+)\s*$", re.M)
_HOMO_LUMO = re.compile(r"HOMO - LUMO gap \[eV\] :\s+([\d.]+)")


def cp2k_result(run_dir: str) -> str:
    """Key numbers of a CP2K run from its cp2k.out: energy, geometry steps, SCF steps, gap, wall time."""
    out = Path(run_dir) / "cp2k.out"
    if not out.is_file():
        return f"no {out}"
    text = out.read_text(errors="replace")
    parts = []
    opt = _OPT.findall(text)
    if opt:
        state = "converged" if "GEOMETRY OPTIMIZATION COMPLETED" in text else "not converged"
        parts.append(f"{opt[-1]} steps {state}")
    e = _E.findall(text)
    if e:
        parts.append(f"E {float(e[-1]):.6f} Ha")
    g = _GRAD.findall(text)
    if g:
        parts.append(f"max gradient {float(g[-1]):.1e}")
    scf = _SCF.findall(text)
    if scf and not opt:
        parts.append(f"SCF {scf[-1]} steps")
    gap = _HOMO_LUMO.findall(text)
    if gap:
        parts.append(f"HOMO-LUMO {float(gap[-1]):.3f} eV")
    t = _TIME.findall(text)
    if t:
        parts.append(hms(float(t[-1])))
    elif "PROGRAM ENDED" not in text:
        parts.append("not finished")
    return ", ".join(parts) or "no result in cp2k.out"


def window_result(path: str) -> str:
    """The MO window of `qdw cp2k window`: KS gap and how many MOs QDEX gets."""
    w = json.loads(Path(path).read_text())
    return (f"KS gap {w['gap_ev']:.3f} eV, {w['occupied_in_window']} occupied + {w['virtual_in_window']} virtual MOs "
            f"({w['e_min_ev']:.1f} to {w['e_max_ev']:.1f} eV)")


# ---------------- QDEX ----------------
def qdex_result(qdex_dir: str) -> str:
    """Gaps and lowest excitons from webapp.json (or qdex.out when it is missing)."""
    w = Path(qdex_dir) / "webapp.json"
    if not w.is_file():
        return f"no {w.name}"
    s = json.loads(w.read_text())["summary"]
    parts = []
    for key, lab in (("spin_free", "SF"), ("soc", "SOC")):
        x = s.get(key) or {}
        if x.get("lowest_ev") is not None:
            parts.append(f"{lab}: QP gap {x.get('qp_gap_ev', float('nan')):.2f} eV, lowest exciton "
                         f"{x['lowest_ev']:.3f} eV (f {x.get('lowest_f', 0):.2f}), {x.get('n_states_stored', '?')} states")
    return "; ".join(parts) or "no excitons"


# ---------------- Orchestr.AI PROPS steps ----------------
def _num(v, n=3):
    return f"{v:.{n}g}" if isinstance(v, float) else str(v)


# step -> (summary key, label, format) shown in the log
PROPS_KEYS = {
    "relax": [("steps", "steps", None), ("fmax_eV_A", "fmax", "{:.3f} eV/Å"), ("relaxation_energy_meV_atom", "ΔE", "{:.0f} meV/atom"),
              ("rmsd_A", "RMSD", "{:.2f} Å")],
    "structure": [("bonds_broken", "bonds broken", None), ("bonds_formed", "formed", None), ("ligands_detached", "ligands detached", None),
                  ("core_strain_pct", "core strain", "{:.1f} %")],
    "hessian": [("n_modes", "modes", None), ("n_imaginary", "imaginary", None), ("lowest_cm1", "lowest", "{:.1f} cm-1"),
                ("zpe_eV", "ZPE", "{:.3f} eV")],
    "stability": [("decomposition_dG_300K_per_MA_eV", "ΔG_dec(300 K) per MX", "{:+.3f} eV"),
                  ("binding_dG_300K_per_unit_eV", "binding ΔG per unit", "{:+.3f} eV")],
    "detachment": [("unit", "unit", None), ("n_steps", "removals", None), ("n_relaxations", "relaxations", None)],
    "field": [("alpha_iso_A3", "α_iso", "{:.0f} Å³"), ("dipole_D", "dipole", "{:.1f} D"), ("sum_rule_born_e", "Born sum rule", "{:.1e} e")],
    "vibspec": [("point_group", "point group", None), ("source", "from", None)],
    "solvation": [("dG_solv_dot_eps2.4_eV", "ΔG_solv(ε 2.4)", "{:+.2f} eV"), ("n_structures", "structures", None)],
    "sites": [("n_classes", "site classes", None), ("n_new_relaxations", "new relaxations", None)],
    "report": [("html", "page", None)],
}


def props_result(record_dir: str, step: str) -> str:
    f = Path(record_dir) / "props" / f"{step}.json"
    if not f.is_file():
        return ""
    try:
        s = json.loads(f.read_text()).get("summary") or {}
    except (ValueError, AttributeError):
        return ""
    out = []
    for key, lab, fmt in PROPS_KEYS.get(step, []):
        if key in s and s[key] is not None and not isinstance(s[key], (dict, list)):
            v = s[key]
            out.append(f"{lab} {fmt.format(v) if fmt else _num(v)}")
    if not out:          # unknown step: the first few scalars
        out = [f"{k} {_num(v)}" for k, v in s.items() if isinstance(v, (int, float, str)) and len(str(v)) < 40][:4]
    return ", ".join(out)


_STEP = re.compile(r"\[qdprops\]\s+(\S+): (done in ([\d.]+) s|cached)\s*$")
_REC = re.compile(r"\[qdprops\] (\S+): (\d+) atoms; steps (.*)")
_FAIL = re.compile(r"\[qdprops\]\s+(.*(?:FAILED|failed).*)")
_BATCH = re.compile(r"\[qdprops\] batch: (.*)")


def follow(slurm_out: str, stage: str, records: list[str], stop: str, poll: float = 10.0) -> None:
    """
    Follow an Orchestr.AI PROPS job's output and log each step as it ends (with its key result),
    each record as it starts, and tracebacks; until the file `stop` exists (then one last read).
    """
    pos, buf = 0, ""
    by_name = {Path(r).name: r for r in records}
    current = records[0] if len(records) == 1 else None
    todo: list[str] = []
    path = Path(slurm_out)
    last = False
    while True:
        last = last or Path(stop).exists()
        if path.is_file():
            with open(path, errors="replace") as f:
                f.seek(pos)
                chunk = f.read()
                pos = f.tell()
            buf += chunk
            *lines, buf = buf.split("\n")
            for ln in lines:
                m = _REC.search(ln)
                if m:
                    current = by_name.get(m.group(1), current)
                    todo = [x.strip() for x in m.group(3).split(",")]
                    write(stage, f"{m.group(1)}: {m.group(2)} atoms, steps {m.group(3)}")
                    continue
                m = _STEP.search(ln)
                if m:
                    step = m.group(1)
                    nxt = todo[todo.index(step) + 1] if step in todo and todo.index(step) + 1 < len(todo) else None
                    if m.group(2).startswith("done"):
                        res = props_result(current, step) if current else ""
                        write(stage, f"{step} done in {hms(float(m.group(3)))}" + (f": {res}" if res else "")
                              + (f" -> {nxt}" if nxt else ""))
                    continue
                m = _FAIL.search(ln)
                if m:
                    write(stage, f"{'ERROR' if 'FAILED' in ln else 'WARN'} {m.group(1).strip()} (read {path})")
                    continue
                if ln.startswith("[props]") and "DFT chains" not in ln:      # piece hand-overs, submit failures
                    write(stage, ln[len("[props]"):].strip())
                    continue
                m = _BATCH.search(ln)
                if m:
                    write(stage, f"batch: {m.group(1)}")
                elif "Traceback (most recent call last)" in ln:
                    write(stage, f"ERROR traceback in {path}")
        if last:
            return
        time.sleep(poll)


# ---------------- status ----------------
def _squeue() -> str:
    try:
        r = subprocess.run(["squeue", "--me", "-o", "%.10i %.24j %.18P %.8T %.10M %.10l %.5D %R"],
                           capture_output=True, text=True, timeout=30)
        return r.stdout.rstrip() or r.stderr.strip()
    except (OSError, subprocess.TimeoutExpired) as e:
        return f"(squeue: {e})"


def _count_done(d: Path) -> tuple[int, int]:
    runs = [x for x in d.iterdir() if (x / "cp2k.inp").is_file()] if d.is_dir() else []
    done = sum("PROGRAM ENDED" in (x / "cp2k.out").read_text(errors="replace")[-3000:]
               for x in runs if (x / "cp2k.out").is_file())
    return done, len(runs)


def chain_progress(cp2k_dir: str) -> list[str]:
    """One line per CP2K chain step: done, running (with its progress) or still to do."""
    c = Path(cp2k_dir)
    if not c.is_dir():
        return ["not started"]
    out = []
    g = c / "geo_opt" / "cp2k.out"
    if g.is_file():
        t = g.read_text(errors="replace")
        steps = _OPT.findall(t)
        grad = _GRAD.findall(t)
        state = "done" if "GEOMETRY OPTIMIZATION COMPLETED" in t else "running"
        out.append(f"geo_opt   {state}: {steps[-1] if steps else 0} steps"
                   + (f", max gradient {float(grad[-1]):.1e}" if grad else "")
                   + f", last write {datetime.fromtimestamp(g.stat().st_mtime):%H:%M}")
    for name in ("field", "esp"):
        d, n = _count_done(c / name)
        if n:
            out.append(f"{name:<9s} {d}/{n} single points finished")
    for name, done in (("pdos", any(c.glob("pdos/*.pdos"))), ("window", (c / "window.json").is_file()),
                       ("trexio", any(c.glob("trexio/orbitals*"))), ("trim", (c / "orbitals.h5").is_file()),
                       ("solvation", (c / "solvation.json").is_file())):
        if done:
            out.append(f"{name:<9s} done")
    for q in ("qdex", "qdex_trap"):
        if (c / q / "webapp.json").is_file():
            out.append(f"{q:<9s} done: {qdex_result(str(c / q))}")
        elif (c / q / "config.yaml").is_file():
            out.append(f"{q:<9s} running")
    if (c / ".done").is_file():
        out.append("chain     finished")
    return out


def status(run_dir: str, n_log: int = 25) -> str:
    run = Path(run_dir)
    log = Path(os.environ.get("QDW_LOG") or run / "pilot.log")
    lines = [f"=== {run} · {datetime.now():%Y-%m-%d %H:%M:%S}", "", "--- Slurm jobs", _squeue(), ""]
    for rec in sorted(p.parent for p in run.glob("rec/*/record.json")):
        props = sorted((rec / "props").glob("*.json"), key=lambda p: p.stat().st_mtime) if (rec / "props").is_dir() else []
        lines.append(f"--- {rec.name}")
        lines.append("  props:  " + (", ".join(p.stem for p in props if p.stem not in ("properties",)) or "none yet"))
        lines += ["  cp2k:   " + x for x in chain_progress(str(rec / "cp2k"))]
        lines.append("")
    lines.append(f"--- {log} (last {n_log})")
    if log.is_file():
        lines += log.read_text(errors="replace").rstrip("\n").split("\n")[-n_log:]
    else:
        lines.append("(empty)")
    return "\n".join(lines)
