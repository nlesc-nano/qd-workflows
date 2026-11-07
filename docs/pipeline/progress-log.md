# Following a run

Every job of a run writes to one progress log, so the whole run can be followed step by step without
reading the Slurm outputs: the GPU stage, each CP2K chain and the CPU stage add a line when a job
starts and ends, and when a step starts or ends, with the step's key result.

Set `QDW_LOG` before preparing the jobs; the job scripts keep it, and the jobs they submit inherit it:

```bash
export QDW_LOG=$RUN/pilot.log
qdw props prepare $RUN/rec/* --out $RUN/batch --cluster nhr
```

A run of one dot then reads like this:

```text
2026-10-18 09:00:12 | gpu   | 9582100@ggpu03          | start: props, NVIDIA A100-SXM4-80GB, limit 04:00:00
2026-10-18 09:00:40 | gpu   | 9582100@ggpu03          | CdSe-Se-Cd176Se147Cl58-clean: 381 atoms, steps relax, structure, hessian, ...
2026-10-18 09:12:51 | gpu   | 9582100@ggpu03          | relax done in 12m11s: steps 140, fmax 0.009 eV/Å, ΔE -98 meV/atom, RMSD 0.61 Å -> structure
...
2026-10-18 11:40:02 | gpu   | 9582100@ggpu03          | submitted the CP2K chain(s) 9582201 and the CPU stage 9582202 (after the chains)
2026-10-18 11:52:30 | chain | 9582201@c[0101-0102]    | start CdSe-...: 336 ranks on 2 node(s), 1290000 MB per node, limit 12:00:00
2026-10-18 11:52:31 | chain | 9582201@c[0101-0102]    | geo_opt started
2026-10-18 15:03:44 | chain | 9582201@c[0101-0102]    | geo_opt done: 61 steps converged, E -10234.567891 Ha, max gradient 4.1e-04, 3h11m
2026-10-18 15:03:50 | chain | 9582201@c[0101-0102]    | field: 15 single points, 6 side by side, 56 ranks, 193500M each
...
2026-10-18 19:20:11 | chain | 9582201@c[0101-0102]    | qdex done: SF: QP gap 2.61 eV, lowest exciton 2.420 eV (f 0.41), 1650 states; SOC: ...
```

Each line: time, stage (`gpu`, `chain`, `cpu`), Slurm job id and nodes, then the event. Lines that
need attention start with `WARN` or `ERROR` and name the file to read. A job that stops early says
in which step and with which exit code.

| Stage | What is logged |
|---|---|
| GPU and CPU stages | job start (GPU model or cores) and end; every Orchestr.AI step as it ends, with its time, its key numbers and the next step; failures and tracebacks; the hand-over to the chains and the CPU stage |
| CP2K chain | job start (ranks, nodes, memory); geo_opt (steps, energy, max gradient, wall time); the field and ESP single points (how many, side by side, ranks and memory each, then the total time); PDOS and TREXIO; the MO window (KS gap, MOs kept); QDEX and the trap-filtered QDEX (QP gap, lowest exciton SF and SOC); the dashboards |

The GPU and CPU stages log their steps through a follower (`qdw log follow`), which reads the job's
own output while the job runs; nothing changes in Orchestr.AI.

## Status

`qdw status <run>` shows your Slurm jobs, the progress of every record (props outputs so far; geo_opt
steps and max gradient, field and ESP single points finished, which chain steps are done) and the end
of the log. `--watch 30` refreshes it every 30 s:

```bash
qdw status $RUN --watch 30
```

The commands behind it (`qdw log event|cp2k|window|qdex|follow`) are in the
[command line reference](../reference/cli.md); without `QDW_LOG` they print the line instead.
