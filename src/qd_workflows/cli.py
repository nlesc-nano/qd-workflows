"""qdw: command line of qd-workflows."""
from __future__ import annotations

import argparse
import json
import sys

from .cluster import load_cluster


def _cp2k(args) -> int:
    from .cp2k import inputs, outputs, trexio
    from .xyz import write_xyz

    if args.action == "prepare":
        window = json.loads(open(args.window).read()) if args.window else None
        info = inputs.prepare(args.xyz, args.step, args.project, args.out, load_cluster(args.cluster),
                              wfn=args.wfn, window=window, charge=args.charge, optimizer=args.optimizer,
                              data_dir=args.data_dir, time=args.time, qos=args.qos)
        print(json.dumps({k: info[k] for k in ("project", "step", "n_atoms", "basis_functions",
                                                "occupied_mos", "cell_A", "resources")}, indent=1))
    elif args.action == "chain":
        info = inputs.chain(args.xyz, args.project, args.out, load_cluster(args.cluster), cluster_ref=args.cluster,
                            record=args.record, charge=args.charge, optimizer=args.optimizer, data_dir=args.data_dir,
                            time=args.time, qos=args.qos)
        print(json.dumps({k: info[k] for k in ("project", "n_atoms", "basis_functions", "resources", "chain")}, indent=1))
    elif args.action == "window":
        win = outputs.mo_window(args.pdos, args.below, args.above)
        text = json.dumps(win, indent=1)
        if args.output:
            open(args.output, "w").write(text + "\n")
        print(text)
    elif args.action == "trim":
        print(json.dumps(trexio.trim(args.src, args.dst, json.loads(open(args.window).read())), indent=1))
    elif args.action == "relaxed":
        symbols, coords = outputs.last_frame(args.pos)
        write_xyz(args.output, symbols, coords, "CP2K geo_opt, last frame")
        print(f"{len(symbols)} atoms -> {args.output}")
    elif args.action == "esp":
        from .cp2k import esp
        if args.esp_action == "prepare":
            print(" ".join(esp.prepare(args.base, args.out, stripped=args.stripped, monomers=args.monomers)))
        elif args.esp_action == "status":
            print(esp.status(args.dir))
        else:
            data = esp.write(args.dir, args.output)
            print(json.dumps({n: {k: round(v, 4) if isinstance(v, float) else v for k, v in s.items()}
                              for n, s in data["species"].items()}, indent=1))
    elif args.action == "label":
        from .xyz import read_xyz
        symbols, _ = read_xyz(args.xyz)
        energy = outputs.read_energy(args.out)
        forces = outputs.read_forces(args.forces, len(symbols))
        print(json.dumps({"energy_ev": energy, "max_force_ev_a": float(abs(forces).max()), "n_atoms": len(symbols)}))
    return 0


def _qdex(args) -> int:
    from . import qdex
    if args.action == "prepare":
        info = qdex.prepare(args.record, args.window, args.out, load_cluster(args.cluster), orbitals=args.orbitals,
                            geom=args.geom, threads=args.threads, max_states=args.max_states)
    elif args.action == "webapp":
        info = qdex.webapp(args.dir, args.window)
        if args.output:
            open(args.output, "w").write(json.dumps(info, separators=(",", ":")) + "\n")
        print(json.dumps(info["summary"], indent=1))
        return 0
    else:
        info = qdex.summary(args.dir, args.bright_f)
        if args.output:
            open(args.output, "w").write(json.dumps(info, indent=1) + "\n")
    print(json.dumps(info, indent=1))
    return 0


def _props(args) -> int:
    from .props import bench, chains, prepare
    if args.action == "chains":
        import yaml
        recs = list(args.records or [])
        if args.config:
            recs += yaml.safe_load(open(args.config))["props"]["records"]
        ids = chains(recs, load_cluster(args.cluster), args.cluster, max_atoms=args.max_atoms, submit=args.submit,
                     time=args.time, qos=args.qos)
        print(" ".join(ids))                   # last line: the job ids (the GPU job's hand-off reads it)
        return 0
    if args.action == "bench":
        print(json.dumps(bench(args.out, load_cluster(args.cluster), sizes=args.sizes, time=args.time,
                               qos=args.qos, extra=args.extra), indent=1))
        return 0
    info = prepare(args.records, args.out, load_cluster(args.cluster), steps=args.steps.split(",") if args.steps else None,
                   time=args.time, qos=args.qos, max_atoms=args.max_atoms, name=args.name, pieces=args.pieces, md=args.md,
                   cpu_time=args.cpu_time)
    print(json.dumps(info, indent=1))
    return 0


def _reaction(args) -> int:
    from .reaction import spec
    if args.action == "check":
        rx = spec.load(args.reaction)
        print(f"{rx['name']}: balanced, {rx['equation']}")
        for n, s in rx["species"].items():
            print(f"  {n:<14s} nu {rx['nu'][n]:>4}  {s['n_atoms']:>5d} atoms  charge {s['charge']:+d}  mult {s['mult']}  "
                  f"E {s['energy_eV']:.6f} eV  {s['phase']}  c {s['c']:g} M")
        de = sum(v * rx["species"][n]["energy_eV"] for n, v in rx["nu"].items())
        print(f"  dE_DFT = {de:+.4f} eV ({de * 23.060548:+.2f} kcal/mol) per reaction as written")
        return 0
    if args.action == "prepare":
        from .reaction.jobs import prepare
        print(json.dumps(prepare(args.reaction, args.out, load_cluster(args.cluster), gpu_time=args.gpu_time,
                                 cpu_time=args.cpu_time, qos=args.qos, data_dir=args.data_dir,
                                 test_queue=args.test_queue), indent=1))
        return 0
    from .reaction.analyse import analyse
    s = analyse(args.run)
    for name, m in s["models"].items():
        print(f"{name:>16s}: dG {m['dG_eV']:+.3f} eV  dG0 {m['dG0_eV']:+.3f}  dH0 {m['dH0_eV']:+.3f}  "
              f"TdS0 {m['TdS0_eV']:+.3f}  (T = {s['T_K']} K)  dG = 0 at {m['T_dG_zero_K'] or '-'} K")
    for lab, v in s.get("dE_alternatives_eV_per", {}).items():
        print(f"  dE {lab}: {'n/a' if v is None else f'{v:+.3f} eV'}")
    for n, f in s["flags"].items():
        for x in f:
            print(f"  ! {n}: {x}")
    print(f"-> {args.run}/results/")
    return 0


def _compare(args) -> int:
    from .compare import compare
    res = compare(args.ref, args.new, args.rtol)
    print(f"{res['compared']} values compared, {res['failed']} outside tolerance")
    if res["path"]:
        worst = max(res["path"], key=lambda p: p["diff"])
        better = sum(p["diff"] < -0.05 for p in res["path"])
        print(f"  desorption path: cumulative dE vs reference at {len(res['path'])} levels, worst {worst['diff']:+.3f} eV "
              f"({worst['key']}), {better} levels better by > 0.05 eV")
    for f in res["failures"][: args.show]:
        print(f"  {f['key']}: ref {f['ref']:.6g}  new {f['new']:.6g}  diff {f['diff']:+.3g}")
    if res["only_in_ref"] or res["only_in_new"]:
        print(f"  only in ref: {len(res['only_in_ref'])}  only in new: {len(res['only_in_new'])}")
    return 1 if res["failed"] else 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="qdw", description=__doc__)
    sub = ap.add_subparsers(dest="group", required=True)
    cp = sub.add_parser("cp2k", help="CP2K PBE track").add_subparsers(dest="action", required=True)

    p = cp.add_parser("prepare", help="input, geometry and job script for one step")
    p.add_argument("xyz", help='structure (XYZ)')
    p.add_argument("--step", required=True, choices=("geo_opt", "pdos", "trexio", "sample", "charges"), help='which CP2K step')
    p.add_argument("--project", required=True, help='CP2K project name (file prefix)')
    p.add_argument("--out", required=True, help='output directory')
    p.add_argument("--cluster", default="default", help='cluster profile (clusters/<name>.yaml or a path)')
    p.add_argument("--wfn", help="geo_opt restart file read by the later steps")
    p.add_argument("--window", help="window.json from `qdw cp2k window` (trexio step)")
    p.add_argument("--charge", type=int, default=0, help='total charge')
    p.add_argument("--optimizer", default="auto", help="geo_opt: auto (LBFGS above 1,000 atoms), BFGS or LBFGS")
    p.add_argument("--data-dir", help="local copy of the CP2K data files (default: $QDW_CP2K_DATA, then the cluster path)")
    p.add_argument("--time", help="wall time, overriding the cluster profile (e.g. 00:10:00)")
    p.add_argument("--qos", help="QoS, overriding the cluster profile (e.g. test)")

    p = cp.add_parser("chain", help="one job for geo_opt, PDOS, window, TREXIO and trim")
    p.add_argument("xyz", help="MACE-relaxed structure")
    p.add_argument("--project", required=True, help='CP2K project name (file prefix)')
    p.add_argument("--out", required=True, help='output directory of the chain')
    p.add_argument("--cluster", default="default", help='cluster profile')
    p.add_argument("--charge", type=int, default=0, help='total charge')
    p.add_argument("--optimizer", default="auto", help='geo_opt optimiser: auto (LBFGS above 1,000 atoms), BFGS or LBFGS')
    p.add_argument("--data-dir", help='local copy of the CP2K data files (default: $QDW_CP2K_DATA, then the cluster path)')
    p.add_argument("--record", help="record.json of the structure; adds QDEX as the last step")
    p.add_argument("--time", help="wall time of the whole chain")
    p.add_argument("--qos", help='QoS, overriding the cluster profile')

    p = cp.add_parser("window", help="MO window from the pdos step")
    p.add_argument("pdos", nargs="+", help="the pdos step's .pdos files")
    p.add_argument("-o", "--output", help='window.json to write')
    p.add_argument("--below", type=float, default=5.0, help="eV below the valence band edge")
    p.add_argument("--above", type=float, default=5.0, help="eV above the conduction band edge")

    p = cp.add_parser("trim", help="keep only the window's MOs in a TREXIO .h5")
    p.add_argument("src", help='TREXIO .h5 from the trexio step')
    p.add_argument("dst", help='trimmed .h5 to write')
    p.add_argument("--window", required=True, help='window.json from qdw cp2k window')

    p = cp.add_parser("relaxed", help="last geo_opt frame as an XYZ file")
    p.add_argument("pos", help='geo_opt trajectory (<project>-pos-1.xyz)')
    p.add_argument("-o", "--output", required=True, help='XYZ to write')

    ep = cp.add_parser("esp", help="solvation cubes: intact and fully stripped dot, monomers -> COSMO (solvation.json)"
                       ).add_subparsers(dest="esp_action", required=True)
    p = ep.add_parser("prepare", help="the runs, from the chain's single-point input (field/base/cp2k.inp)")
    p.add_argument("base", help="the chain's single-point input (field/base/cp2k.inp)")
    p.add_argument("out", help='directory of the esp runs')
    p.add_argument("--stripped", help="fully stripped dot: props/detach_<K>.xyz")
    p.add_argument("--monomers", help="props/monomers.xyz (frames MA, MX)")
    p = ep.add_parser("status", help='"all done" or the runs still to do')
    p.add_argument("dir", help='directory of the esp runs')
    p = ep.add_parser("collect", help="COSMO energy and isosurface area of every run")
    p.add_argument("dir", help='directory of the esp runs')
    p.add_argument("-o", "--output", required=True, help='solvation.json to write')

    p = cp.add_parser("label", help="energy and forces of a sample single point")
    p.add_argument("xyz", help='structure of the single point')
    p.add_argument("--out", required=True, help="cp2k.out")
    p.add_argument("--forces", required=True, help="forces.xyz")

    qp = sub.add_parser("qdex", help="QDEX on the TREXIO orbitals").add_subparsers(dest="action", required=True)
    p = qp.add_parser("prepare", help="config.yaml from the record and the MO window")
    p.add_argument("record", help='record.json of the dot')
    p.add_argument("--window", required=True, help='window.json from qdw cp2k window')
    p.add_argument("--out", required=True, help='QDEX run directory')
    p.add_argument("--cluster", default="default", help='cluster profile')
    p.add_argument("--orbitals", default="../orbitals.h5", help="TREXIO file, relative to --out")
    p.add_argument("--geom", default="../relaxed.xyz", help="geometry, relative to --out")
    p.add_argument("--threads", type=int, default=8, help='OpenMP threads for QDEX')
    p.add_argument("--max-states", type=int, help="occupied and virtual MOs in the sBSE, at most (default 1000)")
    p = qp.add_parser("summary", help="frontier levels, gaps and excitons (spin-free and SOC) as JSON")
    p.add_argument("dir", help="QDEX run directory (qdex_electronic.h5, qdex_excitations.h5)")
    p.add_argument("-o", "--output", help='JSON file to write (default: print)')
    p.add_argument("--bright-f", type=float, default=0.05, help='oscillator strength of the first bright exciton')
    p = qp.add_parser("webapp", help="webapp.json: summary, band-edge excitons and absorption spectra")
    p.add_argument("dir", help="QDEX run directory")
    p.add_argument("-o", "--output", help='webapp.json to write')
    p.add_argument("--window", type=float, default=0.5, help="band-edge states up to this far above the lowest (eV)")

    pp = sub.add_parser("props", help="Orchestr.AI PROPS runs on a GPU cluster").add_subparsers(dest="action", required=True)
    p = pp.add_parser("prepare", help="config.yaml and GPU job script for a set of records")
    p.add_argument("records", nargs="+", help="record directories (record.json + start.xyz)")
    p.add_argument("--out", required=True, help='job directory')
    p.add_argument("--cluster", default="default", help='cluster profile')
    p.add_argument("--steps", help="comma-separated steps (default: all)")
    p.add_argument("--max-atoms", type=int, help='records above this many atoms are left out (default 2,000)')
    p.add_argument("--time", help='wall time of the GPU job')
    p.add_argument("--qos", help='QoS, overriding the cluster profile')
    p.add_argument("--name", default="props", help='job name')
    p.add_argument("--pieces", action="store_true", help="resubmit the job until the run finishes (short QoS)")
    p.add_argument("--md", action="store_true", help="also run MD (constant-T and ramp replicas): the MD subset")
    p.add_argument("--cpu-time", default=None, help="wall time of the CPU (assembly) stage, if the profile splits the run")

    p = pp.add_parser("chains", help="one CP2K chain (geo_opt, field, PDOS, TREXIO, QDEX) per record at <record>/cp2k")
    p.add_argument("records", nargs="*", help="record directories")
    p.add_argument("--config", help="a PROPS config.yaml whose props.records are used")
    p.add_argument("--cluster", default="default", help='cluster profile')
    p.add_argument("--max-atoms", type=int, default=2000, help='records above this many atoms get no chain')
    p.add_argument("--submit", action="store_true", help="sbatch each chain; prints the job ids")
    p.add_argument("--time", help='wall time of each chain')
    p.add_argument("--qos", help='QoS, overriding the cluster profile')
    p = pp.add_parser("bench", help="GPU job timing force calls and Hessians against the atom count")
    p.add_argument("--out", required=True, help='job directory')
    p.add_argument("--cluster", default="default", help='cluster profile')
    p.add_argument("--sizes", default="100,300,1000,2000,5000", help='atom counts of the benchmark spheres')
    p.add_argument("--time", default="01:00:00", help='wall time')
    p.add_argument("--qos", help='QoS, overriding the cluster profile')
    p.add_argument("--extra", default="", help="more orchestr_ai.qd.bench options, e.g. '--analytic-max 1000'")

    rp = sub.add_parser("reaction", help="reaction free energies vs T, eps and concentrations").add_subparsers(
        dest="action", required=True)
    p = rp.add_parser("check", help="read the reaction file: balance, units, dE_DFT")
    p.add_argument("reaction", help="reaction.yaml")
    p = rp.add_parser("prepare", help="run directory: MACE (GPU), CP2K (CPU) and analysis jobs, submit.sh")
    p.add_argument("reaction", help="reaction.yaml")
    p.add_argument("--out", required=True, help='run directory')
    p.add_argument("--cluster", default="default", help='cluster profile')
    p.add_argument("--gpu-time", default="01:00:00", help='wall time of the MACE (GPU) job')
    p.add_argument("--cpu-time", help="wall time of the CP2K job (default: estimated from the largest species)")
    p.add_argument("--qos", help='QoS, overriding the cluster profile')
    p.add_argument("--data-dir", help="local copy of the CP2K data files (default: $QDW_CP2K_DATA, then the cluster path)")
    p.add_argument("--test-queue", action="store_true",
                   help="all jobs in the test partitions (1 h); CP2K in pieces that resubmit themselves")
    p = rp.add_parser("analyse", help="results/: results.json, summary.json, plots, reaction.html")
    p.add_argument("run", help="run directory")

    p = sub.add_parser("compare", help="compare two props/properties.json (Phase 0 gate)")
    p.add_argument("ref", help='reference properties.json')
    p.add_argument("new", help='properties.json to compare')
    p.add_argument("--rtol", type=float, default=1e-3, help='relative tolerance for values without an absolute one')
    p.add_argument("--show", type=int, default=30, help='show at most this many differences')

    args = ap.parse_args(argv)
    return {"cp2k": _cp2k, "qdex": _qdex, "props": _props, "reaction": _reaction,
            "compare": _compare}[args.group](args)


if __name__ == "__main__":
    sys.exit(main())
