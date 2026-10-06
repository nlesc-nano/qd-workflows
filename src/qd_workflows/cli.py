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
    from .props import bench, prepare
    if args.action == "bench":
        print(json.dumps(bench(args.out, load_cluster(args.cluster), sizes=args.sizes, time=args.time,
                               qos=args.qos, extra=args.extra), indent=1))
        return 0
    info = prepare(args.records, args.out, load_cluster(args.cluster), steps=args.steps.split(",") if args.steps else None,
                   time=args.time, qos=args.qos, max_atoms=args.max_atoms, name=args.name, pieces=args.pieces)
    print(json.dumps(info, indent=1))
    return 0


def _compare(args) -> int:
    from .compare import compare
    res = compare(args.ref, args.new, args.rtol)
    print(f"{res['compared']} values compared, {res['failed']} outside tolerance")
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
    p.add_argument("xyz")
    p.add_argument("--step", required=True, choices=("geo_opt", "pdos", "trexio", "sample"))
    p.add_argument("--project", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--cluster", default="default")
    p.add_argument("--wfn", help="geo_opt restart file read by the later steps")
    p.add_argument("--window", help="window.json from `qdw cp2k window` (trexio step)")
    p.add_argument("--charge", type=int, default=0)
    p.add_argument("--optimizer", default="auto", help="geo_opt: auto (LBFGS above 1,000 atoms), BFGS or LBFGS")
    p.add_argument("--data-dir", help="local copy of the CP2K data files (default: $QDW_CP2K_DATA, then the cluster path)")
    p.add_argument("--time", help="wall time, overriding the cluster profile (e.g. 00:10:00)")
    p.add_argument("--qos", help="QoS, overriding the cluster profile (e.g. test)")

    p = cp.add_parser("chain", help="one job for geo_opt, PDOS, window, TREXIO and trim")
    p.add_argument("xyz", help="MACE-relaxed structure")
    p.add_argument("--project", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--cluster", default="default")
    p.add_argument("--charge", type=int, default=0)
    p.add_argument("--optimizer", default="auto")
    p.add_argument("--data-dir")
    p.add_argument("--record", help="record.json of the structure; adds QDEX as the last step")
    p.add_argument("--time", help="wall time of the whole chain")
    p.add_argument("--qos")

    p = cp.add_parser("window", help="MO window from the pdos step")
    p.add_argument("pdos", nargs="+")
    p.add_argument("-o", "--output")
    p.add_argument("--below", type=float, default=5.0, help="eV below the valence band edge")
    p.add_argument("--above", type=float, default=5.0, help="eV above the conduction band edge")

    p = cp.add_parser("trim", help="keep only the window's MOs in a TREXIO .h5")
    p.add_argument("src")
    p.add_argument("dst")
    p.add_argument("--window", required=True)

    p = cp.add_parser("relaxed", help="last geo_opt frame as an XYZ file")
    p.add_argument("pos")
    p.add_argument("-o", "--output", required=True)

    p = cp.add_parser("label", help="energy and forces of a sample single point")
    p.add_argument("xyz")
    p.add_argument("--out", required=True, help="cp2k.out")
    p.add_argument("--forces", required=True, help="forces.xyz")

    qp = sub.add_parser("qdex", help="QDEX on the TREXIO orbitals").add_subparsers(dest="action", required=True)
    p = qp.add_parser("prepare", help="config.yaml from the record and the MO window")
    p.add_argument("record")
    p.add_argument("--window", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--cluster", default="default")
    p.add_argument("--orbitals", default="../orbitals.h5", help="TREXIO file, relative to --out")
    p.add_argument("--geom", default="../relaxed.xyz", help="geometry, relative to --out")
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--max-states", type=int, help="occupied and virtual MOs in the sBSE, at most (default 1000)")
    p = qp.add_parser("summary", help="frontier levels, gaps and excitons (spin-free and SOC) as JSON")
    p.add_argument("dir", help="QDEX run directory (qdex_electronic.h5, qdex_excitations.h5)")
    p.add_argument("-o", "--output")
    p.add_argument("--bright-f", type=float, default=0.05)
    p = qp.add_parser("webapp", help="webapp.json: summary, band-edge excitons and absorption spectra")
    p.add_argument("dir", help="QDEX run directory")
    p.add_argument("-o", "--output")
    p.add_argument("--window", type=float, default=0.5, help="band-edge states up to this far above the lowest (eV)")

    pp = sub.add_parser("props", help="Orchestr.AI PROPS runs on a GPU cluster").add_subparsers(dest="action", required=True)
    p = pp.add_parser("prepare", help="config.yaml and GPU job script for a set of records")
    p.add_argument("records", nargs="+", help="record directories (record.json + start.xyz)")
    p.add_argument("--out", required=True)
    p.add_argument("--cluster", default="default")
    p.add_argument("--steps", help="comma-separated steps (default: all)")
    p.add_argument("--max-atoms", type=int)
    p.add_argument("--time")
    p.add_argument("--qos")
    p.add_argument("--name", default="props")
    p.add_argument("--pieces", action="store_true", help="resubmit the job until the run finishes (short QoS)")

    p = pp.add_parser("bench", help="GPU job timing force calls and Hessians against the atom count")
    p.add_argument("--out", required=True)
    p.add_argument("--cluster", default="default")
    p.add_argument("--sizes", default="100,300,1000,2000,5000")
    p.add_argument("--time", default="01:00:00")
    p.add_argument("--qos")
    p.add_argument("--extra", default="", help="more orchestr_ai.qd.bench options, e.g. '--analytic-max 1000'")

    p = sub.add_parser("compare", help="compare two props/properties.json (Phase 0 gate)")
    p.add_argument("ref")
    p.add_argument("new")
    p.add_argument("--rtol", type=float, default=1e-3)
    p.add_argument("--show", type=int, default=30)

    args = ap.parse_args(argv)
    return {"cp2k": _cp2k, "qdex": _qdex, "props": _props, "compare": _compare}[args.group](args)


if __name__ == "__main__":
    sys.exit(main())
