"""Step 5, neb: the monomer rate constants per role from the same site states as K_X (bare CdSe; the growth model
puts in the actual mu_CdSe).

  k_disp  CdSe + dot-X -> P -> [TS] -> dot-CdSe + X      reference G(dot-X)  (states: ref -> M)
  k_dir   CdSe + dot*  -> P -> [TS] -> dot-CdSe          reference G(dot*)   (states: empty -> M)

Endpoints (endpoints()), atoms context + unit + incoming for k_disp, context + incoming for k_dir:
  k_disp IS    relaxed ref, CdSe gap A off the dot above the unit, its Se on the side of the Cl that becomes Se
         FS_A  the site Cd stays (Cd at the M Cd, incoming Se at the M Se); incoming Cd + both Cl leave as CdCl2
         FS_B  the whole site unit leaves; incoming Cd at the M Cd, incoming Se at the M Se
         the leaving CdCl2 goes where the CdSe came in (same centre and axis)
  k_dir  IS    relaxed dot*, CdSe at least h_min above its place in dot-CdSe and gap off the dot
         FS_A  relaxed M
mobile: context atoms within r_mob of the site, or that move > 0.05 A between the two relaxed states; the rest frozen.
Run (run_neb()):
  1. relax the endpoints (float32, L-BFGS to 0.02 eV/A, frozen atoms fixed); free CdSe (float64, 3N - 5 modes, 1 M)
  2. each path: `images` images, IDPP, improved tangent, spring k, FIRE to fmax_loose; every image in one batched call
  3. the path with the lower highest image, if that image is above both ends (else: no barrier, collision-limited):
     a window from the loose images on either side of it (IDPP, 3 images per side and the highest one), loose again,
     then the climbing image, FIRE with 0.05 A steps to fmax_ci; stopped, and the window's highest loose image kept,
     if any atom moves more than max_drift from the window's loose path
  4. float64 energies on the path; partial Hessian at the TS and at the start P on the site's local set + unit +
     incoming; the lowest 8 modes kept with their overlap with the path tangent at the TS. TS check (TS_ok): the
     lowest mode below ts_mode_cm and along the path (overlap >= ts_overlap), the next above ts_second_cm (noise);
     a soft TS mode can sit above imag_cm, so the count n_imag alone is not the check
     dG_act = G(TS) - G(P), k1 = (kT/h) e^(-dG_act/kT); dG_P and dG_TS against G(site) + G(CdSe, 1 M)
A deep minimum before the highest image (no_barrier set by the ends only) needs a look at the path."""
from __future__ import annotations

import json
import math
import os
import time
from pathlib import Path

import numpy as np

from . import engine as en
from . import recipe as rc

KINDS = ("disp", "dir")


def _atoms(symbols, x, tags=None):
    from ase import Atoms
    return Atoms(symbols, positions=x, tags=tags)


def endpoints(kind: str, site: dict, relaxed: dict, p: dict, cation="Cd", anion="Se", cl="Cl") -> dict:
    """Start and end structures of one NEB (see the module text). site: the sites-step entry; relaxed: that site's
    relaxed states (ref, empty, M). Returns {structures: {IS, FS_A[, FS_B]}, meta}."""
    n0 = int(site["n_context"])
    xm = np.array(relaxed["M"]["positions"], float)
    sy_ctx = site["states"]["empty"]["symbols"]
    r_mob, gap = float(p["r_mob"]), float(p["gap"])
    if kind == "disp":
        xr = np.array(relaxed["ref"]["positions"], float)
        ctx = xr[:n0]
        cd_s, cl_a, cl_b = xr[n0], xr[n0 + 1], xr[n0 + 2]          # ref unit order: Cd, Cl (-> Se), Cl (removed)
        near = np.linalg.norm(ctx - cd_s, axis=1) < r_mob
        moved = np.linalg.norm(xm[:n0] - ctx, axis=1) > 0.05
        unit_c = xr[n0:].mean(0)
        nrm = unit_c - ctx[near].mean(0)
        nrm /= np.linalg.norm(nrm)
        ax = (cl_a - cl_b) - ((cl_a - cl_b) @ nrm) * nrm
        ax /= np.linalg.norm(ax)
        half = float(p["d_cdse"]) / 2
        for h in np.arange(2.0, 12.0, 0.1):
            c_in = unit_c + h * nrm
            if min(np.linalg.norm(xr - q, axis=1).min() for q in (c_in + half * ax, c_in - half * ax)) >= gap:
                break
        else:
            raise RuntimeError("no room above the site")
        se_in, cd_in = c_in + half * ax, c_in - half * ax
        dcl = float(p["d_cdcl"])
        leave = [c_in, c_in + dcl * ax, c_in - dcl * ax]
        fs_ctx = ctx.copy()
        mobile = np.flatnonzero(near | moved)
        fs_ctx[mobile] = xm[:n0][mobile]
        sy = list(sy_ctx) + [cation, cl, cl, cation, anion]
        st = {"IS": np.vstack([xr, cd_in, se_in]),
              "FS_A": np.vstack([fs_ctx, xm[n0], leave[1], leave[2], leave[0], xm[n0 + 1]]),
              "FS_B": np.vstack([fs_ctx, leave[0], leave[1], leave[2], xm[n0], xm[n0 + 1]])}
        unit, incoming, routes, centre = [n0, n0 + 1, n0 + 2], [n0 + 3, n0 + 4], ["A", "B"], unit_c
    elif kind == "dir":
        xe = np.array(relaxed["empty"]["positions"], float)
        cd_m, se_m = xm[n0], xm[n0 + 1]
        near = np.linalg.norm(xe - cd_m, axis=1) < r_mob
        moved = np.linalg.norm(xm[:n0] - xe, axis=1) > 0.05
        centre = 0.5 * (cd_m + se_m)
        nrm = centre - xe[near].mean(0)
        nrm /= np.linalg.norm(nrm)
        ax = (se_m - cd_m) - ((se_m - cd_m) @ nrm) * nrm
        ax /= np.linalg.norm(ax)
        half = float(p["d_cdse"]) / 2
        for h in np.arange(float(p["h_min"]), 12.0, 0.1):
            c_in = centre + h * nrm
            pts = [c_in + q * half * ax for q in (1, -1)]
            if min(np.linalg.norm(xe - q, axis=1).min() for q in pts) >= gap:
                break
        else:
            raise RuntimeError("no room above the site")
        se_in, cd_in = pts
        fs_ctx = xe.copy()
        mobile = np.flatnonzero(near | moved)
        fs_ctx[mobile] = xm[:n0][mobile]
        sy = list(sy_ctx) + [cation, anion]
        st = {"IS": np.vstack([xe, cd_in, se_in]), "FS_A": np.vstack([fs_ctx, cd_m, se_m])}
        unit, incoming, routes = [], [n0, n0 + 1], ["A"]
    else:
        raise ValueError(kind)
    meta = {"kind": kind, "n_context": n0, "mobile": [int(i) for i in mobile], "local_phva": list(site["local"]),
            "unit": unit, "incoming": incoming, "routes": routes, "normal": nrm.tolist(), "symbols": sy,
            "start_height_A": float(np.linalg.norm(c_in - centre))}
    return {"structures": {k: v.round(5).tolist() for k, v in st.items()}, "meta": meta}


def marked(meta: dict) -> list[str]:
    """Symbols with the moving atoms marked for viewing: site Cd -> Hg, site Cl -> Br, incoming Cd -> Zn, Se -> Te."""
    s = list(meta["symbols"])
    for i, e in zip(meta["unit"], ("Hg", "Br", "Br")):
        s[i] = e
    for i, e in zip(meta["incoming"], ("Zn", "Te")):
        s[i] = e
    return s


def write_xyz(path: Path, symbols, x, comment: str):
    path.write_text(f"{len(symbols)}\n{comment}\n" + "".join(f"{e:2s} {q[0]:12.5f} {q[1]:12.5f} {q[2]:12.5f}\n"
                                                           for e, q in zip(symbols, np.asarray(x, float))))


def run_neb(engine, ep: dict, g_site: float, r: dict, workdir: Path, log=print) -> dict:
    """One NEB from endpoints() output; resumable through workdir/result.json."""
    from ase.calculators.calculator import Calculator, all_changes
    from ase.calculators.singlepoint import SinglePointCalculator
    from ase.constraints import FixAtoms
    from ase.io import read, write
    from ase.mep import NEB
    from ase.optimize import FIRE, LBFGS

    p, T = r["neb"], float(r["T_K"])
    ph = r["phva"]
    meta = ep["meta"]
    workdir.mkdir(parents=True, exist_ok=True)
    out = workdir / "result.json"
    R = json.loads(out.read_text()) if out.exists() else {"kind": meta["kind"], "G_site": g_site}
    t0 = time.time()
    n0, sy = meta["n_context"], meta["symbols"]
    fixed = sorted(set(range(n0)) - set(meta["mobile"]))
    phva_set = meta["local_phva"] + meta["unit"] + meta["incoming"]
    routes = meta["routes"]
    inc = meta["incoming"]

    def lg(m):
        log(f"[{time.time() - t0:6.0f} s] {m}")

    def save():
        tmp = out.with_suffix(".tmp")
        tmp.write_text(json.dumps(R))
        os.replace(tmp, out)

    class B32(Calculator):
        implemented_properties = ["energy", "forces"]

        def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
            super().calculate(atoms, properties, system_changes)
            e, f = engine.ev32([(atoms.get_chemical_symbols(), atoms.positions.copy())])[0]
            self.results = {"energy": float(e), "forces": np.asarray(f, float)}

    class BatchNEB(NEB):
        """every interior image in one call, as single-point results; off while ASE's IDPP runs on this object"""
        batch = True

        def get_forces(self):
            if not self.batch:
                return super().get_forces()
            inner = self.images[1:-1]
            res = engine.ev32([(im.get_chemical_symbols(), im.positions.copy()) for im in inner])
            for im, (e, f) in zip(inner, res):
                im.calc = SinglePointCalculator(im, energy=float(e), forces=np.asarray(f, float))
            return super().get_forces()

    def frozen(at):
        at.set_constraint(FixAtoms(indices=fixed))
        return at

    def spc(im):
        e, f = engine.ev32([(im.get_chemical_symbols(), im.positions.copy())])[0]
        im.calc = SinglePointCalculator(im, energy=float(e), forces=np.asarray(f, float))

    # 1. endpoints and free CdSe
    ends = {}
    for name in ["IS"] + [f"FS_{q}" for q in routes]:
        f = workdir / f"{name}_relaxed.xyz"
        if name in R.get("ends", {}) and f.exists():
            ends[name] = frozen(read(f))
            continue
        at = frozen(_atoms(sy, np.array(ep["structures"][name], float)))
        at.calc = B32()
        opt = LBFGS(at, logfile=None)
        ok = opt.run(fmax=0.02, steps=1500)
        write(f, at, format="extxyz")
        e64 = float(engine.ev64([(sy, at.positions.copy())])[0][0])
        R.setdefault("ends", {})[name] = {"E": e64, "steps": opt.nsteps, "converged": bool(ok),
                                          "d_incoming": float(np.linalg.norm(at.positions[inc[1]] - at.positions[inc[0]]))}
        ends[name] = at
        save()
        lg(f"{name}: relaxed in {opt.nsteps} steps (converged {ok}), E64 {e64:.4f} eV")
    if "CdSe" not in R:
        R["CdSe"] = en.free_molecule(engine, [sy[inc[0]], sy[inc[1]]], [[0, 0, 0], [float(p["d_cdse"]), 0, 0]], T,
                                     float(ph["step_A"]), float(ph["floor_cm"]))
        save()
        lg(f"free CdSe: G({T:.0f} K, 1 M) {R['CdSe']['G']:.4f} eV")

    n_img, k_spring = int(p["images"]), float(p["spring"])

    def band(fs):
        images = [ends["IS"].copy()] + [ends["IS"].copy() for _ in range(n_img)] + [ends[fs].copy()]
        for im in images:
            frozen(im)
        spc(images[0])
        spc(images[-1])
        neb = BatchNEB(images, k=k_spring, climb=False, method="improvedtangent", allow_shared_calculator=True)
        neb.batch = False
        neb.interpolate(method="idpp", apply_constraint=False)
        neb.batch = True
        return neb

    # 2. loose bands
    for q in routes:
        key = f"loose_{q}"
        if key in R:
            continue
        neb = band(f"FS_{q}")
        opt = FIRE(neb, logfile=None)
        ok = opt.run(fmax=float(p["fmax_loose"]), steps=400)
        e = [im.get_potential_energy() for im in neb.images]
        R[key] = {"E32": e, "steps": opt.nsteps, "converged": bool(ok), "barrier32": max(e[1:-1]) - e[0]}
        write(workdir / f"path_loose_{q}.xyz", neb.images, format="extxyz")
        save()
        lg(f"path {q}: loose band in {opt.nsteps} steps, highest image {R[key]['barrier32']:+.3f} eV above IS, profile "
           + " ".join(f"{v - e[0]:+.2f}" for v in e))

    # 3. window climbing image on the lower path
    route = min(routes, key=lambda q: R[f"loose_{q}"]["barrier32"])
    e_lo = R[f"loose_{route}"]["E32"]
    i_hi = 1 + int(np.argmax(e_lo[1:-1]))
    no_barrier = bool(e_lo[i_hi] <= max(e_lo[0], e_lo[-1]) + 0.02)
    if not no_barrier and "ci" not in R:
        lo = read(workdir / f"path_loose_{route}.xyz", ":")
        a, b, c = (frozen(lo[i].copy()) for i in (i_hi - 1, i_hi, i_hi + 1))
        images = [a]
        for u, v in ((a, b), (b, c)):
            seg = [u.copy()] + [u.copy() for _ in range(3)] + [v.copy()]
            NEB(seg).interpolate(method="idpp", apply_constraint=False)
            images += [frozen(im) for im in seg[1:]]
        spc(images[0])
        spc(images[-1])
        neb = BatchNEB(images, k=k_spring, climb=False, method="improvedtangent", allow_shared_calculator=True)
        ok_lo = FIRE(neb, logfile=None, maxstep=0.05).run(fmax=float(p["fmax_loose"]), steps=200)
        win0 = [im.positions.copy() for im in neb.images]
        write(workdir / f"path_win_{route}.xyz", neb.images, format="extxyz")
        neb.climb = True
        opt = FIRE(neb, logfile=None, maxstep=0.05)
        ok, drift, steps = False, 0.0, 0
        while steps < 600:
            ok = opt.run(fmax=float(p["fmax_ci"]), steps=25)
            steps += 25
            drift = max(float(np.linalg.norm(im.positions - x0, axis=1).max()) for im, x0 in zip(neb.images, win0))
            if ok or drift > float(p["max_drift"]):
                break
        used = drift <= float(p["max_drift"])
        if not used:
            for im, x0 in zip(neb.images[1:-1], win0[1:-1]):
                im.calc = None
                im.positions = x0
        write(workdir / f"path_ci_{route}.xyz", neb.images, format="extxyz")
        R["ci"] = {"route": route, "window": [i_hi - 1, i_hi + 1], "loose_window_converged": bool(ok_lo),
                   "steps": steps, "converged": bool(ok), "drift": drift, "ci_used": used}
        save()
        lg(f"path {route}: window around loose image {i_hi}, climbing image {steps} steps (converged {ok}), "
           f"largest move {drift:.2f} A" + ("" if used else " -> stopped, the window's loose path kept"))

    # 4. float64 path, PHVA at the TS and at P
    if no_barrier:
        names = (("IS", ends["IS"]),)
    else:
        path = read(workdir / f"path_ci_{route}.xyz", ":")
        if "E64" not in R["ci"]:
            R["ci"]["E64"] = [float(e) for e, _f in engine.ev64([(im.get_chemical_symbols(), im.positions.copy())
                                                                  for im in path])]
            save()
        e64 = R["ci"]["E64"]
        i_ts = 1 + int(np.argmax(e64[1:-1]))
        names = (("TS", path[i_ts]), ("IS", ends["IS"]))
        write(workdir / "TS.xyz", path[i_ts], format="extxyz")
    for name, at in names:
        if name in R.get("phva", {}):
            continue
        fr, vec = en.phva(engine, at.get_chemical_symbols(), at.positions.copy(), phva_set, float(ph["step_A"]), vectors=True)
        R.setdefault("phva", {})[name] = {
            "F_vib": en.vib_free(fr, T, 0, float(ph["floor_cm"])), "n_imag": int(sum(f < float(p["imag_cm"]) for f in fr)),
            "lowest": sorted(fr)[:8], "E64": float(engine.ev64([(at.get_chemical_symbols(), at.positions.copy())])[0][0])}
        if name == "TS":
            # how much each low mode lies along the path: mass-weighted tangent between the TS's neighbours
            tan = (path[i_ts + 1].positions - path[i_ts - 1].positions)[phva_set]
            tan = (tan * np.sqrt(en.masses(at.get_chemical_symbols())[phva_set])[:, None]).reshape(-1)
            tan /= np.linalg.norm(tan)
            R["phva"]["TS"]["overlap"] = [float(abs(vec[:, j] @ tan)) for j in np.argsort(fr)[:8]]
        save()
        lg(f"PHVA {name}: F_vib {R['phva'][name]['F_vib']:.4f} eV, lowest modes "
           + ", ".join(f"{f:.0f}" for f in R["phva"][name]["lowest"][:3]) + " cm^-1"
           + (f", overlap with the path {R['phva'][name]['overlap'][0]:.2f}" if name == "TS" else ""))
    kt = en.KB * T
    g_p = R["phva"]["IS"]["E64"] + R["phva"]["IS"]["F_vib"]
    e_is = R["ends"]["IS"]["E"]
    res = {"kind": meta["kind"], "route": route, "no_barrier": no_barrier, "dG_P_1M": g_p - g_site - R["CdSe"]["G"],
           "dE_end_minus_start": R["ends"][f"FS_{route}"]["E"] - e_is, "G_P": g_p}
    if no_barrier:
        res["dE_loose_max_from_IS"] = e_lo[i_hi] - e_lo[0]
    else:
        g_ts = R["phva"]["TS"]["E64"] + R["phva"]["TS"]["F_vib"]
        dg = g_ts - g_p
        res.update({"dE_act": e64[i_ts] - e_is, "G_TS": g_ts, "dG_act": dg,
                    "k1_s": 1.380649e-23 * T / 6.62607015e-34 * math.exp(-dg / kt),
                    "dG_TS_1M": g_ts - g_site - R["CdSe"]["G"], "TS_lowest_cm": R["phva"]["TS"]["lowest"],
                    "TS_overlap": R["phva"]["TS"]["overlap"], "TS_n_imag": R["phva"]["TS"]["n_imag"],
                    # the check: one negative mode, and it is the one along the path (a soft one may sit above imag_cm)
                    "TS_ok": bool(R["phva"]["TS"]["lowest"][0] < float(p.get("ts_mode_cm", -5.0))
                                  and R["phva"]["TS"]["lowest"][1] > float(p.get("ts_second_cm", -10.0))
                                  and R["phva"]["TS"]["overlap"][0] >= float(p.get("ts_overlap", 0.5))),
                    "converged_CI": R["ci"]["ci_used"] and R["ci"]["converged"]})
    R["result"] = res
    save()
    if no_barrier:
        lg(f"{meta['kind']}: no barrier (highest loose image {res['dE_loose_max_from_IS']:+.3f} eV); P {res['dG_P_1M']:+.3f} eV vs 1 M")
    else:
        lg(f"{meta['kind']} path {route}: dE_act {res['dE_act']:+.3f}, dG_act {res['dG_act']:+.3f} eV, k1 {res['k1_s']:.2e} /s; "
           f"P {res['dG_P_1M']:+.3f}, TS {res['dG_TS_1M']:+.3f} eV vs CdSe 1 M; TS lowest {res['TS_lowest_cm'][0]:.0f} cm^-1, "
           f"overlap {res['TS_overlap'][0]:.2f} -> TS check {'ok' if res['TS_ok'] else 'FAILED'}")
    return R


def choose_sites(r: dict, summary: list) -> dict:
    """Per role, the site of its largest environment (the protocol's one site per role)."""
    out = {}
    for row in summary:
        envs = [v for v in row["environments"] if not v["no_whole_unit"]]
        if envs:
            out[row["role"]] = max(envs, key=lambda v: v["n"])["site_name"]
    return out


def prepare(r: dict, out: Path, sizes=None, roles=None, log=print) -> list[Path]:
    """Endpoints and marked xyz for every chosen site and kind: neb/<tag>/<site>/<kind>/endpoints.json."""
    out = Path(out)
    paths = []
    p, mat = r["neb"], r["material"]
    for size in (sizes or [rc.tag(s) for s in p["sizes"]]):
        t = size if isinstance(size, str) else rc.tag(size)
        job = json.loads((out / "sites" / f"{t}.json").read_text())
        summary = json.loads((out / "sites" / f"{t}_summary.json").read_text())
        states = json.loads((out / "states" / f"{t}.json").read_text())
        for role, name in choose_sites(r, summary).items():
            if roles and role not in roles:
                continue
            st = states["sites"][name]
            for kind in KINDS:
                d = out / "neb" / t / name / kind
                if (d / "result.json").exists():
                    log(f"{t} {role:26s} {kind:4s} {name}: has a result, kept (remove the directory to redo it)")
                    continue
                d.mkdir(parents=True, exist_ok=True)
                ep = endpoints(kind, job["sites"][name], st["relaxed"], p, mat["cation"], mat["anion"], r["ligand"]["anion"])
                ep["meta"].update({"site": name, "role": role, "size_tag": t,
                                   "G_site": st["G"]["ref" if kind == "disp" else "empty"]})
                (d / "endpoints.json").write_text(json.dumps(ep))
                ms = marked(ep["meta"])
                for k, x in ep["structures"].items():
                    write_xyz(d / f"{k}_marked.xyz", ms, x, f"{kind} {name} {k} (unrelaxed); Hg = site Cd, Br = site Cl, "
                              "Zn = incoming Cd, Te = incoming Se")
                log(f"{t} {role:26s} {kind:4s} {name}: {len(ms)} atoms, {len(ep['meta']['mobile'])} mobile, "
                    f"CdSe start {ep['meta']['start_height_A']:.2f} A above the site")
                paths.append(d)
    return paths


def run(r: dict, out: Path, engine, dirs=None, log=print) -> list[Path]:
    """Every prepared NEB under out/neb (or the given dirs) without a result yet."""
    out = Path(out)
    todo = [Path(d) for d in dirs] if dirs else sorted(q.parent for q in (out / "neb").glob("*/*/*/endpoints.json"))
    done = []
    for d in todo:
        ep = json.loads((d / "endpoints.json").read_text())
        res = d / "result.json"
        if res.exists() and "result" in json.loads(res.read_text()):
            continue
        log(f"== {d.relative_to(out) if d.is_relative_to(out) else d}")
        run_neb(engine, ep, float(ep["meta"]["G_site"]), r, d, log)
        done.append(d)
    return done
