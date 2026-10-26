"""
Contact geometry and transport from atomistic structures, for the nucleation flux and the growth of dots.

  contact radius   R_se: radius of the sphere with the solvent-excluded volume of a cluster (union of van der Waals
                   spheres closed by a rolling probe of radius `probe`, on a voxel grid); the monomer-cluster contact
                   distance of the Smoluchowski kernel is R_se(1) + R_se(m)
  capacitance      C of any body (in A) by walk-on-spheres: walkers launched uniformly on a sphere of radius R that
                   encloses the body hit it with probability C / R. A walker outside the launch sphere returns to it
                   with probability R / |x|, at a point drawn exactly from the exterior harmonic measure,
                   cos(t) = (r^2 + 1 - 1/s^2) / (2 r), s = 2 xi / (r^2 - 1) + 1 / (r + 1)   (r = |x| / R, xi uniform)
                   Bodies: a union of atomic spheres (exact distance) or a convex polyhedron given by its planes
                   (max_i (n_i.x - h_i) is a distance that never overshoots the body)
  hydrodynamics    R_h = C of the atomic body (Hubbard-Douglas, the approximation of ZENO: within ~1-2 % for compact
                   bodies), D = kT / (c pi eta R_h), c = 6 (stick) or 4 (slip); eta: the solvent viscosity at T (input)

    python -m qd_workflows.transport geometry <thermo.json> [--probe 2.0]   (adds R_se_A, R_h_A to every state)
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

KB = 1.380649e-23          # J/K
PROBE = 2.0                # A, solvent probe radius (an alkane CH2/CH3 group)


def radii(symbols):
    from .reaction.solvation import VDW
    return np.array([VDW.get(s, 2.0) for s in symbols], float)


# ---------------- solvent-excluded radius ----------------
def ses_volume(symbols, positions, probe: float = PROBE, spacing: float = 0.25) -> float:
    """Solvent-excluded volume (A^3): the vdW union closed by a probe sphere (dilate by `probe`, erode by `probe`)."""
    from scipy import ndimage
    p = np.asarray(positions, float)
    r = radii(symbols)
    lo = p.min(0) - r.max() - 2 * probe - 1.0
    hi = p.max(0) + r.max() + 2 * probe + 1.0
    n = np.ceil((hi - lo) / spacing).astype(int) + 1
    g = [lo[i] + spacing * np.arange(n[i]) for i in range(3)]
    X = np.stack(np.meshgrid(*g, indexing="ij"), -1).reshape(-1, 3)
    from scipy.spatial import cKDTree
    t = cKDTree(p)
    # distance to the vdW surface of the union: min_i |x - c_i| - r_i (the 8 nearest centres decide it)
    d, j = t.query(X, k=min(8, len(p)))
    d = np.atleast_2d(d.T).T if d.ndim == 1 else d
    j = np.atleast_2d(j.T).T if j.ndim == 1 else j
    surf = (d - r[j]).min(1).reshape(n)
    accessible = surf <= probe                                   # within probe of the vdW union (SAS body)
    # erode by the probe: keep the points farther than `probe` from outside the accessible body
    dist_out = ndimage.distance_transform_edt(accessible, sampling=spacing)
    excluded = dist_out > probe
    return float(excluded.sum() * spacing ** 3)


def ses_radius(symbols, positions, probe: float = PROBE, spacing: float = 0.25) -> float:
    return (3.0 * ses_volume(symbols, positions, probe, spacing) / (4.0 * math.pi)) ** (1.0 / 3.0)


# ---------------- capacitance by walk-on-spheres ----------------
def _wos(dist, R: float, n: int, eps: float, rng, max_steps: int = 100000) -> tuple[float, float]:
    """Hit probability of a body (dist: (k,3) -> (k,) safe distance, > 0 outside) for walkers launched uniformly on
    the sphere of radius R about the origin; returns (C, standard error)."""
    def on_sphere(k):
        v = rng.normal(size=(k, 3))
        return v / np.linalg.norm(v, axis=1)[:, None]

    x = R * on_sphere(n)
    alive = np.ones(n, bool)
    hit = np.zeros(n, bool)
    for _ in range(max_steps):
        idx = np.flatnonzero(alive)
        if not len(idx):
            break
        xi = x[idx]
        d = dist(xi)
        h = d < eps
        hit[idx[h]] = True
        alive[idx[h]] = False
        idx, xi, d = idx[~h], xi[~h], d[~h]
        xi = xi + d[:, None] * on_sphere(len(idx))                 # jump to the sphere of safe radius d
        r = np.linalg.norm(xi, axis=1)
        out = r > R
        if out.any():
            ro = r[out] / R
            back = rng.random(out.sum()) < 1.0 / ro                 # returns with probability R / |x|
            esc = np.flatnonzero(out)[~back]
            alive[idx[esc]] = False
            ret = np.flatnonzero(out)[back]
            if len(ret):
                rr = ro[back]
                xi_ = rng.random(len(ret))
                s = 2.0 * xi_ / (rr ** 2 - 1.0) + 1.0 / (rr + 1.0)
                cu = np.clip((rr ** 2 + 1.0 - 1.0 / s ** 2) / (2.0 * rr), -1.0, 1.0)
                ax = xi[ret] / np.linalg.norm(xi[ret], axis=1)[:, None]
                # a unit vector at angle acos(cu) from ax, uniform azimuth
                tmp = on_sphere(len(ret))
                perp = tmp - (tmp * ax).sum(1)[:, None] * ax
                perp /= np.linalg.norm(perp, axis=1)[:, None]
                su = np.sqrt(1.0 - cu ** 2)
                xi[ret] = R * (cu[:, None] * ax + su[:, None] * perp)
        x[idx] = xi
    p = hit.mean()
    return R * p, R * math.sqrt(max(p * (1 - p), 1e-12) / n)


def capacitance_spheres(symbols, positions, n: int = 40000, seed: int = 0, scale: float = 1.0) -> tuple[float, float]:
    """Capacitance (A) of a union of atomic spheres (vdW radii x scale)."""
    from scipy.spatial import cKDTree
    p = np.asarray(positions, float)
    p = p - p.mean(0)
    r = radii(symbols) * scale
    t = cKDTree(p)
    k = min(8, len(p))
    R = float(np.max(np.linalg.norm(p, axis=1) + r)) * 1.05

    def dist(x):
        d, j = t.query(x, k=k)
        if k == 1:
            d, j = d[:, None], j[:, None]
        return (d - r[j]).min(1)
    return _wos(dist, R, n, 1e-3 * R, np.random.default_rng(seed))


def capacitance_polyhedron(normals, distances, n: int = 40000, seed: int = 0, dilate: float = 0.0) -> tuple[float, float]:
    """Capacitance (A) of the convex polyhedron P = {x: n_i.x <= h_i} (origin inside), or with dilate = R > 0 of its
    parallel body P + B(R) (every point within R of P: the contact surface of a monomer of radius R on the dot).
    The parallel body needs the exact distance to P: the plane distance where x projects into the face of the
    farthest plane, else the distance to the nearest edge."""
    N = np.asarray(normals, float)
    N = N / np.linalg.norm(N, axis=1)[:, None]
    h = np.asarray(distances, float)
    from scipy.spatial import HalfspaceIntersection
    verts = HalfspaceIntersection(np.hstack([N, -h[:, None]]), np.zeros(3)).intersections
    R = (float(np.linalg.norm(verts, axis=1).max()) + dilate) * 1.05
    if dilate <= 0:
        def dist(x):
            return (x @ N.T - h[None, :]).max(1)
        return _wos(dist, R, n, 1e-3 * R, np.random.default_rng(seed))
    # vertices (merged) and edges: two vertices on two common planes
    tol = 1e-6 * R
    vs = []
    for v in verts:
        if not any(np.linalg.norm(v - w) < 10 * tol for w in vs):
            vs.append(v)
    V = np.array(vs)
    act = np.abs(V @ N.T - h[None, :]) < 1e-5 * R
    ea, eb = [], []
    for a in range(len(V)):
        for b in range(a + 1, len(V)):
            if (act[a] & act[b]).sum() >= 2:
                ea.append(V[a])
                eb.append(V[b])
    ea, eb = np.array(ea), np.array(eb)
    ed = eb - ea
    el2 = np.maximum((ed * ed).sum(1), 1e-30)

    def dist(x):
        s = x @ N.T - h[None, :]
        i = s.argmax(1)
        si = s[np.arange(len(x)), i]
        p = x - si[:, None] * N[i]                                  # projection on the farthest plane
        inside = ((p @ N.T - h[None, :]) <= 1e-9 * R).all(1)
        d = si.copy()
        out = ~inside & (si > 0)
        if out.any():
            xo = x[out]
            t = np.clip(((xo[:, None, :] - ea[None]) * ed[None]).sum(2) / el2[None], 0.0, 1.0)
            q = ea[None] + t[..., None] * ed[None]
            d[out] = np.sqrt(((xo[:, None, :] - q) ** 2).sum(2)).min(1)
        return d - dilate
    return _wos(dist, R, n, 1e-3 * R, np.random.default_rng(seed))


def stokes_d(r_h_A: float, T: float, viscosity_mPa_s: float, c: float = 6.0) -> float:
    """Stokes-Einstein diffusion coefficient (A^2/s) of a body of hydrodynamic radius r_h (A)."""
    return KB * T / (c * math.pi * viscosity_mPa_s * 1e-3 * r_h_A * 1e-10) * 1e20


# ---------------- thermo annotation ----------------
def add_geometry(thermo: dict, base: Path, probe: float = PROBE, log=print) -> int:
    """R_se_A (solvent-excluded radius) and R_h_A (capacitance of the vdW body) of every state."""
    from .xyz import read_xyz
    done = 0
    for nd in thermo["nodes"]:
        for s in nd["states"]:
            if "R_se_A" in s and "R_h_A" in s:
                continue
            sym, pts = read_xyz(base / s["file"])
            s["R_se_A"] = round(ses_radius(list(sym), pts, probe), 3)
            s["R_h_A"] = round(capacitance_spheres(list(sym), pts, n=20000)[0], 3)
            done += 1
    thermo["transport"] = {"R_se": f"solvent-excluded radius, probe {probe} A", "R_h": "capacitance of the vdW body"}
    log(f"[transport] geometry of {done} states")
    return done


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python -m qd_workflows.transport")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("geometry", help="add R_se_A and R_h_A to every state of a thermo file")
    p.add_argument("thermo")
    p.add_argument("--probe", type=float, default=PROBE)
    a = ap.parse_args(argv)
    path = Path(a.thermo)
    th = json.loads(path.read_text())
    add_geometry(th, path.resolve().parent, a.probe)
    path.write_text(json.dumps(th))
    return 0


if __name__ == "__main__":
    sys.exit(main())
