"""
Facet geometry of a growing dot from QD_builder, built on the fly (the facet-area size law, growth.py).

A dot of the growth model is its lattice-site plane distances h_f and its CdSe units m. The builder makes the
charge-balanced dot whose lattice-site planes (manifest facet_geometry.lattice_fitted: each plane through the
outermost atom on a lattice site of the Wulff cut, whatever its species; adsorbed ligands excluded) sit at h_f:

    construction distances d_f (QD_builder facet_distances), start d_f = h_f + 1 A, then
    d_f <- d_f + (h_f - d_lat,f)   until every lattice plane is within half a layer of h_f (a few builds)

and returns the built dot's facet areas A_f and edge lengths L_e (per edge type "100|111A", ...) of its lattice-site
polyhedron, scaled to the model's m as A (m/m_b)^(2/3) and L (m/m_b)^(1/3) (m_b the built dot's anion count). Builds
are cached on h_f to `step_A` (0.5 A): the geometry follows the dot as it grows, one build every half angstrom.

Needs a QD_builder whose manifest has facet_geometry.lattice_fitted (QD_builder feat/lattice-site-planes).
"""
from __future__ import annotations

import os
import time
from collections import Counter

LATTICE_NEEDED = ("QD_builder without facet_geometry.lattice_fitted in its manifest: the facet-area law needs the "
                  "lattice-site report (QD_builder branch feat/lattice-site-planes)")


class BuilderGeometry:
    def __init__(self, spec: dict, families, step_A: float = 0.5, tol_A: float = 0.45, iters: int = 5):
        """spec: cif ($VARS expanded), charges, passivation, symmetry, centre (construction origin species), anion
        (the species counted as m), positive_q_mode."""
        self.spec = dict(spec)
        self.spec["cif"] = os.path.expandvars(str(spec["cif"]))
        self.fams = list(families)
        self.step, self.tol, self.iters = float(step_A), float(tol_A), int(iters)
        self.cache: dict = {}
        self.builds = 0
        self.seconds = 0.0

    def recipe(self, dc: dict) -> dict:
        s = self.spec
        return {"charges": s["charges"], "passivation": s["passivation"], "symmetry": s.get("symmetry", {}),
                "construction_origin": {"center_on_species": s["centre"]},
                "facet_distances": {f: round(float(dc[f]), 4) for f in self.fams}}

    def build(self, dc: dict) -> dict:
        """One builder run at construction distances dc: the lattice-site polyhedron and the composition."""
        from .sitekin.build import run_builder
        t0 = time.time()
        s, _x, man = run_builder(self.spec["cif"], self.recipe(dc), self.spec.get("positive_q_mode", "add"))
        self.builds += 1
        self.seconds += time.time() - t0
        geo = man.get("facet_geometry") or {}
        lf = geo.get("lattice_fitted")
        if lf is None:
            raise RuntimeError(LATTICE_NEEDED)
        if "error" in lf:
            raise RuntimeError(f"lattice_fitted: {lf['error']}")
        d = {}
        for p in lf["planes"]:
            d[p["label"]] = max(d.get(p["label"], -1e9), float(p["d_A"]))
        comp = dict(Counter(s.tolist()))
        return {"dc": dict(dc), "d_lat": d, "areas": {f: float(lf["areas_A2"].get(f, 0.0)) for f in self.fams},
                "edges": {k: float(v["length_A"]) for k, v in lf["edges"].items()},
                "volume": float(lf["volume_A3"]), "m": int(comp.get(self.spec["anion"], 0)), "composition": comp}

    def key(self, h: dict) -> tuple:
        return tuple(int(round(float(h[f]) / self.step)) for f in self.fams)

    def dot(self, h: dict) -> dict:
        """The built dot whose lattice-site planes are closest to h (cached on h to step_A)."""
        k = self.key(h)
        if k in self.cache:
            return self.cache[k]
        target = {f: k[i] * self.step for i, f in enumerate(self.fams)}
        dc = {f: max(target[f] + 1.0, 2.0) for f in self.fams}
        best, seen, last = None, set(), ""
        for _ in range(self.iters):
            sig = tuple(round(dc[f], 3) for f in self.fams)
            if sig in seen:
                break
            seen.add(sig)
            try:
                b = self.build(dc)
            except Exception as exc:                     # a cut that holds nothing, a failed charge balance, ...
                if "lattice_fitted" in str(exc) and "QD_builder without" in str(exc):
                    raise
                if best is not None:
                    break
                dc = {f: dc[f] + 1.0 for f in self.fams}  # nothing yet: a larger cut
                last = f"{type(exc).__name__}: {exc}"
                continue
            # a family without a face whose target plane lies outside the dot is hidden in the model too: its plane
            # cannot be met (the other planes set its support) and is not chased
            dev = {f: 0.0 if (b["areas"][f] <= 0.0 and target[f] >= b["d_lat"].get(f, target[f]) - self.tol)
                   else target[f] - b["d_lat"].get(f, target[f]) for f in self.fams}
            err = max(abs(x) for x in dev.values())
            if best is None or err < best[0]:
                best = (err, b)
            if err <= self.tol:
                break
            dc = {f: max(dc[f] + dev[f], 2.0) for f in self.fams}
        if best is None:
            raise RuntimeError(f"no dot built for h = {target} ({last})")
        best[1]["h_err_A"] = best[0]
        self.cache[k] = best[1]
        return best[1]

    def get(self, h: dict, m: float) -> tuple:
        """Facet areas (A^2) and edge lengths (A) for a dot of m units with lattice-site planes h: the built dot's,
        scaled from its own anion count to m. Also returns the built dot's record."""
        b = self.dot(h)
        r = (float(m) / b["m"]) if b["m"] > 0 else 1.0
        area = {f: b["areas"][f] * r ** (2.0 / 3.0) for f in self.fams}
        edges = {e: L * r ** (1.0 / 3.0) for e, L in b["edges"].items()}
        return area, edges, b
