"""
Even Z-type placement for QD_builder's z_type_target (the MX_q added or displaced after the charge balance).

The builder's dot and its charge-balance ligands are untouched; only where z_type_target puts (or takes) the units
changes, so that every facet of a family receives the same number on equivalent sites:

  equivalent sites  core sites (the empty cation site above a free surface anion for addition, the unit's cation for
                    displacement) that are images of each other under the crystal point group
                    (builder.polyhedron._point_group_rotations) about the dot's centre
  instance          the facet (edge, vertex) of the family a site lies on: its nearest construction plane(s) of that
                    family's label, each plane labelled by the site roles of the core atoms on it
  counts            the builder's split of the units over the families (its divisor rule on the per_family weights),
                    each family's count rounded to a whole number per facet
  order             families in the order of per_family (the caller lists them strongest-bound first), each filled
                    round-robin over its facets, class by class; families share the anion sites, so a family never
                    takes the sites of one listed before it
Everything else (host sites, crowding by ligands, the anions of each unit, clashes, the displacement rule) is the
builder's.

    with even_placement():
        builder_main(...)        # or generate_library._build(...)
"""
from __future__ import annotations

import contextlib
from collections import Counter, defaultdict

import numpy as np

CTX: dict = {}


def _fam_labels(fam: str):
    if fam.startswith("edge:"):
        return fam[5:].split("|")
    if fam.startswith("vertex:"):
        return fam[7:].split("+")
    if fam in ("edge", "vertex", "all"):
        return None
    return [fam[6:] if fam.startswith("facet:") else fam]


def instance(fam: str, x) -> tuple:
    """Plane ids of the facet / edge / vertex of family fam nearest to x."""
    x = np.asarray(x, float)
    depth = [(d - n @ x, p, lab) for p, (n, d, lab) in enumerate(CTX["planes"])]
    labs = _fam_labels(fam)
    if labs is None:
        k = 2 if fam == "edge" else 3 if fam == "vertex" else 1
        return tuple(sorted(p for _dp, p, _l in sorted(depth)[:k]))
    out = []
    for lab in labs:
        cand = [(dp, p) for dp, p, l2 in depth if l2 == lab and p not in out]
        if cand:
            out.append(min(cand)[1])
    return tuple(sorted(out))


def n_instances(fam: str) -> int:
    labs = _fam_labels(fam)
    if labs is None or len(labs) != 1:
        return 1
    return max(1, sum(1 for _n, _d, lab in CTX["planes"] if lab == labs[0]))


def klass(x) -> tuple:
    """Canonical position of the equivalence class of a site at x (images about the dot's centre)."""
    y = np.asarray(x, float) - CTX["centre"]
    return min(tuple(np.round(R @ y * 2) / 2) for R in CTX["ops"])


def _even_order(fam, pos):
    cls = defaultdict(list)
    for i, x in enumerate(pos):
        cls[klass(x)].append(i)
    order = []
    for g in sorted(cls.values(), key=lambda g: (-len(g), min(g))):
        by = defaultdict(list)
        for i in sorted(g):
            by[instance(fam, pos[i])].append(i)
        lanes = [by[k] for k in sorted(by)]
        for r in range(max(len(v) for v in lanes)):
            order += [v[r] for v in lanes if r < len(v)]
    return order


def _lanes(fam, idxs, sites):
    pos = np.asarray([sites[i]["pos"] for i in idxs], float)
    rank = {j: r for r, j in enumerate(_even_order(fam, pos))}
    by = defaultdict(list)
    for j, i in enumerate(idxs):
        by[instance(fam, pos[j])].append((rank[j], i))
    return [[i for _r, i in sorted(v)] for _k, v in sorted(by.items())]


def _even_units(fam, plan, sites):
    """Units a plan puts evenly: the family's facets times the fewest units on any of them."""
    if not plan:
        return 0
    c = Counter(instance(fam, sites[i]["pos"]) for i, _t in plan)
    n = n_instances(fam)
    return n * (min(c.values()) if len(c) >= n else 0)


def _in_order(families, counts, exhausted):
    for fam, wt in families:
        if fam not in exhausted and counts[fam] < wt:
            return fam
    return None


def _install():
    """Patch z_type_target (and the planner it uses) for the duration of a build; returns the undo."""
    import builder.z_type_ligand_posttreat as ZL
    import builder.z_type_target as ZT
    from builder.polyhedron import _point_group_rotations

    orig = {"run": ZT.run_z_type_target, "add": ZT._add_units, "remove": ZT._remove_units,
            "next": ZT._next_family, "plans": ZL._family_plans}

    def even_plans(sites, cands, n, family_sites, distribution, seed):
        users = ZL._site_users(cands)
        used, placed = set(), set()
        base = orig["plans"](sites, cands, n, family_sites, distribution, seed)
        plans = {}
        for fam, idxs in family_sites.items():
            if not idxs:
                plans[fam] = []
                continue
            lanes = _lanes(fam, idxs, sites)
            ptr = [0] * len(lanes)
            u_loc, p_loc = set(used), set(placed)
            plan = []
            while True:
                moved = False
                for li, lane in enumerate(lanes):
                    while ptr[li] < len(lane):
                        i = lane[ptr[li]]
                        ptr[li] += 1
                        if i in p_loc:
                            continue
                        free = ZL._free_sites(i, cands, users, u_loc, p_loc)
                        if len(free) < n:
                            continue
                        plan.append((i, free[:n]))
                        p_loc.add(i)
                        u_loc.update(c["id"] for c in free[:n])
                        moved = True
                        break
                if not moved:
                    break
            want = int(CTX["counts"].get(fam, len(idxs)))
            own = [(i, t) for i, t in base[fam] if i not in placed and not any(c["id"] in used for c in t)]
            if own:                              # the builder's own compatible set, round-robin, when it is more even
                take_of = dict(own)
                ln = _lanes(fam, [i for i, _t in own], sites)
                alt = [(lane[r], take_of[lane[r]]) for r in range(max(map(len, ln))) for lane in ln if r < len(lane)]
                if _even_units(fam, alt[:want], sites) > _even_units(fam, plan[:want], sites):
                    plan = alt
            for i, take in plan[:want]:
                placed.add(i)
                used.update(c["id"] for c in take)
            plans[fam] = plan
        return plans

    def counts_for(families, n):
        shares = {f: 0 for f, _w in families}
        for _ in range(n):
            f = orig["next"](families, shares, set())
            if f is None:
                break
            shares[f] += 1
        out = {}
        for f, _w in families:
            step = n_instances(f)
            out[f] = int(round(shares[f] / step)) * step
        return shares, out

    def per_instance(fam, xs):
        return sorted(Counter(instance(fam, x) for x in xs).values())

    def add_units(syms, pts, cfg, struct, planes, roles, families, n_add, complex_, distribution, seed):
        shares, counts = counts_for(families, n_add)
        fam_list = [(f, float(counts[f])) for f, _w in families if counts[f] > 0]
        CTX["counts"] = {f: c for f, c in counts.items()}
        n0 = len(syms)
        if fam_list:
            ZT._next_family, ZL._family_plans = _in_order, even_plans
            try:
                out = orig["add"](syms, pts, cfg, struct, planes, roles, fam_list, sum(counts.values()), complex_,
                                  distribution, seed)
            finally:
                ZT._next_family, ZL._family_plans = orig["next"], orig["plans"]
            s2, p2, per_family, order, entry = out
        else:
            s2, p2, per_family, order, entry = list(syms), np.asarray(pts, float), {}, [], None
        new_metal = [p for s, p in zip(s2[n0:], p2[n0:]) if CTX["metal"] and s == CTX["metal"]]
        full = {}
        for f, w in families:
            rec = dict(per_family.get(f, {"weight": w, "added": 0}))
            rec.update(weight=w, share=shares[f], even_count=counts[f], facets=n_instances(f))
            if rec.get("added"):
                rec["per_facet"] = per_instance(f, [x for x, fam in zip(new_metal, order) if fam == f])
            full[f] = rec
        return s2, p2, full, order, entry

    def remove_units(syms, pts, cfg, roles, core, families, n_remove, metal, ligand, q, cif_path):
        shares, counts = counts_for(families, n_remove)
        fam_list = [(f, float(counts[f])) for f, _w in families if counts[f] > 0]
        cuts = ZT.derive_pair_cuts_from_cif(cif_path, cfg.charges, safety=1.00)
        work_syms, work_pts = list(syms), np.asarray(pts, float).copy()
        work_roles, work_core = list(roles), np.asarray(core, bool).copy()
        native_anions = {s for s, c in zip(syms, core) if c and cfg.charges.get(s, 0) < 0}
        done_inst, removed_at = defaultdict(int), defaultdict(list)
        cnt = {f: 0 for f, _w in fam_list}
        initial, exhausted, order = {}, set(), []
        while sum(cnt.values()) < sum(counts.values()):
            bonds = ZT._bonds(work_syms, work_pts, cfg.charges, cuts)
            nbr = defaultdict(list)
            for i, j in bonds:
                nbr[i].append(j)
                nbr[j].append(i)
            n_comp = ZT._n_components(len(work_syms), bonds, set())
            centre = work_pts[work_core].mean(axis=0)

            def candidates(fam):
                out = []
                for c, s in enumerate(work_syms):
                    if s != metal or not work_core[c] or not ZT._matches(fam, work_roles[c]):
                        continue
                    n_own = sum(1 for x in nbr[c] if work_syms[x] == ligand)
                    if n_own == 0:
                        continue
                    cn_core = sum(1 for x in nbr[c] if work_syms[x] in native_anions)
                    out.append((done_inst[(fam, instance(fam, work_pts[c]))], klass(work_pts[c]), cn_core, -n_own,
                                -float(np.linalg.norm(work_pts[c] - centre)), c))
                return [o[-1] for o in sorted(out)]

            if not initial:
                initial = {f: len(candidates(f)) for f, _w in fam_list}
            fam = _in_order(fam_list, cnt, exhausted)
            if fam is None:
                break
            ok = False
            for c in candidates(fam):
                unit = ZT._removal_unit(c, work_syms, work_pts, nbr, metal, ligand, q, work_core)
                if unit is None:
                    continue
                removed = {c, *unit}
                if ZT._n_components(len(work_syms), bonds, removed) > n_comp:
                    continue
                inst = instance(fam, work_pts[c])
                done_inst[(fam, inst)] += 1
                removed_at[fam].append(inst)
                keep = [i for i in range(len(work_syms)) if i not in removed]
                work_syms = [work_syms[i] for i in keep]
                work_pts = work_pts[keep]
                work_roles = [work_roles[i] for i in keep]
                work_core = work_core[keep]
                ok = True
                break
            if not ok:
                exhausted.add(fam)
                continue
            cnt[fam] += 1
            order.append(fam)
        per_family = {}
        for f, w in families:
            rec = {"weight": w, "candidates": initial.get(f, 0), "removed": cnt.get(f, 0), "share": shares[f],
                   "even_count": counts[f], "facets": n_instances(f)}
            if removed_at[f]:
                rec["per_facet"] = sorted(Counter(removed_at[f]).values())
            per_family[f] = rec
        return work_syms, work_pts, per_family, order

    def run(syms, pts, cfg, struct, planes, cif_path, *, roles, core, ligand):
        P = np.asarray(pts, float)
        core_m = np.asarray(core, bool)
        N = np.array([np.asarray(n, float) for n, _d in planes])
        Dd = np.array([float(d) for _n, d in planes])
        ops = _point_group_rotations(struct)
        # label of a plane: its Miller family (from the normal in the lattice, as the builder labels facets); for a
        # polar family (n and -n in different orbits) A or B from the facet roles of the atoms on its orbit
        digits = []
        for n in N:
            g = struct.lattice.matrix @ n
            g = g / min(abs(v) for v in g if abs(v) > 1e-6)
            digits.append("".join(str(int(round(abs(v)))) for v in sorted(g, key=lambda v: -abs(v))))
        orbit_of, orbits = {}, []
        for p in range(len(planes)):
            if p in orbit_of:
                continue
            orb = [q for q in range(len(planes)) if any(np.dot(R @ N[p], N[q]) > 0.999 for R in ops)]
            for q in orb:
                orbit_of[q] = len(orbits)
            orbits.append(orb)
        lab = {}
        for oi, orb in enumerate(orbits):
            fam = digits[orb[0]]
            polar = not any(np.dot(-N[orb[0]], N[q]) > 0.999 for q in orb)
            if not polar:
                lab.update({q: fam for q in orb})
                continue
            tally = Counter()
            for x, r, c in zip(P, roles or [], core_m):
                if c and (r or "").startswith(f"facet:{fam}"):
                    same = [q for q in range(len(planes)) if digits[q] == fam]
                    q = min(same, key=lambda q: Dd[q] - N[q] @ x)
                    if orbit_of[q] == oi:
                        tally[r[6:]] += 1
            lab.update({q: tally.most_common(1)[0][0] if tally else None for q in orb})
        labelled = [(N[p], Dd[p], lab[p]) for p in range(len(planes))]
        cations = {s for s, c in zip(syms, core_m) if c and cfg.charges.get(s, 0) > 0}
        CTX.clear()
        CTX.update(planes=labelled, ops=ops, centre=P[core_m].mean(axis=0),
                   counts={}, metal=next(iter(cations)) if len(cations) == 1 else None)
        return orig["run"](syms, pts, cfg, struct, planes, cif_path, roles=roles, core=core, ligand=ligand)

    ZT.run_z_type_target, ZT._add_units, ZT._remove_units = run, add_units, remove_units

    def undo():
        ZT.run_z_type_target, ZT._add_units, ZT._remove_units = orig["run"], orig["add"], orig["remove"]
        ZT._next_family, ZL._family_plans = orig["next"], orig["plans"]
        CTX.clear()
    return undo


@contextlib.contextmanager
def even_placement(enabled: bool = True):
    """Z-type target units placed evenly over the facets (see the module docstring) inside the block."""
    if not enabled:
        yield
        return
    undo = _install()
    try:
        yield
    finally:
        undo()
