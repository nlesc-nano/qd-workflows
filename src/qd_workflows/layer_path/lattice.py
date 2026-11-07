"""
Zinc-blende lattice of the layer path: sites, plane cuts, symmetry classes and dot keys.

Grid: a/4 units. Anion sites (Se, or Cl the charge balance put on an anion site) on fcc points (even coordinates,
sum = 0 mod 4), Cd sites = anion + (1, 1, 1). Plane levels are d = n.r with the unnormalised normal n; one layer is
SPACING[family] in d: {100} 2 (a Cd and an anion plane), {111} 4 (one double layer), {110} 2.

Identities (both up to the 24 Td operations and a lattice translation):
- class key (canon): the lattice-site set of a layer state, anions species-blind ("A", "Cd");
- dot key (dot_key): every atom of a built dot on the a/16 grid, Z-type Cl included, species kept.
The registry stores both as repr() strings, as the scripts that wrote it did; eval(key, {"np": np}) reads them back.
"""
from __future__ import annotations

import itertools
from pathlib import Path

import numpy as np

A = 6.14054240                                   # CdSe zinc-blende lattice constant (A), as the scan's CIF
Q = A / 4.0
TD_OPS = [np.diag(sg) @ np.eye(3, dtype=int)[list(perm)]
          for perm in itertools.permutations(range(3))
          for sg in itertools.product((1, -1), repeat=3) if np.prod(sg) > 0]
_SIGNS = list(itertools.product((1, -1), repeat=3))
NORMALS = {
    "100": [tuple(int(v) for v in s * row) for row in np.eye(3, dtype=int) for s in (1, -1)],
    "111A": [s for s in _SIGNS if np.prod(s) < 0],
    "111B": [s for s in _SIGNS if np.prod(s) > 0],
    "110": sorted({tuple(sg * v for sg, v in zip(s3, perm)) for perm in itertools.permutations((1, 1, 0))
                   for s3 in itertools.product((1, -1), repeat=3)}),
}
SPACING = {"100": 2, "111A": 4, "111B": 4, "110": 2}
FAMS = ("100", "111A", "111B")
PLANES = [(f, n) for f in FAMS for n in NORMALS[f]]                  # the 14 {100}/{111} planes
PLANES26 = PLANES + [("110", n) for n in NORMALS["110"]]             # with the 12 {110} planes
MARGIN, FAR = 0.3, 60.0          # A above the outermost lattice plane kept; distance of a plane that carves nothing
OFF = 10 ** 6                    # plane level of a plane that carves nothing


# --- Cd/Se framework (overlay of scan dots) ------------------------------------------------------------------------

def grid(xyz):
    """Cd and Se of a scan structure (rows [symbol, role, x, y, z]) on the a/4 grid from the first Se."""
    at = [a for a in xyz if a[0] in ("Cd", "Se")]
    sym = np.array([a[0] for a in at])
    P = np.array([a[2:5] for a in at], float)
    X = (P - P[list(sym).index("Se")]) / Q
    G = np.rint(X).astype(int)
    assert np.abs(X - G).max() < 0.05, np.abs(X - G).max()
    return sym, G


def keyset(sym, G):
    return {(s, *g) for s, g in zip(sym, map(tuple, G))}


def embeddings(ps, pG, cs, cG, anion="Se"):
    """All (M, t, parent in the child's frame) with M a Td op and t a lattice translation that put every parent site
    on a child site of the same species."""
    ck = keyset(cs, cG)
    out = []
    for M in TD_OPS:
        PM = pG @ M.T
        anchor = PM[ps == anion][0]
        for t in cG[cs == anion] - anchor:
            if all((s, *g) in ck for s, g in zip(ps, map(tuple, PM + t))):
                out.append((M, t, PM + t))
    return out


def canon_cdse(sym, G):
    """Canonical Cd/Se framework (Se the anchor species) up to Td ops and lattice translation."""
    best = None
    se = sym == "Se"
    for M in TD_OPS:
        Y = G @ M.T
        t = Y[se][np.lexsort(Y[se].T[::-1])[0]]
        rep = tuple(sorted((s, *map(int, g)) for s, g in zip(sym, Y - t)))
        best = rep if best is None or rep < best else best
    return best


# --- lattice sites ---------------------------------------------------------------------------------------------------

def is_anion_site(g):
    return bool(((g % 2) == 0).all() and g.sum() % 4 == 0)


def is_cd_site(g):
    return is_anion_site(g - 1)


def sites(sym, P, ref):
    """Lattice sites of a dot in a/4 units from the Se `ref` (Cd and anion sites), with what sits on each."""
    out = {}
    for s, p in zip(sym, P):
        x = (p - ref) / Q
        g = np.rint(x).astype(int)
        if np.abs(x - g).max() > 0.12:
            continue
        if s == "Cd" and is_cd_site(g):
            out[("Cd", *g)] = "Cd"
        elif s in ("Se", "Cl") and is_anion_site(g):
            out[("A", *g)] = s
    return out


def from_struct(st):
    """Sites of a scan structure (data.json structs entry), from its first Se."""
    sym = np.array([a[0] for a in st["xyz"]])
    P = np.array([a[2:5] for a in st["xyz"]], float)
    return sites(sym, P, P[list(sym).index("Se")])


def arr(S):
    keys = list(S)
    return np.array([k[0] for k in keys]), np.array([k[1:] for k in keys], int)


def inside(G, planes):
    """Mask of grid points with n.r <= d for every (n, d)."""
    m = np.ones(len(G), bool)
    for n, d in planes:
        m &= G @ np.array(n) <= d
    return m


def canon(sym, G):
    """Class key: the site set up to Td ops + lattice translation (anchor: the lowest anion site)."""
    best = None
    for M in TD_OPS:
        Y = G @ M.T
        a = Y[sym == "A"]
        t = a[np.lexsort(a.T[::-1])[0]]
        rep = tuple(sorted((s, *map(int, g)) for s, g in zip(sym, Y - t)))
        best = rep if best is None or rep < best else best
    return best


def lattice_box(G, pad=10):
    """Every Cd and anion site of the box around G (pad a/4 on each side): (positions, species 'A'/'Cd')."""
    lo, hi = G.min(0) - pad, G.max(0) + pad
    box = np.array(list(itertools.product(*[range(x, y + 1) for x, y in zip(lo, hi)])))
    an = box[((box % 2) == 0).all(1) & (box.sum(1) % 4 == 0)]
    return np.vstack([an, an + 1]), np.array(["A"] * len(an) + ["Cd"] * len(an))


def embed_sites(ps, pG, cs, cG, use=None):
    """First Td op + translation that puts the parent's sites (those in `use`) on the child's; parent in the child's
    frame, or None."""
    use = np.ones(len(ps), bool) if use is None else use
    ck = keyset(cs, cG)
    for M in TD_OPS:
        PM = pG @ M.T
        anchor = PM[ps == "A"][0]
        for t in cG[cs == "A"] - anchor:
            if all((s, *g) in ck for s, g, u in zip(ps, map(tuple, PM + t), use) if u):
                return PM + t
    return None


def central_anion(G, sym):
    """The anion site nearest the centroid of the sites G (the builder's construction origin)."""
    a = G[sym == "A"]
    return a[np.argmin(np.linalg.norm(a - G.mean(0), axis=1))]


def planes_A(LG, levels, O, planes=PLANES):
    """Builder planes (normal, distance in A from the origin O) of the cut LG by `levels` (normal -> d): planes that
    carve nothing go far out (greedy, in the fixed order of `planes`)."""
    m = inside(LG, list(levels.items()))
    active = dict(levels)
    for _f, n in planes:
        trial = {k: (OFF if k == n else v) for k, v in active.items()}
        if inside(LG, list(trial.items())).sum() == m.sum():
            active = trial
    out = []
    for _f, n in planes:
        nv = np.array(n, float)
        h = FAR if active[n] >= OFF else (active[n] - nv @ O) * Q / np.linalg.norm(nv) + MARGIN
        out.append((tuple(n), float(h)))
    return out, active


def key_sites(rk):
    """(species, positions) of a class key string."""
    key = eval(rk, {"np": np})
    return np.array([x[0] for x in key]), np.array([x[1:] for x in key], int)


def class_planes(rk):
    """Builder planes of a class key: its 14 tight levels, origin on its central anion, faceless planes far out."""
    _sym, G = key_sites(rk)
    LG, Ls = lattice_box(G)
    lev = {n: int((G @ np.array(n)).max()) for _f, n in PLANES}
    msk = inside(LG, list(lev.items()))
    return planes_A(LG, lev, central_anion(LG[msk], Ls[msk]))[0]


# --- built dots ----------------------------------------------------------------------------------------------------

def read_xyz(p):
    lines = Path(p).read_text().split("\n")
    at = [l.split() for l in lines[2:2 + int(lines[0])]]
    return np.array([a[0] for a in at]), np.array([[float(v) for v in a[1:4]] for a in at])


def write_xyz(p, sym, x, comment=""):
    Path(p).write_text(f"{len(sym)}\n{comment}\n" +
                       "".join(f"{s} {v[0]:.6f} {v[1]:.6f} {v[2]:.6f}\n" for s, v in zip(sym, x)))


def dot_key(sy, x):
    """Dot identity: every atom on the a/16 grid (Z-type Cl off the lattice included) up to Td ops + translation."""
    G = x / Q
    G = G - G[list(sy).index("Se")]
    G4 = np.rint(G * 4).astype(int)
    best = None
    for M in TD_OPS:
        Y = G4 @ M.T
        se = Y[sy == "Se"]
        t = se[np.lexsort(se.T[::-1])[0]]
        rep = tuple(sorted((s, *v) for s, v in zip(sy, map(tuple, (Y - t).tolist()))))
        best = rep if best is None or rep < best else best
    return best


def se_key(sy, x):
    """Se framework up to Td ops + lattice translation (the lineage's identity; brothers share it)."""
    se = np.asarray(x, float)[np.asarray(sy) == "Se"] / Q
    G = np.rint(se - se[0]).astype(int)
    best = None
    for M in TD_OPS:
        Y = G @ M.T
        t = Y[np.lexsort(Y.T[::-1])[0]]
        rep = tuple(sorted(map(tuple, (Y - t).tolist())))
        best = rep if best is None or rep < best else best
    return best


def formula(sy):
    return "".join(f"{e}{int((np.asarray(sy) == e).sum())}" for e in ("Cd", "Se", "Cl"))
