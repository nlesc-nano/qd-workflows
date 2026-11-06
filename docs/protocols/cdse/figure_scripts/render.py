"""Structure figures for the CdSe protocol docs: depth-sorted circle projections (Cd orange, Se blue as on the
dashboard), the atoms a step adds outlined in black and everything else faded."""
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from ase.io import read

OUT = Path(sys.argv[1])
OUT.mkdir(parents=True, exist_ok=True)
SCR = Path("/private/tmp/claude-501/-Users-ivaninfante-Documents-University-escience-qd-workflows--claude-worktrees-confident-davinci-e3556a/c7860e0d-ca2f-4ef8-9d28-226f755fcca2/scratchpad")
GS = Path("/private/tmp/claude-501/-Users-ivaninfante-Documents-University-escience-qd-workflows--claude-worktrees-growth-models-discussion-cf78cd/1f382c86-e331-4984-af42-73bef9e1580d/scratchpad")
NOTES1 = Path("/Users/ivaninfante/Library/Application Support/Claude/side-session-notes/local_61fddce9-c61a-4268-b693-c539c6f85e24")
GRID = Path("/Users/ivaninfante/Documents/University/escience/qd-workflows/.claude/worktrees/cdse-grid/cd16_lattice_grid/nucleation_builder_main/full_run2")
ISL = Path("/Users/ivaninfante/Documents/University/escience/network_growth_data/island_builds_cdcl2/builds")
BUILDS = GS / "overlay/net/builds"

COL = {"Cd": "#f0a030", "Se": "#2f6fdf", "Cl": "#3aa655", "O": "#d9443a", "C": "#6b6b6b", "H": "#e8e8e8"}
RAD = {"Cd": 1.45, "Se": 1.55, "Cl": 1.15, "O": 0.85, "C": 0.8, "H": 0.45}
ROLE = {"facet": "#3f8fd2", "edge": "#8e44ad", "vertex": "#d0453a", "interior": "#9a9a9a"}


def rot(view, up=(0, 0, 1)):
    z = np.asarray(view, float); z /= np.linalg.norm(z)
    u = np.asarray(up, float)
    if abs(np.dot(u, z)) > 0.95:
        u = np.array([0, 1, 0.0]) if abs(z[1]) < 0.9 else np.array([1, 0, 0.0])
    x = np.cross(u, z); x /= np.linalg.norm(x)
    y = np.cross(z, x)
    return np.vstack([x, y, z])


def mix(c, w):
    c = np.array(matplotlib.colors.to_rgb(c))
    return tuple(c + (1 - c) * w)


def draw(ax, sym, pos, view=(1, 1, 1), up=(0, 0, 1), new=None, colors=None, title=None, center=None, span=None, ghost=None):
    R = rot(view, up)
    c = pos.mean(0) if center is None else center
    p = (pos - c) @ R.T
    order = np.argsort(p[:, 2])
    zmin, zmax = p[:, 2].min(), p[:, 2].max()
    for i in order:
        s = sym[i]
        depth = 0 if zmax == zmin else (zmax - p[i, 2]) / (zmax - zmin)
        base = colors[i] if colors is not None else COL.get(s, "#aaaaaa")
        if new is not None and not new[i]:
            fc, ec, lw = mix(base, 0.62 + 0.15 * depth), mix("#444444", 0.6), 0.4
        else:
            fc, ec, lw = mix(base, 0.25 * depth), ("#111111" if new is not None else mix("#222222", 0.3)), (1.3 if new is not None else 0.5)
        ax.add_patch(Circle((p[i, 0], p[i, 1]), RAD.get(s, 1.0) * 0.82, fc=fc, ec=ec, lw=lw))
    if ghost is not None and len(ghost[0]):
        gq = (ghost[1] - c) @ R.T
        for gs, gp in zip(ghost[0], gq):
            ax.add_patch(Circle((gp[0], gp[1]), RAD.get(gs, 1.0) * 0.82, fc="none", ec="#c0392b", lw=1.2, ls="--"))
    L = span or (np.abs(p[:, :2]).max() + 2.0)
    ax.set_xlim(-L, L); ax.set_ylim(-L, L); ax.set_aspect("equal"); ax.axis("off")
    if title:
        ax.set_title(title, fontsize=9)
    return L


def is_new(sym_a, pos_a, sym_b, pos_b, tol=0.35):
    new = np.ones(len(sym_b), bool)
    for i, (s, x) in enumerate(zip(sym_b, pos_b)):
        d = np.linalg.norm(pos_a - x, axis=1)
        j = np.where((d < tol) & (np.array(sym_a) == s))[0]
        new[i] = len(j) == 0
    return new


def kabsch(mov, ref):
    a, b = mov - mov.mean(0), ref - ref.mean(0)
    u, _, vt = np.linalg.svd(a.T @ b)
    d = np.sign(np.linalg.det(u @ vt))
    r = u @ np.diag([1, 1, d]) @ vt
    return (mov - mov.mean(0)) @ r + ref.mean(0)


def cubic_rotations():
    import itertools
    out = []
    for perm in itertools.permutations(range(3)):
        for sg in itertools.product((1, -1), repeat=3):
            R = np.zeros((3, 3))
            for i, j in enumerate(perm):
                R[i, j] = sg[i]
            if np.linalg.det(R) > 0:
                out.append(R)
    return out


def align_lattice(ref, mov, el="Se"):
    """rotate/translate mov (sym, pos) so its Se sites coincide with ref's (24 proper cubic rotations)"""
    ra = ref[1][np.array(ref[0]) == el]; mb = mov[1][np.array(mov[0]) == el]
    best = (-1, None)
    for R in cubic_rotations():
        q = mb @ R.T
        for j in range(len(ra)):
            t = ra[j] - q[0]
            d = np.min(np.linalg.norm((q + t)[:, None, :] - ra[None, :, :], axis=2), axis=1)
            n = int(np.sum(d < 0.4))
            if n > best[0]:
                best = (n, (R, t))
    R, t = best[1]
    return mov[0], mov[1] @ R.T + t


def xyz(path, index=-1):
    a = read(str(path), index=index)
    return list(a.get_chemical_symbols()), a.get_positions()


def strip(formula):
    return formula


def save(fig, name):
    fig.savefig(OUT / name, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("wrote", name)


def legend(fig, extra=None, y=0.0):
    items = [("Cd", COL["Cd"]), ("Se", COL["Se"]), ("Cl", COL["Cl"])] + (extra or [])
    hs = [plt.Line2D([], [], marker="o", ls="", mfc=c, mec="#333", ms=8, label=l) for l, c in items]
    fig.legend(handles=hs, loc="lower center", ncol=len(hs), frameon=False, fontsize=8, bbox_to_anchor=(0.5, y))


def grid_rows(path):
    """frames and comment lines of a grid row file m??.xyz"""
    lines = Path(path).read_text().splitlines()
    out, i = [], 0
    while i < len(lines):
        n = int(lines[i]); com = lines[i + 1]
        at = [l.split() for l in lines[i + 2:i + 2 + n]]
        sym = [a[0] for a in at]; pos = np.array([[float(v) for v in a[1:4]] for a in at])
        keep = [j for j, s in enumerate(sym) if s != "X" and s != "Xx"]
        out.append((com, [sym[j] for j in keep], pos[keep]))
        i += 2 + n
    return out


def kv(com, key):
    for t in com.split():
        if t.startswith(key + "="):
            return t.split("=", 1)[1]


# 1. start dot and seeds -------------------------------------------------------------------------------------------
def fig_seeds():
    fig, axs = plt.subplots(1, 4, figsize=(12, 3.4))
    s, p = xyz(GRID / "start/start.xyz"); draw(axs[0], s, p, title="start dot, builder\nCd28Se13Cl30 (13, 15)")
    s, p = xyz(GRID / "start/seed_x.xyz"); draw(axs[1], s, p, title="X-type seed, builder\nCd16Se13Cl6 (13, 3)")
    row = {int(kv(c, "k")): (sy, po) for c, sy, po in grid_rows(GRID / "xyz/m13.xyz")}
    draw(axs[2], *row[3], title="grid node (13, 3)\nrelaxed lowest state")
    draw(axs[3], *row[15], title="grid node (13, 15)\nrelaxed (annealed) lowest state")
    legend(fig, y=-0.04)
    save(fig, "grid_seeds.png")


def fig_grid_row():
    rows = grid_rows(GRID / "xyz/m13.xyz")
    ks = [0, 3, 6, 9, 12, 15]
    row = {int(kv(c, "k")): (c, sy, po) for c, sy, po in rows}
    ks = [k for k in ks if k in row]
    fig, axs = plt.subplots(1, len(ks), figsize=(2.3 * len(ks), 2.8))
    for ax, k in zip(axs, ks):
        c, sy, po = row[k]
        draw(ax, sy, po, title=f"(13, {k})  {kv(c, 'formula')}\n{kv(c, 'source')}")
    legend(fig, y=-0.06)
    save(fig, "grid_row13.png")


def fig_grid_descent():
    ms = [13, 11, 9, 7, 5, 3, 2, 1]
    fig, axs = plt.subplots(1, len(ms), figsize=(2.0 * len(ms), 2.6))
    for ax, m in zip(axs, ms):
        rows = grid_rows(GRID / f"xyz/m{m:02d}.xyz")
        best = min(rows, key=lambda r: abs(int(kv(r[0], "k")) - 3))
        draw(ax, best[1], best[2], title=f"({m}, {kv(best[0], 'k')})\n{kv(best[0], 'formula')}", span=9.5)
    legend(fig, y=-0.08)
    save(fig, "grid_descent.png")


# 2. builder stages and site roles ---------------------------------------------------------------------------------
def fig_builder():
    d = BUILDS / "c000015"
    fig, axs = plt.subplots(1, 4, figsize=(12, 3.4))
    for ax, f, t in zip(axs, ["out_cut.xyz", "out_02_after_prepass.xyz", "out.xyz"],
                         ["1. Wulff cut (lattice sites)", "2. after the prepass\n(low-CN Se swapped to Cl)", "3. charge balanced, Q = 0\nCd16Se13Cl6"]):
        s, p = xyz(d / f); draw(ax, s, p, title=t)
    s, p = xyz(d / "out.xyz")
    man = json.loads((d / "out.json").read_text())["site_roles"]
    roles = man["roles"]
    cols = [ROLE[r.split(":")[0]] for r in roles]
    draw(axs[3], s, p, colors=cols, title="4. builder site roles\n(manifest site_roles)")
    hs = [plt.Line2D([], [], marker="o", ls="", mfc=c, mec="#333", ms=8, label=l) for l, c in ROLE.items()]
    fig.legend(handles=[plt.Line2D([], [], marker="o", ls="", mfc=COL[e], mec="#333", ms=8, label=e) for e in ("Cd", "Se", "Cl")] + hs,
               loc="lower center", ncol=7, frameon=False, fontsize=8, bbox_to_anchor=(0.5, -0.04))
    save(fig, "builder_stages.png")


# 3. layer path: seed -> Cd38Se28Cl20 ------------------------------------------------------------------------------
def fig_layer_path():
    B = NOTES1 / "layer_path_Cd38Se28Cl20/builds"
    st = {x["i"]: x for x in json.loads((B / "states.json").read_text())["states"]}
    path = [0, 1, 5, 10, 18, 21, 26, 27]
    fig, axs = plt.subplots(2, 4, figsize=(12, 6.6))
    prev = None
    ref = xyz(B / "s27/out.xyz")
    cen = ref[1].mean(0)
    for ax, i in zip(axs.flat, path):
        s, p = xyz(B / f"s{i:02d}/out.xyz")
        new = None if prev is None else is_new(prev[0], prev[1], s, p)
        mv = "seed" if i == 0 else "+ " + (set(st[i]["moved"]) - set(st[prev_i]["moved"])).pop()
        draw(ax, s, p, view=(1, 0.35, 0.55), new=new, center=cen, span=9.0, title=f"s{i:02d}  {st[i]['formula']}\n{mv}")
        prev, prev_i = (s, p), i
    legend(fig, [("added by this step (outlined)", "#ffffff")], y=0.0)
    save(fig, "layer_path.png")


def fig_parent_child():
    B = NOTES1 / "layer_path_Cd38Se28Cl20/builds"
    a = xyz(B / "s00/out.xyz"); b = xyz(B / "s27/out.xyz")
    cen = b[1].mean(0)
    fig, axs = plt.subplots(1, 3, figsize=(10, 3.6))
    draw(axs[0], *a, view=(1, 0.35, 0.55), center=cen, span=9.0, title="parent: seed Cd16Se13Cl6 (13, 3)")
    draw(axs[1], *b, view=(1, 0.35, 0.55), center=cen, span=9.0, title="child: Cd38Se28Cl20 (28, 10)")
    draw(axs[2], *b, view=(1, 0.35, 0.55), center=cen, span=9.0, new=is_new(a[0], a[1], *b),
         title="child, with the parent's atoms faded\n7 layers added on 7 facet planes")
    legend(fig, [("not in the parent (outlined)", "#ffffff")], y=-0.03)
    save(fig, "parent_child.png")


# 4. islands -------------------------------------------------------------------------------------------------------
def fig_islands(parent_dir, fam, n0, tags, name, up=(0, 0, 1)):
    pa = xyz(BUILDS / parent_dir / "out.xyz")
    n = len(tags) + 1
    fig, axs = plt.subplots(2, n, figsize=(2.4 * n, 5.2))
    cen = pa[1].mean(0)
    side = np.cross(n0, up if abs(np.dot(np.array(n0) / np.linalg.norm(n0), up)) < 0.9 else (1, 0, 0))
    span = None
    for j, (lab, path) in enumerate([("parent " + pa_formula(parent_dir), None)] + tags):
        b = pa if path is None else xyz(ISL / path / "out.xyz")
        new = None if path is None else is_new(pa[0], pa[1], *b)
        span = draw(axs[0, j], *b, view=n0, up=up, new=new, center=cen, span=span, title=lab)
        draw(axs[1, j], *b, view=side, up=n0, new=new, center=cen, span=span)
    axs[0, 0].text(-0.12, 0.5, f"top view\n(along {fam} normal)", transform=axs[0, 0].transAxes, rotation=90, va="center", ha="center", fontsize=8)
    axs[1, 0].text(-0.12, 0.5, "side view\n(facet up)", transform=axs[1, 0].transAxes, rotation=90, va="center", ha="center", fontsize=8)
    legend(fig, [("island atoms (outlined)", "#ffffff")], y=0.0)
    save(fig, name)


REG = {}
for l in (SCR / "reg_dots.txt").read_text().splitlines():
    t = l.split()
    if len(t) >= 3:
        REG[t[2]] = t[1]


def pa_formula(d):
    return REG.get(d, d)


# 5. Cd shells -----------------------------------------------------------------------------------------------------
def fig_shells():
    base = xyz(BUILDS / "c000036/out.xyz")
    items = [("k0 − 2", "fp00033_k-2"), ("k0 − 1", "fp00033_k-1"), ("k0 (builder)", None), ("k0 + 1", "fp00033_k+1"), ("k0 + 2", "fp00033_k+2")]
    fig, axs = plt.subplots(1, 5, figsize=(12, 3.0))
    cen = base[1].mean(0)
    for ax, (lab, d) in zip(axs, items):
        b = base if d is None else align_lattice(base, xyz(BUILDS / d / "out.xyz"))
        new = None if d is None else is_new(base[0], base[1], *b)
        gone = None
        if d is not None:
            miss = is_new(b[0], b[1], *base)
            gone = ([x for x, f_ in zip(base[0], miss) if f_], base[1][miss])
        f = REG.get(d or "c000036")
        draw(ax, *b, new=new, center=cen, span=11, title=f"{lab}\n{f}", ghost=gone)
    legend(fig, [("added to the k0 build (outlined)", "#ffffff")], y=-0.06)
    fig.text(0.5, -0.1, "dashed red: atoms of the k0 build removed (a CdCl2 unit = one Cd and its two Cl)", ha="center", fontsize=8, color="#c0392b")
    save(fig, "cd_shells.png")


# 6. acetate dressing ----------------------------------------------------------------------------------------------
def fig_acetate():
    ds = json.loads((SCR / "doc_structs.json").read_text())
    fig, axs = plt.subplots(2, 3, figsize=(10, 6.6))
    for r, (name, d) in enumerate([("c000015", "c000015"), ("c000036", "c000036")]):
        s0, p0 = xyz(BUILDS / d / "out.xyz")
        s = ds[name]["symbols"]; st = np.array(ds[name]["start_positions"]); pr = kabsch(np.array(ds[name]["positions"]), st)
        cen = p0.mean(0)
        L = draw(axs[r, 0], s0, p0, center=cen, title=f"CdCl2 builder dot\n{REG[d]}")
        draw(axs[r, 1], s, st, center=st[[i for i, x in enumerate(s) if x in ("Cd", "Se")]].mean(0), span=L + 2,
             title="acetate dressed (builder\nX-type exchange), before relaxing")
        draw(axs[r, 2], s, pr, center=pr[[i for i, x in enumerate(s) if x in ("Cd", "Se")]].mean(0), span=L + 2,
             title="relaxed, MACE-MH-1 + D3(BJ)")
    legend(fig, [("O", COL["O"]), ("C", COL["C"]), ("H", COL["H"])], y=0.0)
    save(fig, "acetate_dressing.png")


if __name__ == "__main__":
    fig_seeds(); fig_grid_row(); fig_grid_descent(); fig_builder(); fig_layer_path(); fig_parent_child()
    fig_islands("c002306", "{111}B", (1, -1, -1),
                [(f"+{d} Cd", t) for d, t in [(3, "p1817_111B_pmm_s9_10_14_12"), (5, "p1817_111B_pmm_s11_10_14_14"),
                                               (8, "p1817_111B_pmm_s12_12_12_16"), (10, "p1817_111B_pmm_s13_12_14_16"),
                                               (12, "p1817_111B_pmm_s14_14_14_16")]], "islands_111B.png")
    fig_islands("c009066", "{111}A", (1, 1, -1),
                [(f"+{d} Cd", t) for d, t in [(1, "p6281_111A_ppm_s1_20_22_26"), (2, "p6281_111A_ppm_s2_20_24_26"),
                                               (4, "p6281_111A_ppm_s4_22_24_26"), (6, "p6281_111A_ppm_s6_24_24_26")]], "islands_111A.png")
    fig_islands("c005070", "{100}", (0, 1, 0),
                [(f"+{d} Cd", t) for d, t in [(5, "p3592_100_ppp_s5_18_10_12_14"), (6, "p3592_100_ppp_s9_19_11_13_16")]], "islands_100.png")
    fig_shells(); fig_acetate()
