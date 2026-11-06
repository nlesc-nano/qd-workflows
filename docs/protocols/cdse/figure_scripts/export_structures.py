"""Write every structure shown in the CdSe protocol docs as .xyz, one folder per set, with structures.json for the
3D viewer (added atoms, removed atoms, the views along the facet of interest) and one VMD script per set.

    python export_structures.py OUTDIR      (qd_builder env; reads the same sources as render.py)
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.argv = [sys.argv[0], sys.argv[1]]
import render as R  # noqa: E402  (paths, readers, alignment, the registry table)

OUT = Path(sys.argv[1])
OUT.mkdir(parents=True, exist_ok=True)
SETS = []

# VMD colour ids for the element colours of the docs (Cd orange, Se blue, Cl green, O red, C grey, H white)
VMD_COL = {"Cd": 3, "Se": 0, "Cl": 7, "O": 1, "C": 2, "H": 8}
VDW = {"Cd": 1.58, "Se": 1.90, "Cl": 1.75, "O": 1.52, "C": 1.70, "H": 1.20}


def view_matrix(view, up=(0, 0, 1)):
    """rows x, y, z of the camera frame: z along the view direction (towards the viewer)"""
    return R.rot(view, up).tolist()


def write_xyz(path, sym, pos, comment):
    with open(path, "w") as fh:
        fh.write(f"{len(sym)}\n{comment}\n")
        for s, p in zip(sym, pos):
            fh.write(f"{s:2s} {p[0]:12.6f} {p[1]:12.6f} {p[2]:12.6f}\n")


def formula(sym):
    order = ["Cd", "Se", "Cl", "C", "O", "H"]
    return "".join(f"{e}{sym.count(e)}" for e in order if sym.count(e))


def add_set(key, title, text, items, views, ref_note):
    """items: (name, label, sym, pos, ref) with ref = (sym, pos) of the reference or None"""
    d = OUT / key
    d.mkdir(exist_ok=True)
    rows = []
    cen = None
    for name, label, sym, pos, ref in items:
        sym = list(sym); pos = np.asarray(pos, float)
        new, gone = [], []
        if ref is not None:
            nm = R.is_new(ref[0], ref[1], sym, pos)
            new = [int(i) for i in np.where(nm)[0]]
            miss = R.is_new(sym, pos, ref[0], ref[1])
            gone = [[ref[0][i]] + [round(float(v), 4) for v in ref[1][i]] for i in np.where(miss)[0]]
        f = formula(sym)
        write_xyz(d / f"{name}.xyz", sym, pos, f"{label} | {f} | added vs reference: {len(new)} atoms, removed: {len(gone)}")
        if cen is None:
            cen = pos.mean(0)
        rows.append({"name": name, "label": label, "formula": f, "file": f"{key}/{name}.xyz", "n": len(sym),
                     "new": new, "gone": gone,
                     "sym": sym, "pos": [[round(float(v), 4) for v in p] for p in pos]})
    SETS.append({"key": key, "title": title, "text": text, "ref_note": ref_note, "center": [float(v) for v in cen],
                 "views": [{"name": n, "R": view_matrix(v, u)} for n, v, u in views], "items": rows})
    vmd(d, key, title, rows, views[0])
    print(f"{key}: {len(rows)} structures")


def vmd(d, key, title, rows, view):
    """open every structure of the set as one molecule each, VDW spheres, element colours, camera along the view"""
    Rm = np.array(view_matrix(view[1], view[2]))
    m4 = [[*Rm[i], 0.0] for i in range(3)] + [[0.0, 0.0, 0.0, 1.0]]
    lines = [f"# {title}: vmd -e {key}.tcl  (run from this folder)", "display projection Orthographic",
             "color Display Background white", "axes location Off"]
    for e, c in VMD_COL.items():
        lines.append(f"color Element {e} {c}")
    lines.append("color change rgb 3 0.94 0.63 0.19")
    lines.append("color change rgb 0 0.18 0.44 0.87")
    for r in rows:
        lines.append(f"mol new {r['name']}.xyz type xyz waitfor all")
        lines.append("mol delrep 0 top")
        lines.append("mol representation VDW 1.0 30")
        lines.append("mol color Element")
        lines.append("mol addrep top")
        lines.append(f"mol rename top {{{r['label']}}}")
        for e, v in VDW.items():
            lines.append(f"set s [atomselect top \"name {e}\"]; $s set radius {v}; $s set element {e}; $s delete")
        if r["new"]:
            idx = " ".join(map(str, r["new"]))
            lines.append(f"mol representation VDW 1.05 30\nmol color ColorID 16\nmol material Transparent\nmol selection {{index {idx}}}\nmol addrep top\nmol selection all\nmol material Opaque")
        lines.append("mol off top")
    lines.append("mol on 0")
    mat = "{" + " ".join("{" + " ".join(f"{v:.6f}" for v in row) + "}" for row in m4) + "}"
    lines.append("foreach m [molinfo list] { molinfo $m set rotate_matrix {" + mat + "} }")
    lines.append("display resetview")
    lines.append("# In VMD: Graphics > Representations; Molecule list: toggle the 'D' column to show one structure at a time."
                 " The black outline (rep 2) marks the atoms added against the set's reference.")
    (d / f"{key}.tcl").write_text("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------------------------------------------
def grid_sets():
    s0, p0 = R.xyz(R.GRID / "start/start.xyz")
    s1, p1 = R.xyz(R.GRID / "start/seed_x.xyz")
    row = {int(R.kv(c, "k")): (c, sy, po) for c, sy, po in R.grid_rows(R.GRID / "xyz/m13.xyz")}
    add_set("grid_seeds", "Grid: start dot and seeds",
            "The builder's start dot (13, 15) and X-type seed (13, 3), and the lowest relaxed state of the same two grid "
            "nodes (CdCl2 grid cdse_full2).",
            [("start_dot_13_15_builder", "start dot (13, 15), builder", s0, p0, None),
             ("seed_x_13_3_builder", "X-type seed (13, 3), builder", s1, p1, None),
             ("grid_13_3_relaxed", "grid node (13, 3), relaxed", *row[3][1:], None),
             ("grid_13_15_relaxed", "grid node (13, 15), relaxed", *row[15][1:], None)],
            [("along [111]", (1, 1, 1), (0, 0, 1)), ("along [100]", (1, 0, 0), (0, 0, 1)), ("along [110]", (1, 1, 0), (0, 0, 1))],
            "no reference")
    items = [(f"m13_k{k:02d}", f"(13, {k}) {R.kv(c, 'source')}", sy, po, None) for k, (c, sy, po) in sorted(row.items())]
    add_set("grid_row13", "Grid: the m = 13 row",
            "Lowest relaxed state of every node of the m = 13 row of the CdCl2 grid.", items,
            [("along [111]", (1, 1, 1), (0, 0, 1)), ("along [100]", (1, 0, 0), (0, 0, 1))], "no reference")
    items = []
    for m in range(13, 0, -1):
        rows = R.grid_rows(R.GRID / f"xyz/m{m:02d}.xyz")
        best = min(rows, key=lambda r: abs(int(R.kv(r[0], "k")) - 3))
        items.append((f"m{m:02d}_k{int(R.kv(best[0], 'k')):02d}", f"({m}, {R.kv(best[0], 'k')})", best[1], best[2], None))
    add_set("grid_descent", "Grid: removing CdSe units (k ≈ 3)",
            "Lowest state of the node nearest to k = 3 in each row, m = 13 down to 1.", items,
            [("along [111]", (1, 1, 1), (0, 0, 1)), ("along [100]", (1, 0, 0), (0, 0, 1))], "no reference")


def builder_set():
    d = R.BUILDS / "c000015"
    cut = R.xyz(d / "out_cut.xyz"); pre = R.xyz(d / "out_02_after_prepass.xyz"); fin = R.xyz(d / "out.xyz")
    add_set("builder_stages", "Builder: the seed through the builder stages",
            "Wulff cut, after the prepass, charge balanced (Q = 0). Added/removed atoms are against the previous stage.",
            [("1_wulff_cut", "1. Wulff cut (lattice sites)", *cut, None),
             ("2_after_prepass", "2. after the prepass", *pre, cut),
             ("3_charge_balanced", "3. charge balanced Cd16Se13Cl6", *fin, pre)],
            [("along [111]", (1, 1, 1), (0, 0, 1)), ("along [100]", (1, 0, 0), (0, 0, 1))], "previous stage")


def layer_set():
    B = R.NOTES1 / "layer_path_Cd38Se28Cl20/builds"
    st = {x["i"]: x for x in json.loads((B / "states.json").read_text())["states"]}
    seed = R.xyz(B / "s00/out.xyz")
    path = [0, 1, 5, 10, 18, 21, 26, 27]
    items, prev = [], None
    for i in path:
        s, p = R.xyz(B / f"s{i:02d}/out.xyz")
        mv = "seed" if prev is None else "+ " + (set(st[i]["moved"]) - set(st[prev[2]]["moved"])).pop()
        items.append((f"path_{len(items)}_s{i:02d}", f"s{i:02d} {mv}", s, p, None if prev is None else prev[:2]))
        prev = (s, p, i)
    add_set("layer_path", "Layer path: seed → Cd38Se28Cl20, one plane per step",
            "One path through the intermediates of the lineage edge seed → Cd38Se28Cl20. Added atoms are against the "
            "previous step.", items,
            [("oblique", (1, 0.35, 0.55), (0, 0, 1)), ("along [100]", (1, 0, 0), (0, 0, 1)),
             ("along [111]", (1, 1, 1), (0, 0, 1)), ("along [1-1-1]", (1, -1, -1), (0, 0, 1))], "previous step")
    items = []
    for i in sorted(st):
        s, p = R.xyz(B / f"s{i:02d}/out.xyz")
        items.append((f"s{i:02d}", f"s{i:02d} {st[i]['formula']} | {len(st[i]['moved'])} layers: {', '.join(st[i]['moved']) or 'seed'}",
                      s, p, None if i == 0 else seed))
    add_set("layer_states", "Layer path: all 28 intermediates (vs the seed)",
            "All 28 symmetry-distinct layer states of the edge; added atoms are against the seed.", items,
            [("oblique", (1, 0.35, 0.55), (0, 0, 1)), ("along [100]", (1, 0, 0), (0, 0, 1)), ("along [111]", (1, 1, 1), (0, 0, 1))],
            "the seed")


def island_sets():
    si = json.loads((R.ISL.parent / "states_index.json").read_text())
    for key, pid, pdir, fam, n0 in [("islands_111B", 1817, "c002306", "{111}B", (1, -1, -1)),
                                    ("islands_111A", 6281, "c009066", "{111}A", (1, 1, -1)),
                                    ("islands_100", 3592, "c005070", "{100}", (0, 1, 0))]:
        pa = R.xyz(R.BUILDS / pdir / "out.xyz")
        rows, seen = [], set()
        for x in sorted(si, key=lambda x: (x.get("dCd") or 0, x["tag"])):
            if x["parent"] != pid or x["family"] not in fam.replace("{", "").replace("}", "") or x.get("Q", 0) != 0:
                continue
            if list(x["n0"]) != list(n0) or not x.get("pure") or x.get("dup_of") is not None or (x.get("dCd") or 0) < 1:
                continue
            if x["tag"] in seen or not (R.ISL / x["tag"] / "out.xyz").exists():
                continue
            seen.add(x["tag"])
            rows.append((x["tag"], f"+{x['dCd']} Cd ({x['formula']})", *R.xyz(R.ISL / x["tag"] / "out.xyz"), pa))
        side = tuple(np.cross(n0, (0, 0, 1) if abs(n0[2]) < 0.9 * np.linalg.norm(n0) else (1, 0, 0)))
        add_set(key, f"Islands on {fam}: parent {R.REG[pdir]}",
                f"Builder islands on one {fam} facet of {R.REG[pdir]} (registry {pid}), sorted by the number of new Cd. "
                "Added atoms are against the parent.",
                [("parent", f"parent {R.REG[pdir]}", *pa, None)] + rows,
                [(f"down the {fam} normal", n0, (0, 0, 1) if abs(n0[2]) < 0.9 * np.linalg.norm(n0) else (1, 0, 0)),
                 ("side view, facet up", side, n0)], "the parent")


def shells_set():
    base = R.xyz(R.BUILDS / "c000036/out.xyz")
    items = [("k0", "k0 (builder) Cd25Se22Cl6", *base, None)]
    for lab, d in [("k0-2", "fp00033_k-2"), ("k0-1", "fp00033_k-1"), ("k0+1", "fp00033_k+1"), ("k0+2", "fp00033_k+2")]:
        b = R.align_lattice(base, R.xyz(R.BUILDS / d / "out.xyz"))
        items.append((lab.replace("+", "p").replace("-", "m"), f"{lab} {R.REG[d]}", *b, base))
    add_set("cd_shells", "Cd shells k0 ± 1, ± 2 (z_type_target)",
            "A fill dot and its shell variants, superposed on the k0 build. Added: Z-type CdCl2 units; removed (k < k0): "
            "displaced units.", items,
            [("along [111]", (1, 1, 1), (0, 0, 1)), ("along [100]", (1, 0, 0), (0, 0, 1))], "the k0 build")


def acetate_set():
    ds = json.loads((R.SCR / "doc_structs.json").read_text())
    items = []
    for name in ("c000015", "c000036"):
        s0, p0 = R.xyz(R.BUILDS / name / "out.xyz")
        s = ds[name]["symbols"]; st = np.array(ds[name]["start_positions"]); pr = R.kabsch(np.array(ds[name]["positions"]), st)
        items += [(f"{name}_cdcl2_builder", f"{R.REG[name]} CdCl2 builder", s0, p0, None),
                  (f"{name}_acetate_dressed", f"{name} acetate dressed (start)", s, st, None),
                  (f"{name}_acetate_relaxed", f"{name} acetate relaxed (MACE-MH-1 + D3)", s, pr, None)]
    add_set("acetate_dressing", "Acetate dressing of two CdCl2 builder dots",
            "CdCl2 builder dot, the same dot dressed with acetate by the builder's X-type exchange, and relaxed.", items,
            [("along [111]", (1, 1, 1), (0, 0, 1)), ("along [100]", (1, 0, 0), (0, 0, 1))], "no reference")


if __name__ == "__main__":
    grid_sets(); builder_set(); layer_set(); island_sets(); shells_set(); acetate_set()
    (OUT / "structures.json").write_text(json.dumps({"sets": SETS}, separators=(",", ":")))
    readme = ["# Structures of the CdSe protocol docs", "",
              "One folder per set; every structure is a plain .xyz (Å). The comment line gives the label, the formula and "
              "how many atoms were added or removed against the set's reference.", "",
              "- Open a set in VMD: `cd <set>; vmd -e <set>.tcl` (VDW spheres, Cd orange, Se blue, camera along the "
              "facet of interest; added atoms get a black transparent shell).",
              "- Or open any .xyz in VESTA, Avogadro or OVITO.",
              "- The interactive viewer (rotate, views along each facet, added atoms highlighted): `viewer.html` next to "
              "this folder, or the Artifact page.", "", "| set | structures | reference for 'added' |", "|---|---|---|"]
    for s in SETS:
        readme.append(f"| `{s['key']}` — {s['title']} | {len(s['items'])} | {s['ref_note']} |")
    (OUT / "README.md").write_text("\n".join(readme) + "\n")
    print("total", sum(len(s["items"]) for s in SETS), "structures")
