"""
Monomer library: the MA and MX_q monomers (CdSe, CdCl2, InP, InCl3, ...) computed once per formula and
shared by every dot, whatever its material or phase (CdCl2 is the same for CdS, CdSe and CdTe dots,
zinc blende or wurtzite). Part of stage 1 of the campaign.

    MACE-MH-1 (GPU job): Orchestr.AI's monomer library, <props.refs>/<head>/monomers/<name>-<key>/:
                         geometry, energy, frequencies, G(T); the training frames (--ensemble).
    CP2K (CPU job):      <cp2k.monomer_library>/<name>/hirshfeld  zero-field single point (Hirshfeld
                                                                  charges, the field step's mono_<label>)
                         <cp2k.monomer_library>/<name>/esp        potential and density cubes (the ESP
                                                                  step's <label>, COSMO solvation)
                         Small molecules: many runs side by side on one node, a few ranks each.

The DFT chain of a dot links its field/mono_<label> and esp/<label> to the library entry
(`qdw monomers link`) when the entry has finished and its input and geometry are the ones the chain
would run; otherwise the chain runs the monomer itself.

    qdw monomers species CdSe InP PbS ...         # the monomer names of these materials
    qdw monomers prepare <dir> --cluster nhr --materials CdSe,InP,...   # gpu.sh, cpu.sh, submit.sh
    qdw monomers status [--cluster nhr]
    qdw monomers link <field dir> <esp dir> --monomers props/monomers.xyz --library <lib>
"""
from __future__ import annotations

import json
import re
import shutil
import tempfile
from pathlib import Path

import numpy as np

# Formal cation charges of the library's families (II-VI, III-V, IV-VI).
CATION_CHARGE = {"Cd": 2, "Zn": 2, "Hg": 2, "Pb": 2, "Sn": 2, "Al": 3, "Ga": 3, "In": 3}
LIGAND = "Cl"
SIDE_BY_SIDE_RANKS = 8          # ranks per monomer run (small molecules: a whole node would be wasted)
GEOMETRY_TOL_A = 1e-5


# ---------------- names ----------------
def split_material(material: str) -> tuple[str, str]:
    """'CdSe' -> ('Cd', 'Se')."""
    parts = re.findall(r"[A-Z][a-z]?", material)
    if len(parts) != 2 or "".join(parts) != material:
        raise ValueError(f"not a binary material name: {material}")
    return parts[0], parts[1]


def species(materials, ligand: str = LIGAND) -> list[dict]:
    """The monomers of these materials, each once: MA (CdSe) and MX_q (CdCl2), with the label the
    DFT chain gives them (MA, MX)."""
    out: dict[str, dict] = {}
    for mat in materials:
        m, a = split_material(mat)
        q = CATION_CHARGE[m]
        out.setdefault(mat, {"name": mat, "label": "MA", "materials": []})["materials"].append(mat)
        mx = f"{m}{ligand}{q if q > 1 else ''}"
        out.setdefault(mx, {"name": mx, "label": "MX", "materials": []})["materials"].append(mat)
    return list(out.values())


# ---------------- inputs ----------------
def read_frames(path) -> list[tuple[str, str, list, np.ndarray]]:
    """(label, name, symbols, positions) of an xyz file of monomers ("MA CdSe" / "MX CdCl2"; a
    comment with one word is both label and name)."""
    lines, out, i = Path(path).read_text().splitlines(), [], 0
    while i < len(lines) and lines[i].strip():
        n = int(lines[i])
        words = lines[i + 1].split()
        rows = [ln.split() for ln in lines[i + 2:i + 2 + n]]
        out.append((words[0], words[1] if len(words) > 1 else words[0], [r[0] for r in rows],
                    np.array([[float(x) for x in r[1:4]] for r in rows])))
        i += n + 2
    return out


_KIND = re.compile(r"^\s*&KIND\s+(\S+).*?^\s*&END(?:\s+KIND)?\s*$\n?", re.M | re.S | re.I)


def normalized(inp: str, elements) -> str:
    """A CP2K input as compared between a dot's chain and the library: no comments or blank lines,
    and only the &KIND blocks of `elements` (a dot's base input has the kinds of all its elements)."""
    keep = set(elements)
    t = _KIND.sub(lambda m: m.group(0) if m.group(1) in keep else "", inp)
    rows = [re.sub(r"\s+", " ", ln.split("#")[0].split("!")[0]).strip() for ln in t.splitlines()]
    return "\n".join(r for r in rows if r)


def _xyz(sym, pts, comment: str) -> str:
    return f"{len(sym)}\n{comment}\n" + "".join(f"{s} {x:.8f} {y:.8f} {z:.8f}\n"
                                               for s, (x, y, z) in zip(sym, np.asarray(pts, float)))


def write_cp2k(library: str, monomers_xyz: str, cluster: dict) -> list[str]:
    """The library's CP2K runs from the MACE geometries (one xyz frame per monomer, "<name> <name>"
    or "<label> <name>"): <library>/<name>/hirshfeld and <library>/<name>/esp, each with cp2k.inp and
    geom.xyz written exactly as a dot's chain writes its monomer runs. Returns the run directories."""
    from orchestr_ai.qd.field import monomer_cell, monomer_input
    from .cp2k import esp
    from .cp2k.inputs import cell_abc, prepare
    lib = Path(library)
    runs = []
    for label, name, sym, pts in read_frames(monomers_xyz):
        lab = label if label in ("MA", "MX") else _label_of(name)
        d = lib / name
        runs += [str(d / "hirshfeld"), str(d / "esp")]
        if run_done(d / "hirshfeld") and run_done(d / "esp"):
            continue                           # finished: never rewritten
        d.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp) / "m.xyz").write_text(_xyz(sym, pts, name))
            # (the restart line is removed from both monomer inputs: they start from an atomic guess)
            prepare(str(Path(tmp) / "m.xyz"), "sample", f"{name}-base", str(Path(tmp) / "base"), cluster,
                    wfn="../../geo_opt/none-RESTART.wfn")
            base = (Path(tmp) / "base" / "cp2k.inp").read_text()
        h = d / "hirshfeld"
        h.mkdir(exist_ok=True)
        (h / "cp2k.inp").write_text(monomer_input(base, f"mono_{lab}", monomer_cell(base, pts)))
        (h / "geom.xyz").write_text(_xyz(sym, pts, lab))
        cell = cell_abc(pts, min_length=esp.MONOMER_MIN_CELL_A)
        esp._write(d / "esp", esp.esp_input(base, lab, scratch=True, cell=cell), sym, pts, f"{lab} monomer (MACE)")
        (d / "monomer.json").write_text(json.dumps({"name": name, "label": lab, "symbols": list(sym),
                                                    "positions": np.asarray(pts, float).tolist()}, indent=1))
    return runs


def _label_of(name: str) -> str:
    """MX for a cation with ligands (CdCl2, InCl3), MA otherwise."""
    m, a = re.findall(r"[A-Z][a-z]?", name)[:2]
    return "MX" if a == LIGAND else "MA"


# ---------------- state ----------------
def run_done(d: Path) -> bool:
    """A library run finished: CP2K ended; the ESP run also wrote both cubes."""
    out = d / "cp2k.out"
    if not out.is_file():
        return False
    t = out.read_text(errors="ignore")
    if "PROGRAM ENDED" not in t[-4000:] or "SCF run NOT converged" in t:
        return False
    return d.name != "esp" or (any(d.glob("v_hartree*")) and any(d.glob("density*")))


def status(library: str) -> dict:
    lib = Path(library)
    res = {}
    for d in sorted(p for p in lib.iterdir() if (p / "monomer.json").is_file()) if lib.is_dir() else []:
        res[d.name] = {k: run_done(d / k) for k in ("hirshfeld", "esp")}
    return {"library": str(lib), "monomers": res,
            "done": sorted(n for n, r in res.items() if all(r.values())),
            "todo": sorted(n for n, r in res.items() if not all(r.values()))}


# ---------------- the chain's side ----------------
def link(field_dir: str | None, esp_dir: str | None, monomers_xyz: str, library: str) -> list[dict]:
    """
    Replace a dot's monomer runs (field/mono_<label>, esp/<label>) by links to the library entry of
    the same monomer, when that entry has finished, has the same geometry and, after normalized(),
    the same input as the run the chain prepared. Returns one row per run: linked, or why not.
    """
    lib = Path(library)
    rows = []
    for label, name, sym, pts in read_frames(monomers_xyz):
        for kind, top, run in (("hirshfeld", field_dir, f"mono_{label}"), ("esp", esp_dir, label)):
            if not top:
                continue
            mine, entry = Path(top) / run, lib / name / kind
            row = {"monomer": name, "run": f"{Path(top).name}/{run}", "library": str(entry)}
            if mine.is_symlink():
                row["status"] = "linked"
            elif not (mine / "cp2k.inp").is_file():
                row["status"] = "not prepared"
            elif not run_done(entry):
                row["status"] = "not in the library" if not (entry / "cp2k.inp").is_file() else "library run not finished"
            else:
                why = _differs(mine, entry, sym, pts)
                if why:
                    row["status"] = f"differs: {why}"
                else:
                    shutil.rmtree(mine)
                    mine.symlink_to(entry.resolve(), target_is_directory=True)
                    row["status"] = "linked"
            rows.append(row)
    return rows


def _differs(mine: Path, entry: Path, sym, pts) -> str:
    from .xyz import read_xyz
    s2, p2 = read_xyz(entry / "geom.xyz")
    if list(s2) != list(sym) or np.abs(np.asarray(p2, float) - np.asarray(pts, float)).max() > GEOMETRY_TOL_A:
        return "geometry"
    a = normalized((mine / "cp2k.inp").read_text(), sym)
    b = normalized((entry / "cp2k.inp").read_text(), sym)
    if a != b:
        diff = next((f"{x!r} vs {y!r}" for x, y in zip(a.splitlines(), b.splitlines()) if x != y), "length")
        return f"input ({diff})"
    return ""


# ---------------- jobs ----------------
def prepare_jobs(out_dir: str, cluster: dict, cluster_ref: str, materials, *, time_gpu: str = "01:00:00",
                 time_cpu: str = "02:00:00", qos: str | None = None, test_queue: bool = False) -> dict:
    """
    gpu.sh (MACE monomers and their training frames), cpu.sh (the CP2K library runs, side by side on
    one node) and submit.sh (the CPU job after the GPU job) in out_dir. test_queue: both jobs on the
    profile's test partitions (1 h).
    """
    import copy
    if test_queue:
        cluster = copy.deepcopy(cluster)
        for block in (cluster["props"], cluster["cp2k"]):
            if not block.get("test_partition"):
                raise KeyError("--test-queue: the profile names no test_partition")
            block["partition"] = block["test_partition"]
        time_gpu = time_cpu = "01:00:00"
    from .cluster import sbatch_account
    from .cp2k.inputs import cp2k_launch, cp2k_partition, cp2k_setup
    from .props import _job
    out = Path(out_dir).resolve()            # the jobs run in out_dir: every path in them is absolute
    out.mkdir(parents=True, exist_ok=True)
    cp = cluster["cp2k"]
    lib = cp.get("monomer_library")
    if not lib:
        raise KeyError("cluster profile needs cp2k.monomer_library (where the CP2K monomer library lives)")
    sp = species(materials)
    names = " ".join(s["name"] for s in sp)
    pp = cluster["props"]
    run = (f"python -m orchestr_ai.qd.references monomers {names} --ensemble --out {out}/monomers.xyz "
           f"--model {pp['model']} --head {pp.get('head', 'omat_pbe')} > {out}/mace.jsonl")
    (out / "gpu.sh").write_text(_job(cluster, "monomers-mace", time_gpu, qos, run))

    per_node = int(cluster["cores_per_node"])
    n = SIDE_BY_SIDE_RANKS
    py, ppath = cluster["qdw_python"], cluster.get("qdw_pythonpath", "")
    mem = f"{int(cluster.get('node_memory_gb', 0))}G" if cluster.get("node_memory_gb") else "0"
    cpu = f"""#!/bin/bash
#SBATCH -J monomers-cp2k
#SBATCH -t {time_cpu}
#SBATCH -N 1
#SBATCH --ntasks-per-node={per_node}
#SBATCH --mem={mem}
#SBATCH --qos={qos or cluster.get('qos', 'regular')}
{cp2k_partition(cp)}{sbatch_account(cluster)}{chr(10).join(cluster.get('setup', []))}
# CP2K monomer library: every monomer's Hirshfeld single point and ESP cubes, {n} ranks each,
# {per_node // n} side by side on one node (small molecules). Finished runs are skipped.
{cp2k_setup(cp)}
export OMP_NUM_THREADS=1
qdw() {{ env -u PYTHONHOME PYTHONPATH={ppath} LD_LIBRARY_PATH=$(dirname "$(dirname "{py}")")/lib "{py}" -m qd_workflows.cli "$@"; }}
LIB={lib}
qdw monomers write {out}/monomers.xyz --library $LIB --cluster {cluster_ref}
STEP_MEM=$(( ${{SLURM_MEM_PER_NODE:-0}} * 9 / 10 / {per_node // n} ))M; [ "${{SLURM_MEM_PER_NODE:-0}}" -gt 0 ] || STEP_MEM=
cp2k_n() {{ {cp2k_launch(cp, "$1", step_mem=True)} {cp['executable']} -i cp2k.inp -o cp2k.out > srun.log 2>&1; }}
k=0
for d in $(qdw monomers todo --library $LIB); do
    (cd "$d" && cp2k_n {n}) &
    k=$((k + 1)); if [ $k -ge {per_node // n} ]; then wait -n; k=$((k - 1)); fi
done
wait
qdw monomers status --library $LIB
"""
    (out / "cpu.sh").write_text(cpu)
    (out / "submit.sh").write_text(
        "#!/bin/bash\nset -e\ncd \"$(dirname \"$0\")\"\n"
        "g=$(sbatch --parsable gpu.sh)\nc=$(sbatch --parsable --dependency=afterok:${g%%;*} cpu.sh)\n"
        "echo \"monomers: GPU job $g, CPU job $c\"\n")
    (out / "species.json").write_text(json.dumps(sp, indent=1))
    return {"species": sp, "library": lib, "jobs": [str(out / f) for f in ("gpu.sh", "cpu.sh", "submit.sh")]}


def todo(library: str) -> list[str]:
    """The library's run directories that have not finished."""
    st = status(library)
    lib = Path(library)
    return [str(lib / n / k) for n, r in st["monomers"].items() for k, ok in r.items() if not ok]
