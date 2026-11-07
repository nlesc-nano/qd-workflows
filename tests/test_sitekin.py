"""qdw sites (site-kinetics protocol): picking and the role table against the 2026-10 LUMI results of the 1.0 dot,
the GPU steps with a toy EMT engine in place of MACE, and the job script."""
import json
import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from qd_workflows.sitekin import engine as en
from qd_workflows.sitekin import jobs, neb, recipe, sites, states, table

DATA = Path(__file__).parent / "data" / "sitekin"


def load(name):
    return json.loads((DATA / name).read_text())


def test_recipe_template():
    r = recipe.load()
    assert [recipe.tag(s["size"]) for s in r["sizes"]] == ["1p0", "1p5", "2p0", "2p5", "3p0"]
    assert recipe.size_entry(r, "1p0")["surface"] == "plain" and recipe.size_entry(r, 2.5)["surface"] == "half_swap"
    assert recipe.kt(r) == pytest.approx(0.0430867, rel=1e-5)


def test_pick_reproduces_the_computed_1p0_job():
    r = recipe.load()
    job, summary = sites.pick(load("dot_1p0.json"), load("reference_1p0.json"), 4.5, 6.0, 7.0)
    old = load("job_1p0_2026-10.json")
    assert sorted(job["sites"]) == sorted(old["sites"])
    for k, v in job["sites"].items():
        o = old["sites"][k]
        for f in ("cd", "cl", "se_at", "cl_removed", "local", "n_context"):
            assert v[f] == o[f], f
        for st in ("ref", "empty", "M"):
            assert v["states"][st] == o["states"][st]
    row = next(x for x in summary if x["role"] == "edge:100|111A")
    assert row["n"] == 12 and [e["n"] for e in row["environments"]] == [12]
    assert all(e["no_whole_unit"] for x in summary if x is not row for e in x["environments"])
    # the default 5.5 A reach borrows no second Cl on this dot: no whole unit at all
    job2, summary2 = sites.pick(load("dot_1p0.json"), load("reference_1p0.json"), 4.5, 6.0, r["sites"]["borrow_A"])
    assert not job2["sites"] and all(e["no_whole_unit"] for x in summary2 for e in x["environments"])


def test_table_reproduces_the_computed_1p0_row():
    r = recipe.load()
    _job, summary = sites.pick(load("dot_1p0.json"), load("reference_1p0.json"), 4.5, 6.0, 7.0)
    t = table.size_table(summary, load("states_1p0_2026-10.json"), r["material"]["g_bulk_eV"], recipe.kt(r))
    (old,) = load("role_table_1p0_2026-10.json")["edge:100|111A"]["envs"]
    (env,) = t["edge:100|111A"]["envs"]
    assert env["n"] == old["n"] == 12
    assert env["bind"] == pytest.approx(old["bind"], abs=1e-9)
    assert env["dmu"] == pytest.approx(old["dmu"], abs=1e-9)
    assert env["K_X"] == pytest.approx(math.exp(-old["bind"] / recipe.kt(r)))
    assert t["edge:100|111A"]["K_X_mean"] == pytest.approx(env["K_X"])


def test_free_energy_pieces_match_the_computed_cdcl2():
    st = load("states_1p0_2026-10.json")["CdCl2"]
    assert en.vib_free(st["freqs"], 500.0, 5) == pytest.approx(st["F_vib"], abs=1e-12)
    assert en.gas_trans_rot(["Cd", "Cl", "Cl"], st["positions"], 500.0) == pytest.approx(st["G_trans_rot_1M"], abs=1e-6)


class ToyEngine:
    """EMT in place of MACE (both precisions), optionally in a periodic cell; relax_many by ASE L-BFGS."""

    def __init__(self, cell=None, pbc=False):
        from ase.calculators.emt import EMT
        self.calc, self.cell, self.pbc = EMT(), cell, pbc

    def _atoms(self, symbols, x):
        from ase import Atoms
        a = Atoms(symbols, positions=np.asarray(x, float), cell=self.cell, pbc=self.pbc)
        a.calc = self.calc
        return a

    def _one(self, symbols, x):
        a = self._atoms(symbols, x)
        return a.get_potential_energy(), a.get_forces()

    def ev32(self, items):
        return [self._one(s, x) for s, x in items]

    ev64 = ev32

    def relax_many(self, items, fmax, max_steps, precision="float32"):
        from ase.optimize import LBFGS
        out = []
        for s, x in items:
            a = self._atoms(s, x)
            opt = LBFGS(a, logfile=None)
            ok = opt.run(fmax=fmax, steps=max_steps)
            out.append(SimpleNamespace(positions=a.positions.copy(), energy=a.get_potential_energy(), steps=opt.nsteps,
                                       reason="fmax" if ok else "steps"))
        return out


def _toy_recipe():
    r = recipe.load()
    r["material"].update({"cation": "Au", "anion": "Ag"})
    r["ligand"]["anion"] = "Cu"
    return r


def test_states_step_runs_and_keeps_the_bookkeeping(tmp_path):
    from ase.build import fcc100
    slab = fcc100("Al", size=(3, 3, 2), vacuum=0.0)
    x = slab.positions - slab.positions.mean(0)
    top = int(np.argmax(x[:, 2]))
    ctx = [i for i in range(len(slab)) if i != top]
    sc, xc = ["Al"] * len(ctx), x[ctx]
    au = x[top]
    cu = [au + [1.8, 0, 1.2], au + [-1.8, 0, 1.2]]
    job = {"dot": "toy", "reference": "relaxed", "sites": {"s": {
        "role": "facet:100", "env": "toy", "cd": 0, "n_context": len(ctx), "local": list(range(len(ctx)))[-4:],
        "states": {"ref": {"symbols": sc + ["Au", "Cu", "Cu"], "positions": np.vstack([xc, au, *cu]).tolist()},
                   "empty": {"symbols": sc, "positions": xc.tolist()},
                   "M": {"symbols": sc + ["Au", "Ag"], "positions": np.vstack([xc, au, cu[0]]).tolist()}}}}}
    S = states.run_job(ToyEngine(), job, tmp_path / "states.json", _toy_recipe(), log=lambda m: None)
    s = S["sites"]["s"]
    assert set(s["vib"]) == {"ref", "empty", "M"} and S["done"]
    assert s["dG_bind_X"] == pytest.approx(s["G"]["ref"] - s["G"]["empty"] - S["CdCl2"]["G"])
    # resumable: a second call changes nothing
    S2 = states.run_job(ToyEngine(), job, tmp_path / "states.json", _toy_recipe(), log=lambda m: None)
    assert S2["sites"]["s"]["G"] == s["G"]


def test_neb_runs_a_toy_hop_with_one_negative_mode(tmp_path):
    """an Au adatom hopping hollow -> hollow on periodic Al(100); the second incoming atom (Al) sits still"""
    from ase.build import add_adsorbate, fcc100
    slab = fcc100("Al", size=(3, 3, 3), vacuum=6.0)
    add_adsorbate(slab, "Au", 1.7, "hollow")
    add_adsorbate(slab, "Al", 1.7, "hollow", offset=(1, 1))
    x = slab.positions.copy()
    n0 = len(slab) - 2
    fs = x.copy()
    fs[n0, 0] += slab.cell[0, 0] / 3
    bottom = [i for i in range(n0) if x[i, 2] < x[:n0, 2].min() + 0.5]
    meta = {"kind": "dir", "n_context": n0, "mobile": [i for i in range(n0) if i not in bottom],
            "local_phva": [i for i in range(n0) if x[i, 2] > x[:n0, 2].max() - 0.5], "unit": [], "incoming": [n0, n0 + 1],
            "routes": ["A"], "symbols": slab.get_chemical_symbols(), "normal": [0, 0, 1], "start_height_A": 0.0}
    ep = {"structures": {"IS": x.tolist(), "FS_A": fs.tolist()}, "meta": meta}
    R = neb.run_neb(ToyEngine(slab.cell, [True, True, False]), ep, 0.0, _toy_recipe(), tmp_path, log=lambda m: None)
    res = R["result"]
    assert not res["no_barrier"] and res["dE_act"] > 0.1
    assert res["TS_lowest_cm"][0] < 0 < res["TS_lowest_cm"][1] and res["TS_ok"]      # soft (heavy Au): above -30
    assert res["TS_overlap"][0] > 0.5
    assert res["k1_s"] == pytest.approx(1.380649e-23 * 500 / 6.62607015e-34 * math.exp(-res["dG_act"] / (en.KB * 500)))
    assert (tmp_path / "TS.xyz").exists() and R["ci"]["ci_used"]


def test_marked_symbols_and_site_choice():
    meta = {"symbols": ["Se", "Cd", "Cl", "Cl", "Cd", "Se"], "unit": [1, 2, 3], "incoming": [4, 5]}
    assert neb.marked(meta) == ["Se", "Hg", "Br", "Br", "Zn", "Te"]
    summary = [{"role": "edge", "environments": [{"n": 3, "no_whole_unit": False, "site_name": "a"},
                                                 {"n": 9, "no_whole_unit": False, "site_name": "b"},
                                                 {"n": 20, "no_whole_unit": True, "site_name": None}]}]
    assert neb.choose_sites(recipe.load(), summary) == {"edge": "b"}


def test_run_recipe_is_copied_and_kept(tmp_path):
    from qd_workflows.sitekin import build
    run = tmp_path / "run"
    r, p = recipe.for_run(run)
    assert p == recipe.template_path()
    build.run(r, run, sizes=["9p9"], recipe_path=p, log=lambda m: None)       # no such size: only the recipe copy
    r2, p2 = recipe.for_run(run)
    assert p2 == run / "recipe.yaml" and r2 == r
    other = tmp_path / "other.yaml"
    other.write_text(recipe.template_path().read_text().replace("T_K: 500.0", "T_K: 450.0"))
    with pytest.raises(ValueError):
        build.run(recipe.load(other), run, sizes=["9p9"], recipe_path=other, log=lambda m: None)


def test_job_script():
    r = recipe.load()
    s = jobs.script(r, "/scratch/run", ["reference", "states"], None, time="01:00:00")
    assert "#SBATCH -A project_465003437" in s and "#SBATCH -t 01:00:00" in s
    assert "sites reference /scratch/run > log_reference.txt 2>&1 && " in s and "--recipe" not in s
    assert '"\\$WITH_CONDA; ' in s and "QDW_OEQ=1" in s
    with pytest.raises(ValueError):
        jobs.script(r, "/x", ["table"])
