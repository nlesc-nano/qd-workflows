"""MA detachment grid: candidate pairs, fragments, job script."""
import numpy as np
import pytest

pytest.importorskip("orchestr_ai.qd.steps.structure")
from qd_workflows import grid  # noqa: E402

CHARGES = {"Cd": 2, "Se": -2, "Cl": -1}
NATIVE = ["Cd", "Se"]
BOND = 2.63


def _cd2se_cl():
    """Se bonded to two Cd; the second Cd carries one Cl."""
    sym = ["Se", "Cd", "Cd", "Cl"]
    pts = np.array([[0, 0, 0], [BOND, 0, 0], [0, BOND, 0], [0, BOND + 2.45, 0]], float)
    return sym, pts


def test_pairs_are_bonded_cation_anion_pairs():
    sym, pts = _cd2se_cl()
    units = grid.enumerate_pairs(sym, pts, CHARGES, NATIVE, BOND)
    assert sorted((u["cation"], tuple(u["ligands"])) for u in units) == [(1, (0,)), (2, (0,))]
    by_cd = {u["cation"]: u for u in units}
    assert by_cd[2]["cation_ligands"] == 1 and by_cd[1]["cation_ligands"] == 0
    assert all(u["atoms"] == frozenset((u["cation"], 0)) and u["ligand_mu"] == [2] for u in units)


def test_fragments_and_cations_without_anion():
    sym, pts = _cd2se_cl()
    assert grid.fragments(sym, pts, CHARGES, NATIVE, BOND) == 1
    far = pts.copy()
    far[3] += [0, 20, 0]                       # the Cl leaves
    assert grid.fragments(sym, far, CHARGES, NATIVE, BOND) == 2
    sym2 = sym + ["Cd"]
    pts2 = np.vstack([pts, pts[3] + [0, 2.45, 0]])  # a Cd held by the Cl only
    assert grid.cations_without_anion(sym2, pts2, CHARGES, NATIVE, BOND) == 1


def test_prepare_writes_a_test_partition_job(tmp_path):
    cluster = {"qos": "normal", "account": "acc", "setup": [],
               "props": {"partition": "gpu:shared", "test_partition": "gpu:test", "gres": "gpu:A100:1",
                         "conda_env": "env", "model": "/m/model", "gxtb": "", "refs": "/r", "head": "omat_pbe"}}
    info = grid.prepare("/rec/Cd16", str(tmp_path), cluster, "/code/src")
    job = (tmp_path / "job.sh").read_text()
    assert info["partition"] == "gpu:test" and "#SBATCH -p gpu:test" in job
    assert "env PYTHONPATH=/code/src:$PYTHONPATH python -u -m qd_workflows.grid run /rec/Cd16" in job
    assert "--model /m/model --head omat_pbe" in job


def test_prepare_pieces_requeue_and_stop_before_the_limit(tmp_path):
    cluster = {"qos": "normal", "account": "acc", "setup": [],
               "props": {"partition": "gpu:shared", "test_partition": "gpu:test", "gres": "gpu:A100:1",
                         "conda_env": "env", "model": "/m/model", "gxtb": "", "refs": "/r"}}
    info = grid.prepare("/rec/Cd16", str(tmp_path), cluster, "/code/src", partition="gpu:test,gpu:shared", pieces=True)
    job = (tmp_path / "job.sh").read_text()
    assert info["partition"] == "gpu:test,gpu:shared" and "#SBATCH -p gpu:test,gpu:shared" in job
    assert "--dependency=afterany:$SLURM_JOB_ID job.sh" in job
    assert "timeout 3540 env PYTHONPATH=/code/src" in job
