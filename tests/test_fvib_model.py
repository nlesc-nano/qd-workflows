"""Learned F_vib (fvib_model) on synthetic zinc-blende fragments: descriptor invariance and additivity, the angle block,
the fit, the domain check, the splits."""
import json

import numpy as np
import pytest

from qd_workflows import fvib_model as fm

A = 6.08


def _zb(n, rng, cl=0):
    """A random compact zinc-blende CdSe fragment of about n atoms, with `cl` Cl on outer anion sites."""
    g = np.array([(i, j, k) for i in range(-4, 5) for j in range(-4, 5) for k in range(-4, 5)])
    se = g[g.sum(axis=1) % 2 == 0] * A / 2                       # fcc anion sites
    cd = se + A / 4
    X = np.r_[se, cd]
    sym = np.array(["Se"] * len(se) + ["Cd"] * len(cd))
    c = rng.normal(size=3)
    keep = np.argsort(np.linalg.norm(X - c, axis=1))[:n]
    X, sym = X[keep], sym[keep]
    out = np.where(sym == "Se")[0][np.argsort(-np.linalg.norm(X[sym == "Se"] - c, axis=1))[:cl]]
    sym[out] = "Cl"
    return list(sym), X + rng.normal(scale=0.15, size=X.shape)


def _rot(rng):
    q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    return q * np.sign(np.linalg.det(q))


def test_descriptor_invariant_and_additive():
    rng = np.random.default_rng(1)
    sym, X = _zb(40, rng, cl=4)
    d = fm.descriptor(sym, X)
    assert d.shape == (len(fm.feature_names()),)
    p = rng.permutation(len(sym))
    d2 = fm.descriptor([sym[i] for i in p], X[p] @ _rot(rng).T + 3.0)
    assert np.allclose(d, d2)
    s2, X2 = _zb(25, rng, cl=2)
    far = fm.descriptor(sym + s2, np.r_[X, X2 + [100.0, 0, 0]])
    assert np.allclose(far, d + fm.descriptor(s2, X2))


def _dataset(rng, n=300):
    temps = [400.0, 500.0, 600.0]
    structs = [_zb(int(rng.integers(20, 90)), rng, cl=int(rng.integers(0, 6))) for _ in range(n)]
    X = fm.descriptors(structs)
    w = rng.normal(scale=1e-3, size=X.shape[1]) * ((X > fm.SEEN).sum(axis=0) > 0)   # features the data have
    Y = np.array([X @ w * t / 500 - 0.11 * X[:, :3].sum(axis=1) * t / 500 for t in temps]).T
    m = X[:, 1].astype(int)
    fw = np.array([f"f{i // 3}" for i in range(n)])
    edge = np.array([f"e{i // 30}" if i % 3 else "" for i in range(n)])
    P = np.array([(i, i + 1) for i in range(n - 1)])
    return {"temps": temps, "X": X, "Y": Y, "m": m, "fw": fw, "edge": edge, "pairs": P,
            "pair_kind": np.array(["growth"] * len(P))}


def test_fit_recovers_linear_target_and_round_trips(tmp_path):
    rng = np.random.default_rng(2)
    data = _dataset(rng)
    model = fm.FvibModel(data["temps"], lam=1e-10, min_rows=1).fit(data["X"], data["Y"])
    assert np.abs(model.predict(data["X"]) - data["Y"]).max() < 1e-6
    f = tmp_path / "model.json"
    f.write_text(json.dumps(model.to_json()))
    back = fm.FvibModel.from_json(json.loads(f.read_text()))
    assert np.allclose(back.predict(data["X"]), model.predict(data["X"]))
    assert back.predict_structure(*_zb(50, rng))[0].shape == (3,)


def test_domain_check():
    rng = np.random.default_rng(3)
    data = _dataset(rng)
    model = fm.FvibModel(data["temps"]).fit(data["X"], data["Y"])
    big = fm.descriptor(*_zb(300, rng, cl=10))
    assert not model.out_of_domain(big)[0]                       # larger, but nothing new
    sym, X = _zb(40, rng)
    X = np.vstack([X, X[sym.index("Se")] + [0.0, 0.0, 2.2]])     # an Se-Cl bond the data never have
    assert model.out_of_domain(fm.descriptor(sym + ["Cl"], X))[0]


def test_splits_and_metrics():
    rng = np.random.default_rng(4)
    data = _dataset(rng)
    sp = fm.splits(data, n_folds=3)
    for c in fm.SIZE_CUTS:
        tr, te = sp[f"size_m<{c}"]
        assert (data["m"][tr] < c).all() and (data["m"][te] >= c).all()
    lab = fm.lineage_labels(data)
    assert all(len(set(lab[data["fw"] == f])) == 1 for f in np.unique(data["fw"]))
    folds = [sp[f"edges_fold{i}"][1] for i in range(3)]
    assert sorted(np.concatenate(folds).tolist()) == list(range(len(data["m"])))
    rep, R = fm.split_metrics(data, *sp["random_80_20"], lam=1e-8, report_t=(500.0,))
    assert R.shape == data["Y"].shape
    assert rep["T"]["500"]["step_growth"]["mae"] < 1e-3


def test_target_matches_thermo():
    pytest.importorskip("orchestr_ai.qd.steps.hessian")
    from orchestr_ai.qd.steps.hessian import thermo
    f = [5.0, 80.0, 150.0, 300.0, -20.0]
    assert np.allclose(fm.fvib_target(f, [500.0]), thermo(np.array(f), [500.0])["F_vib_eV"])


def test_angular_block_of_a_tetrahedron():
    t = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]]) * 2.63 / np.sqrt(3)
    d = dict(zip(fm.feature_names(), fm.descriptor(["Cd"] + ["Se"] * 4, np.r_[[[0.0, 0, 0]], t])))
    assert d["coord_Cd_4"] == 1 and d["coord_Se_1"] == 4
    ang = np.array([d[f"ang_Cd_{c:+.2f}"] for c in fm.COS_C])
    assert np.allclose(ang, 6 * np.exp(-((fm.COS_C + 1 / 3) / fm.COS_W) ** 2))     # six angles at cos = -1/3


def test_species_parameter_for_acetate_dots():
    rng = np.random.default_rng(5)
    sym, X = _zb(80, rng)
    ac = [("O", [1.2, 0, 0]), ("O", [-1.2, 0, 0]), ("C", [0, 0.7, 0]), ("C", [0, 2.2, 0]), ("H", [0.9, 2.6, 0.4])]
    sy2 = sym + [s for s, _ in ac]
    X2 = np.r_[X, np.array([v for _, v in ac]) + [20.0, 0, 0]]
    with pytest.raises(ValueError):
        fm.descriptor(sy2, X2)                                   # O, C, H are not in the CdCl2 species
    d = fm.descriptor(sy2, X2, fm.SPECIES_ACETATE)
    names = fm.feature_names(fm.SPECIES_ACETATE)
    assert d.shape == (len(names),)
    assert dict(zip(names, d))["N_H"] == 1
    data = {"temps": [300.0, 500.0], "X": np.array([fm.descriptor(*_zb(60 + i, rng), fm.SPECIES_ACETATE)
                                                    for i in range(30)])}
    data["Y"] = np.c_[data["X"][:, 0], 2 * data["X"][:, 1]] * 0.01
    model = fm.FvibModel(data["temps"], species=fm.SPECIES_ACETATE).fit(data["X"], data["Y"])
    back = fm.FvibModel.from_json(json.loads(json.dumps(model.to_json())))
    assert back.species == fm.SPECIES_ACETATE
    assert np.allclose(back.predict(data["X"]), model.predict(data["X"]))
    assert back.out_of_domain(d)[0]                              # an O-O pair at 2.4 A: never seen in training
