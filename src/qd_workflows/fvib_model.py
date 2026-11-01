"""
Learned vibrational free energy of relaxed builder dots (CdSe/CdCl2), for states of the network growth model that are
relaxed with MACE but get no finite-difference Hessian.

The target is the quasi-RRHO F_vib(T) of orchestr_ai.qd.steps.hessian.thermo (harmonic U, Grimme entropy, modes
floored at 10 cm-1), the vibrational part of the network's G (references.ideal_gas_g). The model is per-atom
additive on the relaxed geometry and linear in its descriptor, one weight set per temperature:

    F_vib(T) = sum_{i<j} sum_n w_{s_i s_j, n}(T) phi_n(r_ij) f_c(r_ij) + sum_i sum_p c_{s_i, p}(T) g_p(rho_i)
               + sum_s b_s(T) N_s + b_0(T)

phi_n: Gaussians of width 0.25 A on 2-7 A every 0.25 A; f_c: cosine cutoff at 7.5 A; rho_i: the atom's smooth
coordination (sum_j exp(-((r_ij - 2.7)/0.3)^2)), g_p in (rho, rho^2, sqrt(rho)). Ridge on standardised features.

Inputs (formats of the tree scripts, as network_export): results_by_key.json (key -> symbols, positions,
frequencies_cm1, n_imaginary), registry.json (dots: key -> id, framework "se"), layer_steps.json (a, b, dSe, dCd,
edge), and the grid's thermo_merged.json (temperatures).

    python -m qd_workflows.fvib_model fit --results results_by_key.json --registry registry.json \
        --steps layer_steps.json --thermo thermo_merged.json --out DIR [--nproc 4] [--lam L]
    python -m qd_workflows.fvib_model predict --model DIR/model.json file.xyz ...

fit writes model.json (trained on every dot), validation.json (size and lineage-edge splits, random as a reference),
residuals.png, and caches the descriptors in dataset.npz (a second fit with the same inputs reuses it).
"""
from __future__ import annotations

import argparse
import collections
import json
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np

SPECIES = ("Cd", "Se", "Cl")
PAIRS = tuple((a, b) for i, a in enumerate(SPECIES) for b in SPECIES[i:])
CENTRES = np.arange(2.0, 7.001, 0.25)
WIDTH = 0.25
R_CUT = 7.5
RHO_R0, RHO_W = 2.7, 0.3
EMBED = ("rho", "rho2", "sqrt_rho")
SIZE_CUTS = (30, 60, 100)
REPORT_T = (400.0, 500.0, 600.0)
SEEN = 0.05                   # a feature value above this counts as present (a twentieth of a pair)
RARE = 0.05                   # present in fewer than this share of the training dots: a rare feature
BOND_R = 3.5                  # pair features up to this distance are the bonding shell (the domain check)


def feature_names() -> list[str]:
    return ([f"N_{s}" for s in SPECIES] + [f"pair_{a}{b}_{r:.2f}" for a, b in PAIRS for r in CENTRES]
            + [f"emb_{s}_{g}" for s in SPECIES for g in EMBED])


BOND = np.array([n.startswith("pair_") and float(n.rsplit("_", 1)[1]) <= BOND_R for n in feature_names()])


def descriptor(symbols, positions) -> np.ndarray:
    """Composition, pair radial sums per species pair, and per-species sums of the embedding functions."""
    from scipy.spatial import cKDTree
    sy = np.asarray(symbols)
    X = np.asarray(positions, float)
    p = cKDTree(X).query_pairs(R_CUT, output_type="ndarray").reshape(-1, 2)
    r = np.linalg.norm(X[p[:, 0]] - X[p[:, 1]], axis=1)
    fc = 0.5 * (1 + np.cos(np.pi * r / R_CUT))
    g = np.exp(-((r[:, None] - CENTRES) / WIDTH) ** 2) * fc[:, None]
    si, sj = sy[p[:, 0]], sy[p[:, 1]]
    pair = [g[((si == a) & (sj == b)) | ((si == b) & (sj == a))].sum(axis=0) for a, b in PAIRS]
    rho = np.zeros(len(sy))
    w = np.exp(-((r - RHO_R0) / RHO_W) ** 2)
    np.add.at(rho, p[:, 0], w)
    np.add.at(rho, p[:, 1], w)
    emb = []
    for s in SPECIES:
        x = rho[sy == s]
        emb += [x.sum(), (x * x).sum(), np.sqrt(x).sum()]
    comp = [float((sy == s).sum()) for s in SPECIES]
    return np.concatenate([comp, np.concatenate(pair), emb])


def _desc_job(job):
    return descriptor(*job)


def descriptors(structs, nproc: int = 1) -> np.ndarray:
    if nproc > 1:
        with Pool(nproc) as pool:
            return np.array(pool.map(_desc_job, structs, chunksize=64))
    return np.array([descriptor(s, x) for s, x in structs])


def fvib_target(freqs_cm1, temps) -> np.ndarray:
    """F_vib (eV) on the temperatures, exactly as the network's G has it."""
    from orchestr_ai.qd.steps.hessian import thermo
    return np.array(thermo(np.asarray(freqs_cm1, float), list(temps))["F_vib_eV"])


class FvibModel:
    """Ridge on standardised features, one weight column per temperature; the intercept is not penalised.
    lam is per training row: (Z'Z + lam n I) w = Z'y. A feature seen in fewer than `min_rows` training dots (a bond
    type the data hardly have, such as Se-Cl at 2.2 A) gets no weight. The domain check looks at the bonding shell
    (pairs up to BOND_R): a structure is out of domain if it has a bond type without weight, or more of a rare one
    (present in < 5 % of the training dots, such as Se-Se at 2.7 A) than any training dot has. A larger dot is not
    out of domain for being larger."""

    def __init__(self, temps, lam: float = 1e-6, min_rows: int = 20):
        self.temps = [float(t) for t in temps]
        self.lam = float(lam)
        self.min_rows = int(min_rows)
        self.mu = self.sd = self.w = self.b = self.active = self.rare = self.xmax = None

    def fit(self, X, Y) -> "FvibModel":
        X, Y = np.asarray(X, float), np.asarray(Y, float)
        seen = (X > SEEN).sum(axis=0)
        self.active = seen >= self.min_rows
        self.rare = self.active & (seen < RARE * len(X))
        self.xmax = X.max(axis=0)
        self.mu = np.where(self.active, X.mean(axis=0), 0.0)
        self.sd = np.where(self.active, X.std(axis=0), 1.0)
        self.sd[self.sd < 1e-12] = 1.0
        Z = ((X - self.mu) / self.sd)[:, self.active]
        ym = Y.mean(axis=0)
        A = Z.T @ Z + self.lam * len(Z) * np.eye(Z.shape[1])
        self.w = np.zeros((X.shape[1], Y.shape[1]))
        self.w[self.active] = np.linalg.solve(A, Z.T @ (Y - ym))
        self.b = ym
        return self

    def predict(self, X) -> np.ndarray:
        return ((np.atleast_2d(X) - self.mu) / self.sd) @ self.w + self.b

    def out_of_domain(self, X) -> np.ndarray:
        """True for structures with a feature the model was not trained on."""
        X = np.atleast_2d(X)
        new = BOND & ~self.active
        rare = BOND & self.rare
        return (X[:, new] > SEEN).any(axis=1) | (X[:, rare] > 1.1 * self.xmax[rare] + SEEN).any(axis=1)

    def predict_structure(self, symbols, positions) -> tuple[np.ndarray, bool]:
        x = descriptor(symbols, positions)
        return self.predict(x)[0], bool(self.out_of_domain(x)[0])

    def at(self, T: float) -> int:
        return int(np.argmin(np.abs(np.array(self.temps) - T)))

    def to_json(self) -> dict:
        return {"kind": "fvib_pair_embedding_ridge", "features": feature_names(), "temps": self.temps,
                "lam": self.lam, "min_rows": self.min_rows, "active": self.active.tolist(), "rare": self.rare.tolist(),
                "xmax": self.xmax.tolist(), "mu": self.mu.tolist(),
                "sd": self.sd.tolist(), "w": self.w.tolist(), "b": self.b.tolist(),
                "descriptor": {"centres": CENTRES.tolist(), "width": WIDTH, "r_cut": R_CUT, "rho_r0": RHO_R0,
                               "rho_w": RHO_W}}

    @classmethod
    def from_json(cls, d: dict) -> "FvibModel":
        if d["features"] != feature_names():
            raise ValueError("model.json was written with a different descriptor")
        m = cls(d["temps"], d["lam"], d["min_rows"])
        m.mu, m.sd, m.w, m.b = (np.array(d[k]) for k in ("mu", "sd", "w", "b"))
        m.active, m.rare = np.array(d["active"], bool), np.array(d["rare"], bool)
        m.xmax = np.array(d["xmax"])
        return m


# ---------------------------------------------------------------- data set and validation

def load_dataset(results, registry, steps, temps, nproc: int = 1, cache: Path | None = None) -> dict:
    """Descriptors, targets and labels of every relaxed dot, the step pairs between them and their lineage edges."""
    res = json.loads(Path(results).read_text())
    keys = list(res)
    if cache is not None and cache.exists():
        z = np.load(cache, allow_pickle=False)
        if list(z["keys"]) == keys and np.allclose(z["temps"], temps):
            data = {k: z[k] for k in z.files}
            data["temps"] = list(temps)
            print(f"[fvib] descriptors from {cache}", flush=True)
            return data
    reg = json.loads(Path(registry).read_text())["dots"]
    row_of_id = {reg[k]["id"]: j for j, k in enumerate(keys) if k in reg}
    fw = np.array([str(reg.get(k, {}).get("se", "")) for k in keys])
    cnt = [collections.Counter(res[k]["symbols"]) for k in keys]
    m = np.array([c["Se"] for c in cnt])
    k_ = np.array([c["Cd"] - c["Se"] for c in cnt])
    print(f"[fvib] {len(keys)} relaxed dots, descriptors on {nproc} processes", flush=True)
    X = descriptors([(res[k]["symbols"], res[k]["positions"]) for k in keys], nproc)
    Y = np.array([fvib_target(res[k]["frequencies_cm1"], temps) for k in keys])
    E = np.array([res[k]["E_eV"] for k in keys])
    nimag = np.array([int(res[k].get("n_imaginary", 0)) for k in keys])
    pairs, kinds, edge_votes = set(), {}, collections.defaultdict(collections.Counter)
    for r in json.loads(Path(steps).read_text()):
        a, b = row_of_id.get(r["a"]), row_of_id.get(r["b"])
        if a is None or b is None or a == b:
            continue
        a, b = (a, b) if (m[a], k_[a]) <= (m[b], k_[b]) else (b, a)
        pairs.add((a, b))
        kinds[(a, b)] = "growth" if m[b] > m[a] else "other"
        edge_votes[a][r["edge"]] += 1
        edge_votes[b][r["edge"]] += 1
    # shell links as the export builds them: lowest-E isomer of each (framework, k) to that of the next k
    low = {}
    for j in range(len(keys)):
        if fw[j] and ((fw[j], k_[j]) not in low or E[j] < E[low[(fw[j], k_[j])]]):
            low[(fw[j], k_[j])] = j
    for (f, kk), j in low.items():
        i2 = low.get((f, kk + 1))
        if i2 is not None:
            pairs.add((j, i2))
            kinds[(j, i2)] = "shell"
    P = np.array(sorted(pairs), dtype=np.int64).reshape(-1, 2)
    edge = np.array([edge_votes[j].most_common(1)[0][0] if edge_votes[j] else "" for j in range(len(keys))])
    data = {"keys": np.array(keys), "temps": list(temps), "X": X, "Y": Y, "m": m, "k": k_, "fw": fw, "edge": edge,
            "n_imag": nimag, "pairs": P, "pair_kind": np.array([kinds[tuple(p)] for p in P])}
    if cache is not None:
        np.savez_compressed(cache, **{k: v for k, v in data.items() if k != "temps"}, temps=np.array(temps))
    return data


def split_metrics(data: dict, train, test, lam: float, report_t=REPORT_T) -> tuple[dict, np.ndarray]:
    """Train on `train`, score on the in-domain part of `test`: absolute, within-size and step errors (both ends of a
    step in it). Out-of-domain test dots (a feature the training set lacks) are counted, not scored."""
    model = FvibModel(data["temps"], lam).fit(data["X"][train], data["Y"][train])
    R = data["Y"] - model.predict(data["X"])                       # residual, (dots, temps)
    m = data["m"]
    test = np.asarray(test)
    ood = model.out_of_domain(data["X"][test])
    test = test[~ood]
    in_test = np.zeros(len(m), bool)
    in_test[test] = True
    P, kind = data["pairs"], data["pair_kind"]
    both = in_test[P[:, 0]] & in_test[P[:, 1]]
    out = {"n_train": int(len(train)), "n_test": int(len(test)), "n_out_of_domain": int(ood.sum()), "T": {}}
    for T in report_t:
        it = model.at(T)
        r = R[test, it]
        mt = m[test]
        centred = r - (np.bincount(mt, r) / np.maximum(np.bincount(mt), 1))[mt]
        row = {"abs_mae": float(np.abs(r).mean()), "abs_rmse": float(np.sqrt((r ** 2).mean())),
               "within_size_mae": float(np.abs(centred).mean()),
               "drift_eV_per_100m": float(np.polyfit(m[test], r, 1)[0] * 100) if np.ptp(m[test]) else 0.0}
        for kk in ("growth", "shell", "all"):
            sel = both if kk == "all" else both & (kind == kk)
            d = np.abs(R[P[sel, 1], it] - R[P[sel, 0], it])
            row[f"step_{kk}"] = ({"n": int(sel.sum()), "mae": float(d.mean()), "p95": float(np.percentile(d, 95)),
                                  "frac_gt_0.05": float((d > 0.05).mean()), "frac_gt_0.1": float((d > 0.1).mean())}
                                 if sel.any() else {"n": 0})
        out["T"][f"{T:g}"] = row
    return out, R


def lineage_labels(data: dict) -> np.ndarray:
    """Each dot's lineage edge; a dot on no layer step (shell variants, fill dots) takes its framework's most common
    edge, else its framework, so that a framework never sits on both sides of a split."""
    by_fw = collections.defaultdict(collections.Counter)
    for e, f in zip(data["edge"], data["fw"]):
        if e:
            by_fw[f][e] += 1
    return np.array([e or (by_fw[f].most_common(1)[0][0] if by_fw[f] else f"fw:{f}")
                     for e, f in zip(data["edge"], data["fw"])])


def splits(data: dict, n_folds: int = 5, seed: int = 0) -> dict:
    """Size cuts (train m < cut, test m >= cut), lineage edges in folds, and a random 80/20 split for reference."""
    m, idx = data["m"], np.arange(len(data["m"]))
    rng = np.random.default_rng(seed)
    out = {f"size_m<{c}": (idx[m < c], idx[m >= c]) for c in SIZE_CUTS}
    label = lineage_labels(data)
    edges = np.unique(label)
    fold_of = dict(zip(rng.permutation(edges), np.arange(len(edges)) % n_folds))
    f = np.array([fold_of[e] for e in label])
    for i in range(n_folds):
        out[f"edges_fold{i}"] = (idx[f != i], idx[f == i])
    r = rng.random(len(m))
    out["random_80_20"] = (idx[r < 0.8], idx[r >= 0.8])
    return out


def choose_lam(data: dict, grid=(1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2), T: float = 500.0) -> tuple[float, dict]:
    """The lam with the lowest growth-step MAE over the lineage-edge folds at T."""
    folds = {k: v for k, v in splits(data).items() if k.startswith("edges_")}
    scan = {}
    for lam in grid:
        maes = [split_metrics(data, tr, te, lam, (T,))[0]["T"][f"{T:g}"]["step_growth"].get("mae", np.nan)
                for tr, te in folds.values()]
        scan[f"{lam:g}"] = float(np.nanmean(maes))
    best = min(scan, key=scan.get)
    return float(best), scan


def validate(data: dict, lam: float) -> tuple[dict, dict]:
    report, resid = {}, {}
    for name, (tr, te) in splits(data).items():
        report[name], R = split_metrics(data, tr, te, lam)
        resid[name] = (te, R)
    report["edges_mean"] = _mean_of([v for k, v in report.items() if k.startswith("edges_")])
    return report, resid


def _mean_of(rows: list):
    """Element-wise mean of nested metric dicts (the lineage-edge folds)."""
    if isinstance(rows[0], dict):
        return {k: _mean_of([r.get(k, np.nan) for r in rows]) for k in max(rows, key=len)}
    return float(np.nanmean(rows))


def plot(data: dict, resid: dict, path: Path, T: float = 500.0) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    it = int(np.argmin(np.abs(np.array(data["temps"]) - T)))
    names = [f"size_m<{c}" for c in SIZE_CUTS]
    fig, ax = plt.subplots(2, len(names), figsize=(4.2 * len(names), 7), squeeze=False)
    P, kind = data["pairs"], data["pair_kind"]
    for j, name in enumerate(names):
        te, R = resid[name]
        ax[0, j].scatter(data["m"][te], R[te, it], s=3, alpha=0.3)
        ax[0, j].axhline(0, color="k", lw=0.6)
        ax[0, j].set(title=f"train {name.split('_')[1]}, test the rest", xlabel="m (CdSe)",
                     ylabel=f"F_vib residual at {T:g} K (eV)")
        in_test = np.zeros(len(data["m"]), bool)
        in_test[te] = True
        sel = in_test[P[:, 0]] & in_test[P[:, 1]] & (kind == "growth")
        d = R[P[sel, 1], it] - R[P[sel, 0], it]
        ax[1, j].hist(d, bins=np.linspace(-0.15, 0.15, 61))
        ax[1, j].set(xlabel="step error (eV)", ylabel="growth steps",
                     title=f"MAE {np.abs(d).mean():.3f} eV, >0.1 eV {np.mean(np.abs(d) > 0.1):.1%}")
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


def read_xyz(path):
    lines = Path(path).read_text().split("\n")
    at = [ln.split() for ln in lines[2:2 + int(lines[0])]]
    return [a[0] for a in at], np.array([[float(v) for v in a[1:4]] for a in at])


def fit(a) -> dict:
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    temps = json.loads(Path(a.thermo).read_text())["temperatures"]
    data = load_dataset(a.results, a.registry, a.steps, temps, a.nproc, out / "dataset.npz")
    print(f"[fvib] {len(data['m'])} dots, {len(data['pairs'])} step pairs "
          f"{dict(collections.Counter(data['pair_kind'].tolist()))}, {len(np.unique(lineage_labels(data)))} lineage "
          f"groups", flush=True)
    scan = None
    lam = a.lam
    if lam is None:
        lam, scan = choose_lam(data)
        print(f"[fvib] lam {lam:g} (edge-fold growth-step MAE at 500 K: {scan})", flush=True)
    report, resid = validate(data, lam)
    model = FvibModel(temps, lam).fit(data["X"], data["Y"])
    (out / "model.json").write_text(json.dumps(model.to_json()))
    (out / "validation.json").write_text(json.dumps({"lam": lam, "lam_scan": scan, "splits": report}, indent=1))
    plot(data, resid, out / "residuals.png")
    for name in [f"size_m<{c}" for c in SIZE_CUTS] + ["edges_mean", "random_80_20"]:
        r = report[name]["T"]["500"]
        print(f"[fvib] {name:14s} 500 K: growth step MAE {r['step_growth']['mae']:.4f} eV "
              f"(>0.1 eV {r['step_growth']['frac_gt_0.1']:.3f}), shell {r['step_shell']['mae']:.4f}, "
              f"within size {r['within_size_mae']:.4f}, abs {r['abs_mae']:.4f}, drift {r['drift_eV_per_100m']:+.3f} "
              f"eV/100 m; {report[name]['n_out_of_domain']:.0f} test dots out of domain", flush=True)
    return report


def predict(a) -> None:
    model = FvibModel.from_json(json.loads(Path(a.model).read_text()))
    it = model.at(a.T)
    for f in a.xyz:
        fv, ood = model.predict_structure(*read_xyz(f))
        print(f"{f}\tF_vib({model.temps[it]:g} K) = {fv[it]:.4f} eV" + ("\tOUT OF DOMAIN" if ood else ""))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python -m qd_workflows.fvib_model")
    sub = ap.add_subparsers(dest="cmd", required=True)
    f = sub.add_parser("fit", help="fit on the relaxed dots, validate, write model.json")
    f.add_argument("--results", required=True)
    f.add_argument("--registry", required=True)
    f.add_argument("--steps", required=True)
    f.add_argument("--thermo", required=True, help="the grid's thermo_merged.json (temperatures)")
    f.add_argument("--out", required=True)
    f.add_argument("--nproc", type=int, default=4)
    f.add_argument("--lam", type=float, default=None, help="ridge strength per row; default: chosen on edge folds")
    p = sub.add_parser("predict", help="F_vib of relaxed structures (xyz)")
    p.add_argument("--model", required=True)
    p.add_argument("--T", type=float, default=500.0)
    p.add_argument("xyz", nargs="+")
    a = ap.parse_args(argv)
    fit(a) if a.cmd == "fit" else predict(a)
    return 0


if __name__ == "__main__":
    sys.exit(main())
