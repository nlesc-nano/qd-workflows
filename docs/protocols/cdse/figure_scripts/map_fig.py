"""Omega(m, k) of the exported acetate network at the pools of the default run: the landscape the growth model walks."""
import json, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from qd_workflows.network_growth import load, _pool
from qd_workflows.grid_thermo import KB_EV
S = sys.argv[1]
D = "/Users/ivaninfante/Documents/University/escience/network_growth_data/cdse_acetate_p2"
net = load(D, f"{D}/thermo_acetate_pool.json", "src/qd_workflows/templates/growth/cdse_zb_acetate.yaml")
run = json.load(open(f"{S}/run_acetate.json"))
t = np.array(run["t"])
T, eps = 500.0, 2.0
g, th, _ = net.G(T, eps)
kt = KB_EV * T
m, k = net.m.astype(int), net.k.astype(int)
seed = int(np.where((m == 13) & (k == 3))[0][0])
snaps = [(5.0, "t = 5 s, before the burst"), (14.0, "t = 14 s, the burst"), (100.0, "t = 100 s, growth phase"), (600.0, "t = 600 s, end of the run")]
fig, axs = plt.subplots(2, 4, figsize=(16, 8.4), gridspec_kw={"height_ratios": [1.25, 1], "hspace": 0.38})
prof = []
for j, (ts, lab) in enumerate(snaps):
    i = int(np.argmin(abs(t - ts)))
    pl = _pool(th, run["c_MA"][i], run["c_X_left"][i])
    om = g - m * pl["mu_MA"] - k * pl["mu_MX"]
    om = om - om[seed]
    M = np.full((m.max() + 1, k.max() + 1), np.nan)
    for a, b, o in zip(m, k, om):
        if np.isnan(M[a, b]) or o < M[a, b]:
            M[a, b] = o
    ax = axs[0, j]
    lim = 40
    im = ax.imshow(M.T, origin="lower", aspect="auto", cmap="RdBu_r", norm=TwoSlopeNorm(0, -lim, lim),
                   extent=[-0.5, m.max() + 0.5, -0.5, k.max() + 0.5], interpolation="nearest")
    ax.plot(13, 3, "k*", ms=9)
    ax.set_title(f"{lab}\nC_M = {run['c_MA'][i]:.2g} M, α = {run['alpha'][i]:.2f}", fontsize=9)
    ax.set_xlabel("m (CdSe units)", labelpad=1); ax.set_ylabel("k (Cd(OAc)2 units)" if j == 0 else "")
    # Omega(m) over all states of each m
    ms = np.arange(m.min(), m.max() + 1)
    Om = np.array([-kt * np.log(np.sum(np.exp(-(om[m == x] - om[m == x].min()) / kt))) + om[m == x].min() if (m == x).any() else np.nan for x in ms])
    prof.append((lab, ms, Om))
cb = fig.colorbar(im, ax=axs[0, :], shrink=0.85, pad=0.01); cb.set_label("Ω(m, k) − Ω(seed) of the lowest state (eV)")
ax = axs[1, 0]
cnt = np.bincount(m)
ax.bar(np.arange(len(cnt)), cnt, width=1, color="#888")
ax.set_xlabel("m"); ax.set_ylabel("relaxed states per m"); ax.set_title("where the map has states", fontsize=9)
ax.axvspan(180, m.max(), color="#d9443a", alpha=.12)
for j, (lab, ms, Om) in enumerate(prof[:3] + [prof[3]]):
    pass
ax = axs[1, 1]
for (lab, ms, Om), col in zip(prof, ("#1b9e77", "#d95f02", "#7570b3", "#e7298a")):
    ax.plot(ms, Om, c=col, lw=1.2, label=lab.split(",")[0])
ax.axhline(0, c="#888", lw=.8); ax.set_ylim(-60, 40)
ax.set_xlabel("m"); ax.set_ylabel("Ω(m) − Ω(seed) (eV)"); ax.set_title("Ω(m) over all k and isomers, relative to the seed", fontsize=9)
ax.legend(fontsize=8)
ax = axs[1, 2]
tt = np.array(run["t"]); okp = tt > 0.05
ax.plot(tt[okp], np.array(run["mu_MA"])[okp] - run["mu_MA"][0], c="#2f6fdf", label="μ_MA (CdSe)")
ax.plot(tt[okp], np.array(run["mu_MX"])[okp] - run["mu_MX"][0], c="#3aa655", label="μ_MX (Cd(OAc)2)")
old = json.load(open(f"{S}/open_run_acetate.json")); to = np.array(old["t"]); oko = to > 0.05
ax.plot(to[oko], np.array(old["mu_MA"])[oko] - old["mu_MA"][0], c="#2f6fdf", ls="--", lw=.8, alpha=.6, label="μ_MA, open-loop conversion (old)")
ax.plot(to[oko], np.array(old["mu_MX"])[oko] - old["mu_MX"][0], c="#3aa655", ls="--", lw=.8, alpha=.6, label="μ_MX, open loop (old)")
for ts, _ in snaps:
    ax.axvline(ts, c="#bbb", lw=.8)
ax.set_xscale("log"); ax.set_xlabel("t (s)"); ax.set_ylabel("μ − μ(t = 1 ms) (eV)")
ax.set_title("pool chemical potentials along the run", fontsize=9); ax.legend(fontsize=8, loc="lower left")
ax2 = ax.twinx()
ax2.plot(tt[okp], np.maximum(np.array(run["c_X_left"])[okp], 1e-17), c="#3aa655", ls=":", lw=1)
ax2.plot(tt[okp], np.maximum(np.array(run["c_MA"])[okp], 1e-17), c="#2f6fdf", ls=":", lw=1)
ax2.set_yscale("log"); ax2.set_ylabel("C_M (blue), C_X (green) in solution, M (dotted)", fontsize=8, labelpad=2)
axs[1, 3].axis("off")
fig.savefig(f"{S}/img/omega_maps.png", dpi=140, bbox_inches="tight", facecolor="white")
print("seed", seed, "states", net.n)
