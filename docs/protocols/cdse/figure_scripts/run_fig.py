import json, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
S = sys.argv[1]
a = json.load(open(f"{S}/run_acetate.json")); c = json.load(open(f"{S}/run_cdcl2.json"))
t = np.array(a["t"]); tc = np.array(c["t"])
fig, axs = plt.subplots(4, 1, figsize=(8, 10), sharex=True)
phases = [(1e-3, 13, "I", "#eef3f8"), (13, 16.5, "II", "#fbe9d8"), (16.5, 200, "III", "#eaf5ea"),
          (200, 330, "IV", "#f3ecf7"), (330, 360, "V", "#f8e1e1")]
for ax in axs:
    for lo, hi, lab, col in phases:
        ax.axvspan(lo, hi, color=col, lw=0)
    ax.set_xscale("log"); ax.grid(alpha=.3, lw=.5)
for lo, hi, lab, col in phases:
    axs[0].text(np.sqrt(lo * hi) if lo > 1e-3 else 0.1, 1.01, lab, transform=axs[0].get_xaxis_transform(), ha="center", va="bottom", fontsize=9, fontweight="bold")
ax = axs[0]
ax.plot(t, a["alpha"], c="#555", label="acetate α (converted)")
ax.plot(tc, c["alpha"], c="#555", ls=":", label="CdCl2 α")
ax.set_ylabel("conversion α"); ax.legend(fontsize=8, loc="center left")
ax = axs[1]
ax.plot(t, np.maximum(a["c_MA"], 1e-14), c="#2f6fdf", label="acetate C_M (CdSe in solution)")
ax.plot(tc, np.maximum(c["c_MA"], 1e-14), c="#2f6fdf", ls=":", label="CdCl2 C_M")
ax.plot(t, np.maximum(a["n_dots"], 1e-14), c="#f0a030", label="acetate dots (m ≥ 13)")
ax.plot(tc, np.maximum(c["n_dots"], 1e-14), c="#f0a030", ls=":", label="CdCl2 dots")
ax.set_yscale("log"); ax.set_ylim(1e-14, 3e-2); ax.set_ylabel("concentration (M)"); ax.legend(fontsize=8, loc="lower right")
ax = axs[2]
m, sd = np.array(a["m_mean"]), np.array(a["m_std"])
ok = np.array(a["n_dots"]) > 1e-15
ax.plot(t[ok], m[ok], c="#a0522d", label="acetate ⟨m⟩")
ax.fill_between(t[ok], (m - sd)[ok], (m + sd)[ok], color="#a0522d", alpha=.2, lw=0)
okc = np.array(c["n_dots"]) > 1e-15
ax.plot(tc[okc], np.array(c["m_mean"])[okc], c="#a0522d", ls=":", label="CdCl2 ⟨m⟩")
ax.axhline(13, c="#888", lw=.8); ax.text(1.5e-3, 16, "seed m = 13", fontsize=7.5, color="#666")
ax.set_ylabel("dot size m (CdSe units)"); ax.legend(fontsize=8, loc="upper left")
ax = axs[3]
ax.plot(t, np.maximum(np.abs(a["nucleation"]), 1e-16), c="#6a3d9a", label="|nucleation flux| (M/s)")
for f, col in (("100", "#1b9e77"), ("111B", "#d95f02"), ("111A", "#7570b3")):
    ax.plot(t, np.maximum(np.abs(a["cdse_flux"][f]), 1e-16), c=col, lw=1, label=f"|CdSe flux by layers on {{{f}}}|")
ax.plot(t, np.array(a["frontier"]) * 1e-3, c="#999", ls="--", label="share of dots on the network edge × 1e-3")
ax.set_yscale("log"); ax.set_ylim(1e-14, 1e-2); ax.set_ylabel("flux (M/s)"); ax.legend(fontsize=7, loc="lower left", ncol=2)
ax.set_xlabel("time (s)")
fig.suptitle("Acetate default case (solid) and CdCl2 (dotted), 500 K, c_P0 0.025 M, c_MX0 0.03 M", fontsize=10, y=0.995)
fig.tight_layout()
fig.savefig(f"{S}/img/run_phases.png", dpi=150, facecolor="white")
print("ok", t[-1], a.get("stopped"))
