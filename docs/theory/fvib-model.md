# Learned vibrational free energy

The network growth model ([network-growth.md](network-growth.md)) needs G for every state. For dots larger than the
relaxed network ($m > 223$) the plan is to keep a real MACE relaxation for every builder structure, so the energy is
computed, and to replace only the finite-difference Hessian, the $N^2$ part of the cost, by a learned vibrational free
energy (`qd_workflows.fvib_model`). Translation, rotation and the solvation $\sigma S$ come from the relaxed geometry
as before:

$$
G_i(T,\varepsilon) = E_i^\mathrm{MACE} + \hat F_\mathrm{vib}(\mathbf R_i;T) + G_\mathrm{tr+rot}(\mathbf R_i;T) + \sigma(\varepsilon)\,S(\mathbf R_i)
$$

## Target and model

The target is the network's own $F_\mathrm{vib}(T)$: the quasi-RRHO free energy of `orchestr_ai.qd.steps.hessian.thermo`
(harmonic $U$, Grimme entropy, modes floored at 10 cm$^{-1}$) from the FD frequencies of the relaxed dots, on the
grid's temperatures. Only this network's data are used (MACE-MH-1 omat_pbe float64, relax and polish).

The model is per-atom additive on the relaxed geometry and linear in its descriptor, with one weight set per
temperature:

$$
\hat F_\mathrm{vib}(T) = \sum_{i<j} \sum_n w_{s_i s_j,n}(T)\,\phi_n(r_{ij})\,f_c(r_{ij}) + \sum_i \sum_p c_{s_i,p}(T)\,g_p(\rho_i) + \sum_s b_s(T)\,N_s + b_0(T)
$$

$\phi_n$ are Gaussians of width 0.25 Å centred every 0.25 Å from 2 to 7 Å, $f_c$ a cosine cutoff at 7.5 Å,
$\rho_i = \sum_j e^{-((r_{ij}-2.7)/0.3)^2}$ the atom's smooth coordination and $g_p \in \{\rho, \rho^2, \sqrt\rho\}$.
The weights come from ridge regression on standardised features; the strength is chosen on the lineage-edge folds.

**Domain.** A feature present in fewer than 20 training dots gets no weight. A structure is out of domain if its
bonding shell (pairs up to 3.5 Å) has such a bond type, or more of a rare one (present in fewer than 5 % of the training
dots) than any training dot. Being larger does not make a dot out of domain. Out-of-domain structures need a Hessian.
In the 17,564 dots these are a few Se–Cl contacts at 2.2 Å (6 dots) and the largest Se–Se dimer counts (Se–Se near
2.5 Å occurs in 94 dots).

## Validation (17,564 relaxed dots, m 3–223)

The error that matters is that of a step free energy $\Delta F_\mathrm{vib}$ between neighbouring dots, since a step
error $\delta$ changes a rate by $e^{\delta/k_BT}$ (0.1 eV is a factor of 10 at 500 K). The steps are those of
`layer_steps.json` between relaxed dots (61,695 that add CdSe) and the shell links of the export (3,274: the lowest-E
isomer of each framework and CdCl$_2$ count to that of the next count). A step is scored when both of its dots are in
the test set. *Within size* is the error after removing the test set's mean residual at each $m$; a residual that
drifts with $m$ acts as a shift of $\mu_\mathrm{CdSe}$.

| split (500 K) | train / test dots | growth step MAE | p95 | > 0.1 eV | shell step MAE | within size | absolute | drift (eV per 100 m) |
|---|---|---|---|---|---|---|---|---|
| train $m<30$, test $m \ge 30$ | 1,845 / 15,285 | 0.035 | 0.100 | 5.1 % | 0.013 | 0.037 | 0.106 | −0.21 |
| train $m<60$, test $m \ge 60$ | 7,752 / 9,809 | 0.027 | 0.074 | 1.5 % | 0.012 | 0.026 | 0.041 | −0.07 |
| train $m<100$, test $m \ge 100$ | 12,126 / 5,438 | 0.029 | 0.080 | 2.0 % | 0.011 | 0.024 | 0.028 | −0.01 |
| lineage edges, 5 folds (mean) | 14,051 / 3,505 | 0.025 | 0.074 | 2.1 % | 0.016 | 0.019 | 0.022 | 0.00 |
| random 80/20 (reference) | 14,042 / 3,513 | 0.021 | 0.061 | 0.5 % | 0.016 | 0.018 | 0.018 | 0.00 |

Errors in eV; out-of-domain test dots are counted and not scored (434 for the $m<30$ cut, at most 9 otherwise). From 400
to 600 K the growth step MAE rises by about 30 % (0.023 to 0.031 eV for the $m<60$ cut).

Trained on $m<60$, the model predicts dots three times larger with a mean step error of 0.027 eV, a factor of 1.9 in a
rate at 500 K. Its absolute error drifts by −0.07 eV per 100 CdSe. Groups of isomers of one framework share their
residual (the vertical bands in `residuals.png`), so the step errors inside a group are smaller than the absolute
ones. Dots with several imaginary modes (quasi-RRHO counts them as $|\nu|$) are among the largest residuals.

## Use

```bash
python -m qd_workflows.fvib_model fit --results results_by_key.json --registry registry.json --steps layer_steps.json --thermo thermo_merged.json --out DIR
```

```bash
python -m qd_workflows.fvib_model predict --model DIR/model.json relaxed.xyz
```

`fit` writes `model.json` (trained on every dot), `validation.json` (every split at 400, 500 and 600 K, and the ridge
scan), `residuals.png` (residual against $m$ and the step-error histogram for the size cuts) and caches the
descriptors in `dataset.npz`. It takes about 30 s on four processes. `predict` prints $\hat F_\mathrm{vib}$ and flags
out-of-domain structures.

## Next

Validation above the training range: about 100 dots at $m$ 250–500 relaxed with a full FD Hessian at the network's
settings (on the go-ahead). The model is then used for the relaxed states beyond $m = 223$, with each state's
$F_\mathrm{vib}$ marked as computed or predicted in the export.
