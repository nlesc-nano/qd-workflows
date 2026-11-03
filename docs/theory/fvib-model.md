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
\begin{aligned}
\hat F_\mathrm{vib}(T) = {} & \sum_{i<j} \sum_n w_{s_i s_j,n}(T)\,\phi_n(r_{ij})\,f_c(r_{ij}) + \sum_i \sum_p c_{s_i,p}(T)\,g_p(\rho_i) \\
& + \sum_i \sum_{j<k \in \mathrm{sh}(i)} \sum_c a_{s_i,c}(T)\,\psi_c(\cos\theta_{jik}) + \sum_i d_{s_i,n_i}(T) + \sum_s b_s(T)\,N_s + b_0(T)
\end{aligned}
$$

$\phi_n$ are Gaussians of width 0.25 Å centred every 0.25 Å from 2 to 7 Å, $f_c$ a cosine cutoff at 7.5 Å,
$\rho_i = \sum_j e^{-((r_{ij}-2.7)/0.3)^2}$ the atom's smooth coordination and $g_p \in \{\rho, \rho^2, \sqrt\rho\}$.
The angle term runs over pairs of first-shell neighbours ($r < 3.2$ Å) of each atom: $\psi_c$ are Gaussians of width
0.12 on the cosine of the bond angle, centred every 1/6 from −1 to 1, and $n_i$ is the atom's first-shell count (0 to
6, with 6 for more), giving 66 features.
The weights come from ridge regression on standardised features; the strength is chosen on the lineage-edge folds.

**Domain.** A feature present in fewer than 20 training dots gets no weight. A structure is out of domain if its
bonding shell (pairs up to 3.5 Å) has such a bond type, or more of a rare one (present in fewer than 5 % of the training
dots) than any training dot. Being larger does not make a dot out of domain. Out-of-domain structures need a Hessian.
In the 17,564 dots these are a few Se–Cl contacts at 2.2 Å (6 dots) and the largest Se–Se dimer counts (Se–Se near
2.5 Å occurs in 94 dots).

## Validation (17,564 relaxed dots, m 3–223, and the 99 P1 dots)

The error that matters is that of a step free energy $\Delta F_\mathrm{vib}$ between neighbouring dots, since a step
error $\delta$ changes a rate by $e^{\delta/k_BT}$ (0.1 eV is a factor of 10 at 500 K). The steps are those of
`layer_steps.json` between relaxed dots (61,695 that add CdSe) and the shell links of the export (3,274: the lowest-E
isomer of each framework and CdCl$_2$ count to that of the next count). A step is scored when both of its dots are in
the test set. *Within size* is the error after removing the test set's mean residual at each $m$; a residual that
drifts with $m$ acts as a shift of $\mu_\mathrm{CdSe}$.

| split (500 K) | train / test dots | growth step MAE | p95 | > 0.1 eV | shell step MAE | within size | absolute | drift (eV per 100 m) |
|---|---|---|---|---|---|---|---|---|
| train $m<30$, test $m \ge 30$ | 1,845 / 15,373 | 0.029 | 0.082 | 2.5 % | 0.013 | 0.031 | 0.096 | −0.20 |
| train $m<60$, test $m \ge 60$ | 7,752 / 9,908 | 0.024 | 0.066 | 0.9 % | 0.011 | 0.024 | 0.040 | −0.07 |
| train $m<100$, test $m \ge 100$ | 12,126 / 5,537 | 0.025 | 0.066 | 0.9 % | 0.010 | 0.022 | 0.027 | 0.00 |
| lineage edges, 5 folds (mean) | 14,130 / 3,525 | 0.021 | 0.059 | 0.5 % | 0.014 | 0.016 | 0.018 | 0.00 |
| random 80/20 (reference) | 14,111 / 3,540 | 0.019 | 0.052 | 0.2 % | 0.013 | 0.015 | 0.016 | 0.00 |

Errors in eV, at λ = 1e-7 chosen on the lineage folds. The P1 dots are in the training and test sets but have no
steps to the network. Out-of-domain test dots are counted and not scored (445 for the $m<30$ cut, at most 12
otherwise). From 400 to 600 K the growth step MAE rises by about 40 % (0.020 to 0.028 eV for the $m<60$ cut). Without
the angle term (the first version of the model) the growth step MAEs were 0.035, 0.027, 0.029, 0.025 and 0.021 eV.

Trained on $m<60$, the model predicts dots three times larger with a mean step error of 0.024 eV, a factor of 1.7 in a
rate at 500 K. Its absolute error drifts by −0.07 eV per 100 CdSe. Groups of isomers of one framework share their
residual (the vertical bands in `residuals.png`), so the step errors inside a group are smaller than the absolute
ones. Dots with several imaginary modes (quasi-RRHO counts them as $|\nu|$) are among the largest residuals.

## Use

```bash
python -m qd_workflows.fvib_model fit --results results_by_key.json --registry registry.json --steps layer_steps.json --thermo thermo_merged.json --out DIR [--extra states.json ...]
```

```bash
python -m qd_workflows.fvib_model predict --model DIR/model.json relaxed.xyz
```

`--extra states.json ...` adds relaxed states from `relax_states.py` runs outside the registry (the P1 dots). `fit`
writes `model.json` (trained on every dot), `validation.json` (every split at 400, 500 and 600 K, and the ridge
scan), `residuals.png` (residual against $m$ and the step-error histogram for the size cuts) and caches the
descriptors in `dataset.npz`. It takes about 40 s on four processes. `predict` prints $\hat F_\mathrm{vib}$ and flags
out-of-domain structures.

## Validation above the training range (P1, 99 dots, m 253–494)

99 builder dots at $m$ 253–494 ($N$ 596–1,228) were relaxed with an FD Hessian at the network's settings on LUMI
(about 44 GCD-h). They comprise 39 dots on chains along the lineage edges below the lowest-$\Omega$ states at
$m$ 150–223, 10 shell variants ($k_0 \pm 1$) and 25 random layer pairs, with 63 steps inside the set. All 99
converged; 25 have one to three soft imaginary modes (−10 to −23 cm$^{-1}$).

The first model (pair + embedding, no angles) trained on every $m \le 223$ dot:

| T | absolute MAE | mean | drift (eV per 100 m) | step MAE | p95 | > 0.1 eV |
|---|---|---|---|---|---|---|
| 400 K | 0.038 | +0.002 | +0.010 | 0.038 | 0.109 | 8 % |
| 500 K | 0.046 | +0.002 | +0.009 | 0.045 | 0.133 | 10 % |
| 600 K | 0.053 | +0.001 | +0.008 | 0.052 | 0.156 | 13 % |

There is no bias, no drift and no out-of-domain dot. The step error is 1.7–2 times the in-range error at the same
step size ($\Delta m$ 3–20: about 0.04 against 0.022 eV), about what a per-atom additive error growing as
$\sqrt N$ gives. The largest errors are on the chain steps of 25–38 CdSe. The shell steps have 0.027 eV, and steps
with an imaginary-mode end have 0.057 eV against 0.037 eV.

**Bond angles.** A histogram of first-shell bond-angle cosines per centre species, together with the first-shell
coordination counts, lowers the in-range growth-step error by 10–20 % on every size cut ($m<60$: 0.027 to
0.024 eV). On P1 it gives 0.043 eV per step, and the share of steps above 0.1 eV drops from 10 % to 5 %. Adding the
P1 dots to the training set (other lineage groups only) lowers the P1 absolute error from 0.045 to 0.039 eV but not
the step error. The module now includes both: the angle features, and `--extra` for relaxed states outside the
registry.

## Use beyond the network

The model is meant for the relaxed states beyond $m = 223$, with each state's $F_\mathrm{vib}$ marked as computed
or predicted in the export. Expect about 0.04 eV per step at 500 K for dots of 600–1,200 atoms. That is a factor of
about 2.5 in a rate, with 5–10 % of steps above 0.1 eV. Where single rates matter (the lowest-$\Omega$ path), an FD
Hessian on those states removes the error.

## Acetate-capped dots

`--species Cd,Se,O,C,H` (`fvib_model.SPECIES_ACETATE`) builds the same descriptor on the acetate elements. The bonds
inside an acetate (C–H, C–O, C–C) lie below the 2 Å pair grid. Their modes are nearly the same in every acetate, so
the composition terms absorb them. Cd–O and the crowding of neighbouring acetates fall in the pair, angle and
coordination features.

Test (NHR, 3 Nov 2026): 10,473 acetate network dots with FD Hessians, from the P0 runs ($m$ 7–13), P1 ($m$ 14–60) and
the finished part of P2 ($m$ 61–223, 2,620 dots), λ 1e-7 from the lineage folds. At 500 K:

| split | growth step MAE | > 0.1 eV | absolute MAE | drift (eV per 100 m) | out of domain |
|---|---|---|---|---|---|
| train $m < 61$ (P0 + P1), test P2 | 0.062 | 19 % | 0.054 | −0.018 | 32 |
| train $m < 100$ | 0.066 | 23 % | 0.052 | +0.051 | 3 |
| lineage folds | 0.050 | 12 % | 0.039 | −0.002 | 9 |
| random 80/20 | 0.047 | 9 % | 0.036 | +0.007 | 9 |

The error is about 2.5 times that for CdCl2 dots, but it hardly grows with $m$: on P2 the absolute MAE goes from
0.044 eV at $m$ 61–89 to 0.064 eV at $m$ 160–223, with a mean of −0.01 to −0.03 eV. It also grows with the number
of imaginary modes, from 0.045 eV with none to 0.069 eV with seven or more. Those soft acetate modes make the FD
reference itself noisy. Steps of 1–2 CdSe have 0.057 eV.

For comparison, the same start relaxed on two machines lands 0.27 eV rms apart in $E$ (acetate network P0). A
step error of 0.06 eV in $F_\mathrm{vib}$ is small against that noise, so the network's acetate states can be relaxed
without an FD Hessian. Out-of-domain states (about 1 % with the $m \le 60$ training set) still need a Hessian.
