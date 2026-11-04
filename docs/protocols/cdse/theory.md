# Theory: free energies and the pool

Everything in the model comes back to one number per state: its **grand potential** at the current pool,

$$
\Omega_i(t)=G_i(T,\varepsilon)-m_i\,\mu_\mathrm{MA}(t)-k_i\,\mu_\mathrm{MX}(t).
$$

This page explains each piece: how $G_i$ is computed from a relaxed structure (§1–4), and how the pool sets
$\mu_\mathrm{MA}$ and $\mu_\mathrm{MX}$ (§5–6). §7 explains what $\Omega$ means.

## 1. Electronic energy: the level of theory

| item | CdCl2 grid and network | acetate grid and network |
|---|---|---|
| model | MACE-MH-1 (`macemh1model`), head `omat_pbe`, float64 | same |
| dispersion | none | D3(BJ), PBE parameters (torch-dftd 0.5.3), inside every evaluation: relaxation, MD, Hessian |
| optimiser | batched L-BFGS | same |
| Hessian, grid | analytic MACE (autograd) | central finite differences, 0.01 Å |
| Hessian, network | finite differences | finite differences where computed; learned F_vib elsewhere |

**Relaxation.** All production runs use the same relaxer:
`batched-lbfgs-0.02-plateau-0.15-20-1e-05-stall-3-1000-polish-0.005-0.3`.

- Converge to a maximum force of 0.02 eV/Å.
- A plateau is accepted below 0.15 eV/Å when the energy changes by less than 1e-5 eV/atom over 20 steps, with 3
  stall windows.
- At most 1000 steps, maximum step 0.2 Å, L-BFGS memory 100, initial Hessian 70 eV/Å².
- **Polish**: every structure within 0.3 eV of the best of its batch is re-relaxed to 0.005 eV/Å.

A `RelaxCache` (Orchestr.AI) stores every result under the exact input geometry, so reruns and restarts are free.
Its key contains the model file's absolute path, so a cache moved between clusters needs its keys rewritten.

Rules (Ivan): MACE-MH-1 throughout; LBFGS only (FIRE gave minima 0.8–2.6 eV higher, FIRE2 diverged); no classical
force fields, metadynamics or xTB; DFT (CP2K) only for COSMO of free species. MACE + D3 energies on LUMI
(OpenEquivariance) and NHR (cuEquivariance) agree to 1e-6 eV at the same geometry.

The bulk reference: MACE gives a = 6.2019 Å and E = −5.49030 eV per CdSe pair. MACE + D3 gives a = 6.1176 Å and
E = −6.10058 eV.

## 2. Free energy of one state in the gas phase

Each state is treated as an ideal-gas molecule at the standard state $c^\circ$ = 1 M, i.e. at the pressure
$p=c^\circ RT$ (Orchestr.AI `references.ideal_gas_g`):

$$
G_i^\circ(T)=E_i+\underbrace{\tfrac52k_BT-TS_\mathrm{trans}}_{\text{translation, incl. }pV}+\underbrace{\tfrac32k_BT-TS_\mathrm{rot}}_{\text{rotation}}+F_\mathrm{vib}(T)
$$

| term | formula | computed from |
|---|---|---|
| $E_i$ | MACE (+ D3) energy of the relaxed state | the relaxation |
| $S_\mathrm{trans}$ | $k_B\left[\ln\!\left(\left(\frac{2\pi Mk_BT}{h^2}\right)^{3/2}\frac{k_BT}{p}\right)+\frac52\right]$ | total mass $M$, $p=c^\circ RT$ |
| $S_\mathrm{rot}$ | $k_B\left[\ln\!\left(\frac{\sqrt\pi}{\sigma}\sqrt{\frac{T^3}{\Theta_A\Theta_B\Theta_C}}\right)+\frac32\right]$, $\Theta_X=\frac{h^2}{8\pi^2I_Xk_B}$ | principal moments $I_X$ of the relaxed geometry; σ |
| σ | rotational symmetry number, pymatgen `PointGroupAnalyzer(tolerance=0.1)` | the relaxed geometry (1 if the analysis fails) |
| $F_\mathrm{vib}$ | quasi-RRHO (below) | the Hessian, or the learned model |

**Vibrations: quasi-RRHO** (Grimme). Low modes are interpolated towards free rotors so that soft modes do not blow up
the entropy:

$$
F_\mathrm{vib}=\sum_j\left[\frac{h\nu_j}{2}+\frac{h\nu_j}{e^{x_j}-1}\right]-T\sum_j\Big[w_jS^\mathrm{harm}_j+(1-w_j)S^\mathrm{rot}_j\Big],\qquad
x_j=\frac{h\nu_j}{k_BT},\quad w_j=\frac{1}{1+(\nu_0/\nu_j)^4}
$$

$$
S^\mathrm{harm}_j=k_B\!\left[\frac{x_j}{e^{x_j}-1}-\ln(1-e^{-x_j})\right],\qquad
S^\mathrm{rot}_j=k_B\!\left[\tfrac12+\ln\sqrt{\frac{8\pi^3\mu'_jk_BT}{h^2}}\right],\qquad
\mu'_j=\frac{\mu_jB_\mathrm{av}}{\mu_j+B_\mathrm{av}},\ \mu_j=\frac{h}{8\pi^2\nu_j}
$$

with $\nu_0$ = 100 cm⁻¹ and $B_\mathrm{av}$ = 1e-44 kg m².

- **Imaginary modes**: counted when ν < −10 cm⁻¹; used as |ν|, and |ν| < 10 cm⁻¹ is raised to 10 cm⁻¹. No state is
  re-optimised along them. The CdCl2 grid has none; 121 of the 242 acetate grid states have 1–3.
- **Temperatures**: G is stored on 31 temperatures, 50 to 800 K in steps of 25 K, and interpolated linearly.
- **PHVA**: a partial Hessian (the core plus the carboxylate head, outer chain frozen) is used only when the inner
  region is smaller than the structure. Acetate has no C_β, so its inner region is the whole structure: **the acetate
  grid and network used full Hessians.** PHVA is for oleate and longer chains.

### 2.1 The learned F_vib

The network has too many large states for a Hessian each (P0 measured 2.28 ms·N² per Hessian). Ivan's option A
(2026-11-01) keeps the exact MACE energy of every state and learns only $F_\mathrm{vib}$:

- a linear pair + embedding ridge model on the relaxed geometry, species Cd, Se, O, C, H;
- translation, rotation, σ and the surface still come from the geometry;
- trained on the FD states; on acetate, trained on m ≤ 60 and tested on m 61–223, the error of a growth step is
  0.062 eV, against 0.27 eV rms between two clusters' minima of the same start;
- final model refitted on all 10,783 FD states.

The export uses FD F_vib wherever a Hessian exists and the learned F_vib only where none does ("Keep FD where it
exists", 2026-11-03); those states carry the flag `fvib: model`, or `model_ood` outside the model's domain.

## 3. A grid node: several isomers

A composition $(m,k)$ of the grid holds several states (isomers). Its free energy is their Boltzmann sum. Configurations
evaluated without a Hessian (the other quenches of the anneal) take the thermal part of the node's best state:

$$
G(m,k;T)=-k_BT\ln\sum_{i\in(m,k)}e^{-G_i(T)/k_BT},\qquad G_i\approx E_i+\big[G_\mathrm{best}(T)-E_\mathrm{best}\big]\ \text{(no Hessian)}
$$

In the network every relaxed state is its own state; isomers of one composition are joined by `iso` steps (see
[the master equation](master-equation.md)).

## 4. Solvation

Implicit solvent only, ε = 2.0 ("about octadecene"). Explicit solvent and ligand shells are not included; Ivan asked
that this remain noted as a later extension.

$$
\Delta G_\mathrm{solv}(\varepsilon)=\underbrace{\Big(1-\frac1\varepsilon\Big)E_\infty}_{\text{electrostatics}}+\underbrace{(\alpha+\gamma)(\varepsilon)\,S}_{\text{cavity + dispersion}},\qquad
(\alpha+\gamma)(\varepsilon)=2.16\ln\varepsilon-14.79\ \mathrm{mN\,m^{-1}}
$$

The cavity term is the Hille et al. (2019) β = 0 fit, refitted on 67 non-aqueous solvents (`grid_solvation` evaluates
it at $\max(\varepsilon, 1.8)$).

**Clusters and dots: σ·S.** S is computed from the geometry; the electrostatic term is taken per area:

- S = the solvent-exposed area of the union of vdW spheres × 1.4 (Cd 1.58, Se 1.90, Cl 1.75, O 1.52, C 1.70,
  H 1.20 Å), about 4 points per Å², KD-tree neighbour search.
- $E_\infty=e_\mathrm{el}S$ with $e_\mathrm{el}$ = −1.06 meV/Å², the mean over 196 CdSe campaign dots (±0.23, no
  size trend).

$$
\sigma(\varepsilon)=\Big(1-\frac1\varepsilon\Big)e_\mathrm{el}+(\alpha+\gamma)(\max(\varepsilon,1.8)),\qquad
\sigma(2)=-1.36\ \mathrm{meV\,\AA^{-2}},\qquad G_i(T,\varepsilon)=G_i^\circ(T)+\sigma(\varepsilon)S_i
$$

**Free species: CP2K COSMO.** A free molecule's charge distribution and surface are not those of a dot, so the
per-area value fitted on dots is replaced by a PBE gas-phase single point with COSMO (conductor energy $E_\infty$,
isosurface ρ₀ = 5e-4 e bohr⁻³):

| species | E_inf (eV) | S (Å²) |
|---|---|---|
| CdSe monomer | −0.327178 | 117.799 |
| CdCl2 | −0.198261 | 135.843 |
| CdSe(CdCl2)₁, own grid geometry | −0.273933 | 208.805 |
| CdSe(CdCl2)₂ | −0.459247 | 287.439 |
| CdSe(CdCl2)₃ | −0.471024 | 358.550 |
| Cd(OAc)₂ | −0.229848 | 210.278 |
| Cd(oleate)₂ | −0.13244 | 647.55 |

The acetate m = 1 complexes have no COSMO of their own and take σ·S. CP2K always runs on one whole node, all cores as
MPI ranks, with the SCF in 100-cycle pieces restarted from the wavefunction.

## 5. Precursor conversion

The Se precursor P reacts with the cadmium salt to give the CdSe supply, one salt per CdSe, second order:

$$
\mathrm{P}+\mathrm{MX}_q\longrightarrow\mathrm{MA}+\ldots,\qquad
Q=-\frac{dc_\mathrm{MX}}{dt}=k_\mathrm{conv}(T)\,c_\mathrm{MX}\,c_\mathrm P,\qquad
k_\mathrm{conv}(T)=k_\mathrm{conv}(500\,\mathrm K)\,e^{-\frac{E_a}{k_B}\left(\frac1T-\frac1{500\,\mathrm K}\right)}
$$

With $a=c_\mathrm{MX,0}$ and $b=c_\mathrm{P,0}$, the converted amount $x(t)$ is closed form, and $\alpha=x/c_\mathrm{P,0}$
is the fraction of precursor converted:

$$
x(t)=ab\,\frac{1-e^{-|a-b|k_\mathrm{conv}t}}{\max(a,b)-\min(a,b)\,e^{-|a-b|k_\mathrm{conv}t}}\quad(a\ne b),\qquad
x(t)=\frac{a^2k_\mathrm{conv}t}{1+ak_\mathrm{conv}t}\quad(a=b)
$$

- Defaults: T = 500 K, $c_\mathrm{P,0}$ = 0.025 M, $c_\mathrm{MX,0}$ = 0.03 M, $k_\mathrm{conv}$ = 1 M⁻¹ s⁻¹, $E_a$ = 0.
  $k_\mathrm{conv}$ is set by hand, not computed (Ivan: "I should be able to set this as I wish").
- $x(t)$ uses $c_\mathrm{MX,0}-x$: it does not see the salt that the dots take into their shells. With the defaults
  the conversion leaves 0.005 M of salt, less than the dots' shells need, so late in a run the salt in solution is
  used up while $x(t)$ still produces CdSe (see [Results](results.md)).

## 6. The pool: complexes and chemical potentials

The CdSe in solution is not a bare molecule: it is spread over the packages CdSe(CdX2)ₙ, n = 0…3, the grid's own m = 1
row. Two unknowns describe the pool: $C_M$, all CdSe in solution, and $C_X$, all salt in solution (free plus bound in
packages).

$$
\mu_\mathrm{MX}=G^\circ_\mathrm{MX}+\Delta G_\mathrm{solv,MX}+k_BT\ln\frac{\gamma_\mathrm{MX}\,c_\mathrm{MX,free}}{c^\circ}
$$

$$
x_n=\frac{e^{-(G_{1,n}-n\mu_\mathrm{MX})/k_BT}}{\sum_{n'}e^{-(G_{1,n'}-n'\mu_\mathrm{MX})/k_BT}},\qquad
c_\mathrm{MX,free}+C_M\langle n\rangle=C_X
$$

$$
\mu_\mathrm{MA}=k_BT\ln\frac{\gamma_\mathrm{MA}C_M}{c^\circ}-k_BT\ln\sum_n e^{-(G_{1,n}-n\mu_\mathrm{MX})/k_BT}
$$

| quantity | meaning | how it is found |
|---|---|---|
| $x_n$ | share of the dissolved CdSe in the package with n salts | Boltzmann weights of the m = 1 row at $\mu_\mathrm{MX}$ |
| $c_\mathrm{MX,free}$ | free salt | root of the salt balance, bisection in ln c |
| $\mu_\mathrm{MX}$ | chemical potential of one salt unit | its free G + COSMO + ideal mixing |
| $\mu_\mathrm{MA}$ | chemical potential of one CdSe unit | the packages' grand partition function |
| $\gamma$ | activity coefficients | 1 by default (page range 0.3–3) |

Each package then has $\Omega_{1,n}=G_{1,n}-\mu_\mathrm{MA}-n\mu_\mathrm{MX}=-k_BT\ln(C_Mx_n/c^\circ)$: the pool is in
equilibrium with itself. When $\mu_\mathrm{MX}$ falls (little free salt), the packages lose their salt and the
log-sum falls, so $\mu_\mathrm{MA}$ **rises** at the same $C_M$: bare CdSe is less stable than packaged CdSe.

Oleate enters only here and as steric coverage: Cd(oleate)2 solubility (assumed, tunable) and a maximum coverage of
about 20–25 Å² per chain. Cd(OAc)2 is nearly insoluble.

## 7. What Ω means

$\Omega_i$ is the free energy of forming state $i$ from the pool: take $m_i$ CdSe and $k_i$ salt units out of solution
(each at its chemical potential) and assemble them. Three consequences are used everywhere:

1. **Equilibrium populations.** At a fixed pool, a population in equilibrium with the solution is
   $N_i=c^\circ e^{-\Omega_i/k_BT}$. Low Ω means populated.
2. **Steps.** For a step $a\to b$ that adds $\Delta m$ CdSe and $\Delta k$ salt from solution,
   $\Delta\Omega=\Omega_b-\Omega_a=(G_b-G_a)-\Delta m\,\mu_\mathrm{MA}-\Delta k\,\mu_\mathrm{MX}$ is its free energy of
   reaction at the current pool. $\Delta\Omega<0$: downhill (growth favoured); $>0$: uphill (dissolution favoured).
3. **The map moves with the pool.** $G_i$ is fixed by the structure; the pool changes $\mu_\mathrm{MA}$ and
   $\mu_\mathrm{MX}$ in time. A high CdSe supply makes the $-m\mu_\mathrm{MA}$ term tilt the whole map towards large
   m; a falling supply tilts it back. This is how the same map describes nucleation, growth and dissolution
   ([master equation](master-equation.md), §8).

Per size, the grid and the network are summarised by

$$
\Omega(m)=-k_BT\ln\sum_{k}\sum_{i\in(m,k)}e^{-\Omega_i/k_BT},\qquad \ln S=\frac{\mu_\mathrm{MA}-g_\mathrm{bulk}}{k_BT}
$$

where S is the supersaturation against bulk CdSe ($g_\mathrm{bulk}$: bulk energy plus harmonic phonons; −0.610 eV
D3 shift for acetate).
