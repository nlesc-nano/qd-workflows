# Network growth model

A growth model that runs on the relaxed structures of the layer-path network and nothing else
(`qd_workflows.network_growth`, page `qd_workflows.network_sim`). It takes no parameter from `growth.py`: no facet rate
constants, island $\beta$, facet $\gamma$, coverage isotherms or NEB rates.

## States and steps

A dot is one relaxed state of the network (`network_export`: the tree registry's dots connected to the seed
Cd$_{16}$Se$_{13}$Cl$_6$), or, below the seed, one $(m, k)$ node of the nucleation grid ($m$ CdSe units, $k$ CdCl$_2$
units; a network state with $m \le 13$ whose composition the grid has is that grid node). The $m = 1$ row of the grid
is the solution: the CdSe(CdCl$_2$)$_n$ complexes.

| step | from the data | arriving unit |
|---|---|---|
| layer | a registry step that starts a new sheet on facet $f$ | the CdSe carrier, over the island barrier |
| attach | a registry step that adds CdSe otherwise (a sheet already started, mixed directions) | the CdSe carrier |
| shell | same CdSe framework, neighbouring $k$ | free CdCl$_2$ |
| isomer | same composition (same framework and $k$, or a registry step with no change) | none |
| grid | $(m,k) \to (m+1, k+n)$ and $(m,k) \to (m, k+1)$ between grid nodes | the complex $n$, free CdCl$_2$ |
| nucleus | two complexes $\to$ $(2, n_1+n_2)$ | the complex $n_2$ |

## Rates

Every state has the grand potential at the current solution,

$$
\Omega_i(t) = G_i(T,\varepsilon) - m_i\,\mu_\mathrm{CdSe}(t) - k_i\,\mu_\mathrm{CdCl_2}(t),
$$

with $G_i$ from the state's own frequencies (ideal gas at 1 M, its rotational symmetry number) plus the grid's cluster
solvation on its surface, and the $\mu$ from the pool (`grid_thermo.pool`). Each step $a \to b$ ($b$ the larger) has
one top, and both directions climb to it:

$$
\Omega^\ddagger_{ab} = \max\left(\Omega_a + w_s,\ \Omega_b,\ \Omega_a + w_s + B_{ab}\right), \qquad w_s = -k_BT \ln\frac{c_s}{c^\circ}
$$

$$
k_{a\to b} = \nu_{ab}\, e^{-(\Omega^\ddagger_{ab}-\Omega_a)/k_BT}, \qquad k_{b\to a} = \nu_{ab}\, e^{-(\Omega^\ddagger_{ab}-\Omega_b)/k_BT}
$$

$w_s$ places the arriving unit $s$ at the dot, so the forward rate is proportional to its concentration $c_s$ and the
backward rate is not. Their ratio is $e^{-(\Omega_b-\Omega_a)/k_BT}$ whatever the top, so at a fixed pool the network
relaxes to $N_i = c^\circ e^{-\Omega_i/k_BT}$ (detailed balance by construction).

**Island barrier** $B_{ab}$, on layer steps only. For a parent $p$ and facet family $f$ with relaxed islands (the
builder's $\{211\}$/$\{210\}$ step cuts) and ladder rungs, the lowest state with $\Delta m$ CdSe added, shell free:

$$
\Omega^*_{p,f}(\Delta m) = \min_{\text{states with }\Delta m}\left[\Delta G - \Delta m\,\mu_\mathrm{CdSe} - \Delta k\,\mu_\mathrm{CdCl_2}\right]
$$

$B$ is the highest $\Omega^*$ before the curve first drops below the parent (the highest point when it never does),
recomputed at every pool. A layer step whose parent has no islands takes the median curve of its facet and parent-size
bin ($m$ < 30, 30–60, 60–100, 100–150, ≥ 150), else of its facet; $\{110\}$ takes the median over all.

**Prefactor.** $\nu_{ab}$ is the collision rate of the unit on the step's sites at $c^\circ = 1$ M in series with
diffusion to the dot (Stokes–Einstein, radii from the computed surfaces), times one scale $s_\nu$:

$$
\nu_{ab} = s_\nu\, c^\circ \left[\frac{1}{n_\mathrm{sites}\, a_\mathrm{site} N_A \sqrt{k_BT/2\pi M_s}} + \frac{1}{4\pi (D_a + D_s)(R_a + R_s) N_A}\right]^{-1}
$$

with $n_\mathrm{sites}$ the new sheet's size for a layer step and the dot's surface over $a_\mathrm{site}$ otherwise.
Isomerisations take $k_BT/h$.

## The flask

$x(t)$ CdSe is converted from the precursor (`grid_thermo.converted`, $k_\mathrm{conv}$ Arrhenius in $T$), and the pool
holds what the dots do not:

$$
C_M = x - \sum_i m_i N_i, \qquad C_X = c_{\mathrm{MX},0} - x - \sum_i k_i N_i, \qquad \frac{dN_i}{dt} = \sum_j \left(k_{j\to i} N_j - k_{i\to j} N_i\right) + s_i
$$

($s_i$: the nuclei formed from two complexes). Each backward-Euler step solves the population (sparse LU) and the pool
together: Newton on $\ln C_M, \ln C_X$, with the exact Jacobian (the rates' derivatives by differences of the assembly,
the population's through the same LU). Cd, Se and Cl are conserved to the Newton tolerance.

**Fast mode** (`prune_eV`, default 2 eV). Within one size $\mu_\mathrm{CdSe}$ cancels, so which states of a size lie
high depends only on $\mu_\mathrm{CdCl_2}$, and that only moves between the CdCl$_2$ of the start and what is left at full
conversion. A run keeps every grid node, the seed, and each dot that comes within `prune_eV` of the lowest state of its
size somewhere in that range (checked at both ends and the middle); the populations are returned on the whole network.
At 2 eV the network keeps 3,200–4,300 of 12,474 states and runs 5–11× faster; against the exact runs (default, 4×
precursor, 560 K) the pool and the dots agree to $10^{-5}$ and the CdSe beyond $m = 13$ to within 1 %. `prune_eV = 0` is
the exact model.

## Use

```bash
python -m qd_workflows.network_export --registry registry.json --results results_by_key.json --g-cache g_cache.json --steps layer_steps.json --ladders ladders.json --islands states_index.json --island-states batch*/states.json --thermo thermo_merged.json --params templates/growth/cdse_zb.yaml --out DIR
```

```bash
python -m qd_workflows.network_sim --network DIR --thermo thermo_merged.json --params templates/growth/cdse_zb.yaml --port 8773
```

`python -m qd_workflows.network_growth run|check ...` runs the model or the detailed-balance check from the command
line.

The page's step 1 is the nucleation simulator's (`simulator.Simulator.omega_live`, run nucleation-only with the nucleus
at the seed's $m = 13$): the grid's $\Delta\Omega(m)$, $P(m,k)$, the barrier along the time and $\Delta\Omega(m,t)$ from
the precursor alone, on the same thermo, solvation and pool (`grid_thermo.pool`) as the network, so its $\mu$ at a time
equals the live map's. After a run it adds the grid's $\Omega$ at the pool the run leaves (dashed). The Becker–Döring
$J$ of the old growth model is left out: the network forms its nuclei from the complexes itself.

## First result (the template's synthesis: 500 K, 0.025 M Se precursor, 0.03 M CdCl$_2$, 600 s)

The 12,474-state network (12,370 relaxed dots, 132 grid nodes, 74,125 steps) runs in about 3 minutes (204 implicit
steps, 278 sparse factorisations). At a fixed pool the stationary state matches $c^\circ e^{-\Omega/k_BT}$ to $10^{-13}$.

Nothing grows past the seed. The small clusters ($m$ 10–12) and the seeds take up the CdSe as it converts and hold
$C_M$ near $7.6\times10^{-4}$ M. At that pool the cheapest route from the seed to $m \ge 60$ climbs 2.0 eV, and its top
is a relaxed state at $m \approx 46$, not an island barrier (2.03 eV with the islands switched off as well). The route
falls to 0.86 eV at $C_M = 3\times10^{-3}$ M and 0.79 eV at $10^{-2}$ M (0.56 and 0.20 eV without islands). The lowest
states at these pools are dots of $m$ 165–198, about 10 eV below the seed, out of reach. Whether the $m$ 20–50 hump
is real (a magic-size trap) or a gap of the network (only the lineage's intermediates were relaxed there) decides
whether this model grows dots at all.

No condition in the page's range changes that: ten times faster conversion, four times the precursor ($c_{P,0} = 0.1$
M, $c_{\mathrm{MX},0} = 0.12$ M) or 560 K all end at $\langle m \rangle = 13$ after 600 s, because the seeds and the
$m$ 10–12 clusters hold $C_M$ below about $10^{-3}$ M whatever the supply. At that pool the lowest state of every $m$
from 14 to 71 lies 0.5–2 eV above the seed (none below it before $m \approx 72$), and the best route climbs 2.15 eV,
about 1.1 eV more than the lowest states alone would need: the low states of different lineage edges are not one
layer apart. The window is also thin in shells: 1,530 of the 1,953 Se frameworks at $m$ 20–50 have one CdCl$_2$
count only.

## Limits

- The network stops at $m \approx 223$ (about 2.9 nm); dots that reach a state with no further growth step stay there
  (the page reports their share as "on the network's edge").
- Island curves exist for 30 (parent, facet) pairs; a median curve predicts a single parent's barrier to about
  0.8–1.1 eV (leave-one-out), so individual layer steps are uncertain and the result must be checked for its
  sensitivity to them.
- Implicit solvent only; no explicit ligands.
