# Thermochemistry

The free energy of every species in the pipeline, computed after the Hessian (GPU stage: steps
`hessian`, `stability`, `detachment`). Every reaction free energy and equilibrium in
[Reactions and equilibria](equilibria.md) is built from these.

## Notation: what MA and MX$_q$ are

A library dot is a charge-balanced binary nanocrystal: a metal **M** (Cd, Zn, Pb, Hg, In, Ga, Cs…),
a core anion **A** (Se, S, Te, P, As…) and a monovalent surface ligand **X** (Cl, Br, I). With the
metal charge $q$, it is split into two kinds of neutral units:

| Unit | What it is | Examples |
|---|---|---|
| **MA** | one formula unit of the bulk semiconductor: a cation–anion pair of the core | CdSe, InAs, PbS |
| **MX$_q$** | one metal-halide unit on the surface (a Z-type ligand) | CdCl$_2$ ($q$ = 2), InCl$_3$ ($q$ = 3), PbBr$_2$ |

so that a dot with $n$ MA units and $m$ MX$_q$ units is written

$$
[\mathrm{MA}]_n(\mathrm{MX}_q)_m
$$

Example: Cd$_{16}$Se$_{13}$Cl$_6$ = $[\mathrm{CdSe}]_{13}(\mathrm{CdCl_2})_3$: 13 of the 16 Cd pair with
the 13 Se, the other 3 Cd carry the 6 Cl. The InAs dot of the reaction example,
In$_{52}$As$_{35}$Cl$_{51}$(CH$_3$NH$_2$)$_{21}$, is $[\mathrm{InAs}]_{35}(\mathrm{InCl_3})_{17}(\mathrm{CH_3NH_2})_{21}$.

Three **reference species** give the energy scale for growing, dissolving or stripping a dot:

| Reference | Physical meaning | Used for |
|---|---|---|
| **bulk MA** (crystal) | what the dot would become if it grew without limit: the thermodynamic sink of the core | decomposition free energy, MA solubility |
| **MA monomer** (diatomic molecule, e.g. CdSe) | the core taken apart into isolated units; not a real solution species | binding (cohesive) free energy, coupled equilibria |
| **MX$_q$ monomer** (molecule, e.g. linear CdCl$_2$, trigonal InCl$_3$) | a Z-type ligand once it has left the surface; it stays in solution as a molecular complex | every desorption and binding reaction |

Energies $E$ are MACE-MH-1 energies at MACE-MH-1 minima. $k_BT$ at 298.15 K is 25.69 meV.

## Free energy of one molecule-like species

The dot, its ligand-stripped states and the two monomers are each treated as an ideal solute at the
1 M standard state ($c^\circ = 1$ mol L⁻¹), with ideal-gas translation, rigid-rotor rotation and
harmonic vibrations from the MACE-MH-1 Hessian:

$$
G_i^\circ(T) = E_i + F_{\mathrm{vib},i}(T) + H_{\mathrm{trans}} + H_{\mathrm{rot}} - T\left(S_{\mathrm{trans}} + S_{\mathrm{rot}}\right)
$$

### Translation

Evaluated at the volume per molecule of a 1 M solution, i.e. at the pressure $p = c^\circ R T$ in
ASE's `IdealGasThermo` ($m$: the molecule's mass):

$$
H_{\mathrm{trans}} = \tfrac52 k_B T, \qquad
S_{\mathrm{trans}} = k_B\left[\ln\!\left(\left(\frac{2\pi m k_B T}{h^2}\right)^{3/2}\frac{k_B T}{p}\right) + \frac52\right]
$$

### Rotation

Rigid rotor with the moments of inertia $I_X$ of the minimum and the rotational symmetry number
$\sigma$ of its point group (12 for a T$_d$ dot, 2 for linear CdCl$_2$):

$$
H_{\mathrm{rot}} = \tfrac32 k_B T, \qquad
S_{\mathrm{rot}} = k_B\left[\ln\!\left(\frac{\sqrt{\pi}}{\sigma}\sqrt{\frac{T^3}{\Theta_A\Theta_B\Theta_C}}\right) + \frac32\right], \qquad
\Theta_X = \frac{h^2}{8\pi^2 I_X k_B}
$$

For linear molecules $H_{\mathrm{rot}} = k_BT$ and $S_{\mathrm{rot}} = k_B[\ln(T/\sigma\Theta) + 1]$.

### Vibration

Harmonic internal energy with the zero-point energy, and Grimme's quasi-RRHO entropy, so that the
soft modes of large dots (ligand torsions, surface rocking) do not inflate $S$. With
$x_j = h\nu_j / k_B T$ for each of the $3N - 6$ frequencies (imaginary modes taken as $|\nu_j|$, modes
below 10 cm⁻¹ raised to 10 cm⁻¹):

$$
F_{\mathrm{vib}} = U_{\mathrm{vib}} - T S_{\mathrm{vib}}, \qquad
U_{\mathrm{vib}} = \sum_j \left[\frac{h\nu_j}{2} + \frac{h\nu_j}{e^{x_j} - 1}\right]
$$

Each mode's entropy switches smoothly from the harmonic oscillator (above ~100 cm⁻¹) to a free
rotor (below):

$$
S_{\mathrm{vib}} = \sum_j \Big[w_j\, S^{\mathrm{harm}}_j + (1 - w_j)\, S^{\mathrm{rot}}_j\Big], \qquad
w_j = \frac{1}{1 + (\nu_0/\nu_j)^4}, \quad \nu_0 = 100\ \mathrm{cm^{-1}}
$$

$$
S^{\mathrm{harm}}_j = k_B\left[\frac{x_j}{e^{x_j} - 1} - \ln\!\left(1 - e^{-x_j}\right)\right], \qquad
S^{\mathrm{rot}}_j = k_B\left[\frac12 + \ln\sqrt{\frac{8\pi^3 \mu'_j k_B T}{h^2}}\right]
$$

with the free-rotor moment of inertia of a mode, damped towards an average molecular value:

$$
\mu_j = \frac{h}{8\pi^2 \nu_j}, \qquad \mu'_j = \frac{\mu_j B_{\mathrm{av}}}{\mu_j + B_{\mathrm{av}}}, \qquad B_{\mathrm{av}} = 10^{-44}\ \mathrm{kg\,m^2}
$$

## Free energy of bulk MA

The bulk crystal neither translates nor rotates: its free energy per formula unit is the
electronic energy plus the phonon free energy. The bulk is built from the record's CIF and
cell-relaxed with MACE-MH-1; the force constants come from finite displacements (0.01 Å) in a
supercell (ASE `Phonons`), and the phonon density of states $g(\varepsilon)$ is sampled on a
24 × 24 × 24 $q$-mesh. Per formula unit ($N_\mathrm{fu}$ units in the cell):

$$
G^{\mathrm{bulk}}_{\mathrm{MA}}(T) = E^{\mathrm{bulk}}_{\mathrm{MA}} + \frac{1}{N_\mathrm{fu}}\int g(\varepsilon)\left[\frac{\varepsilon}{2} + k_BT\ln\!\left(1 - e^{-\varepsilon/k_BT}\right)\right] d\varepsilon
$$

The reference is computed once per material and MACE head and cached (`props.refs`), so every
dot of a material uses the same bulk.

## From free energies to chemical potentials

In solution each molecule-like species gains its solvation free energy
([Solvation](solvation.md)), and its chemical potential depends on its concentration $c_i$
(ideal dilute solution, $c^\circ = 1$ M):

$$
G_i^{\mathrm{sol}}(T, \varepsilon) = G_i^\circ(T) + \Delta G_{\mathrm{solv},i}(\varepsilon), \qquad
\mu_i(T, \varepsilon, c_i) = G_i^{\mathrm{sol}}(T, \varepsilon) + k_B T \ln\frac{c_i}{c^\circ}
$$

One decade of concentration moves $\mu_i$ by $k_BT\ln 10$ = 59.2 meV at 298 K.

Two references are treated differently:

- **Bulk MA is a solid.** A pure solid has activity 1 whatever its amount, so it gets neither a
  concentration term nor a solvation term: $\mu^{\mathrm{bulk}}_{\mathrm{MA}} = G^{\mathrm{bulk}}_{\mathrm{MA}}(T)$.
- **The MA monomer can be shifted.** The diatomic MA molecule is a computational reference: in a
  real synthesis MA comes from molecular precursors (e.g. cadmium oleate and a selenium phosphine),
  which are far more stable than an isolated CdSe molecule. The dashboards therefore add an
  adjustable shift $\Delta\mu_{\mathrm{prec}}$, negative for a stabilised precursor:

$$
\mu_{\mathrm{MA}} = G^{\mathrm{sol}}_{\mathrm{MA}}(T, \varepsilon) + \Delta\mu_{\mathrm{prec}} + k_B T \ln\frac{c_{\mathrm{MA}}}{c^\circ}
$$

These chemical potentials enter the reactions of the dot: decomposition into bulk, binding from
the monomers, ligand desorption and the coupled equilibria ([Reactions and equilibria](equilibria.md)).
