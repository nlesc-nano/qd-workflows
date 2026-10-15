# Thermochemistry

A dot is a charge-balanced binary $[\mathrm{MA}]_n(\mathrm{MX}_q)_m$: $n$ MA units (the core) and $m$ Z-type ligand units $\mathrm{MX}_q$. Energies $E$ are MACE-MH-1 energies at MACE-MH-1 minima; $k_BT$ at 298.15 K is 25.69 meV.

Computed after the Hessian (GPU stage, `hessian`, `stability`, `detachment`) and used by every free energy in the pipeline.

## Free energy of one species

Every molecule-like species (the dot and its stripped states, the MA and $\mathrm{MX}_q$ monomers) is
an ideal solute at the 1 M standard state, with rigid-rotor rotation and harmonic vibrations from the
MACE-MH-1 Hessian:

$$
G_i^\circ(T) = E_i + F_{\mathrm{vib},i}(T) + H_{\mathrm{trans}} + H_{\mathrm{rot}} - T\left(S_{\mathrm{trans}} + S_{\mathrm{rot}}\right)
$$

### Translation

 at $c^\circ = 1$ M, i.e. at the pressure $p = c^\circ R T$ (ASE `IdealGasThermo`):

$$
H_{\mathrm{trans}} = \tfrac52 k_B T, \qquad
S_{\mathrm{trans}} = k_B\left[\ln\!\left(\left(\frac{2\pi m k_B T}{h^2}\right)^{3/2}\frac{k_B T}{p}\right) + \frac52\right]
$$

### Rotation

 (rigid rotor, rotational symmetry number $\sigma$ from the point group):

$$
H_{\mathrm{rot}} = \tfrac32 k_B T, \qquad
S_{\mathrm{rot}} = k_B\left[\ln\!\left(\frac{\sqrt{\pi}}{\sigma}\sqrt{\frac{T^3}{\Theta_A\Theta_B\Theta_C}}\right) + \frac32\right], \qquad
\Theta_X = \frac{h^2}{8\pi^2 I_X k_B}
$$

(linear molecules: $H_{\mathrm{rot}} = k_BT$, $S_{\mathrm{rot}} = k_B[\ln(T/\sigma\Theta) + 1]$).

### Vibration

 harmonic $U$, quasi-RRHO entropy (Grimme), so that soft modes of large dots do not
inflate $S$. With $x_j = h\nu_j / k_B T$ (imaginary modes as $|\nu_j|$, modes below 10 cm⁻¹ raised to 10):

$$
F_{\mathrm{vib}} = U_{\mathrm{vib}} - T S_{\mathrm{vib}}, \qquad
U_{\mathrm{vib}} = \sum_j \left[\frac{h\nu_j}{2} + \frac{h\nu_j}{e^{x_j} - 1}\right]
$$

$$
S_{\mathrm{vib}} = \sum_j \Big[w_j\, S^{\mathrm{harm}}_j + (1 - w_j)\, S^{\mathrm{rot}}_j\Big], \qquad
w_j = \frac{1}{1 + (\nu_0/\nu_j)^4}, \quad \nu_0 = 100\ \mathrm{cm^{-1}}
$$

$$
S^{\mathrm{harm}}_j = k_B\left[\frac{x_j}{e^{x_j} - 1} - \ln\!\left(1 - e^{-x_j}\right)\right], \qquad
S^{\mathrm{rot}}_j = k_B\left[\frac12 + \ln\sqrt{\frac{8\pi^3 \mu'_j k_B T}{h^2}}\right]
$$

$$
\mu_j = \frac{h}{8\pi^2 \nu_j}, \qquad \mu'_j = \frac{\mu_j B_{\mathrm{av}}}{\mu_j + B_{\mathrm{av}}}, \qquad B_{\mathrm{av}} = 10^{-44}\ \mathrm{kg\,m^2}
$$

### Bulk MA

The growth sink: cell-relaxed bulk from the record's CIF, phonons by finite displacements
in a supercell, per formula unit:

$$
G_{\mathrm{MA}}^{\mathrm{bulk}}(T) = E^{\mathrm{bulk}}_{\mathrm{MA}} + F^{\mathrm{phonon}}_{\mathrm{vib}}(T)
$$

## Chemical potentials

$$
G_i^{\mathrm{sol}}(T, \varepsilon) = G_i^\circ(T) + \Delta G_{\mathrm{solv},i}(\varepsilon), \qquad
\mu_i = G_i^{\mathrm{sol}} + k_B T \ln\frac{c_i}{c^\circ}
$$

The MA monomer can be shifted by $\Delta\mu_{\mathrm{prec}}$ (a stabilised molecular precursor):
$\mu_{\mathrm{MA}} = G^{\mathrm{sol}}_{\mathrm{MA}} + \Delta\mu_{\mathrm{prec}} + k_BT\ln c_{\mathrm{MA}}$. Bulk MA is a
solid: no solvation, no concentration term.

