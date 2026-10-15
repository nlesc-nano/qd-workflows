# Solvation

The cubes are computed in the CP2K chain (`esp`) and turned into energies in the CPU stage
(`solvation`). It is the same model as `qdw reaction` ([Reaction free energies](../reactions.md)).

$$
\Delta G_{\mathrm{solv}}(\varepsilon) = \underbrace{\left(1 - \frac{1}{\varepsilon}\right) E_\infty}_{\text{electrostatics}}
\;+\; \underbrace{(\alpha + \gamma)(\varepsilon)\, S}_{\text{cavity + dispersion}}
$$

Both terms are continuous in $\varepsilon$: one number per species ($E_\infty$, $S$) serves every
dielectric constant on the dashboards' sliders.

## Electrostatics: COSMO on the DFT potential

The solvent is a conductor screened by $(1 - 1/\varepsilon)$ (COSMO). The cavity is the
electron-density isosurface

$$
\rho(\mathbf r) = \rho_0 = e^{-7.6} = 5\times10^{-4}\ e\,\mathrm{bohr}^{-3},
$$

the centre of the SCCS switching function of Andreussi et al. with the parameters of Hille et al.
The surface is found by marching cubes on the CP2K density cube and divided into elements $a$ on
the mesh vertices (each with a third of the area of its triangles, about one element per
0.3 Å²; coarser for large dots so that there are at most 20,000). The screening charges
$\sigma_a$ cancel the potential of the solute on the surface, and $E_\infty$ is their interaction
energy with it:

$$
\sum_b A_{ab}\,\sigma_b = -\phi_a, \qquad
A_{ab} = \frac{1}{|\mathbf r_a - \mathbf r_b|}\ (a \ne b), \quad A_{aa} = 1.07\sqrt{\frac{4\pi}{s_a}}, \qquad
E_\infty = \tfrac12 \sum_a \sigma_a\, \phi_a
$$

($s_a$: the element's area). The potential $\phi$ is the CP2K PBE electrostatic (Hartree +
nuclear) potential of the gas-phase single point, interpolated on the elements; CP2K's
`V_HARTREE_CUBE` has the opposite sign of the usual convention, which the reader corrects. No
atomic charges are involved.

### Why not atomic charges or spheres

- Atom-centred spheres are the wrong boundary for an ionic nanocrystal surface: they leave
  crevices between surface atoms where the potential is strong. The conductor energy of an InAs
  dot changed 4× between 1.1 and 1.4 × the van der Waals radii.
- CP2K's default Hirshfeld charges overestimate dipoles 2–30× (AsH$_3$: 5.3 D against 0.19 D), so
  any charge model built on them (Generalized Born) inherits large errors.
- The density isosurface needs no radii and follows the reconstructed surface.

## Cavity and dispersion

Hille et al., *J. Chem. Phys.* **150**, 041710 (2019), fitted the non-electrostatic terms of SCCS
on about 6,000 experimental partition coefficients in 67 non-aqueous solvents. Their fit with
$\beta = 0$ (surface term only; their main fit has a volume term that gives −18 eV for a
5,700 Å³ dot, so it is for molecules only), refitted here over their Table S2:

$$
(\alpha + \gamma)(\varepsilon) = 2.16\,\ln\varepsilon - 14.79\ \ \mathrm{mN\,m^{-1}}
\qquad (\varepsilon \ge 1.8,\ R^2 = 0.71)
$$

times the area $S$ of the same isosurface (1 mN m⁻¹ = $6.24\times10^{-5}$ eV Å⁻²). The fit has no
solvent below $\varepsilon = 1.8$; there the term goes to zero with $(1 - 1/\varepsilon)$:

$$
(\alpha+\gamma)\,S\Big|_{\varepsilon < 1.8} = (\alpha+\gamma)(1.8)\,S\;\frac{1 - 1/\varepsilon}{1 - 1/1.8}
$$

## Which structures are computed

| Species | $E_\infty$ and $S$ |
|---|---|
| intact dot, $k = 0$ | DFT cubes (restart of the geo_opt wavefunction) |
| fully stripped dot, $k = K$ | DFT cubes at the MACE geometry of the last ladder state |
| stripped states $0 < k < K$ | linear in $k$ between the two |
| MA, MX$_q$ monomers | DFT cubes |
| site products (dot − one unit) | the $k = 1$ value for every site |

$$
X_k = \left(1 - \frac{k}{K}\right) X_0 + \frac{k}{K}\, X_K, \qquad X \in \{E_\infty, S\}
$$

## Checks

On the InAs reaction example (`qdw reaction`, 285-atom dots and small molecules):

- The isosurface area agrees with the SCCS quantum surface within 1–4 % for every species, dots
  included.
- At $\varepsilon = 2$ the COSMO electrostatics agree with CP2K SCCS for the small molecules (HCl
  −0.039 vs −0.038 eV, CH$_3$NH$_2$ −0.048 vs −0.042 eV).
- Two independent surface meshes agree within 2 %.
- For one dot the SCCS value violates the variational bound
  $\Delta G^\mathrm{SCCS}_\mathrm{el} \le \Delta G^\mathrm{frozen}_\mathrm{el}$ (SCCS relaxes the density,
  so it can only be more negative), which marks the SCCS number, not COSMO, as unreliable there.

## Limits

- No electronic polarisation of the solute: COSMO on the frozen gas-phase density gives about
  60–65 % of SCCS at $\varepsilon = 80$ for small polar molecules, but agrees at $\varepsilon \approx 2$,
  the range of the nonpolar solvents used in synthesis (toluene, octadecene).
- The stripped states between $k = 0$ and $K$ are interpolated, not computed.
- One gas-phase geometry per species (no re-optimisation in solvent).

## Implementation

| Piece | Where |
|---|---|
| Cube runs and `solvation.json` | qd-workflows `cp2k/esp.py` (`qdw cp2k esp prepare / status / collect`) |
| Cube reader, isosurface, COSMO matrix, β = 0 fit | qd-workflows `reaction/solvation.py` |
| Interpolation in $k$, cavity term, exports | Orchestr.AI `qd/steps/solvation.py`, `qd/solution.py`, `qd/dashboards.py` |
