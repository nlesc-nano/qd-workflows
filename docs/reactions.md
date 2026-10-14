# Reaction free energies (`qdw reaction`)

On demand: a balanced reaction whose species come with DFT-relaxed structures and DFT
energies. qd-workflows adds the thermal corrections (MACE-MH-1 Hessians), solvation free
energies (CP2K single points, post-processed) and evaluates the reaction free energy against
temperature, dielectric constant, the electronic-energy level and the chemical potential of
every species, on an interactive page.

```bash
qdw reaction check   reaction.yaml
qdw reaction prepare reaction.yaml --cluster nhr --out runs/x --test-queue
(cd runs/x && bash submit.sh)
qdw reaction analyse runs/x          # results/: results.json, summary.json, plots, reaction.html
```

Examples: `examples/reaction_InAs_HCl/` (InAs dot + 12 HCl) and `examples/reaction_InAs_adduct/`
(the same with the CH3NH2...HCl adduct as the HCl source).

## Input

```yaml
name: InAs285_HCl
energy_method: HLE17                 # what the given energies are
equation: "origin + 12 HCl -> QD_Cl + 4 AsH3 + 8 CH3NH2"
energy_unit: hartree                 # hartree | eV, required
T: [250, 800, 10]                    # K: start, stop, step (298.15 is always added)
per: HCl                             # also per HCl (a species, or a number)
conditions: {eps: 2.38, electrostatics: cosmo, cavity: beta0}
scan: {species: HCl, log10_c: [-8, 1]}
ligand: CH3NH2                       # dG_lig slider per released ligand
species:
  origin: {xyz: origin_285.xyz, energy: xyz, charge: 0, mult: 1, c: 1.0e-6}
  AsH3:   {xyz: AsH3.xyz, energy: xyz, c: 1.0e-6, p_bar: 1.0e-4}
  CH3NH2: {xyz: CH3NH2.xyz, energy: xyz, phase: solvent, c_liquid: 3.04}
  # energy: a number, or xyz (the "E = ..." of the comment line)
  # phase: solute (c, M) | gas (p_bar) | solvent (c_liquid, M); optional symmetry_number, charges
settings:
  sp_functional: HLE17               # GAPW single point at the given geometry (charges, .wfn)
  energy_check_functional: PBE0      # isolated PBE0 (ADMM2) single points: an energy set
  sccs: false                        # SCCS runs (optional check; see below)
```

`qdw reaction check` verifies element and charge balance and suggests the coefficients.

## Model

    mu_i  = E_i + G_i^corr(T) + dG_solv,i(eps) + kT ln a_i
    dG_r  = sum_i nu_i mu_i + nu_L dG_lig = dG_r^0 + kT ln Q

* **E_i**, one of four sets (page selector): the given energies (e.g. HLE17), PBE0@given
  geometries, MACE-MH1 (at its minima), PBE0@MACE-MH1 geometries. Each set is internally
  consistent; a set appears only when every species has it.
* **G^corr**: ideal gas at 1 M (Sackur-Tetrode at 1/(N_A c0), rigid rotor with symmetry number,
  harmonic ZPE and U_vib, quasi-RRHO entropy (Grimme, w0 = 100 cm-1), k ln(2S+1)), from the
  MACE-MH-1 Hessian at the MACE minimum reached from the given structure. Identical to 1 µeV to
  Orchestr.AI `qd.references.ideal_gas_g`.
* **a_i**: c / 1 M (dissolved), c_liquid / 1 M (pure liquid, e.g. neat oleylamine 3.04 M), or
  p / (c0 R T) (a dissolved gas at equilibrium with its headspace: no dG_solv term).
* **dG_lig**: free energy per released ligand that a small surrogate leaves out (chain packing,
  conformations); calibrate on a measured bound/free ratio.

### Solvation

    dG_solv(eps) = (1 - 1/eps) E_inf + (alpha + gamma)(eps) S

* **Electrostatics (default `cosmo`)**: conductor-like screening of the DFT electrostatic
  potential of the gas-phase single point, with the cavity on the electron-density isosurface
  rho_0 = e^-7.6 = 5e-4 e/bohr^3 (the centre of the SCCS switching of Hille et al.). The ESP run
  is the HLE17 single point restarted from its converged .wfn with the potential and density
  cubes printed (one SCF cycle); marching cubes on the density, surface elements on the mesh
  vertices (~0.55 Å). No electronic polarisation of the solute.
  Checks on the InAs example: cavity area within 1-4 % of the SCCS quantum surface for every
  species (dot included); dG_el(eps 2) HCl -0.039 vs SCCS -0.038 eV, CH3NH2 -0.048 vs -0.042,
  QD_Cl -1.33 vs -1.79; two independent surface meshes agree within 2 %. The density level
  matters for ionic surfaces (QD_Cl: -9.2 eV at 1e-3 instead of 5e-4); spheres around atoms are
  worse (4x change between 1.1 and 1.4 x vdW).
* **Cavity + dispersion (default `beta0`)**: Hille et al., J. Chem. Phys. 150, 041710 (2019),
  their beta = 0 fit (SI Table S2, in `src/qd_workflows/reaction/data/`), with
  (alpha + gamma)(eps) = 2.16 ln eps - 14.79 mN/m fitted here over the 67 non-aqueous solvents
  (R^2 0.71), times the isosurface area. Their main fit (beta = -0.5 GPa) has a volume term
  that gives -18 eV for a 5,700 Å^3 dot: molecules only. Defined for eps >= 1.8.
* **SCCS (optional)**: CP2K SCCS single points. Only usable with GPW (no SCCS with GAPW), the
  periodic Poisson solver (MT, ANALYTIC and MULTIPOLE give wrong energies with SCCS in CP2K
  2026.2), one OpenMP thread per rank (with threads the polarisation loop diverges or gives
  wrong values), restarted along an eps ladder from the single point's .wfn, FULL_KINETIC, OT/CG
  for large solutes. Reliable for small molecules; for the InAs dots it converged at eps 2 only.
* **GB**: Generalized Born on RESP charges. CP2K's default Hirshfeld charges (Gaussian shape
  functions) overestimate dipoles 2-30x (AsH3 5.3 D vs 0.19 D) and must not be used for this.

## Recipe: what is put in the flask

`examples/reaction_InAs_toluene/`: washed dots in toluene treated with HCl, oleylamine (CH3NH2 as
stand-in) added or not. The reaction and the side equilibria (`equilibria:`) are solved together
from the amounts put in (`recipe:`): dot concentration, reagent equivalents per dot, additive
concentration (0 allowed; the etched dots release their own amines).

    c_i      = c_i0 + sum_r nu_ri xi_r                (dissolved: mass balance)
    a_i      fixed                                    (gas at its p: open flask; pure liquid)
    0        = dG_r^0 / kT + sum_i nu_ri ln a_i       (every reaction r)
    c_dot0   = c_dot, c_reagent0 = n_eq c_dot, c_additive0 = c_add, all else 0
    theta    = xi_0 / c_dot                           (fraction of the dots converted)

G is convex in the extents, so each residual rises with its extent: bisection on the main extent,
and for each trial the side extent (closed form for one A + B -> C: (c + y) = K (a - y)(b - y)).
The page has sliders for c_dot, log10 equivalents, the additive and a dG correction of the side
equilibrium (ion pairing beyond the computed H-bonded complex). All 12 HCl per dot are one step:
theta is the fraction of fully converted dots.

## On NHR

* Jobs are charged to `bep00157` (cluster profile `account:`; `SBATCH_ACCOUNT` in `~/.bashrc`).
* `--test-queue`: 1 h pieces that resubmit themselves; CP2K stops cleanly 150 s before the limit
  (EXIT file) and the next piece continues from the .wfn; the last piece submits the analysis.
* Run CP2K pure MPI (192 x 1 per node): with OpenMP threads every rank ended up on one core.
  Small molecules with exact exchange: 2 ranks (8-24 ranks were 10-30x slower per SCF step).
* The analysis needs scikit-image (marching cubes) and scipy.
