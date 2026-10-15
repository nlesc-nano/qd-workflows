"""
Reaction free energies on demand: a balanced reaction sum_i nu_i X_i = 0 whose species
come with DFT-relaxed structures and DFT energies.

    G_i(T, eps, c_i) = E_i^DFT + G_i^corr(T) + dG_solv,i(eps) + kT ln(c_i / 1 M)

G^corr: ideal gas at 1 M (translation, rigid rotor, quasi-RRHO vibrations) from the MACE-MH-1
Hessian at the MACE minimum reached from the DFT structure (mace.py, thermo.py).
dG_solv: SCCS electrostatics (CP2K GPW, a few eps, interpolated in 1 - 1/eps) or Generalized
Born on CP2K PBE Hirshfeld charges, plus the cavity/dispersion term of Hille et al. (2019),
which depends on eps only (solvation.py). Formulas of the reaction: model.py.

    qdw reaction check    reaction.yaml                  balance, units, what will run
    qdw reaction prepare  reaction.yaml --cluster nhr --out runs/x
    (cd runs/x && bash submit.sh)                         MACE (GPU), CP2K (CPU), then analyse
    qdw reaction analyse  runs/x                          results/ (json, plots, reaction.html)
"""
