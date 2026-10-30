"""Site-kinetics protocol: surface sites of a QD_builder dot, their CdCl2-type binding K_X, the facet chemical
potential Delta mu_f and the NEB rate constants k_dir, k_disp that feed the per-facet size law

    dm/dt = N_A sum_f A_f (c*_M - C_eq e^(Delta mu_f/kT)) / (delta_f/D + 1/k_eff^(f)),
    k_eff^(f) = (k_dir + k_disp K_X,edge c_X) / (1 + K_X,eff c_X).

Steps (qdw sites <step>), one directory per run, every file JSON:
  build      QD_builder dot per size, the builder's own site roles per atom (manifest site_roles), carried by atom
             through the surface treatment (half swap)                                   -> dots/<tag>.json
  reference  relax each dot, no anneal (GPU)                                             -> reference/<tag>.json
  pick       one site per symmetry-distinct environment of each role, the CdCl2 unit and the three states
             ref (dot-X), empty (dot*), M (dot-CdSe)                                     -> sites/<tag>.json
  states     relax every state, partial-Hessian free energies, free CdCl2 (GPU)          -> states/<tag>.json
  table      dG_bind,X, K_X and Delta mu_f per size, role and environment               -> role_table.json
  neb        k_disp and k_dir endpoints from the states, CI-NEB, PHVA at the TS (GPU)    -> neb/<tag>/<site>/<kind>/
  inputs     role table + NEB results in one file, and marked xyz of every NEB endpoint  -> kinetics_inputs.json
  job        a batch script for a GPU step (cluster block of the recipe)
Computed values only: nothing here reads the growth template or the webapp's placeholders."""
