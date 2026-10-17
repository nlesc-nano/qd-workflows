# Campaign decisions

Open items need Ivan's answer; the 2-hourly brief lists them under "Needs you". When he answers, move the
item to "Decided" with his words.

## Open

- **2026-10-17 · Hyperion home files:** moving CsPbX3 frees none of the 282,572 home files (97 % are in `.conda`:
  170,530 in envs, 104,047 in pkgs). Recommendation: `conda clean --packages` (pkgs, 5.5 GB, safe: rebuilt on demand)
  after your OK; envs stay. Asked by: campaign-storage.
- **2026-10-17 · AlSb monomer does not converge (both runs of job 9582047; 27 of 28 monomers done).** OT DIIS,
  EPS_SCF 1e-6, MAX_SCF 200: energy flat to ~1e-10 Ha (-7.42238626 Ha) after step 50, gradient 1.1-2.9e-6, just above
  EPS_SCF; likely near-degenerate HOMO/LUMO. Options: OT MINIMIZER CG + LINESEARCH 2PNT, more MAX_SCF, or OUTER_SCF,
  for the monomer inputs only. Any change goes into the shared monomer input code (Orchestr.AI `field.monomer_input`,
  qd-workflows esp/monomers), or the chains' monomer runs stop matching the library. Recommendation: CG + 2PNT for all
  monomers (same inputs everywhere), rerun AlSb first as a test. Asked by: campaign session (CdSe workflow agents planning).
- **2026-10-17 · Spin state of the MA diatomics.** Group 13-15 diatomics (AlP, AlAs, AlSb, GaP, GaAs, InP, ...) likely
  have triplet ground states, while every monomer run is closed-shell and the MACE references have no spin. Options:
  keep closed-shell (consistent with the dots) or run MA as triplets (UKS, MULTIPLICITY 3). Affects the MA monomer
  reference (Hirshfeld charges, ESP solvation), not the dots. Recommendation: keep closed-shell, and check the energy
  gap singlet vs triplet for one case (InP) before deciding. Asked by: campaign session.
- **2026-10-17 · Webapp for dots without a ladder:** show the stability section from stability.json alone
  (today the Thermodynamics tab needs detachment.json). Asked by: campaign session.
- **2026-10-17 · Start date of the bep00157 yearly grant** (budget period assumed 2026-10-01 to 2027-09-30).
  Asked by: campaign session.

## Decided

- 2026-10-17 · Desorption/adsorption ladder extended to dots up to 2,000 atoms: "The ladder do it up to 2000 atoms" (Ivan).
  Estimate ~760 GPU-h (~115k core-h) for the 447 dots of 600-2,000 atoms at 12 steps. Order as proposed: ladders
  of dots up to 600 atoms first; then a calibration batch of 4-6 CdSe/InP dots across 600-2,000 atoms; then the rest,
  small dots first, in batches under 20,000 core-h a day. Come back to Ivan only if the measured cost is above ~2x the
  estimate. Not started: stage 1 is not finished.
- 2026-10-17 · Hyperion: CsPbX3 data (it is on /scratch/iinfante, not in home; 30 paths, 1.32 TiB, 46,385 files, list in
  `proposals/2026-10-17-hyperion-cspbx3.txt`) goes to the archive: "You should just move the files in /scratch/iinfante
  in the archive in hyperion as well where you have 10 TB available" (Ivan). Archive = /data/ivan.infante/cspbx3
  (/data/iinfante does not exist). Copy running since 23:21 (rsync, log ~/cspbx3_copy.log on Hyperion). The originals
  are removed only after the copy is verified and Ivan confirms. Still open: home /dipc/iinfante at 94 % of files
  (97 % are in .conda envs/pkgs; not touched by this move).
- 2026-10-17 · Move the campaign's data to the NHR project space `/scratch/projects/bep00157/qdw_campaign`:
  "Settle the project-space move first, and go with your recommendation." (Ivan). Plan: copy (not move) code,
  references, monomers and ledger once monomer job 9582047 has ended, then switch `nhr.yaml`/`nhr_campaign.yaml`
  and redeploy; the old copy is deleted only after Ivan approves a list. The pilot (`qdw_tests/cd176_pilot`) stays
  where it is. Status: copy done 2026-10-17 (job 9582047 ended, counts and sizes match); paths switched in the checkout, old copy kept.
- 2026-10-17 · Scope: all II–VI, III–V, IV–VI library dots up to 2,000 atoms, CdSe/InP first. (Ivan)
- 2026-10-17 · III–V mostly swap; II–VI vacancy and swap; IV–VI clean. (Ivan)
- 2026-10-17 · DFT labelling of ML frames waits until after the GPU campaign. (Ivan)
- 2026-10-17 · CP2K on whole nodes, all cores as MPI ranks; monomers may use fewer cores. (Ivan)
- 2026-10-17 · Monomers once per formula, shared by every dot, in stage 1. (Ivan)
- 2026-10-17 · CP2K cells on one common grid, 6 Å vacuum per side for monomers too. (Ivan)
