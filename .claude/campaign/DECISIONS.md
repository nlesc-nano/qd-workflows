# Campaign decisions

Open items need Ivan's answer; the 2-hourly brief lists them under "Needs you". When he answers, move the
item to "Decided" with his words.

## Open

- **2026-10-17 · Webapp for dots without a ladder:** show the stability section from stability.json alone
  (today the Thermodynamics tab needs detachment.json). Asked by: campaign session.
- **2026-10-17 · Start date of the bep00157 yearly grant** (budget period assumed 2026-10-01 to 2027-09-30).
  Asked by: campaign session.

## Decided

- 2026-10-17 · AlSb monomer (Hirshfeld run stops at 200 SCF steps after 3 identical runs; its ESP run converged, in exactly 200):
  to change the monomer input for all monomers (CG + 2PNT): "No." (Ivan), who asked whether Hirshfeld is needed at all.
  Parked, no input changed. What uses it: the field step's zero-field monomer run (a dot's chain links it from the
  library or runs it itself); the MACE side and the ESP/RESP charges do not. Matters only for AlSb dots in stage 2 bis.
  Open question for Ivan: drop the Hirshfeld monomer runs from the pipeline?
- 2026-10-17 · The 11 earlier InP dots whose DFT geometry carries extra organic atoms (C, N, H in 1:1:5, methylamine): "Yes strip"
  (Ivan). `qdw campaign legacy-records --strip-extra` (qd-workflows 9c445af) removed them, re-checked the In/P/Cl
  composition against each id (all 11 neutral, 92 to 1,012 atoms), recorded the stripped atoms in record.json
  (`origin.stripped_atoms`). Stage 1 relax of the 11: whole node, 4 groups, job 9582087 (<= 600 core-h).
- 2026-10-17 · Stage 1 of the earlier dots: "For the legacy records I am ok with code change . Legacy structures should
  always be re relaxed with the correct input etc . If there is no geometry forget about them and rerun them from
  scratch . I am ok with whole node batch ." (Ivan). Done: `qdw campaign legacy-records` (qd-workflows 6f3707a, tested)
  wrote 12 records (11 CdSe, 1 InP) in /scratch/projects/bep00157/qdw_campaign/stage1/legacy_records, listed in the
  ledger. Stage 1 relax on one whole gpu-a100 node (4 groups, one A100 each, <= 600 core-h): job 9582085. Without a usable
  geometry (rebuilt from the QD_builder library, not yet added to the ledger): CdSe Cd177Se140Cl174 (file has 391 atoms,
  id says 491), Cd301Se224Cl154, Cd768Se632Cl272, Se-Cd508Se421Cl174, Se-Cd600Se489Cl222. Not in scope: the 2 CdSe dots above
  2,000 atoms.
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
