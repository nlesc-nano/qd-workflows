# Campaign decisions

Open items need Ivan's answer; the 2-hourly brief lists them under "Needs you". When he answers, move the
item to "Decided" with his words.

## Open

- **2026-10-17 · Move the campaign's data to the NHR project space?** The user scratch is at 237,001 of
  250,000 files (95 %); the campaign will create hundreds of thousands of files. Proposal: campaign directory
  `/scratch/projects/bep00157/qdw_campaign` (12 TiB, 1 M files soft); paths in `nhr.yaml` (campaign) and
  `nhr_campaign.yaml` (code, references, monomer library) change accordingly. Asked by: campaign session.
- **2026-10-17 · Hyperion home at 282,572 of 300,000 files (94 %).** Clean up or move files out of
  `/dipc/iinfante`; the campaign-storage agent can propose a list. Asked by: campaign session.
- **2026-10-17 · Extend the desorption/adsorption ladder to 2,000 atoms?** Estimate: ~760 GPU-h (~115k core-h)
  for the 447 dots of 600–2,000 atoms with 12 steps; 3–5x more if the steps grow with the shell. Then every
  dot has thermodynamics in the webapp. Asked by: campaign session.
- **2026-10-17 · Webapp for dots without a ladder:** show the stability section from stability.json alone
  (today the Thermodynamics tab needs detachment.json). Asked by: campaign session.
- **2026-10-17 · Start date of the bep00157 yearly grant** (budget period assumed 2026-10-01 to 2027-09-30).
  Asked by: campaign session.

## Decided

- 2026-10-17 · Scope: all II–VI, III–V, IV–VI library dots up to 2,000 atoms, CdSe/InP first. (Ivan)
- 2026-10-17 · III–V mostly swap; II–VI vacancy and swap; IV–VI clean. (Ivan)
- 2026-10-17 · DFT labelling of ML frames waits until after the GPU campaign. (Ivan)
- 2026-10-17 · CP2K on whole nodes, all cores as MPI ranks; monomers may use fewer cores. (Ivan)
- 2026-10-17 · Monomers once per formula, shared by every dot, in stage 1. (Ivan)
- 2026-10-17 · CP2K cells on one common grid, 6 Å vacuum per side for monomers too. (Ivan)
