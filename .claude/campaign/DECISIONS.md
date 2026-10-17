# Campaign decisions

Open items need Ivan's answer; the 2-hourly brief lists them under "Needs you". When he answers, move the
item to "Decided" with his words.

**Standing rule (2026-10-18, Ivan):** "After you do stage 1 , run stage 2 on gpu for mace and stage 2 bis on cpu for the dft
runs . I will tell you when I am going to sleep but until tomorrow you have to take decisions . You are asking too much to
me . Many things you can solve by yourself . I just want to be informed". The supervisor decides batches, retries, sizing
and order within the plan and the 20,000 core-h/day guardrail (split: stage 1 4,000, stage 2 4,000, stage 2 bis 12,000),
records each decision below as "(supervisor)", and informs Ivan. Deleting data still needs Ivan.

## Open

- **2026-10-18 · Agents blocked by the auto-mode permission check.** The DFT agent's `sbatch` on NHR (and then a read-only
  ssh) was denied ("Shared Cluster Mutation"); nothing submitted. The supervisor does not submit in its place without Ivan.
  Ivan: allow the submissions (main session or a permission rule for the agents)?
- **2026-10-17 · Hyperion home /dipc/iinfante at 94 % of its 300,000 files** (282,572; 97 % in `.conda`: 170,530 in envs,
  104,047 in pkgs; the CsPbX3 move does not touch home). Recommendation: `conda clean --packages` (pkgs, 5.5 GB, rebuilt on
  demand); envs stay. Needs Ivan (deletion). Not blocking the campaign (it runs on NHR).

## Decided

- 2026-10-18 (supervisor) · Stage 2 batch 1 (job 9582120, 163 core-h): structure and Hessian done for the 23 earlier dots
  (6 with one small imaginary mode, -5.9 to -13.5 cm-1, noise); stability and everything after failed: the legacy records
  had no bulk CIF. Bug of `legacy-records`, fixed (qd-workflows: origin.cif = <material>_zb.cif, tested), the 23 records
  patched, the same batch resubmitted: job 9582140 (whole node, 6 h cap, <= 3,600 core-h, estimate ~3,000).
- 2026-10-18 (supervisor) · 8 InP earlier dots show many bond changes in the MACE structure check (4-9 broken, 8-50 formed:
  In407P336Cl201F12, In407P336Cl213, In415P324Cl273_swap, In149P108Cl123, P-In204P149Cl165_swap, P-In204P153Cl153,
  P-In52P31Cl63, P-In40P23Cl51; mostly the ones whose methylamine was stripped). Not held: they go on through stages 2 and
  2 bis like the others; flagged for Ivan to look at in the webapp.
- 2026-10-18 (supervisor) · Stage 2 bis order and the six chains above 2,000 core-h: approved. Chain scripts exist for all 23
  relaxed earlier dots (<record>/cp2k/chain.sh). Estimates (pilot Cd176, 429 core-h, scaled n_bf^1.72): 17 chains of 96 to
  1,728 core-h, about 8,350 in all, first; then InP-In407P336Cl201F12 2,003, InP-In407P336Cl213 2,003, InP-In415P324Cl273_swap
  2,171, CdSe-Cd446Se352Cl188 2,178; then CdSe-Cd663Se560Cl206 4,171 and CdSe-Se-Cd856Se731Cl250 6,460, over the next days
  within 12,000 core-h/day. CdSe-Cd373Se304Cl138 may hit its 12 h cap: one resubmit. Waiting on the permission item above.
- 2026-10-18 (supervisor) · Spin state of the MA diatomics: kept closed-shell, as every dot and the MACE references; the MA
  monomers now only enter the ESP/COSMO reference. Revisit only if a III-V ladder looks wrong.
- 2026-10-18 (supervisor) · Webapp for dots without a ladder: no longer needed, every dot up to 2,000 atoms gets the ladder.
- 2026-10-18 (supervisor) · bep00157 budget period kept as 2026-10-01 to 2027-09-30 until the grant letter says otherwise.
- 2026-10-18 · Library selection for the ledger: the 904 records of the scope (541 II-VI, 306 III-V, 57 IV-VI, reproduced by
  series) are 847 distinct structures: 57 records are the same structure in two series (54 III-V clean in builder and
  builder_100; 3 II-VI identical in a vacancy and a swap series), 14 ids carry two different structures (vacancy and swap).
  Proposal: relax each structure once, key the 14 swap variants apart; Ivan: "Yes go ahead". Key used: `<id>_swap` (the id
  is also the CP2K project and job name; `_swap` is the convention of the earlier InP dots), record.json id rewritten in the
  copy, `origin.campaign_key` says so. Ledger (backup campaign.yaml.before-library-2026-10-18): 862 dots = 831 library +
  4 earlier CdSe rebuilt from their library record + 23 relaxed earlier dots + pilot + 3 without a start (Cd177Se140Cl174
  not in the library; the 2 earlier CdSe above 2,000 atoms). Duplicates' other series paths in `also_in` (results to be
  copied to both webapp paths). Started the same day, in parallel as the plan says: stage 1 of the 835 new dots, stage 2
  (GPU, no CPU steps) and stage 2 bis (CP2K chains) of the 23 relaxed dots, each through its role agent.
- 2026-10-18 · NHR user-scratch cleanup (proposal `proposals/2026-10-18-nhr-scratch.txt`): Ivan: "can't you run the calculations
  in the project folder that has more space ? Didn't we decide this before already ?" Yes: the campaign runs in
  /scratch/projects/bep00157/qdw_campaign, so the 95 % file warning on /scratch/usr/becivinf does not block it. Proposal
  shelved, nothing approved or deleted. Note: the campaign still READS from /scratch/usr/becivinf (MACE model, CIFs, QDEX
  build, basis data, mostly under qdw_tests): that stays out of any cleanup until moved to the project space.
- 2026-10-18 · Hirshfeld monomer runs: Ivan asked whether the webapp needs Hirshfeld; supervisor checked (the webapp never
  reads or shows it; only a GB fallback used it) and recommended dropping the runs; Ivan: "Yes". Done: Orchestr.AI
  3960395 (field step without monomers; solvation without COSMO data = vacuum, with a warning and a model line; sites
  and dashboard note likewise), qd-workflows e15c4bd (library and chain: ESP only; `monomers link --esp`),
  QDSpaceWebApp e86c9c2d (a "vacuum" model is shown as vacuum). Deployed on NHR (code/COMMITS). AlSb is no longer a problem:
  the library is 28 of 28 (ESP). The 27 old Hirshfeld run directories stay in the library, unused. The ESP step now matters
  more: a dot whose ESP runs fail gets vacuum thermodynamics.
- 2026-10-17 · AlSb monomer (Hirshfeld run stops at 200 SCF steps after 3 identical runs; its ESP run converged, in exactly 200):
  to change the monomer input for all monomers (CG + 2PNT): "No." (Ivan). Parked, no input changed.
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
