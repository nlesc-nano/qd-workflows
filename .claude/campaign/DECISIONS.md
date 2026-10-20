# Campaign decisions

Open items need Ivan's answer; the 2-hourly brief lists them under "Needs you". When he answers, move the
item to "Decided" with his words.

**Standing rule (2026-10-18, Ivan):** "After you do stage 1 , run stage 2 on gpu for mace and stage 2 bis on cpu for the dft
runs . I will tell you when I am going to sleep but until tomorrow you have to take decisions . You are asking too much to
me . Many things you can solve by yourself . I just want to be informed". The supervisor decides batches, retries, sizing
and order within the plan and the yearly allocation (no daily cap since 2026-10-18 evening),
records each decision below as "(supervisor)", and informs Ivan. Deleting data still needs Ivan.

## Open

## Decided

- 2026-10-20 06:25 (supervisor) · The watchdog replans GPU (stage 2) leftovers itself after a time-out too (logs/gpu_replan.sh:
  dots not done and not in any queued/running job, 12 processes per node, tag s2auto<date>). s2rest02 TIMEOUT: 1 dot replanned
  (9589822). More geo_opt OOMs of 1,400-1,500-atom dots were resubmitted on 8 nodes by the watchdog.
- 2026-10-20 06:00 (supervisor) · Dots of 1,400-1,950 atoms were OOM-killed in geo_opt on 4 CLX nodes (MaxRSS ~3.5 GB per rank,
  ~340 GB per node > 288 GB): dots >= 1,300 atoms run on 8 CLX nodes (768 cores; memory spread, not speed); 37 queued chains
  moved; watchdog and replan use the rule. Alternative left for Ivan: 1 Genoa node (192 cores, 1.5 TB) needs the Genoa build in
  the chain scripts (only 4 Genoa nodes idle now). Also fixed: a broken line in replan_chains.py (would have crashed the next
  automatic replan) and the watchdog's failure-reason text.
- 2026-10-19 23:55 · Cores: "Large memory is only needed to run trexio single point ... All others can run in standard genoa
  nodes with 192 cores or two nodes of clx ... Below 6-800 atoms you can easily use 192 core ... For the largest ones 384 cores
  is super enough" (Ivan). Done: < 800 atoms on 192 cores (2 CLX nodes), >= 800 atoms on 384 cores (4 CLX nodes), only TREXIO
  on a large-memory node (by the chain itself). The 338 queued bundles/chains resubmitted with these sizes; the watchdog and
  the replan use them. Kept on CLX: the campaign profile and every chain.sh use the CLX build (Genoa would need the Genoa
  build's environment).
- 2026-10-19 23:40 (supervisor, after Ivan: "you always restart the relaxation, so the 12 hours limit is not a problem. If you put
  in a long queue that it may enter very late") · All queued jobs wait on the 32-jobs-per-user limit; cpu-clx has ~260 idle
  nodes, cpu-clx:large (32 nodes) is full. Large dots (>= 800 atoms, 212 chains) now run on standard nodes (2-4 nodes, 12 h,
  continuing from the last geometry) up to the MO window; the chain itself submits TREXIO/trim/QDEX on 1 large-memory node
  (2 above 1,500 atoms) - marker file TREXIO_ON_LARGE, qd-workflows 7bc80e1. The 212 queued ones were moved. Relaxation status
  (geo_opt) at 23:30: done/finished chains 243 + 27 past geo_opt; in progress 23; not started 453 (115 < 400 atoms, 130 of
  400-800, 208 >= 800).
- 2026-10-19 23:10 (supervisor) · DFT bundle packing was 3x too optimistic (my 0.5 h + 0.004 h/atom came from the earlier dots,
  whose geo_opt started from a DFT geometry): library chains take ~4-6 h at 250-380 atoms and >12 h on one node at 630-700
  atoms. 132 queued bundles cancelled; the 222 chains not done and not running replanned by the cost law (stage2/
  replan_chains.py): 81 bundles (chcNNN) + 98 own jobs on 1-7 nodes (time per node <= 10 h), ~190k core-h by the cost law.
  The watchdog replans and submits automatically after a bundle time-out (logs/replan_submit.sh, tag chr<date>).
- 2026-10-19 22:30 (supervisor) · Planner bug of mine: groups.txt of the 8/12-group jobs listed only groups 0-3, and the
  exclusion lists were built from it, so 5 dots ended up in two active batches (s2rest_01/02 and s2rest2_01): the process of
  s2rest2_01/gpu10 (duplicates only) stopped, 3 other duplicates run twice (identical results, atomic writes). Fixed: groups.txt
  lists every group and the planner takes busy dots from the configs of all queued/running stage 2 jobs. s2gt600c_01/03 TIMEOUT
  (their dots are in active jobs); the one uncovered dot (InP-In407P336Cl201F12) in s2rest3_01.
- 2026-10-19 22:00 (supervisor) · Webapp (local): 53 more finished dots published (51 library + Cd333, Cd373); 64 dots with MACE
  + DFT results in the index. Ingest fix (QDSpaceWebApp commit): builder records are merged only when their start geometry is
  identical - a vacancy and a swap cut can share the fingerprint (composition + binned radial distances), and the swap page
  was dropped. Index 1,785 -> 1,816 structures; a second record with the same id gets "-v2" (e.g.
  InP-In-In31P16Cl45-reconstructed-v2 is the swap dot). Stage 3 of the 10 re-relaxed dots needs their relax settings (fmax
  0.005, relax_continue, L-BFGS): batch_02r. Watchdog: any first chain failure resubmitted once (reason logged), >= 800 atoms
  back to cpu-clx:large. InP-In863P794Cl207: CP2K SIGABRT at exit on 6 nodes (as In365), resubmitted.
- 2026-10-19 17:40 (supervisor) · NHR closed for an import update ~11:00-17:28 (jobs kept running). After: the 32-job limit was
  filled by DFT bundles, so stage 3 batch_02 (queued since 10:40) and the GPU leftovers waited: the 317 pending bundles/chains
  niced (1000) so GPU and stage 3 jobs get free slots first. s2gt600b_10/11 TIMEOUT, s2gt600c_06 "FAILED" only for empty-group
  folders (planner fixed: loop over existing groups); the 21 dots not done and not in a job: s2rest2_01..02 (9587038-39).
  Totals: stage 2 829/858, chains complete 159, stage 3 10 (+63 queued); budget 219,858 core-h used (4.9 %), last 24 h 174k.
- 2026-10-19 (Ivan, in chat: "Approve conda clean --packages on Hyperion") · Hyperion home /dipc/iinfante at 94 % of its
  300,000 files: `conda clean --packages --yes` run on hyperion-login-01 after a dry run (322 packages, 3.81 GB). pkgs 5.5 GB
  -> 1.7 GB; envs (13 GB) untouched. Home .conda files 274k -> 164k by our count (find); the quota figure was not re-read.

- 2026-10-19 10:40 (supervisor) · Stage 2 bis had nearly stalled (89 chains started for 742 dots ready): NHR allows 500 queued and
  32 running jobs per user, so chains are now bundled: 192 single-node bundles (1-6 chains one after another on a whole node,
  12 h, CdSe/InP first; stage2bis/chbNNN.sh, packing time max(0.6 x cost law, 0.5 h + 0.004 h/atom)) and 161 multi-node chains
  (all >= 800 atoms: cpu-clx:large, --mem=700G for the TREXIO step; 32 large nodes, so these take days). Upper bound 727k core-h.
  Stage 3 batch_02: 63 dots (job 9585763). Stage 2 leftovers: 78 dots in 4 jobs s2rest_01..04. Budget: 161,696 core-h used
  (3.6 %), last 24 h 125,962 (the GPU burst); whole campaign projected well under 1 M core-h of the 4.5 M.
- 2026-10-19 · Ladder: "I am ok with what you did on 3 processes if it makes the calc more efficient . Let's not change the ladder
  yet . Let's finish the calcs first then we fix the ladder" (Ivan). 3 processes per GPU kept for the remaining large dots; the
  ladder code (steps/detachment.py) is not touched until the campaign's calculations are finished; then profile and speed it up
  (proposed session "Profile and speed up the large-dot MACE ladder").
- 2026-10-19 08:00 (supervisor) · The ladder of large dots is single-threaded Python for hours (one core per process, GPUs at 0 %,
  10-21 GB GPU memory each). Not optimised now: any change to steps/detachment.py changes its cache key and would invalidate
  the 725 finished ladders (and stage 3 of those dots) - a code task for after the campaign. Instead the 105 remaining dots of
  600-2,000 atoms run 12 processes per node (3 per GPU): 6 jobs s2gt600c_01..06 (9585090-95), replacing the 8 queued
  s2gt600b jobs.
- 2026-10-19 · Soft modes: "if after the second round of relaxation they stay with soft mode let's keep them" (Ivan); imaginary
  modes in the spectra: left out of IR/Raman, noted on the page (Ivan: "Ok"). Thermochemistry unchanged (|nu|, floor 10 cm-1,
  Grimme quasi-RRHO). Done: Orchestr.AI 3818d5a (vibspec peak lists and report spectra without imaginary modes; summary
  n_imaginary_left_out), QDSpaceWebApp commit (Structure view: spectra from real modes only; the Hessian card says "soft,
  lowest X cm-1; not in the IR/Raman spectra"), checked on CdSe Cd95. The static report pages refresh at each dot's next
  stage 3 run.
- 2026-10-19 07:40 (supervisor) · Re-relax results (21 dots, fmax 0.005, L-BFGS after the thread fix: 1-160 s per relax):
  imaginary modes gone in 13; 8 keep small ones (-11 to -35 cm-1: Cd477, Cd44Se37, In52P31, Cd183Se150, Cd640Se561-clean,
  Cd55, Cd31, Cd95), soft modes rather than loose convergence. Stage 2 so far: 725 of 858 dots done (446 up to 600 atoms, 279
  above); 133 left (12 jobs). Ivan's re-relax rule applied to the 103 further dots with imaginary modes (38 below -20 cm-1,
  worst AlP-P-Al116P99Cl51-clean -193 cm-1): jobs 9585077-79 (reimag2_01..03, 8 groups per node).
- 2026-10-19 · Field runs stalling above EPS_SCF 1e-7 (16+ runs in 9 dots): "Use 10-6 as threshold" (Ivan). Done: field runs use
  EPS_SCF 1e-6 (Orchestr.AI 36f7430, deployed); the 31 unfinished field inputs switched (20 stalled ones restart from their own
  .wfn; old inputs/outputs kept as *.eps1e-7); runs already converged at 1e-7 kept; the 9 held chains resubmitted
  (9585068-76); the watchdog holds cleared.
- 2026-10-19 04:30 (supervisor) · Bug of mine in the chain template: the rerun of a field/esp run restarted from its own .wfn
  had no `|| true`, so a second failure killed the chain silently under `set -e`. Fixed (qd-workflows a90494d), template and
  the 68 existing chain.sh patched; the watchdog holds any failed chain with a non-converged field run.
- 2026-10-19 01:30 (supervisor) · My error: the 42 library chains resubmitted on 2026-10-18 got a forced 3 h wall time (and the
  watchdog 4 h), overriding the cost law (4-8 h for 400-530 atoms): many timed out and were resubmitted. The watchdog now uses
  each chain's own wall time (cost law), at least 3 h.
- 2026-10-18 22:50 (supervisor) · Watchdog now fixes the routine cases itself (NHR script logs/watch_fix_remote.sh, run every
  30 min from the supervisor session): a chain that timed out or whose geo_opt stopped unconverged is resubmitted once (4 h;
  logs/watchdog_retries.txt); a running geo_opt with converged forces under the old criteria is restarted (forces-only);
  every fix is a campaign-log line. It wakes the supervisor for anything else: other failures, a second failure of a chain,
  overloaded GPU nodes, idle GPUs (except the known CPU phases of the 4-group ladder jobs), geo_opt > 250 steps unconverged.
- 2026-10-18 21:50 (supervisor) · Watchdog findings: (1) the ladder of large dots is CPU-bound in phases (Python site search,
  GPUs at 0 % for long stretches; 7-12 GB GPU memory per process): the 342 queued dots of 600-2,000 atoms replanned with 8
  groups per node (2 processes per GPU), 21 jobs s2gt600b_01..21; the 7 running s2gt600 jobs (61 dots) continue. (2)
  CdSe-Cd-Cd201Se128Cl146 geo_opt plateau under the old criteria (RMS gradient 3-5e-5): resubmitted with forces-only (9583695).
- 2026-10-18 · DFT geo_opt plateaus: "when dft plateu and you think has converged you can accept the structure" (Ivan). Done:
  geo_opt converges on the forces only (MAX_FORCE 4.5e-4, RMS_FORCE 3.0e-4 Ha/bohr kept; MAX_DR, RMS_DR 1.0 bohr, i.e. off):
  new chains (qd-workflows inputs.geo_motion) and the 55 unfinished geo_opt inputs (old kept as cp2k.inp.before-forces-only);
  a resubmitted chain continues from its last geometry with it. InP-P-In40P23Cl51-reconstructed (forces converged for 20
  steps, still stepping) resubmitted (9583602); InP-In149P108Cl123 picks it up when its queued job starts.
- 2026-10-18 20:40 (supervisor) · Efficiency problem found after Ivan asked why a re-relax took hours: the stage 2 job.sh (copied
  from library_01) set OMP_NUM_THREADS=144 for each of the 4 GPU processes of a node: load 61-341 on 144 cores, GPUs 0-76 %;
  BFGS (used up to 500 atoms) at 47 s per step for Cd181 (131 steps, no oscillation). Fixed: 8 threads per GPU process
  (OMP/MKL/OpenBLAS) in every new stage 2 job; the re-relax of the small dots with L-BFGS for every size, resubmitted (9583568);
  the 44 queued s2gt600 jobs resubmitted with the fixed script; the 8 running s2le600 jobs finish as they are. A watchdog in
  the supervisor session checks every 30 min (idle GPUs, overloaded nodes, failed/timed-out jobs, geo_opt > 150 steps) and
  wakes the supervisor to fix what it finds. Ivan: "what are the agents doing ? they should check this stuff and weird
  behavior and fix it".
- 2026-10-18 20:15 (supervisor) · Ladder calibration above 600 atoms (re-relax job 9582815): 1,420-5,860 GPU-s per dot for
  1,087-1,715 atoms (0.4-1.6 GPU-h), below the ~1.7 GPU-h/dot of the estimate Ivan approved; so no need to come back to him.
  Stage 2 with the ladder for all 403 remaining dots of 600-2,000 atoms (library + earlier): 50 whole-node jobs s2gt600_01..50
  (9583457-9583506), 12 h caps (partition max 24 h), upper bound ~190k core-h. Re-relax with fmax 0.005 removed the imaginary
  modes of Cd815, Cd640Se477, Cd605, Cd534, Cd181; Cd477 (-34.8 cm-1) and Cd640Se561-clean (8 modes, -13.9) keep them; the
  re-relax of the smaller dots is slow (Cd181, 393 atoms: 1.7 h of BFGS).
- 2026-10-18 19:30 (supervisor) · QDEX reruns with fuzzy: the 10 published earlier dots done in 1-3 min each (2 fuzzy dashboards
  each); republished in the local webapp (dft/ re-imported, index rebuilt, still 1,785 structures); Electronic structure tab
  checked (fuzzy bands with bulk-band overlay). Cd333, Cd373, In149P104_swap (and In365, In204P153 when their jobs end) follow.
- 2026-10-18 19:00 (supervisor) · NHR QOS normal allows 32 running jobs per user (the bottleneck now, not nodes; scontrol top is
  not allowed): QDEX reruns kept at normal priority, pending library chains niced (Nice=500) so the reruns start first.
  In365 completed on the large node (TREXIO OK) but before the CIF fix: QDEX rerun 9583360. In149P108 geo_opt still converging
  after 5.4 h: continues (9583361). 6 library chains that hit the old 1 h limit resubmitted with 3 h. To consider: several
  small-dot chains per job (each on its own whole node, or one after another) to use fewer of the 32 job slots.
- 2026-10-18 18:30 (supervisor) · No fuzzy dashboards for the earlier dots: each chain copies record.json when it is created,
  and the 23 earlier chains were created before the legacy records got origin.cif, so QDEX ran with run_fuzzy false. Fixed:
  CIF added to the 23 chain records (running chains pick it up at their QDEX step); the 13 finished chains rerun only QDEX
  (old qdex/ and qdex_trap/ kept as *.nofuzzy-20261018, .done as .done.nofuzzy-20261018): jobs 9583286-98, 2 h. The webapp
  pages are republished when they finish. Library chains carry the builder CIF (not affected).
- 2026-10-18 18:00 (supervisor, no daily cap any more) · Released: the re-relax jobs 9582814/15 start now. Stage 2 for all 423
  library dots up to 600 atoms (21 materials; CdSe/InP first; re-relax dots excluded): 10 whole-node jobs s2le600_01..10
  (9583164-73), fresh configs with the ladder limit 2000, predicted ~4.5 h each (~27k core-h; time model 1.3 x 7.1 n^0.93
  GPU-s per dot from today's runs). Stage 2 bis: chains of the 66 library dots whose ladder is done (new chain.sh with the
  own-.wfn restart), 66 whole-node jobs. Ladders above 600 atoms wait for the calibration in 9582815.
- 2026-10-18 · Budget: "Dont worry about that. Do not self-impose anything, just control it" (Ivan). The 20,000 core-h/day
  guardrail and the 2,000 core-h per-submission gate are removed (SUPERVISOR.md, RULES.md); the supervisor controls the used
  total and the projected year-end use against the 4.5 M core-h allocation and tells Ivan if it goes above.
- 2026-10-18 17:30 (supervisor) · Chains: Cd333 completed (its field fxzpp converged with the outer loop before the policy
  change) and InP-In149P104Cl135_swap completed. InP-In365P288Cl231 failed again: TREXIO out of memory (884 atoms), as Cd373;
  truncated orbitals moved aside, resubmitted on cpu-clx:large, 2 nodes, --mem=700G (740G is not accepted): job 9583159.
  InP-P-In204P149Cl165_swap timed out at geo_opt step 181 (converging), resubmitted once 6 h: job 9583158. InP-P-In204P153Cl153
  second retry running (9582792), InP-In149P108Cl123 retry running (9582615). Rule for the library chains: dots of ~800 atoms
  and more run their TREXIO step on cpu-clx:large (--mem=700G): to build into the chain generator before those chains.
- 2026-10-18 · Webapp: "yes, in the webapp locally" (Ivan). The 10 finished earlier dots are in the local QDSpaceWebApp
  (not committed, not pushed): the 6 CdSe in public/II-VI/CdSe/builder/<id>/ (props/ + dft/ next to the builder record; their
  MACE relax started from the earlier DFT geometry), the 4 InP in a new series public/III-V/InP/builder_legacy/<id>/ (legacy
  record + fingerprint); props/solvation_cubes.json = the chain's solvation.json; dft/ by import_dft_results.py;
  library_index.json rebuilt by ingest_library.py (1,781 -> 1,785 structures: 4 new, 6 updated, nothing else). Pages checked in
  the dev server (CdSe Cd156, InP P-In52: MACE + DFT PBE + QDEX, excitons, stability with DFT COSMO). Cd55 and Cd95 will be
  republished after their re-relax.
- 2026-10-18 · Imaginary modes: "Those with imaginary freqs I would try to relax again from the last relaxed point and tighten a
  bit more the force threshold" (Ivan). Done: Orchestr.AI setting relax_continue (relax from props/relaxed.xyz; in run.py, so
  steps/relax.py and every other dot's cache are unchanged), deployed. The 21 dots with imaginary modes (of 125 with a Hessian,
  -10 to -36 cm-1; list $C/stage2/imaginary_dots.txt) re-relax with fmax 0.005 eV/Å (was 0.01), then Hessian and stage 2 again:
  jobs 9582814 (15 dots <= 600 atoms) and 9582815 (6 dots of 1,087-1,715 atoms, which also calibrate the ladder above 600
  atoms), whole nodes, 6 h caps, held by Slurm until 2026-10-19 00:05 (today's budget hold).
- 2026-10-18 15:00 (supervisor) · Stage 3 batch_01 (job 9582805, 9 min, 0 failed): the 10 earlier dots CdSe Cd156, Cd139, Cd95,
  Cd92, Cd55, Cd60 and InP In31, In31_bare, P-In40, P-In52 are complete (report.json; solvation COSMO on the DFT potential).
  With the pilot: 11 finished dots. Webapp copy waits for Ivan (SUPERVISOR: copies only as Ivan approves). Note: Cd55 and Cd95
  are in the re-relax list and their stage 2/3 will be redone.
- 2026-10-18 14:45 (supervisor) · Stage 2 so far: library_02 (job 9582506) TIMEOUT at 4 h, 26 of 68 dots done, 0 failed (rest
  resumes from caches); library_03 (9582507) 32 of 32 done but without the ladder (built with the old 600 limit); ladder rerun
  of the 10 earlier dots ready for stage 3: job 9582803 (running, ~0.6 h). Ladder cost measured: 1,300-3,700 GPU-s per dot for
  CdSe of 500-560 atoms (~1 h), far above the 390 s of the small dots; the 600-2,000-atom ladder is not measured yet and may
  exceed the 2026-10-17 estimate (~115k core-h). Decided: a calibration batch of 4 dots across 600-2,000 atoms (one per GPU)
  first, tomorrow inside the 8,000 core-h stage 2 budget; if the measured cost is above ~2x the estimate, report to Ivan (as he
  decided). No new GPU submissions today: the rolling 24 h use is above 20,000 core-h.
  Agents were again blocked by the auto-mode check (a read-only squeue/tail of the GPU agent) although Ivan's allow rule for
  `ssh -x nhr` is in ~/.claude/settings.json; the supervisor submits and checks on NHR itself from now on, agents analyse only.
  Structure flags (for Ivan, no action): ligands detached in 8 library dots (CdSe-Se-Cd256Se213Cl86, -Cd247Se209Cl76,
  -Se-Cd615Se545Cl140, -Cd813Se730Cl166, -Cd-Cd815Se730Cl170, -Cd605Se535Cl140, -Se-Cd534Se472Cl124, -Cd-Cd533Se471Cl124);
  many bonds broken (>= 10) in 24; Hessian below -20 cm-1 in CdSe-Cd-Cd477Se344Cl266 (-36.0) and CdSe-Cd-Cd815Se730Cl170 (-28.6).
- 2026-10-18 · SCF that stops at MAX_SCF 200 (field and ESP runs): "dont use outer OT as it is not very good. If it hits 200
  cycles just restart it from that .wfn and then it will converge in few cycles" and "Do not change preconditioner, no looser
  and tighter thresholds, etc. Just as I told you . After hitting 200 cycles, the you restart it." (Ivan). Done: outer loop
  removed; the field runs (Orchestr.AI 7db1fa8) and the dot ESP runs (qd-workflows 86f7c9c) write their own .wfn every 50
  cycles; the chain restarts a run that ended "SCF run NOT converged" once from it (`qdw cp2k scf-restart`: SCF_GUESS restart
  + its own WFN_RESTART_FILE_NAME, nothing else). The .wfn are kept: "I would keep the .wfn and remove the wfn backups only.
  We check how much disk is taken afterwards and decide what to do" (Ivan); no backups are written (BACKUP_COPIES 0),
  qd-workflows 50f2f21. Monomer inputs unchanged. The 3 ESP
  runs that had already converged with the outer loop (Cd55 dot_12, Cd95 dot_add12, In31 dot_12: 106, 115, 87 steps) are
  kept; Cd333 field/fxzpp is running with it. Chains generated from now on use the new chain.sh.
- 2026-10-18 ~14:00 (supervisor) · DFT problems decided:
  - Cd373 (815 atoms): out of memory in the TREXIO single point (288 GB/node). Truncated trexio/orbitals.h5 moved aside
    (orbitals.h5.truncated-oom-9582423), resubmit on cpu-clx:large --mem=740G, 4 h. Large dots will need the large-memory
    partition for TREXIO: to build into the chain generator (memory by atom count) before the library chains above ~800 atoms.
  - In365: heap corruption at CP2K exit after a converged geo_opt: one identical resubmit.
  - InP-P-In204P153Cl153: a second resubmit allowed if it times out while geo_opt still converges.
  - SCF stalls at MAX_SCF 200 (Cd333 field fxzpp; ESP of the stripped/added dots of Cd55, Cd95, In31): OT outer loop
    (inner 50 x outer 10, same EPS_SCF) added to the field runs (Orchestr.AI e41634e) and to the stripped/added-dot ESP runs
    (qd-workflows 6cb662c); monomer inputs unchanged (the library still links). The 4 run inputs rewritten, chains resubmitted.
  - Stage 3 batch_01 failed "detachment is not cached": the ladder-limit change (600 -> 2000) changes the detachment key for
    every dot. Detachment reruns on the GPU (identical results) for the 10 dots ready for stage 3 first, then the others.
- 2026-10-18 13:20 (supervisor) · The ladder had been skipped above 600 atoms (Orchestr.AI default detach_max_atoms 600; the
  campaign never raised it), against Ivan's decision of 2026-10-17. Fixed: qd-workflows 0f99e8d (profile props.settings go into
  every job; nhr_campaign sets detach_max_atoms 2000), tested, deployed. Dots above 600 atoms already through stage 2
  (library_02/03, 9 earlier dots) get the ladder in follow-up runs (other steps stay cached). Budget: last 24 h 22,621 core-h,
  above the 20,000 guardrail (stage 1 + chain timeouts + stage 2 together); held back for the rest of today: stage 2 within its
  5,000, DFT/CPU <= 4,000 more. DFT failures (2 TIMEOUT at 4 h, 2 FAILED with ESP/field aborts) and stage 3 for the 10 dots in
  phase cpu delegated to the DFT agent.
- 2026-10-18 (supervisor) · Stage 2 validated: batch1 (job 9582434) 23 earlier dots, 0 failed, 58 min; library_01 (job 9582435)
  40 CdSe/InP 35-217 atoms, 0 failed, 1 h 59 min (~680 GPU-s per dot: detachment 390, md 210, detach_ensemble 70; ~30 core-h per
  small dot). Library stage 2 continues through the GPU agent, CdSe/InP first, small first; budget split changed: stage 2 <=
  5,000 core-h for the rest of 2026-10-18, <= 8,000 core-h/day from 2026-10-19 (stage 1 is finished); DFT keeps 12,000.
  DFT chains: 5 of 17 completed, 5 hit the 1 h limit (2 still in geo_opt) and were resubmitted once with 4 h (jobs
  9582488-92, <= 1,920 core-h). Campaign log $C/logs/campaign.log for every new job (RULES.md).
- 2026-10-18 (supervisor) · Stage 2 batch 1, second retry (job 9582407) failed again for the same 14 dots: the structure step
  had been cached before the CIF was added (its cache key ignored the CIF), so bulk_bond was None and the detachment step
  failed on float(None). Fixed: Orchestr.AI 345cfcc (the bulk CIF enters the cache key of structure, vibspec, stability,
  detachment, solvation, sites, report; relax unaffected, so stage 1 stays cached), tested and deployed. Resubmitted batch 1:
  job 9582434 (all 23 redo structure and the following steps). First library batch of stage 2: job 9582435, the 40 smallest
  CdSe/InP dots (35-217 atoms), whole node, 3 h cap (<= 1,800 core-h). Today so far: stage 1 ~1,930, stage 2 ~750 + up to
  5,400 (two jobs), DFT ~8,350 (17 chains): under 20,000.
- 2026-10-18 · NHR submissions: "Submit to nhr and also grant permission to the agents to submit too" (Ivan). Submitted by the
  supervisor: stage 2 batch 1 retry, job 9582407 (14 dots left, with the reference lock fix); the 17 smaller DFT chains, jobs
  9582408-9582424 (~8,350 core-h estimated). The supervisor cannot change permission settings itself (blocked as
  self-modification): Ivan adds the allow rule for the agents. Next: stage 2 of the 835 library dots once 9582407 shows the
  lock fix works; the six large chains on the following days within the DFT budget.
- 2026-10-18 (supervisor) · Stage 1 library relax moved atoms by more than 2 Å in 577 of 835 dots (mostly Cl: 7,505 atoms; Hg
  1,240; largest 7.5 Å in InP-In-In165P139Cl78-clean; relaxation energies -56 to -267 meV/atom in batch 1). Expected for
  unrelaxed builder geometries; not held. Stage 2's structure step (ligands detached, bonds broken) decides per dot; flagged
  dots will be listed for Ivan. Snapshots: $C/stage1/status_before_relax.json, status_after_relax.json.
- 2026-10-18 (supervisor) · Stage 1 complete: all 835 library dots relaxed and converged (fmax <= 0.05 eV/Å), 6 whole-node
  batches (jobs 9582129, 9582147-51, ~3.2 node-hours, ~1,900 core-h). With the 23 earlier dots: 858 relaxed of 862 in the ledger.
- 2026-10-18 (supervisor) · Stage 2 batch 1 retry (job 9582140): 9 of 23 dots finished every GPU step (the 9 largest); 14 failed
  in stability/detachment because the 4 GPU runs computed the same CdSe/InP bulk references at once (no lock, non-atomic
  write). Fixed: Orchestr.AI e88847d (file lock + atomic write, tested; flock works on the project space), deployed.
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
  (/data/iinfante does not exist). Done 2026-10-18 00:56: copy ended cleanly (1.45 TB), all 30 paths matched in
  file count and bytes, the 30 originals removed (Ivan approved the proposal by name). /scratch/iinfante: 1.38 TiB -> 58 GiB. Still open: home /dipc/iinfante at 94 % of files
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
