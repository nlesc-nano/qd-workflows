# Campaign rules (every agent reads this first)

The structure database of quantum dots: the whole QD_builder library of II–VI, III–V and IV–VI dots up to
2,000 atoms (~904 dots, 21 materials) plus 28 earlier CdSe/InP DFT dots. Owner: Ivan Infante.
Code: this checkout (branch `feat/campaign`), deployed on NHR in `/scratch/usr/becivinf/qdw_campaign/code`.

## Stages

1. **MACE-MH-1 relaxation of every dot** (GPU), the 28 earlier dots first (from their DFT geometries), and the
   **monomer library**: MA and MX_q once per formula (`qdw monomers`), shared by every dot.
2. **GPU campaign** (MACE): Hessian, stability, Wigner, MD, desorption/adsorption ladders (dots ≤ 600 atoms
   for now), small dots first.
2 bis. **DFT campaign** (CP2K PBE + QDEX): one chain per dot, small dots first, independent of stage 2.
3. **CPU assembly** per dot once its GPU results and its chain are both finished → QDSpace webapp.
- DFT labels of the ML frames: **not planned**. Ivan decides after the GPU campaign. Never propose or start it.
- Dots above 2,000 atoms: not in the campaign.

## Selection

II–VI: clean, vacancy- and swap-reconstructed. III–V: clean and swap (vacancy only for the earlier InP dots).
IV–VI: clean. The NHR ledger (`/scratch/usr/becivinf/qdw_campaign/campaign.yaml`, `selection`) is the reference.

## Where to look

- State of everything: `qdw campaign poll --json` from the Mac (all clusters), or on a cluster
  `qdw campaign status <dir> --json`. Read its JSON; do not grep CP2K outputs unless a problem points to a file.
- Run logs (`pilot.log`-style): one line per step, ERROR/WARN lines name the file to read.
- Monomer library: `qdw monomers status --cluster nhr_campaign`.

## Cluster rules (NHR)

- `ssh -x nhr` and `scp/rsync -e "ssh -x"` (X11 forwarding hangs).
- Every job on account **bep00157**. Check `sinfo -o "%.18P %.16F"` before choosing a partition; test
  partitions (1 h) share one running job per user, regular partitions often have idle nodes.
- **CP2K: whole nodes, every core an MPI rank, OMP_NUM_THREADS=1.** No side-by-side low-rank runs of dots.
  Only small molecules (monomers) may run side by side with few ranks.
- Memory per node (`--mem`), never `--mem-per-cpu` on Genoa; a standard CLX node gives a job at most 330 GB.
- CP2K cells come from `cell_abc` (common grid, 6 Å vacuum each side): never write cells by hand.
- Profile for campaign jobs: `nhr_campaign` (CLX build, whole nodes, campaign code and references).
- MACE runs on NHR GPUs, never on the Mac.

## Storage

- `qdw campaign status` reports every quota area (space and files) with warnings from 80 % of the limit.
- File counts fill first: NHR user scratch 250,000 files (soft), Hyperion home 300,000. The campaign's growing
  data belongs in the NHR project space `/scratch/projects/bep00157` (12 TiB, 1 M files soft) once Ivan agrees.
- Cleanup goes through the `campaign-storage` agent: proposal file first, deletion only after Ivan approves
  that file by name.

## What an agent may do without asking

- Read anything; run status commands; `squeue`, `sacct`, `sinfo`.
- Resubmit a failed or timed-out step **once**, with identical settings, when the failure is clearly transient
  (node failure, time limit with the step restartable). Say so in the report.

## What needs Ivan's approval (return it to the supervisor as a question)

- Any change of settings, inputs, cutoffs, cells, code or profiles.
- A submission estimated above **2,000 core-h**, or more than **20,000 core-h** in one day for the campaign.
- Cancelling a job; deleting or moving any file; anything in another user's or the pilot's directories.
- git push, merges, pull requests.

## Decisions and notifications

- Anything that needs Ivan goes into `.claude/campaign/DECISIONS.md` (this checkout) as an open item: date,
  question, options, recommendation, who asked. The supervisor also asks him directly in its session and
  sends a push notification starting with "ACTION:".
- The 2-hourly brief (scheduled task "QD campaign brief") lists the open items of DECISIONS.md under
  "Needs you", next to failures, budget and storage warnings.
- When Ivan answers, the item moves to "Decided" with his answer, word for word.

## Reports

Short and factual: what you checked, what you found (dot, job, file), what you did, what needs a decision.
Numbers with units. Never claim a job succeeded without reading its state from `sacct` or the status JSON.
No co-author or "Generated with" lines in commits.
