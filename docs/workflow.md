# The qd-workflows pipeline

One batch of library records goes through three Slurm stages: MACE-MH-1 on a GPU, one CP2K
chain per dot on CPU nodes, then the CPU stage that assembles the properties and the report.
Every step is cached: a job that stops continues where it stopped when resubmitted.

```mermaid
flowchart TB
    QB["QD_builder<br/>record.json + start.xyz"] -->|"qdw props prepare"| GPU

    subgraph GPU["1 · GPU stage: MACE-MH-1"]
        direction LR
        RELAX["relax"] --> STRUCT["structure<br/>core / surface / ligand"] --> HESS["hessian<br/>analytic ≤ 2,000 atoms"]
        HESS --> STAB["stability<br/>monomers.xyz"]
        HESS --> DET["detachment ≤ 600 atoms<br/>detach_k.xyz"]
        HESS --> WIG["wigner (+ md)"]
    end

    GPU -->|"one chain per dot ≤ 2,000 atoms"| CHAIN

    subgraph CHAIN["2 · CP2K chain per dot: PBE, GAPW"]
        direction LR
        GEO["geo_opt"] --> FIELD["field<br/>13 finite-field single points<br/>+ monomers"] --> ESP["esp<br/>cubes: intact + stripped dot,<br/>monomers → solvation.json"]
        ESP --> PDOS["pdos → MO window"] --> TRX["TREXIO → trim"] --> QDEX["QDEX<br/>QP gap, excitons, SOC,<br/>fuzzy bands, PDOS, COOP"]
    end

    CHAIN -->|"after every chain of the batch"| CPU

    subgraph CPU["3 · CPU stage: assembly"]
        direction LR
        F2["field<br/>Born charges, Raman tensors"] --> VIB["vibspec<br/>IR / Raman"]
        F2 --> SOLV["solvation<br/>COSMO + cavity"] --> SITES["sites<br/>binding map"]
        VIB --> REP["report<br/>properties.json + dashboards"]
        SITES --> REP
    end

    CPU --> WEB["QDSpace webapp"]
    CHAIN -->|"QDEX: webapp.json, h5, MO cubes"| WEB

    GPU -.->|"Wigner / MD samples →<br/>CP2K sample single points"| FT["fine-tuning set<br/>DFT energies + forces"]
    FT -.->|"fine-tuned MACE model"| GPU

    RX["qdw reaction<br/>on demand, outside the batch"]:::side
    classDef side fill:#f4f4f1,stroke:#8a8984,stroke-dasharray: 4 3
```

## What each engine provides

| Quantity | Source (dots ≤ 2,000 atoms) |
|---|---|
| Structure, frequencies, thermochemistry, desorption ladder, Wigner / MD ensembles | MACE-MH-1 |
| Electronic structure: levels, QP gap, excitons, absorption, PDOS, COOP | CP2K PBE → QDEX |
| IR / Raman intensities | DFT Born charges and Raman tensors (finite field) on the MACE modes |
| Solvation | COSMO on the DFT potential + cavity/dispersion |

xTB is not used in the batch. Orchestr.AI keeps its xTB steps (electronic levels, g-xTB IR/Raman,
GB on xTB charges) for runs without a DFT chain (`vibspec_source gxtb`, `solvation_charges xtb`).

## Size routing

| Atoms | What runs |
|---|---|
| ≤ 600 | everything, including the desorption ladder and its thermochemistry |
| ≤ 2,000 | analytic Hessian, the CP2K chain and QDEX; no desorption ladder above 600 |
| > 2,000 | MACE relaxation and structure analysis only: every step that needs the Hessian (and the report) is skipped |

## Commands

```bash
qdw props prepare <records> --cluster nhr --out jobs/batch   # GPU job; it submits the chains and the CPU job
(cd jobs/batch && sbatch job.sh)
qdw cp2k chain relaxed.xyz --project <id> --record record.json --out cp2k   # one chain by hand
qdw reaction prepare reaction.yaml --cluster nhr --out runs/x                # on demand
```

Dashed in the diagram: the fine-tuning loop. Its pieces exist (`qdw cp2k prepare --step sample`,
`qdw cp2k label`), but nothing runs it automatically yet. `qdw reaction` stands outside the batch.
