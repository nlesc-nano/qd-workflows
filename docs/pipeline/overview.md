# Pipeline overview

One batch of library records goes through three Slurm stages. MACE-MH-1 runs on a GPU; each dot
gets one CP2K chain on CPU nodes; then a CPU job assembles the properties and the report. The GPU
job submits the chains and the CPU job itself, so a batch is one `sbatch`.

```mermaid
flowchart TB
    QB["QD_builder<br/>record.json + start.xyz"] -->|"qdw props prepare"| GPU

    subgraph GPU["1 · GPU stage: MACE-MH-1"]
        direction LR
        RELAX["relax"] --> STRUCT["structure<br/>core / surface / ligand"] --> HESS["hessian<br/>frequencies → G°(T), 1 M, qRRHO"]
        HESS --> STAB["stability<br/>bulk + monomer references<br/>ΔG_dec, ΔG_bind, T_diss (gas)"]
        HESS --> DET["detachment ≤ 600 atoms<br/>ladder ΔG_k, ⟨k⟩(T, μ)"]
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
        F2 --> SOLV["solvation<br/>ΔG_solv(ε) per species"] --> SITES["sites<br/>ΔG_site per unit"]
        SOLV --> EQ["solution thermochemistry<br/>μ_i(T, ε, c), ΔG_k, ΔG_dec, ΔG_bind,<br/>coupled equilibria"]
        VIB --> REP["report<br/>properties.json + dashboards"]
        SITES --> REP
        EQ --> REP
    end

    CPU --> WEB["QDSpace webapp"]
    CHAIN -->|"QDEX: webapp.json, h5, MO cubes"| WEB

    GPU -.->|"Wigner / MD samples →<br/>CP2K sample single points"| FT["fine-tuning set<br/>DFT energies + forces"]
    FT -.->|"fine-tuned MACE model"| GPU

    RX["qdw reaction<br/>on demand, outside the batch"]:::side
    classDef side stroke-dasharray: 4 3
```

| Stage | Where | Steps | Page |
|---|---|---|---|
| 1 · GPU stage | one GPU (A100) | relax, structure, hessian, stability, detachment, wigner (md) | [GPU stage](gpu-stage.md) |
| 2 · CP2K chain | CPU nodes, one job per dot | geo_opt, field, esp, pdos, window, TREXIO, trim, QDEX | [CP2K chain](cp2k-chain.md) |
| 3 · CPU stage | one CPU node | field, vibspec, solvation, sites, report | [CPU stage](cpu-stage.md) |

## Thermochemistry and equilibria

What happens to the energies after the Hessian. Every box is a formula in the
[Theory](../theory/thermochemistry.md) pages; the dashboards evaluate the last rows live, with
sliders for $T$, $\varepsilon$, the concentrations and the precursor shift.

```mermaid
flowchart TB
    E["E_i: MACE-MH-1 energy at the minimum"] --> G0
    NU["ν_j: MACE-MH-1 Hessian"] --> G0
    G0["G°_i(T): ideal solute at 1 M<br/>translation + rotation (σ) + vibrations (quasi-RRHO)"]
    BULK["bulk MA: phonons → G_bulk(T)"]
    G0 --> GS["G_sol,i(T, ε) = G°_i + ΔG_solv,i(ε)<br/>COSMO on the DFT potential + cavity"]
    GS --> MU["μ_i = G_sol,i + kT ln c_i<br/>(+ Δμ_prec for the MA monomer)"]

    MU --> R1["decomposition<br/>dot → n MA(bulk) + m MX_q<br/>ΔG_dec per MA"]
    BULK --> R1
    MU --> R2["binding<br/>dot → n MA + m MX_q<br/>ΔG_bind per unit, T_diss"]
    MU --> R3["desorption ladder<br/>dot_k−1 → dot_k + MX_q<br/>ΔG_k(T, ε, c)"]
    MU --> R4["sites<br/>dot → dot − u_i + MX_q<br/>ΔG_site,i"]

    R3 --> K["⟨k⟩(T, μ_MX): Boltzmann over every<br/>configuration the search relaxed"]
    R3 --> CEQ["coupled equilibria<br/>c_i = e^(−ΔG°_i/kT) c_MA^n_i c_MX^m_i<br/>+ mass balance for C_MA, C_MX<br/>(+ bulk precipitation above c_sat)"]
    R2 --> CEQ
    K --> DASH["report and dashboards"]
    CEQ --> DASH
    R1 --> DASH
    R4 --> DASH
```

## What each engine provides

| Quantity | Source (dots ≤ 2,000 atoms) |
|---|---|
| Structure, frequencies, thermochemistry, desorption ladder, Wigner / MD ensembles | MACE-MH-1 |
| Electronic structure: levels, QP gap, excitons, absorption, PDOS, COOP | CP2K PBE → QDEX |
| IR / Raman intensities | DFT Born charges and Raman tensors (finite field) on the MACE modes |
| Solvation | COSMO on the DFT potential + cavity/dispersion |

No xTB: `qdw props prepare` sets `vibspec_source: dft` and `solvation_charges: dft`, so IR/Raman
and solvation come only from the DFT chain. (Orchestr.AI keeps its xTB steps for standalone runs.)

## Size routing

| Atoms | What runs |
|---|---|
| ≤ 600 | everything, including the desorption ladder and its thermochemistry |
| ≤ 2,000 | analytic Hessian, the CP2K chain and QDEX; no desorption ladder above 600 |
| > 2,000 | not in the batch: `qdw props prepare` leaves these records out (to be decided) |

## The fine-tuning loop

Dashed in the diagram. The Wigner and MD ensembles of the GPU stage are labelled with CP2K
single points (`qdw cp2k prepare --step sample`, energies and forces, `qdw cp2k label`), which
make a fine-tuning set for MACE; a fine-tuned model reruns the GPU stage. The pieces exist;
nothing runs the loop automatically yet.
