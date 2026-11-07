# Decisions and their reasons

Ivan's rules, with the reason he gave and the date (2026). Dates are 2026, month-day.

| topic | rule | reason / Ivan's words | date |
|---|---|---|---|
| grid | Detach CdSe backwards from the smallest Wulff dot; lattice + relaxed descent, then annealed MD; no corpus structures | generic, applies to other systems | 10-23/24 |
| grid | k range from the builder's moves and the Se deficit, never a fixed k_max | "You do not have to reimplement anything … already available in the QD_builder" | 10-23 |
| grid | One workflow, no mixing with older grids | "once we go to a new system, we cannot run all stuff we have done so far. It should be consistent" | 10-25 |
| grid | Acetate: anneal every node up to m 13, then pool the two runs' lowest states; never touch the grid again | anneal noise ±0.3–0.7 eV per node | 10-28/30, 11-01 |
| grid → growth | Hand-over at m = 13; its k comes from μ | (13, 3) is the builder's first Wulff dot | 10-29, 11-02 |
| builder | Build with the builder, then relax; never graft onto relaxed dots; native Cl untouched | "The best model is always the one that is built with the builder and then relaxed" | 10-25/28 |
| builder | Site roles from the builder's manifest, carried by atom; one site per distinct environment | "why every time you have to reinvent a new way to classify" | 10-29 |
| builder | Facet areas through the outermost lattice sites | area is geometric | 10-30 |
| network | Grow by full layers per direction, charge-balanced by the builder | – | 10-30 |
| network | The growth model uses only the network map, not NEB/γ/size-law parameters | "a completely new growth model … Only on the new map data" | 10-31 |
| network | Islands built by the builder (plane out one layer, vicinal planes cut back) | 2D nucleation on facets | 10-31 |
| network | No anneal for m > 13; at most one short MD | "already a crystal"; "otherwise it becomes a mess" | 11-01/02 |
| network | Acetate = CdCl2 builder dots dressed by the builder's X-type exchange | no need to rerun the builder | 11-01 |
| network | Acetate only from 2026-11-04 | – | 11-04 |
| theory | MACE-MH-1 omat_pbe; D3 for carboxylates; LBFGS; no FF, metadynamics or xTB | – | 10-25/27, 11-01 |
| thermo | Free energies with Hessians; σ·S for dots, COSMO for free species | vary T and concentrations | 10-23/25 |
| thermo | Learned F_vib for relax-only states; keep FD where it exists; no extra FD checks | step error 0.06 eV ≪ 0.27 eV relaxation noise; "Waste of computer time" | 11-01/03/04 |
| kinetics | Packages react like bare CdSe on the total C_M | "to keep it clean we study the reactivity with only cdse" | 10-30 |
| kinetics | Computed values only, never the webapp placeholders | "Lets not mix stuff" | 10-29 |
| model | Qualitative and tunable (activities, concentrations, k_conv) to experiments | "what we do is qualitative and could be fine tuned later" | 10-25 |
| model | CdCl2 is an architecture test only; never push it to grow | "If you force it to grow you are probably breaking some physics" | 10-29 … 11-02 |
| model | Runs must stop cleanly; speed-ups must not change the physics | "You cannot wait forever" | 10-30, 11-03 |
| compute | NHR account bep00157, pack 4 GPU tasks per A100 node; prefer LUMI; nothing heavy on the laptop; no calculation without a usable output | budget | 10-13 … 10-29 |

Superseded along the way: the continuous size law with Δμ per builder role (plain Cd-weighted mean), K_X from site
states and NEB k_dir/k_disp (10-29/30), and the facet-area law (10-30), all on `feat/growth-unified`.
