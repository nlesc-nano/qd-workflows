# Classical nucleation-and-growth parameters for II-VI, III-V and IV-VI colloidal nanocrystals: literature survey

Research only; no code changed. Date: 2026-10-26.

## How the numbers were checked

- Publisher sites (pubs.acs.org, pubs.rsc.org, science.org) were behind bot verification and could not be read. NSF PAR and OSTI refused connections. Numbers therefore come from open copies: PMC / Europe PMC (with supplementary files), arXiv, eScholarship, a Hamburg PhD thesis, and repository abstracts.
- Verification tags used in the tables:
  - **FT**: read in the full text or SI.
  - **Abs**: read in an abstract or repository record.
  - **Snip**: seen only in a search-engine snippet. Treat as unconfirmed.
  - **derived**: computed here from reported numbers. Not a literature value.
- Main gaps:
  - Rempel/Bawendi/Jensen JACS 2009: the parameter table could not be read.
  - Abé/Hens ACS Nano 2012 and 2013, and De Nolf/Hens JACS 2015: the SI parameter sets could not be read.
  - Lazzari/Abolhasani/Jensen RCE 2017 and Mantzaris CES 2005: could not be read.
  - Owen JACS 2010 rate constants and Hendricks Science 2015 per-compound k: could not be read.
  - McMurtry/Owen Chem. Mater. 2020 (InP) full text: could not be read.
  - Wang/Segets Chem. Eng. Sci. 2023 (InP PBE fit of γ, k, Ea): paywalled.
- **Correction to the brief:** doi 10.1021/acsomega.0c05223 is **not by Lazzari**. It is Wang, Xing, Wei, Wang, Wang, Shen, "Modeling of nucleation and growth in the synthesis of PbS colloidal quantum dots under variable temperatures", ACS Omega 2021, 6, 3701 (PMC7876681).

---

## (a) Parameter tables

### 1. Surface energy / interfacial tension γ

There are three kinds of number here, and they should not be mixed:

- **Kinetic-model γ for capped dots in solution.** These are assumed or fitted values of about 0.1 to 0.45 J/m², and exist only for CdSe.
- **Bare-surface values.** DFT gives about 0.2 to 1.6 J/m²; calorimetry gives about 1 to 1.8 J/m².
- **Calorimetric surface enthalpies of oleate-capped dots.** These are **negative** (InP −3.7, ZnS −1.4 J/m²).

| Material | Value | Units | Kind | Conditions / facet | Source | Check |
|---|---|---|---|---|---|---|
| CdSe | 0.125 (scanned 5e-3 to 1) | J/m² | assumed (model input) | MC Ostwald-ripening simulations; wurtzite V_m | Talapin et al., J. Phys. Chem. B 105, 12278 (2001), doi:10.1021/jp012229m; parameter table in D. V. Talapin, PhD thesis, Univ. Hamburg 2002 (Table in ch. 6) | FT (thesis) |
| CdSe | 0.17 (printed "J/cm²"; must be J/m²) | J/m² | fitted, capped (oleic acid / TOP) | TOPSe + Cd oleate in ODE, ~265 °C | Bullen & Mulvaney, Nano Lett. 4, 2303 (2004), doi:10.1021/nl0496724 | Abs |
| CdSe | 0.44 | J/m² | model input (from Peng 1998 conditions) | 573 K, fit to Peng/Wickham/Alivisatos 1998 data | Fanelli, Cregan, Font, Myers, Int. J. Heat Mass Transf. 165, 120643 (2021); arXiv:1901.08990 (Table 1) | FT |
| CdSe (wurtzite) | 1.31 ± 0.26 (hydrated); 1.65 ± 0.27 (anhydrous) | J/m² | measured, calorimetry (bare / water-covered) | 20–60 nm particles | Xu, Zhou, Navrotsky, J. Mater. Res. 26, 720 (2011), doi:10.1557/jmr.2011.20 | Abs |
| CdSe | numbers not found | – | DFT, per facet, bare and passivated | (10-10), (11-20), (0001), (000-1) | Manna et al., J. Phys. Chem. B 109, 6183 (2005); Csik, Russo, Mulvaney, Chem. Phys. Lett. 414, 322 (2005); Rempel, Trout, Bawendi, Jensen, J. Phys. Chem. B 109, 19320 (2005) | Abs only |
| ZnS sphalerite | 1.76 ± 0.32 (anhydrous); 1.25 ± 0.21 (hydrated) | J/m² | measured, calorimetry | – | Subramani, …, Navrotsky, Geochim. Cosmochim. Acta 340, 99 (2023), doi:10.1016/j.gca.2022.11.003; arXiv:2601.01316 | FT |
| ZnS wurtzite | ~1.50 ± 0.50 (anhydrous, estimated); 0.99 ± 0.32 (hydrated) | J/m² | measured, calorimetry | – | same | FT |
| ZnS sphalerite / wurtzite | 0.86 / 0.57 | J/m² | computed, classical MD, bare | 300 K, average | Zhang, Huang, Gilbert, Banfield, J. Phys. Chem. B 107, 13051 (2003), doi:10.1021/jp036108t | Snip |
| ZnS (ZB) | 1.41 (111), 1.42 (-1-1-1), 0.38 (110) | J/m² | DFT-GGA, bare, unreconstructed | Zn-rich limit; 88.2 / 88.6 / 23.9 meV/Å² | Zhang, Zhang, Tse, Chan, Deng, Zhu, Sci. Rep. 6, 20055 (2016); arXiv:1506.05570 (Table I) | FT |
| ZnS | −1.41 ± 0.03 | J/m² | calorimetric surface enthalpy, oleate-capped | 298 K | Calvin et al., PNAS 121, e2307633121 (2024); arXiv:2210.16616 | quoted in a secondary table |
| InP | −3.71 ± 0.18 | J/m² | calorimetric surface enthalpy, oleate-capped | 298 K | Calvin, Brewer, Crook, Kaufman, Alivisatos, PNAS 121, e2307633121 (2024) | FT |
| InP | 0.58 (110), 1.17 (111) | J/m² | DFT, bare (?) | – | original paper not identified (slide deck) | Snip, **unverified** |
| InAs | 0.65 (40.5 meV/Å²) | J/m² | DFT-LDA, bare | relaxed (110) | Leitsmann & Bechstedt, arXiv:cond-mat/0611521 (2006) | FT |
| GaAs | 0.83 (110); 0.72–0.87 (111), (-1-1-1) depending on reconstruction and chemical potential (range 0.69–1.11 over all facets) | J/m² | DFT-LDA, bare | 45–69 meV/Å² | Moll, Kley, Pehlke, Scheffler, Phys. Rev. B 54, 8844 (1996); arXiv:mtrl-th/9607005 | FT |
| GaAs | 0.86 ± 0.15 | J/m² | measured (fracture), {110} | – | Messmer & Bilello, J. Appl. Phys. 52, 4623 (1981), as cited by Moll 1996 | FT (citation) |
| GaP | 1.40 (111), 1.61 (-1-1-1), 0.74 (110) | J/m² | DFT-GGA, bare | Ga-rich limit; 87.5 / 100.8 / 46.3 meV/Å² | Zhang et al., Sci. Rep. 6, 20055 (2016) | FT |
| PbS | 0.38 (100), 0.45 (110), 0.29 (111, reconstructed) | J/m² | DFT-PBE, bare, vacuum | 24 / 28 / 18 meV/Å²; solvation lowers by ≤7 meV/Å² | Mathew et al., J. Chem. Phys. 140, 084106 (2014), arXiv:1310.4242; numbers from MPInterfaces, arXiv:1602.07784 | FT |
| PbS | −0.618, −0.537 | J/m² | surface **stress** (XRD), oleate-capped | not γ | Bertolotti et al., Nat. Mater. 15, 987 (2016), as quoted by Calvin 2024 | FT (citation) |
| PbSe | 0.184 (100), 0.318 (110), 0.328 (111)-Pb reconstructed | J/m² | DFT, bare | – | Fang, van Huis, Vanmaekelbergh, Zandbergen, ACS Nano 4, 211 (2010), doi:10.1021/nn9013406 | Abs |
| CdS, CdTe, ZnSe, ZnTe, HgTe | **not found** | – | – | CdS (oleate) γ is negative in Calvin 2024, but the value is only in a figure | – | – |

### 2. Bulk solubility C_inf (and solubility proxies)

| Material | Value | Units | Kind | Monomer definition / conditions | Source | Check |
|---|---|---|---|---|---|---|
| CdSe | 1e-2 mol/m³ (= 1e-5 M) | mol/m³ | assumed (C_flat) | generic "monomer" | Talapin 2001 / thesis 2002 | FT (thesis) |
| CdSe | 5.53e-2 mol/m³ (= 5.5e-5 M); a variant gives 5.53e-3 mol/m³ | mol/m³ | fitted | 573 K, Peng 1998 data | Fanelli et al. 2021 (Table 1) | FT (variant: Snip) |
| CdSe | 6.4e-5 M² (ion product, **for 2.2 nm particles**, not bulk) | M² | measured | residual TOPSe × Cd oleate at onset of ripening, ODE, 265 °C | Bullen & Mulvaney 2004 | Abs |
| CdSe | ~0.1–0.2 mM selenourea, below which no nucleation occurs | mM | threshold read by the authors as a solubility limit | Cd oleate 12 mM + oleic acid 24 mM, 240 °C | Hamachi et al., Chem. Sci. 10, 6539 (2019), doi:10.1039/c9sc00989b | FT |
| CdS | ~0.9–1 mM thiourea (nucleation threshold) | mM | same | same | Hamachi 2019 | FT |
| PbS | 2e-10 | mol/L | fitted c_1,∞ | improved kinetic-rate-equation model, fitted to Shrestha PbS data | Wang et al., ACS Omega 6, 3701 (2021), PMC7876681 | FT |
| PbS | 1.0–2.9 mM ([PbS]_ind, 11–32 % of S) | mM | measured solute at nucleation onset (a critical-supersaturation proxy, **not** C_inf) | Pb oleate 10.8 mM, thioureas 9 mM, hexadecane, 80–120 °C | Campos et al., Chem. Sci. 13, 4555 (2022), PMC9019910 | FT |
| PbS | ~4 mM peak solute (~45 % of S) | mM | measured | 110 °C | Abécassis et al., Chem. Sci. 13, 4977 (2022), PMC9067564 | FT |
| PbSe | 0.15–0.58 mM ([PbSe]_ind, 1–7 % of Se) | mM | same proxy | selenoureas, 80–120 °C | Campos 2022 | FT |
| Others | **not found** | – | – | – | – | – |

The "monomer" in these papers is a solute-phase species that is not structurally identified, i.e. "MX solute" formed by precursor conversion (Owen, Hens). Hens' ACS Nano 2013 (doi:10.1021/nn3059168) shows that free carboxylic acid acts like a higher solute solubility (Abs); no number is given.

### 3. Monomer diffusion coefficient D

| Material | Value | Units | Kind | Conditions | Source | Check |
|---|---|---|---|---|---|---|
| CdSe | 1e-12 m²/s (= 1e-8 cm²/s) | m²/s | assumed | model | Talapin thesis 2002 / JPCB 2001 | FT (thesis) |
| CdSe | 3.01e-18 m²/s (variant 3.62e-18) | m²/s | fitted, **physically implausible** (artefact of the fit) | 573 K | Fanelli et al. 2021 | FT |
| PbS (Pb oleate) | 3.6e-10 m²/s (= 3.6e-6 cm²/s) | m²/s | **measured**, DOSY NMR | Pb oleate in hexadecane, 110 °C | Abécassis et al. 2022 | FT |
| PbS | 3.46e-11 m²/s | m²/s | fitted k_G in G = v_m k_G C/r (α = 1); 10× below DOSY | 110 °C | Abécassis et al. 2022 | FT |
| (ODE solvent, reference only) | ~1.1e-9 m²/s | m²/s | DOSY of ODE in CDCl3, room temperature | – | Dhaene et al., "The trouble with ODE", ChemRxiv 2019 (Chem. Mater. 2019) | Snip |
| Cd oleate, Zn, In, Ga, Hg precursors | **not found** | – | – | – | – | – |

De Nolf, …, Hens, JACS 137, 2495 (2015), doi:10.1021/ja509941g, controls CdSe size by changing the carboxylic-acid chain length. They interpret the effect as changing the solute's D and solubility (Abs); the D values are not visible.

### 4. Surface reaction rate constant k_r / Damköhler-type ratio, and regime

| Material | Value | Units | Kind | Conditions | Regime concluded | Source | Check |
|---|---|---|---|---|---|---|---|
| CdSe | 2.2e-6 cm/s (= 220 Å/s, derived conversion) | cm/s | measured, steady-state deposition rate constant | TOPSe + Cd oleate, ODE, ~265 °C | **reaction-limited** ("far slower than the diffusion limit") | Bullen & Mulvaney 2004 | Abs |
| CdSe | 7.97e-10 m/s (= 7.97 Å/s); variant 9.56e-10 | m/s | fitted, coupled to the implausible D | 573 K | D/k ≈ 3.8 nm | Fanelli 2021 | FT |
| CdSe | K = (RT/2γV_m)·D/k_g^flat scanned 1e-3 to 1e3; K < 0.01 means diffusion, K > 100 means reaction | – | dimensionless regime parameter | MC simulations | diffusion control gives better focusing | Talapin 2001 / thesis | FT |
| CdSe | H = D/(k r) = 0.00049 (early), 0.95 (later) | – | fitted | CdO/ODPA + TBPSe, ODE/ODA/TOPO, 280 °C injection, 230 °C growth | diffusion early, then mixed | Su et al., Nanoscale Res. Lett. 5, 823 (2010), PMC2894121 | FT (D, k themselves only in a figure) |
| CdSe | κ_a σ_s / D_1 = 5.03e-2 (tends to 0 for activation-limited growth); κ_a ρ_1,∞ = 7.68e-3 nm⁻² s⁻¹ | – | fitted to Chen et al. 2008 CdSe data | – | leans **reaction** | Kim et al., PNAS 122, e2424950122 (2025), PMC12167951 | FT |
| CdSe | – | – | population balance model | – | **strongly reaction-limited** | van Embden, Sader, Davidson, Mulvaney, J. Phys. Chem. C 113, 16342 (2009), doi:10.1021/jp9027673 | Abs |
| CdSe | – | – | model | – | focusing is possible under pure reaction control | Rempel, Bawendi, Jensen, JACS 131, 4479 (2009), doi:10.1021/ja809156t | Abs |
| CdSe | – | – | in-situ SAXS | Cd oleate / TOPSe hot injection | **reaction**: "superfocusing", surface reactivity drops steeply with size | Prins, …, Hens, Nano Lett. 21, 2487 (2021), doi:10.1021/acs.nanolett.0c04813 | Abs |
| CdSe | – | – | SAXS | Cd stearate + Se, solvothermal | **mixed** (diffusion plus reaction terms) | Biswas, Varghese, Rao, J. Mater. Sci. Technol. 24, 615 (2008) | Abs |
| PbS | k_G = 3.46e-11 m²/s with α = 1 (G = v_m k_G C/r) | m²/s | fitted | Pb oleate / thioureas, hexadecane, 110 °C | **reaction / ligand-shell permeation**, not bulk diffusion (D_fit ≪ DOSY; growth insensitive to a 7× viscosity change) | Abécassis 2022 | FT |
| PbS, PbSe | – | – | – | 80–140 °C | size-dependent, self-limiting attachment; diffusion limitation unlikely | Campos 2022 | FT |
| PbS | k_a' = 4.05e3 s⁻¹ | s⁻¹ | fitted attachment constant | variable T | attachment treated as reaction-limited | Wang et al., ACS Omega 2021 | FT |
| InP | – | – | PBE with homogeneous nucleation | In(My)3 + P(TMS)3 | **reaction-limited** growth assumed; fitted γ, k, Ea not readable | Wang, Traoré, …, Segets, Chem. Eng. Sci. (2023), doi:10.1016/j.ces.2023.119062 | Abs |
| InP | – | – | – | aminophosphines, 180 °C | continuous nucleation (~half the reaction); focusing by size-dependent growth (reaction/surface-limited in practice) | McMurtry, …, Owen, Chem. Mater. 32, 4358 (2020), doi:10.1021/acs.chemmater.0c01561 | Abs/Snip |
| InAs | k = 6 D m_j V_m fitted as "35" (units not given); D not reported | – | diffusion-limited growth law with continuous cluster supply | InAs clusters dripped at 300 °C | **diffusion** (model assumption) | Kim, Park, Jeong, Nat. Commun. 12, 3013 (2021), PMC8140152 | FT |
| CdS (uncapped) | – | – | coarsening r³ ∝ t | – | diffusion-limited LSW coarsening | Viswanatha & Sarma, J. Phys. Chem. Lett. (2010), doi:10.1021/jz9001339 | Abs |
| Others | **not found** | – | – | – | – | – | – |

### 5. Molar volume and prefactors

V_m below is **derived** from room-temperature lattice constants: zinc blende and rock salt both have 4 formula units per cubic cell, so V = a³/4. The lattice constants are standard handbook values (Madelung, *Semiconductors: Data Handbook*, Springer 2004). The literature V_m values that the papers actually used are listed alongside.

| Material | a (Å) | V_fu (Å³) | V_m (cm³/mol) | Literature V_m used |
|---|---|---|---|---|
| CdSe (ZB) | 6.052 | 55.42 | 33.37 | 32.9 (wurtzite; Talapin; Fanelli) = 54.6 Å³ |
| CdSe (WZ) | a 4.300, c 7.011 | 56.13 | 33.80 | – |
| CdS (ZB) | 5.818 | 49.23 | 29.65 | – |
| CdTe | 6.481 | 68.06 | 40.98 | – |
| ZnSe | 5.668 | 45.52 | 27.41 | – |
| ZnS (ZB) | 5.409 | 39.56 | 23.83 | – |
| ZnTe | 6.101 | 56.77 | 34.19 | – |
| HgTe | 6.461 | 67.43 | 40.61 | – |
| InP | 5.869 | 50.54 | 30.44 | – |
| InAs | 6.058 | 55.58 | 33.47 | – |
| GaAs | 5.653 | 45.16 | 27.20 | – |
| GaP | 5.451 | 40.49 | 24.38 | – |
| PbS | 5.936 | 52.29 | 31.49 | 31.5 (Abécassis 2022) |
| PbSe | 6.124 | 57.42 | 34.58 | – |

**Prefactors and other reported kinetic constants:**

- Wang 2021 (PbS):
  - k_a' = 4.05e3 s⁻¹ (fitted).
  - Monomer radius 0.23 nm; monomer number density 19.12 nm⁻³.
  - Scan ranges: k_f 1e-2 to 1e1 s⁻¹; k_a 1e2 to 1e5 s⁻¹; k_D 1e-1 to 1e2 s⁻¹ (FT).
- Kim PNAS 2025 (CdSe): κ_a ρ_1,∞ = 7.68e-3 nm⁻² s⁻¹; Δε_f = 1.24 kT; Δε_e = −1.66 kT (FT).
- Fanelli 2021 (CdSe): N0 = 8.04e21 m⁻³; capillary length 6e-9 m (FT).
- Su 2010 (CdSe): N0 = 3.92e16 cm⁻³; LSW ripening constant 8.32e-5 nm³/s (FT).
- Bullen 2004 (CdSe): nucleus radius 0.8–1.0 nm; 2–8 % of Cd consumed in nucleation (Abs).
- **No classical CNT pre-exponential factor (J0) was found** in any accessible paper.

### 6. Precursor conversion kinetics

| Material | Precursor / system | k (s⁻¹) or rate | Ea | Conditions | Source | Check |
|---|---|---|---|---|---|---|
| CdS | 13 N,N'-substituted thioureas + Cd oleate | 2.5e-3 (S-Im(H,Ph)) down to 1.3e-5 (S-Pym(iPr2)); ±25 %; "reactivity exponent", not strictly first order | not reported | thiourea 10 mM, Cd oleate 12 mM, OA 24 mM, 240 °C | Hamachi et al., Chem. Sci. 10, 6539 (2019) | FT |
| CdSe | 12 selenoureas + Cd oleate | 2.0e-1 (Se-Im(H,Ph)) down to 2.1e-4 (Se-Im(iPr2)) | not reported | same, 240 °C | Hamachi 2019 | FT |
| CdSe | R3P=Se + Cd oleate | initial rate Q0 = 3.0–36 µM/s, linear in [R3PSe]; conversion limits nucleation and growth | not found | – | Owen, Chan, Liu, Alivisatos, JACS 132, 18206 (2010), doi:10.1021/ja106777j | Abs/Snip |
| CdSe | TOPSe/Cd-oleate (microfluidic) | first-order conversion; growth more T-dependent than conversion | not found | – | Maceiczyk, Bezinge, deMello, React. Chem. Eng. 1, 261 (2016) | Abs |
| CdSe (wurtzite) | growth / dissolution | – | Ea(growth) 14 ± 6 kJ/mol; Ea(dissolution) 27 ± 8 kJ/mol | – | Brauser et al., Chem. Mater. 28, 3831 (2016) | Snip |
| PbS | 13 thioureas + Pb oleate | 1e-4 to 1e-1 overall, e.g. TU(Ph,C12): 6.56e-4 (100 °C), 3.5–4.8e-3 (110 °C), 1.02e-2 (120 °C) | not reported; two-point estimates from the SI table ≈ 100–120 kJ/mol (range 43–173), **derived** | Pb oleate 10.8 mM, TU 9 mM, hexadecane, 80–140 °C | Campos et al., Chem. Sci. 13, 4555 (2022), SI §8 | FT |
| PbSe | selenoureas + Pb oleate | SeU(Me2,Cy) 1.3–1.7e-3 (80 °C), 2.7–4.5e-3 (90 °C); SeU(Pip,Cy) 2.41e-2 (80 °C) | derived ≈ 104 kJ/mol | same | Campos 2022 | FT |
| PbS | N-aryl-N'-dodecylthioureas | k_NMR 8.8–13.3e-3; k_UV 2.1–7.2e-3; k_SAXS 2.4–4.3e-3 | – | 110 °C | Abécassis 2022 | FT |
| PbS | thiourea library | spans more than 5 orders of magnitude; Hammett ρ = 1.3 | – | 90–150 °C | Hendricks, Campos, Cleveland, Jen-La Plante, Owen, Science 348, 1226 (2015) | Snip/Abs |
| InP | InCl3 + P(NEt2)3, oleylamine | ~70 % aminophosphine consumed between 2 and 15 min, i.e. k ≈ 1.5e-3 s⁻¹ (**derived**, assuming first order); slower with water (≈ 3.5e-4 and 1e-4 s⁻¹ derived) | – | 180 °C | Vikram et al., Chem. Mater. 32, 3577 (2020) | FT |
| InP | In carboxylate + P(SiMe3)3 | P precursor consumed at injection; computed In–P bond-formation barriers 13 (inter-) / 21 (intra-complex) kcal/mol | AIMD | – | Allen, Walker, Bawendi, Angew. 2010; Xie, Zhao, Jensen, Kulik, arXiv:1512.08565 | Abs |
| InAs | (TMGe)3As vs (TMS)3As | lower reactivity gives a narrower distribution; conversion faster than growth in single injection | – | 295 °C | Harris & Bawendi, JACS 134, 20211 (2012); Franke et al., Nat. Commun. 7, 12749 (2016) | Abs/FT |
| ZnS | Zn ricinoleate, OLA/DDT | – | coarsening Ea 36.2 kJ/mol (weak estimate) | 240–300 °C | Mensah et al., J. Chem. 2022, 9384592 | FT |
| CdTe, ZnSe, ZnTe, HgTe, GaAs, GaP | – | **not found** | – | – | – | – |

---

## (b) Diffusion vs reaction control in the literature

**CdSe.** The early hot-injection literature interpreted size focusing with the diffusion-limited Sugimoto size law:

- Peng, Wickham, Alivisatos, JACS 120, 5343 (1998), on CdSe and InAs.
- Talapin et al. 2001, whose simulations show the best focusing under diffusion control.

The later quantitative work on Cd oleate / phosphine-selenide hot injection consistently finds growth is **reaction-limited, or mixed and dominated by the surface reaction**:

- Bullen & Mulvaney 2004 measured a deposition constant far below the diffusion limit.
- van Embden et al. 2009 concluded growth is "strongly reaction-limited".
- Rempel et al. 2009 showed focusing is possible under pure reaction control.
- Prins/Hens 2021 attribute "superfocusing" to surface reactivity falling with size.
- Kim et al. 2025 fit a ratio that leans towards activation control.

The Hens group (Abé 2012/2013, De Nolf 2015) works with a mixed growth law in which both D and solubility matter. Some solvothermal or TOPO studies report mixed control, or diffusion control early on (Biswas 2008; Su 2010).

**CdS / CdSe from chalcogenoureas (heat-up-like, conversion-limited).** In Owen's chemistry, precursor conversion sets the number of nuclei (Hendricks 2015, Hamachi 2019). Diffusion vs reaction control of growth was not tested explicitly.

**PbS / PbSe (hot injection at low temperature, 80–140 °C).** These give the clearest experimental test:

- Growth is **not** diffusion-limited: the fitted D is 10× below the DOSY value, growth does not respond to a 7× change in viscosity, and growth is self-limiting and size-dependent (Abécassis 2022, Campos 2022).
- Nucleation is persistent, more than 50 % of the reaction time at low T, rather than a LaMer burst.

**InP.**

- With P(SiMe3)3, the precursor is consumed instantly. Growth proceeds by ripening of magic-sized clusters (Cossairt: two-step, via In37P20 clusters; Allen/Bawendi 2010: "non-molecular ripening").
- With aminophosphines, conversion is slow and rate-limiting. Nucleation is continuous, and focusing comes from size-dependent (reaction or surface-limited) growth (McMurtry/Owen 2020).
- The only PBE fit found assumes reaction-limited growth (Wang/Segets 2023).
- No source concludes that InP growth is diffusion-limited.

**InAs.**

- Single injection: conversion is far faster than growth, so growth is ripening-dominated.
- Size focusing needs continuous precursor or cluster supply (Franke/Bawendi 2016; Kim/Jeong 2021). The Jeong group explicitly models it as **diffusion-limited** growth fed by clusters.

**ZnS, ZnSe, ZnTe, CdTe, HgTe, GaAs, GaP.** No quantitative test was found.

- ZnSe, ZnTe and CdTe literature points to reactivity-limited growth, often through magic-sized clusters.
- HgTe growth is anisotropic and involves aggregation.
- GaAs and GaP have no kinetic studies at all.

**Summary.**

- Modern evidence favours **reaction (surface or ligand-shell) limited growth** for the oleate-based hot-injection chemistries of CdSe, PbS, PbSe and InP.
- Diffusion-limited growth survives mainly as the assumption behind the focusing theory, and in continuous-supply InAs models.
- In heat-up and slow-conversion syntheses, the **precursor conversion rate** is the controlling input for the number of particles and their final size.

---

## (c) Recommended literature presets

Same content as `literature_presets.json`.

- **Units:** V_m in Å³ per formula unit; k_r in Å/s, using the growth law dr/dt = V_m k_r (C − C_r); D in cm²/s; C_inf in M.
- **Flags:**
  - "assumed" means the source itself assumed the value rather than measuring or fitting it.
  - "derived" means computed here.
  - "null" means not found.

| Material | γ (J/m²) | C_inf (M) | D (cm²/s) | k_r (Å/s) | control | V_m (Å³) | Notes |
|---|---|---|---|---|---|---|---|
| CdSe | 0.17 (fitted, capped; Bullen 2004) | 1e-5 (**assumed**, Talapin 2001) | 1e-8 (**assumed**, Talapin 2001; likely 100–1000× too low against DOSY-type values) | 220 (Bullen 2004, 2.2e-6 cm/s, 265 °C) | reaction | 55.42 (ZB, derived) | alternatives: γ 0.125 (Talapin, assumed), 0.44 (Fanelli); C_inf 5.5e-5 M (Fanelli fit); calorimetric bare γ 1.31–1.65 |
| CdS | null | null (≈1 mM nucleation threshold, Hamachi 2019, is not C_inf) | null | null | reaction (**assumed**, by analogy) | 49.23 (derived) | conversion k from Hamachi 2019 |
| CdTe | null | null | null | null | reaction (**assumed**) | 68.06 | – |
| ZnSe | null | null | null | null | reaction (**assumed**) | 45.52 | – |
| ZnS | 0.38 (DFT-GGA bare (110); Zhang 2016) | null | null | null | reaction (**assumed**) | 39.56 | calorimetry 1.76 bare / 1.25 hydrated; MD 0.86 |
| ZnTe | null | null | null | null | reaction (**assumed**) | 56.77 | – |
| HgTe | null | null | null | null | mixed (**assumed**; classical model questionable) | 67.43 | – |
| InP | null (0.58 for (110) unverified; −3.71 capped enthalpy) | null | null | null | reaction (inferred: McMurtry 2020; Wang/Segets 2023) | 50.54 | – |
| InAs | 0.65 (DFT-LDA bare (110); Leitsmann 2006) | null | null | null | diffusion (model assumption in Kim/Jeong 2021) | 55.58 | – |
| GaAs | 0.83 (DFT-LDA bare (110); Moll 1996) | null | null | null | mixed (**assumed**; no data) | 45.16 | – |
| GaP | 0.74 (DFT-GGA bare (110); Zhang 2016) | null | null | null | mixed (**assumed**; no data) | 40.49 | – |
| PbS | 0.38 (DFT-PBE bare (100); Mathew 2014) | 2e-10 (fitted, Wang ACS Omega 2021) | 3.6e-6 (DOSY Pb oleate, 110 °C; Abécassis 2022) | null (k_G = 3.46e-11 m²/s with α = 1 instead) | reaction | 52.29 | – |
| PbSe | 0.184 (DFT bare (100); Fang 2010) | null | 3.6e-6 (**assumed**, Pb oleate value from PbS) | null | reaction (Campos 2022) | 57.42 | – |

**Caveat on γ.** CdSe's preset is a fitted value for ligand-capped dots in solution. The other materials' presets are bare-facet DFT values for the low-energy facet, which are upper bounds for capped dots. Calorimetry shows oleate capping can make the surface enthalpy negative. A consistent cross-material model therefore needs either one effective capped γ per material (only CdSe has one in the literature) or a common scaling, which would be an assumption.
