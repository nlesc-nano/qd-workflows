# Results and experiments

## 1. CdCl2: the architecture test

On the CdCl2 network nothing grows past the seed in the default case. The dots stay at m = 13 (about 1e-5 M) and the
pool at about 1e-4 M. A 600 s run takes 109 s on 4 threads. **This is a valid result**: CdCl2 tests the architecture
only, and Ivan asked four times never to push it to grow ([decisions](decisions.md)).

## 2. Acetate default case

P1 + P2 network, 2026-11-04, NHR job 9644764, 500 K, $c_\mathrm{P,0}$ 0.025 M, $c_\mathrm{MX,0}$ 0.03 M.

| t (s) | α | $C_M$ (M) | $C_X$ (M) | dots (M) | ⟨m⟩ ± σ | phase |
|---|---|---|---|---|---|---|
| 5 | 0.14 | 3.5e-3 | 2.7e-2 | ≈ 0 | – | I |
| 14 | 0.30 | 7.2e-3 | 2.2e-2 | 2.5e-6 | 146 | II |
| 32 | 0.51 | 4.5e-3 | 1.4e-2 | 5.3e-5 | 150 ± 35 | III |
| 99 | 0.79 | 1.0e-3 | 3.5e-3 | 1.2e-4 | 150 ± 35 | III |
| 180 | 0.90 | 1.2e-5 | 1.0e-4 | 1.5e-4 | 151 | III |
| 254 | 0.94 | 2.4e-12 | 2e-16 | 1.5e-4 | 157 | IV |
| 300 | 0.95 | 1.9e-11 | 5e-16 | 1.5e-4 | 159 ± 28 | IV |
| 356 | 0.97 | 1.6e-8 | 9e-15 | 1.8e-4 | 136 ± 59 | V, solver stall |

The phases are explained with the maps in [the master equation](master-equation.md) §8:

- **I**: the pool builds.
- **II**: a burst carries new nuclei straight to m ≈ 150, where the map thins.
- **III**: nucleation-limited; the dots do not grow further.
- **IV**: the salt in solution is used up and the dots ripen between facets.
- **V**: a second burst, an artefact of the conversion ignoring salt depletion.

The run is fast up to the end of the pool (78 s in under 4 min).

## 3. Comparison with experiments

| study | system | nucleation | conversion |
|---|---|---|---|
| Prins, Hens et al., Nano Lett. 2021 (in situ SAXS, SI) | Cd stearate + TOP-Se, 245/230 °C | number of nanocrystals constant after 0.31 min; 12.5 % of the reaction (15–20 % in the main text) | 98 % of the yield in 2.5 min |
| Abécassis et al., Nano Lett. 2015 (in situ SAXS, heat-up) | Cd carboxylate + Se, 240 °C | a 30 s burst | limited by Se activation |
| Owen, Chan, Liu, Alivisatos, JACS 2010 | Cd phosphonate + phosphine selenides | number of nuclei set by the conversion rate | initial rate 3–36 µM/s |
| Hamachi, Owen et al., Chem. Sci. 2019 | Cd oleate + seleno/thioureas, 240 °C | number of nuclei set by the conversion rate | rate constants 1.3e-5 to 0.2 s⁻¹ |

- **Conversion**: half converted at about 30 s, a few times slower than the TOP-Se recipe and within the measured
  range.
- **Nucleation**: starts late (13 s, while the pool builds) and lasts about half the reaction, against 12–20 %
  measured.
- **Growth**: in experiments a fixed number of dots keeps growing. In the model the dots stop at the edge of the map
  (m ≈ 150–180) and the remaining CdSe goes into new dots.

The first two points can be tuned (conversion rate, activities, $s_\nu$), consistent with the model being qualitative
and tunable. The third needs more map above m ≈ 150 (the wave-4 dots and islands on large parents,
[protocol](protocol.md)).
