# Results and experiments

## 1. CdCl2: the architecture test

On the CdCl2 network nothing grows past the seed in the default case. The dots stay at m = 13 (about 1e-5 M) and the
pool at about 1e-4 M. A 600 s run takes 109 s on 4 threads. **This is a valid result**: CdCl2 tests the architecture
only, and Ivan asked four times never to push it to grow ([decisions](decisions.md)).

## 2. Acetate default case

P1 + P2 network, salt-limited conversion (commit 5af3e3f), NHR job 9646722, 500 K, $c_\mathrm{P,0}$ 0.025 M,
$c_\mathrm{MX,0}$ 0.03 M. 600 s in 447 s of wall time.

| t (s) | α | $C_M$ (M) | $C_X$ (M) | dots (M) | ⟨m⟩ ± σ | edge share | phase |
|---|---|---|---|---|---|---|---|
| 5 | 0.14 | 3.5e-3 | 2.7e-2 | ≈ 0 | – | – | I |
| 14 | 0.30 | 7.2e-3 | 2.2e-2 | 2.5e-6 | 146 ± 38 | 0.15 | II |
| 32 | 0.50 | 4.7e-3 | 1.5e-2 | 5.0e-5 | 150 ± 35 | 0.18 | III |
| 99 | 0.73 | 1.7e-3 | 5.8e-3 | 1.0e-4 | 150 ± 35 | 0.18 | III |
| 194 | 0.81 | 7.7e-4 | 2.7e-3 | 1.3e-4 | 150 ± 35 | 0.18 | III |
| 305 | 0.85 | 4.5e-4 | 1.7e-3 | 1.4e-4 | 149 ± 36 | 0.18 | III |
| 600 | 0.89 | 9.8e-5 | 4.7e-4 | 1.5e-4 | 150 ± 35 | 0.18 | III |

The phases are explained with the maps in [the master equation](master-equation.md) §8:

- **I**: the pool builds.
- **II**: a burst carries new nuclei straight to m ≈ 150, where the map thins.
- **III**: nucleation-limited to the end; the dots do not grow further, the remaining CdSe goes into new dots.

**CdCl2 with the salt-limited conversion**: α 0.979 at 600 s (0.991 before), dots up to 20 % more late in the run,
still at m = 13; 600 s in 123 s.

**The old open-loop run** (NHR job 9644764) agreed to 1e-4 up to 13 s. After about 200 s it ran out of salt in
solution and stalled at 356 s after a spurious second burst ([master equation](master-equation.md) §8.2). Its files
are kept next to the new ones in `network_growth_data/default_salt_vs_open/`.

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
