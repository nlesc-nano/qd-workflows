# Stages D–E: the acetate network and the export

## 1. The acetate network

The acetate network uses the same registry and the same steps as CdCl2. Only the structures and their free energies
change.

- **Inputs**: the CdCl2 **builder** `out.xyz` of every dot (checked by md5), named by the registry `dir`, or by the
  island `tag` for islands.
- **Pipeline** (`p1.py relax GRID_WD MODEL LIST OUT`), per dot:
  1. dress with acetate ([builder](builder.md) §4);
  2. relax, then polish (MACE-MH-1 + D3(BJ), float64);
  3. full FD Hessian (the inner region is the whole acetate dot);
  4. G on the grid's temperatures.

  With `P1_NO_HESSIAN=1` it stops after the polish, and F_vib comes from the learned model.
- **No anneal above m = 13**; at most one short MD per dot, with one seed. Ivan: "I would just use a very short MD to
  accommodate better the acetates on the surface not to reconstruct the full structure. Otherwise it becomes a mess."
  P0 showed no MD variant beats relaxation alone (−0.08 and +0.09 eV mean, ±0.6 eV scatter), so P1 and P2 were
  relaxation only.

| phase | m | where | inputs | result | cost |
|---|---|---|---|---|---|
| P0 | pilot: m ≤ 13 plus tests | LUMI small-g | 41 + 34 | pipeline works for 71–1,544 atoms; network dots with m ≤ 13 lie 1.24 eV above the grid's lowest | small |
| P1 | 14–60 + islands | NHR gpu-a100, 4 batches per node (LUMI was reserved) | 8,199 | 8,169 states (30 dressing failures), FD Hessians | 45.5 node-h |
| P2, NHR | 61–223 + islands, batches b000–b063 and b287–b319 | NHR gpu-a100 + test queue | 3,027 | 3,026 states, FD Hessians | about 32k core-h |
| P2, LUMI | batches b064–b286 | LUMI standard-g, 8 GCDs per node | 6,983 | relaxation only, learned F_vib | 48.3 node-h = 193 GPU-h |

- **Dot selection**: every CdCl2-relaxed registry dot in the m range (7,849 at m 14–60, 9,681 at m 61–223). The
  script that chose the island inputs was not found.
- **Cost model from P0**: Hessian 2.28 ms·N², relax + polish 0.30 s per atom, MD 0.19 s per atom per ps. One A100 ≈
  3.9 MI250X GCDs.
- **Relaxation noise**: the same start relaxed on NHR and on LUMI ends in different minima, 0.27 eV rms.
- **Total**: 34 + 7,819 + 9,680 = **17,533 acetate dots**.

## 2. Export

The export runs on NHR cpu-clx, one node with 96 cores (`job_net_p2.sh`): about 1 h, or 7 min when the G cache
exists.

```
FVIB_PRED=predictions.json python ac_build.py net states_index.json thermo_acetate_pool.json cdse_zb_acetate.yaml TREE 96 \
    p0_out_A_states.json P1/b*/states.json P2/b???/states.json P2/b???q?/states.json
cp layer_steps_cdcl2_fill.json TREE/layer_steps.json
python attach_dg.py net TREE/results_by_key.json TREE/g_cache.json TREE/layer_steps.json
python ladders.py TREE/layer_steps.json TREE/ladders.json
python -m qd_workflows.network_export --registry net/registry.json --results TREE/results_by_key.json \
    --g-cache TREE/g_cache.json --steps TREE/layer_steps.json --ladders TREE/ladders.json \
    --islands states_index.json --island-states TREE/island_states.json \
    --thermo thermo_acetate_pool.json --params cdse_zb_acetate.yaml --out OUT --nproc 96
```

| script | does |
|---|---|
| `ac_build.py` | maps each relaxed state to its registry key (dots) or island name; writes `results_by_key.json`, `island_states.json` and `g_cache.json` (G at T = 500 K, ε = 2 including σ·S; the learned F_vib where there is no Hessian) |
| `attach_dg.py` | adds $\Delta G=G_b-G_a$ to every registry step whose ends both have G |
| `ladders.py` | builds the ladders ([islands](islands.md) §3) |
| `network_export` | assembles the network (below) |

`network_export` does the following:

- **Steps**: steps between relaxed dots, deduplicated, oriented so that (m, k) grows, and typed:
  - `layer`: Δm > 0 along one facet family, and not the continuation of a started sheet;
  - `attach`: any other step with Δm > 0;
  - `shell`: Δm = 0, Δk ≠ 0;
  - `iso`: the same composition.

  Extra shell links join, per Se framework, the lowest isomer of each k to the lowest of the next k.
- **What is kept**: the connected component of the seed (registry lineage 0 = dot 12 = Cd16Se13Cl6).
- **Islands**: the island rows with status todo or relaxed, pure, 1–10 new Cd and parent in the component, plus the
  ladder rungs of those (parent, facet) pairs.
- **Free energy**: `ideal_gas_g` with σ = 1 on the thermo temperatures, plus the vdW surface. A state without
  frequencies gets translation + rotation + pV + the learned F_vib. σ is recovered from the cached G and added as
  $k_BT\ln\sigma$.
- **Files**:
  - `network.json`: states {id, formula, m, k, fw, n_imag, fvib}, steps {a, b, dm, dk, kind, facet, sheet}, islands,
    temperatures, seed, element codes Cd Se Cl C O H;
  - `arrays.npz`: G_gas on the temperatures, S, σ;
  - `geometry.npz`: every structure.

| export | states | steps | island rows |
|---|---|---|---|
| CdCl2 tree + fill | 15,160, m 7–223 | 78,838 | 578 |
| acetate P1 | 7,491, m 7–60 | 34,731 | 179 |
| acetate P1 + P2 | 15,131, m 7–223 | 78,870: layer 26,432, attach 37,835, shell 10,518, iso 4,085 | 576 |

In the acetate P1 + P2 export, 5,595 of 15,707 structures use the learned F_vib. σ is 1 for 15,126 states, 2 for 4
and 3 for 1. No cached G is off the recomputed one by more than 2 meV. Laptop copy:
`network_growth_data/cdse_acetate_p2`.
