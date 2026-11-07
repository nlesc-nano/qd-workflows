# Structures of the CdSe protocol docs

One folder per set; every structure is a plain .xyz (Å). The comment line gives the label, the formula and how many atoms were added or removed against the set's reference.

- Open a set in VMD: in the Tk Console (Extensions menu), `cd <set>` then `source <set>.tcl` (or `vmd -e <set>.tcl` from that folder) (VDW spheres, Cd orange, Se blue, camera along the facet of interest; added atoms get a black transparent shell).
- Or open any .xyz in VESTA, Avogadro or OVITO.
- The interactive viewer (rotate, views along each facet, added atoms highlighted): `viewer.html` next to this folder, or the Artifact page.

| set | structures | reference for 'added' |
|---|---|---|
| `grid_seeds` — Grid: start dot and seeds | 4 | no reference |
| `grid_row13` — Grid: the m = 13 row | 10 | no reference |
| `grid_descent` — Grid: removing CdSe units (k ≈ 3) | 13 | no reference |
| `builder_stages` — Builder: the seed through the builder stages | 3 | previous stage |
| `layer_path` — Layer path: seed → Cd38Se28Cl20, one plane per step | 8 | previous step |
| `layer_states` — Layer path: all 28 intermediates (vs the seed) | 28 | the seed |
| `islands_111B` — Islands on {111}B: parent Cd54Se43Cl22 | 16 | the parent |
| `islands_111A` — Islands on {111}A: parent Cd99Se84Cl30 | 10 | the parent |
| `islands_100` — Islands on {100}: parent Cd42Se26Cl32 | 3 | the parent |
| `cd_shells` — Cd shells k0 ± 1, ± 2 (z_type_target) | 5 | the k0 build |
| `acetate_dressing` — Acetate dressing of two CdCl2 builder dots | 6 | no reference |
