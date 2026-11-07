# Acetate dressing of two CdCl2 builder dots: vmd -e acetate_dressing.tcl  (run from this folder)
display projection Orthographic
color Display Background white
axes location Off
color Element Cd 3
color Element Se 0
color Element Cl 7
color Element O 1
color Element C 2
color Element H 8
color change rgb 3 0.94 0.63 0.19
color change rgb 0 0.18 0.44 0.87
mol new c000015_cdcl2_builder.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {Cd16Se13Cl6 CdCl2 builder}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol off top
mol new c000015_acetate_dressed.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {c000015 acetate dressed (start)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol off top
mol new c000015_acetate_relaxed.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {c000015 acetate relaxed (MACE-MH-1 + D3)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol off top
mol new c000036_cdcl2_builder.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {Cd25Se22Cl6 CdCl2 builder}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol off top
mol new c000036_acetate_dressed.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {c000036 acetate dressed (start)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol off top
mol new c000036_acetate_relaxed.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {c000036 acetate relaxed (MACE-MH-1 + D3)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol off top
mol on 0
foreach m [molinfo list] { molinfo $m set rotate_matrix {{{-0.707107 0.707107 0.000000 0.000000} {-0.408248 -0.408248 0.816497 0.000000} {0.577350 0.577350 0.577350 0.000000} {0.000000 0.000000 0.000000 1.000000}}} }
display resetview
# In VMD: Graphics > Representations; Molecule list: toggle the 'D' column to show one structure at a time. The black outline (rep 2) marks the atoms added against the set's reference.
