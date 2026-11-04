# Grid: removing CdSe units (k ≈ 3): vmd -e grid_descent.tcl  (run from this folder)
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
mol new m13_k03.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {(13, 3)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol off top
mol new m12_k03.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {(12, 3)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol off top
mol new m11_k03.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {(11, 3)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol off top
mol new m10_k03.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {(10, 3)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol off top
mol new m09_k03.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {(9, 3)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol off top
mol new m08_k03.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {(8, 3)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol off top
mol new m07_k03.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {(7, 3)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol off top
mol new m06_k03.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {(6, 3)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol off top
mol new m05_k03.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {(5, 3)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol off top
mol new m04_k03.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {(4, 3)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol off top
mol new m03_k03.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {(3, 3)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol off top
mol new m02_k03.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {(2, 3)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol off top
mol new m01_k03.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {(1, 3)}
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
