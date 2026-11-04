# Layer path: seed → Cd38Se28Cl20, one plane per step: vmd -e layer_path.tcl  (run from this folder)
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
mol new path_0_s00.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {s00 seed}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol off top
mol new path_1_s01.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {s01 + 111B(-1, 1, -1)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol representation VDW 1.05 30
mol color ColorID 16
mol material Transparent
mol selection {index 0 3 14 32 34 35 41 42 43}
mol addrep top
mol selection all
mol material Opaque
mol off top
mol new path_2_s05.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {s05 + 100(0, 0, -1)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol representation VDW 1.05 30
mol color ColorID 16
mol material Transparent
mol selection {index 6 21 28 42}
mol addrep top
mol selection all
mol material Opaque
mol off top
mol new path_3_s10.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {s10 + 111A(1, 1, -1)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol representation VDW 1.05 30
mol color ColorID 16
mol material Transparent
mol selection {index 15 21 31 42 44 46 48 49}
mol addrep top
mol selection all
mol material Opaque
mol off top
mol new path_4_s18.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {s18 + 100(0, 1, 0)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol representation VDW 1.05 30
mol color ColorID 16
mol material Transparent
mol selection {index 17 19 29 41 42 43 49}
mol addrep top
mol selection all
mol material Opaque
mol off top
mol new path_5_s21.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {s21 + 111B(1, -1, -1)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol representation VDW 1.05 30
mol color ColorID 16
mol material Transparent
mol selection {index 4 7 23 24 25 35 53 55 58 59 62 63 64 65}
mol addrep top
mol selection all
mol material Opaque
mol off top
mol new path_6_s26.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {s26 + 100(1, 0, 0)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol representation VDW 1.05 30
mol color ColorID 16
mol material Transparent
mol selection {index 26 27 52 53 54 55 64}
mol addrep top
mol selection all
mol material Opaque
mol off top
mol new path_7_s27.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {s27 + 111B(1, 1, 1)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol representation VDW 1.05 30
mol color ColorID 16
mol material Transparent
mol selection {index 23 24 32 33 35 36 37 59 64 65 73 75 76 81 82 83 84 85}
mol addrep top
mol selection all
mol material Opaque
mol off top
mol on 0
foreach m [molinfo list] { molinfo $m set rotate_matrix {{{-0.330350 0.943858 0.000000 0.000000} {-0.434873 -0.152205 0.887536 0.000000} {0.837708 0.293198 0.460739 0.000000} {0.000000 0.000000 0.000000 1.000000}}} }
display resetview
# In VMD: Graphics > Representations; Molecule list: toggle the 'D' column to show one structure at a time. The black outline (rep 2) marks the atoms added against the set's reference.
