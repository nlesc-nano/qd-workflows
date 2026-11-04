# Layer path: all 28 intermediates (vs the seed): vmd -e layer_states.tcl  (run from this folder)
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
mol new s00.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {s00 Cd16Se13Cl6 | 0 layers: seed}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol off top
mol new s01.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {s01 Cd19Se13Cl12 | 1 layers: 111B(-1, 1, -1)}
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
mol new s02.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {s02 Cd17Se12Cl10 | 1 layers: 111A(1, 1, -1)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol representation VDW 1.05 30
mol color ColorID 16
mol material Transparent
mol selection {index 16 30 33 36 38}
mol addrep top
mol selection all
mol material Opaque
mol off top
mol new s03.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {s03 Cd22Se13Cl18 | 2 layers: 111B(1, -1, -1), 111B(-1, 1, -1)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol representation VDW 1.05 30
mol color ColorID 16
mol material Transparent
mol selection {index 0 3 4 15 18 19 35 37 38 39 41 45 47 48 49 50 51 52}
mol addrep top
mol selection all
mol material Opaque
mol off top
mol new s04.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {s04 Cd20Se13Cl14 | 2 layers: 111A(1, 1, -1), 111B(-1, 1, -1)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol representation VDW 1.05 30
mol color ColorID 16
mol material Transparent
mol selection {index 0 3 14 19 33 35 36 39 42 44 45 46}
mol addrep top
mol selection all
mol material Opaque
mol off top
mol new s05.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {s05 Cd20Se15Cl10 | 2 layers: 100(0, 0, -1), 111B(-1, 1, -1)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol representation VDW 1.05 30
mol color ColorID 16
mol material Transparent
mol selection {index 0 3 6 15 21 28 36 37 42 43 44}
mol addrep top
mol selection all
mol material Opaque
mol off top
mol new s06.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {s06 Cd25Se13Cl24 | 3 layers: 111B(1, 1, 1), 111B(1, -1, -1), 111B(-1, 1, -1)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol representation VDW 1.05 30
mol color ColorID 16
mol material Transparent
mol selection {index 0 3 4 15 18 19 20 23 24 38 40 41 42 44 47 49 50 52 53 54 55 56 57 58 59 60 61}
mol addrep top
mol selection all
mol material Opaque
mol off top
mol new s07.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {s07 Cd23Se13Cl20 | 3 layers: 111A(1, 1, -1), 111B(1, -1, -1), 111B(-1, 1, -1)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol representation VDW 1.05 30
mol color ColorID 16
mol material Transparent
mol selection {index 0 3 4 15 18 19 22 36 38 39 40 42 44 47 48 50 51 52 53 54 55}
mol addrep top
mol selection all
mol material Opaque
mol off top
mol new s08.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {s08 Cd24Se16Cl16 | 3 layers: 100(0, 0, -1), 111B(1, -1, -1), 111B(-1, 1, -1)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol representation VDW 1.05 30
mol color ColorID 16
mol material Transparent
mol selection {index 0 3 4 7 8 17 20 21 25 29 33 41 42 44 47 49 50 51 52 53 54 55}
mol addrep top
mol selection all
mol material Opaque
mol off top
mol new s09.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {s09 Cd23Se15Cl16 | 3 layers: 100(0, 0, -1), 111B(1, 1, 1), 111B(-1, 1, -1)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol representation VDW 1.05 30
mol color ColorID 16
mol material Transparent
mol selection {index 0 3 6 15 18 21 22 24 31 39 40 43 45 47 48 49 50 51 52 53}
mol addrep top
mol selection all
mol material Opaque
mol off top
mol new s10.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {s10 Cd22Se16Cl12 | 3 layers: 100(0, 0, -1), 111A(1, 1, -1), 111B(-1, 1, -1)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol representation VDW 1.05 30
mol color ColorID 16
mol material Transparent
mol selection {index 0 3 6 15 16 21 23 30 31 39 40 42 44 46 48 49}
mol addrep top
mol selection all
mol material Opaque
mol off top
mol new s11.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {s11 Cd21Se17Cl8 | 3 layers: 100(0, 1, 0), 100(0, 0, -1), 111B(-1, 1, -1)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol representation VDW 1.05 30
mol color ColorID 16
mol material Transparent
mol selection {index 0 3 6 15 17 22 26 30 37 39 43 44 45}
mol addrep top
mol selection all
mol material Opaque
mol off top
mol new s12.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {s12 Cd26Se13Cl26 | 4 layers: 111A(1, 1, -1), 111B(1, 1, 1), 111B(1, -1, -1), 111B(-1, 1, -1)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol representation VDW 1.05 30
mol color ColorID 16
mol material Transparent
mol selection {index 0 3 4 15 18 19 20 23 24 25 39 41 42 43 45 47 49 51 52 53 55 56 57 58 59 60 61 62 63 64}
mol addrep top
mol selection all
mol material Opaque
mol off top
mol new s13.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {s13 Cd27Se16Cl22 | 4 layers: 100(0, 0, -1), 111B(1, 1, 1), 111B(1, -1, -1), 111B(-1, 1, -1)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol representation VDW 1.05 30
mol color ColorID 16
mol material Transparent
mol selection {index 0 3 4 7 8 17 20 21 22 25 26 28 32 36 44 45 47 49 51 52 54 55 56 57 58 59 60 61 62 63 64}
mol addrep top
mol selection all
mol material Opaque
mol off top
mol new s14.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {s14 Cd27Se17Cl20 | 4 layers: 100(0, 0, -1), 111A(1, 1, -1), 111B(1, -1, -1), 111B(-1, 1, -1)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol representation VDW 1.05 30
mol color ColorID 16
mol material Transparent
mol selection {index 0 3 4 7 8 17 18 21 22 23 26 28 32 36 37 45 46 48 49 50 52 54 55 56 58 59 60 61 62 63}
mol addrep top
mol selection all
mol material Opaque
mol off top
mol new s15.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {s15 Cd25Se16Cl18 | 4 layers: 100(0, 0, -1), 111A(1, 1, -1), 111B(1, 1, 1), 111B(-1, 1, -1)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol representation VDW 1.05 30
mol color ColorID 16
mol material Transparent
mol selection {index 0 3 6 15 16 19 22 23 24 26 33 34 42 43 45 47 48 50 51 53 54 55 56 57 58}
mol addrep top
mol selection all
mol material Opaque
mol off top
mol new s16.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {s16 Cd25Se18Cl14 | 4 layers: 100(0, 1, 0), 100(0, 0, -1), 111B(1, -1, -1), 111B(-1, 1, -1)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol representation VDW 1.05 30
mol color ColorID 16
mol material Transparent
mol selection {index 0 3 4 7 8 17 19 21 22 26 30 31 35 42 44 46 48 50 51 52 53 54 55 56}
mol addrep top
mol selection all
mol material Opaque
mol off top
mol new s17.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {s17 Cd24Se17Cl14 | 4 layers: 100(0, 1, 0), 100(0, 0, -1), 111B(1, 1, 1), 111B(1, -1, -1)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol representation VDW 1.05 30
mol color ColorID 16
mol material Transparent
mol selection {index 2 5 16 17 18 19 22 23 28 32 39 40 43 45 46 48 49 50 51 52 53 54}
mol addrep top
mol selection all
mol material Opaque
mol off top
mol new s18.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {s18 Cd24Se20Cl8 | 4 layers: 100(0, 1, 0), 100(0, 0, -1), 111A(1, 1, -1), 111B(-1, 1, -1)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol representation VDW 1.05 30
mol color ColorID 16
mol material Transparent
mol selection {index 0 3 6 15 16 17 19 23 25 29 33 34 41 42 43 45 47 49 51}
mol addrep top
mol selection all
mol material Opaque
mol off top
mol new s19.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {s19 Cd30Se17Cl26 | 5 layers: 100(0, 0, -1), 111A(1, 1, -1), 111B(1, 1, 1), 111B(1, -1, -1), 111B(-1, 1, -1)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol representation VDW 1.05 30
mol color ColorID 16
mol material Transparent
mol selection {index 0 3 4 7 8 17 18 21 22 23 24 27 28 29 31 35 39 40 48 49 51 52 53 55 56 58 59 60 61 63 64 65 66 67 68 69 70 71 72}
mol addrep top
mol selection all
mol material Opaque
mol off top
mol new s20.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {s20 Cd29Se19Cl20 | 5 layers: 100(0, 1, 0), 100(0, 0, -1), 111B(1, 1, 1), 111B(1, -1, -1), 111B(-1, 1, -1)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol representation VDW 1.05 30
mol color ColorID 16
mol material Transparent
mol selection {index 0 3 4 7 8 17 19 21 22 23 24 27 28 30 34 35 39 46 47 49 51 53 54 56 57 58 59 60 61 62 63 64 65 66 67}
mol addrep top
mol selection all
mol material Opaque
mol off top
mol new s21.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {s21 Cd29Se21Cl16 | 5 layers: 100(0, 1, 0), 100(0, 0, -1), 111A(1, 1, -1), 111B(1, -1, -1), 111B(-1, 1, -1)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol representation VDW 1.05 30
mol color ColorID 16
mol material Transparent
mol selection {index 0 3 4 7 8 17 18 19 21 23 24 25 28 30 34 35 39 40 47 48 49 51 53 54 55 57 58 59 61 62 63 64 65}
mol addrep top
mol selection all
mol material Opaque
mol off top
mol new s22.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {s22 Cd27Se19Cl16 | 5 layers: 100(0, 1, 0), 100(0, 0, -1), 111A(1, 1, -1), 111B(1, 1, 1), 111B(1, -1, -1)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol representation VDW 1.05 30
mol color ColorID 16
mol material Transparent
mol selection {index 2 5 16 17 18 19 20 23 24 25 26 31 35 36 43 44 45 48 49 51 52 53 54 56 57 58 59 60 61}
mol addrep top
mol selection all
mol material Opaque
mol off top
mol new s23.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {s23 Cd26Se20Cl12 | 5 layers: 100(1, 0, 0), 100(0, 1, 0), 100(0, 0, -1), 111B(1, -1, -1), 111B(-1, 1, -1)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol representation VDW 1.05 30
mol color ColorID 16
mol material Transparent
mol selection {index 0 3 4 7 8 17 19 21 22 23 27 31 32 36 43 44 45 47 49 51 52 53 54 55 56 57}
mol addrep top
mol selection all
mol material Opaque
mol off top
mol new s24.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {s24 Cd34Se22Cl24 | 6 layers: 100(0, 1, 0), 100(0, 0, -1), 111A(1, 1, -1), 111B(1, 1, 1), 111B(1, -1, -1), 111B(-1, 1, -1)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol representation VDW 1.05 30
mol color ColorID 16
mol material Transparent
mol selection {index 0 3 4 7 8 17 18 19 21 23 24 25 26 27 30 31 32 33 35 39 40 44 45 52 53 54 55 57 59 60 61 63 64 65 66 67 69 70 71 72 73 74 75 76 77 78 79}
mol addrep top
mol selection all
mol material Opaque
mol off top
mol new s25.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {s25 Cd31Se22Cl18 | 6 layers: 100(1, 0, 0), 100(0, 1, 0), 100(0, 0, -1), 111B(1, 1, 1), 111B(1, -1, -1), 111B(-1, 1, -1)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol representation VDW 1.05 30
mol color ColorID 16
mol material Transparent
mol selection {index 0 3 4 7 8 17 19 21 22 23 24 25 28 29 30 32 36 37 41 48 49 50 51 52 54 56 58 59 60 61 62 63 64 65 66 67 68 69 70}
mol addrep top
mol selection all
mol material Opaque
mol off top
mol new s26.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {s26 Cd31Se25Cl12 | 6 layers: 100(1, 0, 0), 100(0, 1, 0), 100(0, 0, -1), 111A(1, 1, -1), 111B(1, -1, -1), 111B(-1, 1, -1)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol representation VDW 1.05 30
mol color ColorID 16
mol material Transparent
mol selection {index 0 3 4 7 8 17 18 19 21 23 24 25 26 27 30 32 36 37 41 42 49 50 51 52 53 54 55 57 59 60 61 63 64 65 66 67}
mol addrep top
mol selection all
mol material Opaque
mol off top
mol new s27.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {s27 Cd38Se28Cl20 | 7 layers: 100(1, 0, 0), 100(0, 1, 0), 100(0, 0, -1), 111A(1, 1, -1), 111B(1, 1, 1), 111B(1, -1, -1), 111B(-1, 1, -1)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol representation VDW 1.05 30
mol color ColorID 16
mol material Transparent
mol selection {index 0 3 4 7 8 17 18 19 21 23 24 25 26 27 28 29 32 33 34 35 36 37 39 43 44 48 49 56 57 58 59 60 61 62 63 64 65 67 69 70 71 73 74 75 76 77 78 79 80 81 82 83 84 85}
mol addrep top
mol selection all
mol material Opaque
mol off top
mol on 0
foreach m [molinfo list] { molinfo $m set rotate_matrix {{{-0.330350 0.943858 0.000000 0.000000} {-0.434873 -0.152205 0.887536 0.000000} {0.837708 0.293198 0.460739 0.000000} {0.000000 0.000000 0.000000 1.000000}}} }
display resetview
# In VMD: Graphics > Representations; Molecule list: toggle the 'D' column to show one structure at a time. The black outline (rep 2) marks the atoms added against the set's reference.
