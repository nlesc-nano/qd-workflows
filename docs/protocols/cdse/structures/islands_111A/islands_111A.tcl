# Islands on {111}A: parent Cd99Se84Cl30: vmd -e islands_111A.tcl  (run from this folder)
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
mol new parent.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {parent Cd99Se84Cl30}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol off top
mol new p6281_111A_ppm_s1_20_22_26.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {+1 Cd (Cd100Se83Cl34)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol representation VDW 1.05 30
mol color ColorID 16
mol material Transparent
mol selection {index 78 166 192 201 202 210 211}
mol addrep top
mol selection all
mol material Opaque
mol off top
mol new p6281_111A_ppm_s1_20_24_24.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {+1 Cd (Cd100Se83Cl34)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol representation VDW 1.05 30
mol color ColorID 16
mol material Transparent
mol selection {index 90 166 192 201 204 210 211}
mol addrep top
mol selection all
mol material Opaque
mol off top
mol new p6281_111A_ppm_s2_20_24_26.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {+2 Cd (Cd101Se83Cl36)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol representation VDW 1.05 30
mol color ColorID 16
mol material Transparent
mol selection {index 78 91 151 168 193 202 205 210 212 213 214}
mol addrep top
mol selection all
mol material Opaque
mol off top
mol new p6281_111A_ppm_s2_22_22_26.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {+2 Cd (Cd101Se83Cl36)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol representation VDW 1.05 30
mol color ColorID 16
mol material Transparent
mol selection {index 78 80 167 176 193 202 203 210 212 213 214}
mol addrep top
mol selection all
mol material Opaque
mol off top
mol new p6281_111A_ppm_s3_22_24_24.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {+3 Cd (Cd102Se83Cl38)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol representation VDW 1.05 30
mol color ColorID 16
mol material Transparent
mol selection {index 79 91 92 168 178 194 203 206 211 213 214 215 216 217}
mol addrep top
mol selection all
mol material Opaque
mol off top
mol new p6281_111A_ppm_s3_24_24_22.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {+3 Cd (Cd102Se83Cl38)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol representation VDW 1.05 30
mol color ColorID 16
mol material Transparent
mol selection {index 90 97 101 167 177 194 201 210 212 213 214 216 220 222}
mol addrep top
mol selection all
mol material Opaque
mol off top
mol new p6281_111A_ppm_s4_22_24_26.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {+4 Cd (Cd103Se84Cl38)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol representation VDW 1.05 30
mol color ColorID 16
mol material Transparent
mol selection {index 78 80 92 93 153 170 180 196 205 208 213 215 216 217 218 219}
mol addrep top
mol selection all
mol material Opaque
mol off top
mol new p6281_111A_ppm_s5_24_24_24.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {+5 Cd (Cd104Se84Cl40)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol representation VDW 1.05 30
mol color ColorID 16
mol material Transparent
mol selection {index 79 91 92 99 103 169 179 180 196 198 207 210 215 217 218 219 221 225 227}
mol addrep top
mol selection all
mol material Opaque
mol off top
mol new p6281_111A_ppm_s6_24_24_26.xyz type xyz waitfor all
mol delrep 0 top
mol representation VDW 1.0 30
mol color Element
mol addrep top
mol rename top {+6 Cd (Cd105Se85Cl40)}
set s [atomselect top "name Cd"]; $s set radius 1.58; $s set element Cd; $s delete
set s [atomselect top "name Se"]; $s set radius 1.9; $s set element Se; $s delete
set s [atomselect top "name Cl"]; $s set radius 1.75; $s set element Cl; $s delete
set s [atomselect top "name O"]; $s set radius 1.52; $s set element O; $s delete
set s [atomselect top "name C"]; $s set radius 1.7; $s set element C; $s delete
set s [atomselect top "name H"]; $s set radius 1.2; $s set element H; $s delete
mol representation VDW 1.05 30
mol color ColorID 16
mol material Transparent
mol selection {index 78 80 92 93 100 104 154 171 181 182 198 200 209 212 217 219 220 221 223 227 229}
mol addrep top
mol selection all
mol material Opaque
mol off top
mol on 0
foreach m [molinfo list] { molinfo $m set rotate_matrix {{{-0.707107 0.707107 0.000000 0.000000} {0.408248 0.408248 0.816497 0.000000} {0.577350 0.577350 -0.577350 0.000000} {0.000000 0.000000 0.000000 1.000000}}} }
display resetview
# In VMD: Graphics > Representations; Molecule list: toggle the 'D' column to show one structure at a time. The black outline (rep 2) marks the atoms added against the set's reference.
