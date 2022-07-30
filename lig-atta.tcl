set sel [lsort -unique [[atomselect top "resname DT"] get residue]]
set sel1 [lrange $sel 0 26]
set sel2 [lrange $sel 27 41]

foreach ind {473 485 455 474 484 407 437 467 372 408 438 468 106 143 184 231 54 86 126 178 55 87 127 125 12 36 77} res $sel1 {
set attachment_pt [measure center [atomselect top "index $ind"]]
set lig_pt1 [measure center [atomselect top "(residue $res) and (name S)"]]
set lig_pt1_ind [[atomselect top "(residue $res) and (name S)"] get index]
set lig_pt2 [measure center [atomselect top "(residue $res) and (name CA)"]]
set cent [measure center [atomselect top "(residue $res) and (type HS)"]]
set vec [vecsub $lig_pt2 $lig_pt1]
set vec_norm [vecnorm $vec]
set lig [atomselect top "residue $res"]
set ang [expr acos(max(-1.0,min(1.0,[vecdot $vec_norm {0 1 0}])))]
set cross [veccross $vec_norm {0 1 0}]
$lig move [eval "trans center {$cent} offset {$attachment_pt} axis {$cross} $ang rad"]
topo addbond $ind $lig_pt1_ind
}



foreach ind {351 312 272 228 313 273 229 177 124 171 212 253 252 211 170} res $sel2 {
set attachment_pt [measure center [atomselect top "index $ind"]]
set lig_pt1 [measure center [atomselect top "(residue $res) and (name S)"]]
set lig_pt1_ind [[atomselect top "(residue $res) and (name S)"] get index]
set lig_pt2 [measure center [atomselect top "(residue $res) and (name CA)"]]
set cent [measure center [atomselect top "(residue $res) and (type HS)"]]
set vec [vecsub $lig_pt2 $lig_pt1]
set vec_norm [vecnorm $vec]
set lig [atomselect top "residue $res"]
set ang [expr acos(max(-1.0,min(1.0,[vecdot $vec_norm {0 -1 0}])))]
set cross [veccross $vec_norm {0 -1 0}]
$lig move [eval "trans center {$cent} offset {$attachment_pt} axis {$cross} $ang rad"]
topo addbond $ind $lig_pt1_ind
}

set need [atomselect top {not type HS}]
$need writepdb t5.pdb
$need writepsf t5.psf

