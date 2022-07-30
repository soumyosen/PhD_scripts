set iron [atomselect top {name Se}]
set iron_ind [$iron get index]

foreach ind $iron_ind {
set selenium [atomselect top "name Se and (within 2.56 of index $ind)"]
set sele_ind [$selenium get index]
foreach s_ind $sele_ind {
if {$ind == $s_ind} {
 continue} else {
topo addbond $ind $s_ind
}
}
}
set all [atomselect top all]
$all writepdb t4.pdb
$all writepsf t4.psf
