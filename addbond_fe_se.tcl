set iron [atomselect top {name Fe}]
set iron_ind [$iron get index]

foreach ind $iron_ind {
set selenium [atomselect top "name Se and (within 2.4 of index $ind)"]
set sele_ind [$selenium get index]
foreach s_ind $sele_ind {
topo addbond $ind $s_ind
}
}

set all [atomselect top all]
$all writepdb t1.pdb
$all writepsf t1.psf
 
