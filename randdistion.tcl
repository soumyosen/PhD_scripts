 set crd {}
 set range 501
 for {set i 1} {$i < $range} {incr i} {
 set numx [expr rand()*310]
 set numy [expr rand()*310]
 set numz [expr rand()*140]
 set coord [concat "$numx $numy $numz"]
 lappend crd $coord
}

 set ion_residue [lsort -unique [[atomselect top "segname Gi"] get residue]]
 foreach residue $ion_residue pt $crd {
 set atoms [atomselect top "residue $residue"]
 set C [atomselect top "residue $residue and type CG30"]
 set Ccrd [measure center $C]
 $atoms moveby [vecscale [vecsub $Ccrd $pt] -1]
}
