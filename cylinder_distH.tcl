set crd {}
set PI 3.1415926535897931
set radius 20.0
for {set i 2} {$i <= 50} {incr i 2} {
 for {set j 1} {$j <= 32} {incr j} {
  set x [expr [expr $PI*$j*360]/[expr 32*180]]
  set coord_x [expr $radius*cos($x)]
  set y [expr [expr $PI*$j*360]/[expr 32*180]]
  set coord_y [expr $radius*sin($y)]
  set coord [concat "$coord_x $coord_y $i"]
  lappend crd $coord
}
}

#puts "$crd"
#puts [llength $crd]
#mol load graphics testing
#foreach pt $crd {
#graphics top point $pt
#}

set sioh_residue [lsort -unique [[atomselect top "segname D"] get residue]]
foreach resid $sioh_residue pt $crd {
set atoms [atomselect top "residue $resid"]
#set C [atomselect top "residue $residue and type CG30"]
set Ccrd [measure center $atoms]
$atoms moveby [vecscale [vecsub $Ccrd $pt] -1]
}
