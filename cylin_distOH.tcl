set crd {}
set PI 3.1415926535897931
set radius 23.0
for {set i 4} {$i <= 156} {incr i 4} {
 for {set j 1} {$j <= 26} {incr j} {
  set x [expr [expr $PI*$j*360]/[expr 26*180]]
  set coord_x [expr $radius*cos($x)]
  set y [expr [expr $PI*$j*360]/[expr 26*180]]
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

set sioh_residue [lsort -unique [[atomselect top "segname E"] get residue]]
foreach resid $sioh_residue pt $crd {
set atoms [atomselect top "residue $resid"]
#set C [atomselect top "residue $residue and type CG30"]
set Ccrd [measure center $atoms]
$atoms moveby [vecscale [vecsub $Ccrd $pt] -1]
set z [lindex $pt 2]
set cent [concat "0.0 0.0 $z"]
set si_crd [measure center [atomselect top "residue $resid and type OHT"]]
set h_crd [measure center [atomselect top "residue $resid and type HHT"]]
set v1 [vecsub $h_crd $si_crd]
set v1norm [vecnorm $v1]
set v2 [vecsub $si_crd $cent]
set v2norm [vecnorm $v2]
set ang [expr acos(max(-1.0,min(1.0,[vecdot [vecnorm $v1norm] [vecnorm $v2norm]])))]
#puts "$ang"
set cross "{[veccross $v1 $v2]}"
#puts "$cross"
$atoms move [eval "trans center {$si_crd} axis $cross $ang rad"]
$atoms delete
}


