 set crd {}
 set range 1000
 set PI 3.1415926535897931
 set radius 23.0
 for {set i 1} {$i <= $range} {incr i} {
 set x [expr $PI*(rand()*360)/180]
 set coord_x [expr $radius*cos($x)]
 set y [expr $PI*(rand()*360)/180]
 set coord_y [expr $radius*sin($y)]
 set coord_z [expr rand()*155]
 set coord [concat "$coord_x $coord_y $coord_z"]
 lappend crd $coord
}

#puts "$crd"

puts "$crd"
puts [llength $crd]
mol load graphics testing
foreach pt $crd {
graphics top point $pt
}




#set sioh_residue [lsort -unique [[atomselect top "segname W"] get residue]]
#foreach resid $sioh_residue pt $crd {
#set atoms [atomselect top "residue $resid"]
##set C [atomselect top "residue $residue and type CG30"]
#set Ccrd [measure center $atoms]
#$atoms moveby [vecscale [vecsub $Ccrd $pt] -1]
##}
##
