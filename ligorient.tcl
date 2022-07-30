 set nf [molinfo top get numframes]
 set llimit 0.0
 set ulimit 180.0
 set step 5.0
 set num_bin [expr ($ulimit-$llimit)/$step]
 set histogram {};set histogram_x {}
 for {set i 1} {$i <= [expr $num_bin+1]} {incr i} {
 lappend histogram 0.0
 lappend histogram_x [expr ($i-1)*$step]
 }

 set au [atomselect top "segname NANO and resname NP"]
 set cau [measure center $au]
 set l {}
 for {set i 0} {$i < $nf} {incr i} {
 molinfo top set frame $i
 set oatom [atomselect top "segname NANO and type CG30"]
 set ocoord [$oatom get {x y z}]
 set catom [atomselect top "segname NANO and type CG3F"]
 set ccoord [$catom get {x y z}]
 
 foreach o1r $ocoord c1r $ccoord {
 set v1 [vecsub $c1r $cau]
 set v1norm [vecnorm $v1]
 set v2 [vecsub $o1r $c1r]
 set v2norm [vecnorm $v2]
 set ang [expr acos(max(-1.0,min(1.0,[vecdot $v2norm $v1norm])))]
 set angdeg [expr $ang*57.296]
 lappend l $angdeg

}
 $oatom delete
 $catom delete
}

 foreach angtav $l {
 set index [expr int(floor(($angtav-$llimit)/$step))]
 set old_occurance [lindex $histogram $index]
 set new_occurance [expr $old_occurance+1]
 set histogram [lreplace $histogram $index $index $new_occurance]
 puts "angle:$angtav, index: $index, o1: $old_occurance, o2: $new_occurance"
}
 
 puts $histogram
 puts $histogram_x


 
$au delete

