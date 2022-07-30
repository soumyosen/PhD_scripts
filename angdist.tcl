# build empty histogram
 set llimit 0.0
 set ulimit 180.0
 set step 10.0
 set num_bin [expr ($ulimit-$llimit)/$step]
 set histogram {};set histogram_x {}
 for {set i 1} {$i <= $num_bin} {incr i} {
 lappend histogram 0.0
 lappend histogram_x [expr ($i-1)*$step+$llimit+$step/2]
 }
 
set au [atomselect top "resname NP7"]
set cau [measure center $au]
set n [atomselect top "resname li7 and type NG"]
set nxyz [$n get {index x y z}]
set s [atomselect top "resname li7 and type SG"]
set sxyz [$s get {index x y z}]

foreach coord $nxyz {
set ptN [lrange $coord 1 3]
set v1 [vecsub $ptN $cau]
set v1nor [vecnorm $v1]
#puts "$v1nor"
set indN [lindex $coord 0]
#puts "thanks"
foreach indexs $sxyz {
set indS [lindex $indexs 0]
if {[expr ($indS-$indN)] == 13} {set ptS [lrange $indexs 1 3]}
}
set v2 [vecsub $ptS $ptN]
set v2nor [vecnorm $v2]
#puts "$v2nor"
set ang [expr acos(max(-1.0,min(1.0,[vecdot $v1nor $v2nor])))]
set angdeg [expr $ang*57.296]
puts "$angdeg"
set index [expr int(floor(($angdeg-$llimit)/$step))]
set old_occurance [lindex $histogram $index]
set new_occurance [expr $old_occurance+1]
set histogram [lreplace $histogram $index $index $new_occurance]
puts "angle:$angdeg, index: $index, o1: $old_occurance, o2: $new_occurance"
puts $histogram
puts $histogram_x

} 

# build empty histogram
#set llimit 0.0
#set ulimit 180.0
#set step 10.0
#set num_bin [expr ($ulimit-$llimit)/$step]
#set histogram {};set histogram_x {}
#for {set i 1} {$i <= $num_bin} {incr i} {
#lappend histogram 0.0
#lappend histogram_x [expr ($i-1)*$step+$llimit+$step/2]
#}

#set angdeg 3.4
#set index [expr int(floor(($angdeg-$llimit)/$step))]
#set old_occurance [lindex $histogram $index]
#set new_occurance [expr $old_occurance+1]
#set histogram [lreplace $histogram $index $index $new_occurance]
#puts "angle:$angdeg, index: $index, o1: $old_occurance, o2: $new_occurance"
#puts $histogram
#puts $histogram_x
