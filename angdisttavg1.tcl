 set nf [molinfo top get numframes]
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
 
set au [atomselect top "resname NP3"]
set cau [measure center $au]
set n [atomselect top "resname li3 and type NG"]
set s [atomselect top "resname li3 and type SG"]
set l {}
for {set i 0} {$i < $nf} {incr i} {
 molinfo top set frame $i
 set l2 {}
 set nxyz [$n get {x y z}]
 set sxyz [$s get {x y z}]
# puts "N $i $nxyz"
# puts "S $i $sxyz"
 foreach ncoor $nxyz scoor $sxyz {
 set v1 [vecnorm [vecsub $ncoor $cau]]
 set v2 [vecnorm [vecsub $scoor $ncoor]]
 set ang [expr acos(max(-1.0,min(1.0,[vecdot $v1 $v2])))] 
 set angdeg [expr $ang*57.296]
# puts "$i $angdeg"
#puts " N=$nxyz S=$sxyz "
#set b [llength $nxyz]
#puts "number $b"
 #puts "$i $angdeg"
 if { $i == 0 } {
  lappend l $angdeg
 } else {
  lappend l2 $angdeg
 }
}
if { $i != 0 } {set l [vecadd $l $l2]}
#puts "print $i l: $l\n"
#puts "print $i l2: $l2\n"
}
set lav [vecscale $l [expr double(1)/$nf]]
#puts "total angle $lav" 

 foreach angtav $lav {
 set index [expr int(floor(($angtav-$llimit)/$step))]
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
