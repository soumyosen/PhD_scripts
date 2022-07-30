set nf [molinfo 0 get numframes]
set conv 57.2958
set tot 0.0

proc segtonum {name} {
set l [split "$name" {}]
set str [join "[lindex $l 1][lindex $l 2]"]
return $str
}

#set all {}

for {set i 0} {$i < $nf} {incr i} {
molinfo top set frame $i
#set a [atomselect 0 "(segname B1 B2 B3 B4 B4 B6 B7 B8 B9 B10 B11 B12 B13 B14 B15 B16 B17 B18 B19 B20 B21 B22 B23 B24 B25 B26 B27 B28 B29 H1 H2 H3 H4 H5 H6 H7 H8 H9 H10 H11 H12 H13 H14 H15 H16 H17 H18 H19 H20 H21 H22 H23 H24 H25 H26 H27 H28 H29) and (within 10 of segname NANO N P)" frame $i]
#set a [atomselect 0 "(segname B1 B2 B3 B4 B4 B6 B7 B8 B9 B10 B11 B12 B13 B14 B15 B16 B17 B18 B19 B20 B21 B22 B23 B24 B25 B26 B27 B28 B29 H1 H2 H3 H4 H5 H6 H7 H8 H9 H10 H11 H12 H13 H14 H15 H16 H17 H18 H19 H20 H21 H22 H23 H24 H25 H26 H27 H28 H29) and (within 10 of segname NANO N P)"]
#set a_seg [lsort -unique [$a get segname]]
#set B {}
#set H {}

#foreach word $a_seg {
#set l1 [split "$word" {}]
#set temp [lindex $l1 0]
#puts "$temp"

#if {[lindex $l1 0] == {B}} {
#	lappend B $word 
#} else {
#	lappend H $word}
#}

#puts "$B"
#puts "$H"

#set u [lsort -dictionary $B]
#puts "$u"

#set e [lindex $u 0]
set e_num [segtonum B4]
#set f [lindex $u end]
set f_num [segtonum B25]
#puts "$e  $f"

set tot_perfr 0.0
for {set j 4} {$j < 25} {incr j} {
set a1 [measure center [atomselect top "segname B$j and resid 18 and backbone and name CA"]]
set a2 [measure center [atomselect top "segname B$j and resid 32 and backbone and name CA"]]
set v1 [vecnorm [vecsub $a2 $a1]]

set k [expr $j+1]
set b1 [measure center [atomselect top "segname B$k and resid 18 and backbone and name CA"]]
set b2 [measure center [atomselect top "segname B$k and resid 32 and backbone and name CA"]]
set v2 [vecnorm [vecsub $b2 $b1]]

set ang [expr acos(max(-1.0,min(1.0,[vecdot $v1 $v2])))]
set ang_deg [expr $ang*$conv]
puts "$ang_deg"
#set tot_perfr [expr $tot_perfr+$ang_deg]

}
#puts "total per frame $tot_perfr"
#set num [expr $f_num-$e_num]
#puts "number $num"
#set ang_deg_perpep [expr $tot_perfr/$num]
#set tot [expr $tot+$ang_deg_perpep]
#puts "$i  $ang_deg_perpep"   
#$a delete
}

#set avg_ang [expr $tot/$nf]
#puts "average angle $avg_ang"
