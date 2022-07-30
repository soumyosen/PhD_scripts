set nf [molinfo 0 get numframes]
set tot 0.0
for {set i 1} {$i <= $nf} {incr i} {
molinfo top set frame $i
set a [atomselect 0 "(segname B1 B2 B3 B4 B5 B6 B7 B8 B9 B10 B11 B12 B13 B14 B15 B16 B17 B18 B19 B20 B21 B22 B23 B24 B25 B26 B27 B28 B29 and resid 12 to 40) and (within 10 of segname NANO N P)" frame $i]
set a_seg [lsort -unique [$a get segname]]
set a_seg_num [llength $a_seg]
set a_sel [atomselect 0 "segname $a_seg and resid 12 to 40" frame $i]
set b_sel [atomselect 1 "segname $a_seg and resid 12 to 40"]
set transform_matrix [measure fit $a_sel $b_sel]
set all [atomselect 0 all frame $i]
$all move $transform_matrix
set a_sel1 [atomselect 0 "segname $a_seg and resid 12 to 40" frame $i]
set rms [measure rmsd $a_sel1 $b_sel]
set rms_seg [expr $rms/$a_seg_num]
puts "$i  $rms_seg"
set tot [expr $tot+$rms_seg]
}
set avg_rms [expr $tot/$nf]
puts "average rmsd $avg_rms"
