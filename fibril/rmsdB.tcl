set nf [molinfo 0 get numframes]

mol load psf last_onlyfib.psf pdb last_onlyfib.pdb

set b_sel [atomselect 1 "segname B1 B2 B3 B4 B5 B6 B7 B8 B9 B10 B11 B12 B13 B14 B15 B16 B17 B18 B19 B20 B21 B22 B23 B24 B25 B26 B27 B28 B29 H1 H2 H3 H4 H5 H6 H7 H8 H9 H10 H11 H12 H13 H14 H15 H16 H17 H18 H19 H20 H21 H22 H23 H24 H25 H26 H27 H28 H29 and resid 12 to 40 and backbone"]

set tot 0.0
for {set i 0} {$i < $nf} {incr i} {
molinfo top set frame $i
#set a [atomselect 0 "(segname B1 B2 B3 B4 B5 B6 B7 B8 B9 B10 B11 B12 B13 B14 B15 B16 B17 B18 B19 B20 B21 B22 B23 B24 B25 B26 B27 B28 B29 H1 H2 H3 H4 H5 H6 H7 H8 H9 H10 H11 H12 H13 H14 H15 H16 H17 H18 H19 H20 H21 H22 H23 H24 H25 H26 H27 H28 H29) and (within 10 of segname NANO P)" frame $i]
set a_seg [concat "B4 B5 B6 B7 B8 B9 B10 B11 B12 B13 B14 B15 B16 B17 B18 B19 B20 B21 B22 B23 B24 B25"]
set a_seg_num [llength $a_seg]
set a_sel [atomselect 0 "segname B1 B2 B3 B4 B5 B6 B7 B8 B9 B10 B11 B12 B13 B14 B15 B16 B17 B18 B19 B20 B21 B22 B23 B24 B25 B26 B27 B28 B29 H1 H2 H3 H4 H5 H6 H7 H8 H9 H10 H11 H12 H13 H14 H15 H16 H17 H18 H19 H20 H21 H22 H23 H24 H25 H26 H27 H28 H29 and resid 12 to 40 and backbone" frame $i]
set transform_matrix [measure fit $a_sel $b_sel]
set all [atomselect 0 all frame $i]
$all move $transform_matrix

set tot_perpep 0.0
foreach seg $a_seg {
set a_sel1 [atomselect 0 "segname $seg and resid 12 to 40 and backbone" frame $i]
set b_sel1 [atomselect 1 "segname $seg and resid 12 to 40 and backbone"]
set rms [measure rmsd $a_sel1 $b_sel1]
set tot_perpep [expr $tot_perpep+$rms]
}

set rms_seg [expr $tot_perpep/$a_seg_num]
#puts "number $a_seg_num"
puts "$i  $rms_seg"
set tot [expr $tot+$rms_seg]
#$a delete
$a_sel delete
$a_sel1 delete
$b_sel1 delete
$all delete
}
set avg_rms [expr $tot/$nf]
puts "average rmsd $avg_rms"
$b_sel delete
