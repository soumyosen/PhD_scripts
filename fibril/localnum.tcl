########### Measuring atom numbers close to the nanoparticle
set nf [molinfo top get numframes]
set tot 0.0
#set total 34684
set total [atomselect top "(segname B4 B5 B6 B7 B8 B9 B10 B11 B12 B13 B14 B15 B16 B17 B18 B19 B20 B21 B22 B23 B24 B25 H4 H5 H6 H7 H8 H9 H10 H11 H12 H13 H14 H15 H16 H17 H18 H19 H20 H21 H22 H23 H24 H25) and resid 12 to 40"]
set tot_num [$total num]

for {set i 0} {$i < $nf} {incr i} {
molinfo top set frame $i
set a [atomselect top "(segname B4 B5 B6 B7 B8 B9 B10 B11 B12 B13 B14 B15 B16 B17 B18 B19 B20 B21 B22 B23 B24 B25 H4 H5 H6 H7 H8 H9 H10 H11 H12 H13 H14 H15 H16 H17 H18 H19 H20 H21 H22 H23 H24 H25) and resid 12 to 40 and (within 10 of segname NANO N P)" frame $i]
set a_num [$a num]
set tot [expr $tot+$a_num]
$a delete
}
set avg_num [expr $tot/$nf]
puts "$avg_num"
set avg_num_out [expr $tot_num-$avg_num]
puts "$avg_num_out"
