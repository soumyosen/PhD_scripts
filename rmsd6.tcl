
#mol load psf last30ns.psf pdb last50_30ns.dcd
#mol load psf last_onlyfib.psf pdb last_onlyfib.pdb

#set sel1 [atomselect 0 "(segname B1 B2 B3 B4 B5 B6 B7 B8 B9 B10 B11 B12 B13 B14 B15 B16 B17 B18 B19 B20 B21 B22 B23 B24 B25 B26 B27 B28 B29 H1 H2 H3 H4 H5 H6 H7 H8 H9 H10 H11 H12 H13 H14 H15 H16 H17 H18 H19 H20 H21 H22 H23 H24 H25 H26 H27 H28 H29) and resid 12 to 40"]
set sel2 [atomselect 1 "(segname B1 B2 B3 B4 B5 B6 B7 B8 B9 B10 B11 B12 B13 B14 B15 B16 B17 B18 B19 B20 B21 B22 B23 B24 B25 B26 B27 B28 B29 H1 H2 H3 H4 H5 H6 H7 H8 H9 H10 H11 H12 H13 H14 H15 H16 H17 H18 H19 H20 H21 H22 H23 H24 H25 H26 H27 H28 H29) and resid 12 to 40"]

set nf [molinfo 0 get numframes]

for {set i 1} {$i <= $nf} {incr i} {
molinfo top set frame $i

set sel1 [atomselect 0 "(segname B1 B2 B3 B4 B5 B6 B7 B8 B9 B10 B11 B12 B13 B14 B15 B16 B17 B18 B19 B20 B21 B22 B23 B24 B25 B26 B27 B28 B29 H1 H2 H3 H4 H5 H6 H7 H8 H9 H10 H11 H12 H13 H14 H15 H16 H17 H18 H19 H20 H21 H22 H23 H24 H25 H26 H27 H28 H29) and resid 12 to 40" frame $i]
set transform_matrix [measure fit $sel1 $sel2]

set all [atomselect 0 all frame $i]
$all move $transform_matrix

set rms [measure rmsd $sel1 $sel2]
puts "RMSD of beta sheet in frame $i is $rms"
}
