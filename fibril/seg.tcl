set nf [molinfo top get numframes]

for {set i 0} {$i < $nf} {incr i} {
	molinfo top set frame $i
	set atoms [atomselect top {segname B1 B2 B3 B4 B4 B6 B7 B8 B9 B10 B11 B12 B13 B14 B15 B16 B17 B18 B19 B20 B21 B22 B23 B24 B25 B26 B27 B28 B29 and (resid 12 to 40) and (within 3 of segname NANO N P)}]
	set seg [lsort -unique [$atoms get segname]]
	puts "closeset peptides $seg"
}
