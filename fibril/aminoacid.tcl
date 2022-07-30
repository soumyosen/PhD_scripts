set nf [molinfo top get numframes]

for {set i 0} {$i < $nf} {incr i} {
	molinfo top set frame $i
	set atoms [atomselect top {not water and not segname ION NANO N P and (within 3 of segname N)}]
	set aa [lsort -unique [$atoms get resname]]
        set aaid [lsort -unique [$atoms get resid]]
	puts "$aa $aaid "
}
