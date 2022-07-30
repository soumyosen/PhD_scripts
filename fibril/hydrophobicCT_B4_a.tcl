set nf [molinfo top get numframes]
set tot_contact 0.0

for {set i 0} {$i < $nf} {incr i} {
molinfo top set frame $i

set contact 0

for {set j 4} {$j < 26} {incr j} {
set atoms_ind [[atomselect top "segname B$j and resid 12 to 40 and noh"] get index]
set k [expr $j+1]
set l [expr $j+2]
set m [expr $j+3]
set n [expr $j+4]
set o [expr $j+5]


foreach index $atoms_ind {
set neighbor_atoms [[atomselect top "segname B$k B$l B$m B$n B$o and resid 12 to 40 and noh and (within 4 of index $index)"] num]
set contact [expr $contact+$neighbor_atoms]
}
}
puts "$i  $contact"
set tot_contact [expr $tot_contact+$contact]
}

set avg_contact [expr $tot_contact/$nf]
puts "$avg_contact"
