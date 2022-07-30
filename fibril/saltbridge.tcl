set nf [molinfo top get numframes]
set tot_contact 0.0

for {set i 0} {$i < $nf} {incr i} {
molinfo top set frame $i

set contact 0

set atoms_ind [[atomselect top "segname B4 B5 B6 B7 B8 B9 B10 B11 B12 B13 B14 B15 B16 B17 B18 B19 B20 B21 B22 B23 B24 B25 H4 H5 H6 H7 H8 H9 H10 H11 H12 H13 H14 H15 H16 H17 H18 H19 H20 H21 H22 H23 H24 H25 and resid 12 to 40 and resname ARG LYS and name NH1 NH2 NZ"] get index]

foreach index $atoms_ind {
set neighbor_atoms [[atomselect top "segname B4 B5 B6 B7 B8 B9 B10 B11 B12 B13 B14 B15 B16 B17 B18 B19 B20 B21 B22 B23 B24 B25 H4 H5 H6 H7 H8 H9 H10 H11 H12 H13 H14 H15 H16 H17 H18 H19 H20 H21 H22 H23 H24 H25 and resid 12 to 40 and resname ASP GLU and name OD1 OD2 OE1 OE2 and (within 4 of index $index)"] num]
set contact [expr $contact+$neighbor_atoms]
}

puts "$i  $contact"
set tot_contact [expr $tot_contact+$contact]
}

set avg_contact [expr $tot_contact/$nf]
puts "$avg_contact"
