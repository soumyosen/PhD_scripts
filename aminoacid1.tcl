set close_aa [lsort -unique [[atomselect top "(segname P1) and (within 10 of segname G)"] get resid]]

set residues {}
foreach residue $close_aa {
set name [lsort -unique [[atomselect top "resid $residue"] get resname]]
lappend residues $name
}

puts "$residues"


proc countwords {str} {
    foreach word [split $str] {incr count($word)}
    foreach word [array names count] {puts "$word $count($word)"}
}


countwords $residues
