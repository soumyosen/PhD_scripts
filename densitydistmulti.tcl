 for {set i 1} {$i < 50} {incr i} {
 set dentot 0 
 set nf [molinfo top get numframes] 
 for {set j 1} {$j<=$nf} {incr j} {
 set pi 3.14
 set vol [expr 4*$pi*(pow($i,3)-pow(($i-1),3))/3]
 set fa [atomselect top "resname FA"]
 set fa_residue [lsort -unique [$fa get residue]]
 set nfa [llength $fa_residue]
 set denmultifa 0
 foreach residue $fa_residue {
 set atom [atomselect top "resname FA and residue $residue"]
 set facent [measure center $atom]
 set x [lindex $facent 0]
 set y [lindex $facent 1]
 set z [lindex $facent 2]
 set sel_string [string map {-- +} "segname MAIN ONE and not resname FA and (sqr(x-$x)+sqr(y-$y)+sqr(z-$z)<sqr($i)) and (not sqr(x-$x)+sqr(y-$y)+sqr(z-$z)<sqr([expr $i-1]))"]
 set sel [atomselect top $sel_string]
 set selind [$sel get index]
 set m 0
 foreach ind $selind {
 set M [atomselect top "index $ind"]
 set mass [$M get mass]
 set m [expr $m+$mass]
$M delete
} 
 set density [expr $m/$vol]
 set denmultifa [expr $denmultifa+$density]
$sel delete
$atom delete
}
 set denmultifaavg [expr $denmultifa/$nfa]
 set dentot [expr $dentot+$denmultifaavg]
$fa delete
}
 set denavg [expr $dentot/$nf]
 puts "$i  $denavg"
}
