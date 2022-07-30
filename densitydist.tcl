 for {set i 1} {$i < 51} {incr i} {
 set dentot 0 
 set nf [molinfo top get numframes] 
 for {set j 1} {$j<=$nf} {incr j} {
 molinfo top set frame $j
 set pi 3.14
 set vol [expr 4*$pi*(pow($i,3)-pow(($i-1),3))/3]
 set fa [atomselect top "resname FA"]
 set facent [measure center $fa]
 set x [lindex $facent 0]
 set y [lindex $facent 1]
 set z [lindex $facent 2]
 set sel_string [string map {-- +} "resname DEE MEE IPA and (sqr(x-$x)+sqr(y-$y)+sqr(z-$z)<sqr($i)) and (not sqr(x-$x)+sqr(y-$y)+sqr(z-$z)<sqr([expr $i-1]))"]
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
 set dentot [expr $dentot+$density]
$fa delete
$sel delete
}
 set denavg [expr $dentot/$nf]
 puts "$i  $denavg"
}
