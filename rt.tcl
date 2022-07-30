menu main on
mol load psf NW1-3.psf dcd t4.dcd
set outfile [open NP3rt.xvg w]
set sel [atomselect top "same residue as (water and within 2.4 of resname li3 and type SG NG)"]
set nf [molinfo top get numframes]
set m 0
for {set i 0} {$i < $nf} {incr i} {
  $sel frame $i
  $sel update
  set residuelist [$sel get residue]
  set num [llength $residuelist]
  set num1 [expr $num/3]
  set m [expr $m+$num1]
  set tim [expr {0.002*($i)}]
  puts $outfile "$tim $residuelist $num1"

}
 set avg [expr $m/$nf]
 puts $outfile "avg $avg"
exit
