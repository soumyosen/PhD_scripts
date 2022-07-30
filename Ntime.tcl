
set outfile [open Ntime.plot w]
set nf [molinfo top get numframes]

# reference frame
set ref_fr 0
molinfo top set frame $ref_fr
set O1 [atomselect top "type OT and (within 4 of (resname li3 and type SG NG))" ]
set ind [$O1 get index]
set Nr [llength $ind]

# find everage number of water in hydration shell
for {set i 1} {$i < $nf} {incr i} {
  molinfo top set frame $i
  $O1 update
  set Nr [expr $Nr+[$O1 num]]
}
set Nr1 [expr $Nr/$nf]
# number of water in ref frame that are still in hydration shell
set O2 [atomselect top "index $ind and (within 4 of (resname li3 and type SG NG))"]
for {set i 0} {$i < $nf} {incr i} {
  molinfo top set frame $i
  $O2 update
  set tim [expr {4*($i)}]
  puts $outfile "$tim [expr double([$O2 num])/$Nr1]"

}

$O1 delete
$O2 delete

