menu main on
mol load psf test_amber.psf dcd equil-e3.dcd
mol addfile equil-e4.dcd type dcd waitfor all
mol addfile equil-e5.dcd type dcd waitfor all
mol addfile equil-e6.dcd type dcd waitfor all

 set outfile [open sasa-prot-dsrna-rbd.xvg w] 
 set outfile2 [open sasa-rna-dsrna-rbd.xvg w] 
 set outfile3 [open sasa-all-dsrna-rbd.xvg w] 
 set outfile4 [open sasa-diff-dsrna-rbd.xvg w] 

 set nf [molinfo top get numframes]
 set p [atomselect top "protein"]
 set n [atomselect top "nucleic"]
 set np [atomselect top "nucleic or protein"]

# for {set i 0} {$i < $nf} {incr i} {
 for {set i 0} {$i < $nf} {incr i} {
      $p frame $i
      $n frame $i
      $np frame $i
   set sasap [measure sasa 1.4 $p -restrict $p]
   set sasan [measure sasa 1.4 $n -restrict $n]
   set sasanp [measure sasa 1.4 $np -restrict $np]
   set tot [expr {0.5*($sasap+$sasan-$sasanp)}] 
   set tim [expr {0.002*($i)}]
   puts $outfile "$tim $sasap"
   puts $outfile2 "$tim $sasan"
   puts $outfile3 "$tim $sasanp"
   puts $outfile4 "$tim $tot"
 } 

exit
