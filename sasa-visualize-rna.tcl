menu main on
mol load psf test_amber.psf dcd equil-e2.dcd

 set outfile [open sasa-prot-dsrna-rbd.xvg w] 
 set outfile2 [open sasa-rna-dsrna-rbd.xvg w] 
 set outfile3 [open sasa-all-dsrna-rbd.xvg w] 
 set outfile4 [open sasa-diff-dsrna-rbd.xvg w] 

 set nf [molinfo top get numframes]
 set all [atomselect top all]
 set p [atomselect top "protein"]
 set n [atomselect top "nucleic"]
 set np [atomselect top "nucleic or protein"]

# for {set i 0} {$i < $nf} {incr i} {
 for {set i 0} {$i < 1} {incr i} {
      $p frame $i
      $n frame $i
      $np frame $i
   set sasap [measure sasa 1.4 $p -restrict $p]
   set sasan [measure sasa 1.4 $n -points pts -restrict $n]
   set sasanp [measure sasa 1.4 $np -restrict $np]
   set tot [expr {0.5*($sasap+$sasan-$sasanp)}] 
   puts $outfile "$i $sasap"
   puts $outfile2 "$i $sasan"
   puts $outfile3 "$i $sasanp"
   puts $outfile4 "$i $tot"
    foreach pt $pts {
    #draw point $pt radius 0.5
     graphics top sphere $pt radius 0.4
    }

 } 

#exit
