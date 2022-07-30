for {set i 2} {$i <= 92} {incr i 3} {
   set a [expr $i+1]
   set b [expr $i+2]
   topo addbond "$i" "$a"
   topo addbond "$i" "$b"
}

topo guessangles
topo guessdihedrals
[atomselect top all] writepsf b.psf


