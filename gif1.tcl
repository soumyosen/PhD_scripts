proc rotate_axis {vec deg {molid top}} {
    ## example: rotate_axis {1 1 0} 20
    # get the current matrix
    lassign [molinfo $molid get rotate_matrix] curr
    # the transformation matrix
    set r [trans axis $vec $deg]
    # get the new matrix
    set m [transmult $r $curr]
    # and apply it to the molecule
    molinfo $molid set rotate_matrix "{ $m }"
}

set frame 0
for {set fr 0} {$fr < [molinfo top get numframes]} {incr fr 2} {
#for {set fr 0} {$fr < 100} {incr fr 1} {}
molinfo top set frame $fr

#    rotate y by 0.5
    incr frame
    set filename [format "%04d" $frame].ppm
    render snapshot $filename
}


exec convert -quality 100 *.ppm movie.mpeg

