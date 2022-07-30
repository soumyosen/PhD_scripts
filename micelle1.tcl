 proc Sphere {num_pts radius} {
    set PI 3.1415926535897931
    set dlong [expr $PI*(3-sqrt(5.0))]
    set dz [expr 2.0/$num_pts]
    set long 0.0
    set z [expr 1-$dz/2]
    set sphere {}
    set num {}
    for {set k 1} {$k<=$num_pts} {incr k} {
        set r [expr sqrt(1-$z*$z)]
        set sp_x [expr $radius*cos($long)*$r]
        set sp_y [expr $radius*sin($long)*$r]
        set sp_z [expr $radius*$z]
        lappend sphere "$sp_x $sp_y $sp_z"
        lappend num "$k"
        set z [expr $z-$dz]
        set long [expr $long+$dlong]
    }
    return $sphere
    return $num
}

 set Sphere1 [Sphere 60 120]
 set midpt [vecscale [eval "vecadd $Sphere1"] [expr 1./60.]]
#puts "points $Sphere1" 
#puts "midpt $midpt"
#foreach pt $Sphere1 {
#draw sphere $pt radius 2
#}
#draw sphere $midpt radius 2

 set lig_residue [lsort -unique [[atomselect top "segname aaaa"] get residue]]

 foreach residue $lig_residue pt $Sphere1 {
 set Nc [atomselect top "residue $residue and resname bz am alp rt"]
 set Ncoord [measure center $Nc]
 set N [atomselect top "residue $residue and resname THY ADE GUA CYT"]
 set Nmid [measure center $N]
 set all [atomselect top "residue $residue"]
 set v1 [vecnorm [vecsub $Nmid $Ncoord]]
 set v2 [vecnorm [vecsub $pt $midpt]]
 set ang [expr acos(max(-1.0,min(1.0,[vecdot $v1 $v2])))]
 set cross "{[veccross $v1 $v2]}"
 $all move [eval "trans center {$Ncoord} offset {$pt} axis $cross $ang rad"]
 $Nc delete
 $N delete
 $all delete
  }
 






 


