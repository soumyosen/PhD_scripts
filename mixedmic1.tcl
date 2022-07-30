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

 set Sphere1 [Sphere 60 50]
 set midpt [vecscale [eval "vecadd $Sphere1"] [expr 1./60.]]
 set Sphere2 [Sphere 12 50]
#puts "total $Sphere1"
  
 set closestpt {}
 foreach pt $Sphere2 {
 set dist1 9999.9
 foreach pt1 $Sphere1 {
 set dist2 [vecdist $pt1 $pt]
 if {$dist2 < $dist1} {set dist1 $dist2; set closest [concat $pt1]}
 }
 lappend closestpt $closest
 }
# puts "2k-fa $closestpt"

set rempt {}
foreach pt2 $Sphere1 {
 if { [ lsearch $closestpt $pt2 ] < 0 } {lappend rempt $pt2} 
}
#puts "$rem"
#puts "0.6k $rempt"
#set num1 [llength $closestpt]
#set num2 [llength $rem]
#puts "2k-fa $num1"
#puts "0.6k $num2"

# puts "all $Sphere1"

#puts "points $Sphere1" 
#
#puts "midpt $midpt"
#foreach pt $Sphere1 {
#draw sphere $pt radius 2
#}
#draw sphere $midpt radius 2
#
#
# 1) 2 sphere pt
# 2) for each pt, arrange 
 
 set lig14_residue [lsort -unique [[atomselect top "segname MAIN"] get residue]]

 foreach residue $lig14_residue pt3 $rempt {
 set Nc [atomselect top "residue $residue and resname ACT DEC TRZ"]
 set Ncoord [measure center $Nc]
 set N [atomselect top "residue $residue and resname IPA and type NG2S"]
 set Nmid [measure center $N]
 set all [atomselect top "residue $residue"]
 set v1 [vecnorm [vecsub $Nmid $Ncoord]]
 set v2 [vecnorm [vecsub $pt3 $midpt]]
 set ang [expr acos(max(-1.0,min(1.0,[vecdot $v1 $v2])))]
 set cross "{[veccross $v1 $v2]}"
 $all move [eval "trans center {$Ncoord} offset {$pt3} axis $cross $ang rad"]
 $Nc delete
 $N delete
 $all delete
  }

 
 set lig44_residue [lsort -unique [[atomselect top "segname ONE"] get residue]]

 foreach residue1 $lig44_residue pt4 $closestpt {
 set Nc [atomselect top "residue $residue1 and resname ACT"]
 set Ncoord [measure center $Nc]
 set N [atomselect top "residue $residue1 and type NG2S"]
 set Nmid [measure center $N]
 set all [atomselect top "residue $residue1"]
 set v1 [vecnorm [vecsub $Nmid $Ncoord]]
 set v2 [vecnorm [vecsub $pt4 $midpt]]
 set ang [expr acos(max(-1.0,min(1.0,[vecdot $v1 $v2])))]
 set cross "{[veccross $v1 $v2]}"
 $all move [eval "trans center {$Ncoord} offset {$pt4} axis $cross $ang rad"]
 $Nc delete
 $N delete
 $all delete
  }  
 




 


