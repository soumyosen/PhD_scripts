proc Sphere {num_pts radius} {
    set PI 3.1415926535897931
    set dlong [expr $PI*(3-sqrt(5.0))]
    set dz [expr 2.0/$num_pts]
    set long 0.0
    set z [expr 1-$dz/2]
    set sphere {}
    for {set k 1} {$k<=$num_pts} {incr k} {
        set r [expr sqrt(1-$z*$z)]
        set sp_x [expr $radius*cos($long)*$r]
        set sp_y [expr $radius*sin($long)*$r]
        set sp_z [expr $radius*$z]
        lappend sphere "$sp_x $sp_y $sp_z"
        set z [expr $z-$dz]
        set long [expr $long+$dlong]
    }
    return $sphere
}

set Sphere1 [Sphere 90 12]
#puts "$Sphere1"

set pt_list1 {}
set pt_list2 {}
foreach pt $Sphere1 {
set z [lindex $pt 2]
if {$z > 0} {
set pt_list1 [lappend pt_list1 $pt]
} else {
set pt_list2 [lappend pt_list2 $pt]
}
}
#puts "$pt_list1"
#puts "$pt_list2"
#puts "[llength $pt_list1]"
#puts "[llength $pt_list2]" 


proc selectattachpt {goldindcoord ptcoord} {
 set Au_indexes {}
 foreach pt $ptcoord {
 set index_of_closest_Au -1
 set dist1 9999.9
 foreach indcoord $goldindcoord {
 set pt2 [lrange $indcoord 1 3]
 set dist2 [vecdist $pt $pt2]
 if {$dist2 < $dist1} {set dist1 $dist2;set index_of_cloest_Au [lindex $indcoord 0]}
 }
 lappend Au_indexes $index_of_cloest_Au
}
 return $Au_indexes
}


set allAu [atomselect top {type Au1 Au2 Au3}]
$allAu moveby [vecscale [measure center $allAu] -1]
set Au_xyz [$allAu get {index x y z}]
$allAu delete


set goldlig1 [selectattachpt $Au_xyz $pt_list1]
set goldlig2 [selectattachpt $Au_xyz $pt_list2]


proc lattach {l_attachIndex l_alignIndex np_index} {
 set sel3 "(index $l_alignIndex)" ;
 set sel1 "(same residue as $sel3)" ;
 set sel2 "(index $l_attachIndex)" ;
 set sel4 "(index $np_index)" ;
 foreach i {1 2 3 4} {set sel$i [atomselect top "[set sel$i]"]}
 foreach i {2 3 4} {set xyz$i "{[measure center [set sel$i]]}"}
 set v1 [vecnorm [eval "vecsub $xyz3 $xyz2"]] ;
 set v2 [vecnorm [eval "vecsub $xyz4 {0 0 0}"]] ;
 set ang [expr acos(max(-1.0,min(1.0,[vecdot [vecnorm $v1] [vecnorm $v2]])))]
 set cross "{[veccross $v1 $v2]}"
 $sel1 move [eval "trans center $xyz2 offset $xyz4 axis $cross $ang rad"]
 $sel1 set segname [$sel4 get segname]
 foreach i {2 3 4} {set xyz$i "{[measure center [set sel$i]]}"}
 $sel1 move [eval "trans bond $xyz4 $xyz3 [expr rand()*360] deg"]
 foreach i {1 2 3 4} {[set sel$i] delete}
 return
}


set lig_residue1 [[atomselect top {segname N and type TYPE}] get residue]
foreach residue $lig_residue1 Au_index $goldlig1 {
set H [atomselect top "residue $residue and type TYPE"]
set O [atomselect top "residue $residue and type SG3O1"]
lattach [$H get index] [$O get index] $Au_index
set S [atomselect top "residue $residue and type SG311"]
topo addbond [$S get index] $Au_index
$S delete
$H delete
$O delete
}

puts "thanks"
set lig_residue2 [[atomselect top {segname P and type TYPE}] get residue]
foreach residue2 $lig_residue2 Au_index $goldlig2 {
set H [atomselect top "residue $residue2 and type TYPE"]
set N [atomselect top "residue $residue2 and type NG3P3"]
lattach [$H get index] [$N get index] $Au_index
set S [atomselect top "residue $residue2 and type SG311"]
topo addbond [$S get index] $Au_index
$S delete
$H delete
$N delete
}


set need [atomselect top {not type TYPE}]
$need writepsf np1.psf
$need writepdb np1.pdb

