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

set Sphere [Sphere 580 27]
set lig1pt {}
for {set i 1} {$i < 500} {incr i} {
set pt [lindex $Sphere [expr {int(rand()*527)}]]
if { [ lsearch $lig1pt $pt ] < 0 } {lappend lig1pt $pt}
set num [llength $lig1pt]
if {$num == 290} {
 break
 }
}

set lig2pt {}
foreach pt2 $Sphere {
 if { [ lsearch $lig1pt $pt2 ] < 0 } {lappend lig2pt $pt2}
}

#puts "lig1pt $lig1pt"
#puts "lig2pt $lig2pt" 

#set coord1 [llength $lig1pt]
#puts "coord1 $coord1"
#set coord2 [llength $lig2pt]
#puts "coord2 $coord2"


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

set goldlig1 [selectattachpt $Au_xyz $lig1pt]
set goldlig2 [selectattachpt $Au_xyz $lig2pt]

#puts "goldlig1:   $goldlig1"
#puts "goldlig2:   $goldlig2"

#set num1 [llength $goldlig1]
#set num2 [llength $goldlig2]
#puts "num1 $num1"
#puts "num2 $num2"


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

set lig_residue1 [[atomselect top {segname A and type TYPE}] get residue]
foreach residue $lig_residue1 Au_index $goldlig1 {
set H [atomselect top "residue $residue and type TYPE"]
set C [atomselect top "residue $residue and type CG33"]
lattach [$H get index] [$C get index] $Au_index 
set S [atomselect top "residue $residue and name S"]
topo addbond [$S get index] $Au_index
$S delete
$H delete
$C delete
}


set lig_residue2 [[atomselect top {segname S and type TYPE}] get residue]
foreach residue2 $lig_residue2 Au_index $goldlig2 {
set H [atomselect top "residue $residue2 and type TYPE"]
set S1 [atomselect top "residue $residue2 and type SG3O"]
lattach [$H get index] [$S1 get index] $Au_index
set S [atomselect top "residue $residue2 and type SG31"]
topo addbond [$S get index] $Au_index
$S delete
$H delete
$S1 delete
}

set need [atomselect top {not type TYPE}]
$need writepsf t1.psf
$need writepdb t1.pdb

