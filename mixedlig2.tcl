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

set Sphere1 [Sphere 204 30]
set Sphere2 [Sphere 4 30]
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
#puts "peptide $closestpt"
set rempt {}
foreach pt2 $Sphere1 {
 if { [ lsearch $closestpt $pt2 ] < 0 } {lappend rempt $pt2}
}
#puts "ligand $rempt"



proc selectattachpt {zincindcoord ptcoord} {
 set Zn_indexes {}
 foreach pt $ptcoord {
 set index_of_closest_Zn -1
 set dist1 9999.9
 foreach indcoord $zincindcoord {
 set pt2 [lrange $indcoord 1 3]
 set dist2 [vecdist $pt $pt2]
 if {$dist2 < $dist1} {set dist1 $dist2;set index_of_cloest_Zn [lindex $indcoord 0]}
 }
 lappend Zn_indexes $index_of_cloest_Zn
}
 return $Zn_indexes
}

set allZn [atomselect top {type Zn1 Zn2 Zn3}]
$allZn moveby [vecscale [measure center $allZn] -1]
set Zn_xyz [$allZn get {index x y z}]
$allZn delete


set lig [selectattachpt $Zn_xyz $rempt]
set pep [selectattachpt $Zn_xyz $closestpt]

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


set ligand [[atomselect top {segname G and name H1}] get residue]
foreach residue $ligand Zn_index $lig {
set H [atomselect top "residue $residue and name H1"]
set C [atomselect top "residue $residue and name C1"]
lattach [$H get index] [$C get index] $Zn_index
set S [atomselect top "residue $residue and name S1"]
topo addbond [$S get index] $Zn_index
$S delete
$H delete
$C delete
}



set peptide [[atomselect top {segname P and name H1}] get residue]
foreach residue2 $peptide Zn_index $pep {
set H [atomselect top "residue $residue2 and name H1"]
set S [atomselect top "residue $residue2 and name S1"]
lattach [$H get index] [$S get index] $Zn_index
set N [atomselect top "residue $residue2 and name N1"]
topo addbond [$N get index] $Zn_index
$H delete
$S delete
$N delete
}

set need [atomselect top {not name H1}]
$need writepsf t2.psf
$need writepdb t2.pdb

