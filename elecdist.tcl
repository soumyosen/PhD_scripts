

set micresidue [atomselect top "resname LYS and type NH3"]
set ncoord [$micresidue get {x y z}] 
set tot 0
foreach coord $ncoord {
#set sel [atomselect top "residue $residuemic and resname LYS"]
#set selc [measure center $sel]
set memresidue [lsort -unique [[atomselect top "same residue as (resname DPPC DPPG and z>-16)"] get residue]]
set num [llength $memresidue]
#set tot 0
foreach residuemem $memresidue {
set p [atomselect top "residue $residuemem and type PL O2L OSLP"]
set pmid [measure center $p]
set dist [vecdist $pmid $coord]
set energy [expr double(1/$dist)]
puts "$energy"
set tot [expr $tot+$energy]
$p delete
}
#set avg [expr $tot/$num]
#puts "$avg"
#$sel delete
}
puts "$tot"
$micresidue delete






