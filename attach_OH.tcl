set oxy [lsort -unique [[atomselect top {segname E}] get residue]]
set hyd [atomselect top {segname D and type HH2}]
set hyd_indxyz [$hyd get {index x y z}]
set H_indexes {}
set Si_indexes {}

foreach res $oxy {
set index_of_closest_H -1
set dist1 9999.0
set atoms [atomselect top "residue $res"]
set cent_oxy [measure center [atomselect top "residue $res and type OHT"]]
set ind_oxy [[atomselect top "residue $res and type OHT"] get index]

foreach indcoord $hyd_indxyz {
set hcoord [lrange $indcoord 1 3]
set dist2 [vecdist $hcoord $cent_oxy]
if {$dist2 < $dist1} {set dist1 $dist2;set index_of_closest_H [lindex $indcoord 0]}
 }
lappend H_indexes $index_of_closest_H 
set cent_H [measure center [atomselect top "index $index_of_closest_H"]]
$atoms move [eval "trans center {$cent_oxy} offset {$cent_H}"]
set attach_ind [[atomselect top "type SIH2 and within 1.5 of index $index_of_closest_H"] get index]
puts "$ind_oxy    $index_of_closest_H    $attach_ind"
topo addbond $ind_oxy $attach_ind
lappend Si_indexes $attach_ind
$atoms delete
} 
set c [lsort -unique $H_indexes]
set c_len [llength $c]
puts "H replaced $c_len"
set d [lsort -unique $Si_indexes]
set d_len [llength $d]
puts "Si attached $d_len"
set e [llength $oxy]
puts "OH groups $e"

topo guessangles
[atomselect top "not index $H_indexes"] writepdb sih2_pil_OH5.pdb 
[atomselect top "not index $H_indexes"] writepsf sih2_pil_OH5.psf 
