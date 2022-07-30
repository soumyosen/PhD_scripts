set nf [molinfo top get numframes]
set prot [atomselect top {segname P3}]
set prot_resid1 [atomselect top {segname P3 and resid 110}]
set prot_resid2 [atomselect top {segname P3 and resid 222}]
set prot_resid3 [atomselect top {segname P3 and resid 39}]
set prot_resid4 [atomselect top {segname P3 and resid 233}]
set part [atomselect top {segname NANO G}]

for {set i 1} {$i <= $nf} {incr i} {
 $prot frame $i
 $prot_resid1 frame $i
 $prot_resid2 frame $i
 $prot_resid3 frame $i
 $prot_resid4 frame $i
 $part frame $i
 set prot_cent [measure center $prot]
 set cent_resid1 [measure center $prot_resid1]
 set cent_resid2 [measure center $prot_resid2]
 set cent_resid3 [measure center $prot_resid3]
 set cent_resid4 [measure center $prot_resid4]
 set part_cent [measure center $part]
 
 set v_pa_pr [vecsub $prot_cent $part_cent]
 set v_pa_pr_n [vecnorm $v_pa_pr]
 set v1 [vecsub $cent_resid1 $cent_resid2]
 set v1_n [vecnorm $v1]
 set v2 [vecsub $cent_resid3 $cent_resid4]
 set v2_n [vecnorm $v2]
 set ang1 [expr acos(max(-1.0,min(1.0,[vecdot $v_pa_pr_n $v1_n])))] 
 set ang1_deg [expr $ang1*57.296]
 set ang2 [expr acos(max(-1.0,min(1.0,[vecdot $v_pa_pr_n $v2_n])))] 
 set ang2_deg [expr $ang2*57.296]
 puts "First $ang1_deg   Second $ang2_deg"
}

$prot delete
$prot_resid1 delete
$prot_resid2 delete
$prot_resid3 delete
$prot_resid4 delete
$part delete
