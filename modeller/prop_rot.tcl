
set first [measure center [atomselect top {index 528}]]
set second [measure center [atomselect top {index 1126}]]
set third [measure center [atomselect top {index 1097}]]
set th_fi [vecsub $third $first]
set th_fi_norm [vecnorm $th_fi]
set th_se [vecsub $third $second]
set th_se_norm [vecnorm $th_se]
set cross [veccross $th_fi_norm $th_se_norm]
set midpt [measure center [atomselect top {index 528 1126}]]

for {set i 0} {$i<=19} {incr i} {
set p1 [atomselect top "resid $i"]
set dist [expr $i*7.1]
$p1 moveby "0.0 0.0 $dist"
set angle [expr $i*1.29]
$p1 move [eval "trans center {$midpt} axis {$cross} $angle deg"]
}

[atomselect top all] writepdb fibril_20_2.pdb
