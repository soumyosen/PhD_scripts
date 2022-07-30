
for {set i 0} {$i<=19} {incr i} {
set p1 [atomselect top "resid $i"]
set dist [expr $i*7.1]
$p1 moveby "0.0 0.0 $dist"
}
[atomselect top all] writepdb fibril_20_1.pdb
