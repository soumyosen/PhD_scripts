mol load pdb fin.pdb
set a [atomselect top all]
set c [atomselect top "segname H2 H4 H6 B2 B4 B6 H28 H26 H24 B28 B26 B24 and name CA and resid 17 to 20 31 to 35"]

$a set beta 0.0 
$c set beta 1.0 

puts "CENTER:"
measure center $c

puts "Put these in the colvars file:"
$c get serial

$a writepdb ref_orient.pdb
exit
