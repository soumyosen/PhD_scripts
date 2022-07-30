#!/usr/bin/tcsh
#Qo profile

set startFr 2500;# start from what frame  #2500 for 1np
set numFrame 3750;#[molinfo top get numframes]  #for how many frames

set step 0.5;#step for histogram
set range1 0;#lower range
set range2 140;#upper range
set filename "dist.dat"
set fileId [open $filename "w"]


#create blank histogram
set histogram {} ; # total height
set histogram2 {}; # occurence
for {set i 0} {$i < ($range2-$range1)/$step} {incr i} {
    lappend histogram 0
    lappend histogram2 0
}

#go through all the frames
for {set fr 0} {$fr< $numFrame} {incr fr} {
    puts "Processing Fr [expr $fr+$startFr] ([expr 100*$fr/$numFrame]% overall)"
    molinfo top set frame $fr+$startFr
    set e1 [atomselect top "name PT and index 9 190 371 552 733 914 1095 1276 1457 1638 1819 2000 2181"];# select center of NP (up to 13NP)
    set e1_num [$e1 num];#num of atom selected
    set xyz_np [$e1 get {x y z}];# extract NP coordinates
    #initiate variables
    set x1 0.0
    set y1 0.0
    set z1 0.0
    #find centroid of cluster (x1,y1,z1)
    for {set i 0} {$i < $e1_num} {incr i} {
    set x1 [expr $x1+[lindex [lindex $xyz_np $i] 0]]
    set y1 [expr $y1+[lindex [lindex $xyz_np $i] 1]]
    set z1 [expr $z1+[lindex [lindex $xyz_np $i] 2]]
    }
    set x1 [expr $x1/$e1_num]
    set y1 [expr $y1/$e1_num]
    set z1 [expr $z1/$e1_num]
    #select lipid polar head Qo base of z-coor of centroid
    #set e2 [atomselect top "type Qo"]
    set e2 [atomselect top "type Qo and z>$z1"]
    #set e2 [atomselect top "(type Qo and within 30 of name PT) and z > $z1"];# top
    #set e3 [atomselect top "type Qo and z<$z1"];# bottom
    #find distances (centroid z-plane to Qo z-coor)
    for {set i 0} {$i < [$e2 num]} {incr i} {
        set xyz [$e2 get {x y z}];# extract NP coordinates
	#obtain top Qo coordinate
	set x2 [lindex [lindex $xyz $i] 0]
	set y2 [lindex [lindex $xyz $i] 1]
	set z2 [lindex [lindex $xyz $i] 2]
	#find r distance
	set dist [expr sqrt(pow(($x1-$x2),2)+pow(($y1-$y2),2))]
	#puts "$x1 $y1  $x2 $y2 $dist $i"
	#mark on histogram
	set index [expr {round($dist/$step)}];# .5-.9 round up, .1-.4 round down
        set histogram [lreplace $histogram $index $index [expr [lindex $histogram $index]+[expr abs($z2-$z1)]]]
	set histogram2 [lreplace $histogram2 $index $index [expr [lindex $histogram2 $index]+1]]
     }
}

# write histogram to file
for {set i 0} {$i < [llength $histogram]} {incr i} {
  if { [lindex $histogram2 $i] == 0} {
    set histogram2 [lreplace $histogram2 $i $i 1];#fix divide by 0 problem
  }
  puts $fileId "[expr $i*$step] [expr [lindex $histogram $i]/([lindex $histogram2 $i])]"
}

close $fileId
puts "Completed: 100% done"
