
proc lookup_electron {} {
global atom_name atom_e
set atom_e 0.0
# remap
if { "$atom_name" == "Cla" } { set atom_name "Cl" }
if { "$atom_name" == "Sod" } { set atom_name "Na" }
if { [lsearch {Oh2 O12 O13 O21 O22} $atom_name] != -1 } { set atom_name "O" }
if { [lsearch {C11 C12 C13 C14 C15 C1 C2 C3 C21 C29 C32 C33 C218} $atom_name] != -1 } { set atom_name "C" }
if { [lsearch {H1 H2 Hs Ha H91 H18r H12a} $atom_name] != -1 } { set atom_name "H" }

# periodic table
set atom_e [expr double([lsearch {H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Uut Fl Uup Lv Uus Uuo} $atom_name]+1)]
}

proc createEmptyBin {} {
  global z_min z_step z_bin bins
  for {set i 1} {$i <= [expr $z_bin]} {incr i} {
  dict set bins b${i} bin_bottom [expr $z_min+($i-1)*$z_step] ; # lower bound
  dict set bins b${i} bin_top [expr $z_min+$i*$z_step] ; # upper bound
  dict set bins b${i} charge 0.0
  dict set bins b${i} electron 0.0
  puts "[format "%14.4f%14.4f%14.4f" $i [dict get $bins b$i bin_bottom] [dict get $bins b$i bin_top]]"
  }
}

proc sliceAtom {} {
 # Assign fraction of charge/electron of atom slices to bin
 while 1 {
  # 1st cap, spherical cap formula (removed factor PI/3)
  set z_cut [dict get $bins b$use_bin bin_bottom]
  if { $z_cut < $atom_bottom_z} { set z_cut $atom_bottom_z }
  set h [expr $atom_top_z-$z_cut]
  set v [expr ($h**2)*(3*$atom_radius-$h)]
  # 2nd cap, spherical cap formula (removed factor PI/3)
  set z_lastcut [dict get $bins b$use_bin bin_top]
  if { $z_lastcut > $atom_top_z} { set z_lastcut $atom_top_z }
  set h [expr $atom_top_z-$z_lastcut]
  set v_remove [expr ($h**2)*(3*$atom_radius-$h)]
  # volume ratio = (1st cap - 2nd cap)/atom_volume
  set ratio [expr ($v-$v_remove)/(4*$atom_radius**3)]
  # assign charge to bin (based on volume ratio)
  dict set bins b$use_bin charge [expr [dict get $bins b$use_bin charge]+$atom_charge*$ratio]
  # assign electron to bin (based on volume ratio)
  dict set bins b$use_bin electron [expr [dict get $bins b$use_bin electron]+$atom_e*$ratio]
  # proceed to next atom slice (the bin below the current bin)
  incr use_bin -1
  #  Out of range check (bottom part of the atom)
  if { $use_bin == 0 } {break}
  if { [dict get $bins b$use_bin bin_top] <= $atom_bottom_z } {break}
 } ; # done with one atom
}

package require pbctools

# SET THESE
set numFrame [molinfo top get numframes]
set startFr [expr [molinfo top get numframes]-$numFrame]
set selection [atomselect top all]
set disable_slicing 1 ; # Disable slicing, set value = 1 to enable
set z_min [expr [lindex [measure minmax $selection] 0 2]-10] ; # lowest z to plot
set z_max [expr [lindex [measure minmax $selection] 1 2]+10] ; # highest z to plot
set z_bin 100    ; # how many bins

# parameters
set start_time [clock clicks -milliseconds]
set next_time $start_time 
set x_width [lindex [pbc get] 0 0] 
set y_width [lindex [pbc get] 0 1] 
set z_min [expr double($z_min)]
set z_max [expr double($z_max)]
set z_step [expr ($z_max-$z_min)/$z_bin] ; # width of each bin
createEmptyBin

# Load force field
if { $disable_slicing != 1 } {
package require ilstools
ILStools::readcharmmparams {/home/irena/Desktop/ions/chgNP.inp /home/irena/Desktop/ions/force.inp}
ILStools::assigncharmmparams top
}

# check atom name
puts "Assigning atom_e by atom name"
foreach atom_name [$selection get {name}] {
set atom_name [string totitle $atom_name]
lookup_electron
if { $atom_e == 0 } { puts "need to define number of electron for name $atom_name";vwait forever }
}



for {set fr 0} {$fr < $numFrame} {incr fr} {
set use_fr [expr $fr+$startFr]; molinfo top set frame $use_fr

# Loop through all atoms
set count 0; set allatom [$selection num]
foreach allinfo [$selection get {z name charge radius}] {

 # timer and progress messages
 incr count 1; set time [clock clicks -milliseconds]
 if {$count == $allatom || $fr == $numFrame || $time >= $next_time} then {
 puts "[expr ($time-$start_time)/1000] secs, Frame [expr $use_fr+1] ([expr $numFrame-$fr-1] more to go): [expr $allatom-$count] atoms left"
 set next_time [expr $time + 10000]
 }

 # Exclude the following atoms
 set atom_name [string totitle [lindex $allinfo 1]]
 if { "$atom_name" == "Au" } { continue }

 # Obtain atom info
 set atom_z [lindex $allinfo 0]
 set atom_charge [lindex $allinfo 2]
 lookup_electron

 # more accurate atom_e
 set atom_e [expr $atom_e-$atom_charge]
 if { $atom_e == 0.0 } {continue}
 
 # find z-Boundary of atom and check if whole atom is in the defined z-range
 if { $disable_slicing == 1 } {
 set atom_top_z $atom_z 
 set atom_bottom_z $atom_z 
 } else {
 set atom_radius [lindex $allinfo 3]
 set atom_top_z [expr $atom_z+$atom_radius]
 set atom_bottom_z [expr $atom_z-$atom_radius]
 }
 if { $atom_top_z <= $z_min } {continue}
 if { $atom_bottom_z >= $z_max } {continue}

 # Search for topmost bin containing the atom and in z-range
 set use_bin [expr int(floor(($atom_top_z-$z_min)/$z_step)+1)]
 if { [expr fmod($atom_top_z-$z_min,$z_step)] == 0 } { incr use_bin -1 } ; # boundary case
 while { $use_bin > $z_bin } { incr use_bin -1 } ; # topmost bino= in z-range

 if { $disable_slicing == 1 } {
 # assign charge to bin
 dict set bins b$use_bin charge [expr [dict get $bins b$use_bin charge]+$atom_charge]
 # assign electron to bin
 dict set bins b$use_bin electron [expr [dict get $bins b$use_bin electron]+$atom_e]
 } else {
 sliceAtom
 }

} ; # done with all atoms

} ; # done with all frames

# electron density in unit "electron charge per cubic A"
set electron_list {}
for {set i 1} {$i <= $z_bin} {incr i} {
lappend electron_list [expr [dict get $bins b$i electron]/($numFrame*$x_width*$y_width*$z_step)]
}

# rho = charge density in unit "electron charge per cubic A"
set rho_list {}
for {set i 1} {$i <= $z_bin} {incr i} {
lappend rho_list [expr [dict get $bins b$i charge]/($numFrame*$x_width*$y_width*$z_step)]
}

# z (mid pt of each bin)
set z_list {}
for {set i 1} {$i <= $z_bin} {incr i} {
lappend z_list [expr ([dict get $bins b$i bin_top]+[dict get $bins b$i bin_bottom])/2]
}

# output electron and charge density
set output [open "output_elec_rho2.txt" "w"]
for {set i 0} {$i < $z_bin} {incr i} {
puts $output "[lindex $z_list $i] [lindex $electron_list $i] [lindex $rho_list $i]"
}
close $output


## find E by integrating rho, unit "electron charge per square A"
#set E_list {} ; set total 0.0
#for {set i 1} {$i < $z_bin} {incr i} {
#set total [expr $total+0.5*$z_step*([lindex $rho_list $i]+[lindex $rho_list [expr $i-1]])]
#lappend E_list $total
#}

set E_list {} ; set total 0.0
for {set i 1} {$i < $z_bin} {incr i} {
set total_bottom 0.0;set total_top 0.0
for {set j 1} {$j < $i} {incr j} {
set total_bottom [expr $total_bottom+0.5*$z_step*([lindex $rho_list $i]+[lindex $rho_list [expr $i-1]])]
}
for {set j $i} {$j < $z_bin} {incr j} {
set total_top [expr $total_top+0.5*$z_step*([lindex $rho_list $i]+[lindex $rho_list [expr $i-1]])]
}
set total [expr $total_bottom-$total_top]
lappend E_list $total
}

# find V by integrating E, unit "electron charge per A"
set V_list {} ; set total 0.0
for {set i 1} {$i < $z_bin} {incr i} {
set total [expr $total+0.5*$z_step*([lindex $E_list $i]+[lindex $E_list [expr $i-1]])]
lappend V_list $total
}

# z (mid pt of mid pt)
set zz_list {}
for {set i 1} {$i < $z_bin} {incr i} {
lappend zz_list [expr ([lindex $z_list $i]+[lindex $z_list [expr $i-1]])/2]
}

# output
set output [open "output_E_V2.txt" "w"]
set eps 180.9313492; # eps is electon charge divided by epsilon (V/A^2)
for {set i 0} {$i < [expr $z_bin-1]} {incr i} {
puts $output "[lindex $zz_list $i] [expr $eps*[lindex $E_list $i]] [expr -$eps*[lindex $V_list $i]]"
}
close $output
