 proc readData {filename} {
    set distance {}
    set chargediff  {}
    set f [open $filename r]
    foreach line [split [read $f] \n] {
        if {$line != ""} {
        lappend distance [lindex $line 0]
        lappend chargediff  [lindex $line 1]
        }
    }
    return [list $distance $chargediff]
}


set data [readData deltacharge-4_20.log]
puts "data $data"
set tot 0.0
foreach dist [lindex $data 0] cdiff [lindex $data 1] {
set charge_dist_ratio [expr double($cdiff)/$dist]
#puts "$charge_dist_ratio"
set tot [expr $tot+$charge_dist_ratio]
}

puts "capacitance $tot" 
