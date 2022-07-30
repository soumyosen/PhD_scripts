set fp [open "twistangleB2_1" r]
set file_data [read $fp]
close $fp

set data [split $file_data "\n"]
set count 0

foreach line $data {
	#puts "$line"
	if {$line >= 20.0} {
		set count [expr $count+1]
	}
}
puts "Total number $count"
set probability [expr $count/10500.0]
puts "probability over 20 or lower than 10 is $probability"
