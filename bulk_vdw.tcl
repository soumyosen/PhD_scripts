# 1cal=4.184J, kunit is in N    (1 kcal/mol*A = 69.478578545 pN)
 set kunit [expr 4184/(6.022*pow(10,23)*pow(10,-10))]
 # force in N = $k*$kunit

 set forcesRecalcFreq  1000
 set cutoff 80.
 set numNP 72

set centers1 {1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 122 123 124 125 138 139 140 141 154 155 156 157 170 171 172 173 186 187 188 189 202 203 204 205 210 211 212 213 214 215 220 221 222 223 224 225 230 231 232 233 234 235 240 241 242 243 244 245 250 251}
 
set centers2 {5937 5938 5939 5940 5941 5942 5943 5944 5945 5946 5947 5948 5949 5950 5951 5952 5953 5954 5955 5956 5957 6058 6059 6060 6061 6074 6075 6076 6077 6090 6091 6092 6093 6106 6107 6108 6109 6122 6123 6124 6125 6138 6139 6140 6141 6146 6147 6148 6149 6150 6151 6156 6157 6158 6159 6160 6161 6166 6167 6168 6169 6170 6171 6176 6177 6178 6179 6180 6181 6186 6187}

set centers3 {11873 11874 11875 11876 11877 11878 11879 11880 11881 11882 11883 11884 11885 11886 11887 11888 11889 11890 11891 11892 11893 11994 11995 11996 11997 12010 12011 12012 12013 12026 12027 12028 12029 12042 12043 12044 12045 12058 12059 12060 12061 12074 12075 12076 12077 12082 12083 12084 12085 12086 12087 12092 12093 12094 12095 12096 12097 12102 12103 12104 12105 12106 12107 12112 12113 12114 12115 12116 12117 12122 12123}

set hamaker [expr 6.022E23*(3/(4*3.14159)*1.95*1.602E-19)/4184] ; # kcal/mol
set r 14.0

foreach serial1 $centers1 serial2 $centers2 serial3 $centers3 {
lappend savedforces {0. 0. 0.}
}

set counter [expr $forcesRecalcFreq-1] 

# calcforce is executed every step
proc calcforces { } {
global centers1 centers2 centers3
global hamaker r cutoff
global forcesRecalcFreq counter savedforces

# request atom coord for the next force evaluation
 if { $counter == [expr $forcesRecalcFreq-1] } {
 foreach serial1 $centers1 { addatom $serial1 }
 foreach serial2 $centers2 { addatom $serial2 }
 foreach serial3 $centers3 { addatom $serial3 }
}  



# force calculation
if { $counter == $forcesRecalcFreq } {
loadcoords c
set cent1 {0.0 0.0 0.0}
set cent2 {0.0 0.0 0.0}
set cent3 {0.0 0.0 0.0}
foreach ind1 $centers1 ind2 $centers2 ind3 $centers3 {
set c1 [vecscale $c($ind1) [expr 1/71]]
set cent1 [vecadd $cent1 $c1]
set c2 [vecscale $c($ind2) [expr 1/71]]
set cent2 [vecadd $cent2 $c2]
set c3 [vecscale $c($ind3) [expr 1/71]]
set cent3 [vecadd $cent3 $c3]
}

set lst [list $cent1 $cent2 $cent3]


set savedforces {}
foreach pt1 $lst {
set forces_from_neighbors {}
foreach pt2 $lst {
if {$pt1 == $pt2} {continue} ; # np itself
set vd [vecsub $pt2 $pt1]
set d [veclength $vd]
if {$d > $cutoff} {continue}
set fa [expr 1./(4.*$r*$d+$d*$d)] ; # 1/a
set fb [expr 1./(4.*$r*$r+4.*$r*$d+$d*$d)] ; # 1/b
set force [expr ($hamaker/3.)*(2.*$r+$d)*($fa-$fb-2*$r*$r*($fa*$fa+$fb*$fb))]
lappend forces_from_neighbors [vecscale $vd [expr $force/$d]] ; # $d is for normalizing
}
set applyforce {0. 0. 0.} ; # reset
foreach force $forces_from_neighbors {
set applyforce [vecadd $applyforce $force]
}
lappend savedforces $applyforce
}
set counter 0
clearconfig
}

# apply force to atoms
foreach serial1 $centers1 {
set f1 [lindex $savedforces 0]
addforce $serial1 $f1
}

foreach serial2 $centers2 {
set f2 [lindex $savedforces 1]
addforce $serial2 $f2
}

foreach serial3 $centers3 {
set f3 [lindex $savedforces 2]
addforce $serial3 $f3
}

incr counter
return
}









































