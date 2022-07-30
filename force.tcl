# 1cal=4.184J, kunit is in N    (1 kcal/mol*A = 69.478578545 pN)
set kunit [expr 4184/(6.022*pow(10,23)*pow(10,-10))]
# force in N = $k*$kunit

# add vertice of NP (Au3)
foreach index $X1 {
addatom [expr $index+1] 
lappend X1_list [expr $index+1]
}

foreach index $X2 {
addatom [expr $index+1] 
lappend X2_list [expr $index+1]
}

proc calcforces {} {
  global X1_list X2_list

  foreach index1 $X1_list {
    set k 0.2
    set force1 [vecscale $k {-1 0 0}]
    addforce $index1 $force1
  }

  foreach index2 $X2_list {
    set k 0.2
    set force2 [vecscale $k {1 0 0}]
    addforce $index2 $force2
  }
}

