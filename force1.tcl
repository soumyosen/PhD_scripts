# 1cal=4.184J, kunit is in N    (1 kcal/mol*A = 69.478578545 pN)
set kunit [expr 4184/(6.022*pow(10,23)*pow(10,-10))]
# force in N = $k*$kunit

# add vertice of NP (Au3)
foreach index $ACT {
addatom [expr $index+1] 
lappend ACT_list [expr $index+1]
}


proc calcforces {} {
  global ACT_list
  loadcoords c

  foreach index $ACT_list {

    set k 0.05
    set force [vecscale [expr -$k/[veclength $c($index)]] $c($index)]
    addforce $index $force
  }

}

