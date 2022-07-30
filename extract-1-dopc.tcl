package forget
package require psfgen 1.6

#topology top_all27_lipid.rtf

mol load pdb dopc.pdb 

set lipid1 [ atomselect top "resid 2 and not water" ]

$lipid1 writepdb dopc1.pdb

mol delete all

#mol pdb dopc.pdb

exit 
