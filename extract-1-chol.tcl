package forget


#topology top_all27_lipid.rtf

mol load pdb chol.pdb 

set MEMB1 [ atomselect top "resid 2 and not water" ]


$MEMB1 writepdb chol2.pdb

mol delete all

#mol load psf dopc.psf pdb dopc.pdb

exit 
