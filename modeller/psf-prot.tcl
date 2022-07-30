
 mol new 2LMO_fill.BL00010001.pdb
 set newfile fibril_mono


resetpsf
package require psfgen 

topology top_all27_prot_na.rtf

pdbalias residue HIS HSE 	 
pdbalias atom ILE CD1 CD 	 

#################################
# split the input file
#################################
  set sel [atomselect top "chain A"]
    foreach selatom $sel {
        $selatom set segid A
    }
  $sel writepdb A.pdb

  set sel [atomselect top "chain G"]
    foreach selatom $sel {
        $selatom set segid B
    }
  $sel writepdb G.pdb


##########################################
##########################################

segment A {
  first NTER
  last  CTER
  pdb A.pdb
} 

segment G {
  first NTER
  last  CTER
  pdb G.pdb
} 


#patch DISU U:22 U:61

##########################################
##########################################

coordpdb A.pdb A 	 
coordpdb G.pdb G 	 
guesscoord 	 

#regenerate angle dihedrals
writepdb $newfile.pdb 	 
writepsf $newfile.psf 
