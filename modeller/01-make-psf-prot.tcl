
 mol new 2BEGf1_model.B99990001.pdb
 set newfile fibril


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

  set sel [atomselect top "chain B"]
    foreach selatom $sel {
        $selatom set segid B
    }
  $sel writepdb B.pdb

  set sel [atomselect top "chain C"]
    foreach selatom $sel {
        $selatom set segid C
    }
  $sel writepdb C.pdb

  set sel [atomselect top "chain D"]
    foreach selatom $sel {
        $selatom set segid D
    }
  $sel writepdb D.pdb

  set sel [atomselect top "chain E"]
    foreach selatom $sel {
        $selatom set segid E
    }
  $sel writepdb E.pdb


##########################################
##########################################

segment A {
  first NTER
  last  CTER
  pdb A.pdb
} 

segment B {
  first NTER
  last  CTER
  pdb B.pdb
} 

segment C {
  first NTER
  last  CTER
  pdb C.pdb
} 

segment D {
  first NTER
  last  CTER
  pdb D.pdb
} 

segment E {
  first NTER
  last  CTER
  pdb E.pdb
} 


#patch DISU U:22 U:61

##########################################
##########################################

coordpdb A.pdb A 	 
coordpdb B.pdb B 	 
coordpdb C.pdb C 	 
coordpdb D.pdb D 	 
coordpdb E.pdb E 	 

guesscoord 	 

#regenerate angle dihedrals
writepdb $newfile.pdb 	 
writepsf $newfile.psf 
