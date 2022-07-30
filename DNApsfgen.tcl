#This script will generate a psf file for DNA using VESTA. Currently VESTA's
#psfgen only works for RNA since the CHARMM topology file top_all27_prot_na.inp
#only works for RNA. This script uses psfgen and then patches each base in the
#DNA sequence. 

#modified from both:
#http://www.ks.uiuc.edu/Research/namd/mailing_list/namd-l.2008-2009/3670.html
#http://www.ks.uiuc.edu/~villa/dna/
#accessed on 6/11/2013
# -thomas

resetpsf

mol load pdb DNA.pdb

set newfile DNA_psfgen

set TOPOFILE ./top_all27_prot_na.inp

proc psfalias {} {
        # Define common aliases
        # Here's for nucleics
        pdbalias residue DG GUA
        pdbalias residue DC CYT
        pdbalias residue DA ADE
        pdbalias residue DT THY

        foreach bp { GUA CYT ADE THY URA } {
                pdbalias atom $bp "O5\*" O5'
                pdbalias atom $bp "C5\*" C5'
                pdbalias atom $bp "O4\*" O4'
                pdbalias atom $bp "C4\*" C4'
                pdbalias atom $bp "C3\*" C3'
                pdbalias atom $bp "O3\*" O3'
                pdbalias atom $bp "C2\*" C2'
                pdbalias atom $bp "O2\*" O2'
                pdbalias atom $bp "C1\*" C1'
        }

     }

package require psfgen
topology $TOPOFILE
psfalias

#seperating the DNA into strands

set sel [atomselect top "chain A"]
set listA [lsort -unique [$sel get {resid resname}] ]
foreach selatom $sel {
$selatom set segid ADNA

}
$sel writepdb ADNA.pdb

set sel [atomselect top "chain B"]
set listB [lsort -unique [$sel get {resid resname}] ]
foreach selatom $sel {
$selatom set segid BDNA

}
$sel writepdb BDNA.pdb

# splitting the input pdb

set segidList {ADNA BDNA}
set patchList [list $listA $listB]

segment ADNA {

     first 5TER
     last 3TER

     pdb ADNA.pdb
}

segment BDNA {

     first 5TER
     last 3TER

     pdb BDNA.pdb
}

puts "Begin to patch to make DEOXYribose"

  # patch to make DEOXYribose
  # purine use patch DEO2
  # pyrimidine use patch DEO1

foreach segid $segidList tmpList $patchList {
puts "segid: $segid tmpList:$tmpList"
     foreach record $tmpList {
        foreach {resid resname} $record { break }
        if {$resname == "THY" || $resname == "CYT" } {
            patch DEO1 ${segid}:$resid
            puts "patch DEO1 ${segid}:$resid"
        }

        if {$resname == "ADE" || $resname == "GUA" } {
            puts "patch DEO2 ${segid}:$resid"
            patch DEO2 ${segid}:$resid
        }
     }
     coordpdb ${segid}.pdb

}

puts "End of DEOXYribose "
guesscoord

writepsf $newfile.psf
writepdb $newfile.pdb 
