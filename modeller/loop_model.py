from modeller import *
from modeller.automodel import *    # Load the automodel class

log.verbose()
env = environ()

# directories for input atom files
env.io.atom_files_directory = ['.', '../atom_files']

a = loopmodel(env, alnfile = '2LMO.aln',
		              knowns = '2LMO', 
			      sequence = '2LMO_fill',
			      loop_assess_methods=assess.DOPE)
a.starting_model= 1
a.ending_model  = 1

#a.loop.starting_model = 1
#a.loop.ending_model   = 10
#a.loop.md_level       = refine.slow

a.make()
