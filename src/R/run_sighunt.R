#!/usr/local/bin/Rscript --vanilla

# Usage - run_sighunt.R - Use the SigHunt R package to find lateral gene transfer (LGT) in eukaryotic organisms
# Note - This is a template file that will fail if run on its own.  ###input_file### and ###output_path### in the script are replaced with proper paths via a wrapper script (run_sighunt.pl) that copies the template and makes a usable script from it.

# Author - Shaun Adkins (sadkins@som.umaryland.edu)

# Load the SigHunt library
library(sighunt)

cutoff <- 8

# function that executes my analysis on one sequence
# NOTE: code copied from https://github.com/KamilSJaron/sighunt/blob/master/README.md
get_candidates <- function(sequence){
	signature <- get_signature(sequence, window=5000, step=1000)
# can choose from global_density, sliding_density or eye_of_storm
	dias <- global_density(signature)
	candidates <- dias[dias > cutoff]
	#plot(candidates)
# Return vector of candidate positions that exceeds cutoff
	return(candidates)
}

# NOTE about input - Sequence must be > 5000 bases
myseq <- read_fasta('###input_file###')
list_of_candidates <- lapply(myseq, get_candidates)
#lapply(list_of_candidates, write, file='###output_path###/candidates.txt', append=TRUE, ncolumns=1000)
lapply(list_of_candidates, cat, "\n", file='###output_path###/candidates.txt', append=TRUE)
