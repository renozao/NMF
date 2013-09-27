#######################################################################
#
# Script to compare the results/performances from the NMF package 
# and the original MATLAB implementation of the algorithm from Brunet et al..
#
# This is for when the RcppOctave package is not available, which
# makes things much manual... 
#
#######################################################################

cat("\n# START\n")

# auxiliary function to save data in a plain text format
save.ascii <- function(X, file){
	
	X <- format(X, digits=11)
	write.table(X, file=file, row.names=FALSE, col.names=FALSE, quote=FALSE)	
	invisible()

}

#load NMF
library(NMF)


##################################################
# STEP 1: Generate the test data
##################################################

# load the target: the Golub data
cat("* Prepare data:\n")
cat("\t- Load Golub data\n")
data(esGolub)
target <- exprs(esGolub)
n <- nrow(target)
p <- ncol(target)
rank <- 3

# save target
cat("\t- Save Golub data in ASCII\n")
save.ascii(target, 'target.txt')

# set the seed for reproducibility
set.seed(654321) 
# generate and save W
cat("\t- Generate initial value for W\n")
save.ascii(matrix(runif(n*rank), n, rank), 'W.txt')
# generate and save H
cat("\t- Generate initial value for H\n")
save.ascii(matrix(runif(rank*p), rank, p), 'H.txt')



###########################################################
# STEP 2: Run Brunet et al. algorithm using the NMF package
###########################################################
# load the test data
cat("\n* Reload the generated data\n")
target <- as.matrix(read.table('target.txt'))
W0 <- as.matrix(read.table('W.txt'))
H0 <- as.matrix(read.table('H.txt'))

# wrap the seed into a NMF object
start <- nmfModel(W=W0, H=H0)

# apply Brunet algorithm
cat("* NMF package: run Brunet [optim]\n")
res <- nmf(target, method='brunet', seed=start) # optimized in C++
cat("* NMF package: run Brunet [plain]\n")
res.R <- nmf(target, method='.R#brunet', seed=start) # plain R


######################################################################
# STEP 3: Switch to MATLAB/GNU Octave and run Brunet et al. algorithm 
# 	  using the orginal MATLAB code
######################################################################
# We adapted the original MATLAB code obtained to be able to specify 
# initial values for the matrices W and H.
# Original MATLAB code: http://www.broadinstitute.org/mpr/publications/projects/NMF/nmf.m
# Adapted algorithm: m-files/brunet.m
# MATLAB code to run the test factorization: m-files/brunet-run.m
#
# The results should be saved in a single file named 'ref.brunet.oct'
#, that contains the final values for W and H, in variables named 'W' and 'H'.

###################################################################
# STEP 4: Compare the results
###################################################################
# load the results obtained with the MATLAB code (from Brunet et al. publication)

cat("* Load results from MATLAB/Octave\n")
if( !file.exists('ref.brunet.oct') )
	stop("Could not find file 'ref.brunet.oct'. [Please run the script 'brunet-run.m' to generate it]")

library(foreign)
ref <- read.octave('ref.brunet.oct')

## Note: we used GNU Octave to run the MATLAB code, 
## if the result was obtained using MATLAB, then the following code should load 
## it correctly (we did not test it though).
#library(R.matlab)
#ref <- readMat('ref.brunet.mat')

cat("* Comparison of results:\n\n")
#Sum of differences in W
cat("\t- Sum of absolute differences in W: ", sum( abs(basis(res) -ref$W) ), "\n")
cat("\t [all.equal = ", all.equal(basis(res), ref$W, check.attributes=FALSE), "]\n")

#Sum of differences in H
cat("\t- Sum of absolute differences in H: ", sum( abs(coef(res) - ref$H) ), "\n")
cat("\t [all.equal = ", all.equal(coef(res), ref$H, check.attributes=FALSE), "]\n")

# compare performances
cat("\t- Speed comparaison: \n")
t.ref <- unlist(ref[c('user', 'sys', 'elapsed')])
rbind(`R optim` = runtime(res)[1:3], `MATLAB/Octave`=t.ref, `R plain`=runtime(res.R)[1:3])
cat("\n# DONE\n")

