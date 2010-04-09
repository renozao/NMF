#######################################################################
#
# Script to compare results obtained with the NMF package 
# and the ones obtained from the original MATLAB implementations
#
#######################################################################

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
data(esGolub)
target <- exprs(esGolub)
n <- nrow(target)
p <- ncol(target)
rank <- 3

# save target
save.ascii(target, 'target.txt')

# set the seed for reproducibility
set.seed(654321) 
# generate and save W
save.ascii(matrix(runif(n*rank), n, rank), 'W.txt')
# generate and save H
save.ascii(matrix(runif(rank*p), rank, p), 'H.txt')



###########################################################
# STEP 2: Run Brunet et al. algorithm using the NMF package
###########################################################
# load the test data
target <- as.matrix(read.table('target.txt'))
W0 <- as.matrix(read.table('W.txt'))
H0 <- as.matrix(read.table('H.txt'))

# wrap the seed into a NMF object
start <- newNMF(W=W0, H=H0)

# apply Brunet algorithm
res <- nmf(target, method='brunet', seed=start, .opt='v')



######################################################################
# STEP 3: Switch to MATLAB/GNU Octave and run Brunet et al. algorithm 
# 	  using the orginal MATLAB code
######################################################################
# We adapted the original MATLAB code to be able to specify the initial point
# Brunet algorithm: matlab/nmf-brunet.m
# MATLAB code to run the test factorization: matlab/test-brunet.m
#
# The results should be saved in a single file named 'ref.brunet.oct'
#, that contains the final values for W and H, in variables named 'W' and 'H'.



###################################################################
# STEP 4: Compare the results
###################################################################
# load the results obtained with the MATLAB code (from Brunet et al. publication)
library(foreign)
ref <- read.octave('ref.brunet.oct')

## Note: we used GNU Octave to run the MATLAB code, 
## if the result was obtained using MATLAB, then the following code should load 
## it correctly (we did not test it though).
#library(R.matlab)
#ref <- readMat('ref.brunet.mat')

#Sum of differences in W
sum( abs(basis(res) -ref$W) )

#Sum of differences in H
sum( abs(coef(res) - ref$H) )

