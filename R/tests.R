# Utility functions to test NMF algorithms
# 
# Author: Renaud Gaujoux
# Created: 29 Nov 2012
###############################################################################

#' @include rmatrix.R
NULL

#' Checking NMF Algorithm
#' 
#' \code{nmfCheck} enables to quickly check that a given NMF algorithm runs 
#' properly, by applying it to some small random data.
#' 
#' @param method name of the NMF algorithm to be tested.
#' @param rank rank of the factorization
#' @param x target data. If \code{NULL}, a random 20 x 10 matrix is generated
#' @param seed specifies a seed or seeding method for the computation.  
#' @param ... other arguments passed to the call to \code{\link{nmf}}.
#' 
#' @return the result of the NMF fit invisibly.
#' 
#' @export
#' @examples 
#' 
#' # test default algorithm
#' nmfCheck()
#' 
#' # test 'lee' algorithm
#' nmfCheck('lee')
#' 
nmfCheck <- function(method=NULL, rank=max(ncol(x)/5, 3), x=NULL, seed=1234, ...){
	
	# seed computation
	if( isNumber(seed) ){
		os <- RNGseed()
		on.exit( RNGseed(os), add=TRUE)
		set.seed(seed)
		seed <- NULL
	}
	if( is.null(x) ){
		x <- rmatrix(20, 10)
	}
	res <- nmf(x, rank, method, seed=seed, ...)
}
