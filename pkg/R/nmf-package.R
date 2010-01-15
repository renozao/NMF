#' Framework to perform Non-negative Matrix Factorization (NMF)
#'
#' \tabular{ll}{
#' Package: \tab nmf\cr
#' Type: \tab Package\cr
#' Version: \tab 0.1\cr
#' Date: \tab 2009-08-01\cr
#' License: \tab GPL (>= 2)\cr
#' LazyLoad: \tab yes\cr
#' }
#'
#' This package provides a framework to perform Non-negative Matrix Factorization (NMF).
#' A implements a set of already plublished algorithms and seeding methods, and provides a framework 
#' to test and develop new algorithms. 
#'
#' \code{\link{nmf}} Run a given NMF algorithm
#'
#' @name nmf-package
#' @aliases nmf
#' @docType package
#' @title Framework to perform Non-negative Matrix Factorization (NMF)
#' @author Renaud Gaujoux \email{renaud@@cbio.uct.ac.za}
#' @references
#' \url{http://www.r-project.org/}
#' @keywords package
#' @seealso \code{\link{nmf}}
#' @examples
#' # create a synthetic matrix
#' V <- syntheticNMF(
#' 
#' # perform a 3-rank NMF using the default algorithm
#' res <- nmf(V, 3)

.onLoad <- function(...){
	
	.init.sequence <- function(){
	
		## 0. INITIALIZE PACKAGE SPECFIC OPTIONS
		.init.nmf.options()
		
		## 1. INITIALIZE THE INTERNAL REGISTRY
		.init.nmf.registry()
		
		## 2. INITIALIZE THE NMF MODELS
		.init.nmf.models()		
		
		## 3. INITIALIZE BIOC LAYER
		if( .init.nmf.bioc() )
			message("NMF :: BioConductor layer loaded\n")
	}
	
	# run intialization sequence suppressing messages or not depending on verbosity options
	if( getOption('verbose') ) .init.sequence()
	else suppressMessages(.init.sequence())
	
	return(invisible())
}

.onAttach <- function(...){
	
	.init.sequence <- function(){
				
		## 1. BUILT-IN NMF METHODS
		# initialize built-in seeding methods
		.load.seed.base() # base: none, random
		.load.seed.nndsvd() # Non-Negative Double SVD
		.load.seed.ica() # Positive part of ICA
			
		# initialize built-in algorithms
		.load.algorithm.NMFStrategyIterative() # iterative schema
		.load.algorithm.snmf() # SNMF methods (R/L)
		.load.algorithm.lnmf() # LNMF
		
		## 2. USER-DEFINED NMF METHODS
		.load.methods.user()
	}
	
	# run intialization sequence suppressing messages or not depending on verbosity options
	if( getOption('verbose') || nmf.getOption('debug') ) .init.sequence()
	else suppressMessages(.init.sequence())	
		
	return(invisible())
}

#' Hook to initialize user-defined methods when the package is loaded 
.load.methods.user <- function(){
	# TODO: load some RData file stored in the user's home R directory	
}
