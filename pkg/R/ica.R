#' Seeding method: Absolute Independent Component Analysis
#'
#' @author Renaud Gaujoux
#' @creation 17 Jul 2009


#' Seeding method for Nonnegative Matrix Factorization (NMF) algorithms.
#' 
#' @param object An instance of class \code{NMF} to seed
#' @param x The target matrix
#' @param method The method parameter passed to \code{fastICA}. Can either be 'R' or 'C' and
#' tells which implementation of fastICA to use (R code or C code).
#' @param ... extra parameters passed to \code{fastICA}
#' 
#' @returnType NMF
#' @return an updated version of \code{object}, where the matrix slots \code{W} and \code{H}
#' are set to the positive part of the IC of \code{x}.
#'  
posICA <- function(object, x, ica.method=c('C', 'R'), ...){
			
	# perform ICA using the fastICA package
	if( !require(fastICA) )
		stop("Seeding method requires package fastICA to be installed")
	ica.method <- match.arg(ica.method)
	res <- fastICA(x, nbasis(object), method=ica.method, ...)
	
	# update the 'NMF' object
	object@W <- ifelse( res$S>=0, res$S, .Machine$double.eps ) * res$S ; 
	object@H <- res$A	
	
	# return the updated object
	object
	
} 

.load.seed.ica <- function(){
	
	# Positive ICA
	nmfRegisterSeed(posICA, 'ica', overwrite=TRUE)
	
}
