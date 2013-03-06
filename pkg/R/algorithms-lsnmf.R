# Implementations of LS-NMF
# 
# Reference:
# LS-NMF: a modified non-negative matrix factorization algorithm utilizing uncertainty estimates.
# ﻿Wang, Guoli, Andrew V Kossenkov, and Michael F Ochs. 
# BMC bioinformatics 7 (January 2006): 175. http://www.ncbi.nlm.nih.gov/pubmed/16569230.
#  ﻿
# Author: Renaud Gaujoux
# Creation: 09 Nov 2011
###############################################################################

#' @include registry-algorithms.R
NULL

#' Multiplicative Updates for LS-NMF
#' 
#' Implementation of the updates for the LS-NMF algorithm from \cite{Wang2006}.  
#' 
#' @param i current iteration
#' @param X target matrix
#' @param object current NMF model
#' @param weight value for \eqn{\Sigma}{S}, i.e. the weights that are applied to each 
#' entry in \code{X} by \code{X * weight} (= entry wise product).
#' Weights are usually specified as a matrix of the same dimension as \code{X} 
#' (e.g. uncertainty estimates for each measurement), but may also be passed as a vector, 
#' in which case the standard rules for entry wise product between matrices and vectors apply 
#' (e.g. recylcing elements).
#' @param eps small number passed to the standard euclidean-based NMF updates 
#' (see \code{\link{nmf_update.euclidean}}).
#' @param ... extra arguments (not used)
#'  
#' @return updated object \code{object}
#' @aliases lsnmf-nmf
#' @rdname lsnmf
nmf_update.lsnmf <- function(i, X, object, weight, eps=10^-9, ...)
{
	if( i == 1 ){# pre-compute weighted target matrix
		staticVar('wX', X * weight, init=TRUE)
	}
	
	# retrieve weighted target matrix
	wX <- staticVar('wX')
	
	# retrieve each factor
	w <- .basis(object); h <- .coef(object);	
	
	# compute the estimate WH
	wh <- fitted(object) * weight
	
	# euclidean-reducing NMF iterations	
	# H_au = H_au (W^T V/sigma)_au / (W^T (W H)/sigma)_au
	h <- nmf_update.euclidean.h_R(wX, w, h, wh=wh, eps=eps)	
	
	# update H and recompute the estimate WH
	.coef(object) <- h;
	wh <- fitted(object) * weight
	
	# W_ia = W_ia (V/sigma H^T)_ia / ((W H)/sigma H^T)_ia and columns are rescaled after each iteration	
	w <- nmf_update.euclidean.w_R(wX, w, h, wh=wh, eps=eps)	
	
	#return the modified data
	.basis(object) <- w	
	return(object)
}

#' \code{wrss} implements the objective function used by the LS-NMF algorithm.
#' 
#' @rdname lsnmf 
wrss <- function(object, X, weight){
	sum( ((X - fitted(object)) * weight)^2 )/2
}

# Registration of LS-NMF
setNMFMethod('ls-nmf', objective=wrss
			, Update=nmf_update.lsnmf
			, Stop='stationary')
	
# Unit test for the LS-NMF algorithm
runit.lsnmf <- function(){
	
	set.seed(12345)
	X <- rmatrix(100,20)
	
	res <- nmf(X, 3, 'ls-nmf', weight=1, seed=1)	
	res2 <- nmf(X, 3, '.R#lee', rescale=FALSE, seed=1, .stop=nmf.stop.stationary)
	tol <- 10^-14
	checkTrue( nmf.equal(res, res2, identical=FALSE, tol=tol ), paste("LS-NMF with weight = 1 and .R#Lee (no scale + stationary) give identical results at tolerance=", tol))	
	
}
