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
#' (e.g. recycling elements).
#' 
#' Null weights can be used to handle missing values in the target matrix.
#' In particular, using \code{weight=NA} cancels out all missing values (see examples).
#'  
#' @param eps small number passed to the standard euclidean-based NMF updates 
#' (see \code{\link{nmf_update.euclidean}}).
#' @param ... extra arguments (not used)
#'  
#' @return updated object \code{object}
#' @rdname lsNMF-nmf
#' @examples 
#' 
#' # Handling missing values in data
#' x <- rmatrix(100, 20)
#' NA_values <- sample(length(x), 5)
#' x[ NA_values ] <- NA
#' 
#' # Classic Lee does not work on this
#' res <- nmf(x, 2, 'lee')
#' anyNA(res) # 3 means NA values in basis and coef matrix
#' 
#' # LS-NMF handles missing values by cancelling them with null weights
#' res <- nmf(x, 2, 'ls-nmf', .opt='v')
#' anyNA(res)
#' 
#' @demo Handling Missing Values in Data
#' 
#' # random data
#' y <- rmatrix(100, 20)
#' # add missing values
#' NA_values <- sample(length(y), 5)
#' y[ NA_values ] <- NA
#' 
#' x <- y
#' # Now a trick: as fixed dummy value (because NA values break other stuffs)
#' x[ NA_values ] <- 123456789
#' 
#' # run ls-nmf using weights that cancel out the missing values
#' w <- matrix(1, nrow(x), ncol(x))
#' w[ NA_values ] <- 0
#' 
#' res <- nmf(x, 3, 'ls-nmf', weight = w)
#' 
#' # The result can be used to input missing values
#' x[ NA_values ] <- fitted(res)[ NA_values ]
#' 
#' # NOTE (to convince yourself that the missing/dummy values are not used)
#' # use a common seed (only fixing RNG does not work here because the range of values in the target matrix affects the initial seed)
#' s <- rnmf(3, y)
#' res <- nmf(x, 3, 'ls-nmf', weight = w, seed = s)
#' 
#' # use another dummy value
#' x2 <- y
#' x2[ NA_values ] <- 987654321
#' res2 <- nmf(x2, 3, 'ls-nmf', weight = w, seed = s)
#' 
#' # results are identical
#' nmf.equal(res, res2)
#' 
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
#' @rdname lsNMF-nmf
wrss <- function(object, X, weight = 1){
	sum( ((X - fitted(object)) * weight)^2 , na.rm = TRUE)/2
}

.ls_nmf <- function(y, x, weight = NA, ...){
    
    verbose <- verbose(x)
    if( !verbose ) message <- function(...) NULL
    
    # check weight size
    message("* Checking LS-NMF weights ... ", appendLF = FALSE)
    if( is.matrix(weight) ){
        if( ncol(weight) != ncol(y) || nrow(weight) != nrow(y) ){
            stop("Invalid weight matrix: dimensions [", str_dim(weight), "] does not match data matrix dimensions [", str_dim(y), "]")
        }else message("OK [", str_dim(weight), "]")
        
    }else if( length(weight) == length(y) ){
        message('NOTE [re-shaped to ', str_dim(y), ']')
        weight <- matrix(head(weight, length(y)), nrow(y))
        
    }else if( !is_NA(weight) ){
        if( length(weight) > length(y) ){
            message('WARNING [truncated and re-shaped to ', str_dim(y), ']')
            warning("LS-NMF: truncated and re-shaped weights to match data matrix dimensions ", str_dim(y))
            weight <- head(weight, length(y))
            
        }else if( length(weight) < length(y) ){
            if( length(y) %% length(weight) == 0 ){
                message('NOTE [recycled and re-shaped to ', str_dim(y), ']')
            }else{
                message('WARNING [incompletely recycled and re-shaped to ', str_dim(y), ']')
                warning("LS-NMF: incompletely recycled and re-shaped weights to match data matrix dimensions ", str_dim(y))
            }
            weight <- rep(weight, length.out = length(y))
        }else message('NOTE [recycled and re-shaped to ', str_dim(y), ']')
        # reshape
        weight <- matrix(head(weight, length(y)), nrow(y))
        
    }else message('SKIP')
    
    # check for NA values
    if( verbose ) message("* Checking for missing values in data matrix ... ", appendLF = FALSE)
    if( anyNA(y) ){
        NA_values <- which(is.na(y))
        
        if( is_NA(weight) ){ # create cancelling weight matrix
            if( verbose ) message("OK [Cancelling ", length(NA_values), " value(s)]")
            weight <- matrix(1, nrow(y), ncol(y))
            
        }else if( nz <- sum(weight[NA_values] != 0, na.rm = TRUE) ){
            if( verbose ) message("WARNING [Cancelling ", nz, " value(s)]")
            warning("LS-NMF: ", nz, " missing value(s) were cancelled out in data matrix by forcing their associated weight to be null.")
        }
        
        weight[NA_values] <- 0
        # using dummy value instead of missing values 
        y[NA_values] <- 10^9
        
    }else if( verbose ) message("OK")
    
    # load plain LS-NMF implementation
    algo <- nmfAlgorithm('.ls-nmf')
    res <- run(algo, y, x, ..., weight = weight)
    
    # store used weights
	res$weight <- weight
    
    # return
    res
}

# Registration of LS-NMF
#' @inheritParams run,NMFStrategyIterative,matrix,NMFfit-method
#' @inheritParams nmf.stop.stationary
#' 
#' @aliases lsNMF-nmf
#' @rdname lsNMF-nmf
nmfAlgorithm._lsNMF <- setNMFMethod('.ls-nmf', objective=wrss
		, Update=nmf_update.lsnmf
		, Stop='stationary')

nmfAlgorithm.lsNMF <- setNMFMethod('ls-nmf', .ls_nmf, objective=wrss)
	
# Unit test for the LS-NMF algorithm
runit.lsnmf <- function(){
	
	set.seed(12345)
	X <- rmatrix(100,20)
	
	res <- nmf(X, 3, 'ls-nmf', weight=1, seed=1)	
	res2 <- nmf(X, 3, '.R#lee', rescale=FALSE, seed=1, .stop=nmf.stop.stationary)
	tol <- 10^-14
	checkTrue( nmf.equal(res, res2, identical=FALSE, tol=tol ), paste("LS-NMF with weight = 1 and .R#Lee (no scale + stationary) give identical results at tolerance=", tol))	
	
}
