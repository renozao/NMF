#' @include NMFstd-class.R
NULL

#' NMF Model - Nonsmooth Nonnegative Matrix Factorization
#' 
#' This class implements the \emph{Nonsmooth Nonnegative Matrix Factorization}
#' (nsNMF) model, required by the Nonsmooth NMF algorithm.
#' 
#' The Nonsmooth NMF algorithm is defined by \cite{Pascual-Montano2006} as a
#' modification of the standard divergence based NMF algorithm (see section
#' Details and references below).  It aims at obtaining sparser factor
#' matrices, by the introduction of a smoothing matrix.
#' 
#' @details
#' The Nonsmooth NMF algorithm is a modification of the standard divergence
#' based NMF algorithm (see \code{\linkS4class{NMF}}).  
#' Given a non-negative \eqn{n \times p}{n x p} matrix \eqn{V} and a 
#' factorization rank \eqn{r}, it fits the following model: 
#' 
#' \deqn{V \equiv W S(\theta) H,}{V &#126; W S(theta) H,} 
#' where: 
#' \itemize{ 
#' 
#' \item \eqn{W} and \eqn{H} are such as in the standard model, i.e. 
#' non-negative matrices of dimension \eqn{n \times r}{n x r}
#' and \eqn{r \times p}{r x p} respectively; 
#' 
#' \item \eqn{S} is a \eqn{r \times r} square matrix whose entries depends on 
#' an extra parameter \eqn{0\leq \theta \leq 1} in the following way: 
#' \deqn{S = (1-\theta)I + \frac{\theta}{r} 11^T ,}
#' where \eqn{I} is the identity matrix and \eqn{1}
#' is a vector of ones.
#'  
#' }
#' 
#' The interpretation of S as a smoothing matrix can be explained as follows:
#' Let \eqn{X} be a positive, nonzero, vector. Consider the transformed vector
#' \eqn{Y = S X}. If \eqn{\theta = 0}, then \eqn{Y = X} and no smoothing on
#' \eqn{X} has occurred.  However, as \eqn{\theta \to 1}{theta tends to 1}, the
#' vector \eqn{Y} tends to the constant vector with all elements almost equal
#' to the average of the elements of \eqn{X}. This is the smoothest possible
#' vector in the sense of non-sparseness because all entries are equal to the
#' same nonzero value, instead of having some values close to zero and others
#' clearly nonzero.
#' 
#' @section Creating objects from the Class:
#' 
#' Object of class \code{NMFns} can be created using the standard way with
#' operator \code{\link{new}}
#' 
#' However, as for all NMF model classes -- that extend class 
#' \code{\linkS4class{NMF}}, objects of class \code{NMFns} should be
#' created using factory method \code{\link{nmfModel}} :
#' 
#' \code{new('NMFns')}
#' 
#' \code{nmfModel(model='NMFns')}
#' 
#' \code{nmfModel(model='NMFns', W=w, theta=0.3}
#' 
#' See \code{\link{nmfModel}} for more details on how to use the factory
#' method.
#' 
#' @section Algorithm:
#' 
#' The Nonsmooth NMF algorithm uses a modified version of the multiplicative
#' update equations in Lee & Seung's method for Kullback-Leibler divergence
#' minimization.  
#' The update equations are modified to take into account the --
#' constant -- smoothing matrix.  
#' The modification reduces to using matrix \eqn{W S} instead of matrix \eqn{W} 
#' in the update of matrix \eqn{H}, and similarly using matrix \eqn{S H} 
#' instead of matrix \eqn{H} in the update of matrix \eqn{W}.
#' 
#' After the matrix \eqn{W} has been updated, each of its columns is scaled so
#' that it sums up to 1.
#' 
#' @export
#' @family NMF-model
#' @examples
#' 
#' # create a completely empty NMFns object
#' new('NMFns')
#' 
#' # create a NMF object based on random (compatible) matrices
#' n <- 50; r <- 3; p <- 20
#' w <- rmatrix(n, r) 
#' h <- rmatrix(r, p)
#' nmfModel(model='NMFns', W=w, H=h)
#' 
#' # apply Nonsmooth NMF algorithm to a random target matrix
#' V <- rmatrix(n, p)
#' \dontrun{nmf(V, r, 'ns')}
#' 
#' # random nonsmooth NMF model  
#' rnmf(3, 10, 5, model='NMFns', theta=0.3)
#' 
setClass('NMFns'
	, representation(
				theta = 'numeric' # smoothing matrix
				)
	, contains = 'NMFstd'
  	, prototype = prototype(
				theta = 0.5
				)
	, validity = function(object){
		if( object@theta < 0 || object@theta > 1 ) 
			return(paste("Invalid value for theta (",object@theta,"): must be between 0 and 1", sep=''))
		TRUE
	}
)

#' Show method for objects of class \code{NMFns}
#' @export
setMethod('show', 'NMFns', 
		function(object)
		{			
			callNextMethod()
			cat("theta:", object@theta, "\n")
		}
)

#' Compute estimate for an NMFns object, according to the Nonsmooth NMF model 
#' (cf. \code{\link{NMFns-class}}).
#' 
#' Extra arguments in \code{...} are passed to method \code{smoothing}, and are 
#' typically used to pass a value for \code{theta}, which is used to compute 
#' the smoothing matrix instead of the one stored in \code{object}.
#' 
#' @param S smoothing matrix to use instead of \code{smoothing(object)}
#' It must be a square matrix compatible with the basis and coefficient matrices 
#' used in the computation.
#' @inline
#' 
setMethod('fitted', signature(object='NMFns'), 
	function(object, W, H, S, ...){
		if( missing(W) ) W <- object@W
		if( missing(H) ) H <- object@H
		if( missing(S) ) S <- smoothing(object, ...)
		W %*% (S %*% H)		
	}
)

#' Smoothing Matrix in Nonsmooth NMF Models
#' 
#' The function \code{smoothing} builds a smoothing matrix for using in Nonsmooth 
#' NMF models.
#'  
#' For a \eqn{r}-rank NMF, the smoothing matrix of parameter \eqn{\theta} is 
#' built as follows:
#' \deqn{S = (1-\theta)I + \frac{\theta}{r} 11^T ,}
#' where \eqn{I} is the identity matrix and \eqn{1} is a vector of ones 
#' (cf. \code{\link{NMFns-class}} for more details).
#' 
#' @param x a object of class \code{NMFns}.
#' @param theta the smoothing parameter (numeric) between 0 and 1.
#' @param ... extra arguments to allow extension (not used)   
#' 
#' @return if \code{x} estimates a \eqn{r}-rank NMF, 
#' then the result is a \eqn{r \times r} square matrix.
#' @export
#' 
#' @examples
#' x <- nmfModel(3, model='NMFns')
#' smoothing(x)
#' smoothing(x, 0.1)
#' 
smoothing <- function(x, theta=x@theta, ...){	
	# check validity of theta
	if( theta < 0 || theta > 1 ) 
		stop("Invalid smoothing parameter theta [",theta,"]: theta must be susch that 0 <= theta <=1")
	diag(1-theta, nbasis(x)) + theta / nbasis(x)		
}

