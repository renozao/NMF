# Class that implements the standard NMF model
# 
# Author: Renaud Gaujoux \email{renaud@@cbio.uct.ac.za}
###############################################################################

#' @include NMF-class.R
NULL

#' NMF Model - Standard model 
#' 
#' This class implements the standard model of Nonnegative Matrix
#' Factorization.
#' It provides a general structure and generic functions to manage
#' factorizations that follow the standard NMF model, as defined by 
#' \cite{Lee2001}.
#' 
#' Let \eqn{V} be a \eqn{n \times m} non-negative matrix and \eqn{r} a positive
#' integer.  In its standard form (see references below), a NMF of \eqn{V} is
#' commonly defined as a pair of matrices \eqn{(W, H)} such that:
#' 
#' \deqn{V \equiv W H,}
#' 
#' where: 
#' \itemize{ 
#' \item \eqn{W} and \eqn{H} are \eqn{n \times r} and \eqn{r
#' \times m} matrices respectively with non-negative entries; 
#' \item \eqn{\equiv} is to be understood with respect to some loss function.  
#' Common choices of loss functions are based on Frobenius norm or Kullback-Leibler
#' divergence.  
#' }
#' 
#' Integer \eqn{r} is called the \emph{factorization rank}.  
#' Depending on the context of application of NMF, the columns of \eqn{W} 
#' and \eqn{H} are given different names: 
#' \describe{ 
#' \item{columns of \code{W}}{basis vector, metagenes, factors, source, image basis}
#' \item{columns of \code{H}}{mixture coefficients, metagene sample expression profiles, weights}
#' \item{rows of \code{H}}{basis profiles, metagene expression profiles}
#' }
#' 
#' NMF approaches have been successfully applied to several fields. 
#' The package NMF was implemented trying to use names as generic as possible 
#' for objects and methods.  
#' 
#' The following terminology is used: 
#' \describe{ 
#' \item{samples}{the columns of the target matrix \eqn{V}} 
#' \item{features}{the rows of the target matrix \eqn{V}}
#' \item{basis matrix}{the first matrix factor \eqn{W}}
#' \item{basis vectors}{the columns of first matrix factor \eqn{W}}
#' \item{mixture matrix}{the second matrix factor \eqn{H}} \item{mixtures
#' coefficients}{the columns of second matrix factor \eqn{H}} 
#' }
#' 
#' However, because the package NMF was primarily implemented to work with gene
#' expression microarray data, it also provides a layer to easily and
#' intuitively work with objects from the Bioconductor base framework.  
#' See \link{bioc-NMF} for more details.
#' 
#' @slot W A \code{matrix} that contains the basis matrix, i.e. the \emph{first} 
#' matrix factor of the factorisation
#' @slot H A \code{matrix} that contains the coefficient matrix, i.e. the 
#' \emph{second} matrix factor of the factorisation
#' @slot bterms a \code{data.frame} that contains the primary data that 
#' define fixed basis terms. See \code{\link{bterms}}.
#' @slot ibterms integer vector that contains the indexes of the basis components
#' that are fixed, i.e. for which only the coefficient are estimated.
#' 
#' IMPORTANT: This slot is set on construction of an NMF model via 
#' \code{\link[=nmfModel,formula,ANY-method]{nmfModel}} and is not recommended to 
#' not be subsequently changed by the end-user.
#' @slot cterms  a \code{data.frame} that contains the primary data that 
#' define fixed coefficient terms. See \code{\link{cterms}}.
#' @slot icterms integer vector that contains the indexes of the basis components
#' that have fixed coefficients, i.e. for which only the basis vectors are estimated.
#' 
#' IMPORTANT: This slot is set on construction of an NMF model via 
#' \code{\link[=nmfModel,formula,ANY-method]{nmfModel}} and is not recommended to 
#' not be subsequently changed by the end-user. 
#' 
#' @export
#' @family NMF-model 
#' @examples 
#' # create a completely empty NMFstd object
#' new('NMFstd')
#' 
#' # create a NMF object based on one random matrix: the missing matrix is deduced
#' # Note this only works when using factory method NMF 
#' n <- 50; r <- 3; 
#' w <- rmatrix(n, r) 
#' nmfModel(W=w)
#' 
#' # create a NMF object based on random (compatible) matrices
#' p <- 20
#' h <- rmatrix(r, p)
#' nmfModel(W=w, H=h)
#' 
#' # create a NMF object based on incompatible matrices: generate an error
#' h <- rmatrix(r+1, p)
#' try( new('NMFstd', W=w, H=h) )
#' try( nmfModel(w, h) )
#' 
#' # Giving target dimensions to the factory method allow for coping with dimension
#' # incompatibilty (a warning is thrown in such case) 
#' nmfModel(r, W=w, H=h)
#'  
setClass('NMFstd' 
		, representation(
			W = 'matrix' # basis matrix
			, H = 'matrix' # mixture coefficients matrix
			, bterms = 'data.frame' # fixed basis terms: nrow(bterms) = nrow(x)
			, ibterms = 'integer' # index of the fixed basis terms
			, cterms = 'data.frame' # fixed coef terms: ncol(cterms) = ncol(x)
			, icterms = 'integer' # index of the fixed coefficient terms
		)
		
		, prototype = prototype(
				W = matrix(as.numeric(NA), 0, 0),
				H = matrix(as.numeric(NA), 0, 0)
		)
		
		, validity = function(object){
			
			# dimension compatibility: W and H must be compatible for matrix multiplication
			if( ncol(object@W) != nrow(object@H) ){
				return(paste('Dimensions of W and H are not compatible [ncol(W)=', ncol(object@W) , '!= nrow(H)=', nrow(object@H), ']'))
			}
			# give a warning if the dimensions look strange: rank greater than the number of samples
			if( !is.empty.nmf(object) && ncol(object@H) && ncol(object@W) > ncol(object@H) ){
				warning(paste('Dimensions of W and H look strange [ncol(W)=', ncol(object@W) , '> ncol(H)=', ncol(object@H), ']'))
			}
			
			# everything went fine: return TRUE
			return(TRUE)
		}
		, contains = 'NMF'
)


#' Get the basis matrix in standard NMF models 
#' 
#' This function returns slot \code{W} of \code{object}.
#' 
#' @examples
#' # random standard NMF model
#' x <- rnmf(3, 10, 5)
#' basis(x)
#' coef(x)
#' 
#' # set matrix factors
#' basis(x) <- matrix(1, nrow(x), nbasis(x))
#' coef(x) <- matrix(1, nbasis(x), ncol(x))
#' # set random factors
#' basis(x) <- rmatrix(basis(x))
#' coef(x) <- rmatrix(coef(x))
#' 
#' # incompatible matrices generate an error:
#' try( coef(x) <- matrix(1, nbasis(x)-1, nrow(x)) )
#' # but the low-level method allow it
#' .coef(x) <- matrix(1, nbasis(x)-1, nrow(x))
#' try( validObject(x) )
#' 
setMethod('.basis', 'NMFstd',
	function(object){ 
		object@W
	}
)
#' Set the basis matrix in standard NMF models 
#' 
#' This function sets slot \code{W} of \code{object}.
setReplaceMethod('.basis', signature(object='NMFstd', value='matrix'), 
	function(object, value){ 
		object@W <- value		
		object
	} 
)

#' Get the mixture coefficient matrix in standard NMF models 
#' 
#' This function returns slot \code{H} of \code{object}.
setMethod('.coef', 'NMFstd',
	function(object){
		object@H
	}
)
#' Set the mixture coefficient matrix in standard NMF models 
#' 
#' This function sets slot \code{H} of \code{object}.
setReplaceMethod('.coef', signature(object='NMFstd', value='matrix'), 
	function(object, value){ 
		object@H <- value			
		object
	}
)

#' Compute the target matrix estimate in \emph{standard NMF models}.
#' 
#' The estimate matrix is computed as the product of the two matrix slots 
#' \code{W} and \code{H}:
#' \deqn{\hat{V} = W H}{V ~ W H} 
#' 
#' @param W a matrix to use in the computation as the basis matrix in place of 
#' \code{basis(object)}. 
#' It must be compatible with the coefficient matrix used 
#' in the computation (i.e. number of columns in \code{W} = number of rows in \code{H}).
#' @param H a matrix to use in the computation as the coefficient matrix in place of 
#' \code{coef(object)}. 
#' It must be compatible with the basis matrix used 
#' in the computation (i.e. number of rows in \code{H} = number of columns in \code{W}).
#'  
#' @export
#' @inline
#' 
#' @examples
#' # random standard NMF model
#' x <- rnmf(3, 10, 5)
#' all.equal(fitted(x), basis(x) %*% coef(x))
#' 
#' 
setMethod('fitted', signature(object='NMFstd'), 
	function(object, W, H, ...){
		if( missing(W) ) W <- object@W
		if( missing(H) ) H <- object@H
		return(W %*% H)
	}
)

