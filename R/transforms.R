# Transformation methods for matrix-like and NMF objects
# 
# Author: Renaud Gaujoux
# Creation: 19 Jan 2012
###############################################################################

#' @include NMF-class.R
NULL

#' Transforming from Mixed-sign to Nonnegative Data
#' 
#' \code{nneg}  is a generic function to transform a data objects that 
#' contains negative values into a similar object that only contains 
#' values that are nonnegative or greater than a given threshold.
#' 
#' @param object The data object to transform
#' @param ... extra arguments to allow extension or passed down to \code{nneg,matrix}
#' or \code{rposneg,matrix} in subsequent calls.
#' 
#' @return an object of the same class as argument \code{object}.
#' @export
#' @inline
#' @family transforms
#' 
setGeneric('nneg', function(object, ...) standardGeneric('nneg'))
#' Transforms a mixed-sign matrix into a nonnegative matrix, optionally apply a
#' lower threshold. 
#' This is the workhorse method, that is eventually called by all other 
#' methods defined in the \code{\link{NMF}} package.
#' 
#' @param method Name of the transformation method to use, that is partially 
#' matched against the following possible methods:
#' \describe{
#' \item{pmax}{Each entry is constrained to be above threshold \code{threshold}.}
#' 
#' \item{posneg}{The matrix is split into its "positive" and "negative" parts, 
#' with the entries of each part constrained to be above threshold \code{threshold}.
#' The result consists in these two parts stacked in rows (i.e. \code{\link{rbind}}-ed)
#' into a single matrix, which has double the number of rows of the input 
#' matrix \code{object}.}
#' 
#' \item{absolute}{The absolute value of each entry is constrained to be above 
#' threshold \code{threshold}.}
#' 
#' \item{min}{Global shift by adding the minimum entry to each entry, only if 
#' it is negative, and then apply threshold.
#' }
#' 
#' }
#' 
#' @param threshold Nonnegative lower threshold value (single numeric). 
#' See argument \code{shit} for details on how the threshold is used and affects
#' the result.
#' @param shift a logical indicating whether the entries below the threshold 
#' value \code{threshold} should be forced (shifted) to 0 (default) or to 
#' the threshold value itself. 
#' In other words, if \code{shift=TRUE} (default) all entries in 
#' the result matrix are either 0 or strictly greater than \code{threshold}.
#' They are all greater or equal than \code{threshold} otherwise.
#' 
#' @seealso \code{\link{pmax}}
#' @examples
#' 
#' # random mixed sign data (normal distribution)
#' set.seed(1)
#' x <- rmatrix(5,5, rnorm, mean=0, sd=5)
#' x
#' 
#' # pmax (default)
#' nneg(x)
#' # using a threshold
#' nneg(x, threshold=2)
#' # without shifting the entries lower than threshold
#' nneg(x, threshold=2, shift=FALSE)
#' 
#' # posneg: split positive and negative part
#' nneg(x, method='posneg')
#' nneg(x, method='pos', threshold=2)
#' 
#' # absolute
#' nneg(x, method='absolute')
#' nneg(x, method='abs', threshold=2)
#' 
#' # min
#' nneg(x, method='min')
#' nneg(x, method='min', threshold=2)
#' 
setMethod('nneg', 'matrix'
, function(object, method=c('pmax', 'posneg', 'absolute', 'min'), threshold=0, shift=TRUE){
	# match argument
	method <- match.arg(method)
	if( !is.numeric(threshold) || length(threshold) != 1L )
		stop("nneg - Invalid threshold value in argument `threshold` [",threshold,"]: must be a single numeric value.")
	if( threshold < 0 )
		stop("nneg - Invalid threshold value in argument `threshold` [",threshold,"]: must be nonnegative.")
	
	# 1. Transform if there is any negative entry
	m <- min(object)			
	if( m < 0 ){
		object <- 
		switch(method
		, pmax = pmax(object, 0)
		, posneg = rbind(pmax(object, 0), pmax(-object, 0))
		, absolute = pmax(abs(object), 0)
		, min = object - m
		, stop("NMF::nneg - Unexpected error: unimplemented transformation method '", method, "'.")
		)
	}

	if( threshold > 0 ){
		# 2. Apply threshold if any
		object <- pmax(object, threshold)
		
		# 3. Shifting: entries under threshold
		if( shift ) object[object<=threshold] <- 0
	}
	
	# return modified object
	object
}
)

#' Apply \code{nneg} to the basis matrix of an \code{\link{NMF}} 
#' object (i.e. \code{basis(object)}).
#' All extra arguments in \code{...} are passed to the method \code{nneg,matrix}.
#' 
#' @examples
#' 
#' # random 
#' M <- nmfModel(x, rmatrix(ncol(x), 3))
#' nnM <- nneg(M) 
#' basis(nnM)
#' # mixture coefficients are not affected
#' identical( coef(M), coef(nnM) )
#' 
setMethod('nneg', 'NMF', 
	function(object, ...){
		basis(object) <- nneg(basis(object), ...)
		object
	}
)

#' \code{posneg} is a shortcut for \code{nneg(..., method='posneg')}, to split 
#' mixed-sign data into its positive and negative part. 
#' See description for method \code{"posneg"}, in \code{\link{nneg}}.
#' 
#' @export
#' @rdname nneg
#' @examples
#' # shortcut for the "posneg" transformation
#' posneg(x)
#' posneg(x, 2)
#' 
posneg <- function(...) nneg(..., method='posneg')

#' Transforming from Nonnegative to Mixed Sign Data
#' 
#' \code{rposneg} performs the "reverse" transformation of the \code{\link{posneg}} function.
#' 
#' @return an object of the same type of \code{object}
#' @rdname nneg
#' @inline
#' 
setGeneric('rposneg', function(object, ...) standardGeneric('rposneg'))
#' @param unstack Logical indicating whether the positive and negative parts 
#' should be unstacked and combined into a matrix as \code{pos - neg}, which contains 
#' half the number of rows of \code{object} (default), or left 
#' stacked as \code{[pos; -neg]}.
#'   
#' @export
#' @examples
#' 
#' # random mixed sign data (normal distribution)
#' set.seed(1)
#' x <- rmatrix(5,5, rnorm, mean=0, sd=5)
#' x
#'  
#' # posneg-transform: split positive and negative part
#' y <- posneg(x)
#' dim(y)
#' # posneg-reverse
#' z <- rposneg(y)
#' identical(x, z)
#' rposneg(y, unstack=FALSE)
#' 
#' # But posneg-transformation with a non zero threshold is not reversible
#' y1 <- posneg(x, 1)
#' identical(rposneg(y1), x)
#' 
setMethod('rposneg', 'matrix'
, function(object, unstack=TRUE){
	
	# check that the number of rows is pair
	if( nrow(object) %% 2 != 0 )
		stop("rposneg - Invalid input matrix: must have a pair number of rows [",nrow(object),"].")
	n2 <- nrow(object)
	n <- n2/2
	if( unstack ) object <- object[1:n,,drop=FALSE] - object[(n+1):n2,,drop=FALSE]
	else object[(n+1):n2,] <- - object[(n+1):n2,,drop=FALSE]
	
	# return modified object
	object
}
)

#' Apply \code{rposneg} to the basis matrix of an \code{\link{NMF}} object.
#' 
#' @examples
#' 
#' # random mixed signed NMF model 
#' M <- nmfModel(rmatrix(10, 3, rnorm), rmatrix(3, 4))
#' # split positive and negative part
#' nnM <- posneg(M)
#' M2 <- rposneg(nnM)
#' identical(M, M2)
setMethod('rposneg', 'NMF'
, function(object, ...){ 
	basis(object) <- rposneg(basis(object), ...)
	object
}
)

#' Transformation NMF Model Objects
#' 
#' \code{t} transpose an NMF model, by transposing and swapping its basis and 
#' coefficient matrices: \eqn{t([W,H]) = [t(H), t(W)]}.
#' 
#' The function \code{t} is a generic defined in the \pkg{base} package.
#' The method \code{t.NMF} defines the trasnformation for the general NMF interface. 
#' This method may need to be overloaded for NMF models, whose structure requires 
#' specific handling.
#' 
#' @param x NMF model object.
#' 
#' @family transforms
#' @S3method t NMF
#' @examples
#' 
#' x <- rnmf(3, 100, 20)
#' x
#' # transpose
#' y <- t(x)
#' y
#' 
#' # factors are swapped-transposed
#' stopifnot( identical(basis(y), t(coef(x))) )
#' stopifnot( identical(coef(y), t(basis(x))) )
#' 
t.NMF <- function(x){
	# transpose and swap factors
	w <- t(basis(x))
	.basis(x) <- t(coef(x))
	.coef(x) <- w
	# return object
	x
}
