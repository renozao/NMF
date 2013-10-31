# Generation of random matrices
# 
# Defines the generic function `rmatrix` and basic methods for it.
#
# Author: Renaud Gaujoux
###############################################################################


#' Generating Random Matrices
#' 
#' The S4 generic \code{rmatrix} generates a random matrix from a given object.  
#' Methods are provided to generate matrices with entries drawn from any 
#' given random distribution function, e.g. \code{\link{runif}} or 
#' \code{\link{rnorm}}.
#' 
#' @param x object from which to generate a random matrix 
#' 
#' @export
setGeneric('rmatrix', function(x, ...) standardGeneric('rmatrix'))
#' Generates a random matrix of given dimensions, whose entries 
#' are drawn using the distribution function \code{dist}.
#' 
#' This is the workhorse method that is eventually called by all other methods.
#' It returns a matrix with:
#' \itemize{
#' \item \code{x} rows and \code{y} columns if \code{y} is not missing and 
#' not \code{NULL};
#' \item dimension \code{x[1]} x \code{x[2]} if \code{x} has at least two elements;
#' \item dimension \code{x} (i.e. a square matrix) otherwise.
#' }
#' 
#' The default is to draw its entries from the standard uniform distribution using
#' the base function \code{\link{runif}}, but any other function that generates 
#' random numeric vectors of a given length may be specified in argument \code{dist}.
#' All arguments in \code{...} are passed to the function specified in \code{dist}.
#' 
#' The only requirement is that the function in \code{dist} is of the following form:
#' 
#' \samp{
#' function(n, ...){
#' # return vector of length n
#' ...
#' }}
#' 
#' This is the case of all base random draw function such as \code{\link{rnorm}}, 
#' \code{\link{rgamma}}, etc\ldots
#'  
#' 
#' @param y optional specification of number of columns
#' @param dist a random distribution function or a numeric seed (see details of method 
#' \code{rmatrix,numeric})
#' @param byrow a logical passed in the internal call to the function 
#' \code{\link{matrix}}
#' @param dimnames \code{NULL} or a \code{list} passed in the internal call to 
#' the function \code{\link{matrix}}
#' @param ... extra arguments passed to the distribution function \code{dist}.
#' 
#' @inline
#' 
#' @examples
#' ## Generate a random matrix of a given size
#' rmatrix(5, 3)
#' \dontshow{ stopifnot( identical(dim(rmatrix(5, 3)), c(5L,3L)) ) }
#' 
#' ## Generate a random matrix of the same dimension of a template matrix
#' a <- matrix(1, 3, 4)
#' rmatrix(a)
#' \dontshow{ stopifnot( identical(dim(rmatrix(a)), c(3L,4L)) ) }
#' 
#' ## Specificy the distribution to use
#' 
#' # the default is uniform
#' a <- rmatrix(1000, 50)
#' \dontrun{ hist(a) }
#' 
#' # use normal ditribution
#' a <- rmatrix(1000, 50, rnorm)
#' \dontrun{ hist(a) }
#' 
#' # extra arguments can be passed to the random variate generation function 
#' a <- rmatrix(1000, 50, rnorm, mean=2, sd=0.5)
#' \dontrun{ hist(a) }
#' 
setMethod('rmatrix', 'numeric', 
		function(x, y=NULL, dist=runif, byrow = FALSE, dimnames = NULL, ...){
			
			x <- as.integer(x)
			# early exit if x has length 0
			if( length(x) == 0L )
				stop("NMF::rmatrix - invalid empty vector in argument `x`.")
			
			# check/ensure that 'dist' is a function.
			if( is.null(dist) ) dist <- runif
			if( isNumber(dist) ){
				os <- RNGseed()
				on.exit( RNGseed(os), add=TRUE)
				set.seed(dist)
				dist <- runif
			}
			if( !is.function(dist) )
				stop("NMF::rmatrix - invalid value for argument 'dist': must be a function [class(dist)='", class(dist), "'].")
			
			# if 'y' is not specified:
			if( is.null(y) ){
				
				if( length(x) == 1L ) y <- x # create a square matrix 
				else{ # assume x contains all dimensions (e.g. returned by dim())
					y <- x[2L]
					x <- x[1L]
				}
				
			}else{
				y <- as.integer(y)
				y <- y[1L] # only use first element
			}
			
			# build the random matrix using the distribution function
			matrix(dist(x*y, ...), x, y, byrow=byrow, dimnames=dimnames)	
		}
)

#' Default method which calls \code{rmatrix,vector} on the dimensions of \code{x}
#' that is assumed to be returned by a suitable \code{dim} method:
#' it is equivalent to \code{rmatrix(dim(x), y=NULL, ...)}.
#' 
#' @examples
#' 
#' # random matrix of the same dimension as another matrix
#' x <- matrix(3,4)
#' dim(rmatrix(x))
#' 
setMethod('rmatrix', 'ANY', 
		function(x, ...){
			rmatrix(x=dim(x), y=NULL, ...)
		}
)

