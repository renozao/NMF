# Generation of Random NMF Models
# 
# Author: Renaud Gaujoux
# Creation: 03 Jul 2012
###############################################################################

#' @include nmfModel.R
NULL

.rnmf_fixed <- oneoffVariable('none')

#' Generates a random NMF model of the same class and rank as another NMF model.
#' 
#' This is the workhorse method that is eventually called by all other methods.
#' It generates an NMF model of the same class and rank as \code{x}, compatible with the 
#' dimensions specified in \code{target}, that can be a single or 2-length 
#' numeric vector, to specify a square or rectangular target matrix respectively.
#' 
#' The second dimension can also be passed via argument \code{ncol}, so that 
#' calling \code{rnmf(x, 20, 10, ...)} is equivalent to \code{rnmf(x, c(20, 10), ...)}, 
#' but easier to write.
#' 
#' The entries are uniformly drawn between \code{0} and \code{max} 
#' (optionally specified in \code{...}) that defaults to 1.
#' 
#' By default the dimnames of \code{x} are set on the returned NMF model. 
#' This behaviour is disabled with argument \code{keep.names=FALSE}. 
#' See \code{\link{nmfModel}}.
#' 
#' @param ncol single numeric value that specifies the number of columns of the
#' coefficient matrix. Only used when \code{target} is a single numeric value.
#' @param keep.names a logical that indicates if the dimension names of the 
#' original NMF object \code{x} should be conserved (\code{TRUE}) or discarded
#' (\code{FALSE}).
#' @param dist specification of the random distribution to use to draw the entries 
#' of the basis and coefficient matrices.
#' It may be specified as:
#' \itemize{
#' 
#' \item a \code{function} which must be a distribution function such as e.g. 
#' \code{\link{runif}} that is used to draw the entries of both the basis and 
#' coefficient matrices. It is passed in the \code{dist} argument of 
#' \code{\link{rmatrix}}.
#' 
#' \item a \code{list} of arguments that are passed internally to \code{\link{rmatrix}}, 
#' via \code{do.call('rmatrix', dist)}.   
#' 
#' \item a \code{character} string that is partially matched to \sQuote{basis} or 
#' \sQuote{coef}, that specifies which matrix in should be drawn randomly, the 
#' other remaining as in \code{x} -- unchanged.
#' 
#' \item a \code{list} with elements \sQuote{basis} and/or \sQuote{coef}, which 
#' specify the \code{dist} argument separately for the basis and coefficient 
#' matrix respectively.
#' 
#' These elements may be either a distribution function, or a list of arguments that  
#' are passed internally to \code{\link{rmatrix}}, via 
#' \code{do.call('rmatrix', dist$basis)} 
#' or \code{do.call('rmatrix', dist$coef)}.   
#' }
#' 
#' @inline
#' @examples 
#' 
#' ## random NMF of same class and rank as another model
#' 
#' x <- nmfModel(3, 10, 5)
#' x
#' rnmf(x, 20) # square
#' rnmf(x, 20, 13)
#' rnmf(x, c(20, 13))
#' 
#' # using another distribution
#' rnmf(x, 20, dist=rnorm) 
#' 
#' # other than standard model
#' y <- rnmf(3, 50, 10, model='NMFns')
#' y
#' \dontshow{ stopifnot( identical(dim(y), c(50L,10L,3L)) ) }
#' \dontshow{ stopifnot( is(y, 'NMFns') ) }
#' 
setMethod('rnmf', signature(x='NMF', target='numeric'), 
	function(x, target, ncol=NULL, keep.names=TRUE, dist=runif){
		
		# store original dimnames
		if( keep.names ) dn <- dimnames(x)
		
		# valid parameter 'target'
		if( length(target) != 1 && length(target) != 2 )
			stop('NMF::rnmf - invalid target dimensions [length must be 1 or 2. Here length = ', length(target) ,']')
		if( any(is.na(target)) ) 
			stop('NMF::rnmf - invalid target dimensions [NA values in element(s): ', paste(which(is.na(target)), collapse=' and '), ']')		
		# shortcut for symetric case: provide only one dimension
		if( length(target) == 1L ){
			ncol <- if( !is.null(ncol) ){
				if( !is.numeric(ncol) || length(ncol) != 1 || is.na(ncol) )
					stop("NMF::rnmf - invalid argument `ncol`: must be a single numeric value")
				ncol
			}else target
			target <- c(target, ncol)
		}
		
		# retrieve dimension of the target matrix
		n <- target[1]; m <- target[2];
		# retrieve the factorization rank					
		r <- nbasis(x)
		
		## draw basis and coef matrices
		# interpret argument dist
		if( length(dist) == 0L ) dist <- runif
		if( is.character(dist) ){
			dist <- match.arg(dist, c('basis', 'coef'))
			dist <- setNames(list(runif), dist)
		}
		if( is.function(dist) ){
			dist <- list(basis = list(x=n, y=r, dist=dist)
					, coef = list(x=r, y=m, dist=dist))
		}else if( is.list(dist) ){
			if( !all(names(dist) %in% c('basis', 'coef')) ){
				dist <- list(basis=c(list(x=n, y=r), dist)
							, coef=c(list(x=r, y=m), dist))
			}else{
				if( !is.null(dist$basis) )
					dist$basis <- c(list(x=n, y=r), dist$basis)
				if( !is.null(dist$coef) )
					dist$coef <- c(list(x=r, y=m), dist$coef)
			}
		}
		
		fixed <- .rnmf_fixed()
		#Vc# Initialize random matrix: W
		# NB: this will keep the values of fixed basis terms
		if( !is.null(dist$basis) && !('basis' %in% fixed) ){
			basis(x) <- do.call('rmatrix', dist$basis);
		}
		#Vc# Initialize random matrix: H
		# NB: this will keep the values of fixed coef terms
		if( !is.null(dist$coef) && !('coef' %in% fixed) ){
			coef(x) <- do.call('rmatrix', dist$coef);
		}
		
		# if one needs to keep the names (possibly or reducing/increasing) 
		if( keep.names && !is.null(dn) )
			dimnames(x) <- list(dn[[1]][1:n], dn[[2]][1:m], dn[[3]][1:r])
		
		# return the modified object
		x
	}
)

#' Generates a random NMF model compatible and consistent with a target matrix.
#' 
#' The entries are uniformly drawn between \code{0} and \code{max(target)}.
#' It is more or less a shortcut for:
#' \samp{ rnmf(x, dim(target), max=max(target), ...)} 
#' 
#' It returns an NMF model of the same class as \code{x}.
#' 
#' @param use.dimnames a logical that indicates whether the dimnames of the 
#' target matrix should be set on the returned NMF model. 
#' 
#' @inline
#' 
#' @examples
#' # random NMF compatible with a target matrix
#' x <- nmfModel(3, 10, 5)
#' y <- rmatrix(20, 13)
#' rnmf(x, y) # rank of x
#' rnmf(2, y) # rank 2
#' 
setMethod('rnmf', signature(x='ANY', target='matrix'), 
	function(x, target, ..., dist=list(max=max(max(target, na.rm=TRUE), 1)), use.dimnames=TRUE){	
				
		# build a random NMF with the dimensions of the target matrix upper-bounded by the target's maximum entry.
		res <- rnmf(x, dim(target), ..., dist=dist)
		# compute the upper-bound of the random entries and enforce it if possible
		no.na <- abs(target[!is.na(target)])
		if( length(no.na) > 0 ){
			m <- max(no.na)
			basis(res) <- pmin(basis(res), m)
			coef(res) <- pmin(coef(res), m)
		}
		
		# set the dimnames from the target matrix if necessary
		if( use.dimnames )
			dimnames(res) <- dimnames(target)
		
		# return result
		res
	}
)
#' Shortcut for \code{rnmf(x, as.matrix(target))}.
setMethod('rnmf', signature(x='ANY', target='data.frame'),
	function(x, target, ...){
		rnmf(x, as.matrix(target), ...)
	}
)
#' Generates a random NMF model of the same dimension as another NMF model.
#' 
#' It is a shortcut for \code{rnmf(x, nrow(x), ncol(x), ...)}, which returns
#' a random NMF model of the same class and dimensions as \code{x}.
#' 
#' @examples
#' ## random NMF from another model
#' 
#' a <- nmfModel(3, 100, 20)
#' b <- rnmf(a)
#' \dontshow{ stopifnot( !nmf.equal(a,b) ) }
#' 
setMethod('rnmf', signature(x='NMF', target='missing'), 
	function(x, target, ...){
		rnmf(x, c(nrow(x),ncol(x)), ...)
	}
)
#' Generates a random NMF model of a given rank, with known basis and/or 
#' coefficient matrices.  
#'
#' This methods allow to easily generate partially random NMF model, where one 
#' or both factors are known.
#' Although the later case might seems strange, it makes sense for NMF models that 
#' have fit extra data, other than the basis and coefficient matrices, that 
#' are drawn by an \code{rnmf} method defined for their own class, which should 
#' internally call \code{rnmf,NMF,numeric} and let it draw the basis and 
#' coefficient matrices.
#' (e.g. see \code{\linkS4class{NMFOffset}} and \code{\link{rnmf,NMFOffset,numeric-method}}).  
#' 
#' Depending on whether arguments \code{W} and/or \code{H} are missing, 
#' this method interprets \code{x} differently:
#' \itemize{
#' 
#' \item \code{W} provided, \code{H} missing: \code{x} is taken as the number of 
#' columns that must be drawn to build a random coefficient matrix 
#' (i.e. the number of columns in the target matrix).
#' 
#' \item \code{W} is missing, \code{H} is provided: \code{x} is taken as the number of 
#' rows that must be drawn to build a random basis matrix 
#' (i.e. the number of rows in the target matrix).
#' 
#' \item both \code{W} and \code{H} are provided: \code{x} is taken as the target 
#' rank of the model to generate.

#' \item Having both \code{W} and \code{H} missing produces an error, as the 
#' dimension of the model cannot be determined in this case. 
#' }
#'
#' The matrices \code{W} and \code{H} are reduced if necessary and possible 
#' to be consistent with this value of the rank, by the internal call to 
#' \code{\link{nmfModel}}.
#'  
#' All arguments in \code{...} are passed to the function \code{\link{nmfModel}} 
#' which is used to build an initial NMF model, that is in turn passed to 
#' \code{rnmf,NMF,numeric} with \code{dist=list(coef=dist)} or 
#' \code{dist=list(basis=dist)} when suitable.
#' The type of NMF model to generate can therefore be specified in argument 
#' \code{model} (see \code{\link{nmfModel}} for other possible arguments). 
#' 
#' The returned NMF model, has a basis matrix equal to \code{W} (if not missing) 
#' and a coefficient matrix equal to \code{H} (if not missing), or drawn 
#' according to the specification provided in argument \code{dist} 
#' (see method \code{rnmf,NMF,numeric} for details on the supported values for \code{dist}).
#' 
#' @examples
#' # random NMF model with known basis matrix
#' x <- rnmf(5, W=matrix(1:18, 6)) # 6 x 5 model with rank=3
#' basis(x) # fixed 
#' coef(x) # random
#' 
#' # random NMF model with known coefficient matrix
#' x <- rnmf(5, H=matrix(1:18, 3)) # 5 x 6 model with rank=3 
#' basis(x) # random
#' coef(x) # fixed
#' 
#' # random model other than standard NMF
#' x <- rnmf(5, H=matrix(1:18, 3), model='NMFOffset')
#' basis(x) # random
#' coef(x) # fixed
#' offset(x) # random
#' 
setMethod('rnmf', signature(x='numeric', target='missing'),
	function(x, target, ..., W, H, dist=runif){
		
		# get fixed matrices to restore on exit:
		# one must enforce honouring the fixed matrices to prevent the call to 
		# rnmf from a sub-class method to change them.  
		of <- .rnmf_fixed()
		on.exit( .rnmf_fixed(of) )
		
		if( !missing(W) && missing(H) ){ # fixed basis matrix: x = n samples
			# one must not change the values H
			.rnmf_fixed('basis')
			x <- nmfModel(ncol(W), nrow(W), x, W=W, ...)
			dist <- list(coef=dist)
		}else if( missing(W) && !missing(H) ){ # fixed coef matrix: x = n features
			# one must not change the values H
			.rnmf_fixed('coef')
			x <- nmfModel(nrow(H), x, ncol(H), H=H, ...)
			dist <- list(basis=dist)
		}else if( !missing(W) && !missing(H) ){ # fixed basis and coef: x = rank
			# one must not change the values of W and H
			.rnmf_fixed(c('basis', 'coef'))
			x <- nmfModel(x, nrow(W), ncol(H), W=W, H=H, ...)
		}else
			stop("NMF::rnmf - Missing both arguments `W` and/or `H`: at least one of them must be specified.")
		
		rnmf(x, dist=dist) 
	}
)
#' Generates a random NMF model with known basis and coefficient matrices.
#' 
#' This method is a shortcut for calling \code{rnmf,numeric,missing} with a 
#' suitable value for \code{x} (the rank), when both factors are known:
#' code{rnmf(min(ncol(W), nrow(H)), ..., W=W, H=H)}.
#' 
#' Arguments \code{W} and \code{H} are required.
#' Note that calling this method only makes sense for NMF models that contains 
#' data to fit other than the basis and coefficient matrices, 
#' e.g. \code{\linkS4class{NMFOffset}}. 
#' 
#' @examples
#' 
#' # random model other than standard NMF
#' x <- rnmf(W=matrix(1:18, 6), H=matrix(21:38, 3), model='NMFOffset')
#' basis(x) # fixed
#' coef(x) # fixed
#' offset(x) # random
#' 
setMethod('rnmf', signature(x='missing', target='missing'),
	function(x, target, ..., W, H){
		rnmf(min(ncol(W), nrow(H)), ..., W=W, H=H)
	}
)

#' Generates a random standard NMF model of given dimensions.
#' 
#' This is a shortcut for \code{rnmf(nmfModel(x, target, ncol, ...)), dist=dist)}.
#' It generates a standard NMF model compatible with the dimensions passed in 
#' \code{target}, that can be a single or 2-length numeric vector, to specify 
#' a square or rectangular target matrix respectively. 
#' See \code{\link{nmfModel}}.
#' 
#' @inheritParams nmfModel,numeric,numeric-method
#' 
#' @examples
#' 
#' ## random standard NMF of given dimensions
#' 
#' # generate a random NMF model with rank 3 that fits a 100x20 matrix  
#' rnmf(3, 100, 20)
#' \dontshow{ stopifnot( identical(dim(rnmf(3, 100, 20)), c(100L,20L,3L)) ) }
#' # generate a random NMF model with rank 3 that fits a 100x100 matrix
#' rnmf(3, 100)
#' \dontshow{ stopifnot( identical(dim(rnmf(3, 100)), c(100L,100L,3L)) ) }
#' 
setMethod('rnmf', signature(x='numeric', target='numeric'), 
		function(x, target, ncol=NULL, ..., dist=runif){		
			rnmf(nmfModel(x, target, ncol, ...), dist=dist)
		}
)
#' Generate a random formula-based NMF model, using the method 
#' \code{\link{nmfModel,formula,ANY-method}}.
setMethod('rnmf', signature(x='formula', target='ANY'), 
	function(x, target, ..., dist=runif){
		# missing target is NULL
		if( missing(target) ) target <- NULL
		rnmf(nmfModel(x, target, ...), dist=dist)
	}
)

