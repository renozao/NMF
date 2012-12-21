#' @include NMFstd-class.R
NULL

#' NMF Model - Nonnegative Matrix Factorization with Offset
#' 
#' This class implements the \emph{Nonnegative Matrix Factorization with
#' Offset} model, required by the NMF with Offset algorithm.
#' 
#' The NMF with Offset algorithm is defined by \cite{Badea2008} as a modification
#' of the euclidean based NMF algorithm from \code{Lee2001} (see section Details and
#' references below).  
#' It aims at obtaining 'cleaner' factor matrices, by the introduction of an 
#' offset matrix, explicitly modelling a feature specific baseline 
#' -- constant across samples.
#' 
#' @section Creating objects from the Class:
#' 
#' Object of class \code{NMFOffset} can be created using the standard way with
#' operator \code{\link{new}}
#' 
#' However, as for all NMF model classes -- that extend class 
#' \code{\linkS4class{NMF}}, objects of class \code{NMFOffset} should be
#' created using factory method \code{\link{nmfModel}} :
#' 
#' \code{new('NMFOffset')}
#' 
#' \code{nmfModel(model='NMFOffset')}
#' 
#' \code{nmfModel(model='NMFOffset', W=w, offset=rep(1, nrow(w)))}
#' 
#' See \code{\link{nmfModel}} for more details on how to use the factory
#' method.
#' 
#' @export
#' @family NMF-model
#' @examples
#'  
#' # create a completely empty NMF object
#' new('NMFOffset')
#' 
#' # create a NMF object based on random (compatible) matrices
#' n <- 50; r <- 3; p <- 20
#' w <- rmatrix(n, r) 
#' h <- rmatrix(r, p)
#' nmfModel(model='NMFOffset', W=w, H=h, offset=rep(0.5, nrow(w)))
#' 
#' # apply Nonsmooth NMF algorithm to a random target matrix
#' V <- rmatrix(n, p)
#' \dontrun{nmf(V, r, 'offset')}
#'
#' # random NMF model with offset  
#' rnmf(3, 10, 5, model='NMFOffset')
#' 
setClass('NMFOffset'
	, representation(
				offset = 'numeric' # offset vector
				)
	, contains = 'NMFstd'
  	, prototype=prototype(
  				offset = numeric()
				)
	
)

#' Show method for objects of class \code{NMFOffset}
#' @export
setMethod('show', 'NMFOffset', 
	function(object)
	{		
		callNextMethod()
		cat("offset: ")
		if( length(object@offset) > 0 ){
			cat('[', head(object@offset, 5)
				, if( length(object@offset) > 5 ) "..." else NULL
				, ']')
		}
		else cat('none')
		cat("\n")
	}
)

#' @section Initialize method:
#' The initialize method for \code{NMFOffset} objects tries to correct the initial 
#' value passed for slot \code{offset}, so that it is consistent with the dimensions 
#' of the \code{NMF} model:
#' it will pad the offset vector with NA values to get the length equal to the 
#' number of rows in the basis matrix.
#' 
#' @param offset optional numeric vector used to initialise slot \sQuote{offset}.
#' 
#' @rdname NMFOffset-class
setMethod("initialize", 'NMFOffset', 
		function(.Object, ..., offset){			
			.Object <- callNextMethod()
			# correct the offset slot if possible
			if( missing(offset) ) offset <- numeric()
			if( !is.numeric(offset) ) stop("Unvalid value for parameter 'offset': a numeric vector is expected")			
			 
			# force length to be consistent with the factorization's dimension
			n <- nrow(.Object)
			if( n > 0 ) .Object@offset <- c( offset, rep(NA, max(0, n - length(offset))) )[1:n]
			
			# return the initialized valid object
			.Object
		}
)

#' @export
setGeneric('offset', package='stats')
#' Offsets in NMF Models with Offset
#' 
#' The function \code{offset} returns the offset vector from an NMF model 
#' that has an offset, e.g. an \code{NMFOffset} model.
#' @param object an instance of class \code{NMFOffset}.
#' 
setMethod('offset', signature(object='NMFOffset'), 
	function(object){
		object@offset
	}
)

#' Computes the target matrix estimate for an NMFOffset object.
#' 
#' The estimate is computed as:
#' \deqn{ W H + offset }
#' 
#' @param offset offset vector
#' @inline
setMethod('fitted', signature(object='NMFOffset'), 
	function(object, W, H, offset=object@offset){
		if( missing(W) ) W <- object@W
		if( missing(H) ) H <- object@H
		object@W %*% object@H + offset
	}
)

#' Generates a random NMF model with offset, from class \code{NMFOffset}.
#' 
#' The offset values are drawn from a uniform distribution between 0 and 
#' the maximum entry of the basis and coefficient matrices, which are drawn 
#' by the next suitable \code{\link{rnmf}} method, which is the workhorse 
#' method \code{rnmf,NMF,numeric}.
#'   
#' @examples
#' 
#' # random NMF model with offset
#' x <- rnmf(2, 3, model='NMFOffset')
#' x
#' offset(x)
#' # from a matrix
#' x <- rnmf(2, rmatrix(5,3, max=10), model='NMFOffset')
#' offset(x)
#' 
setMethod('rnmf', signature(x='NMFOffset', target='numeric'), 
function(x, target, ...){	
	
	# call the parent's 'rnmf' method to build a standard random NMF factorization
	res <- callNextMethod()
		
	#Vc# Initialize a random offset of length the number of genes
	res@offset <- runif(nrow(res), min=0, max=max(basis(res), coef(res)));
	
	# return the initialized NMFOffset object
	res
})
