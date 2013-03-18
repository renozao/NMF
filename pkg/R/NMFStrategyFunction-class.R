#' @include NMFStrategy-class.R
NULL

#' Interface for Single Function NMF Strategies 
#'
#' This class implements the virtual interface \code{\link{NMFStrategy}} for 
#' NMF algorithms that are implemented by a single workhorse R function. 
#' 
#' @slot algorithm a function that implements an NMF algorithm.
#' It must have signature \code{(y='matrix', x='NMFfit')}, where \code{y} is the 
#' target matrix to approximate and \code{x} is the NMF model assumed to be  
#' seeded with an appropriate initial value -- as it is done internally by 
#' function \code{\link{nmf}}.
#' 
#' Note that argument names currently do not matter, but it is recommended to 
#' name them as specified above.
#' 
setClass('NMFStrategyFunction'
	, representation(
		algorithm = 'function' # the function that implements the algorithm				
	)
	, contains = 'NMFStrategy'
)
#' Runs the NMF algorithms implemented by the single R function -- and stored in slot \code{'algorithm'} 
#' of \code{object}, on the data object \code{y}, using \code{x} as starting point.
#' It is equivalent to calling \code{object@@algorithm(y, x, ...)}.
#' 
#' This method is usually not called directly, but only via the function \code{\link{nmf}}, which 
#' takes care of many other details such as seeding the computation, handling RNG settings, or setting up 
#' parallelisation.
#' 
#' @rdname NMFStrategy
setMethod('run', signature(object='NMFStrategyFunction', y='matrix', x='NMFfit'),
	function(object, y, x, ...){
		if( !is.function(fun <- algorithm(object)) )  
			stop("NMFStrategyFunction '", name(object), "': algorithm is not defined.")
		
		# run the function that defines the algorithm and return the result
		fun(y, x, ...)
	}
)

#' @S3method nmfFormals NMFStrategyFunction
nmfFormals.NMFStrategyFunction <- function(x, ...){
	args <- formals(x@algorithm)
	args[-(1:2)]
}
	

#' Returns the single R function that implements the NMF algorithm -- as stored in 
#' slot \code{algorithm}.
setMethod('algorithm', signature(object='NMFStrategyFunction'),
	function(object){
		slot(object, 'algorithm')
	}
)
#setReplaceMethod('algorithm', signature(object='NMFStrategyFunction', value='character'),
#	function(object, value){
#		slot(object, 'algorithm') <- value
#		object
#	}
#)
#' Sets the function that implements the NMF algorithm, stored in slot \code{algorithm}. 
setReplaceMethod('algorithm', signature(object='NMFStrategyFunction', value='function'),
	function(object, value){
		slot(object, 'algorithm') <- value
		object
	}
)
