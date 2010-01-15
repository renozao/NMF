#' @include NMFStrategy-class.R
#'
#' Class to define NMF algorithms with a single function.
#'
setClass('NMFStrategyFunction'
	, representation(
		algorithm = 'function' # the function that implements the algorithm				
	)
	, contains = 'NMFStrategy'
)

setMethod('run', signature(object='NMFStrategyFunction', target='matrix', start='NMFfit'),
	function(object, target, start, ...){
		if( !is.function(fun <- algorithm(object)) )  
			stop("NMFFunction '", name(object), "': algorithm not defined.")
		
		# run the function that defines the algorithm and return the result
		fun(target, start, ...)
	}
)

#' Accessor methods to slot \code{algorithm}
if ( is.null(getGeneric('algorithm')) ) setGeneric('algorithm', function(object, ...) standardGeneric('algorithm'))
setMethod('algorithm', signature(object='NMFStrategyFunction'),
	function(object){
		slot(object, 'algorithm')
	}
)
if ( is.null(getGeneric('algorithm<-')) ) setGeneric('algorithm<-', function(object, ..., value) standardGeneric('algorithm<-'))
setReplaceMethod('algorithm', signature(object='NMFStrategyFunction', value='character'),
	function(object, value){
		slot(object, 'algorithm') <- value
		object
	}
)
setReplaceMethod('algorithm', signature(object='NMFStrategyFunction', value='function'),
	function(object, value){
		slot(object, 'algorithm') <- value
		object
	}
)
