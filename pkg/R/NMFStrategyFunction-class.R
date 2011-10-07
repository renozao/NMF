###% Class to define NMF algorithms with a single function.
###%
###% @include NMFStrategy-class.R
setClass('NMFStrategyFunction'
	, representation(
		algorithm = 'function' # the function that implements the algorithm				
	)
	, contains = 'NMFStrategy'
)

setMethod('run', signature(method='NMFStrategyFunction', x='matrix', seed='NMFfit'),
	function(method, x, seed, ...){
		if( !is.function(fun <- algorithm(method)) )  
			stop("NMFStrategyFunction '", name(method), "': algorithm not defined.")
		
		# run the function that defines the algorithm and return the result
		fun(x, seed, ...)
	}
)

###% Accessor methods to slot \code{algorithm}
setMethod('algorithm', signature(object='NMFStrategyFunction'),
	function(object){
		slot(object, 'algorithm')
	}
)
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
