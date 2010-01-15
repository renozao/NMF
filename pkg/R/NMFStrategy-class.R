#' Base abstract class that defines the interface for NMF algorithms.
#'
#' @slot name character string giving the name of the strategy
#'
#' @slot objective the objective function associated with the algorithm (Frobenius, Kullback-Leibler, etc...). 
#' It is either a character string as a key registered by \code{nmfRegisterDistance} or a function definition. 
#' In the latter case, the given function must have the following signature (x=matrix, y=matrix) and return a nonnegative real value.
#'
#' @slot model a character string giving either the (sub)class name of the NMF-class instance used and returned by the strategy, or a function name.
#' In the latter case, the given function must have the following signature \code{(v=matrix, r=integer, ...)}, where the \code{v} is the target
#' matrix to approximate and \code{r} is the rank of factorization to achieve.
#'
setClass('NMFStrategy'
	, representation(
				name = 'character' # name of the method (also key)
				, objective = '.functionSlot' # the objective function used to compute the error (defined by name or function)
				, model = 'character' # NMF model to use
				, mixed = 'logical' # can the input data be negative?
	)
	, prototype=prototype(name='', objective='euclidean', model='NMFstd', mixed=FALSE)
	, validity=function(object){
		
		# slot 'name' must be a non-empty character string
		obj <- name(object)
		if( !is.character(obj) || length(obj)!=1 || obj=='' )
			return("Slot 'name' must be a non-empty character string.")
			
		# slot 'objective' must either be a non-empty character string or a function
		obj <- objective(object)
		if( is.character(obj) && obj == '' )
			return("Slot 'objective' must either be a non-empty character string or a function definition.")
			
		# slot 'model' must be the name of a class that extends class 'NMF'
		obj <- model(object)
		if( obj == 'NMF' || !extends(obj, 'NMF') )
			return("Slot 'model' must be the name of a class that STRICTLY extends class 'NMF'.")
		
		# slot 'mixed' must be a single logical		
		obj <- slot(object, 'mixed')
		if( length(obj) != 1 )
			return( paste("Slot 'mixed' must be a single logical [length=", length(obj), "]", sep='') )
	}
	, contains = 'VIRTUAL'
)

setMethod('show', 'NMFStrategy',
		function(object){
			
			cat('<object of class: ', class(object), ">\n")
			cat("name:\t", name(object), "\n")
			svalue <- objective(object)
			svalue <- if( is.function(svalue) ) '<function>' else paste("'", svalue,"'", sep='')
			cat("objective:\t", svalue, "\n")
			cat("NMF model:\t", model(object), "\n")
			return(invisible())
		}
)

#' Main interface to run the algorithm
if ( is.null(getGeneric('run')) ) setGeneric('run', function(object, target, start, ...) standardGeneric('run'))
setMethod('run', signature(object='NMFStrategy', target='matrix', start='NMFfit'),
	function(object, target, start, ...){
		stop("NMFStrategy::run is a pure virtual method that should be overloaded in class '", class(object),"'.")
	}
)

#' Accessor methods to slot \code{name}
if ( is.null(getGeneric('name')) ) setGeneric('name', function(object, ...) standardGeneric('name'))
setMethod('name', signature(object='NMFStrategy'),
	function(object){
		slot(object, 'name')
	}
)
if ( is.null(getGeneric('name<-')) ) setGeneric('name<-', function(object, ..., value) standardGeneric('name<-'))
setReplaceMethod('name', signature(object='NMFStrategy', value='character'),
	function(object, value){
		slot(object, 'name') <- value
		validObject(object)
		object
	}
)

#' Accessor methods to slot \code{objective}
if ( is.null(getGeneric('objective')) ) setGeneric('objective', function(object, ...) standardGeneric('objective'))
setMethod('objective', signature(object='NMFStrategy'),
	function(object, x, y){
	
		# when both x and y are missing then returns slot objective
		if( missing(x) && missing(y) ) return(slot(object, 'objective'))
		
		# return the distance computed using the strategy's objective function
		distance(x, y, method=slot(object, 'objective'))
		
	}
)
if ( is.null(getGeneric('objective<-')) ) setGeneric('objective<-', function(object, ..., value) standardGeneric('objective<-'))
setReplaceMethod('objective', signature(object='NMFStrategy', value='character'),
	function(object, value){
		#TODO: test for the existence of objective method
		slot(object, 'objective') <- value
		validObject(object)
		object
	}
)
setReplaceMethod('objective', signature(object='NMFStrategy', value='function'),
	function(object, value){
		slot(object, 'objective') <- value
		validObject(object)
		object
	}
)

#' Accessor methods to slot \code{model}
if ( !isGeneric('model') ) setGeneric('model', function(object, ...) standardGeneric('model'))
setMethod('model', signature(object='NMFStrategy'),
	function(object){
		slot(object, 'model')
	}
)
if ( is.null(getGeneric('model<-')) ) setGeneric('model<-', function(object, ..., value) standardGeneric('model<-'))
setReplaceMethod('model', signature(object='NMFStrategy', value='character'),
	function(object, value){
		slot(object, 'model') <- value
		validObject(object)
		object
	}
)

#' Accessor methods to slot \code{mixed}
if ( is.null(getGeneric('is.mixed')) ) setGeneric('is.mixed', function(object, ...) standardGeneric('is.mixed'))
setMethod('is.mixed', signature(object='NMFStrategy'),
		function(object){
			return( slot(object, 'mixed') )
		}
)

###########################################################################
# REGISTRY METHODS FOR ALGORITHMS
###########################################################################

#' Register a new algorithm into the NMF registry.
#'
if ( is.null(getGeneric('nmfRegisterAlgorithm')) ) setGeneric('nmfRegisterAlgorithm', function(method, key, ...) standardGeneric('nmfRegisterAlgorithm') )
setMethod('nmfRegisterAlgorithm', signature(method='ANY', key='character'), 
		function(method, key, ...){	
			nmfRegister(method, key, registry.name='algorithm', ...)
		}
)
setMethod('nmfRegisterAlgorithm', signature(method='NMFStrategy', key='missing'), 
		function(method, key, ...){
			
			# get the strategy name
			key <- method@name
			
			# register the method
			nmfRegisterAlgorithm(method, key, ...)		
		}
)
setMethod('nmfRegisterAlgorithm', signature(method='function', key='character'), 
	function(method, key, model, objective, ...){
			
		# build a NMFStrategyFunction object on the fly to wrap function 'method'
		strategy.params <- list('NMFStrategyFunction', name=key, algorithm=method)
		if( !missing(model) ) strategy.params <- c(strategy.params, model=model)
		if( !missing(objective) ) strategy.params <- c(strategy.params, objective=objective)
		strategy <- do.call('new', strategy.params)
		
		# valid the new strategy
		validObject(strategy)
		
		# register the method
		nmfRegisterAlgorithm(strategy, key, ...)
	}
)

#' Factory method to create NMFStrategy objects.
#'
#' Create predefined NMFStrategy objects that implement algorithms from different papers.
nmfAlgorithm <- function(name=NULL, ...){	
	
	nmfGet(name, registry.name='algorithm', ...)
	
}

#' Returns TRUE if the algorithm is registered FALSE otherwise
existsNMFAlgorithm <- function(name, exact=TRUE){	
	
	res <- !is.null( nmfGet(name, registry.name='algorithm', error=FALSE, exact=exact) )
	return(res)
	
}
