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
		if( !is.character(obj) )
			return("Slot 'model' must be a character vector")
		invalid.class <- function(cl){ !extends(cl, 'NMF') }
		if( any( inv <- sapply(obj,invalid.class) ) )
			return(paste("Slot 'model' must contain only names of a class that extends class 'NMF' [failure on class(es) "
					, paste( paste("'", obj[inv], "'", sep=''), collapse=', ')  
					,"]"
					, sep=''))
		
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
if ( is.null(getGeneric('run')) ) setGeneric('run', function(method, x, seed, ...) standardGeneric('run'))
setMethod('run', signature(method='NMFStrategy', x='matrix', seed='NMFfit'),
	function(method, x, seed, ...){
		stop("NMFStrategy::run is a pure virtual method that should be overloaded in class '", class(method),"'.")
	}
)

#' Accessor methods to slot \code{name}
if ( is.null(getGeneric('name')) ) setGeneric('name', function(object, ...) standardGeneric('name'))
setMethod('name', signature(object='NMFStrategy'),
	function(object, all=FALSE){
		if( !all ) slot(object, 'name')[1] else slot(object, 'name')
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
	function(object, x, y, ...){
	
		obj.fun <- slot(object, 'objective')
		
		# when both x and y are missing then returns slot objective
		if( missing(x) && missing(y) ) return(obj.fun)
		
		# return the distance computed using the strategy's objective function
		if( !is.function(obj.fun) )
			distance(x, y, method=obj.fun, ...)
		else # directly compute the objective function
			obj.fun(x, y, ...)
		
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
	function(method, key, overwrite=FALSE, save=FALSE, ...){
		
		# build the NMFStrategy
		strategy <- newNMFStrategy(method, key, ...)
		
		# register the method
		nmfRegisterAlgorithm(strategy, key, overwrite, save, ...)
	}
)

#' Factory method to create NMFStrategy objects.
#'
#' Create predefined NMFStrategy objects that implement algorithms from different papers.
if ( !isGeneric('newNMFStrategy') ) setGeneric('newNMFStrategy', function(method, key, ...) standardGeneric('newNMFStrategy') )
setMethod('newNMFStrategy', signature(method='function', key='character'), 
	function(method, key, ...){
			
		# build a NMFStrategyFunction object on the fly to wrap function 'method'
		strategy.params <- list('NMFStrategyFunction', name=key, algorithm=method)
		strategy <- do.call('new', c(strategy.params, list(...)))
		
		# valid the new strategy
		validObject(strategy)
		
		# register the method
		strategy
	}
)

#' Access to registered algorithms
nmfAlgorithm <- function(name=NULL, model, ...){	
	
	algo <- nmfGet(name, registry.name='algorithm', ...)
	if( missing(model) )
		return(algo)
	else{ # lookup for an algorithm suitable for the given NMF model
		if( !is.character(model) || nchar(model) == 0 || !extends(model, 'NMF') )
			stop("argument 'model' must be a the name of class that extends class 'NMF'")
		
		# if the algo was defined then say if it is defined for the given model
		if( inherits(algo, 'NMFStrategy') )
			return( is.element(model, model(algo)) )
		
		# start lookup
		algo.ok <- NULL
		for( name in algo ){
			algo.test <- nmfAlgorithm(name)
			if( is.element(model, model(algo.test)) )
				algo.ok <- c(algo.ok, name)
		}
		
		# return the vector of algorithm names
		algo.ok
	}
}

#' Returns TRUE if the algorithm is registered FALSE otherwise
existsNMFAlgorithm <- function(name, exact=TRUE){	
	
	res <- !is.null( nmfGet(name, registry.name='algorithm', error=FALSE, exact=exact) )
	return(res)
	
}
