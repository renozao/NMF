#' Seeding methods for NMF package.
#' 
#' @author Renaud Gaujoux
#' @created 22 Jul 2009 
#'
#' @include registry.R 
#' @include nndsvd.R  
NA

#' Base class that defines the interface for NMF seeding methods.
#'
#' @slot name character string giving the name of the strategy
#'
#'
setClass('NMFSeed'
		, representation(
				name = 'character' # name of the method (also key)
				, method = 'function' # the method actual definition
		)
		, prototype=prototype(name='')
		, validity=function(object){
			
			# slot 'name' must be a non-empty character string
			obj <- name(object)
			if( !is.character(obj) || length(obj)!=1 || obj=='' )
				return("Slot 'name' must be a non-empty character string.")
			
			# slot 'method' must either be a non-empty character string or a function
#			obj <- method(object)
#			if( !(is.character(obj) && obj != '') && !is.function(obj) )
#				return("Slot 'method' must either be a non-empty character string or a defined function.")

			return(TRUE)
		}		
)

setMethod('show', 'NMFSeed',
		function(object){			
			cat('<object of class: ', class(object), ">\n")
			cat("name:\t", name(object), "\n")
			svalue <- method(object)
			svalue <- if( is.function(svalue) ) '<function>' else paste("'", svalue,"'", sep='')
			cat("method:\t", svalue, "\n")
			return(invisible())
		}
)

#' Accessor methods to slot \code{name}
if ( !isGeneric('name') ) setGeneric('name', function(object, ...) standardGeneric('name'))
setMethod('name', signature(object='NMFSeed'),
		function(object){
			slot(object, 'name')
		}
)
if ( !isGeneric('name<-') ) setGeneric('name<-', function(object, ..., value) standardGeneric('name<-'))
setReplaceMethod('name', signature(object='NMFSeed', value='character'),
		function(object, value){
			slot(object, 'name') <- value
			validObject(object)
			object
		}
)

#' Accessor methods to slot \code{objective}
if ( !isGeneric('method') ) setGeneric('method', function(object, ...) standardGeneric('method'))
setMethod('method', signature(object='NMFSeed'),
	function(object){						
		slot(object, 'method')
	}
)

if ( !isGeneric('method<-') ) setGeneric('method<-', function(object, ..., value) standardGeneric('method<-'))
setReplaceMethod('method', signature(object='NMFSeed', value='function'),
	function(object, value){
		slot(object, 'method') <- value
		validObject(object)
		object
	}
)

###########################################################################
# REGISTRY METHODS FOR SEEDING METHODS
###########################################################################

#' Register a new seeding method into the NMF registry.
#'
if ( is.null(getGeneric('nmfRegisterSeed')) ) setGeneric('nmfRegisterSeed', function(method, key, ...) standardGeneric('nmfRegisterSeed') )
setMethod('nmfRegisterSeed', signature(method='ANY', key='character'), 
		function(method, key, ...){	
			
			# wrap function method into a new NMFSeed object
			seedObj <- new('NMFSeed', name=key, method=method)
			# register the newly created object
			nmfRegister(seedObj, key, registry.name='seed', ...)
			
		}
)

#' Factory method to retrieve seeding methods from the NMF registry.
#'
nmfSeed <- function(name=NULL, ...){
	
	nmfGet(name, registry.name='seed', ...)
	
}

#' Returns TRUE if the algorithm is registered FALSE otherwise
existsNMFSeed <- function(name, exact=TRUE){	
	
	res <- !is.null( nmfGet(name, registry.name='seed', error=FALSE, exact=exact) )
	return(res)
	
}

###########################################################################
# INITIALIZATION FUNCTIONS
###########################################################################

#' Hook to initialize base seeding methods when the package is loaded  
.load.seed.base <- function(){
	
	# None [do nothing]
	nmfRegisterSeed(function(object, x, ...){object}, 'none', overwrite=TRUE)
	
	# Random
	nmfRegisterSeed(random, 'random', overwrite=TRUE)
		
}