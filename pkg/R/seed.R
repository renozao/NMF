#' @include registry.R
NULL

# create sub-registry for seeding methods
setPackageRegistry('seed', "NMFSeed"
				, description="Seeding methods for NMF algorithms")

#' Base class that defines the interface for NMF seeding methods.
#' 
#' This class implements a simple wrapper strategy object that defines a unified
#' interface to seeding methods, that are used to initialise NMF models before 
#' fitting them with any NMF algorithm. 
#' 
#' @slot name character string giving the name of the seeding strategy
#' @slot method workhorse function that implements the seeding strategy.
#' It must have signature \code{(object="NMF", x="matrix", ...)} and initialise
#' the NMF model \code{object} with suitable values for fitting the target 
#' matrix \code{x}.
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
			
			return(TRUE)
		}		
)

#' Show method for objects of class \code{NMFSeed} 
setMethod('show', 'NMFSeed',
		function(object){			
			cat('<object of class: ', class(object), ">\n")
			cat("name:\t", name(object), "\n")
			svalue <- algorithm(object)
			svalue <- if( is.function(svalue) ) '<function>' else paste("'", svalue,"'", sep='')
			cat("method:\t", svalue, "\n")
			return(invisible())
		}
)

#' Returns the name of the seeding method described by \code{object}.
setMethod('name', signature(object='NMFSeed'),
		function(object){
			slot(object, 'name')
		}
)
#' Sets the name of the seeding method described by \code{object}.
#' 
#' This method is not meant to be called directly
setReplaceMethod('name', signature(object='NMFSeed', value='character'),
	function(object, value){
		slot(object, 'name') <- value
		validObject(object)
		object
	}
)

#' Returns the workhorse function of the seeding method described by \code{object}. 
setMethod('algorithm', signature(object='NMFSeed'),
	function(object){						
		slot(object, 'method')
	}
)
#' Sets the workhorse function of the seeding method described by \code{object}.
setReplaceMethod('algorithm', signature(object='NMFSeed', value='function'),
	function(object, value){
		slot(object, 'method') <- value
		validObject(object)
		object
	}
)

###########################################################################
# REGISTRY METHODS FOR SEEDING METHODS
###########################################################################

###% Register a new seeding method into the NMF registry.
###%
setGeneric('NMFSeed', function(key, method, ...) standardGeneric('NMFSeed') )
setMethod('NMFSeed', signature(key='character', method='ANY'), 
	function(key, method, ...){
		# wrap function method into a new NMFSeed object
		new('NMFSeed', name=key, method=method, ...)
	}
)

setNMFSeed <- function(..., overwrite=FALSE, verbose=nmf.getOption('verbose')){
		
	# development/tracking trick 
	if( !isNamespaceLoaded('NMF') ) overwrite <- TRUE 
	lverbose <- # if not specified: always when loading or in dev mode
			if( missing(verbose) ) isLoadingNamespace() || !isNamespaceLoaded('NMF')
			else verbose 
	
	# wrap function method into a new NMFSeed object
	method <- NMFSeed(...)
	parent.method <- attr(method, 'parent')
	key <- name(method)[1]
	
	if( lverbose ){
		tmpl <- if( !is.null(parent.method) && parent.method != key )
					stringr::str_c(" based on template '", parent.method, "'")
		
		message("Registering NMF seeding method '", key,"'", tmpl,"... ", appendLF=FALSE)
	}
	
	# register the newly created object
	res <- nmfRegister(method, key, registry.name='seed'
					, overwrite=overwrite, verbose=verbose)
	
	if( !is.null(res) && res > 0L ){
		if( lverbose ) message( if(res == 1L) "OK" else "UPDATED" )
		method
	}else{
		if( lverbose ) message( "ERROR" )
		NULL
	}
	
}

nmfRegisterSeed <- setNMFSeed

#' Seeding Strategies for NMF Algorithms
#' 
#' \code{nmfSeed} lists and retrieves NMF seeding methods.
#' 
#' Currently the internal registry contains the following seeding methods, 
#' which may be specified to the function \code{\link{nmf}} via its argument 
#' \code{seed} using their access keys:
#' 
#' \describe{
#' \item{random}{ The entries of each factors are drawn from a uniform 
#' distribution over \eqn{[0, max(x)]}, where $x$ is the target matrix.}
#' \item{nndsvd}{ Nonnegative Double Singular Value Decomposition.
#' 
#' The basic algorithm contains no randomization and is based on two SVD processes, 
#' one approximating the data matrix, the other approximating positive sections 
#' of the resulting partial SVD factors utilising an algebraic property of 
#' unit rank matrices.
#' 
#' It is well suited to initialise NMF algorithms with sparse factors.
#' Simple practical variants of the algorithm allows to generate dense factors.
#' 
#' \strong{Reference:} \cite{Boutsidis2008}}
#' \item{ica}{ Uses the result of an Independent Component Analysis (ICA) 
#' (from the \code{fastICA} package).
#' Only the positive part of the result are used to initialise the factors.}
#' \item{none}{ Fixed seed.
#' 
#' This method allows the user to manually provide initial values for 
#' both matrix factors.}
#' }
#' 
#' @param name access key of a seeding method stored in registry.
#' If missing, \code{nmfSeed} returns the list of all available seeding methods.
#' @param ... extra arguments used for internal calls
#'  
#' @export
#' 
#' @examples
#' 
#' # list all registered seeding methods
#' nmfSeed()
#' # retrieve one of the methods
#' nmfSeed('ica') 
#' 
nmfSeed <- function(name=NULL, ...){
	
	nmfGet(name, registry.name='seed', ...)
	
}

#' \code{existsNMFSeed} tells if a given seeding method exists in the registry.
#' 
#' @param exact a logical that indicates if the access key should be matched 
#' exactly or partially.
#'  
#' @rdname nmfSeed
#' @export
existsNMFSeed <- function(name, exact=TRUE){	
	
	res <- !is.null( nmfGet(name, registry.name='seed', error=FALSE, exact=exact) )
	return(res)
	
}

###########################################################################
# REGISTRATION
###########################################################################

## Register base seeding methods
# None: do nothing and return object unchanged
setNMFSeed('none', function(object, x, ...){object}, overwrite=TRUE)
# Random: use function rnmf
setNMFSeed('random', rnmf, overwrite=TRUE)
