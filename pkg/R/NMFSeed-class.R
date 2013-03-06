#' @include registry.R
#' @include NMFStrategy-class.R
NULL

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
			method = 'function' # the method actual definition
		)
		, contains = 'Strategy'
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

#' \code{NMFSeed} is a constructor method that instantiate 
#' \code{\linkS4class{NMFSeed}} objects. 
#'
#' @param key access key as a single character string
#' @param method specification of the seeding method, as a function that takes 
#' at least the following arguments:
#' \describe{
#' \item{object}{uninitialised/empty NMF model, i.e. that it has 0 rows and 
#' columns, but has already the rank requested in the call to \code{\link{nmf}} 
#' or \code{\link{seed}}.}
#' \item{x}{target matrix}
#' \item{...}{extra arguments}
#' }
#'  
#' @export
#' @rdname setNMFSeed
#' @inline
setGeneric('NMFSeed', function(key, method, ...) standardGeneric('NMFSeed') )
#' Default method simply calls \code{\link{new}} with the same arguments. 
setMethod('NMFSeed', signature(key='character', method='ANY'), 
		function(key, method, ...){
			# wrap function method into a new NMFSeed object
			new('NMFSeed', name=key, method=method, ..., package=topns_name())
		}
)

#' Creates an \code{NMFSeed} based on a template object (Constructor-Copy), 
#' in particular it uses the \strong{same} name.
setMethod('NMFSeed', signature(key='NMFSeed', method='ANY'), 
		function(key, method, ...){
			
			# do not change the object if single argument
			if( nargs() == 1L ) return(key)
			
			# build an object based on template object
			new(class(method), key, method=method, ..., package=topns_name())
			
		}
)

