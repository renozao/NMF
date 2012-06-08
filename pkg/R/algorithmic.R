# Definition of a generic interface for algorithms
# 
# Author: Renaud Gaujoux
###############################################################################


#' Generic Interface for Algorithms
#' 
#' @description
#' The functions documented here are S4 generics that define an general interface for 
#' -- optimisation -- algorithms.
#' 
#' This interface builds upon the broad definition of an algorithm as a workhorse function
#' to which is associated auxiliary objects such as an underlying model or an objective function 
#' that measures the adequation of the model with observed data.
#' It aims at complementing the interface provided by the \code{\link{stats}} package. 
#' 
#' @section Interface fo NMF algorithms:
#' This interface is implemented for NMF algorithms by the classes \code{\link{NMFfit}}, 
#' \code{\link{NMFfitX}} and \code{\link{NMFStrategy}}, and their respective sub-classes.
#' The examples given in this documentation page are mainly based on this implementation.
#' 
#' @param object an object computed using some algorithm, or that describes an algorithm
#' itself.
#' @param value replacement value
#' @param ... extra arguments to allow extension 
#' 
#' @name algorithmic-NMF
#' @rdname algorithmic
NULL

#' @details 
#' \code{algorithm} and \code{algorithm<-} get/set an object that describes the 
#' algorithm used to compute another object, or with which it is associated.
#' It may be a simple character string that gives the algorithm's names, or an object that
#' includes the algorithm's definition itself (e.g. an \code{\link{NMFStrategy}} object).
#' 
#' @export
#' @rdname algorithmic
setGeneric('algorithm', function(object, ...) standardGeneric('algorithm') )
#' @export
#' @rdname algorithmic
setGeneric('algorithm<-', function(object, ..., value) standardGeneric('algorithm<-') )

#' @details
#' \code{seeding} get/set the seeding method used to initialise the computation of an object, 
#' i.e. usually the function that sets the starting point of an algorithm.
#' 
#' @export
#' @rdname algorithmic
setGeneric('seeding', function(object, ...) standardGeneric('seeding') )
#' @export
#' @rdname algorithmic
setGeneric('seeding<-', function(object, ..., value) standardGeneric('seeding<-') )

#' @details
#' \code{niter} and \code{niter<-} get/set the number of iterations performed 
#' to compute an object.
#' The function \code{niter<-} would usually be called just before returning the result 
#' of an algorithm, when putting together data about the fit. 
#' 
#' @export
#' @rdname algorithmic
setGeneric('niter', function(object, ...) standardGeneric('niter'))
#' @rdname algorithmic
#' @export
setGeneric('niter<-', function(object, ..., value) standardGeneric('niter<-'))


#' @details
#' \code{nrun} returns the number of times the algorithm has been run to compute
#' an object.
#' Usually this will be 1, but may be be more if the algorithm involves multiple
#' starting points.
#'  
#' @export
#' @rdname algorithmic
setGeneric('nrun', function(object, ...) standardGeneric('nrun') )
#' Default method that returns the value of attribute \sQuote{nrun}.
#' 
#' Such an attribute my be attached to objects to keep track of data about 
#' the parent fit object (e.g. by method \code{\link{consensus}}), which 
#' can be used by subsequent function calls such as plot functions 
#' (e.g. see \code{\link{consensusmap}}). 
#' This method returns \code{NULL} if no suitable data was found.
setMethod('nrun', 'ANY', 
	function(object){
		attr(object, 'nrun')
	}
)

#' @details
#' \code{objective} and \code{objective<-} get/set the objective function associated 
#' with an object. 
#' Some methods for \code{objective} may also compute the objective value with respect to 
#' some target/observed data.
#' 
#' @export
#' @rdname algorithmic 
setGeneric('objective', function(object, ...) standardGeneric('objective'))
#' @export
#' @rdname algorithmic
setGeneric('objective<-', function(object, ..., value) standardGeneric('objective<-'))

#' @details
#' \code{runtime} returns the CPU time required to compute an object.
#' This would generally be an object of class \code{\link[=proc.time]{proc_time}}.
#' 
#' @export
#' @rdname algorithmic
setGeneric('runtime', function(object, ...) standardGeneric('runtime') )
#' @details
#' \code{runtime.all} returns the CPU time required to compute a collection of 
#' objects, e.g. a sequence of independent fits.   
#' 
#' @export
#' @rdname algorithmic
setGeneric('runtime.all', function(object, ...) standardGeneric('runtime.all') )

#' @details
#' \code{seqtime} returns the sequential CPU time -- that would be -- required 
#' to compute a collection of objects.
#' It would differ from \code{runtime.all} if the computations were performed 
#' in parallel.   
#' 
#' @export
#' @rdname algorithmic
setGeneric('seqtime', function(object, ...) standardGeneric('seqtime') )

#' @details
#' \code{modelname} returns a the type of model associated with an object.
#' 
#' @rdname algorithmic
#' @export
setGeneric('modelname', function(object, ...) standardGeneric('modelname'))

#' Default method which returns the class name(s) of \code{object}.
#' This should work for objects representing models on their own. 
#' 
#' For NMF objects, this is the type of NMF model, that corresponds to the 
#' name of the S4 sub-class of \code{\linkS4class{NMF}}, inherited by \code{object}.
#' 
#' @examples
#' # get the type of an NMF model
#' modelname(nmfModel(3))
#' modelname(nmfModel(3, model='NMFns'))
#' modelname(nmfModel(3, model='NMFOffset'))
#' 
setMethod('modelname', 'ANY', 
	function(object)
	{
		as.character(class(object))
	}
)

#' @details
#' \code{run} calls the workhorse function that actually implements a strategy/algorithm,  
#' and run it on some data object. 
#'
#' @param y data object, e.g. a target matrix
#' @param x a model object used as a starting point by the algorithm, 
#' e.g. a non-empty NMF model. 
#' 
#' @export
#' @rdname algorithmic
setGeneric('run', function(object, y, x, ...) standardGeneric('run'))

#' @details
#' \code{logs} returns the log messages output during the computation of an 
#' object.
#' @export
#' @rdname algorithmic
setGeneric('logs', function(object, ...) standardGeneric('logs'))
#' Default method that returns the value of attribute/slot \code{'logs'} or, if this latter  
#' does not exists, the value of element \code{'logs'} if \code{object} is a \code{list}.
#' It returns \code{NULL} if no logging data was found.
setMethod('logs', 'ANY', 
	function(object)
	{
		res <- attr(object, 'logs')
		if( !is.null(res) ) res 
		else if( is.list(object) ) object$logs
	}
)

#' @details
#' \code{compare} compares objects obtained from running separate algorithms.
#' 
#' @export
#' @rdname algorithmic
setGeneric('compare', function(object, ...) standardGeneric('compare') )
