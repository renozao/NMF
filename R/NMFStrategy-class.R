# S4 class for NMG algorithms
# 
# Author: Renaud Gaujoux
###############################################################################


#' @include algorithmic.R
#' @include NMFSet-class.R
NULL

#' Generic Strategy Class
#' 
#' This class defines a common interface for generic algorithm strategies 
#' (e.g., \code{\linkS4class{NMFStrategy}}).
#' 
#' @slot name character string giving the name of the algorithm
#' @slot package name of the package that defined the strategy.
#' @slot defaults default values for some of the algorithm's arguments.
#' 
#' @keywords internal
setClass('Strategy'
	, contains = 'VIRTUAL'
	, representation = representation(
		name = 'character' # the strategy name
		, package = 'character' # the package that defines the strategy
		, defaults = 'list'
	)
	, prototype = prototype(
		package = character()
		, name = character()
	)
	, validity=function(object){
		
		# slot 'name' must be a non-empty character string
		obj <- name(object)
		if( !length(obj) || (length(obj)>1L || obj=='') )
			return(str_c("Slot 'name' must be a single non-empty character string [", obj, ']'))
		TRUE
	}
)

#' Accessing Strategy Names 
#' 
#' \code{name} and \code{name<-} gets and sets the name associated with an object.
#' In the case of \code{Strategy} objects it is the the name of the algorithm.
#' 
#' @param object an R object with a defined \code{name} method
#' @param ... extra arguments to allow extension
#' @param value replacement value 
#' 
#' @export
#' @inline
#' @rdname Strategy-class
setGeneric('name', function(object, ...) standardGeneric('name'))
#' Returns the name of an algorithm
#' @param all a logical that indicates if all the names associated with a strategy 
#' should be returned (\code{TRUE}), or only the first (primary) one (\code{FALSE}).
setMethod('name', signature(object='Strategy'),
	function(object, all=FALSE){
		n <- slot(object, 'name')
		if( length(n) && !all ) n[1L] else n
	}
)
#' @export
#' @inline
#' @rdname Strategy-class 
setGeneric('name<-', function(object, ..., value) standardGeneric('name<-'))
#' Sets the name(s) of an NMF algorithm
setReplaceMethod('name', signature(object='Strategy', value='character'),
	function(object, value){
		slot(object, 'name') <- value
		validObject(object)
		object
	}
)

defaultArgument <- function(name, object, value, force=FALSE){

	# taken from methods::hasArg
	aname <- as.character(substitute(name))
	miss <- eval(substitute(missing(name)), sys.frame(sys.parent()))
	defaults <- attr(object, 'defaults')

	if( !miss && !force ) eval(substitute(name), sys.frame(sys.parent()))
	else if( aname %in% names(defaults) ) defaults[[aname]] 
	else value
}

#' Virtual Interface for NMF Algorithms
#' 
#' This class partially implements the generic interface defined for general 
#' algorithms defined in the \pkg{NMF} package (see \code{\link{algorithmic-NMF}}).
#'  
#' @slot objective the objective function associated with the algorithm (Frobenius, Kullback-Leibler, etc...). 
#'  It is either an access key of a registered objective function or a function definition. 
#'  In the latter case, the given function must have the following signature \code{(x="NMF", y="matrix")}
#'  and return a nonnegative real value.
#' 
#' @slot model a character string giving either the (sub)class name of the NMF-class instance used 
#' and returned by the strategy, or a function name.
#' 
#' @slot mixed a logical that indicates if the algorithm works on mixed-sign data. 
#' 
#' @keywords internal
setClass('NMFStrategy'
	, representation(
				objective = '.functionSlot' # the objective function used to compute the error (defined by name or function)
				, model = 'character' # NMF model to use
				, mixed = 'logical' # can the input data be negative?
	)
	, prototype=prototype(objective='euclidean', model='NMFstd', mixed=FALSE)
	, validity=function(object){
		
		# slot 'objective' must either be a non-empty character string or a function
		obj <- objective(object)
		if( is.character(obj) && obj == '' )
			return("Slot 'objective' must either be a non-empty character string or a function definition.")
			
		# slot 'model' must be the name of a class that extends class 'NMF'
		obj <- modelname(object)
		if( !is.character(obj) )
			return("Slot 'model' must be a character vector")
		if( any(inv <- !sapply(obj, isNMFclass)) )
			return(paste("Slot 'model' must contain only names of a class that extends class 'NMF' [failure on class(es) "
					, paste( paste("'", obj[inv], "'", sep=''), collapse=', ')  
					,"]"
					, sep=''))
		
		# slot 'mixed' must be a single logical		
		obj <- slot(object, 'mixed')
		if( length(obj) != 1 )
			return( paste("Slot 'mixed' must be a single logical [length=", length(obj), "]", sep='') )
	}
	, contains = c('VIRTUAL', 'Strategy')
)

#' @export
#' @rdname NMFStrategy-class
setMethod('show', 'NMFStrategy',
		function(object){			
			cat('<object of class: ', class(object), ">\n", sep='')
			cat(" name: ", name(object), " [", packageSlot(object), "]\n", sep='')
			svalue <- objective(object)
			svalue <- if( is.function(svalue) ) str_args(svalue, exdent=10) else paste("'", svalue,"'", sep='')
			cat(" objective:", svalue, "\n")
			cat(" model:", modelname(object), "\n")
			if( length(object@defaults) > 0L ){
				cat(" defaults:", str_desc(object@defaults, exdent=10L), "\n")
			}
			return(invisible())
		}
)

# Coerce method for 'NMFStrategy' objects into 'character': give the main name
setAs('NMFStrategy', 'character'
	, def = function(from) name(from)	
) 

#' Factory Method for NMFStrategy Objects
#' 
#' Creates NMFStrategy objects that wraps implementation of NMF algorithms into 
#' a unified interface.
#' 
#' @param name name/key of an NMF algorithm.
#' @param method definition of the algorithm
#' @param ... extra arguments passed to \code{\link{new}}.
#' 
#' @export
#' @inline
setGeneric('NMFStrategy', function(name, method, ...) standardGeneric('NMFStrategy') )
#' Creates an \code{NMFStrategyFunction} object that wraps the function \code{method}
#' into a unified interface.
#' 
#' \code{method} must be a function with signature \code{(y="matrix", x="NMFfit", ...)}, 
#' and return an object of class \code{\linkS4class{NMFfit}}.
setMethod('NMFStrategy', signature(name='character', method='function'), 
		function(name, method, ...){
			
			# build a NMFStrategyFunction object on the fly to wrap function 'method'
			NMFStrategy(name=name, algorithm=method, ...)
			
		}
)

#' Creates an \code{NMFStrategy} object based on a template object (Constructor-Copy).
setMethod('NMFStrategy', signature(name='character', method='NMFStrategy'), 
		function(name, method, ...){
			
			package <- topns_name()
			# build an NMFStrategy object based on template object
			strategy <- new(class(method), method, name=name, ..., package=package)
			
			# valid the new strategy
			validObject(strategy)
			
			# add trace of inheritance from parent NMF algorithm
			attr(strategy, 'parent') <- name(method)[1]
			
			# return new object
			strategy
		}
)

#' Creates an \code{NMFStrategy} based on a template object (Constructor-Copy), 
#' in particular it uses the \strong{same} name.
setMethod('NMFStrategy', signature(name='NMFStrategy', method='missing'), 
		function(name, method, ...){
			
			# do not change the object if single argument
			if( nargs() == 1L ) return(name)
			
			# use the name as a key
			# NB: need special trick to avoid conflict between argument and function 
			mname <- match.fun('name')(name)

			NMFStrategy(name=mname, method=name, ...)
		}
)

#' Creates an \code{NMFStrategy} based on a registered NMF algorithm that is used 
#' as a template (Constructor-Copy), in particular it uses the \strong{same} name.
#' 
#' It is a shortcut for \code{NMFStrategy(nmfAlgorithm(method, exact=TRUE), ...)}.
setMethod('NMFStrategy', signature(name='missing', method='character'), 
		function(name, method, ...){
			NMFStrategy(nmfAlgorithm(method, exact=TRUE), ...)
		}
)


#' Creates an \code{NMFStrategy} based on a template object (Constructor-Copy) 
#' but using a randomly generated name.
setMethod('NMFStrategy', signature(name='NULL', method='NMFStrategy'), 
		function(name, method, ...){
			
			# use the name as a key
			# NB: need special trick to avoid conflict between argument and function 
			mname <- match.fun('name')(method)
			mname <- basename(tempfile(str_c(mname, '_')))
			
			NMFStrategy(name=mname, method=method, ...)
		}
)

#' Creates an \code{NMFStrategy} based on a registered NMF algorithm that is used 
#' as a template.
setMethod('NMFStrategy', signature(name='character', method='character'), 
		function(name, method, ...){
			NMFStrategy(name=name, method=nmfAlgorithm(method, exact=TRUE), ...) 
		}
)
#' Creates an \code{NMFStrategy} based on a registered NMF algorithm (Constructor-Copy) 
#' using a randomly generated name.
#' 
#' It is a shortcut for \code{NMFStrategy(NULL, nmfAlgorithm(method), ...)}.
setMethod('NMFStrategy', signature(name='NULL', method='character'), 
		function(name, method, ...){
			NMFStrategy(NULL, method=nmfAlgorithm(method, exact=TRUE), ...) 
		}
)
#' Creates an NMFStrategy, determining its type from the extra arguments passed 
#' in \code{...}: if there is an argument named \code{Update} then an 
#' \code{NMFStrategyIterative} is created, or if there is an argument 
#' named \code{algorithm} then an \code{NMFStrategyFunction} is created.
#' Calls other than these generates an error.
#'  
setMethod('NMFStrategy', signature(name='character', method='missing'), 
		function(name, method, ...){
			
			package <- topns_name()
			# check iterative strategy
			if( hasArg2('Update') ){ # create a new NMFStrategyIterative object
				new('NMFStrategyIterative', name=name, ..., package=package)
			}else if( hasArg2('mcode') ){
				new('NMFStrategyOctave', name=name, ..., package=package)
			}else if( hasArg2('algorithm') ){
				new('NMFStrategyFunction', name=name, ..., package=package)
			}else{
				stop('NMFStrategy - Could not infer the type of NMF strategy to instantiate.')
			}
			
		}
)

#' Pure virtual method defined for all NMF algorithms to ensure 
#' that a method \code{run} is defined by sub-classes of \code{NMFStrategy}.
#' 
#' It throws an error if called directly.
#' @rdname NMFStrategy
setMethod('run', signature(object='NMFStrategy', y='matrix', x='NMFfit'),
	function(object, y, x, ...){
		stop("NMFStrategy::run is a pure virtual method that should be overloaded in class '", class(object),"'.")
	}
)
#' Method to run an NMF algorithm directly starting from a given NMF model.
#' @rdname NMFStrategy
setMethod('run', signature(object='NMFStrategy', y='matrix', x='NMF'),
	function(object, y, x, ...){
		run(object, y, NMFfit(fit=x, seed='none', method=name(object)), ...)
	}
)

#' Computes the value of the objective function between the estimate \code{x}
#' and the target \code{y}.
#' 
#' @param x an NMF model that estimates \code{y}.
#' 
#' @inline
setMethod('deviance', 'NMFStrategy',
	function(object, x, y, ...){
		
		obj.fun <- slot(object, 'objective')
		
		# return the distance computed using the strategy's objective function
		if( !is.function(obj.fun) )
			deviance(x, y, method=obj.fun, ...)
		else # directly compute the objective function
			obj.fun(x, y, ...)
		
	}
)
		
#' Gets the objective function associated with an NMF algorithm.
#'  
#' It is used in \code{\link[=deviance,NMFStrategy-method]{deviance}} 
#' to compute the objective value for an NMF model with respect to 
#' a given target matrix. 
#' 
#' @export
#' @rdname NMFStrategy-class
setMethod('objective', 'NMFStrategy',
	function(object){
		slot(object, 'objective')
	}
)
#' Sets the objective function associated with an NMF algorithm, with a character string
#' that must be a registered objective function.
#' @export
#' @rdname NMFStrategy-class
setReplaceMethod('objective', signature(object='NMFStrategy', value='character'),
	function(object, value){
		#TODO: test for the existence of objective method
		slot(object, 'objective') <- value
		validObject(object)
		object
	}
)
#' Sets the objective function associated with an NMF algorithm, with a function
#' that computes the approximation error between an NMF model and a target matrix.
#' @export
#' @rdname NMFStrategy-class
setReplaceMethod('objective', signature(object='NMFStrategy', value='function'),
	function(object, value){
		slot(object, 'objective') <- value
		validObject(object)
		object
	}
)

#' Returns the model(s) that an NMF algorithm can fit.
#' 
#' @examples
#' # get the type of model(s) associated with an NMF algorithm
#' modelname( nmfAlgorithm('brunet') )
#' modelname( nmfAlgorithm('nsNMF') )
#' modelname( nmfAlgorithm('offset') )
#' 
setMethod('modelname', signature(object='NMFStrategy'),
	function(object){
		slot(object, 'model')
	}
)
#' \code{is.mixed} tells if an NMF algorithm works on mixed-sign data. 
#' @export
#' @rdname NMFStrategy-class
is.mixed <-	function(object){
	return( slot(object, 'mixed') )
}

#' Showing Arguments of NMF Algorithms
#' 
#' This function returns the extra arguments that can be passed
#' to a given NMF algorithm in call to \code{\link{nmf}}.
#' 
#' @param x algorithm specification
#' @param ... extra argument to allow extension
#' 
#' @export
nmfFormals <- function(x, ...){
	UseMethod('nmfFormals')
}

#' @S3method nmfFormals character
nmfFormals.character <- function(x, ...){
	s <- nmfAlgorithm(x)
	nmfFormals(s, ...)
}

#' @S3method nmfFormals NMFStrategy
nmfFormals.NMFStrategy <- function(x, ...){
	m <- getMethod('run', signature(object='NMFStrategy', y='matrix', x='NMFfit'))
	args <- allFormals(m)
	# prepend registered default arguments
	expand_list(x@defaults, args)
}

#' \code{nmfArgs} is a shortcut for \code{args(nmfWrapper(x))}, to 
#' display the arguments of a given NMF algorithm.
#' 
#' @rdname nmfFormals
#' @export
#' @examples 
#' 
#' # show arguments of an NMF algorithm
#' nmfArgs('brunet')
#' nmfArgs('snmf/r')
nmfArgs <- function(x){
	args(nmfWrapper(x))
}
