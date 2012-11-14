#' @include registry.R
#' @include algorithmic.R
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
#' should be return (\code{TRUE}), or only the first (primary) one (\code{FALSE}).
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
#' algorithms defined in the package NMF (see \code{\link{algorithmic-NMF}}).
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
			
			# build an NMFStrategy object based on template object
			strategy <- new(class(method), method, name=name, ..., package=topns_name())
			
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
			
			# check iterative strategy
			if( hasArg2('Update') ){ # create a new NMFStrategyIterative object
				new('NMFStrategyIterative', name=name, ..., package=topns_name())
			}else if( hasArg2('algorithm') ){
				new('NMFStrategyFunction', name=name, ..., package=topns_name())
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

###########################################################################
# REGISTRY METHODS FOR ALGORITHMS
###########################################################################

# create sub-registry for NMF algorithm
setPackageRegistry('algorithm', "NMFStrategy", description="NMF algorithms") 

#' Registering NMF Algorithms
#' 
#' Adds a new algorithm to the registry of algorithms that perform 
#' Nonnegative Matrix Factorization.
#'  
#' @inheritParams NMFStrategy
#' @param overwrite logical that indicates if any existing NMF method with the 
#' same name should be overwritten (\code{TRUE}) or not (\code{FALSE}), 
#' in which case an error is thrown.
#' @param verbose a logical that indicates if information about the registration 
#' should be printed (\code{TRUE}) or not (\code{FALSE}).
#' @param ... arguments passed to the factory function \code{\link{NMFStrategy}},
#' which instantiate the \code{\linkS4class{NMFStrategy}} object that is stored
#' in registry. 
#' 
#' @export
setNMFMethod <- function(name, method, ..., overwrite=isLoadingNamespace(), verbose=TRUE){
	
	library(pkgmaker)
	
	# build call to constructor
	call_const <- match.call(NMFStrategy)
	call_const[[1]] <- as.name('NMFStrategy')
	call_const$verbose <- NULL
	call_const$overwrite <- NULL
	# swap name and method if method is missing and name is a registered method
	if( missing(method) && !missing(name) && is.character(name) && existsNMFMethod(name) ){
		call_const$method <- name
		call_const$name <- NULL
	}
	# build the NMFStrategy object (in the parent frame to get the package slot right)
	e <- parent.frame()
	method <- eval(call_const, envir=e)

	parent.method <- attr(method, 'parent')
	key <- match.fun('name')(method)[1]

	if( verbose ){
		tmpl <- 
		if( !is.null(parent.method) && parent.method != key )
			stringr::str_c(" based on template '", parent.method, "'")
	
		pkg <- packageSlot(method)
		message("Registering NMF algorithm '", pkg, '::', key,"'", tmpl,"... ", appendLF=FALSE)
	}
	
	# add to the algorithm registry
	res <- nmfRegister(method, key, registry.name='algorithm'
					, overwrite=overwrite, verbose=verbose>1L)
	
	if( !is.null(res) && res > 0L ){
		if( verbose ) message( if(res == 1L) "OK" else "UPDATED" )
		invisible(method)
	}else{
		if( verbose ) message( "ERROR" )
		invisible(NULL)
	}
}

#' \code{nmfRegisterAlgorithm} is an alias to \code{setNMFMethod} for backward
#' compatibility.
#' 
#' @export 
#' @rdname setNMFMethod
nmfRegisterAlgorithm <- setNMFMethod


#' Registry for NMF Algorithms 
#' 
#' @name methods-NMF
#' @rdname registry-algorithm
#' @family regalgo Registry for NMF algorithms
NULL

#' Testing Compatibility of Algorithm and Models
#' 
#' \code{canFit} is an S4 generic that tests if an algorithm can 
#' fit a particular model.
#' 
#' @param x an object that describes an algorithm
#' @param y an object that describes a model
#' @param ... extra arguments to allow extension
#' 
#' @export
#' @inline
#' @family regalgo
setGeneric('canFit', function(x, y, ...) standardGeneric('canFit') )
#' Tells if an NMF algorithm can fit a given class of NMF models
#' 
#' @param exact for logical that indicates if an algorithm is considered able to fit 
#' only the models that it explicitly declares (\code{TRUE}), or if it should be
#' considered able to also fit models that extend models that it explicitly fits. 
#'    
setMethod('canFit', signature(x='NMFStrategy', y='character'),
	function(x, y, exact=FALSE){
		
		if( !exact ){
			
			# check for one model amongst all the models fittable by the strategy
			can <- if( length(mo <- modelname(x)) > 1 )
						sapply(mo, function(m) extends(y, m))
					else extends(y, mo)
			any(can)
			
		}else
			is.element(y, modelname(x))
	}
)
#' Tells if an NMF algorithm can fit the same class of models as \code{y}
setMethod('canFit', signature(x='NMFStrategy', y='NMF'),
		function(x, y, ...){
			canFit(x, modelname(y), ...)
		}
)
#' Tells if a registered NMF algorithm can fit a given NMF model
setMethod('canFit', signature(x='character', y='ANY'),
		function(x, y, ...){
			canFit(nmfAlgorithm(x), y, ...)
		}
)

#' \code{selectNMFMethod} tries to select an appropriate NMF algorithm that is 
#' able to fit a given the NMF model.
#' 
#' @param name name of a registered NMF algorithm
#' @param model class name of an NMF model, i.e. a class that inherits from class 
#' \code{\linkS4class{NMF}}.
#' @param load a logical that indicates if the selected algorithms should be loaded
#' into \code{NMFStrategy} objects
#' @param all a logical that indicates if all algorithms that can fit \code{model}
#' should be returned or only the default or first found.
#' @param quiet a logical that indicates if warnings or errors should be thrown 
#' in case of the selected algorithm is not the default algorithm.
#' 
#' @return \code{selectNMFMethod} returns a character vector or \code{NMFStrategy} objects, 
#' or NULL if no suitable algorithm was found.
#' 
#' @rdname registry-algorithm
#' 
selectNMFMethod <- function(name, model, load=FALSE, exact=FALSE, all=FALSE, quiet=FALSE){
	
	# lookup for an algorithm suitable for the given NMF model
	if( !isNMFclass(model) )
		stop("argument 'model' must be the name of a class that extends class 'NMF'")
	
	
	algo_list <- if( !missing(name) ){
			algo <- nmfAlgorithm(name)
			name(algo) 
		}else nmfAlgorithm()
	
	# lookup for all the algorithms that can fit the given model
	#NB: if only one model needs to be selected then first look for an exact fit as 
	# this would need to be done with exact=FALSE and TRUE anyways
	w <- sapply(algo_list, canFit, model, exact= if(all) exact else TRUE)	
	algo <- algo_list[w]
	
	# if no suitable algorithm was found, and an exact match is not required 
	# then look for other potential non-exact algorithms
	if( !all && !exact && length(algo) == 0 ){
		w <- sapply(algo_list, canFit, model, exact=FALSE)
		algo <- algo_list[w]
	}
	
	# return NULL if no algorithm was found
	if( length(algo) == 0L ){
		if( !quiet ) 
			stop("Could not find an NMF algorithm to fit model '", model, "'"
				, if( !missing(name) ) paste(" amongst ", str_out(algo_list, Inf)))
		return(NULL)
	}
			
	# if all=FALSE then try to choose the default algorithm if present in the list, or the first one
	res <- if( !all && length(algo) > 1L ){
		
		idx <- which( algo == nmf.getOption('default.algorithm') ) 
		if( !length(idx) ) idx <- 1L
		
		res <- algo[idx]
		if( !quiet ) 
			warning("Selected NMF algorithm '", res, "' amongst other possible algorithm(s): "
					, paste(paste("'", algo[-idx], "'", sep=''), collapse=", "))
		res
	}else # otherwise return all the algorithms
		algo
	
	# load the methods if required
	if( load ){
		if( length(res) > 1 ) sapply(res, nmfAlgorithm) else nmfAlgorithm(res)
	}
	else
		res	
}


#' \code{getNMFMethod} retrieves NMF algorithm objects from the registry.
#' 
#' @param ... extra arguments passed to \code{\link[pkgmaker]{regfetch}}.
#' 
#' @export
#' @rdname registry-algorithm
getNMFMethod <- function(...) nmfGet(..., registry.name='algorithm', msg='NMF algorithm')

#' Listing and Retrieving NMF Algorithms
#' 
#' \code{nmfAlgorithm} lists access keys or retrieves NMF algorithms that are 
#' stored in registry.
#' It allows to list 
#'  
#' @param name Access key. 
#' If not missing, it must be a single character string that is partially matched 
#' against the available algorithms in the registry.
#' In this case, if \code{all=FALSE} (default), then the algorithm is returned 
#' as an \code{NMFStrategy} object that can be directly passed to \code{\link{nmf}}.
#' An error is thrown if no matching algorithm is found.
#' 
#' If missing or \code{NULL}, then access keys of algorithms -- that 
#' match the criteria \code{version}, are returned.
#' This argument is assumed to be regular expression if \code{all=TRUE} or 
#' \code{version} is not \code{NULL}.
#' @param version version of the algorithm(s) to retrieve. 
#' Currently only value \code{'R'} is supported, which searched for plain R 
#' implementations. 
#' @param all a logical that indicates if all algorithm keys should be returned, 
#' including the ones from alternative algorithm versions (e.g. plain R 
#' implementations of algorithms, for which a version based on optimised 
#' C updates is used by default). 
#' @param ... extra arguments passed to \code{\link{getNMFMethod}} when \code{name} 
#' is not \code{NULL} and \code{all=FALSE}. It is not used otherwise.
#' 
#' @return an \code{\linkS4class{NMFStrategy}} object if \code{name} is not 
#' \code{NULL} and \code{all=FALSE}, or a named character vector that contains 
#' the access keys of the matching algorithms.
#' The names correspond to the access key of the primary algorithm: e.g. 
#' algorithm \sQuote{lee} has two registered versions, one plain R (\sQuote{.R#lee}) 
#' and the other uses optimised C updates (\sQuote{lee}), which will all get 
#' named \sQuote{lee}.
#' 
#' @export
#' @family regalgo
#' 
#' @examples 
#' 
#' # list all main algorithms 
#' nmfAlgorithm()
#' # list all versions of algorithms 
#' nmfAlgorithm(all=TRUE)
#' # list all plain R versions 
#' nmfAlgorithm(version='R')
#'  
nmfAlgorithm <- function(name=NULL, version=NULL, all=FALSE, ...){	
	
	# if one passes an NMFStrategy just returns it
	if( is(name, 'NMFStrategy') ) return(name)
	
	# force all=TRUE if type is provided
	if( !is.null(version) ) all <- TRUE
	
	# directly return the algorithm object if a key is supplied and all=FALSE
	if( !is.null(name) && !all ) return( getNMFMethod(name, ...) )
	
	# get all algorithms
	algo <- getNMFMethod(all=TRUE)
	# set names to match the primary key
	algo <- setNames(algo, sub("^\\.(.+#)?", '', algo))	
	# filter out hidden methods
	if( !all ) algo <- algo[!grepl("^\\.", algo)]
	# filter out methods not from the requested algorithm
	if( !is.null(name) ) algo <- algo[grepl(str_c("^", name), names(algo))]
	# filter out types
	if( !is.null(version)  ){
		type <- match.arg(version, c('R'))
		algo <- Filter( function(x) grepl(str_c("^\\.", version, '#'), x), algo)
	}
	
	# return the selected algorithm(s)
	algo
}

#' \code{existsNMFMethod} tells if an NMF algorithm is registered under the
#' 
#' @param exact a logical that indicates if the access key should be matched 
#' exactly (\code{TRUE}) or partially (\code{FALSE}).
#' 
#' @export
#' @rdname registry-algorithm  
existsNMFMethod <- function(name, exact=TRUE){	
	
	!is.null( getNMFMethod(name, error=FALSE, exact=exact) )
	
}

#' Wrapping NMF Algorithms
#' 
#' This function creates a wrapper function for calling the function \code{\link{nmf}} 
#' with a given NMF algorithm.
#' 
#' @param method Name of the NMF algorithm to be wrapped. 
#' It should be the name of a registered algorithm as returned by \code{\link{nmfAlgorithm}}, 
#' or an NMF algorithm object (i.e. an instance of \code{\linkS4class{NMFStrategy}}). 
#' @param ... extra named arguments that define default values for any arguments 
#' of \code{\link{nmf}} or the algorithm itself. 
#' @param .FIXED a logical that indicates if the default arguments defined in \code{...}
#' must be considered as fixed, i.e. that they are forced to have the defined values and cannot
#' be used in a call to the wrapper function, in which case, a warning about discarding them 
#' is thrown if they are used.
#' Non fixed arguments may have their value changed at call time, in which case it is honoured and 
#' passed to the \code{nmf} call.
#' 
#' \code{.FIXED} may also be a character vector that specifies which argument amongst \code{...}
#' should be considered as fixed.
#' @return a function with argument \code{...} and a set of default arguments defined 
#' in \code{...} in the call to \code{nmfWrapper}.
#' 
#' @seealso \code{\link{nmfAlgorithm}}, \code{\link{nmf}}
#' @keywords internal
#' @export
#' 
#' @examples 
#' 
#' # wrap Lee & Seung algorithm into a function
#' lee <- nmfWrapper('lee')
#' lee
#' 
#' # test on random data
#' x <- rmatrix(100,20)
#' res <- nmf(x, 3, 'lee', seed=12345)
#' res2 <- lee(x, 3, seed=12345)
#' nmf.equal(res, res2)
#' 
#' \dontshow{ stopifnot(nmf.equal(res, res2)) }
#' 
nmfWrapper <- function(method, ..., .FIXED=FALSE){
	
	# store original call
	.call <- match.call()
	
	# check that all arguments are named
	if( nargs() > 1L && any(names(.call)[-(1:2)]=='') )
		stop("Invalid call: all arguments must be named.")
	
	# store fixed arguments from default arguments
	.fixedargs <- 'method'
	.defaults <- names(.call)[-1L]
	if( length(.defaults) ){
		e <- parent.frame()
		for(n in .defaults ){
			.call[[n]] <- eval(.call[[n]], envir=e)
		}
		if( isTRUE(.FIXED) ) .fixedargs <- c(.fixedargs, .defaults)
		else if( is.character(.FIXED) ){
			.FIXED <- .FIXED[.FIXED %in% .defaults]
			.fixedargs <- c(.fixedargs, .FIXED)	
		}
	}

	# define wrapper function
	fwrap <- function(...){
		
		# check for fixed arguments passed in the call that need
		# to be discarded
		ca <- match.call()
		nm <- names(ca)[-1L]
		if( any(fnm <- nm %in% .fixedargs) ){
			warning("Discarding fixed arguments from wrapped call to ", .call[1L]
					, " [", str_out(nm[fnm], Inf), '].', immediate.=TRUE)
			ca <- ca[!c(FALSE, fnm)]
		}
		#
		
		# set values of default arguments
		defaults <- formals()[-1L] # skip first argument `...`
		.call <- expand_list(ca, defaults, .exact=FALSE)
		# change into a call to nmf
		.call[[1L]] <- as.name('nmf')
		.call <- as.call(.call)
		# eval
		e <- parent.frame()
		eval(.call, envir=e)
	}
	# add default arguments to signature
	if( length(.defaults) )
		formals(fwrap) <- c(formals(fwrap), as.list(.call[.defaults]))
	
	return( fwrap )
	
}
