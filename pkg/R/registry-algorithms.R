# NMF algorithm registry access methods
# 
# Author: Renaud Gaujoux
###############################################################################

#' @include registry.R
#' @include NMFStrategy-class.R
#' @include NMFStrategyFunction-class.R
#' @include NMFStrategyIterative-class.R
#' @include NMFStrategyOctave-class.R
NULL

# create sub-registry for NMF algorithm
.registryAlgorithm <- setPackageRegistry('algorithm', "NMFStrategy"
		, description = "Algorithms to solve MF optimisation problems"
		, entrydesc = "NMF algorithm") 

nmfAlgorithmInfo <- function(show=TRUE){
    obj <- .registryAlgorithm
    if( show ) print(obj)
    invisible(obj)
}

# specific register method for registering NMFStrategy objects
setMethod('nmfRegister', signature(key='NMFStrategy', method='missing'), 
		function(key, method, ...){
			nmfRegister(name(key), key, ..., regname='algorithm')
		}
)

#' Registering NMF Algorithms
#' 
#' Adds a new algorithm to the registry of algorithms that perform 
#' Nonnegative Matrix Factorization.
#'  
#' @inheritParams NMFStrategy
#' @param ... arguments passed to the factory function \code{\link{NMFStrategy}},
#' which instantiate the \code{\linkS4class{NMFStrategy}} object that is stored
#' in registry. 
#' @param overwrite logical that indicates if any existing NMF method with the 
#' same name should be overwritten (\code{TRUE}) or not (\code{FALSE}), 
#' in which case an error is thrown.
#' @param verbose a logical that indicates if information about the registration 
#' should be printed (\code{TRUE}) or not (\code{FALSE}).
#' 
#' @export
#' @examples 
#' 
#' # define/regsiter a new -- dummy -- NMF algorithm with the minimum arguments
#' # y: target matrix
#' # x: initial NMF model (i.e. the seed)
#' # NB: this algorithm simply return the seed unchanged 
#' setNMFMethod('mynmf', function(y, x, ...){ x })
#' 
#' # check algorithm on toy data
#' res <- nmfCheck('mynmf')
#' # the NMF seed is not changed
#' stopifnot( nmf.equal(res, nmfCheck('mynmf', seed=res)) ) 
#' 
setNMFMethod <- function(name, method, ..., overwrite=isLoadingNamespace(), verbose=TRUE){
		
	# build call to NMFStrategy constructor
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
	# add to the algorithm registry
	res <- nmfRegister(method, overwrite=overwrite, verbose=verbose)
	# return wrapper function invisibly
	wrap <- nmfWrapper(method)
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
#' @param quiet a logical that indicates if the operation should be performed quietly, 
#' without throwing errors or warnings.
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
#' @param ... extra arguments passed to \code{\link[pkgmaker]{pkgreg_fetch}}
#' or \code{\link[pkgmaker]{pkgreg_remove}}.
#' 
#' @export
#' @rdname registry-algorithm
getNMFMethod <- function(...) nmfGet('algorithm', ...)

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
	
	# remove names if no arguments
	if( is.null(version) ) algo <- setNames(algo, NULL)
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

#' \code{removeNMFMethod} removes an NMF algorithm from the registry.
#' 
#' @export
#' @rdname registry-algorithm
removeNMFMethod <- function(name, ...){
	pkgreg_remove('algorithm', key=name, ...)
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
#' lee <- nmfWrapper('lee', seed=12345)
#' args(lee)
#' 
#' # test on random data
#' x <- rmatrix(100,20)
#' res <- nmf(x, 3, 'lee', seed=12345)
#' res2 <- lee(x, 3)
#' nmf.equal(res, res2)
#' res3 <- lee(x, 3, seed=123)
#' nmf.equal(res, res3)
#' 
#' \dontshow{ 
#' stopifnot(nmf.equal(res, res2))
#' stopifnot( !nmf.equal(res, res3)) 
#' }
#' 
#' # argument 'method' has no effect
#' res4 <- lee(x, 3, method='brunet')
#' nmf.equal(res, res4)
#' 
#' \dontshow{ 
#' stopifnot(nmf.equal(res, res4))
#' }
#' 
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
	.defaults <- .defaults[!.defaults %in% 'method']
	if( length(.defaults) ){
#		e <- parent.frame()
#		for(n in .defaults){
#			.call[[n]] <- eval(.call[[n]], envir=e)
#		}
		if( isTRUE(.FIXED) ) .fixedargs <- c(.fixedargs, .defaults)
		else if( is.character(.FIXED) ){
			.FIXED <- .FIXED[.FIXED %in% .defaults]
			.fixedargs <- c(.fixedargs, .FIXED)	
		}
	}
	# store in local environment
	.method <- method
	
	.checkArgs <- function(ca, args){
		# check for fixed arguments passed in the call that need
		# to be discarded
		nm <- names(ca)[-1L]
		if( any(fnm <- !is.na(pmatch(nm, .fixedargs))) ){
			warning("Discarding fixed arguments from wrapped call to ", .call[1L]
					, " [", str_out(nm[fnm], Inf), '].', immediate.=TRUE)
			ca <- ca[!c(FALSE, fnm)]
		}
		#
		
		# start with complete call
		.call <- ca
		# set values of wrapper default arguments if any
		if( length(.defaults) ){
			defaults <- args[.defaults]
			.call <- expand_list(ca, defaults, .exact=FALSE)
		}
		# change into a call to nmf
		.call[[1L]] <- as.name('nmf')
		.call[['method']] <- force(.method)
		as.call(.call)
	}
	
	# define wrapper function
	fwrap <- function(...){
		ca <- match.call()
		args <- formals()
		.call <- .checkArgs(ca, args)
		# eval in parent environment
		e <- parent.frame()
		eval(.call, envir=e)
	}
	
	# add default arguments to signature
	if( length(.defaults) ){
		formals(fwrap) <- expand_list(formals(fwrap), as.list(.call[.defaults]))
	}
	# add arguments from the NMF algorithm
	if( length(meth <- nmfFormals(.method)) ){
		formals(fwrap) <- expand_list(formals(fwrap), meth)
	}
	
	return( fwrap )
	
}

