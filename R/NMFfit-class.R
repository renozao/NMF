# Implementation of class NMFfit 
# 
# This class manages the result of a single run of a NMF algorithm.
# 
# Author: Renaud Gaujoux
###############################################################################

#' @include fixed-terms.R
#' @include nmfModel.R
NULL

#' Base Class for to store Nonnegative Matrix Factorisation results
#' 
#' Base class to handle the results of general \strong{Nonnegative Matrix
#' Factorisation} algorithms (NMF).
#' 
#' It provides a general structure and generic functions to manage the results
#' of NMF algorithms.  It contains a slot with the fitted NMF model (see slot
#' \code{fit}) as well as data about the methods and parameters used to compute
#' the factorization.
#' 
#' The purpose of this class is to handle in a generic way the results of NMF
#' algorithms. Its slot \code{fit} contains the fitted NMF model as an object
#' of class \code{\linkS4class{NMF}}.
#' 
#' Other slots contains data about how the factorization has been computed,
#' such as the algorithm and seeding method, the computation time, the final
#' residuals, etc\dots{}
#' 
#' Class \code{NMFfit} acts as a wrapper class for its slot \code{fit}.  It
#' inherits from interface class \code{\linkS4class{NMF}} defined for generic
#' NMF models.  Therefore, all the methods defined by this interface can be
#' called directly on objects of class \code{NMFfit}. The calls are simply
#' dispatched on slot \code{fit}, i.e.  the results are the same as if calling
#' the methods directly on slot \code{fit}.
#' 
#' @slot fit An object that inherits from class \code{\linkS4class{NMF}}, and 
#' contains the fitted NMF model.
#' 
#' NB: class \code{NMF} is a virtual class. The default class for this
#' slot is \code{NMFstd}, that implements the standard NMF model.
#' 
#' @slot residuals A \code{numeric} vector that contains the final
#' residuals or the residuals track between the target matrix and its NMF
#' estimate(s).  Default value is \code{numeric()}.
#' 
#' See method \code{\link{residuals}} for details on accessor methods and main
#' interface \code{\link{nmf}} for details on how to compute NMF with residuals
#' tracking.
#' 
#' @slot method a single \code{character} string that contains the
#' name of the algorithm used to fit the model. 
#' Default value is \code{''}.
#' 
#' @slot seed a single \code{character} string that contains the
#' name of the seeding method used to seed the algorithm that fitted the NMF 
#' model.
#' Default value is \code{''}.  See \code{\link{nmf}} for more details.
#' 
#' @slot rng an object that contains the RNG settings used for the
#' fit.  
#' Currently the settings are stored as an integer vector, the value of
#' \code{\link{.Random.seed}} at the time the object is created.  
#' It is initialized by the \code{initialized} method.  
#' See \code{\link{getRNG}} for more details.
#' 
#' @slot distance either a single \code{"character"} string that
#' contains the name of the built-in objective function, or a \code{function}
#' that measures the residuals between the target matrix and its NMF estimate.
#' See \code{\link{objective}} and \code{\link{deviance,NMF-method}}.
#' 
#' @slot parameters a \code{list} that contains the extra parameters
#' -- usually specific to the algorithm -- that were used to fit the model.
#' 
#' @slot runtime object of class \code{"proc_time"} that contains
#' various measures of the time spent to fit the model.  
#' See \code{\link[base]{system.time}}
#' 
#' @slot options a \code{list} that contains the options used to
#' compute the object.
#' 
#' @slot extra a \code{list} that contains extra miscellaneous data
#' for internal usage only.  
#' For example it can be used to store extra parameters or temporary data, 
#' without the need to explicitly extend the \code{NMFfit} class. 
#' Currently built-in algorithms only use this slot to
#' store the number of iterations performed to fit the object.
#'   
#' Data that need to be easily accessible by the end-user should rather be set
#' using the methods \code{$<-} that sets elements in the \code{list} slot 
#' \code{misc} -- that is inherited from class \code{\linkS4class{NMF}}.
#' 
#' @slot call stored call to the last \code{nmf} method that generated the
#' object. 
#' 
#' @export
#' @examples 
#' # run default NMF algorithm on a random matrix
#' n <- 50; r <- 3; p <- 20
#' V <- rmatrix(n, p)  
#' res <- nmf(V, r)							
#' 
#' # result class is NMFfit
#' class(res)
#' isNMFfit(res)
#' 
#' # show result
#' res
#' 
#' # compute summary measures
#' summary(res, target=V)
#' 
setClass('NMFfit'
	, representation(
			fit = 'NMF', # NMF model
			residuals = 'numeric', # residuals from the target matrix
			method = 'character', # method used to compute the factorization
			seed = 'character', # seeding method used to compute the factorization
			rng = 'ANY', # numerical random seed
			distance = '.functionSlotNULL', # method used to compute the distance between the target matrix and its NMF estimate
			parameters = 'list', # method used to compute the factorization
			runtime = 'proc_time', # running time to perform the NMF
			options = 'list', # run options
			extra = 'list' # extra list of results output by the method
			, call = 'call' # store last call to nmf()
	)
	
	, prototype = prototype(
			residuals = numeric(),
			method = '',
			seed = '',
			parameters = list(),
			extra = list()
	)
	
	, validity = function(object){
		
		# slot 'objective' must either be a non-empty character string or a function
		obj <- objective(object)
		if( is.character(obj) && obj == '')
			return(paste("Slot 'objective' must either be a non-empty character string or a function definition", sep=''))
		
		
		# everything went fine: return TRUE
		TRUE
	}
	, contains = 'NMF'
)

#' The function \code{NMFfit} is a factory method for NMFfit objects, that should
#' not need to be called by the user.
#' It is used internally by the functions \code{\link{nmf}} and \code{seed} to 
#' instantiate the starting point of NMF algorithms.
#' 
#' @param fit an NMF model
#' @param ... extra argument used to initialise slots in the instantiating 
#' \code{NMFfit} object.
#' @param rng RNG settings specification (typically a suitable value for 
#' \code{\link{.Random.seed}}).   
#' 
#' @rdname NMFfit-class
NMFfit <- function(fit=nmfModel(), ..., rng=NULL){
				
		# use current RNG settings if not otherwise provided
		if( is.null(rng) )
			rng <- getRNG()
		
		new('NMFfit', fit=fit, ..., rng=rng)
}

#' Computes and return the estimated target matrix from an NMF model fitted with 
#' function \code{\link{nmf}}.
#' 
#' It is a shortcut for \code{fitted(fit(object), ...)}, dispatching the call to 
#' the \code{fitted} method of the actual NMF model.  
setMethod('fitted', signature(object='NMFfit'),
	function(object, ...){
		fitted(fit(object), ...)
	}
)

#' Returns the basis matrix from an NMF model fitted with 
#' function \code{\link{nmf}}.
#' 
#' It is a shortcut for \code{.basis(fit(object), ...)}, dispatching the call to 
#' the \code{.basis} method of the actual NMF model.
setMethod('.basis', signature(object='NMFfit'),
	function(object, ...){
		.basis(fit(object), ...)
	}
)
#' Sets the the basis matrix of an NMF model fitted with 
#' function \code{\link{nmf}}.
#' 
#' It is a shortcut for \code{.basis(fit(object)) <- value}, dispatching the call to 
#' the \code{.basis<-} method of the actual NMF model.
#' It is not meant to be used by the user, except when developing 
#' NMF algorithms, to update the basis matrix of the seed object before 
#' returning it.
#' 
setReplaceMethod('.basis', signature(object='NMFfit', value='matrix'), 
	function(object, value){ 
		.basis(fit(object)) <- value
		object
	} 
)

#' Returns the the coefficient matrix from an NMF model fitted with 
#' function \code{\link{nmf}}.
#' 
#' It is a shortcut for \code{.coef(fit(object), ...)}, dispatching the call to 
#' the \code{.coef} method of the actual NMF model.
setMethod('.coef', signature(object='NMFfit'),
	function(object, ...){
		.coef(fit(object), ...)
	}
)
#' Sets the the coefficient matrix of an NMF model fitted with 
#' function \code{\link{nmf}}.
#' 
#' It is a shortcut for \code{.coef(fit(object)) <- value}, dispatching the call to 
#' the \code{.coef<-} method of the actual NMF model.
#' It is not meant to be used by the user, except when developing 
#' NMF algorithms, to update the coefficient matrix in the seed object before
#' returning it.
#' 
setReplaceMethod('.coef', signature(object='NMFfit', value='matrix'), 
	function(object, value){ 
		.coef(fit(object)) <- value
		object
	} 
)

#' Method for single NMF fit objects, which returns the indexes of fixed
#' basis terms from the fitted model.  
setMethod('ibterms', 'NMFfit', 
	function(object){
		ibterms(fit(object))
	}
)
#' Method for single NMF fit objects, which returns the indexes of fixed 
#' coefficient terms from the fitted model.
setMethod('icterms', 'NMFfit', 
	function(object){
		icterms(fit(object))
	}
)


#' Returns the offset from the fitted model. 
setMethod('offset', signature(object='NMFfit'), 
	function(object){
		offset(fit(object))
	}
)

#' Returns the number of iteration performed to fit an NMF model, typically 
#' with function \code{\link{nmf}}.
#' 
#' Currently this data is stored in slot \code{'extra'}, but this might change 
#' in the future.
setMethod('niter', signature(object='NMFfit'),
	function(object, ...){
		object@extra$iteration
	}
)
#' Sets the number of iteration performed to fit an NMF model.
#' 
#' This function is used internally by the function \code{\link{nmf}}.
#' It is not meant to be called by the user, except when developing 
#' new NMF algorithms implemented as single function, to set the number 
#' of iterations performed by the algorithm on the seed, before returning it 
#' (see \code{\linkS4class{NMFStrategyFunction}}).
#' 
setReplaceMethod('niter', signature(object='NMFfit', value='numeric'), 
	function(object, value){
		if( (length(value) != 1) || value < 0  ) 
			stop("NMF::niter - invalid value for 'niter': single non-negative value is required.", call.=FALSE) 
		object@extra$iteration <- value
		object
	} 
)

#' Show method for objects of class \code{NMFfit}
setMethod('show', 'NMFfit', 
	function(object)
	{		
		cat("<Object of class: ", class(object), ">\n", sep='')
		cat(" # Model:\n  ")
		s <- capture.output(show(fit(object)))
		cat(s, sep="\n  ")
		cat(" # Details:\n  ")
		.local <- function(){
			if( algorithm(object) != '' ) cat("algorithm: ", algorithm(object), "\n")
			if( seeding(object) != '' ) cat("seed: ",  seeding(object), "\n")
			
			# initial RNG stream			
			cat("RNG: ", RNGstr(object), "\n", sep='')
	
			# distance/objective function
			svalue <- objective(object)
			svalue <- if( is.function(svalue) ) '<function>' else paste("'", svalue,"'", sep='')
			cat("distance metric: ", svalue, "\n")			
			if( length(residuals(object)) !=0  ) cat("residuals: ",  residuals(object), "\n");
			# show the miscellaneous result values
			if( length(object@misc) > 0L )
				cat("miscellaneous:", str_desc(object@misc, exdent=12L), ". (use 'misc(object)')\n")
			# show the parameters specific to the method		
			if( length(object@parameters) > 0 ){
				cat("parameters:", str_desc(object@parameters, exdent=12L), "\n")
#				p <- sapply(object@parameters, function(x){
#					if( is.vector(x) && length(x) == 1L ) x
#					else paste("<", class(x), ">", sep='')
#				})
#				cat(str_wrap(str_out(p, NA, use.names=TRUE, quote=FALSE), exdent=12), "\n")
			}
			# show number of iterations if present
			if( !is.null(i <- niter(object)) ) cat("Iterations:", i, "\n")
			# show elapsed time if present
			if( length(runtime(object)) > 0 ){ cat("Timing:\n"); show(runtime(object));}
		}
		s <- capture.output(.local())
		cat(s, sep="\n  ")
	}
)



#' Extracting Fitted Models
#' 
#' The functions \code{fit} and \code{minfit} are S4 genetics that extract 
#' the best model object and the best fit object respectively, from a collection 
#' of models or from a wrapper object. 
#' 
#' @details
#' A fit object differs from a model object in that it contains data about the 
#' fit, such as the initial RNG settings, the CPU time used, etc\ldots, while 
#' a model object only contains the actual modelling data such as regression 
#' coefficients, loadings, etc\ldots  
#' 
#' That best model is generally defined as the one that achieves the 
#' maximum/minimum some quantitative measure, amongst all models in a collection.
#' 
#' In the case of NMF models, the best model is the one that achieves the best
#' approximation error, according to the objective function associated with the 
#' algorithm that performed the fit(s).
#' 
#' @param object an object fitted by some algorithm, e.g. as returned by the 
#' function \code{\link{nmf}}.
#' @param value replacement value
#' @param ... extra arguments to allow extension
#'
#' @rdname fit
#' @export 
setGeneric('fit', function(object, ...) standardGeneric('fit'))
#' Returns the NMF model object stored in slot \code{'fit'}. 
setMethod('fit', 'NMFfit', function(object) slot(object, 'fit'))

#' \code{fit<-} sets the fitted model in a fit object.
#' It is meant to be called only when developing new NMF algorithms, e.g. to update 
#' the value of the model stored in the starting point. 
#' 
#' @rdname fit
#' @export
setGeneric('fit<-', function(object, value) standardGeneric('fit<-'))
#' Updates the NMF model object stored in slot \code{'fit'} with a new value.
setReplaceMethod('fit', signature(object='NMFfit', value='NMF'), 
		function(object, value){ 
			slot(object, 'fit') <- value		
			object # TODO: valid object before returning it (+param check=TRUE or FALSE)
		} 
)

#' @rdname fit
#' @export
setGeneric('minfit', function(object, ...) standardGeneric('minfit') )
#' Returns the object its self, since there it is the result of a single NMF run.
setMethod('minfit', 'NMFfit', function(object) object)


#' Returns the type of a fitted NMF model.
#' It is a shortcut for \code{modelname(fit(object)}. 
setMethod('modelname', signature(object='NMFfit'), 
	function(object)
	{
		modelname(fit(object))
	}
)

#' Residuals in NMF Models
#' 
#' The package NMF defines methods for the function \code{\link[stats]{residuals}}
#' that returns the final residuals of an NMF fit or the track of the residuals
#' along the fit process, computed according to the objective function 
#' associated with the algorithm that fitted the model.
#' 
#' When called with \code{track=TRUE}, the whole residuals track is returned, 
#' if available.
#' Note that method \code{\link{nmf}} does not compute the residuals track, 
#' unless explicitly required.
#' 
#' It is a S4 methods defined for the associated generic functions from package
#' \code{stats} (See \link[stats]{residuals}).
#' 
#' @note Stricly speaking, the method \code{residuals,NMFfit} does not fulfill 
#' its contract as defined by the package \code{stats}, but rather acts as function
#' \code{deviance}.  
#' The might be changed in a later release to make it behave as it should.
#' 
#' @param object an \code{NMFfit} object as fitted by function \code{\link{nmf}}, 
#' in single run mode.
#' @param ... extra parameters (not used) 
#' 
#' @return \code{residuals} returns a single numeric value if \code{track=FALSE} 
#' or a numeric vector containing the residual values at some iterations.
#' The names correspond to the iterations at which the residuals were computed.
#' 
#' @family stats
#' @inline
#' @rdname residuals
#' @export 
#' 
setGeneric('residuals', package='stats')
#' Returns the residuals -- track -- between the target matrix and the NMF 
#' fit \code{object}. 
#' 
#' @param track a logical that indicates if the complete track of residuals 
#' should be returned (if it has been computed during the fit), or only the last
#' value.
#' 
#' @param niter specifies the iteration number for which one wants 
#' to get/set/test a residual value. This argument is used only if not \code{NULL}
#' 
setMethod('residuals', 'NMFfit', 
	function(object, track=FALSE, niter=NULL, ...){ 
		## IMPORTANT: keep this '...' and do not add a 'method' argument as this
		## one is passed by NMFfitX::fit (see bug #159) and is not supposed to be 
		## used
		res <- slot(object, 'residuals')
		if( track ) res 
		else if( is.null(niter) ) tail(res, n=1)
		else res[as.character(niter)]
	} 
)

#' \code{residuals<-} sets the value of the last residuals, or, optionally, 
#' of the complete residual track.
#' 
#' @param value residual value
#' 
#' @export
#' @inline
#' @rdname residuals 
setGeneric('residuals<-', function(object, ..., value) standardGeneric('residuals<-') )
#' @inline
setReplaceMethod('residuals', 'NMFfit',
	function(object, ..., niter=NULL, track=FALSE, value){
		if( track ) slot(object, 'residuals') <- value
		else{
			if( !is.null(niter) ) value <- setNames(value, niter)
			slot(object, 'residuals') <- c(slot(object, 'residuals'), value)
		}
		object
	}
)

#' Tells if an \code{NMFfit} object contains a recorded residual track.
#' 
#' @export
#' @rdname residuals
hasTrack <- function(object, niter=NULL){
	if( is.null(niter) ) length( slot(object, 'residuals') ) > 1
	else !is.na(slot(object, 'residuals')[as.character(niter)])
}

#' \code{trackError} adds a residual value to the track of residuals.
#' 
#' @param force logical that indicates if the value should be added to the track
#' even if there already is a value for this iteration number or if the iteration 
#' does not conform to the tracking interval \code{nmf.getOption('track.interval')}.
#' 
#' @rdname residuals
#' @export
trackError <- function(object, value, niter, force=FALSE){	
	track <- run.options(object, 'error.track')
	track.interval <- run.options(object, 'track.interval')
	
	if( force || (track && niter %% track.interval == 0) ){
		# add the new value to the error track
		last.iter <- names(residuals(object))
		duplicate <- if( !is.null(last.iter) ) niter == last.iter else FALSE
		if( !duplicate ){
			iter <- if( niter >= 0 ) niter
			residuals(object, niter=iter) <- value
		}
	}
	object
}

#' Returns the deviance of a fitted NMF model.
#' 
#' This method returns the final residual value if the target matrix \code{y} is
#' not supplied, or the approximation error between the fitted NMF model stored 
#' in \code{object} and \code{y}.
#' In this case, the computation is performed using the objective function 
#' \code{method} if not missing, or the objective of the algorithm that 
#' fitted the model (stored in slot \code{'distance'}).
#' 
#' If not computed by the NMF algorithm itself, the value is automatically
#' computed at the end of the fitting process by the function \code{\link{nmf}}, 
#' using the objective function associated with the NMF algorithm, so that it 
#' should always be available.     
#' 
#' @inline 
setMethod('deviance', 'NMFfit',
	function(object, y, method, ...){
		
		if( missing(y) ) setNames(residuals(object), NULL)
		else{
			# if missing retrieve the actual distance measure from the NMF object
			if( missing(method) ) method = object@distance
			
			# compute the distance between the target and the fitted NMF model
			deviance(fit(object), y, method=method, ...)
		}
	}
)

#' Returns the name of the algorithm that fitted the NMF model \code{object}.
setMethod('algorithm', 'NMFfit', function(object){ object@method } )
#' @inline
setReplaceMethod('algorithm', 'NMFfit',
	function(object, value){
		object@method <- value
		object
	}
)

#' Returns the name of the seeding method that generated the starting point
#' for the NMF algorithm that fitted the NMF model \code{object}.
setMethod('seeding', 'NMFfit', function(object){ object@seed } )
#' @inline
setReplaceMethod('seeding', 'NMFfit',
	function(object, value){
		object@seed <- value
		object
	}
)

#' Returns the objective function associated with the algorithm that computed the 
#' fitted NMF model \code{object}, or the objective value with respect to a given 
#' target matrix \code{y} if it is supplied.
#' 
#' @param y optional target matrix used to compute the objective value.
#' 
setMethod('objective', signature(object='NMFfit'),
	function(object, y){
		
		# when both x and y are missing then returns slot objective
		if( missing(y) ) return(slot(object, 'distance'))
		
		# return the distance computed using the strategy's objective function
		deviance(fit(object), y, method=slot(object, 'distance'))
		
	}
)
#' @inline
setReplaceMethod('objective', signature(object='NMFfit', value='ANY'),
	function(object, value){
		slot(object, 'distance') <- value
		validObject(object)
		object
	}
)

#' Returns the CPU time required to compute a single NMF fit.
setMethod('runtime', 'NMFfit', 
	function(object, ...){ 
		object@runtime
	}
)

#' Identical to \code{runtime}, since their is a single fit. 
setMethod('runtime.all', 'NMFfit', getMethod('runtime', 'NMFfit'))

###% Access methods to run options.
setGeneric('run.options', function(object, ...) standardGeneric('run.options') )
setMethod('run.options', 'NMFfit', 
	function(object, name){
		if( missing(name) ) object@options
		else object@options[[name]]
	}
)
setGeneric('run.options<-', function(object, ..., value) standardGeneric('run.options<-') )
setReplaceMethod('run.options', 'NMFfit', 
	function(object, ..., value){
		
		params <- list(...)
		baseError <- 'Setting NMF runtime options: ' 
		if ( length(params) == 0 ){
			if( !is.list(value) ) stop(baseError, 'options must be given as a list')
			object@options <- value
			return(object)
		}
		else if ( length(params) > 1 ) stop(baseError, 'options cannot set more than one option at a time')
		name <- params[[1]]
		if( !is.character(name) ) stop(baseError, 'option name must be given as a character string')
		# check if the option exists
		#if( !is.element(name, names(nmf.options.runtime())) ) stop(baseError, "option '", name, "' is not defined.")
		
		object@options[[name]] <- value
		object
	}
)
setGeneric('verbose', function(object, ...) standardGeneric('verbose') )
setMethod('verbose', 'NMFfit', 
	function(object){
		return(run.options(object, 'verbose') || nmf.getOption('debug'))
	}
)

setGeneric('plot', package='graphics' )
#' Plots the residual track computed at regular interval during the fit of 
#' the NMF model \code{x}.
#' 
#' @param skip an integer that indicates the number of points to skip/remove from the beginning
#' of the curve.
#' If \code{skip=1L} (default) only the initial residual -- that is computed before any iteration, is
#' skipped, if present in the track (it associated with iteration 0).
#'
#' @export 
setMethod('plot', signature(x='NMFfit', y='missing'),
	function(x, y, skip=-1L, ...){
		
		# retrieve the residuals track
		track <- residuals(x, track=TRUE)
		if( length(track) <= 1 ){
			warning(class(x), ' object has no residuals track')
			return(invisible())
		}
		# skip part of the track
		if( skip == -1L && !is.null(names(track)) ) track <- track[names(track)!='0'] # remove initial residual
		else if( skip > 0 ) track <- track[-(1:skip)]
		
		# set default graphical parameters (those can be overriden by the user)
		params <- .set.list.defaults(list(...)
				, xlab='Iterations'
				, ylab=paste('Objective value ('
							, if( is.character(x@distance) ) x@distance else algorithm(x), ')'
							, sep='' )
				, main=paste("NMF Residuals\nMethod: ", algorithm(x), " - Rank: ", nbasis(x), sep='')
				, cex.main = 1
				, col='#5555ff', lwd=1.4, type='l', cex=0.5)
		
		do.call('plot', c(list(names(track), track), params))
		points(names(track), track, type='p', cex=0.6, col=params$col)
	}
)

#' Computes summary measures for a single fit from \code{\link{nmf}}. 
#' 
#' This method adds the following measures to the measures computed by the method 
#' \code{summary,NMF}:
#' 
#' \describe{
#' \item{residuals}{Residual error as measured by the objective function associated
#' to the algorithm used to fit the model.}
#' \item{niter}{Number of iterations performed to achieve convergence of the algorithm.}
#' \item{cpu}{Total CPU time required for the fit.}
#' \item{cpu.all}{Total CPU time required for the fit. For \code{NMFfit} objects, this element is 
#' always equal to the value in \dQuote{cpu}, but will be different for multiple-run fits.}
#' \item{nrun}{Number of runs performed to fit the model. This is always equal to 1 for 
#' \code{NMFfit} objects, but will vary for multiple-run fits.}
#' }
#' 
#' @inline
#' 
#' @examples
#' # generate a synthetic dataset with known classes: 50 features, 18 samples (5+5+8)
#' n <- 50; counts <- c(5, 5, 8);
#' V <- syntheticNMF(n, counts)
#' cl <- unlist(mapply(rep, 1:3, counts))
#' 
#' # perform default NMF with rank=2
#' x2 <- nmf(V, 2)
#' summary(x2, cl, V)
#' # perform default NMF with rank=2
#' x3 <- nmf(V, 3)
#' summary(x2, cl, V)
#' 
setMethod('summary', signature(object='NMFfit'), 
	function(object, ...){
		
		res <- summary(fit(object), ...)
		
		## IMPORTANT: if adding a summary measure also add it in the sorting 
		## schema of method NMFfitX::compare to allow ordering on it
		
		# retreive final residuals
		res <- c(res, residuals=as.numeric(residuals(object)))
		# nb of iterations
		res <- c(res, niter=as.integer(niter(object)) )
		# runtime
		t <- runtime(object)
		utime <- as.numeric(t['user.self'] + t['user.child'])
		res <- c(res, cpu=utime, cpu.all=utime, nrun=1)		
		
		# return result
		return(res)
	}
)

#' Compares two NMF models when at least one comes from a NMFfit object, 
#' i.e. an object returned by a single run of \code{\link{nmf}}. 
setMethod('nmf.equal', signature(x='NMFfit', y='NMF'), 
		function(x, y, ...){
			nmf.equal(fit(x), y, ...)
		}
)
#' Compares two NMF models when at least one comes from a NMFfit object, 
#' i.e. an object returned by a single run of \code{\link{nmf}}.
setMethod('nmf.equal', signature(x='NMF', y='NMFfit'), 
		function(x, y, ...){
			nmf.equal(x, fit(y), ...)
		}
)
#' Compares two fitted NMF models, i.e. objects returned by single runs of 
#' \code{\link{nmf}}.
setMethod('nmf.equal', signature(x='NMFfit', y='NMFfit'), 
		function(x, y, ...){
			nmf.equal(fit(x), fit(y), ...)
		}
)

