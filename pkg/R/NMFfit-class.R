# Implementatino of class NMFfit 
# 
# This class manages the result of a single run of a NMF algorithm.
# 
# Author: Renaud Gaujoux \email{renaud@@cbio.uct.ac.za}
###############################################################################

#' Class to store the results of \strong{Non-negative Matrix Factorization} algorithms (NMF).
#'
#' Let \eqn{V} be a \eqn{n \times m} non-negative matrix and \eqn{r} a positive integer. A NMF of \eqn{V} is commonly
#' defined as two matrices \eqn{W} and \eqn{H} such that:
#' \deqn{V \equiv W H,}      
#' where:
#' - \eqn{W} and \eqn{H} are \eqn{n \times r} and \eqn{r \times m} non-negative matrices respectivelly;
#' - \eqn{\equiv} is to be understood with respect to some error function (eg. Frobenius norm, Kullbach-Leibler divergence, ...). 
#'
#' Integer \eqn{r} is called the \emph{factorization rank}.
#' Depending on the context of application of NMF, the columns of \eqn{W} and \eqn{H} take different names:
#' - columns of \eqn{W}: metagenes, factors, source, image basis
#' - columns of \eqn{H}: metaprofiles, mixture coefficients, weights
#' 
#' Because package \code{NMF} was primilary intended to microarray data, the following terminology 
#' is used
#' \describe{
#' \item{samples}{the columns of the target matrix \eqn{V}}
#' \item{genes}{the rows of the target matrix \eqn{V}}
#' \item{metagenes}{the columns of matrix \eqn{W}} 
#' \item{metaprofiles}{the columns of matrix \eqn{H}}
#' }
#' 
#' Class \code{NMF} is a class to handle both matrices \eqn{W} and \eqn{H} within a single object,
#' together with data about the methods and parameters used to compute them.
#' 
#' \section{Validity checks}{ The validity method for class \code{NMF} checks for compatibility of slots
#' \code{W} and \code{H}, as those matrices must be compatible with respect to the matrix product. 
#' It also checks the relevance of the factorization, emmiting a warning when the factorization rank is
#' greater than the number of columns in \code{H}.
#' }
#' 
#' @seealso nmf
#' 
#' @slot W a \eqn{n \times r} \code{matrix}, the first matrix factor of the NMF. 
#' @slot H a \eqn{r \times n} \code{matrix}, the second matrix factor of the NMF.
#' @slot residuals a \code{numeric} value of the final residuals. That is a measure of how far the estimate \eqn{WH} is from the target \eqn{V}
#' @slot method a \code{character} string giving the name of the algorithm used to compute the NMF 
#' @slot seed a \code{character} string giving the name of the seeding method used to compute the NMF
#' @slot distance a \code{character} string giving the name of the loss function the algorithm is based on
#' @slot parameters a \code{list} of the algorithm's parameters used as input (in addition to \eqn{V} and \eqn{r}) 
#' @slot runtime a \code{proc_time} object giving the duration of the computation (as returned by \code{system.time}) 
#' @slot extra a \code{list} of extra values set by the algorithm for specific report/tracking 
#' @author Renaud Gaujoux \email{renaud@@cbio.uct.ac.za}
#' @export
setClass('NMFfit'
	, representation(
			fit = 'NMF', # NMF model
			residuals = 'numeric', # residuals from the target matrix
			method = 'character', # method used to compute the factorization
			seed = 'character', # seeding method used to compute the factorization
			distance = '.functionSlot.null', # method used to compute the distance between the target matrix and its NMF estimate
			parameters = 'list', # method used to compute the factorization
			runtime = 'proc_time', # running time to perform the NMF
			options = 'list', # run options
			extra = 'list' # extra list of results output by the method				
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

#' Estimate the target matrix (wrapper function) -> call method on fit
setMethod('fitted', signature(object='NMFfit'),
	function(object, ...){
		fitted(fit(object), ...)
	}
)

#' Get/Set the basis matrix -> call method on fit
setMethod('basis', signature(object='NMFfit'),
	function(object, ...){
		basis(fit(object), ...)
	}
)
setReplaceMethod('basis', signature(object='NMFfit', value='matrix'), 
	function(object, value){ 
		basis(fit(object)) <- value
		object
	} 
)

#' Get/Set the mixture coefficients matrix -> call method on fit
setMethod('coef', signature(object='NMFfit'),
	function(object, ...){
		coef(fit(object), ...)
	}
)
setReplaceMethod('coef', signature(object='NMFfit', value='matrix'), 
	function(object, value){ 
		coef(fit(object)) <- value
		object
	} 
)

setMethod('show', signature(object='NMFfit'), 
	function(object)
	{
		cat("<Object of class:", class(object), ">\n")
		cat(" # Model:\n  ")
		s <- capture.output(show(fit(object)))
		cat(s, sep="\n  ")
		cat(" # Details:\n  ")
		.local <- function(){
			if( algorithm(object) != '' ) cat("algorithm: ", algorithm(object), "\n")
			if( seeding(object) != '' ) cat("seed: ",  seeding(object), "\n")
			# distance/objective function
			svalue <- objective(object)
			svalue <- if( is.function(svalue) ) '<function>' else paste("'", svalue,"'", sep='')
			cat("distance metric: ", svalue, "\n")			
			if( length(residuals(object)) !=0  ) cat("residuals: ",  residuals(object), "\n");
			# show the parameters specific to the method		
			if( length(object@parameters) > 0 ){
				cat("parameters:")
				cat("\n")
				print(object@parameters)
			}
			# show number of iterations if present
			if( !is.null(object$iteration) ) cat("Iterations:", object$iteration, "\n")
			# show elapsed time if present
			if( length(runtime(object)) > 0 ){ cat("Timing:\n"); show(runtime(object));}
		}
		s <- capture.output(.local())
		cat(s, sep="\n  ")
	}
)


#' Returns the fit (i.e. NMF model)
if ( !isGeneric('fit') ) setGeneric('fit', function(object, ...) standardGeneric('fit'))
setMethod('fit', signature(object='NMFfit'), 
	function(object)
	{
		return(slot(object, 'fit'))
	}
)
if ( !isGeneric('fit<-') ) setGeneric('fit<-', function(object, value) standardGeneric('fit<-'))
setReplaceMethod('fit', signature(object='NMFfit', value='NMF'), 
	function(object, value){ 
		slot(object, 'fit') <- value		
		object # TODO: valid object before returning it (+param check=TRUE or FALSE)
	} 
)

#' Returns the NMF model's name
setMethod('model', signature(object='NMFfit'), 
	function(object)
	{
		return(class(fit(object)))
	}
)

if( !isGeneric('residuals') ) setGeneric('residuals', package='stats')
setMethod('residuals', 'NMFfit', 
	function(object, track=FALSE, ...){ 
		## IMPORTANT: keep this '...' and do not add a 'method' argument as this
		## one is passed by NMFSet::fit (see bug #159) and is not supposed to be 
		## used
		res <- slot(object, 'residuals')
		if( track ) res else tail(res, n=1)
	} 
)
if( !isGeneric('residuals<-') ) setGeneric('residuals<-', function(object, value) standardGeneric('residuals<-') )
setReplaceMethod('residuals', 'NMFfit',
	function(object, value){ 		
		slot(object, 'residuals') <- value 
		object
	}
)

#' Track error 
trackError <- function(object, value, iter, force=FALSE){	
	track <- run.options(object, 'error.track')
	track.interval <- run.options(object, 'track.interval')
	# add the new value to the error track
	last.iter <- names(residuals(object))
	duplicate <- if( !is.null(last.iter) ) iter == last.iter else FALSE 
	if( !duplicate && (force || (track && iter %% track.interval == 0)) ){		
		res <- c(residuals(object, track=TRUE), value)
		if( iter >= 0 ) names(res)[length(res)] <- iter
		residuals(object) <- res
	}
	object
}

if (is.null(getGeneric('algorithm'))) setGeneric('algorithm', function(object, ...) standardGeneric('algorithm') )
setMethod('algorithm', 'NMFfit', function(object){ object@method } )
if (is.null(getGeneric('algorithm<-'))) setGeneric('algorithm<-', function(object, ..., value) standardGeneric('algorithm<-') )
setReplaceMethod('algorithm', 'NMFfit',
	function(object, value){
		object@method <- value
		object
	}
)

if (is.null(getGeneric('seeding'))) setGeneric('seeding', function(object, ...) standardGeneric('seeding') )
setMethod('seeding', 'NMFfit', function(object){ object@seed } )
if (is.null(getGeneric('seeding<-'))) setGeneric('seeding<-', function(object, ..., value) standardGeneric('seeding<-') )
setReplaceMethod('seeding', 'NMFfit',
	function(object, value){
		object@seed <- value
		object
	}
)

#' Accessor methods to slot \code{objective}
if ( !isGeneric('objective') ) setGeneric('objective', function(object, ...) standardGeneric('objective'))
setMethod('objective', signature(object='NMFfit'),
	function(object, x){
		
		# when both x and y are missing then returns slot objective
		if( missing(x) ) return(slot(object, 'distance'))
		
		# return the distance computed using the strategy's objective function
		distance(x, fit(object), method=slot(object, 'distance'))
		
	}
)
if ( is.null(getGeneric('objective<-')) ) setGeneric('objective<-', function(object, ..., value) standardGeneric('objective<-'))
setReplaceMethod('objective', signature(object='NMFfit', value='character'),
	function(object, value){
		#TODO: test for the existence of objective method
		slot(object, 'distance') <- value
		validObject(object)
		object
	}
)
setReplaceMethod('objective', signature(object='NMFfit', value='function'),
	function(object, value){
		slot(object, 'distance') <- value
		validObject(object)
		object
	}
)

#' Returns slot \code{runtime}.
#' 
#' @return a numeric vector of class \code{proc_time}
#' 
if (is.null(getGeneric("runtime"))) setGeneric('runtime', function(object, ...) standardGeneric('runtime') )
setMethod('runtime', 'NMFfit', 
	function(object, ...){ 
		object@runtime; 
	}
)

#' Access methods to run options.
if (!isGeneric("run.options")) setGeneric('run.options', function(object, ...) standardGeneric('run.options') )
setMethod('run.options', 'NMFfit', 
	function(object, name){
		if( missing(name) ) object@options
		else object@options[[name]]
	}
)
if (!isGeneric("run.options<-")) setGeneric('run.options<-', function(object, ..., value) standardGeneric('run.options<-') )
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
		if( !is.element(name, names(nmf.options.runtime())) ) stop(baseError, "option '", name, "' is not defined.")
		
		object@options[[name]] <- value
		object
	}
)
if (!isGeneric("verbose")) setGeneric('verbose', function(object, ...) standardGeneric('verbose') )
setMethod('verbose', 'NMFfit', 
	function(object){
		return(run.options(object, 'verbose') || nmf.getOption('debug'))
	}
)

#' Returns an extra slot from the NMF object.
#' An returns an error if \code{name} is not a valid element of the slot \code{extra}.
if (is.null(getGeneric("extra"))) setGeneric("extra", function(object, name) standardGeneric("extra"))
setMethod('extra', 'NMFfit', 
	function(object, name){
		.Defunct("$' or '$<-")	
	}
)

#' Get/Set methods for slot 'extra'
setMethod('$', 'NMFfit', 
	function(x, name){ 
		x@extra[[name, exact=FALSE]]; 
	} 
)

setReplaceMethod('$', 'NMFfit',
	function(x, name, value) {
		x@extra[[name]] <- value
		x
	}
)

#' Plot the residuals track of a NMF result.
#'
#' When slot \code{residuals} of a NMF object contains is not a single value, this function plots
#' the curve of the objective value against the number of iterations.
#'
#' @param x a NMF object
#' @param extra graphical parameters passed to function \code{plot}
#' @return this function is used for its side effect of plotting.
#'
if ( !isGeneric('errorPlot') ) setGeneric('errorPlot', function(x, ...) standardGeneric('errorPlot') )
setMethod('errorPlot', signature(x='NMFfit'), 
	function(x, ...){
		
		# retrieve the residuals track
		track <- residuals(x, track=TRUE)
		if( length(track) <= 1 ){
			warning(class(x), ' object has no residuals track')
			return(invisible())
		}
		
		# set default graphical parameters (those can be overriden by the user)
		params <- .set.list.defaults(list(...)
				, xlab='Iterations'
				, ylab=paste('Objective value ('
							, if( is.character(x@distance) ) x@distance else algorithm(x), ')'
							, sep='' )
				, main=paste("NMF Residuals plot\nrank=", nbasis(x), sep='')
				, col='#5555ff', lwd=1.4, type='l', cex=0.5)
		
		do.call('plot', c(list(names(track), track), params))
		points(names(track), track, type='p', cex=0.6, col=params$col)
	}
)

setMethod('summary', signature(object='NMFfit'), 
	function(object, ...){
		
		res <- summary(fit(object), ...)
		
		## IMPORTANT: if adding a summary measure also add it in the sorting 
		## schema of method NMFSet::compare to allow ordering on it
		
		# nb of iterations
		res <- c(res, niter=as.integer(object$iteration) )
		# runtime
		res <- c(res, time=as.numeric(runtime(object)['user.self']))
		# retreive final residuals
		res <- c(res, residuals=as.numeric(residuals(object)))		
		
		# return result
		return(res)
	}
)

setMethod('distance', signature(target='matrix', x='NMFfit'), 
		function(target, x, method, ...){
			
			# if missing retrieve the actual distance measure from the NMF object
			if( missing(method) ) method = x@distance
			
			# compute the distance between the target and the fitted NMF model
			return(distance(target, fit(x), method=method, ...))
		}
)
