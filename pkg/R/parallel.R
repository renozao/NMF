# Definitions used in the parallel computations of NMF
#
# - reproducible backend
# - reproducible %dopar% operator: %dorng%
# 
# Author: Renaud Gaujoux
# Creation: 08-Feb-2011
###############################################################################

#' @include utils.R
#' @import foreach
#' @import doParallel
NULL

# returns the number of cores to use in all NMF computation when no number is
# specified by the user
getMaxCores <- function(limit=TRUE){
	#ceiling(parallel::detectCores()/2)
	nt <- n <- parallel::detectCores()
	# limit to number of cores specified in options if asked for
	if( limit ){
		if( !is.null(nc <- getOption('cores')) ) n <- nc # global option
		else if( !is.null(nc <- nmf.getOption('cores')) ) n <- nc # NMF-specific option
		else if( n > 2 ) n <- n - 1L # leave one core free if possible
	}
	# forces limiting maximum number of cores to 2 during CRAN checks
	if( n > 2 && isCHECK() ){
		message("# NOTE - CRAN check detected: limiting maximum number of cores [2/", nt, "]")
		n <- 2L
	}
	n
}

#' Utilities and Extensions for Foreach Loops
#' 
#' \code{registerDoBackend} is a unified register function for foreach backends.
#' 
#' @param object specification of a foreach backend, e.g. \sQuote{SEQ}, 
#' \sQuote{PAR} (for doParallel), \sQuote{MPI}, etc\ldots
#' @param ... extra arguments passed to the backend own registration function. 
#' 
#' @keywords internal
#' @rdname foreach
registerDoBackend <- function(object, ...){

	# restore old backend data in case of an error
	old <- getDoBackend()
	on.exit( setDoBackend(old) )
	
	# get old foreach backend object
	ob <- ForeachBackend()
	
	# register new backend: call the register method
	b <- ForeachBackend(object, ...)
	res <- register(b)
	
	# cancel backend restoration
	on.exit()
	# call old backend cleanup method
	doBackendCleanup(ob)
	
	# return old backend
	invisible(ob)
}

#' \code{getDoBackend} returns the internal data of the currently registered foreach \%dopar\% backend.
#' @rdname foreach
#' @export
getDoBackend <- function(){
	fe <- foreach:::.foreachGlobals
	if( !exists("fun", where = fe, inherits = FALSE) )
		return(NULL)
	
	c(foreach:::getDoPar() # this returns the registered %dopar% function + associated data
		# -> add info function from foreach internal environment
		, info= if( exists("info", where = fe, inherits = FALSE) ){
					get('info', fe, inherits=FALSE) 
				}else{
					function(data, item) NULL
				}
		, cleanup = if( exists("cleanup", where = fe, inherits = FALSE) ){
			get('cleanup', fe, inherits=FALSE)
		}
	)
}
#' \code{setDoBackend} is identical to \code{\link[foreach]{setDoPar}}, but 
#' returns the internal of the previously registered backend.
#' 
#' @param data internal data of a foreach \%dopar\% backend.
#' @param cleanup logical that indicates if the previous
#' backend's cleanup procedure should be run, \strong{before} 
#' setting the new backend.
#' 
#' @export
#' @rdname foreach
setDoBackend <- function(data, cleanup=FALSE){
	
	# get old backend data
	ob <- getDoBackend()
	ofb <- ForeachBackend()
	# cleanup old backend if requested
	if( cleanup ){
		doBackendCleanup(ofb)
	}
	
	if( !is.null(data) ){
		bdata <- data
		if( is.backend(data) )	data <- data[!names(data) %in% c('name', 'cleanup')]
		do.call('setDoPar', data)
		setBackendCleanup(bdata)
	}else{
		do.call('setDoPar', list(NULL))
		fe <- foreach:::.foreachGlobals
		if (exists("fun", envir = fe, inherits = FALSE))
			remove("fun", envir = fe)
		setBackendCleanup(NULL)
	}
	# return old backend
	invisible(ob)
}

# setup cleanup procedure for the current backend
setBackendCleanup <- function(object, fun, verbose=FALSE){
	
	fe <- foreach:::.foreachGlobals
	name <- getDoParName()
	if( !is.null(fun <- object$cleanup) ){
		if( verbose ) message("# Registering cleaning up function for '", name, "'... ", appendLF=FALSE)
		assign('cleanup', fun, fe)
		if( verbose ) message("OK")
	}else if (exists("cleanup", envir = fe, inherits = FALSE)){
		if( verbose ) message("# Removing cleaning up function for '", name, "'... ", appendLF=FALSE)
		remove("cleanup", envir = fe)
		if( verbose ) message("OK")
	}
	invisible(object)
}

# run cleanup procedure for a given backend object
doBackendCleanup <- function(object, ..., run=TRUE, verbose=FALSE){
	
	name <- object$name
	if( !is.null(fun <- object$cleanup) ){
		if( verbose ) message("# Cleaning up '", name, "'... ", appendLF=FALSE)
		res <- try(fun(), silent=TRUE) 
		if( verbose ) message(if( is(res, 'try-error') ) 'ERROR' else 'OK')
		if( isTRUE(res) ) object$cleanup <- NULL
		if( verbose ) message('OK', if( !is.null(res) ) str_c(' [', res,']'))
	}
	invisible(object)
}

#' \code{register} is a generic function that register objects.
#' It is used to as a unified interface to register foreach backends.
#' 
#' @param x specification of a foreach backend
#' 
#' @rdname foreach
#' @export
register <- function(x, ...){
	UseMethod('register', x)
}
#' @S3method register foreach_backend
register.foreach_backend <- function(x, ...){
	
	be <- x$name
	# For everything except doSEQ:
	# require definition package (it is safer to re-check)
	if( be != 'doSEQ' ){
		if( !require.quiet(be, character.only=TRUE) )
			stop("Package '", be, "' is required to use foreach backend '", be, "'")
	}
	
	regfun <- .foreach_regfun(x$name)
	res <- 
	if( length(formals(regfun)) > 0L ) do.call(regfun, c(x$data, ...))
	else regfun()
	# throw an error if not successful (foreach::setDoPar do not throw errors!!)
	if( is(res, 'simpleError') ) stop(res)
	# set cleanup procedure if any
	setBackendCleanup(x)
	# return result
	invisible(res)
}

#' \code{ForeachBackend} is a factory method for foreach backend objects.
#' 
#' @export
#' @inline
#' @rdname foreach
setGeneric('ForeachBackend', function(object, ...) standardGeneric('ForeachBackend'))
#' Default method defined to throw an informative error message, when no other
#' method was found.
setMethod('ForeachBackend', 'ANY', 
	function(object, ...){
		if( is.backend(object) ){
			# update arg list if necessary
			if( nargs() > 1L )	object$data <- list(...)
			object
		}else if( is(object, 'cluster') )
			selectMethod('ForeachBackend', 'cluster')(object, ...)
		else
			stop("Could not create foreach backend object with a specification of class '", class(object)[1L], "'")
	}
)

formatDoName <- function(x){
	
	# numeric values are resolved as doParallel
	if( is.numeric(x) ) x <- 'PAR'
	if( is.character(x) ){
		# use upper case if not already specified as 'do*'
		if( !grepl("^do", x) ){
			x <- toupper(x)
			# special treatment for doParallel
			if( x %in% c('PAR', 'PARALLEL') ) x <- 'Parallel'
		}
		# stick prefix 'do' (removing leading 'do' if necessary)
		str_c('do', sub('^do', '', x))
	}else 
		''
}
#' Creates a foreach backend object based on its name.
setMethod('ForeachBackend', 'character', 
	function(object, ...){
		
		object <- formatDoName(object)
		
		# build S3 class name
		s3class <- str_c(object, "_backend")
		
		# create empty S3 object
		obj <- structure(list(name=object, data=list(...))
						, class=c(s3class, 'foreach_backend'))

		# give a chance to a backend-specific ForeachBackend factory method
		# => this will generally fill the object with the elements suitable
		# to be used in a call to foreach::setDoPar: fun, data, info
		# and possibly change the name or the object class, e.g. to allow 
		# subsequent argument-dependent dispatch.
		obj <- ForeachBackend(obj, ...)
		
		# check the registration routine is available
		.foreach_regfun(obj$name)
		
		# set data slot if not already set by the backend-specific method
		if( is.null(obj$data) || (length(obj$data) == 0L && nargs()>1L) ) 
			obj$data <- list(...)
		
		# return object
		obj
	}
)
#' Creates a foreach backend object for the currently registered backend.
setMethod('ForeachBackend', 'missing', 
	function(object, ...){
		be <- getDoParName()
		data <- getDoBackend()
		bdata <- data$data
		res <- if( !is.null(bdata) ) do.call(ForeachBackend, c(list(be, bdata), ...))
		else ForeachBackend(be, ...)
		if( !is.null(data$cleanup) ) res$cleanup <- data$cleanup
		res
	}
)
#' Dummy method that returns \code{NULL}, defined for correct dispatch.
setMethod('ForeachBackend', 'NULL', function(object, ...){ NULL })

setOldClass('cluster')
#' Creates a doParallel foreach backend that uses the cluster described in 
#' \code{object}.
setMethod('ForeachBackend', 'cluster', 
	function(object, ...){
		ForeachBackend('doParallel', cl=object)
	}
)
#' Creates a doParallel foreach backend with \code{object} processes.
setMethod('ForeachBackend', 'numeric', 
	function(object, ...){
		# check numeric specification
		if( length(object) == 0L )
			stop("invalid number of cores specified as a backend [empty]")
		object <- object[1]
		if( object <= 0 )
			stop("invalid negative number of cores [", object, "] specified for backend 'doParallel'")
		
		ForeachBackend('doParallel', cl=object, ...)
	}
)
###############
# doParallel
###############
setOldClass('doParallel_backend')
#' doParallel-specific backend factory
#' 
#' @param cl cluster specification: a cluster object or a numeric that indicates the 
#' number of nodes to use. 
#' @param type type of cluster, See \code{\link[parallel]{makeCluster}}.
setMethod('ForeachBackend', 'doParallel_backend',
	function(object, cl, type=NULL){
		
		# use all available cores if not otherwise specified
		register_cleanup <-
		if( missing(cl) ){
			!length(object$data) || isNumber(object$data) || isNumber(object$data[[1L]])
		}else{
			isNumber(cl)
		}

		# On Windows doParallel::registerDoParallel(numeric) will create a 
		# SOCKcluster with `object` cores.
		# On non-Windows machines registerDoParallel(numeric) will use 
		# parallel::mclapply with `object` cores.
		# => Windows needs a cleanup function that will stop the cluster 
		# when another backend is registered.
		#
		# Fortunately doParallel::registerDoParallel assign the cluster object 
		# to the global variable `.revoDoParCluster`
#		if ( register_cleanup ) {
#			doBackendCleanup(object, function(x, force=FALSE){
#				# do not cleanup if current backend is the same (TODO: improve this)
#				if( !force && getDoParName() == 'doParallel') return()
#				# on Windows: stop cluster stored in global variable `.revoDoParCluster`
#				if( .Platform$OS.type == "windows" ){
#					.cl_cleanup(".revoDoParCluster")
#				}
#				# on all Platforms: try to cleanup PSOCK
#				.cl_cleanup('.doParPSOCKCluster')
#			})
#		}
		
		# required registration data
		# TODO: a function doParallel:::doParallel should exist and do the same 
		# thing as parallel::registerDoParallel without registering the backend
		#object$fun <- doParallel:::doParallel
		object$info <- doParallel:::info
		
		# set type of cluster if explicitly provided
		if( !is.null(type) ) object$data$type <- type
		
		# return object
		object
	}
)

######################################################
# doPSOCK
# Default snow-like cluster from parallel on Windows 
# but works on Unix as well
######################################################

setOldClass('doPSOCK_backend')
#' doSNOW-specific backend factory
setMethod('ForeachBackend', 'doPSOCK_backend',
		function(object, cl){
			
			# use all available cores if not otherwise specified
			if( missing(cl) ) cl <- getMaxCores()
			
			# return equivalent doParallel object
			ForeachBackend('doParallel', cl, type='PSOCK')
		}
)

.cl_cleanup <- function(gvar, envir=.GlobalEnv){
	if( !exists(gvar, envir = envir) ) return()
	cl <- get(gvar, envir = envir)
	try( parallel::stopCluster(cl), silent=TRUE)
	rm(list=gvar, envir = envir)
	TRUE
} 

cleanupCluster <- function(x, cl, stopFun=NULL){
	
	function(){
		
		if( is(x, 'doParallel_backend') ){
			
			# on Windows: stop cluster stored in global variable `.revoDoParCluster`
			if( .Platform$OS.type == "windows" ){
				.cl_cleanup(".revoDoParCluster")
			}
		}
		
		if( is.null(stopFun) ) stopFun <- parallel::stopCluster 
		# stop cluster
		stopFun(cl)
		TRUE
	}
}

#' @S3method register doParallel_backend
register.doParallel_backend <- function(x, ...){
	
	# start cluster if numeric specification and type is defined
	cl <- x$data[[1]]
  if( is.numeric(cl) && (.Platform$OS.type == 'windows' || !is.null(x$data$type)) ){
		names(x$data)[1L] <- 'spec'
		# start cluster
		clObj <- do.call(parallel::makeCluster, x$data)
		x$data <- list(clObj)
		# setup cleanup procedure
		x$cleanup <- cleanupCluster(x, clObj)
	}
	# register
	register.foreach_backend(x, ...)
}

###############
# doMPI 
###############

isMPIBackend <- function(x, ...){
	b <- if( missing(x) ) ForeachBackend(...) else ForeachBackend(object=x, ...)
	if( is.null(b) ) FALSE
	else if( identical(b$name, 'doMPI') ) TRUE 
	else if( length(b$data) ){
		is(b$data[[1]], 'MPIcluster') || is(b$data[[1]], 'mpicluster')
	}else FALSE
}

#' @S3method register doMPI_backend
register.doMPI_backend <- function(x, ...){
	
	if( length(x$data) && isNumber(cl <- x$data[[1]]) ){
		clObj <- doMPI::startMPIcluster(cl)
		x$data[[1]] <- clObj
		# setup cleanup procedure
		x$cleanup <- cleanupCluster(x, clObj, doMPI::closeCluster)
	}
	# register
	register.foreach_backend(x, ...)
}

setOldClass('mpicluster')
#' Creates a doMPI foreach backend that uses the MPI cluster described in 
#' \code{object}.
setMethod('ForeachBackend', 'mpicluster', 
	function(object, ...){
		ForeachBackend('doMPI', cl=object)
	}
)

setOldClass('doMPI_backend')
#' doMPI-specific backend factory
setMethod('ForeachBackend', 'doMPI_backend',
	function(object, cl){
		
		# use all available cores if not otherwise specified
		if( missing(cl) ) cl <- getMaxCores()
				
		# required registration data
		object$fun <- doMPI:::doMPI
		object$info <- doMPI:::info
		
		# return object
		object
	}
)

#as.foreach_backend <- function(x, ...){
#	
#	args <- list(...)
#	if( is.backend(x) ){
#		# update arg list if necessary
#		if( length(args) > 0L )	x$args <- args
#		return(x)
#	}
#	
#	be <-
#	if( is.null(x) ){
#		getDoParName()
#	} else if( is(x, 'cluster') || is.numeric(x) ){
#		# check numeric specification
#		if( is.numeric(x) ){
#			if( length(x) == 0L )
#				stop("invalid number of cores specified as a backend [empty]")
#			x <- x[1]
#			if( x <= 0 )
#				stop("invalid negative number of cores [", x, "] specified for backend 'doParallel'")
#		}
#		
#		args$spec <- x
#		'Parallel'
#	} else if( is(x, 'mpicluster') ){
#		args$spec <- x
#		'MPI'
#	} else if( is.character(x) ){
#		toupper(x)
#	} else 
#		stop("invalid backend specification: must be NULL, a valid backend name, a numeric value or a cluster object [", class(x)[1L], "]")
#
#	if( be %in% c('PAR', 'PARALLEL') ) be <- 'Parallel'
#	# remove leading 'do'
#	be <- str_c('do', sub('^do', '', be))
#	# build S3 class name
#	s3class <- str_c(be, "_backend")
#	
#	# check the registration routine is available
#	regfun <- .foreach_regfun(be)
#	
#	structure(list(name=be, args=args), class=c(s3class, 'foreach_backend'))
#}

is.backend <- function(x) is(x, 'foreach_backend')

#' @S3method print foreach_backend
print.foreach_backend <- function(x, ...){
	cat("<foreach backend:", x$name, ">\n", sep='')
	if( length(x$data) ){
		cat("Specifications:\n")
		str(x$data)
	}
}

.foreach_regfun <- function(name){
	
	# early exit for doSEQ
	if( name == 'doSEQ' ) return( registerDoSEQ )
	
	# build name of registration function
	s <- str_c(toupper(substring(name, 1,1)), substring(name, 2))
	funname <- str_c('register', s)
	s3class <- str_c(name, "_backend")
	
	# require definition package
	if( !require.quiet(name, character.only=TRUE) )
		stop("could not find package for foreach backend '", name, "'")
	# check for registering function or generic
	if( is.null(regfun <- getFunction(funname, mustFind=FALSE, where=asNamespace(name))) ){
		if( is.null(regfun <- getS3method('register', s3class, optional=TRUE)) )
			stop("could not find registration routine for foreach backend '", name, "'")
		#			stop("backend '", name,"' is not supported: function "
		#							,"`", regfun, "` and S3 method `register.", s3class, "` not found.")
	}
	regfun
}


#' \code{getDoParHosts} is a generic function that returns the hostname of the worker nodes used by a backend.
#' 
#' @export
#' @rdname foreach
#' @inline
setGeneric('getDoParHosts', function(object, ...) standardGeneric('getDoParHosts'))
setOldClass('foreach_backend')
#' Default method that tries to heuristaically infer the number of hosts and in last 
#' resort temporarly register the backend and performs a foreach loop, to retrieve the 
#' nodename from each worker.
setMethod('getDoParHosts', 'ANY',
	function(object, ...){
		
		be <- if( missing(object) ) ForeachBackend(...) else ForeachBackend(object, ...)
		if( existsMethod('getDoParHosts', class(be)[1L]) ) return( callGeneric(object) )
		
		# default behaviour
		nodename <- setNames(Sys.info()['nodename'], NULL)
			
		if( is.null(be) || is.null(be$data) ) return( NULL )
		# doSEQ
		if( be$name == 'doSEQ' ) 
			return( nodename )
		if( isNumber(be$data) ) 
			return( rep(nodename, be$data) )
		if( length(be$data) && isNumber(be$data[[1]]) ) 
			return( rep(nodename, be$data[[1]]) )
		if( length(be$data) && be$name == 'doParallel' ) 
			return( sapply(be$data[[1L]], '[[', 'host') )
		
		if( !missing(object) ){ # backend passed: register temporarly
			ob <- getDoBackend()
			on.exit( setDoBackend(ob) )
			registerDoBackend(be)
		}
		setNames(unlist(times(getDoParWorkers()) %dopar% { Sys.info()['nodename'] }), NULL)
	}
)

#' \code{getDoParNHosts} returns the number of hosts used by a backend.
#' 
#' @export
#' @rdname foreach
getDoParNHosts <- function(object){
	if( missing(object) ) foreach::getDoParWorkers()
	else{
		length(getDoParHosts(object))
	}
}

# add new option: limit.cores indicates if the number of cores used in parallel 
# computation can exceed the detected number of CPUs on the host. 
#.OPTIONS$newOptions(limit.cores=TRUE)

#' Computational Setup Functions
#' 
#' @description
#' Functions used internally to setup the computational environment.
#' 
#' \code{setupBackend} sets up a foreach backend given some specifications.
#' 
#' @param spec target parallel specification: either \code{TRUE} or \code{FALSE},
#' or a single numeric value that specifies the number of cores to setup. 
#' @param backend value from argument \code{.pbackend} of \code{nmf}.
#' @param optional a logical that indicates if the specification must be fully 
#' satisfied, throwing an error if it is not, or if one can switch back to 
#' sequential, only outputting a verbose message.
#' @param verbose logical or integer level of verbosity for message outputs.
#' 
#' @return Returns \code{FALSE} if no foreach backend is to be used, \code{NA} if the currently 
#' registered backend is to be used, or, if this function call registered a new backend, 
#' the previously registered backend as a \code{foreach} object, so that it can be restored 
#' after the computation is over.
#' @keywords internals
#' @rdname setup
setupBackend <- function(spec, backend, optional=FALSE, verbose=FALSE){

	pbackend <- backend
	str_backend <- quick_str(pbackend)
	# early exit: FALSE specification or NA backend means not using foreach at all
	if( isFALSE(spec) || is_NA(pbackend) ) return(FALSE)
	# use doParallel with number of cores if specified in backend
	if( is.numeric(pbackend) ){
		spec <- pbackend
		pbackend <- 'PAR'
	}
	# identify doSEQ calls
	doSEQ <- formatDoName(pbackend) == 'doSEQ'
		
	# custom error function
	pcomp <- is.numeric(spec) && !identical(spec[1], 1)
	errorFun <- function(value=FALSE, stop=FALSE, level=1){
		function(e, ...){
			if( !is(e, 'error') ) e <- list(message=str_c(e, ...))
			
			pref <- if( pcomp ) "Parallel" else "Foreach"
			if( !optional || stop ){
				if( verbose >= level ) message('ERROR')
				stop(pref, " computation aborted: ", e$message, call.=FALSE)
			}else if( verbose >= level ){
				message('NOTE')
				message("# NOTE: ", pref, " computation disabled: ", e$message)
			}
			value
		}
	}
	
	# check current backend if backend is NULL
	if( is.null(pbackend) ){
		if( verbose > 1 ){
			message("# Using current backend ... ", appendLF=FALSE)
		}
		ok <- tryCatch({
			if( is.null(parname <- getDoParName()) )
				stop("argument '.pbackend' is NULL but there is no registered backend")
			if( verbose > 1 ) message('OK [', parname, ']')
			TRUE
		}, error = errorFun())
		if( !ok ) return(FALSE)
		# exit now since there is nothing to setup, nothing should change
		# return NULL so that the backend is not restored on.exit of the parent call.
		return(NA)
	}
	##
	
	# test if requested number of cores is actually available
	NCORES <- getMaxCores(limit=FALSE)
	if( verbose > 2 ) message("# Check available cores ... [", NCORES, ']')
	if( verbose > 2 ) message("# Check requested cores ... ", appendLF=FALSE)
	ncores <- if( doSEQ ) 1L
		else{
			ncores <- tryCatch({
					if( is.numeric(spec) ){
						if( length(spec) == 0L )
							stop("no number of cores specified for backend '", str_backend, "'")
						spec <- spec[1]
						if( spec <= 0L )
							stop("invalid negative number of cores [", spec, "] specified for backend '", str_backend, "'")
						spec
					}else # by default use the 'cores' option or half the number of cores
						getMaxCores() #getOption('cores', ceiling(NCORES/2))
				}, error = errorFun(stop=TRUE))
			if( isFALSE(ncores) ) return(FALSE)
			ncores
		}
	if( verbose > 2 ) message('[', ncores, ']')
	
	# create backend object
	if( verbose > 2 ) message("# Loading backend for specification `", str_backend, "` ... ", appendLF=FALSE)
	newBackend <- tryCatch({
			# NB: limit to the number of cores available on the host 
			if( !doSEQ ) ForeachBackend(pbackend, min(ncores, NCORES))
			else ForeachBackend(pbackend)
		}, error = errorFun(level=3))
	if( isFALSE(newBackend) ) return(FALSE)
	if( verbose > 2 ) message('OK')
	
	if( verbose > 1 ) message("# Check host compatibility ... ", appendLF=FALSE)
	ok <- tryCatch({
		# check if we're not running on MAC from GUI
		if( is.Mac(check.gui=TRUE) && (newBackend$name == 'doMC' || (newBackend$name == 'doParallel' && is.numeric(newBackend$data[[1]]))) ){
			# error only if the parallel computation was explicitly asked by the user
			stop("multicore parallel computations are not safe from R.app on Mac OS X."
					, "\n  -> Use a terminal session, starting R from the command line.")
		}
		TRUE
		}, error = errorFun())
	if( !ok ) return(FALSE)
	if( verbose > 1 ) message('OK')
	
	if( verbose > 1 ) message("# Registering backend `", newBackend$name, "` ... ", appendLF=FALSE)
	# try registering the backend
	oldBackend <- getDoBackend()
	# setup retoration of backend in case of an error
	# NB: the new backend cleanup will happens only 
	# if regsitration succeeds, since the cleanup routine is 
	# setup after the registration by the suitable register S3 method. 
	on.exit( setDoBackend(oldBackend, cleanup=TRUE) )
	
	ov <- lverbose(verbose)
	ok <- tryCatch({
			registerDoBackend(newBackend)
			TRUE
		}
		, error ={
			lverbose(ov)
			errorFun()
		})
	lverbose(ov)
	if( !ok ) return(FALSE)
	if( verbose > 1 ) message('OK')
	
	# check allocated cores if not doSEQ backend
	if( newBackend$name != 'doSEQ' ){
		# test allocated number of cores
		if( verbose > 2 ) message("# Check allocated cores ... ", appendLF=FALSE)
		wcores <- getDoParWorkers()
		if( ncores > 0L && wcores < ncores ){
			if( !optional ){
				errorFun(level=3)("only ", wcores, " core(s) available [requested ", ncores ," core(s)]")
			}else if( verbose > 2 ){
				message('NOTE [', wcores, '/', ncores, ']')
				message("# NOTE: using only ", wcores,
						" core(s) [requested ", ncores ," core(s)]")
			}
		}
		else if( verbose > 2 ){
			message('OK [', wcores, '/', ncores
					, if(ncores != NCORES ) str_c(' out of ', NCORES)
					, ']')
		}
	}
	
	# cancel backend restoration
	on.exit()
	# return old backend
	oldBackend
}


# add extra package bigmemory and synchronicity on Unix platforms
if( .Platform$OS.type != 'windows' ){
	setPackageExtra('install.packages', 'bigmemory', pkgs='bigmemory')
	setPackageExtra('install.packages', 'synchronicity', pkgs='synchronicity')
}
# add new option: shared.memory that indicates if one should try using shared memory
# to speed-up parallel computations.
.OPTIONS$newOptions(shared.memory = (.Platform$OS.type != 'windows' && !is.Mac()))


#' \code{setupSharedMemory} checks if one can use the packages \emph{bigmemory} and \emph{sychronicity}
#' to speed-up parallel computations when not keeping all the fits.
#' When both these packages are available, only one result per host is written on disk,
#' with its achieved deviance stored in shared memory, that is accessible to all cores on 
#' a same host.
#' It returns \code{TRUE} if both packages are available and NMF option \code{'shared'} is 
#' toggled on. 
#' 
#' @rdname setup 
setupSharedMemory <- function(verbose){
	
	if( verbose > 1 ) message("# Check shared memory capability ... ", appendLF=FALSE)
	# early exit if option shared is off
	if( !nmf.getOption('shared.memory') ){
		if( verbose > 1 ) message('SKIP [disabled]')
		return(FALSE)
	}
	# early exit if foreach backend is doMPI: it is not working, not sure why
	if( isMPIBackend() ){
		if( verbose > 1 ) message('SKIP [MPI cluster]')
		return(FALSE)
	}
	# not on Windows
	if( .Platform$OS.type == 'windows' ){
		if( verbose > 1 ) message('SKIP [Windows OS]')
		return(FALSE)
	}
	
	if( !require.quiet('bigmemory', character.only=TRUE) ){
		if( verbose > 1 ){
			message('NO', if( verbose > 2 ) ' [Package `bigmemory` required]')
		}
		return(FALSE)
	}
	if( !require.quiet('synchronicity', character.only=TRUE) ){
		if( verbose > 1 ){
			message('NO', if( verbose > 2 ) ' [Package `synchronicity` required]')
		}
		return(FALSE)
	}
	if( verbose > 1 ) message('YES', if( verbose > 2 ) ' [synchronicity]')
	TRUE
}

is.doSEQ <- function(){
	dn <- getDoParName()
	is.null(dn) || dn == 'doSEQ'
}

#' \code{setupTempDirectory} creates a temporary directory to store the best fits computed on each host.
#' It ensures each worker process has access to it.
#' 
#' @rdname setup
setupTempDirectory <- function(verbose){
	
	# - Create a temporary directory to store the best fits computed on each host
	NMF_TMPDIR <- tempfile('NMF_', getwd())
	if( verbose > 2 ) message("# Setup temporary directory: '", NMF_TMPDIR, "' ... ", appendLF=FALSE)
	dir.create(NMF_TMPDIR)
	if( !is.dir(NMF_TMPDIR) ){
		if( verbose > 2 ) message('ERROR')
		nmf_stop('nmf', "could not create temporary result directory '", NMF_TMPDIR, "'")
	}
	
	on.exit( unlink(NMF_TMPDIR, recursive=TRUE) )
	# ensure that all workers can see the temporary directory
	wd <- times(getDoParWorkers()) %dopar% {
		if( !file_test('-d', NMF_TMPDIR) )
			dir.create(NMF_TMPDIR, recursive=TRUE)
		file_test('-d', NMF_TMPDIR)
	}
	# check it worked
	if( any(!wd) ){
		if( verbose > 2 ) message('ERROR')
		nmf_stop('nmf', "could not create/see temporary result directory '", NMF_TMPDIR, "' on worker nodes ", str_out(which(!wd), Inf))	
	}
	if( verbose > 2 ) message('OK')
	on.exit()
	NMF_TMPDIR
}

#' Utilities for Parallel Computations
#'
#' 
#' @rdname parallel
#' @name parallel-NMF 
NULL

#' \code{ts_eval} generates a thread safe version of \code{\link{eval}}.
#' It uses boost mutexes provided by the \code{\link[synchronicity]{synchronicity}}
#' package.
#' The generated function has arguments \code{expr} and \code{envir}, which are passed
#' to \code{\link{eval}}.
#' 
#' @param mutex a mutex or a mutex descriptor.
#' If missing, a new mutex is created via the function \code{\link[synchronicity]{boost.mutex}}.
#' @param verbose a logical that indicates if messages should be printed when 
#' locking and unlocking the mutex.
#' 
#' @rdname parallel
#' @export
ts_eval <- function(mutex = synchronicity::boost.mutex(), verbose=FALSE){
	
	
	library(bigmemory)
	library(synchronicity)
	# describe mutex if necessary
	.MUTEX_DESC <- 
			if( is(mutex, 'boost.mutex') ) synchronicity::describe(mutex)
			else mutex
	
	loadpkg <- TRUE
	function(expr, envir=parent.frame()){
		
		# load packages once
		if( loadpkg ){
			library(bigmemory)
			library(synchronicity)
			loadpkg <<- FALSE
		}
		MUTEX <- synchronicity::attach.mutex(.MUTEX_DESC)
		synchronicity::lock(MUTEX)
		if( verbose )
			message('#', Sys.getpid(), " - START mutex: ", .MUTEX_DESC@description$shared.name)
		ERROR <- "### <Error in mutex expression> ###\n"
		on.exit({
			if( verbose ){
				message(ERROR, '#', Sys.getpid(), " - END mutex: ", .MUTEX_DESC@description$shared.name)
			}
			synchronicity::unlock(MUTEX)
		})
		
		eval(expr, envir=envir)
		
		ERROR <- NULL
	}	
}

#' \code{ts_tempfile} generates a \emph{unique} temporary filename 
#' that includes the name of the host machine and/or the caller's process id, 
#' so that it is thread safe.
#' 
#' @inheritParams base::tempfile
#' @param ... extra arguments passed to \code{\link[base]{tempfile}}.
#' @param host logical that indicates if the host machine name should 
#' be appear in the filename.
#' @param pid logical that indicates if the current process id 
#' be appear in the filename.
#' 
#' @rdname parallel
#' @export
ts_tempfile <- function(pattern = "file", ..., host=TRUE, pid=TRUE){
	if( host ) pattern <- c(pattern, Sys.info()['nodename'])
	if( pid ) pattern <- c(pattern, Sys.getpid())
	tempfile(paste(pattern, collapse='_'), ...)
}

#' \code{hostfile} generates a temporary filename composed with  
#' the name of the host machine and/or the current process id.
#' 
#' @inheritParams base::tempfile
#' @inheritParams ts_tempfile
#' 
#' @rdname parallel
#' @export
hostfile <- function(pattern = "file", tmpdir=tempdir(), fileext='', host=TRUE, pid=TRUE){
	if( host ) pattern <- c(pattern, Sys.info()['nodename'])
	if( pid ) pattern <- c(pattern, Sys.getpid())
	file.path(tmpdir, str_c(paste(pattern, collapse='.'), fileext))
}

#' \code{gVariable} generates a function that access a global static variable, 
#' possibly in shared memory (only for numeric matrix-coercible data in this case).
#' It is used primarily in parallel computations, to preserve data accross 
#' computations that are performed by the same process.
#' 
#' @param init initial value
#' @param shared a logical that indicates if the variable should be stored in shared 
#' memory or in a local environment.
#'  
#' @rdname parallel
#' @export
gVariable <- function(init, shared=FALSE){
	
	if( shared ){ # use bigmemory shared matrices
		if( !is.matrix(init) )
			init <- as.matrix(init)
		library(bigmemory)
		DATA <- bigmemory::as.big.matrix(init, type='double', shared=TRUE)
		DATA_DESC <- bigmemory::describe(DATA)
	}else{ # use variables assigned to .GlobalEnv
		DATA_DESC <- basename(tempfile('.gVariable_'))
	}
	
	.VALUE <- NULL
	.loadpkg <- TRUE
	function(value){
		
		# load packages once
		if( shared && .loadpkg ){
			library(bigmemory)
			.loadpkg <<- FALSE	
		}
		
		# if shared: attach bigmemory matrix from its descriptor object
		if( shared ){
			DATA <- bigmemory::attach.big.matrix(DATA_DESC)
		}
		
		if( missing(value) ){# READ ACCESS
			if( !shared ){
				# initialise on first call if necessary
				if( is.null(.VALUE) ) .VALUE <<- init
				# return variable
				.VALUE
			}else 
				DATA[]
			
		}else{# WRITE ACCESS
			if( !shared ) .VALUE <<- value 
			else DATA[] <- value
			
		}
	}
}

#' \code{setupLibPaths} add the path to the NMF package to each workers' libPaths. 
#' 
#' @param pkg package name whose path should be exported the workers.
#' 
#' @rdname setup
setupLibPaths <- function(pkg='NMF', verbose=FALSE){
	
	# do nothing in sequential mode
	if( is.doSEQ() ) return( character() )
	
	if( verbose ){
		message("# Setting up libpath on workers for package(s) "
			, str_out(pkg, Inf), ' ... ', appendLF=FALSE)
	}
	p <- path.package(pkg)
	if( is.null(p) ) return()
	
	if( !isDevNamespace(pkg) ){ # not a dev package
		plibs <- dirname(p)
		libs <- times(getDoParWorkers()) %dopar% {
			.libPaths(c(.libPaths(), plibs))
		}
		libs <- unique(unlist(libs))
		if( verbose ){
			message("OK\n# libPaths:\n", paste('  ', libs, collapse="\n"))
		}
		libs
		pkg
	}else if( getDoParName() != 'doParallel' || !isNumber(getDoBackend()$data) ){ 
		# devmode: load the package + depends
		if( verbose ){ message("[devtools::load_all] ", appendLF=FALSE) }
		times(getDoParWorkers()) %dopar% {
			capture.output({
				suppressMessages({
					library(devtools)
					library(bigmemory)
					library(rngtools)
					load_all(p)
				})
			})
		}
		if( verbose ){ message("OK") }
		c('bigmemory', 'rngtools')
	}
	else if( verbose ){
		message("OK")
	}
}

#StaticWorkspace <- function(..., .SHARED=FALSE){
#		
#	# create environment
#	e <- new.env(parent=.GlobalEnv)
#	# fill with initial data
#	vars <- list(...)
#	if( .SHARED ){
#		lapply(names(vars), function(x){
#			bm <- bigmemory::as.big.matrix(vars[[x]], type='double', shared=TRUE)
#			e[[x]] <- bigmemory::describe(bm)
#		})
#	}else
#		list2env(vars, envir=e)
#	
#	structure(e, shared=.SHARED, class=c("static_wsp", 'environment'))
#}
#
#`[[.static_wsp` <- function(x, ..., exact = TRUE){
#	if( attr(x, 'shared') ){
#		var <- bigmemory::attach.big.matrix(NextMethod())
#		var[]
#	}else
#		NextMethod()
#}
#
#`[[.static_wsp<-` <- function(x, i, value){
#	
#	if( attr(x, 'shared') ){
#		var <- bigmemory::attach.big.matrix(x[[i]])
#		var[] <- value
#	}else
#		x[[i]] <- value
#	x
#}


isRNGseed <- function(x){
	is.numeric(x) || 
			( is.list(x) 
				&& is.null(names(x)) 
				&& all(sapply(x, is.numeric)) )
}

#' \code{setupRNG} sets the RNG for use by the function nmf.
#' It returns the old RNG as an rstream object or the result of set.seed 
#'  if the RNG is not changed due to one of the following reason:
#'  - the settings are not compatible with rstream  
#' 
#' @param seed initial RNG seed specification
#' @param n number of RNG seeds to generate
#' 
#' @rdname setup
setupRNG <- function(seed, n, verbose=FALSE){
	
	if( verbose == 2 ){
		message("# Setting up RNG ... ", appendLF=FALSE)
		on.exit( if( verbose == 2 ) message("OK") )
	}else if( verbose > 2 ) message("# Setting up RNG ... ")
	
	if( verbose > 3 ){
		message("# ** Original RNG settings:")
		showRNG()
	}
	
	# for multiple runs one always uses RNGstreams
	if( n > 1 ){
		
		# seeding with numeric values only
		if( is.list(seed) && isRNGseed(seed) ){
			if( length(seed) != n )
				stop("Invalid list of RNG seeds: must be of length ", n)
			
			if( verbose > 2 ) message("# Using supplied list of RNG seeds")
			return(seed)
			
		}else if( is.numeric(seed) ){
			
			if( verbose > 2 ){
				message("# Generate RNGStream sequence using seed ("
						, RNGstr(seed), ") ... "
						, appendLF=FALSE)
			}
			res <- RNGseq(n, seed)
			if( verbose > 2 ) message("OK")
			return(res)
			
		}else{ # create a sequence of RNGstream using a random seed
			if( verbose > 2 ){
				message("# Generate RNGStream sequence using a random seed ... "
						, appendLF=FALSE)
			}
			res <- RNGseq(n, NULL)
			if( verbose > 2 ) message("OK")
			return(res)
		}
	}else if( is.numeric(seed) ){ 
		# for single runs: 1-length seeds are used to set the current RNG
		# 6-length seeds are used to set RNGstream
		
		if( !is.vector(seed) ){
			message('ERROR')
			stop("NMF::nmf - Invalid numeric seed: expects a numeric vector.")
		}
		
		# convert to an integer vector
		seed <- as.integer(seed)
		# immediately setup the RNG in the standard way		
		if( length(seed) == 1L ){
			if( verbose > 2 ){
				message("# RNG setup: standard [seeding current RNG]")
				message("# Seeding current RNG with seed (", seed, ") ... "
						, appendLF=FALSE)
			}
			set.seed(seed)
			if( verbose > 2 ) message("OK")				
			return( getRNG() )
		}else if( length(seed) == 6L ){
			if( verbose > 2 ){
				message("# RNG setup: reproducible [using RNGstream]")
				message("# Generate RNGStream sequence using seed ("
						, RNGstr(seed), ") ... "
						, appendLF=FALSE)
			}
			res <- RNGseq(1, seed)
			setRNG(res)
			if( verbose > 2 ) message("OK")
			return( res )
		}else{
			if( verbose > 2 ){
				message("# RNG setup: directly setting RNG")
				message("# Setting RNG with .Random.seed= ("
						, RNGstr(seed), ") ... "
						, appendLF=FALSE)
			}
			setRNG(seed, verbose > 2)
			if( verbose > 2 ) message("OK")
			return( getRNG() )
		}
		stop("NMF::nmf - Invalid numeric seed: unexpected error.")
	}else{
		if( verbose > 2 ) message("# RNG setup: standard [using current RNG]")
		NULL
	}
} 

##################################################################
## END
##################################################################
