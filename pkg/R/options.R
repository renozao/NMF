#' Options management
#' 
#' 
#' @author Renaud Gaujoux \email{renaud@@cbio.uct.ac.za}
#' 

.nmf.Options.Runtime <- character()
#' Get/set package specific options
nmf.options <- function (..., runtime=FALSE){
	
	# no argument: return options as a list	
	current <- as.list(.nmf.Options, all=TRUE)	
	if (nargs() == 0){
		if( length(current) == 0 ) return(NULL)
		return( current )
	}	
		
	if (nargs() == 1 && is.list(...) ) params <- c(...)		
	else params <- list(...)
	
	# no names passed: parameters should be the names of the options to return
	if ( is.null(names(params)) ){
		if( !is.character(c(...)) )
			stop('character strings expected for option names')
		
		cparams <- c(...)
		not.options <- !(cparams %in% names(current))				
		if( !is.null(current[['debug']]) && current[['debug']] && any(not.options) ){
			message("DEBUG:: NMF options: unknown option(s) "
							, paste(paste("'", cparams[not.options], "'", sep=''), collapse=', '))
		}
		
		# retrieve options as a list (use sapply to get non set options also named)
		res <- sapply(cparams, function(n) current[[n]], simplify=FALSE)		
		return(res)
	}
	
	# override the options
	old <- sapply(names(params), 
			function(name){
				# debug message when adding a new option
				if( !is.null(current[['debug']]) && current[['debug']] 
						&& !is.element(name, names(current) ) )
					message("DEBUG:: Adding NMF option '", name,"'")
				if( runtime ){
					if( !is.element(name, .nmf.Options.Runtime) )
						.nmf.Options.Runtime <<- c(.nmf.Options.Runtime, name)					
				}
				
				# assign the new value into the options environment
				val <- params[[name]]
				if( is.null(val) ){
					if( is.element(name,.nmf.Options.Builtin) )
						stop("NMF option '", name, "' is a built-in options and cannot be set to NULL")
					if( !is.null(current[[name]]) ) remove(list=name, envir=.nmf.Options)		
				}
				else assign(name, val, envir=.nmf.Options)
				# return the option's old value
				current[[name]]
			}
			, simplify = FALSE
	)	
	
	# return old values of the modified options
	return(invisible(old))	
}

#' Get the value of a single option
nmf.getOption <- function(name){
	o <- nmf.options(name)
	stopifnot(length(o) == 1)
	return(o[[1]])
}

nmf.options.runtime <- function(){
	nmf.options(.nmf.Options.Runtime)	
}

# initialze the variables that will hold the options' data (IMPORTANT when loaded in a namespace)
.nmf.Options <- NULL
.nmf.Options.Builtin <- NULL
.init.nmf.options <- function(){
	
	# initialize main repository
	envName <- '.nmf.Options'
	parentEnv <- parent.env(environment())
	
	msg.addon <- if( exists(envName) ) ' [reset]' else NULL
	sreg <- utils::capture.output(print(parentEnv))
	message("Create NMF options in ", sreg, msg.addon)		
	.nmf.Options <<- new.env(.nmf.Options)
	
	# populate with default options
	nmf.options.reset()
	
	# store built-in options
	.nmf.Options.Builtin <<- names(nmf.options())
	
	return(invisible(TRUE))
}

#' Define default options
nmf.options.reset <- function(){
	# default algorithm
	nmf.options(default.algorithm='brunet')
	# default seeding method
	nmf.options(default.seed='random')
	# track error during NMF updates
	nmf.options(error.track=FALSE, runtime=TRUE)
	# define the tracking interval
	nmf.options(track.interval=30, runtime=TRUE)
	# define default parallel backend 
	nmf.options(parallel.backend='mc', runtime=TRUE)
	# toogle verbosity
	nmf.options(verbose=FALSE, runtime=TRUE)
	# toogle debug mode
	nmf.options(debug=FALSE)
}

# debugging utility
nmf.debug <- function(fun, ...){
	if( nmf.getOption('debug') ){
		call.stack <- sys.calls()
		n <- length(call.stack)
		if( is.null(fun) ) fun <- as.character(call.stack[[n-1]]) 
		message('DEBUG::', fun, ' -> ', ...)
	}
	return(invisible())
}
