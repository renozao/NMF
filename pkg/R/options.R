###% Options management
###% 
###% 
###% @author Renaud Gaujoux \email{renaud@@cbio.uct.ac.za}
###% 

#.nmf.Options.Runtime <- character()

library(pkgmaker)

# define functions nmf.options and nmf.getOptions
.OPTIONS <- setupPackageOptions(
	# default algorithm
	default.algorithm='brunet'
	# default seeding method
	, default.seed='random'
	# track error during NMF updates
	, error.track=FALSE
	# define the tracking interval
	, track.interval=30
	# define garbage collection interval
	, gc=50
	# define default parallel backend 
	, parallel.backend= option_symlink('backend') # for backward compatibility
	, backend= if( parallel::detectCores() > 1 ) 'par' else 'seq'
	# define default RNG mode
	, reproducible=TRUE
	# toogle verbosity
	, verbose=FALSE
	# toogle debug mode
	, debug=FALSE
, RESET=TRUE)

#' Package Specific Options
#' 
#' \code{nmf.options} sets/get single or multiple options, that are specific
#' to the NMF package. 
#' It behaves in the same way as \code{\link[base]{options}}.
#' 
#' @inheritParams base::options
#' 
#' @export
#' @rdname options
nmf.options <- .OPTIONS$options

#' \code{nmf.getOption} returns the value of a single option, that is specific 
#' to the NMF package.
#' It behaves in the same way as \code{\link[base]{getOption}}.
#' 
#' @inheritParams base::getOption
#' 
#' @export
#' @rdname options
nmf.getOption <- .OPTIONS$getOption

#nmf.options.runtime <- function(){
#	nmf.options(.nmf.Options.Runtime)	
#}

##TODO: store options as an environment in .GlobalEnv => allow restoration with workspace
## initialze the variables that will hold the options' data (IMPORTANT when loaded in a namespace)
#.nmf.Options <- NULL
#.nmf.Options.Builtin <- NULL
#.init.nmf.options <- function(){
#	
#	# initialize main repository
#	envName <- '.nmf.Options'
#	parentEnv <- parent.env(environment())
#	
#	msg.addon <- if( exists(envName) ) ' [reset]' else NULL
#	sreg <- utils::capture.output(print(parentEnv))
#	message("Create NMF options in ", sreg, msg.addon)		
#	.nmf.Options <<- new.env(.nmf.Options)
#	
#	# populate with default options
#	nmf.options.reset()
#	
#	# store built-in options
#	.nmf.Options.Builtin <<- names(nmf.options())
#	
#	return(invisible(TRUE))
#}
#
####% Define default options
#nmf.options.reset <- function(){
#	# default algorithm
#	nmf.options(default.algorithm='brunet')
#	# default seeding method
#	nmf.options(default.seed='random')
#	# track error during NMF updates
#	nmf.options(error.track=FALSE, runtime=TRUE)
#	# define the tracking interval
#	nmf.options(track.interval=30, runtime=TRUE)
#	# define default parallel backend 
#	nmf.options(parallel.backend='mc', runtime=TRUE)
#	# define default RNG mode
#	nmf.options(reproducible=TRUE)
#	# toogle verbosity
#	nmf.options(verbose=FALSE, runtime=TRUE)
#	# toogle debug mode
#	nmf.options(debug=FALSE)
#}

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
