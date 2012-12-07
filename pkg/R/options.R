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

#' NMF Package Specific Options
#' 
#' \code{nmf.options} sets/get single or multiple options, that are specific
#' to the NMF package. 
#' It behaves in the same way as \code{\link[base]{options}}.
#' 
#' @inheritParams base::options
#' @param ... option specifications. For \code{nmf.options} this can be named arguments or 
#' a single unnamed argument that is a named list (see \code{\link{options}}.
#' 
#' For \code{nmf.resetOptions}, this must be the names of the options to reset.
#' Note that \pkg{pkgmaker} version >= 0.9.1 is required for this to work correctly,
#' when options other than the default ones have been set after the package is loaded.
#' 
#' @export
#' @rdname options
#' @examples
#' 
#' # show all NMF specific options
#' nmf.printOptions()
#' 
#' # get some options
#' nmf.getOption('verbose')
#' nmf.getOption('backend')
#' # set new values
#' nmf.options(verbose=TRUE)
#' nmf.options(backend='mc', default.algorithm='lee')
#' nmf.printOptions()
#' 
#' # reset to default
#' nmf.resetOptions()
#' nmf.printOptions()
#' 
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

#' \code{nmf.resetOptions} reset all NMF specific options to their default values.
#' 
#' @param ALL logical that indicates if options that are not part of the default set 
#' of options should be removed.
#' Note that in \pkg{pkgmaker <= 0.9} this argument is only taken into account when
#' no other argument is present. This is fixed in version 0.9.1. 
#' 
#' @export
#' @rdname options
nmf.resetOptions <- .OPTIONS$resetOptions

#' \code{nmf.printOptions} prints all NMF specific options along with their default values, 
#' in a relatively compact way.
#' @export
#' @rdname options
nmf.printOptions <- .OPTIONS$printOptions

#nmf.options.runtime <- function(){
#	nmf.options(.nmf.Options.Runtime)	
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
