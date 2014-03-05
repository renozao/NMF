###% Options management
###% 
###% 
###% @author Renaud Gaujoux \email{renaud@@cbio.uct.ac.za}
###% 

#.nmf.Options.Runtime <- character()

# define functions nmf.options and nmf.getOptions
#' NMF Package Specific Options
#'
#' @section Available options:
#' \describe{
#' 
#' \item{cores}{Default number of cores to use to perform parallel NMF computations.
#' Note that this option is effectively used only if the global option \code{'cores'} is 
#' not set.
#' Moreover, the number of cores can also be set at runtime, in the call to \code{\link{nmf}}, 
#' via arguments \code{.pbackend} or \code{.options} (see \code{\link{nmf}} for more details).}
#' 
#' \item{default.algorithm}{Default NMF algorithm used by the \code{nmf} function when argument 
#' \code{method} is missing. 
#' The value should the key of one of the registered NMF algorithms or a valid specification of an NMF algorithm.
#' See \code{?nmfAlgorithm}.}
#' 
#' \item{default.seed}{Default seeding method used by the \code{nmf} function when argument \code{seed} is missing.
#' The value should the key of one of the registered seeding methods or a vallid specification of a seeding method. 
#' See \code{?nmfSeed}.}
#' 
#' \item{track}{Toggle default residual tracking. 
#' When \code{TRUE}, the \code{nmf} function compute and store the residual track in the result -- if not otherwise specified in argument \code{.options}.
#' Note that tracking may significantly slow down the computations.}
#' 
#' \item{track.interval}{Number of iterations between two points in the residual track. 
#' This option is relevant only when residual tracking is enabled. 
#' See \code{?nmf}.}
#' 
#' \item{error.track}{this is a symbolic link to option \code{track} for backward compatibility.}
#' 
#' \item{pbackend}{Default loop/parallel foreach backend used by the \code{nmf} function when 
#' argument \code{.pbackend} is missing.
#' Currently the following values are supported: \code{'par'} for multicore, 
#' \code{'seq'} for sequential, \code{NA} for standard \code{sapply} (i.e. do not use a foreach loop), 
#' \code{NULL} for using the currently registered foreach backend.}
#' 
#' \item{parallel.backend}{this is a symbolic link to option \code{pbackend} for backward compatibility.}
#' 
#' \item{gc}{Interval/frequency (in number of runs) at which garbage collection is performed.}
#' 
#' \item{verbose}{Default level of verbosity.}
#' 
#' \item{debug}{Toogles debug mode.
#' In this mode the console output may be very -- very -- messy, and is aimed at debugging only.}
#' 
#' \item{maxIter}{ Default maximum number of iteration to use (default NULL).
#' This option is for internal/technical usage only, to globally speed up examples or tests
#' of NMF algorithms. To be used with care at one's own risk...
#' It is documented here so that advanced users are aware of its existence, and can avoid possible 
#' conflict with their own custom options.
#' }
#' } % end description
#' 
#' 
#' @rdname options
#' @name options-NMF
NULL
.OPTIONS <- setupPackageOptions(
	# default algorithm
	default.algorithm='brunet'
	# default seeding method
	, default.seed='random'
	# track error during NMF updates
	, error.track = option_symlink('track') # for backward compatibility
	, track=FALSE
	# define the tracking interval
	, track.interval=30
	# define garbage collection interval
	, gc=50
	# define default parallel backend 
	, parallel.backend= option_symlink('pbackend') # for backward compatibility
	, pbackend= if( parallel::detectCores() > 1 ) 'par' else 'seq'
	# toogle verbosity
	, verbose=FALSE
	# toogle debug mode
	, debug=FALSE
, RESET=TRUE)

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
#' nmf.getOption('pbackend')
#' # set new values
#' nmf.options(verbose=TRUE)
#' nmf.options(pbackend='mc', default.algorithm='lee')
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
