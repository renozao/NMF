# S4 class for NMF strategies implemented in Octave/Matlab 
# Algorithms are run via RcppOctave
#
# Author: Renaud Gaujoux
# Created: 23 Nov 2012
###########################################################
#' @include NMFStrategy-class.R
NULL

#' S4 Interface for Octave-Matlab NMF Algorithms 
#'
#' This class implements the virtual interface \code{\linkS4class{NMFStrategy}} 
#' for NMF algorithms that are implemented in Octave/Matlab, and provided as
#' a set of .m files or as plain code.
#' 
#' The \code{run} method for this class runs the algorithms via the 
#' \code{\link{RcppOctave}} package.  
#' 
#' @slot algorithm character string that gives the name of the main Octave/Matlab 
#' function that implements the algorithm.
#' The function must take at least two arguments: the target matrix and the initial 
#' NMF model, converted into an Octave list object, with elements corresponding to 
#' slots of the corresponding S4 class.  
#' @slot mfiles character vector that contains a set of path to .m files.
#' These files are (re-)sourced every time the strategy is called, and must be 
#' present at runtime in the current directory or in a directory from Octave path. 
#' 
#' @export
setClass('NMFStrategyOctave'
	, representation(
		algorithm = 'character' # the function that implements the algorithm
		, mfiles = 'character'
		, onReturn = 'function' # method called just before returning the resulting NMF object
	)
	, prototype(
		onReturn = function(object, x){
			fit(x) <- new2(modelname(x), object)
			if( !is.null(object$runtime) )
				x@runtime <- structure(unlist(object$runtime), class='proc_time')
			x
		}
	)
	, contains = 'NMFStrategy'
)

setMethod('initialize', 'NMFStrategyOctave',
	function(.Object, ..., mfiles){
		# initialize parent
		.Object <- callNextMethod(.Object, ...)
		# process mfiles
		.Object@mfiles <- mfiles(mfiles)
		# return object
		.Object
	}
)

#' Runs the NMF algorithms implemented by the Octave/Matlab function associated with the 
#' strategy -- and stored in slot \code{'algorithm'} of \code{object}.
#' 
#' This method is usually not called directly, but only via the function \code{\link{nmf}}, which 
#' takes care of many other details such as seeding the computation, handling RNG settings, 
#' or setting up parallel computations.
#' 
#' @rdname NMFStrategy
setMethod('run', signature(object='NMFStrategyOctave', y='matrix', x='NMFfit'),
	function(object, y, x, ...){
		
		fstop <- function(...) stop("NMFStrategyOctave[", name(object), "]: ", ...)
		
		# first thing check for RcppOctave
		if( !require(RcppOctave) )
			fstop("could not load required package RcppOctave")
		
		main <- algorithm(object)
		if( !length(main) || !nchar(main) )
			fstop("main algorithm function is not defined.")
		
		# add path to all mfiles
		pdir <- packagePath('matlab', package=packageSlot(object))
		o_addpath(pdir)
		tdir <- tempdir()
		o_addpath(tdir)
		on.exit({
			rmpath <- RcppOctave::.O$rmpath
			rmpath(pdir); rmpath(tdir) 
		})
		
		# call main function
		res <- .CallOctave(main, y, list(W=basis(x), H=coef(x)), ...)
		# wrap result
		object@onReturn(res, x)
	}
)

#' M Files
#' 
#' \code{mfiles} convert mfile specifications in to real paths to .m files
#' that can be sourced with \code{\link{o_source}}.
#' 
#' @param ... specification of a .m files as character arguments.
#' The elements of the vector can be either file paths or plain Octave/Matlab code, 
#' which are then written to disk in -- temporary -- .m files. 
#' Note that the paths do not need to correspond to existing files.
#' @inheritParams base::tempfile
#' @param dir existing directory where to write the .m files generated from 
#' the plain code elements of \var{x}.
#' 
#' @export
mfiles <- function(..., pattern='mfile_', dir=tempdir()){
	
	if( missing(dir) && !is.null(ns <- getLoadingNamespace()) ){
		dir <- packagePath('matlab', package=ns)
	}
	
	# get args
	x <- unlist(list(...))
	if( !is.character(x) )
		stop("All arguments must be character strings")
	
	# detect type of input
	isfile <- !grepl("\n", x) | !grepl(" ", x)
	# add names if needed
	if( is.null(names(x)) ) names(x) <- rep('', length(x))
	
	code <- x[!isfile]
	if( length(code) ){
		x[!isfile] <- mapply(function(f, x){
			
			# create directory if it does not exist
			if( !file.exists(dir) )	dir.create(dir, recursive=TRUE)
	
			# build file path
			f <- 
			if( nchar(f) ) str_c(file.path(dir, f), ".m")
			else tempfile(pattern, tmpdir=dir, fileext=".m")
	
			# write file
			cat(x, file=f)
			
			# return filepath
			f
		}, names(code), code)
	}
	x
}

#' Returns the name of the Octave/Matlab function that implements the NMF algorithm -- as stored in 
#' slot \code{algorithm}.
setMethod('algorithm', signature(object='NMFStrategyOctave'),
		function(object){
			slot(object, 'algorithm')
		}
)
#' Sets the name of the Octave/Matlab function that implements the NMF algorithm.
#' It is stored in slot \code{algorithm}.
setReplaceMethod('algorithm', signature(object='NMFStrategyOctave', value='character'),
	function(object, value){
		slot(object, 'algorithm') <- head(value, 1L)
		object
	}
)

#' @export
#' @rdname NMFStrategyOctave-class
setMethod('show', 'NMFStrategyOctave', function(object){
		callNextMethod()
		cat(" main: ", algorithm(object), "\n", sep='')
		cat(" mfiles: ", str_out(object@mfiles, Inf), "\n", sep='')
	}
)
