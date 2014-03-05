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
#' \code{\link[RcppOctave]{RcppOctave}} package.  
#' 
#' @slot algorithm character string that gives the name of the main Octave/Matlab 
#' function that implements the algorithm.
#' The function must take at least two arguments: the target matrix and the initial 
#' NMF model, converted into an Octave list object, with elements corresponding to 
#' slots of the corresponding S4 class.  
#' @slot mcode character vector that contains a set of path to .m files.
#' These files are (re-)sourced every time the strategy is called, and must be 
#' present at runtime in the current directory or in a directory from Octave path. 
#' 
setClass('NMFStrategyOctave'
	, representation(
		algorithm = '.functionSlot' # the function that implements the algorithm
		, mcode = 'character'
		, onReturn = 'function' # method called just before returning the resulting NMF object
	)
	, prototype(
		onReturn = function(object, x){
            if( !isNMFfit(object, recursive = FALSE) ){
                
                if( !is.nmf(object) && !is.list(object) ){
                    stop("Result object should be an NMF object or a list")
                }
            	fit(x) <- new2(modelname(x), object)
    			if( !is.null(object$runtime) )
    				x@runtime <- structure(unlist(object$runtime), class='proc_time')
            }
			x
		}
	)
	, contains = 'NMFStrategy'
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
		if( !require.quiet('RcppOctave', character.only=TRUE) )
			fstop("The package RcppOctave is required to run this algorithm.\n"
				, "  Try installing it with: install.packages('RcppOctave')")
		
        # add path to all mfiles
        mdirs <- character()
        ## add package mfiles directory if possible
        if( nzchar(pkg <- packageSlot(object)) ){
            if( nzchar(pkg_mfiles <- RcppOctave::system.mfile(package=pkg)) )
                mdirs <- c(mdirs, pkg_mfiles)
        }
		## add path to specified mfiles
		mfiles <- object@mcode
		if( length(mfiles) && any(nzchar(mfiles)) ){
			mfiles <- RcppOctave::as.mfile(mfiles)
    		mdirs <- c(mdirs, dirname(mfiles))
        }
        ## add to path
        if( length(mdirs) ){
            mdirs <- unique(mdirs)
            # check which dirs were already in Octave path
    		in_path <- sapply(mdirs, RcppOctave::o_inpath)
    		sapply(mdirs[!in_path], RcppOctave::o_addpath)
            # on exit: cleanup Octave path 
    		on.exit({
    			rmpath <- RcppOctave::.O$rmpath
    			sapply(mdirs[!in_path], rmpath)
    		})
        }
        #
        
		# load algorithm
		main <- algorithm(object, load=TRUE)

		# convert matrix storage mode if necessary
		if( storage.mode(y) != 'double' ){
			storage.mode(y) <- 'double'
		}
		# call main function
		res <- main(y, x, ...)
		# wrap result
		object@onReturn(res, x)
	}
)

#' Returns the name of the Octave/Matlab function that implements the NMF algorithm -- as stored in 
#' slot \code{algorithm}.
#' 
#' @param load logical that indicates if the algorithm should be loaded as an 
#' R function. 
#' 
setMethod('algorithm', signature(object='NMFStrategyOctave'),
		function(object, load=FALSE){
			f <- slot(object, 'algorithm')
			if( !load || is.function(f) ) return(f)
			
			if( !length(f) || !nchar(f) )
				fstop("Main function is not defined for NMF algorithm '", name(object), "'.")
			
			# return wrapped into a function
			.main <- RcppOctave::o_get(f)
			function(y, x, ...){
				.main(y, r=as.numeric(nbasis(x)), W=basis(x), H=coef(x), ...)
			}
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
		f <- algorithm(object)
		cat(" main: "
			, if( is.function(f) ) str_fun(f) else str_c(f, ' <Octave function>')
			, "\n", sep='')
		cat(" mcode: ", str_out(object@mcode, Inf), "\n", sep='')
	}
)
