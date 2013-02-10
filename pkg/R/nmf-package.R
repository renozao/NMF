#' @import rngtools
#' @import digest
#' @import stringr
#' @import stats
NULL
library(digest)

#' Defunct Functions and Classes in the NMF Package
#' 
#' @name NMF-defunct
#' @rdname NMF-defunct
NULL

#' Deprecated Functions in the Package NMF
#' 
#' @param object an R object
#' @param ... extra arguments 
#' 
#' @name NMF-deprecated
#' @rdname NMF-deprecated
NULL



#' Algorithms and framework for Nonnegative Matrix Factorization (NMF).
#' 
#' This package provides a framework to perform Non-negative Matrix Factorization (NMF).
#' It implements a set of already published algorithms and seeding methods, and provides a framework 
#' to test, develop and plug new/custom algorithms. 
#' Most of the built-in algorithms have been optimized in C++, and the main interface function provides 
#' an easy way of performing parallel computations on multicore machines.
#' 
#' \code{\link{nmf}} Run a given NMF algorithm
#' 
#' @author Renaud Gaujoux \email{renaud@@cbio.uct.ac.za}
#' @name NMF-package
#' @aliases NMF
#' @docType package
#' @useDynLib NMF
#' 
#' @bibliography ~/Documents/articles/library.bib
#' @references
#' \url{http://www.r-project.org/}
#' @keywords package
#' @seealso \code{\link{nmf}}
#' @examples
#' # generate a synthetic dataset with known classes
#' n <- 50; counts <- c(5, 5, 8);
#' V <- syntheticNMF(n, counts)
#' 
#' # perform a 3-rank NMF using the default algorithm
#' res <- nmf(V, 3)
#' 
#' basismap(res)
#' coefmap(res)
#' 
NA

devnmf <- function(){
	.LOCAL_PKG_NAME <- 'NMF'
	load_all(.LOCAL_PKG_NAME)
	compile_src(.LOCAL_PKG_NAME)
}

.onLoad <- function(libname, pkgname=NULL) {
	
	pkgEnv <- pkgmaker::packageEnv()
	
	pkgmaker::runPostponedAction(group='extra')
	
	.init.sequence <- function(){
	
		## 0. INITIALIZE PACKAGE SPECFIC OPTIONS
		#.init.nmf.options()
				
		## 1. INITIALIZE THE NMF MODELS
		.init.nmf.models()		
		
		## 2. INITIALIZE BIOC LAYER
		b <- body(.onLoad.nmf.bioc)
		bioc.loaded <- eval(b, envir=pkgEnv)
		if( is(bioc.loaded, 'try-error') )
			packageStartupMessage("NMF - Loading BioConductor layer ... ERROR")
		else if ( isTRUE(bioc.loaded) )
			packageStartupMessage("NMF - Loading BioConductor layer ... OK")
		else{
			packageStartupMessage("NMF - Loading BioConductor layer ... NO [missing Biobase]")
			packageStartupMessage("  To enable, try: install.extras('NMF') [with Bioconductor repository enabled]")
		}
		
		# 3. SHARED MEMORY
		if( .Platform$OS.type != 'windows' ){
			packageStartupMessage("NMF - Checking shared memory capabilities ... ", appendLF=FALSE)
			msg <- if( !require.quiet('bigmemory', character.only=TRUE) ) 'bigmemory'
					else if( !require.quiet('synchronicity', character.only=TRUE) ) 'synchronicity'
			
			if( is.null(msg) ) packageStartupMessage('OK')
			else{
				packageStartupMessage(paste('NO [missing ', msg, ']', sep=''))
				packageStartupMessage("  To enable, try: install.extras('NMF')")
			}
		}
		#
	}
		
	# run intialization sequence suppressing messages or not depending on verbosity options
	.init.sequence()
#	if( getOption('verbose') ) .init.sequence()
#	else suppressMessages(.init.sequence())
	
	
	return(invisible())
}

.onUnload <- function(libpath) {
	
	# TODO: pkgmaker::onUnload(libpath)
	# unload compiled library
	dlls <- names(base::getLoadedDLLs())
	if ( 'NMF' %in%  dlls )
		library.dynam.unload("NMF", libpath);	
}

.onAttach <- function(libname, pkgname){
}

