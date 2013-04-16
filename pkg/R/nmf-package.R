#' @import rngtools
#' @import digest
#' @import stringr
#' @import stats
NULL
library(digest)
library(pkgmaker)

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

# local config info
nmfConfig <- mkoptions()

.onLoad <- function(libname, pkgname=NULL) {
	
	pkgEnv <- pkgmaker::packageEnv()
		
	# set default number of cores
	if( pkgmaker::isCHECK() ){
		options(cores=2)
	}else{
		if( nchar(nc <- Sys.getenv('_R_NMF_CORES_')) > 0 ){
			try({
				nmf.options(cores=as.numeric(nc))
			})
		}
	}
	
	.init.sequence <- function(){
	
		## 0. INITIALIZE PACKAGE SPECFIC OPTIONS
		#.init.nmf.options()
				
		## 1. INITIALIZE THE NMF MODELS
		.init.nmf.models()		
		
		## 2. INITIALIZE BIOC LAYER
		b <- body(.onLoad.nmf.bioc)
		bioc.loaded <- eval(b, envir=pkgEnv)
		nmfConfig(bioc=bioc.loaded)
#		if( is(bioc.loaded, 'try-error') )
#			message("NMF - Loading BioConductor layer ... ERROR")
#		else if ( isTRUE(bioc.loaded) )
#			message("NMF - Loading BioConductor layer ... OK")
#		else{
#			message("NMF - Loading BioConductor layer ... NO [missing Biobase]")
#			message("  To enable, try: install.extras('NMF') [with Bioconductor repository enabled]")
#		}
		
		# 3. SHARED MEMORY
		if( .Platform$OS.type != 'windows' ){
#			message("NMF - Checking shared memory capabilities ... ", appendLF=FALSE)
			msg <- if( !require.quiet('bigmemory', character.only=TRUE) ) 'bigmemory'
					else if( !require.quiet('synchronicity', character.only=TRUE) ) 'synchronicity'
					else TRUE
			
			nmfConfig(shared.memory=msg)
#			if( isTRUE(msg) ) message('OK')
#			else{
#				message(paste('NO [missing ', msg, ']', sep=''))
#				message("  To enable, try: install.extras('NMF')")
#			}
		}
		#
	}
		
	# run intialization sequence suppressing messages or not depending on verbosity options
	.init.sequence()
	if( getOption('verbose') ) .init.sequence()
	else suppressMessages(.init.sequence())
	
	
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
	
	# build startup message
	msg <- NULL
	details <- NULL
	## 1. CHECK BIOC LAYER
	bioc.loaded <- nmfConfig('bioc')[[1L]]
	msg <- paste0(msg, 'BioConductor layer')
	if( is(bioc.loaded, 'try-error') ) msg <- paste0(msg, ' [ERROR]')
	else if ( isTRUE(bioc.loaded) ) msg <- paste0(msg, ' [OK]')
	else{
		msg <- paste0(msg, ' [NO: missing Biobase]')
		details <- c(details, "  To enable the Bioconductor layer, try: install.extras('NMF') [with Bioconductor repository enabled]")
	}
	
	# 2. SHARED MEMORY
	msg <- paste0(msg, ' | Shared memory capabilities')
	if( .Platform$OS.type != 'windows' ){
		conf <- nmfConfig('shared.memory')[[1L]]
		if( isTRUE(conf) ) msg <- paste0(msg, ' [OK]')
		else{
			msg <- paste0(msg, ' [NO: ', conf, ']')
			details <- c(details, "  To enable shared memory capabilities, try: install.extras('NMF')")
		}
	}else msg <- paste0(msg, ' [NO: windows]')
	#
	
	# 3. NUMBER OF CORES
	msg <- paste0(msg, ' | Cores ', getMaxCores(), '/', getMaxCores(limit=FALSE))
	#
	
	# FINAL. CRAN FLAG
	if( pkgmaker::isCHECK() ){
		msg <- paste0(msg, ' | CRAN check')
	}
	#
	packageStartupMessage('NMF - ', msg)
	if( !is.null(details) ){
		packageStartupMessage(paste(details, collapse="\n"))
	}
}

