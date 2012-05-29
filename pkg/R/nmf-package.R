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
	
	.init.sequence <- function(){
	
		## 0. INITIALIZE PACKAGE SPECFIC OPTIONS
		#.init.nmf.options()
		
		## 1. INITIALIZE THE INTERNAL REGISTRY
		.init.nmf.registry()
		
		## 2. INITIALIZE THE NMF MODELS
		.init.nmf.models()		
		
		## 3. INITIALIZE BIOC LAYER
		b <- body(.onLoad.nmf.bioc)
		bioc.loaded <- eval(b, envir=pkgEnv)
		if( is(bioc.loaded, 'try-error') )
			packageStartupMessage("NMF - loading BioConductor layer ... ERROR")
		else if ( isTRUE(bioc.loaded) )
			packageStartupMessage("NMF - loading BioConductor layer ... OK")
		else
			packageStartupMessage("NMF - loading BioConductor layer ... SKIPPED")
	}
	
	## LOAD compiled library if one is loading the package
	if( !missing(libname) ) library.dynam('NMF', pkgname, libname)
	else compile_src('NMF')	
	
	# run intialization sequence suppressing messages or not depending on verbosity options
	if( getOption('verbose') ) .init.sequence()
	else suppressMessages(.init.sequence())
	
	# Initialize the package: load NMF algorithms, seeding methods, etc...
	.init.package()
	
	return(invisible())
}

.onUnload <- function(libpath) {
	
	# unload compiled library
	dlls <- names(base::getLoadedDLLs())
	if ( 'NMF' %in%  dlls )
		library.dynam.unload("NMF", libpath);	
}

.onAttach <- function(libname, pkgname){
}

.init.package <- function(){
	
	.init.sequence <- function(){
				
		## 1. BUILT-IN NMF SEED METHODS
		# initialize built-in seeding methods
		#.load.seed.base() # base: none, random
		#.load.seed.nndsvd() # Non-Negative Double SVD
		#.load.seed.ica() # Positive part of ICA
			
		## 2. POPULATE THE REGISTRY WITH BUILT-IN NMF METHODS: ALGORITHMS
		## NB: this has been disactivated in version 0.7.1 
#		.init.nmf.plugin.builtin()
		
		## 2. USER-DEFINED NMF METHODS
#		.init.nmf.plugin.user()
	}
	
	# run intialization sequence suppressing messages or not depending on verbosity options
	if( getOption('verbose') || nmf.getOption('debug') ) .init.sequence()
	else suppressMessages(.init.sequence())	
		
	return(invisible())
}

# Define a super class for all the NMF pluggable elements: strategies, seeding methods, etc...
#setClassUnion('NMFPlugin', c('NMFStrategy', 'NMFSeed'))
isNMFPlugin <- function(x){
	#print(class(x))
	is(x, 'NMFStrategy') || is(x, 'NMFSeed')	
}

###% Internal function to populate the registry with the built-in methods
.init.nmf.plugin.builtin <- function(){
	
	message('NMF: Init built-in plugins')
		
	# the search is performed within the package's namespace/dev environment
	where <- pkgmaker::packageEnv()
	
	# reset algorithm registry (to be sure)
	return( .init.nmf.plugins(where) )
}

# Loads NMF plugins defined in a given environment
.init.nmf.plugins <- function(envir=topenv()){
	
	# the search is performed within the package's namespace/dev environment
	where <- envir
	
	# lookup for all objects that match the correct pattern
	prefix <- "\\.nmf\\.plugin\\."
	pattern <- paste("^", prefix, ".*", sep='')
	load.fun <- ls(where, all.names=TRUE, pattern=pattern)
	load.plugin.name <- sub(paste("^", prefix, "(.*)", sep=''), "\\U\\1", load.fun, perl=TRUE)
	
	# execute all the loading functions
	strat.list <- NULL
	mapply(function(funname, plugin.name){
				# do something only if the name corresponds to a function
				if( !is.null( fun <- getFunction(funname, mustFind=FALSE, where=where)) ){
					
					message("# loading object(s) from ", plugin.name, ' ... ', appendLF=FALSE)
					strats <- try( fun(), silent=TRUE) 
					
					# wrap the result into a list
					if( isNMFPlugin(strats) )
						strats <- list(strats)
					
					# if the plugin returns NULL then do nothing
					if( is.null(strats) ){
						message('DISABLED')
					}
					# otherwise one should have a list of NMFStrategy objects
					else if( !is.list(strats) || !all(sapply(strats, isNMFPlugin)) ){
						warning("NMF package: unable to load built-in plugin ", plugin.name, " [error: invalid result returned by '", funname,"']", call.=FALSE)
						print(strats)
						message('ERROR')
					}
					else{
						# add the strategies to the list to register
						strat.list <<- c(strat.list, strats)
						message('OK')
					}
					
				}##END if
			}
			, load.fun, load.plugin.name)
	
	# register all the strategies defined in the list
	message("# registering all plugin objects ... ")
	lapply(strat.list, 
			function(s){
				# For the moment only NMFStrategies can be plugged
				if( !is(s, 'NMFStrategy') ){
					warning("NMF-package: seeding method plugin is not yet implemented [object '",name(s),"' SKIPPED]", call.=FALSE)
					return()
				}
				err <- try( nmfRegisterAlgorithm(s, overwrite=TRUE), silent=TRUE )
				if( is(err, 'try-error') ){
					warning("NMF package: unable to register built-in strategy : ", name(s)," [error: ", err,"]", call.=FALSE)
					message('ERROR')
				}
			}
	)
	message('DONE')
	
	invisible()
}

###% Hook to initialize user-defined methods when the package is loaded 
.init.nmf.plugin.user <- function(){
	# TODO: load some RData file stored in the user's home R directory	
}
