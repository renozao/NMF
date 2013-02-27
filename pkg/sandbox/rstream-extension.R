# rstream class for user-supplied RNGs
# 
# Author: Renaud Gaujoux
# Creation: 07 Jul 2011
###############################################################################

# Define packageVersion from Biobase if it does not exist
if( !existsFunction('packageVersion') ){
	packageVersion <- function (pkg, lib.loc = NULL) 
	{
		curWarn <- getOption("warn")
		on.exit(options(warn = curWarn), add = TRUE)
		options(warn = -1)
		desc <- packageDescription(pkg, lib.loc, "Version")
		if (is.na(desc)) 
			stop(paste("Package", pkg, "does not exist"))
		desc
	}
	
}

###% Initialize the RNG wrapper library
###% 
###% If an RNG library is already loaded then the current provider is set to 
###% point to its hooks. Otherwise it will point to `rstream`.  
###% 
.init.RNG <- function(pkgname, libname, libspath){
		
	# enable rngtools
	options(rngtools=TRUE)
	# setup safety command on.exit
	on.exit(options(rngtools=FALSE))

	# initialize .Random.seed if necessary
	if (!exists(".Random.seed", envir=.GlobalEnv)) 
		sample(NA)
	
	# identify the RNG wrapper library name
	message("NMF::.init.RNG - Detect RNG layer library ... ", appendLF=FALSE)
	libfile <- list.files(dirname(libspath), pattern=paste("\\", .Platform$dynlib.ext, "$", sep=''), full.names=TRUE) 
	libfile <- libfile[basename(libfile)!=basename(libspath)]
	f <- basename(libfile) 
	if( length(f) != 1 ){
		message("SKIP [not found]")
		return()
	}
	f <- sub(paste("(.*)\\", .Platform$dynlib.ext, "$", sep=''), "\\1", f)		
	message("OK [", f, "]")
		
	# load the library (in devmode the library is loaded at compile time)
	message("NMF::.init.RNG - Initialise RNG wrapper library '", f, "' ...")
	
	# set static values for the Libname and libpath 
	RNGwrapper.Libname(f)
	RNGwrapper.Libpath(libfile)
	# ensure that the RNGwrapper library is not loaded
	RNGwrapper.unload()
	# call RNGlib to try freezing the RNG lib to rstream on Windows
	RNGlib()
	
	# enable rngtools depending on the OS
	#options(rngtools=(.Platform$OS.type != 'windows'))
		
	getRNG()
	# remove safety command from on.exit
	on.exit()
		
	message("### DONE ###")
}

###% Returns the name of the current internally stored RNG provider
###% 
###% It returns the name stored in the C++ static variable _current_provider 
###% 
rngtools.get.provider <- function(fixup=FALSE){
	
	if( !getOption('rngtools.lib') )
		return( '' )
	
	p <- .Call('rngtools_getProvider')
	
	# fixup if necessary: 
	if( fixup && p == '' ){
		p <- RNGlibs(1)
		rngtools.set.provider(p)
		warning("rngtools - Fixing unspecified RNG provider: using `", p, "`")
	}
	p
}

###% Tells if the current internal RNG provider can be directly restored 
###% via .Random.seed.
###% 
###% This is TRUE for base RNG and user-supplied RNG that provides a hook for 
###% unif_seedloc and unif_nseed.
###% 
rstream.is.restorable <- function(){
	if( baseRNGkind()[1] != "user-supplied" )
		TRUE
	else{
		prov <- rngtools.get.provider()
		if( prov != '' )			
			(!is(try(getNativeSymbolInfo("user_unif_nseed", PACKAGE=prov), silent=TRUE), 'try-error')
			&& !is(try(getNativeSymbolInfo("user_unif_seedloc", PACKAGE=prov), silent=TRUE), 'try-error'))
		else
			FALSE
	}
}

###% Sets the current internally stored RNG provider
###% 
###% It populates the C++ static variable _current_provider with hook pointers 
###% obtained from provider \code{name}. 
###% 
rngtools.set.provider <- function(name, reload.hooks=FALSE, verbose=FALSE){
	
	if( !getOption('rngtools.lib') )
		return( '' )
	
	# use the stream user-RNG provider if name is an rstream object
	if( is(name, 'rstream') )
		name <- rstream.user.provider(name)
	
	if( name == '' )
		stop("cannot set an empty current RNG provider ")
	
	if( isTRUE(verbose) )
		verbose <- "Setting user-supplied RNG provider to"
	if( is.character(verbose) )
		message("# ", verbose, " '", name, "' ... ", appendLF=FALSE)
		
	res <- .Call('rngtools_setProvider', as.character(name)[1], reload.hooks)
	
	if( is.character(verbose) )
		message("OK")
	
	res
}

###% Returns the name of the next current internally stored RNG provider
###% 
###% It returns the name stored in the C++ static variable _next_provider 
###% 
rstream.get.nextprovider <- function(){
	
	if( !getOption('rngtools.lib') )
		return( '' )
	
	.Call('rngtools_getNextProvider')
}

###% Sets the next current internally stored RNG provider
###% 
###% It populates the C++ static variable _next_provider with hook pointers 
###% obtained from provider \code{name}. 
###% 
rngtools.set.nextprovider <- function(name){
	
	if( !getOption('rngtools.lib') )
		return( '' )
	
	# use the stream user-RNG provider if name is an rstream object
	if( is(name, 'rstream') )
		name <- rstream.user.provider(name)
	
	.Call('rngtools_setNextProvider', as.character(name)[1])
}

###% Extracts and Updates the Current RNG as an rstream Object
###% 
###% This S4 generic function aims at substituting the non-exported S4 function
###% .rstream.getRNG from the package rstream.
###% It allows to define classes that can fully integrate into the framework 
###% of rstream.
###% 
###% @param stream the current registered rstream object
###% 
setGeneric("rstream.GetRNGstate", function(stream, ...) standardGeneric("rstream.GetRNGstate"))

###% Default method for rstream objects: clone the current stream
###% 
###% This method should work for all rstream classes that store the state 
###% information into an external pointer, that is dynamically updated on the C 
###% side, and copied by the method rstream.clone 
###% 
setMethod("rstream.GetRNGstate", "rstream", 
	function(stream, ...){
		rstream.clone(stream)
	}
)

###% Method for class rstream.runif that return the current R built-in RNG.
###% 
###% The default method that clones the object does not work here, as its state is
###% NOT dynamically updated, together with .Random.seed.
###% This is why one needs to return a new object of rstream.runif based on the 
###% current value of .Random.seed. 
###%  
setMethod("rstream.GetRNGstate", "rstream.runif", 
	function(stream, dest, ...){
		res <- new('rstream.runif')
		# add the user-supplied provider info
		rstream.user.provider(res) <- rstream.user.provider(stream) 
		res
	}
)

###% Set Current RNG from an rstream Object
###% 
###% This S4 generic function aims at substituting the non-exported S4 function
###% rstream.setRNG from the package rstream.
###% It allows to define classes that can fully integrate into the framework 
###% of rstream.
###% 
###% @param stream the rstream object to register as the current RNG
###% 
setGeneric("rstream.PutRNGState", function(stream, ...) standardGeneric("rstream.PutRNGState"))

###% Method for setting an rstream.runif object as the RNG
###% 
###% 
###% 
setMethod("rstream.PutRNGState", "rstream.runif", 
	function(stream, ...){
		
		if ( rstream.antithetic(stream)) {
			warning ("NMF::rstream.PutRNGState(",class(stream)[1],") - antithetic DISABLED")
			rstream.antithetic(stream) <- FALSE 
		}
		# do NOT use RNGkind as it would draw once from the current RNG, 
		# all the required RNG information is already stored in .Random.seed  
		# RNGkind(kind=stream@kind)
		state <- get("state", envir=stream@xstate)
		assign(".Random.seed", state, envir=.GlobalEnv)
		stream		  
	}
)

###% Method for setting an rstream.mrg32k3a object as the RNG
###%
setMethod("rstream.PutRNGState", "rstream.mrg32k3a",
	function(stream, ...){
		# do NOT use RNGkind as it would draw once from the current RNG, 
		# it has been already called before rstream.PutRNGState
		.Call("R_RngStreams_setRNG", stream@stream, PACKAGE="rstream")
		stream
	}
)
		
# backup original base RNGkind
if( !exists('baseRNGkind', topenv()) ){
	baseRNGkind <- base::RNGkind
}

###% Returns the name of the RNG hook wrapper library
###% 
###% The correct name is detected at load time in .onLoad.
###% The actual name is defined in Makevars.
###% 
RNGwrapper.Libname <- function(){
	.libname <- NULL
	function(libname){
		if( !missing(libname) ){
			.libname <<- libname
			dlls <- getLoadedDLLs()
			RNGwrapper.Libpath(dlls[[libname]][['path']])
		}
		
		.libname 
	}
}
RNGwrapper.Libname <- RNGwrapper.Libname()

RNGwrapper.Libpath <- function(){
	.libpath <- NULL
	function(libpath){
		if( !missing(libpath) )
			.libpath <<- libpath
		
		.libpath 
	}
}
RNGwrapper.Libpath <- RNGwrapper.Libpath()

isRNGwrapper <- function(libname){
	if( missing(libname) )
		libname <- RNGlib()
	# test equality
	libname == RNGwrapper.Libname()
}
RNGwrapper.isActive <- function() isRNGwrapper()
 

RNGwrapper.unload <- function(){
	
	res <- NULL
	if ( RNGwrapper.Libname() %in%  names(base::getLoadedDLLs()) )
		res <- dyn.unload(RNGwrapper.Libpath())	
	options(rngtools.lib=FALSE)
	invisible(res)
}

RNGwrapper.load <- function(){
		
	res <- dyn.load(RNGwrapper.Libpath())
	options(rngtools.lib=TRUE)
	# call RNGlib to possibly freeze the cache on RNGwrapper
	RNGlib()
	# return DLLInfo
	invisible(res)
}

#RNGwrapper.updateProvider <- function(provider){
#	
#	update.current <- FALSE
#	verbose <- FALSE
#	if( missing(provider) ){
#		update.current <- TRUE
#		provider <- rngtools.get.provider()
#	}
#	
#	# restore/fix stored RNG provider
#	libs <- RNGlibs()
#	# fixup if the current provider is not loaded anymore
#	if( !is.element(provider, libs) ){	
#		nextprov <- RNGlibs(1)
#		if( provider != '' )
#			warning("RNGwrapper - RNG provider '", provider, "' not found: "
#					, "wrapping last loaded RNG provider '", nextprov, "' instead")			
#		provider <- nextprov
#		if( update.current ) verbose <- TRUE 
#	}
#	# reset the provider
#	if( provider != '' )	
#		rngtools.set.provider(provider, verbose=verbose)	
#}


# Define path.package from R-2.13.1 if not defined already
if( !existsFunction('path.package') ){
	
	path.package <- function (package = NULL, quiet = FALSE) 
	{
		if (is.null(package)) 
			package <- .packages()
		if (length(package) == 0L) 
			return(character())
		s <- search()
		searchpaths <- lapply(seq_along(s), function(i) attr(as.environment(i), 
							"path"))
		searchpaths[[length(s)]] <- system.file()
		pkgs <- paste("package", package, sep = ":")
		pos <- match(pkgs, s)
		if (any(m <- is.na(pos))) {
			if (!quiet) {
				if (all(m)) 
					stop("none of the packages are loaded")
				else warning(sprintf(ngettext(as.integer(sum(m)), 
											"package %s is not loaded", "packages %s are not loaded"), 
									paste(package[m], collapse = ", ")), domain = NA)
			}
			pos <- pos[!m]
		}
		unlist(searchpaths[pos], use.names = FALSE)
	}
	
}

RNGwrapper.enforce <- function(){
	
	# do not do anything if RNGwrapper is already active
	if( isRNGwrapper() )
		return( TRUE )	
	
	if( ! RNGwrapper.Libname() %in%  names(base::getLoadedDLLs()) ){
		
		# if RNGwrapper was not loaded then we'll synchronise it with the current  
		# RNG settings by wrapping the last loaded RNG lib
		RNGwrapper.load()
		#provider <- RNGlibs(1)
		
	}else{ # unload/reload the RNG wrapper library so that it takes precedence over 
		# the other RNG libraries		
		
		# store currently wrapped RNG provider for restoration after reloading
		#provider <- rngtools.get.provider()
		
		RNGwrapper.unload()
		RNGwrapper.load()
		
	}	
	
	# check if the RNGwrapper library is a valid RNG library and unload it if not 	
	if( !is.RNGlib(RNGwrapper.Libname()) ){
		options(rngtools=FALSE)
		warning("rngtools::RNGwrap - RNG wrapper is not a valid RNG library (probably due to too old gcc version): wrapping is now disabled.")		
		RNGwrapper.unload()
		# return FALSE for failure
		return( FALSE )
	}
	
	# On non Windows systems this should be enough to enforce
	if( .Platform$OS.type != 'windows' && !isRNGwrapper() )
		warning("rngtools::RNGwrap - Enforcing RNGwrapper library for non Windows OS did not work as expected.")
	
	if( !isRNGwrapper() ){
		# on Windows there is a DLL caching issue
		# One try to enforce RNGwrapper by unloading the libraries that 
		# could take precedence on it. Currently only the rstream 
		# library is handled.
		
		# get all loaded RNG libraries (including RNGwrapper)
		libs <- RNGlibs(all=TRUE)
		stopifnot( tail(libs, 1) == RNGwrapper.Libname() )
		# if rstream is not initialised yet then try to enforce the 
		# RNGwrapper library by unloading/reloading rstream
		if( 'rstream' %in% libs && !is.rstream.initialised() ){
			
			# unload rstream
			libpath <- path.package('rstream')
			message("\n# Unloading rstream library from '", libpath, "' ... ", appendLF=FALSE)					
			library.dynam.unload('rstream', libpath)
			message("OK")
			
			# call RNGlib to hopefully freeze the cache on RNGwrapper
			RNGlib()
			
			# reload rstream
			message("# Reloading rstream library ... ", appendLF=FALSE)
			library.dynam('rstream', 'rstream')
			message("OK\n")
			
		}
		
	}	
	
	# if the current active RNG library is still not RNGwrapper then exit with a warning
	if( !isRNGwrapper() ){		
		# one cannot enforce rngtools 
#		options(rngtools=FALSE)
#		warning("rngtools::RNGwrap - Could not enforce RNG wrapper library: wrapping is now disabled.")
		return( FALSE )		
	}
		
	# return TRUE for success
	TRUE
}

###% Returns the library that provides the current user-supplied RNG hooks.
###% 
###% This is the library that is first called by runif when using setting RNG 
###% kind to "user-supplied".
###% In general this will be rstream, except if a package providing the RNG hook 
###% 'user_unif_rand' is loaded after rstream, and no call to RNGkind or getRNG 
###% were done thereafter.
###% 
###% @return an object of class NativeSymbolInfo or NULL if no hook were found
###% 
RNGlib <- function(PACKAGE='', full=FALSE, hook="user_unif_rand", ...){
	
	if( !missing(PACKAGE) )
		full = TRUE
	if( !missing(hook) )
		hook <- match.arg(hook, c('user_unif_rand', 'user_unif_init', 'user_unif_nseed', 'user_unif_seedloc'))
	
	# lookup for the hook "user_unif_rand" in all the loaded libraries
	symb.unif_rand <- try( getNativeSymbolInfo(hook, PACKAGE=PACKAGE, ...), silent=TRUE)
	if( is(symb.unif_rand, 'try-error') ){
		
		if( !full ) '' else NULL
		
	}else if( PACKAGE=='' && is.null(symb.unif_rand$package) ){ 
		#special case for MS Windows when PACKAGE is not specified: if two 
		# RNGlibs are loaded, the first one is seen, not the last one as on Unix
		libs <- RNGlibs(all=TRUE, full=TRUE, unlist=FALSE, hook=hook)
		w <- which(sapply(libs, function(l) identical(l$address, symb.unif_rand$address)))
		
		# returns full info or just the name
		if( full ) libs[[w]]
		else names(libs)[w]
		
	}else if( full ) symb.unif_rand
	else symb.unif_rand$package[['name']]
}

###% Returns all the libraries that provides a user-supplied RNG
###% 
###% The library that provides the wrapper hooks for the management multiple 
###% user-supplied RNG is removed from the output list.
###% 
RNGlibs <- function(n=0, all=FALSE, full=FALSE, hook="user_unif_rand", unlist=TRUE){
	dlls <- getLoadedDLLs()
	res <- lapply(dlls, function(d){
			dname <- d[['name']]
			if( dname=='' || (!all && dname == RNGwrapper.Libname()) )
				return(NA)
			
			symb.unif_rand <- RNGlib(PACKAGE=dname, hook=hook)
			if( is.null(symb.unif_rand) )
				NA
			else
				symb.unif_rand
		})
		
	res <- res[!is.na(res)]
	if( !full )
		res <- names(res)	
	
	# limit the results if requested
	if( n>0 )
		res <- tail(res, n)
	
	# return result
	if( unlist && length(res) == 1 )
		res[[1]]
	else
		res
}

###% Tells if a library or package provides RNG hooks
###% 
###% 
is.RNGlib <- function(lib){
	
	# lookup for the hook "user_unif_rand" in all the loaded libraries
	!is.null( RNGlib(PACKAGE=lib) )	
}

###% Returns the package that provides the current RNG managed by rstream
###% 
###% It returns the name of the package to which are currently passed the RNG 
###% calls (runif, set.seed).
###% This is either 'base' if core RNG is in use (e.g. Mersenne-Twister, Marsaglia-Multicarry, etc...) 
###% or the package that provides the actual RNG hooks called by the rstream 
###% wrapper hooks. This one was set either explicitly via RNGkind or implicitly 
###% when rstream was first loaded. In this latter case, the provider was identified 
###% at loading time as 'base' if core RNGs were in use or as the package that was 
###% providing the RNG hook 'user_unif_rand' if the RNG in used was "user-supplied".       
###%
RNGprovider <- function(user.supplied=FALSE){
	
	kind <- baseRNGkind()
	if( kind[1] == 'user-supplied' || user.supplied ){
		rlib <- RNGlib() 
		if( rlib == RNGwrapper.Libname() ) 
			rngtools.get.provider()
		else
			rlib
	}
	else
		'base'
}

RNGproviderInfo <- function(object, user.supplied=FALSE, version=FALSE){
	
	# get the RNG provider's name
	pkg <- if( missing(object) ) RNGprovider(user.supplied = user.supplied)
		   else if( user.supplied ) rstream.user.provider(object)
   		 	else rstream.provider(object)
			   
		
	
	# add version to the package name if requested
	if( version && length(pkg) > 0 ){
		vers <- try( packageVersion(pkg), silent=TRUE )
		if( is(vers, 'try-error') )
			vers <- NULL
		pkg <- c(pkg, as.character(vers))
		pkg <- paste(pkg, collapse='_')
	}
	
	# return info
	pkg
}

RNGinfo <- function(object, prefix=''){
	if( missing(object) ){
		cat(prefix, "RNG kind: ", paste(baseRNGkind(), collapse=" / "), " [", RNGproviderInfo(user.supplied=TRUE, version=TRUE), "]\n")								
		cat(prefix, "RNG state: ", RNGdesc(), "\n")
	}else{
		object <- getRNG(object)
		kinds <- c("Wichmann-Hill", "Marsaglia-Multicarry", "Super-Duper", 
				"Mersenne-Twister", "Knuth-TAOCP", "user-supplied", "Knuth-TAOCP-2002", 
				"default")
		 
		prov <- rstream.provider(object)
		prov <- if(prov == "base") kinds[object@seed[1] %% 100 + 1]
				else paste('pkg:', prov, sep='')
		cat(prefix, "RNG kind: ", paste(prov, baseRNGkind()[2], sep=" / "), " [", RNGproviderInfo(object, user.supplied=TRUE, version=TRUE), "]\n")								
		cat(prefix, "RNG state:", RNGdesc(object), "\n")
	}
} 


RNGrecovery <- function(){
	
	s <- as.integer(c(401,0,0))
	assign(".Random.seed", s, envir=.GlobalEnv)
	RNGkind("default")
	
}

RNGcurrent <- function(){
	if( exists(".rstream.current", envir=rstream:::.rstream.envir) )
		get(".rstream.current", envir=rstream:::.rstream.envir)
	else
		NULL	
}

# check if the rstream library is initialised
is.rstream.initialised <- function(){
	s <- try(.Call("R_RngStreams_GetPackageSeed", PACKAGE="rstream"), silent=TRUE)
	!is(s, 'try-error') && !identical(as.numeric(s), rep(12345, 6))
}

###% Reloads the RNGwrapper library
###% 
###% The default behaviour is to reset the wrapped RNG provider to the last loaded
###% RNG provider if the current RNG is provided by R-base, or to the previously 
###% wrapped RNG provider if RNGwrapper was already used to provide the current 
###% user-supplied RNG, or to the last loaded
###% 
#RNGwrap <- function(provider, force=FALSE){
#		
#	# check if the wrapper hook library is loaded
#	isLoaded <- RNGwrapper.Libname() %in% names(getLoadedDLLs())
#	# get the current RNG library
#	rlib <- RNGlib()
#	
#	# if the RNGwrapper library is not valid or disabled: skip everything
#	if( (!force && !getOption('rngtools')) || (isLoaded && !is.RNGlib(RNGwrapper.Libname())) ){
#		options(rngtools=FALSE)
#		return( rlib )	
#	}
#		
#	# if the last loaded RNG library is not RNGwrapper then (re-)load it
#	if( rlib != RNGwrapper.Libname() || force ){
#		
#		# Reload the RNGwrapper library
#		message("# Reloading the RNG wrapper library ... ", appendLF=FALSE)			
#		if( RNGwrapper.reload() ) message("OK")
#		else{			
#			options(rngtools=FALSE)
#			message("STOP")
#			return( rlib )
#		}
#		# From here one knows that RNGwrapper is the active or last loaded RNG library		
#		
#		# try to enforce the RNGwrapper library if it is not already active		
#		if( !RNGwrapper.isActive() ){
#
#			# On nix systems the RNGlibrary should be already enforced
#			if( .Platform$OS.type != 'windows' ){				
#				stop("rngtools::RNGwrap - Unexpected error for non Windows OS: RNGwrapper library should be the active RNG library")				
#			}else{ 
#				# on Windows there is a DLL caching issue
#				# One try to enforce RNGwarpper by unloading the libraries that 
#				# could take precedence on it. Currently only the rstream 
#				# library is handled.
#				
#				# get all loaded RNG libraries (including RNGwrapper)
#				libs <- RNGlibs(all=TRUE)
#				stopifnot( tail(libs, 1) == RNGwrapper.Libname() )
#				# if rstream is not initialised yet 
#				# AND RNGwrapper is the second one, then one enforces the 
#				# RNGwrapper library by unloading/reloading rstream		
#				if( 'rstream' %in% names(libs) ){
##					[1]=='rstream' && !is.rstream.initialised()
##					&& names(libs)[2]==RNGwrapper.Libname() ){
#					
#					# unload rstream
#					libpath <- find.package('rstream')
#					message("# Unloading rstream library from '", libpath, "' ... ", appendLF=FALSE)					
#					library.dynam.unload('rstream', libpath)
#					message("OK")
#					
#					# call RNGlib to hopefully freeze the cache on RNGwrapper
#					RNGlib()
#					
#					# reload rstream
#					message("# Reloading rstream library ... ", appendLF=FALSE)
#					library.dynam('rstream', 'rstream')
#					message("OK")
#										
#				}
#			}
#			
#			# if the current active RNG library is still not RNGwrapper then 
#			# exit with a warning
#			if( !RNGwrapper.isActive() ){
#				
#				# one cannot enforce rngtools 
#				options(rngtools=FALSE)
#				warning("rngtools::RNGwrap - Could not enforce RNG wrapper library: wrapping is now disabled.")
#				return( RNGlib() )
#				
#			}
#			
#		}
#		
#		# From here: all static variables in RNGwrapper are reset to there default
#		# values, in particular _current_provider and _next_provider are empty
#		# strings
#		
#		# if the current RNG is a core RNG then simply set the current user RNG 
#		# provider to the last loaded RNG library: this keeps the RNGwrapper library 
#		# synchronised with the current settings, and RNGkind will act as normally 
#		# expected
#		if( baseRNGkind()[1] != 'user-supplied' ){
#			# If no specific provider is requested, the library to use 
#			# is presumably the last loaded RNG library so we update the internal
#			# provider to point to it.
#			if( missing(provider) )
#				provider <- RNGlibs(1)
#						
#			rngtools.set.provider(provider, TRUE, verbose=TRUE)			
#		}else{	
#			
#			if( missing(provider) ){
#				# Fixes the potential issue of the cached pointers to the RNG hooks.
#				# A special care is taken for the case where the current active RNG 
#				# (a user-supplied RNG) was called through RNGwrapper.
#				# rstream_postReloadFix returns the detected active RNG provider
#				message("# Auto-wrapping the active user-supplied RNG provider ... ", appendLF=FALSE)
#				prevprov <- .Call('rngtools_detectProvider')
#				if( is.null(prevprov) )
#					message("FAILED")
#				else
#					message("OK [", prevprov, "]")
#				
#			}else{
#				rngtools.set.provider(provider, TRUE, verbose="Switching to RNG provider")
#			}
#		}
#					
#	}else{ # RNGwrapper is the active RNG library
#		
#		# if a specific provider is requested, then wrap it 
#		if( !missing(provider) ){
#			vmsg <- if( baseRNGkind()[1] != 'user-supplied' ) "Setting user-supplied RNG provider to"				
#					else "Switching to RNG provider"
#			
#			# set the RNG provider
#			rngtools.set.provider(provider, TRUE, verbose=vmsg)			
#		}else
#			RNGwrapper.updateProvider()
#	}
#	
#	# everything is working: toggle ON rngtools 
#	options(rngtools=TRUE)
#	
#	# return currently wrapped provider
#	rngtools.get.provider()
#}

RNGhooks <- function(provider='', who=FALSE){
	if( missing(provider) && who ){
		# user_norm_rand
		sapply(c('user_unif_rand', 'user_unif_init', 'user_unif_nseed', 'user_unif_seedloc') 
				, function(s) RNGlib(hook=s))
	}else
		sapply(c('user_unif_rand', 'user_unif_init', 'user_unif_nseed', 'user_unif_seedloc'), is.loaded, provider)
}

checkRNGhooks <- function(provider, wrapper=TRUE, warn=FALSE){
	
	if( missing(provider) )
		provider <- RNGlib()
	
	# check that the RNG hooks are consistent, i.e. from the same provider	
	whooks <- RNGhooks(who=TRUE)
	# hooks that are ok are the ones provided by provider or missing from all RNG libs
	ok <- c('', provider)
	# if wrapping is allowed then resolve hooks provided by the RNGwrapper	
	if( isTRUE(wrapper) )
		whooks[whooks==RNGwrapper.Libname()] <- rngtools.get.provider()
	else if( is.na(wrapper) ) # if NA: do not resolve but allow the hooks from RNGwrapper
		ok <- c(ok , RNGwrapper.Libname())
	# compute missing hooks
	missing.hook <- whooks[ ! (whooks %in% ok) ]
	
	if( length(missing.hook)>0 ){		
		warning("rngtools::checkRNGhooks - inconsistent RNG hooks: usage may crash the R session."
				, " [hooks:", paste(names(whooks)[missing.hook], ':', whooks[missing.hook], collapse=", "),"]")
		return( FALSE )
	}
	
	TRUE	
}

###% Wrapping an RNG library
###% 
###% This function ensures that a given RNG library can be used, even when it is 
###% not the last loaded library.
###% 
RNGwrap <- function(provider, verbose=TRUE, warn=TRUE){
	
	finalCheck <- function(provider, warn){
		# check that all the hooks are correctly resolved
		if( verbose ) message("# Check if all RNG hooks are correctly resolved ... ", appendLF=FALSE)
		if( !checkRNGhooks(provider, warn=warn) ){
			if( verbose ) message("NO")
			return( FALSE )	
		}
		if( verbose ) message("DONE")
		return( TRUE )
	}
	
	## CHECK if the provider could be used directly
	if( verbose ) message("# Check if library '", provider ,"' is loaded ... ", appendLF=FALSE)
	internally.loaded <- FALSE
	if( (! provider %in% names(getLoadedDLLs())) ){
		if( verbose ) message("NO")
		if( verbose ) message("# Load package '", provider ,"' ... ", appendLF=FALSE)
		if( !require.quiet(provider, character.only = TRUE) ){
			if( verbose ) message("ERROR")
			stop("rngtools::RNGwrap - Could not find package '", provider, "'. Check if it is in the search path.")
			return( FALSE )
		}
		if( verbose ) message("OK")
		internally.loaded <- TRUE
	}else if( verbose ) message("YES")
	
	if( verbose ) message("# Check if library '", provider ,"' is a valid RNG library ... ", appendLF=FALSE)
	if( !is.RNGlib(provider) ){
		if( verbose ) message("NO")		
		if( internally.loaded )
			detach(paste('package:', provider, sep=''), unload=TRUE, character.only=TRUE)
		stop("rngtools::RNGwrap - Library '", provider, "' is not a valid RNG provider.")
		return( FALSE )
	}
	if( verbose ) message("YES")
	
	# get active lib
	rlib <- RNGlib()
	# get all other libs
	libs <- RNGlibs()
	olibs <- libs[libs != provider]	
	
	# It is OK to directly use the RNG library only if it provides all the 
	# other RNG hooks or none of the other libraries provide any of them
	can.use.alone <- TRUE
	if( length(olibs) > 0 ){
		
		if( verbose ) message("# Detected multiple RNG libraries: ", paste(libs, collapse=", "))
		
		# define function that checks the presence of hooks accross the libraries
		checkHooks <- function(provider){
			all(sapply(c('user_unif_init', 'user_unif_nseed', 'user_unif_seedloc')
							, function(h){
								is.loaded(h, provider) || !any(sapply(olibs, function(l) is.loaded(h, l)) )
							}))
		}
		
		if( verbose ) message("# Check if RNG library's hooks allow for direct use ... ", appendLF=FALSE)		
		hooks.ok <- checkHooks(provider)
		if( !hooks.ok ){
			if( verbose ) message("NO")
			can.use.alone <- FALSE
		}else if( verbose ) message("YES")
		
	}else if( verbose ) message("# RNG library '", provider, "' is suitable for direct use")
	
	if( verbose ) message("# Check if RNG wrapper library is loaded ... ", appendLF=FALSE)
	wrapperIsLoaded <- is.RNGlib(RNGwrapper.Libname())
	if( verbose ) message(if(wrapperIsLoaded) "YES" else "NO")
	
	if( verbose ) message("# Check if RNG library '", provider ,"' is active ... ", appendLF=FALSE)
	libIsActive <- rlib == provider
	
	if( libIsActive  ){
		
		# unload the wrapper library if it is not necessary
		if( can.use.alone ){
			
			if( wrapperIsLoaded ){
				if( verbose ) message("YES")
				if( verbose ) message("# Unloading uneeded RNG wrapper library ... ", appendLF=FALSE)
				RNGwrapper.unload()
			}
			if( verbose ) message("YES")
			return( finalCheck(provider, warn) )
			
		}else{
			
			if( verbose ) message("YES")
			if( wrapperIsLoaded ){		
				# the requested provider is already the active RNG library: we only
				# point RNGwrapper to the provider so that missing hooks are properly
				# redirected (NB: do not reload the hooks, as we want to directly use the 
				# pointers provided by the RNG library)		
				if( verbose ) message("# Wrapping missing hooks ... ", appendLF=FALSE)
				rngtools.set.provider(provider)
				if( verbose ) message("OK")
				
				return( finalCheck(provider, warn) )				
			}
		}
		# At this point: RNG lib is active but cannot be used alone:
		# the RNG wrapper library needs to be loaded so that the hooks are 
		# correctly resolved
	
	}else{
		
		if( verbose ) message("NO")
		if( can.use.alone && wrapperIsLoaded ){
		
			# If RNGwrapper is the active library and provider is next in the queue
			# => unload RNGwrapper and make sure it activated the correct library
			# This is not OK on Windows where the activation is not guaranteed: 
			# in this case we unload only if these two libs are the only RNG libraries
			if( tail(libs, 1) == provider 
				&& ( (.Platform$OS.type != 'windows' && isRNGwrapper()) 
					|| (.Platform$OS.type == 'windows' && length(libs)==1) ) ){
				
				if( verbose ) message("# Unloading the RNG wrapper library should activate '", provider,"' ... ", appendLF=FALSE)				
				RNGwrapper.unload()
				if( RNGlib() == provider ){
					if( verbose ) message("OK")
					
					return( finalCheck(provider, warn) )
				}
				if( verbose ) message("FAILED")
			}
		}			
		# At this point: the RNG library is not active and although it could be 
		# used alone, one needs the RNG wrapper library to get to it 
	}
	
	## ENFORCE RNGwrapper if necessary
	if( verbose ) message("# Check if RNG wrapper library is active ... ", appendLF=FALSE)
	if( !isRNGwrapper() ){ 
		if( verbose ) message("NO")
		
		# Reload the RNGwrapper library
		if( verbose ) message("# Enforcing RNG wrapper library ... ", appendLF=FALSE)
		if( !RNGwrapper.enforce() ){
			if( verbose ) message("FAILED")
			
			# this check is for windows where caching can produce strange things 
			if( verbose ) message("# Check if RNG library '",provider,"' is now active ... ", appendLF=FALSE)
			rlib <- RNGlib()
			if( rlib == provider ){
								
				# it it can be used alone: we're done 
				if( can.use.alone ){
					
					if( is.RNGlib(RNGwrapper.Libname()) ){
						if( verbose ) message("YES")
						if( verbose ) message("# Unloading uneeded RNG wrapper library ... ", appendLF=FALSE)
						RNGwrapper.unload()
					}
					if( verbose ) message("OK")					
					# final check and return
					return( finalCheck(provider, warn) )
				}else if( verbose ) message("YES")
			
			}else{			
				if( verbose ) message("NO [", rlib, "]")
				return( FALSE )
			}
		}else if( verbose ) message("OK")	
	}else if( verbose ) message("YES")
	##
	
	## WRAP
	# check that all the hooks are possibly correctly resolved
	if( verbose ) message("# Check if all RNG hooks can be resolved ... ", appendLF=FALSE)
	if( !checkRNGhooks(provider, wrapper=NA, warn=warn) ){
		if( verbose ) message("NO")
		return( FALSE )	
	}
	if( verbose ) message("YES")
	
	if( verbose ) message("# Wrapping RNG library '", provider, "' ... ", appendLF=FALSE)
	rngtools.set.provider(provider, TRUE)
	uprovider <- RNGprovider(user.supplied = TRUE)
	if( uprovider == provider ){
		if( verbose ) message("OK")		
		# final check and return 
		return( finalCheck(provider, warn) )
	}
	else if( verbose ) message("FAILED [user-supplied RNG is '", uprovider ,"']")
	
	return( FALSE )	
	##	
}

#testwrap <- function(){
#	source('../pkg/R/package.R', chd=TRUE)
#	try(dyn.unload("../pkg/tmp/urand.so"))
#	try(detach('package:rlecuyer', unload=TRUE))
#	RNGrecovery()
#	getRNG()
#	cat("\n\n\n")	
#	dyn.load("../pkg/tmp/urand.so")
#	#RNGkind("user")
#	#RNGwrap()
#	
#	s <- new('rstream.mrg32k3a', seed=rep(123,6), force=TRUE)
#		
#	setRNG(s)
#	
#	message("Test wrap by missing unif_init")
#	a <- .Random.seed
#	print(head(a))
#	library(rlecuyer)
#	.lec.CreateStream('a')
#	.lec.CurrentStream('a')
#	message("Test runif")
#	print( runif(1) )
#		
#	message("Test restore seed on fake_set_seed")
#	RNGrecovery()
#	a <- .Random.seed
#	print(head(a))
#	RNGwrap('rstream')
#	stopifnot( identical(a,.Random.seed ))
#	print(head(a))
#	message("Test runif")
#	print( runif(1) )
#	
#	message("Test wrap by missing RNGwrap")
#	try(dyn.unload("../pkg/tmp/urand.so"))
#	dyn.load("../pkg/tmp/urand.so")
#	try(detach('package:rlecuyer', unload=TRUE))
#	library(rlecuyer)
#	.lec.CreateStream('a')
#	.lec.CurrentStream('a')	
#	message("RNGwrap NOW")
#	RNGwrap()
#	message("Test runif")
#	print( runif(1) )
#	
#}

###% Fix up the currently stored RNG if it does not exists.
###% 
fixupRNG <- function(silent=FALSE){
	
	current <- RNGcurrent()
	if( is.null(current) # current stream not defined
		# current stream does not match active provider
		|| rstream.provider(current) != RNGprovider()
		# current stream user-supplied RNG does not match current RNG settings 
		|| rstream.user.provider(current) != RNGprovider(user.supplied=TRUE) 
		){
		#message("### Fixing up current RNG")
		if( !silent ) 
			message("Fixing up current RNG")
		assign(".rstream.current", newRNG(), envir=rstream:::.rstream.envir)
		#message("### DONE")
	}
}

###% Directly get/set all internal objects at once
###% 
RNGinternals <- function(value, verbose=FALSE){
	
	if( missing(value) ){ # return the list of all internal objects
		list(.rstream.current = get(".rstream.current", envir=rstream:::.rstream.envir)
			, provider = rngtools.get.provider()
			, nextprovider = rstream.get.nextprovider()
			, .Random.seed = get(".Random.seed", envir=.GlobalEnv)
			)
	}else{
		
		# restore objects
		if( !is.null(value$.Random.seed) ){
			if( verbose )
				message("# Restoring .Random.seed")
			assign(".Random.seed", value$.Random.seed, envir=.GlobalEnv)
		}
		if( !is.null(value$.rstream.current) ){
			if( verbose )
				message("# Restoring current RNG object")
			assign(".rstream.current", value$.rstream.current, envir=rstream:::.rstream.envir)
		}
		if( !is.null(value$provider) ){
			if( verbose )
				message("# Restoring RNG provider")
			rngtools.set.provider(value$provider)
		}
		if( !is.null(value$nextprovider) ){			
			if( verbose )
				message("# Restoring next RNG provider")
			rngtools.set.nextprovider(value$nextprovider)
		}
	}
}

RNGdigest <- function(x){
	
	object <- if( missing(x) )	getRNG() else getRNG(x)
	
	# exit if no RNG was extracted
	if( is.null(object) )
		return(digest(NULL))
	
	# pack the rstream.object if necessary
	if( !rstream.packed(object) ){
		#work on a clone as packing may free allocated ressources (e.g. rstream.mrgk32a3)
		object <- rstream.clone(object)
		rstream.packed(object) <- TRUE
	}
	
	# only digest the invariant part
	p <- object@pack
	# make sure the name does not make the difference
	p$name <- ''
	# add the object's type and the user-supplied RNG provider in the comparison
	cmp <- c(p, ZZZtype=object@type, ZZZprovider=rstream.user.provider(object))
	# order in a consistent way
	cmp <- cmp[order(names(cmp))]
	digest(cmp)
}

rng.equal <- function(x, y){
	if( missing(y) )
		y <- getRNG()
	identical(RNGdigest(x), RNGdigest(y))
}

rng1.equal <- function(x, y){
	if( is(x, 'NMFfitX') )
		x <- getRNG1(x)
	if( is(y, 'NMFfitX') )
		y <- getRNG1(y)
	
	identical(RNGdigest(x), RNGdigest(y))
}

# TODO: remove this function?
#rstream.rand <- function(n=1, provider){
#	
#	if( !missing(provider) ){
#		oldProv <- rngtools.set.provider(provider)
#		on.exit({rngtools.set.provider(oldProv, verbose="Restoring old provider")})
#	}
#	runif(n)	
#}

###% Creates an rstream object for R Core Random Number Generators
###% 
###%  
newRNGbase <- function(seed=NULL, kind){
	
	# if rngtools is disabled then directly use rstream.runif
	if( !getOption('rngtools') ){
		if( missing(kind) )
			kind <- "current"
		return(new('rstream.runif', seed=seed, kind=kind))
	}
	
	res <- if( RNGprovider() != 'base' ){
				
				if( missing(kind) )
					kind <- 'default'
				
				# get the user supplied RNG and restore it on exit
				urng <- .getRNG()
				on.exit(.setRNG(urng))
				# set temporary fake default RNG
				#RNGrecovery()
				# create rstream.runif object
				new('rstream.runif', seed=seed, kind=kind)
								
			}else{
				# create rstream.runif object using current kind if not otherwise
				# specified.
				if( missing(kind) )
					kind <- "current"				
				new('rstream.runif', seed=seed, kind=kind)
			}
	
	# add provider information
	rstream.user.provider(res) <- RNGprovider(user.supplied=TRUE)
	
	# return result
	res
	
}

###% Draw a Random Sample from the Default Random Number Generator
###% 
###% This function draws a random sample from R's default random number generator
###% for a given seed. It does not change the state of the active RNG.
###%  
rbase <- function(n, ...){
	
	# if the active RNG is a core-RNG then generate the seed using the 
	# current RNG (without changing its state), to get the same result
	# as when set.seed is called before calling nmf
	rng <- newRNGbase(...)
	r(rng, n)
	
}

###% This function finds out what rstream class to instantiate from the 
###% current RNG settings and set it as the currently stored RNG in the variable
###% .rstream.current in rstream global environment.
###%  
###% 
newRNG <- function(){
	
	# get current kind	
	ckind <- baseRNGkind()[1]
	
	# Core RNG: return a rstream.runif object based on the current settings
	if( ckind != "user-supplied" ){
		newRNGbase(NULL)
	}else{# User-supplied: finds out from the provider		
	
		# get current provider (for user-supplied RNG)
		provider <- RNGprovider()
		
		rclass <- paste('rstream', provider, sep='.')
		if( provider == 'rstream' )
			new('rstream.mrg32k3a')
		else if( extends(rclass, 'rstream2') )
			new(rclass)
		else if( extends(rclass, 'rstream') ){
			res <- new(rclass)
			rstream.user.provider(res) <- RNGprovider(user.supplied=TRUE)
			res
		}else
			new('rstream.user', provider=provider)
		
	}
	
	
}

# Get the RNG Settings 
setGeneric('getRNG', function(object, ...) standardGeneric('getRNG') ) 
setMethod('getRNG', 'missing',
	function(object, packed=FALSE){
		rng <- if( getOption('rngtools') ) rstream.GetRNGstate(.getRNG())
			   else rstream.RNG()

		# pack the rstream object if requested: this is a fresh clone so there is
		# no issue with freeing allocated ressources
		if( packed )
			rstream.packed(rng) <- TRUE
		
		rng
	}
)
setMethod('getRNG', 'ANY',
	function(object){
		if( .hasSlot(object, 'rng') ) slot(object, 'rng')
		else if( .hasSlot(object, 'rng.seed') ) slot(object, 'rng.seed') # for back compatibility
		else attr(object, 'rng')
	}
)
setMethod('getRNG', 'list',
	function(object){
		if( !is.null(object$rng) ) object$rng  
		else if( is.list(object$noise) ) object$noise$rng
		else attr(object, 'rng')
	}
)
setMethod('getRNG', 'rstream',
	function(object){
		object	
	}
)

.getRNG <- function(fixup=TRUE){
	
	# local pointer to the rstream environment
	.rstream.envir <- rstream:::.rstream.envir
	
	# possibly fixup the stored current RNG
	if( fixup )
		fixupRNG()
	
	## make a new object of the current Rstream object 		
	get(".rstream.current", envir=.rstream.envir)
	
}



###% Generic function that sets the Random Number Generator
###% 
###% 
setGeneric('setRNG', function(object, ...) standardGeneric('setRNG') )
setMethod('setRNG', 'character',
	function(object, ..., verbose=FALSE){
	
		# correct wrong input
		if( length(object) == 0 || !nchar(object) )
			object <- ''

		types <- c('rstream', 'lecuyer', 'base')
		s <- switch(object
				, rstream = new('rstream.mrg32k3a', ...)
				, lecuyer = new('rstream.lecuyer', ...)
				, base = newRNGbase(...)
				, stop("rngtools::setRNG - Invalid RNG type '", object, "'. Should be one of: ", paste(types, collapse=', '), ".")
			)
			
		setRNG(s, verbose=verbose)
	}
)
###% Set the current RNG with an rstream object
###% 
###% This function substitutes the function rstream.RNG (when called with a stream
###% argument). It uses the function getRNG instead of rstream.RNG (with no argument) 
###% to handle the case of user-supplied RNGs, and be able to return the old RNG 
###% in this case.
###% 
setMethod('setRNG', 'rstream',
	function(object, ...){
			
		# unpack and repack if necessary
		was.packed <- rstream.packed(object)
		# unpack if necessary and re-pack on exit
		if( was.packed )
			rstream.packed(object) <- FALSE		
		on.exit({
			if( was.packed )
				rstream.packed(object) <- TRUE
		})

		res <- if( getOption('rngtools') ){			
			## install a clone of the stream object
			.setRNG(rstream.clone(object), ...)
		}else{
			# set RNG directly using rstream interface: this installs a clone 
		 	# of object
			rstream.RNG(object)			
		}
		
		# return RNG object
		invisible(res)
	}
)


.RNGgenSeed <- function(seed){
	
	ru <- if( !missing(seed) ){
				
				# only use first element
				seed <- seed[1]
				rng <- new('rstream.runif', kind='default', seed=seed)
				r(rng, 6)
		}else
				runif(6)
	
	# between 0 and 999999
	ceiling(ru * 999999)
}

# set the seed of rstream 
RNGseed <- function(seed, verbose=FALSE){
	
	# retrieve current seed
	oldseed <- .rstream.get.seed()
	# return current value if missing seed
	if( missing(seed) ) return(oldseed)
	# generate seed if necessary
	if( is.null(seed) ){
		if( verbose ) message("# Generate RNGstream random seed ... ", appendLF=FALSE)
		seed <- .RNGgenSeed()
		if( verbose ) message("OK")
	}else if( is.numeric(seed) ){
		if( length(seed) == 1 ){
			if( verbose ) message("# Generate RNGstream random seed from ", seed, " ... ", appendLF=FALSE)
			seed <- .RNGgenSeed(seed)
			if( verbose ) message("OK")
		}
		else if( length(seed) != 6 )
			stop("RNGseed - Invalid numeric seed: should be a numeric of length 1 or 6")		
	}else if( !is_NA(seed) )
		stop("RNGseed - Invalid seed value: should be a single numeric, NULL or NA")
	
	if( verbose ) message("# Setting RNGstream random seed to: ", paste(seed, collapse=', '), " ... ", appendLF=FALSE)
	.rstream.set.seed(seed)
	oldseed

}

setMethod('setRNG', 'numeric',
	function(object, test=FALSE, verbose=FALSE, ...){
		
		# convert into an integer
		object <- as.integer(object)		
		
		rng <- if( length(object) == 1 ){
			# need special care if a user-supplied RNG is in use:
			# rstream.runif calls RNGkind, which draw once from the active RNG
			newRNGbase(seed=object)
			
		}else if( length(object) == 6 ){
			if( test ) RNGseq(1, seed, packed=FALSE)
			else new('rstream.mrg32k3a', seed=object, force=TRUE)
		}
		else stop("NMF::setRNG - invalid seed length ", length(object)," [expected 1 or 6]")
		
		# if testing then return the current seed
		if( test )
			return(rng)
		
		# set the new RNG, storing the old one
		old <- getRNG()
		setRNG(rng, verbose=verbose)
		# return old RNG as invisible
		invisible(old)
	}
)
setMethod('setRNG', 'ANY',
	function(object, ...){
		rng <- getRNG(object)
		if( is.null(rng) )
			stop("setRNG - could not extract RNG settings from object [class:", class(object), "]")
		setRNG(rng, ...)
	}
)

.setRNG <- function(object, verbose=FALSE, ...){
		
	stream <- object
	## input must be Rstream object
	if (!is(stream,"rstream")) stop ("Invalid argument `object`: must be an 'rstream' object")
	
	## store a clone of the current (future old) RNG as an rstream object to return
	current <- getRNG()
			
	# setup restoration in case of errors
	internals <- RNGinternals()
	on.exit({ 
			if( verbose ) message("# Restoring RNG internal settings")
			RNGinternals(internals, verbose=verbose)
		})
	
	## register the next current RNG provider to use:
	uprovider <- rstream.user.provider(stream)
	#rngtools.set.nextprovider( uprovider )

	# Wrap user-supplied RNG library: this will wrap only if really necessary
	if( verbose ) message("## Setting up user-supplied RNG library '", uprovider, "'")	
	wrapOK <- RNGwrap(uprovider, verbose, warn=FALSE)
	if( verbose ) message("## DONE ##")
	
	## set Rstream object as R generator
	if( rstream.provider(stream) != 'base' ){
		#wrapOK <- FALSE
		if( !wrapOK ){
			h <- RNGhooks(who=TRUE)
			h <- h[h!='']
			stop("NMF::setRNG - Could not setup provider '", uprovider, "' correctly, using RNG functions may give incorrect results or even crash the R session."
				, "\n\tRNG hook(s) provided: ", paste(h, '::', names(h), sep='', collapse=", "))
		}
		
		# for user-suuplied RNG: change the RNG kind to user-supplied, 
		# taking care of not using the current RNG if it is user-supplied itself
		# because RNGkind would draw once from it and it is probably incorrect 
		# anyway since one is already pointing to a possibly different RNG provider
		# If the current is a base RNG, drawing once from it does not have 
		# any side effect other than changing .Random.seed which will be overwritten
		# in any case.
		if( RNGprovider() != 'base' )
			RNGrecovery()
		
		# change kind to user-supplied: this will call the correct user 
		# initialisation hooks (user_unif_init, user_unif_nseed, user_unif_seedloc) 
		# if defined in the provider library
		RNGkind("user-supplied")		
	}
	
	# extract and set the state from the stream  
	stream <- rstream.PutRNGState(stream, ...)
	
#	## force the current provider for base RNGs
#	if( rstream.provider(stream) == 'base' ){
#		rngtools.set.provider(uprovider)
#		rngtools.set.nextprovider('')
#	}
	
	## check that the RNG providers are correct and consistent	
	stopifnot( rstream.user.provider(stream) == uprovider )	
	stopifnot( RNGprovider() == rstream.provider(stream) )
	stopifnot( RNGprovider(user.supplied=TRUE) == uprovider )
	stopifnot( rstream.get.nextprovider() == '' )
	
	## update current Rstream object
	assign(".rstream.current", stream, envir=rstream:::.rstream.envir)
	
	# cancel the restoration commands
	on.exit({})
	
	## return clone of old Rstream object
	invisible(current)	
}


#############################################################################
##                                                                         ##
##   Class: rstream2                                                       ##
##                                                                         ##
#############################################################################
##                                                                         ##
##   Interface to substitute the original class rstream and allow handling ## 
##   multiple user-supplied RNGs.                                          ##
##                                                                         ##
#############################################################################


## Initialize global variables ----------------------------------------------

.rstream2.init <- function () {}
## nothing to do


# Class ---------------------------------------------------------------------
setClass( "rstream2", 
		## additional slots	
		representation(				
			count = 'integer' # counter value from the rstream internal counter 
			, provider = "character" # package providing the user-supplied RNG
		),
		## extends rstream virtual class
		contains = c("rstream", "VIRTUAL") )

setMethod( "initialize", "rstream2",  
	function(.Object, ..., provider=NULL) {
		
		# local pointer to rstream environment
		.rstream.envir <- rstream:::.rstream.envir
		
		#str(.Object)
		# initialise parent slots: named arguments corresponding to slots will
		# overwrite slots in .Object
		.Object <- callNextMethod(.Object, ...)
		#str(.Object)
		
		## check input			
		if( is.null(provider) )
			provider <- RNGprovider(user.supplied=TRUE)
		
		.Object@provider <- provider
		
		## first increment counter for Rstreams
		.Object@count <- as.integer(get(".rstream.Count", envir=.rstream.envir) + 1)
		assign(".rstream.Count", .Object@count, envir=.rstream.envir)		
		
		## Create an environment for storing state of the generator
		.Object@xstate <- new.env(parent=.rstream.envir)
		
		## name of the Rstream object.
		#TODO: could the name be handled in this class?
		# this might be tricky for RNGs that rely on the name, although on could
		# imagine using an internal key (md5-like) that is actually used to 
		# identify the RNG on the C side, and allwo for a more user-friendly name
		# on the R side
		
		## type of Rstream object
		.Object@type <- "rstream2"
		
		## add info about Rstream type
		.Object@info <- "Generic class for uniform random number generator (v2)"
		
		## at creation a Rstream is never packed
		.Object@is.packed <- FALSE
		
		## return new Rstream object
		.Object
	} )

## rstream.count
##    get counter value for an rstream object
setGeneric('rstream.count', function(stream, ...) standardGeneric('rstream.count') )
###% Default rstream.count throws an error
setMethod("rstream.count", "rstream", 
	function(stream, ...){ 
		stop("rstream.count(",class(stream)[1],") - Unimplemented method: should have been overloaded") 
	}
)

###% Default rstream.provider throws an error
setMethod("rstream.count", "rstream2", 
	function(stream, ...){ 
		slot(stream, 'count') 
	}
)

## rstream.provider
##    get the provider of the RNG C level functions
setGeneric('rstream.provider', function(stream, ...) standardGeneric('rstream.provider') )

###% Default rstream.provider returns attribute 'provider' if present or ''
setMethod("rstream.provider", "rstream", 
	function(stream){
		p <- attr(stream, 'provider')
		if( !is.null(p) ) p else ''
	}
)

###% Method rstream.provider for rstream.provider objects: return the RNG provider
###% as stored in slot \code{provider}.
###% 
setMethod("rstream.provider", "rstream2", 
	function(stream, ...){
		slot(stream, 'provider')
	}
)

# implement method `rstream.provider` for known rstream classes
setMethod("rstream.provider", "rstream.runif", function(stream, ...){ 'base' })
setMethod("rstream.provider", "rstream.mrg32k3a", function(stream, ...){ 'rstream' })
setMethod("rstream.provider", "rstream.lecuyer", function(stream, ...){ 'rstream' })

setGeneric('rstream.provider<-', function(stream, ..., value) standardGeneric('rstream.provider<-') )
setReplaceMethod("rstream.provider", "rstream", 
	function(stream, ..., value){
		attr(stream, 'provider') <- value
		stream
	}
)
setReplaceMethod("rstream.provider", "rstream2", 
	function(stream, ..., value){
		slot(stream, 'provider') <- value
		stream
	}
)

###% Returns the provider for the user-supplied RNG of an rstream object. 
###% 
###% This is generally the same as the object's provider returned by 
###% \code{\link{rstream.provider}}, except for objects from class 
###% \code{\linkS4class{rstream.runif}}, for which it is the provider that would 
###% be used if the RNG kind was set to "user-supplied" with \code{\link{RNGkind}}.
###% 
setGeneric('rstream.user.provider', function(stream, ...) standardGeneric('rstream.user.provider') )
setMethod("rstream.user.provider", "rstream", 
	function(stream){
		rstream.provider(stream)
	}
)

###% Specific method for rstream.runif in order to return the attribute 'provider'
###% as the method rstream.provider(rstream.runif) always returns 'base'
###% 
setMethod("rstream.user.provider", "rstream.runif", 
	function(stream){
		p <- attr(stream, 'provider')
		if( !is.null(p) ) p else ''
	}
)
setGeneric('rstream.user.provider<-', function(stream, ..., value) standardGeneric('rstream.user.provider<-') )
setReplaceMethod("rstream.user.provider", "rstream", 
	function(stream, ..., value){
		rstream.provider(stream) <- value
		stream
	}
)

#############################################################################
##                                                                         ##
##   Class: rstream.user                                                   ##
##                                                                         ##
#############################################################################
##                                                                         ##
##   Interface to user-supplied uniform random number generators           ##
##                                                                         ##
#############################################################################



## Initialize global variables ----------------------------------------------

.rstream.user.init <- function () {}
## nothing to do


# Class ---------------------------------------------------------------------

setClass( "rstream.user", 
		## additional slots	
		representation(
				kind	= "character",	  # kind of RNG
				seed	= "integer" ),	  # seed of generator
		## extends rstream virtual class
		contains = "rstream2" )

## Initialize ---------------------------------------------------------------

setMethod( "initialize", "rstream.user",  
		function(.Object, ..., name=NULL, seed=NULL, antithetic=FALSE) {
						
			# local pointer to rstream environment
			.rstream.envir <- rstream:::.rstream.envir
			
			#str(.Object)
			# initialise parent slots: named arguments corresponding to slots will
			# overwrite slots in .Object
			.Object <- callNextMethod(.Object, ...)
			#str(.Object)
			
			if( rstream.provider(.Object) == 'base' )
				stop("Cannot create an rstream.user object for provider 'base': use class rstream.runif instead")
			
			## name of the Rstream object.
			## by default we use provider + number
			if (is.null(name)) name <- paste("user.", rstream.provider(.Object), rstream.count(.Object), sep="")
			assign("name", as.character(name), envir=.Object@xstate)
			
			## overwrite type of Rstream object
			.Object@type <- "user"
			
			## overwrite info about Rstream type
			.Object@info <- "Generic class for user-supplied uniform random number generator"
			
			## set antithetic flag
			assign("anti", as.logical(antithetic), envir=.Object@xstate)
			
			# use RNGinternals to directly get/restore all internal cached objects
			## save current state of runif.
			## however, .Random.seed might not exist.
			## then we run sample(NA) one time to initialize it.
			if (!exists(".Random.seed", envir=.GlobalEnv)) 
				sample(NA)
						
			# setup restoration of all RNG settings
			internals <- RNGinternals()
			on.exit({ RNGinternals(internals, TRUE) })
			
			## register the next current RNG provider to use:
			rngtools.set.nextprovider(.Object)
			
			## Set user-supplied RNG: this will initialize the user-suplied RNG
			# with a random seed (if the provider has hooks user_unif_init, etc...)			
			if( rstream.is.restorable() ){
				RNGkind("user-supplied")
				stopifnot( RNGprovider() == rstream.provider(.Object) )
				
				## set the seed if one is provided by user
				if(!is.null(seed))
					set.seed(seed)
				
				## store the new .Random.seed
				.Object@seed <- get(".Random.seed", envir=.GlobalEnv)
				
			}else{
				warning("Creating an rstream.user the object: current RNG settings do not allow restoration.\n"
						, "The object's seed and state were left uninitialized.")
			}
			
			## store kind of stream
			.Object@kind <- "user-supplied"
			
			## store state of generator in its own environment
			assign("state", .Object@seed, envir=.Object@xstate)
			
			## return new Rstream object
			.Object
		} )


## Validity -----------------------------------------------------------------


## Methods ------------------------------------------------------------------

## Access and Replacement methods ............................................

## rstream.name
##    get and set name of Rstream object
setMethod("rstream.name", "rstream.user", 
		function(stream) { 
			if (stream@is.packed) 
				return (stream@pack$name)
			else
				return (get("name", envir=stream@xstate))
		} )

setReplaceMethod("rstream.name", "rstream.user", 
		function(stream, value) {
			if (stream@is.packed) stop("Cannot change name for PACKED Rstream") 
			assign("name", as.character(value), envir=stream@xstate)
			stream
		} )


## rstream.antithetic
##   get and set flag for antithetic random numbers:  
setMethod("rstream.antithetic", "rstream.user", 
		function(stream) { 
			if (stream@is.packed) 
				return (stream@pack$anti)
			else
				return (get("anti", envir=stream@xstate))
		} )

setReplaceMethod("rstream.antithetic", "rstream.user",
		function(stream, value) { 
			if (stream@is.packed) stop("Cannot change antithetic flag for PACKED Rstream") 
			assign("anti", as.logical(value), envir=stream@xstate)
			stream
		} )


## Sampling methods .........................................................

## rstream.sample
##    make a random sample
setMethod("rstream.sample", "rstream.user",
		function(stream,n=1) .rstream.provider.sample(stream,n) )

setMethod("r", "rstream.user",
		function(stream,n=1) .rstream.provider.sample(stream,n) )

.rstream.provider.sample <- 
		function(stream,n=1) { 
	if (stream@is.packed) stop("Cannot sample from PACKED Rstream") 
	
	## store current state of R generator for restoration on.exit
	oldRNG <- getRNG()
	on.exit({.setRNG(oldRNG)}, add=TRUE)
	
	## temporary set the RNG
	# NB: a clone is actually installed
	.setRNG(stream)
	
	## sample	
	x <- runif(n)
	
	## the original RNG's state should have changed
	
	## return result (depending on the antithetic variable flag)
	if ( rstream.antithetic(stream) )
		return (1-x)
	else
		return (x)
}


## Reset and copy methods ...................................................

## rstream.reset
##   reset Rstream object
setMethod("rstream.reset", "rstream.user", 
		function(stream) { 
			if (stream@is.packed) stop("Cannot reset PACKED Rstream") 
			## copy seed into state variable
			assign("state", stream@seed, envir=stream@xstate) } )


## rstream.clone
##    clone (copy) Rstream object
setMethod("rstream.clone", "rstream.user", 
		function(stream) { 
			if (stream@is.packed) stop("Cannot clone PACKED Rstream") 
			## copy R object
			clone <- stream
			## copy data from old to a new environment
			clone@xstate <- copyEnv(stream@xstate)
			## add a period '.' to name to mark this as clone
			rstream.name(clone) <- paste(rstream.name(clone), ".", sep="")
			## clone the external pointer
			clone@stream <- rstream.cloneptr(stream, clone)
			
			## return clone of Rstream object
			clone
		} )


## rstream.pack, rstream.unpack
##    pack and unpack Rstream object such that all data are contained in R object
##    (and can be easily copied within R)
setReplaceMethod("rstream.packed", "rstream.user", 
		function(stream, value) {
			value <- as.logical(value)
			## check whether there is something to do
			if (value && stream@is.packed)   return(stream)
			if (!value && !stream@is.packed) return(stream)
			
			if (value) {
				## pack
				stream@is.packed <- TRUE
				stream@pack <- as.list(stream@xstate)
			}
			else {
				## unpack
				stream@is.packed <- FALSE
				## copy data into new environment
				stream@xstate <- list2env(stream@pack, parent=rstream:::.rstream.envir)				
			}
			## return un/packed Rstream object
			stream
		} )


## Printing et al. ..........................................................

## print:
##    print all data of a Rstream object
setMethod( "print", "rstream.user",
		function(x, ...) {
			rstream:::.rstream.PrintData(x) 
			cat(	"\tRNGkind = ", x@kind, "\n", sep="")
			cat(	"\tProvider = ", x@provider, "\n", sep="")
		} )

## Rstream objects <-> R generators -----------------------------------------

# USE DEFAULTS FROM rstream

## Cloning of the internal pointer

###% Generic function that clones the external pointer from one rstream object to 
###% another
###%
###% 
setGeneric('rstream.cloneptr', function(stream, dest, ...) standardGeneric('rstream.cloneptr') )

###% Default method to clone the external pointer in rstream objects
###% 
###% Throws an error.
###% 
setMethod("rstream.cloneptr", "rstream", 
		function(stream, dest) { 
			stop("rstream.cloneptr(", class(stream)[1], ") - Unimplemented virtual method: should have been defined")
		}
)

###% Method to clone the external pointer in rstream.runif objects
###% 
###% Returns a Nil pointer.
###% 
setMethod("rstream.cloneptr", "rstream.user", 
	function(stream, dest) {
		if( ptr_isnil(stream@stream) )
			new('externalptr')
		else
			stop("rstream.cloneptr(", class(stream)[1], ") - Unimplemented virtual method: could not clone non nil external pointer")
	}
)

###% Method to clone the external pointer of rstream.mrg32k3a
###% 
###% It calls the C level method 'R_RngStreams_Clone' from the package rstream, 
###% which manages the cloning operation.
###% 
setMethod("rstream.cloneptr", "rstream.mrg32k3a", 
		function(stream, dest){
			if (stream@is.packed) stop("Cannot clone external pointer of PACKED Rstream") 
			.Call("R_RngStreams_Clone", clone, stream@stream, rstream.name(dest), PACKAGE="rstream")		
		} 
)

###% Method to clone the external pointer of rstream.lecuyer
###% 
###% It calls the C level method 'R_RngStreams_Clone' from the package rstream, 
###% which manages the cloning operation.
###% 
setMethod("rstream.cloneptr", "rstream.lecuyer", 
		function(stream, dest){
			if (stream@is.packed) stop("Cannot clone external pointer of PACKED Rstream") 
			.Call("R_RngStreams_Clone", clone, stream@stream, rstream.name(dest), PACKAGE="rstream")		
		} 
)

###% Method for loading the state of the current RNG into an rstream.user object
###% 
###% This method should work for all rstream classes that store the state 
###% information into an external pointer, that is dynamically updated on the C 
###% side, and copied by the method rstream.cloneptr 
###% 
setMethod("rstream.GetRNGstate", "rstream.user", 
	function(stream, ...){
		stream <- rstream.clone(stream)
		
		# Update environment data from .Random.seed
		rseed <- get(".Random.seed", envir=.GlobalEnv)
		state <- assign("state", rseed, envir=stream@xstate)		
		
		stream
	}
)

###% Method for setting an rstream.provider object as the RNG
###% 
###% This method should work for any rstream class that implements a user-supplied 
###% RNG and supplies the C level functions user_unif_rand, user_unif_seedns 
###% and user_unif_seedloc.
###%  
###% 
setMethod("rstream.PutRNGState", "rstream.user", 
	function(stream, allow.anti=FALSE, ...){
		
		if ( !allow.anti && rstream.antithetic(stream)) {
			warning ("rstream::PutRNGState(",class(stream)[1],") - antithetic DISABLED")
			rstream.antithetic(stream) <- FALSE 
		}
		
		# do NOT use RNGkind as it would draw once from the current RNG, 
		# it has been already called before rstream.PutRNGState		
		state <- get("state", envir=stream@xstate)
		assign(".Random.seed", state, envir=.GlobalEnv)
		stream		  
	}
) 


## End ----------------------------------------------------------------------
