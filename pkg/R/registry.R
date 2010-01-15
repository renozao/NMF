#' Define access/setup methods for NMF package registry.
#' 
#' The registry is used to provide a common interface to NMF methods (algorithms, seeding methods, distance, ...).
#' It enables the user to add custom methods which will be accessible in the same way as the built-in ones.
#' 
#' @author Renaud Gaujoux
#' @created 22 Jul 2009 


###########################################################################
# COMMON REGISTRY
###########################################################################

# initialze registry variable (IMPORTANT when loaded in a namespace)
.NMFRegistry <- NULL
.init.nmf.registry <- function(){
	
	# initialize main repository
	registryName <- '.NMFRegistry'
	registryEnv <- parent.env(environment())
	
	msg.addon <- if( exists(registryName) ) ' [reset]' else NULL
	sreg <- utils::capture.output(print(registryEnv))
	message("Create NMF main registry in ", sreg, msg.addon)		
	.NMFRegistry <<- new.env(.NMFRegistry)
	
	# initialize all sub-registries
	sapply(nmfSubRegistry(), function(n) nmfSubRegistry(n, create=TRUE))
	
	return(invisible(TRUE))
}

#' Return the required NMF (sub-)registry.
#'
#' The NMF registry stores built-in and custom NMF methods (algorithms, seeding methods, distance).
#' It provides
nmfRegistry <- function(name, reset=FALSE){
	
	# if no name is provided return the complete NMF registry
	if( missing(name) ){							
		# return the full registry
		return(.NMFRegistry)
	}
	
	# look for the subregistry entry
	sub.registry <- nmfSubRegistry(name)
	
	# return the subregistry
	return(sub.registry)	
}

nmfSubRegistry <- function(name, error=TRUE, create=FALSE){
	
	# local mapping of sub-registries' names
	.RegistryNames <- c(algorithm='.NMFStrategies', seed='.NMFSeed', distance='.NMFDistance')
	
	# return the mapping if no name was supplied
	if( missing(name) ) return(names(.RegistryNames))
	else if( !is.character(name) 
			|| length(name) != 1  
			|| nchar(name) == 0)
		stop('Invalid name for sub-registry: must be a single non-empty character string')
	
	# check if the sub-registry exists
	if( !is.element(name, names(.RegistryNames)) ){
		if( error ) stop("NMF registry: sub-registry '", name, "' does not exist")
		else return(invisible(NULL))
	}
	actual.name <- .RegistryNames[name]
	
	# check the registry has been initialized
	registryEnv <- nmfRegistry()
	if( !exists(actual.name, envir=registryEnv, inherits=FALSE) ){
		if( create ){
			message("Create NMF sub-registry '", name,"'")
			assign(.RegistryNames[name], new.env(registryEnv), envir=registryEnv)
		}
		else stop("NMF registry: sub-registry '", name, "' not initialized")
	}
	sub.registry <- get(actual.name, envir=registryEnv, inherits=FALSE)
	
	# return the sub-registry
	return(sub.registry)
}

#' Return a method stored in the NMF registry.
#' 
#' @param name the key (a character string) of the method to be retrieved
#' @param registry.name the name of the sub-registry where to look for the \code{key}
#' @param exact a boolean. When set to \code{TRUE} the key is searched exactly, otherwise (default) the key
#' is matched partially with the keys registered in the registry.
#' @param error a boolean. When set to \code{TRUE} (default) the function will raise an error if the key is not found.
#' Otherwise it will not raise any error and return \code{NULL}.
#'
nmfGet <- function(name, registry.name, exact=FALSE, error=TRUE){
	
	# force the user to provide the registry's name 
	if( missing(registry.name) ) stop("Parameter 'registry.name' is required")
	
	# retrieve the required registry	
	registry <- nmfRegistry(registry.name)
	
	# if no name is provided then return the complete registry
	if( missing(name) || is.null(name) ) return(ls(registry))
	
	# if the registry is empty: send an error/warning if not running in silent mode
	if( length(registry) == 0 ){
		msg <- paste("Could not find NMF method '", name, "' in registry '", registry.name, "' [sub-registry is empty]", sep='')
		if( error ) stop(msg)			
		return(NULL);
	}
	
	# resolve the strategy's fullname	
	fullname <- if( exact ) name else try( match.arg(name, ls(registry)), silent=TRUE)
	
	# if the strategy was not found throw an error or return NULL
	if ( inherits(fullname, 'try-error') ){
		msg <- paste("Could not find NMF method '", name, "' in registry '", registry.name, "'.\n  -> Registered methods: ", paste(paste("'", ls(registry), "'", sep=''), collapse=', '), '.', sep='')
		if( error ) stop(msg)
		if( getOption('verbose') ) warning(msg)
		return(NULL)
	}
	
	# add the name attribute if necessary (it is not for methods defined by an NMFStrategy object
	res <- registry[[fullname]]
	if( !is.null(res) && is.null(attr(res, 'name')) ) attr(res, 'name') <- fullname
	
	# return the NMF strategy
	res	
}

#' Register a NMF method so that it is accessible via the common interface defined by the \code{nmf} function.
#' TODO: rewrite the doc (obsolete)
#' @param method an NMFStrategy object or a function that defines the method
#' @param key a non-empty character string that will be used as an identifier to access the method
#' @param overwrite a boolean that specify if an existing method (i.e. with exactly the same \code{key}) should be overwritten or not.
#' If \code{FALSE} and a method with the same key exists, an error will be thrown.
#' @param save [Not used] a boolean that if set to \code{TRUE} will save in database so that it is available in other R sessions.
#' @param ... [Not used]
#'
#' @return \code{TRUE} invisibly in case of success.
#'
#' @seealso nmf
#'
if ( is.null(getGeneric('nmfRegister')) ) setGeneric('nmfRegister', function(method, key, ...) standardGeneric('nmfRegister') )
setMethod('nmfRegister', signature(method='ANY', key='character'), 
		function(method, key, registry.name, overwrite=FALSE, save=FALSE, ...){		
			#TODO: add functionnality to save the registered strategy into a file for use is other R sessions
			
			# check if the name provided is not empty
			if( nchar(key) == 0 ) stop('parameter <key> cannot be an empty string')
			
			# for a method given as a NMFStrategy, if the name is empty then use the key as a name
			if( inherits(method, 'NMFStrategy') ){
				if( is.null(method@name) || method@name=='' ) method@name <- key
				#else if( method@name != key ) stop("NMFStrategy slot <name> does not match parameter <key> ['", method@name, "' != '", key, "']")
			}
			# check if the object is already registered
			reg.method <- nmfGet(key, registry.name, exact=TRUE, error=FALSE)
			
			# add the strategy to the registry		
			if( is.null(reg.method) || overwrite ){ # add or overwrite if necessary
				# retrieve the NMFStrategy registry
				registry <- nmfRegistry(registry.name)
				# add the method to the registry
				action <- if( is.null(reg.method) ) 'Register' else 'Overwrite'
				message("NMF:", registry.name, ": ", action, " method '", key, "'")
				assign(key, method, envir=registry)
			}
			else stop("Cannot register NMF method '", key, "' [a method with the same key is already registered].")
			
			# return nothing invisibly
			invisible(TRUE)
		}
)

#' Unregister a NMF method.
#'
#' @param name the key of the method to unregister [character string]
#'
nmfUnregister <- function(name, registry.name){				
	
	# add the strategy to the registry
	strategy <- nmfGet(name, registry.name, exact=TRUE, error=FALSE)
	if( !is.null(strategy) ){
		
		# get the method registry and the method's fullname
		registry <- nmfRegistry(registry.name)
		name <- attr(strategy, 'name')
		message("NMF: Remove method '", name, "' from registry '", registry.name, "'")
		rm(list=c(name), envir=registry)
		
	}
	
	# return TRUE invisibly
	invisible(TRUE)
}
nmfUnregister <- Vectorize(nmfUnregister, 'name')
