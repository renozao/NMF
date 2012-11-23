###% Define access/setup methods for NMF package registry.
###% 
###% The registry is used to provide a common interface to NMF methods (algorithms, seeding methods, distance, ...).
###% It enables the user to add custom methods which will be accessible in the same way as the built-in ones.
###% 
###% @author Renaud Gaujoux
###% @created 22 Jul 2009 


###########################################################################
# COMMON REGISTRY
###########################################################################

#' @import pkgmaker
#' @import registry
nmfRegistry <- function(...) pkgmaker::packageRegistry(...)

.init.nmf.registry <- function(){
	
#	# initialize main repository
#	registryName <- '.NMFRegistry'
#	registryEnv <- parent.env(environment())
#	
#	msg.addon <- if( exists(registryName) ) ' [reset]' else NULL
#	sreg <- utils::capture.output(print(registryEnv))
#	message("Create NMF main registry in ", sreg, msg.addon)		
#	.NMFRegistry <<- new.env(.NMFRegistry)
#	
#	# initialize all sub-registries
#	sapply(nmfSubRegistry(), function(n) nmfSubRegistry(n, create=TRUE))
	
	return(invisible(TRUE))
}

###% Return a method stored in the NMF registry.
###% 
###% @param name the key (a character string) of the method to be retrieved
###% @param registry.name the name of the sub-registry where to look for the \code{key}
###% @param exact a boolean. When set to \code{TRUE} the key is searched exactly, otherwise (default) the key
###% is matched partially with the keys registered in the registry.
###% @param error a boolean. When set to \code{TRUE} (default) the function will raise an error if the key is not found.
###% Otherwise it will not raise any error and return \code{NULL}.
###%
nmfGet <- function(registry.name, name=NULL, ...){
	
	# retrieve the required sub-registry
	registry <- nmfRegistry(registry.name)
	pkgmaker::regfetch(registry, key=name, ...)
	
}

###% Register a NMF method so that it is accessible via the common interface defined by the \code{nmf} function.
###% TODO: rewrite the doc (obsolete)
###% @param method an NMFStrategy object or a function that defines the method
###% @param key a non-empty character string that will be used as an identifier to access the method
###% @param overwrite a boolean that specify if an existing method (i.e. with exactly the same \code{key}) should be overwritten or not.
###% If \code{FALSE} and a method with the same key exists, an error will be thrown.
###% @param save [Not used] a boolean that if set to \code{TRUE} will save in database so that it is available in other R sessions.
###% @param ... [Not used]
###%
###% @return \code{TRUE} invisibly in case of success.
###%
###% @seealso nmf
###%
setGeneric('nmfRegister', function(method, key, ...) standardGeneric('nmfRegister') )
setMethod('nmfRegister', signature(method='ANY', key='character'), 
		function(method, key, registry.name, overwrite=FALSE, verbose=FALSE){		
			#TODO: add functionality to save the registered strategy into a file for use is other R sessions
						
			# check if the name provided is not empty
			if( nchar(key) == 0 ) stop('parameter <key> cannot be an empty string')
			
			# for a method given as a NMFStrategy, if the name is empty then use the key as a name
			if( inherits(method, 'NMFStrategy') ){
				if( is.null(method@name) || method@name=='' ){
					stop("Unexpected error: nmfRegister(method=NMFStrategy) with no method name")
					method@name <- key
				}
				#else if( method@name != key ) stop("NMFStrategy slot <name> does not match parameter <key> ['", method@name, "' != '", key, "']")
			}
			# check if the object is already registered
			reg.method <- nmfGet(key, registry.name=registry.name, exact=TRUE, error=FALSE)
			
			# add the strategy to the registry		
			if( is.null(reg.method) || overwrite ){ # add or overwrite if necessary
				# retrieve the NMFStrategy registry
				registry <- nmfRegistry(registry.name)
				# add the method to the registry
				action <- if( is.null(reg.method) ) 'Register' 
				else{
					registry$delete_entry(key)
					'Overwrite'
				}
				if( verbose ) message("NMF:", registry.name, ": ", action, " method '", key, "'")
				registry$set_entry(key=key, object=method)
				#assign(key, method, envir=registry)
			}
			else stop("Cannot register NMF method '", key, "' [a method with the same key is already registered].")
			
			# setup deferred registering: when loading a namespace other than NMF
			packageNMFObject(key, method)
			#
			
			# return set/update code invisibly
			invisible( if( !is.null(reg.method) )  2L else 1L )
		}
)

###% Unregister a NMF method.
###%
###% @param name the key of the method to unregister [character string]
###%
nmfUnregister <- function(name, registry.name){				
	
	# add the strategy to the registry
	strategy <- nmfGet(name, registry.name=registry.name, exact=TRUE, error=FALSE)
	if( !is.null(strategy) ){
		
		# get the method registry and the method's fullname
		registry <- nmfRegistry(registry.name)
		name <- attr(strategy, 'name')
		message("NMF: Remove method '", name, "' from registry '", registry.name, "'")
		registry$delete_entry(name)
		
		# cancel deferred registering: when loading a namespace other than NMF
		packageNMFObject(name, NULL)
		#
	}
	
	# return TRUE invisibly
	invisible(TRUE)
}
nmfUnregister <- Vectorize(nmfUnregister, 'name')

packageNMFObject <- function(key, method){

	library(pkgmaker)
	# do nothing if 
	# - not loading a namespace
	# - loading the NMF namespace
	if( !isLoadingNamespace(nodev=TRUE) || isLoadingNamespace('NMF') ){
		return()
	}
	# do nothing if the namespace is already completely loaded
	ns <- getLoadingNamespace(env=TRUE)
	ns_name <- getLoadingNamespace()
#	ca <- match.call()
#	message('Call: ', capture.output(print(ca)), ' - ', ns_name , ' - ', capture.output(ns)
#			, paste(pkgmaker:::pkg_calls(), collapse=', '))
	
	# auxiliary function to access the package's registered methods 
	pkgMethods <- simpleRegistry('.__NMFmethods__', ns)
	
	# retrieve method
	if( missing(method) ){
		return( pkgMethods$get(key) )
	}

	# build access key
	lkey <- str_c(packageSlot(method), '::', key)
	
	# removing method
	if( is.null(method) ){
		pkgMethods$set(lkey, NULL)
		# cancel deferred loading
		postponeAction(NULL, lkey, group='NMF')
		return()
	}
	
	# adding a new method
	#message("Adding local NMF method '", lkey,"' to package '", ns_name, "'")
	pkgMethods$set(lkey, method)
	# defer loading (will be executed in the NMF namespace)
	f <- function(...){
		setNMFObject(packageNMFObject(lkey), verbose=getOption('verbose'))
	}
	e <- topenv()
	postponeAction(f, lkey, envir=e, group='NMF')
	#
	
	return()
}

setNMFObject <- function(object, ...){
	
	# use this as common interface for registering NMF algorithms, seeding methods, etc..
	if( is(object, 'NMFSeed') ){
		return( setNMFSeed(object, ...) )
	}else if( is(object, 'NMFStrategy') ){
		return( setNMFMethod(object, ...) )
	}else
		stop("Unexpected error: objects of class '", class(object), "' are not supported.")
	
}
