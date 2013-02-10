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
	
	# retrieve from the given package's sub-registry
	pkgmaker::pkgregfetch(registry.name, key=name, ...)
	
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
setGeneric('nmfRegister', function(key, method, ...) standardGeneric('nmfRegister') )
setMethod('nmfRegister', signature(key='character'), 
	function(key, method, registry.name, ...){		
		#TODO: add functionality to save the registered strategy into a file for use is other R sessions
		
		parent.method <- attr(method, 'parent')
		tmpl <- if( !is.null(parent.method) && parent.method != key ){
			str_c(" based on template '", parent.method, "'")
		}
		pkgmaker:::setPackageRegistryEntry(key, registry.name, method, ..., where='NMF', msg=tmpl)
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
