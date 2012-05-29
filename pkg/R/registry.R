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
				if( is.null(method@name) || method@name=='' ) method@name <- key
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
		
	}
	
	# return TRUE invisibly
	invisible(TRUE)
}
nmfUnregister <- Vectorize(nmfUnregister, 'name')


#' Get Stopping Criterion
#' 
#' Convert a key into a stopping criterion (i.e. a \code{function}).
#' 
#' Functions are returned unchanged, integer values are used to create a stopping 
#' criterion of that number of iterations, numeric values are used to create a 
#' stopping criterion of that stationary threshold.
#' 
#' @param val access key that can be a character string, a single integer or 
#' numeric, or a function.
#' 
#' @return a function
#' @aliases stop-nmf stop-NMF
#' @export
NMFStop <- function(val){
	
	key <- val
	if( is.integer(key) )	nmf.stop.iteration(key)
	else if( is.numeric(key) ) nmf.stop.threshold(key)
	else if( is.function(key) ) key
	else if( is.character(key) ){
		# update .stop for back compatibility:
		if( key == 'nmf.stop.consensus') key <- 'connectivity'
		
		# first lookup for a `nmf.stop.*` function
		key2 <- paste('nmf.stop.', key, sep='')
		e <- pkgmaker::packageEnv()
		sfun <- getFunction(key2, mustFind=FALSE, where = e)
		if( is.null(sfun) ) # lookup for the function as such
			sfun <- getFunction(key, mustFind = FALSE, where = e)			
		if( is.null(sfun) )
			stop("Invalid key ['", key,"']: could not find functions '",key2, "' or '", key, "'")
		sfun
	}else if( identical(val, FALSE) ) # create a function that does not stop 
		function(strategy, i, target, data, ...){FALSE}
	else
		stop("Invalid key: should be a function, a character string or a single integer/numeric value. See ?NMFStop.")	
}
