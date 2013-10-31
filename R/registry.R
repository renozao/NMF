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

# Returns the names of all the packages that contibute to all or a given
# package's primary registry  
registryContributors <- function(package, regname = NULL){
    regs <- packageRegistries(regname = regname, package = package, primary = TRUE)
    if( length(regs) ) unique(names(unlist(lapply(paste0(package, '::', regs), packageRegistries))))
}

###% Return a method stored in the NMF registry.
###% 
###% @param name the key (a character string) of the method to be retrieved
###% @param regname the name of the sub-registry where to look for the \code{key}
###% @param exact a boolean. When set to \code{TRUE} the key is searched exactly, otherwise (default) the key
###% is matched partially with the keys registered in the registry.
###% @param error a boolean. When set to \code{TRUE} (default) the function will raise an error if the key is not found.
###% Otherwise it will not raise any error and return \code{NULL}.
###%
nmfGet <- function(regname, name=NULL, ...){
	
	# retrieve from the given package's sub-registry
	pkgmaker::pkgreg_fetch(regname, key=name, ...)
	
}

###% Register a NMF method so that it is accessible via the common interface defined by the \code{nmf} function.
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
	function(key, method, regname, ...){		
		#TODO: add functionality to save the registered strategy into a file for use is other R sessions
		
		parent.method <- attr(method, 'parent')
		tmpl <- if( !is.null(parent.method) && parent.method != key ){
			str_c(" based on template '", parent.method, "'")
		}
		setPackageRegistryEntry(regname, key, method, ..., where='NMF', msg=tmpl)
	}
)

####% Unregister a NMF method.
####%
####% @param name the key of the method to unregister [character string]
####%
#nmfUnregister <- function(name, regname, quiet=FALSE){				
#	
#	return( pkgreg_remove(regname, key=name, quiet=quiet) )
#	# add the strategy to the registry
#	obj <- nmfGet(name, exact=TRUE, error=FALSE, regname=regname)
#	regentry <- nmfRegistry(regname, entry=TRUE)
#	registry <- regentry$regobj
#	objtype <- regentry$entrydesc
#	
#	if( !is.null(obj) ){
#		# get the method registry and the method's fullname
#		name <- attr(strategy, 'name')
#		
#		if( !quiet ){
#			msg <- paste0("Removing ", objtype," '", name, "' from registry '", regname, "' [", class(obj), ']')
#			message(msg, ' ... ', appendLF=FALSE)
#		}
#		# delete from registry
#		registry$delete_entry(name)
#		if( !quiet ) message('OK')
#		TRUE
#	}else{
#		if( !quiet )
#			warning("Could not remove ", objtype, " '", name, "': no matching registry entry.", call.=FALSE)
#		FALSE
#	}
#}
