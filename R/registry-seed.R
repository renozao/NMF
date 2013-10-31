# Registry for NMF seeding method
# 
# Author: Renaud Gaujoux
###############################################################################

#' @include registry.R
#' @include NMFSeed-class.R
NULL

# create sub-registry for seeding methods
.registrySeed <- setPackageRegistry('seed', "NMFSeed"
		, description = "Initialization methods for NMF algorithms"
		, entrydesc = 'NMF seeding method')

nmfSeedInfo <- function(show=TRUE){
    obj <- .registrySeed
    if( show ) print(obj)
    invisible(obj)
}

#' Seeding Strategies for NMF Algorithms
#' 
#' \code{nmfSeed} lists and retrieves NMF seeding methods.
#' 
#' Currently the internal registry contains the following seeding methods, 
#' which may be specified to the function \code{\link{nmf}} via its argument 
#' \code{seed} using their access keys:
#' 
#' \describe{
#' \item{random}{ The entries of each factors are drawn from a uniform 
#' distribution over \eqn{[0, max(x)]}, where $x$ is the target matrix.}
#' \item{nndsvd}{ Nonnegative Double Singular Value Decomposition.
#' 
#' The basic algorithm contains no randomization and is based on two SVD processes, 
#' one approximating the data matrix, the other approximating positive sections 
#' of the resulting partial SVD factors utilising an algebraic property of 
#' unit rank matrices.
#' 
#' It is well suited to initialise NMF algorithms with sparse factors.
#' Simple practical variants of the algorithm allows to generate dense factors.
#' 
#' \strong{Reference:} \cite{Boutsidis2008}}
#' \item{ica}{ Uses the result of an Independent Component Analysis (ICA) 
#' (from the \code{fastICA} package).
#' Only the positive part of the result are used to initialise the factors.}
#' \item{none}{ Fixed seed.
#' 
#' This method allows the user to manually provide initial values for 
#' both matrix factors.}
#' }
#' 
#' @param name access key of a seeding method stored in registry.
#' If missing, \code{nmfSeed} returns the list of all available seeding methods.
#' @param ... extra arguments used for internal calls
#'  
#' @export
#' 
#' @examples
#' 
#' # list all registered seeding methods
#' nmfSeed()
#' # retrieve one of the methods
#' nmfSeed('ica') 
#' 
nmfSeed <- function(name=NULL, ...){
	
	nmfGet('seed', name, ...)
	
}

#' \code{getNMFSeed} is an alias for \code{nmfSeed}.
#' @rdname nmfSeed
#' @export
getNMFSeed <- nmfSeed

#' \code{existsNMFSeed} tells if a given seeding method exists in the registry.
#' 
#' @param exact a logical that indicates if the access key should be matched 
#' exactly or partially.
#'  
#' @rdname nmfSeed
#' @export
existsNMFSeed <- function(name, exact=TRUE){	
	
	res <- !is.null( getNMFSeed(name, error=FALSE, exact=exact) )
	return(res)
	
}

# specific register method for registering NMFSeed objects
setMethod('nmfRegister', signature(key='NMFSeed', method='missing'), 
		function(key, method, ...){
			nmfRegister(name(key), key, ..., regname='seed')
		}
)

#' Registering NMF Seeding Methods
#' 
#' NMF seeding methods are registered via the function \code{setNMFSeed}, which
#' stores them as \code{\linkS4class{NMFSeed}} objects in a dedicated registry.
#' 
#' @param ... arguments passed to \code{NMFSeed} and used to initialise slots
#' in the \code{\linkS4class{NMFSeed}} object, or to \code{\link[pkgmaker]{pkgreg_remove}}.
#' @inheritParams setNMFMethod
#' 
#' @export
setNMFSeed <- function(..., overwrite=isLoadingNamespace(), verbose=TRUE){
	
	# wrap function method into a new NMFSeed object
	method <- NMFSeed(...)
	# register the newly created object
	res <- nmfRegister(method, overwrite=overwrite, verbose=verbose)	
}

nmfRegisterSeed <- setNMFSeed


#' \code{removeNMFSeed} removes an NMF seeding method from the registry.
#' 
#' @param name name of the seeding method.
#' 
#' @export
#' @rdname setNMFSeed
removeNMFSeed <- function(name, ...){
	pkgreg_remove('seed', key=name, ...)
}

