# TODO: Add comment
# 
# Author: Renaud Gaujoux
# Creation: 08 Nov 2011
###############################################################################

#' RNG Settings
#' 
#' The functions documented here provide a unified interface to work with 
#' RNG settings.
#' 
#' @rdname rng
#' @name RNG-NMF
NULL

#' \code{RNGrecovery} allows to recover from a broken state of \code{.Random.seed}.
#' 
#' @export
#' @rdname rng
RNGrecovery <- function(){
	s <- as.integer(c(401,0,0))
	assign(".Random.seed", s, envir=.GlobalEnv)
	RNGkind("default")
}

###% Returns all the libraries that provides a user-supplied RNG
###% 
###% The library that provides the wrapper hooks for the management multiple 
###% user-supplied RNG is removed from the output list.
###% 
RNGlibs <- function(n=0, full=FALSE, hook="user_unif_rand", unlist=TRUE){
	dlls <- getLoadedDLLs()
	res <- lapply(dlls, function(d){
				dname <- d[['name']]
				if( dname=='' )
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
		libs <- RNGlibs(full=TRUE, unlist=FALSE, hook=hook)
		w <- which(sapply(libs, function(l) identical(l$address, symb.unif_rand$address)))
		
		# returns full info or just the name
		if( full ) libs[[w]]
		else names(libs)[w]
		
	}else if( full ) symb.unif_rand
	else symb.unif_rand$package[['name']]
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
	
	kind <- RNGkind()
	if( kind[1] == 'user-supplied' || user.supplied ) RNGlib()		
	else 'base'
}

RNGscope <- function(seed){
	
	res <- if( missing(seed) ){
				if( exists('.Random.seed', where = .GlobalEnv) )
					get('.Random.seed', .GlobalEnv)
			}else if( is.null(seed) ){
				if( exists('.Random.seed', where = .GlobalEnv) )
					rm('.Random.seed', envir = .GlobalEnv)
			}else{
				old <- RNGscope()
				assign('.Random.seed', seed, .GlobalEnv)
				old
			}
	invisible(res)
}

###% Returns a single string that describes the given RNG state
###% 
RNGdesc <- function(seed){
	
	if( missing(seed) ){
		rp <- RNGprovider()
		rs <- getRNG()
		if( rp == 'base' || length(rs) > 1L )
			seed <- rs
		else 
			return( "Unknown" )		
	}
	
	if( is.null(seed) ) 'NULL'
	else if( is.numeric(seed) ){
		if( length(seed) > 7L )
			paste(str_out(seed, 3),  str_c('[', digest(seed), ']'))
		else 
			str_out(seed, Inf)
	}
	else
		paste(class(seed), ' [', digest(seed), ']', sep='')
}

#' \code{RNGtype} extract the kinds of RNG and Normal RNG.
#'  
#' \code{RNGtype} returns the same type of values as \code{RNGkind()}, except that 
#' it can extract the RNG settings from an object.
#' If \code{object} is missing it returns the kinds of the current RNG settings, 
#' i.e. it is identical to \code{RNGkind()}.
#' 
#' @export
#' @rdname rng
RNGtype <- function(object){
	
	if( missing(object) ){
		RNGkind()
	}else{
		# extract RNG
		object <- getRNG(object)
		# get kind
		kinds <- c("Wichmann-Hill", "Marsaglia-Multicarry", "Super-Duper", 
				"Mersenne-Twister", "Knuth-TAOCP", "user-supplied", "Knuth-TAOCP-2002", 
				"L'Ecuyer-CMRG", "default")
		K <- kinds[object[1L] %% 100 + 1L]
		# get normal kind
		n.kinds <- c("Buggy Kinderman-Ramage", "Ahrens-Dieter", "Box-Muller", 
				"user-supplied", "Inversion", "Kinderman-Ramage", "default")
		NK <- n.kinds[floor(object[1L]/100)+1L]
		
		# return both kinds
		c(K, NK)
	}
	
}

#' \code{RNGinfo} shows displays human readable information about RNG settings.
#' If \code{object} is missing it displays information about the current RNG.
#' 
#' @param indent character string to use as indentation prefix in the output 
#' from \code{RNGinfo}.
#' 
#' @export
#' @rdname rng
RNGinfo <- function(object=getRNG(), indent=''){
	
	# get kind
	kind <- RNGtype(object)
	# determine provider
	prov <- RNGprovider()
	prov <- 
	if( prov == 'base')	kind[1L]
	else paste('package:', prov, sep='')

	# show information
	cat(indent, "RNG kind: ", paste(kind, collapse=" / "), "\n")
	cat(indent, "RNG state:", RNGdesc(object), "\n")

} 


.getRNGattribute <- function(object){
	if( .hasSlot(object, 'rng') ) slot(object, 'rng')
	else if( .hasSlot(object, 'rng.seed') ) slot(object, 'rng.seed') # for back compatibility
	else attr(object, 'rng')
}

#' RNG Settings
#' 
#' The functions documented here provide a unified interface to work with 
#' RNG settings.
#' 
#' \code{getRNG} is an S4 generic that returns the Random Number Generator (RNG) 
#' settings used for computing an object.
#' For example, in the case of results from multiple NMF runs, it returns the 
#' RNG settings used to compute the best fit.
#' 
#' @param object an R object from which RNG settings can be extracted, e.g. an 
#' integer vector containing a suitable value for \code{.Random.seed} or an 
#' object returned by the function \code{\link{nmf}}.
#' @param ... extra arguments passed to a suitable method \code{.getRNG}.
#' 
#' @return the RNG settings as a single integer vector as \code{\link{.Random.seed}} 
#' or \code{NULL} if no RNG data was found.
#' 
#' @rdname rng
#' @export
getRNG <- function(object, ...){
	
	if( missing(object) || is.null(object) ) return( .getRNG() )
	
	# use RNG data from object if available
	rng <- .getRNGattribute(object)
	if( !is.null(rng) ) getRNG(rng, ...)
	else if( isNumber(object)  ){
		nextRNG(object, ...) # return RNG as if after setting seed
	}else .getRNG(object, ...) # call S4 method on object
	
}

#' \code{getRNG} is an S4 generic that returns the Random Number Generator (RNG) 
#' settings used for computing an object.
#' For example, in the case of results from multiple NMF runs, it returns the 
#' RNG settings used to compute the best fit.
#' 
#' @param object an R object from which RNG settings can be extracted, e.g. an 
#' integer vector containing a suitable value for \code{.Random.seed} or an 
#' object returned by the function \code{\link{nmf}}.
#' @param ... extra arguments to allow extension
#' 
#' @return the RNG settings as a single integer vector as \code{\link{.Random.seed}} 
#' or \code{NULL} if no RNG data was found.
#' 
#' @rdname rng
#' @export
setGeneric('.getRNG', function(object, ...) standardGeneric('.getRNG') )
#' Returns the current RNG settings.
#' 
#' @examples 
#' # get current RNG settings
#' head(getRNG())
#' 
setMethod('.getRNG', 'missing',
	function(object){
		
		# return current value of .Random.seed
		# ensuring it exists first 
		if( !exists('.Random.seed', .GlobalEnv) ) 
			sample(NA)
		
		return( get('.Random.seed', .GlobalEnv) )
		
	}
)

#' Default method that tries to extract RNG information from \code{object}, by 
#' looking sequentially to a slot named \code{'rng'}, a slot named \code{'rng.seed'}
#' or an attribute names \code{'rng'}.
#' 
#' It returns \code{NULL} if no RNG data was found.
setMethod('.getRNG', 'ANY',
	function(object, ...){
		.getRNGattribute(object)
	}
)
#' Method for S3 objects, that aims at reproducing the behaviour of the function 
#' \code{getRNG} of the package \code{getRNG}. 
#' 
#' It sequentially looks for RNG data in elements \code{'rng'}, \code{noise$rng} 
#' if element \code{'noise'} exists and is a \code{list}, or in attribute \code{'rng'}.
#'  
setMethod('.getRNG', 'list',
	function(object){
		# lookup for some specific elements
		if( !is.null(object$rng) ) object$rng  
		else if( is.list(object$noise) ) object$noise$rng
		else attr(object, 'rng')
	}
)
#setMethod('.getRNG', 'rstream',
#		function(object){
#			object	
#		}
#)
#' Methods for numeric object, which returns the object itself, if it has more than one 
#' element, coerced into an integer vector if necessary, as it is assumed to 
#' already represent a value for \code{\link{.Random.seed}}.
#' 
#' Or if \code{object} has a single element, the value of \code{.Random.seed} as 
#' it would be after calling \code{set.seed(object, ...)}
#' In this case, all arguments in \code{...} are passed to \code{\link{set.seed}}.
#'    
#' @return \code{getRNG}, \code{getRNG1}, \code{nextRNG} and \code{setRNG} return 
#' an integer vector of length greater than 3 (see \code{\link{.Random.seed}}.
#'  
setMethod('.getRNG', 'numeric',
	function(object, ...){
		as.integer(object)
	}
)

#' \code{getRNG1} is an S4 generic that returns the \strong{initial} RNG settings 
#' used for computing an object.
#' For example, in the case of results from multiple NMF runs, it returns the 
#' RNG settings used to compute the \emph{first} fit.
#' 
#' \code{getRNG1} is defined to provide separate access to the RNG settings as 
#' they were at the very beginning of a whole computation, which might differ 
#' from the RNG settings returned by \code{getRNG}, that allows to reproduce the  
#' result only.
#' 
#' Think of a sequence of separate computations, from which only one result is 
#' used for the result (e.g. the one that maximise a likelihood): 
#' \code{getRNG1} would return the RNG settings to reproduce the complete sequence
#' of computations, while \code{getRNG} would return the RNG settings necessary to 
#' reproduce only the computation whose result has maximum likelihood.  
#' 
#' @rdname rng
#' @export
#' 
setGeneric('getRNG1', function(object, ...) standardGeneric('getRNG1') )
#' Default method that is identical to \code{getRNG(object, ...)}.
setMethod('getRNG1', 'ANY',
	function(object, ...){
		getRNG(object, ...)
	}
)


#' \code{nextRNG} returns the RNG settings as they would be after seeding with 
#' \code{seed}.
#' 
#' @rdname rng
#' @export
nextRNG <- function(object, ...){

	# get/restore .Random.seed on.exit
	orseed <- RNGscope()
	on.exit(RNGscope(orseed))
	
	# return next state of current RNG if object is missing
	if( missing(object) ){
		runif(1)
		return( getRNG() )
	}
	
	# extract RNG from object
	rng <- .getRNGattribute(object)
	if( !is.null(rng) ){
		on.exit()
		return( nextRNG(rng, ...) )
	}
	
	# only work for numeric seeds
	if( !is.numeric(object) )
		stop("Invalid seed: expecting a numeric seed.")
	
	# set RNG 
	.setRNG(object, ...)
	
	# return new RNG settings
	RNGscope()
}

.collapse <- function(x, sep=', ', n){
	
	res <- paste(if( missing(n) ) x else head(x, n), collapse=', ')
	if( length(x) > n )
		res <- paste(res, '...', sep=', ')
	res
}

.showRNG <- function(...){	
	message(paste(capture.output( RNGinfo(..., indent="#") ), collapse="\n"))	
}

RNGshow <- .showRNG

#' \code{setRNG} tries to extract RNG settings from \code{object} and use a 
#' suitable \code{.setRNG} method to set these settings.
#' All arguments are passed to the next call to \code{.setRNG}.
#'
#' @return \code{setRNG} invisibly returns the old RNG settings as 
#' they were before changing them.
#' 
#' @export 
#' @rdname rng
#' @examples 
#' 
#' obj <- list(x=10, rng=123)
#' setRNG(obj)
#' rng <- getRNG()
#' runif(10)
#' set.seed(123)
#' rng.equal(rng)
#' 
setRNG <- function(object, ..., verbose=FALSE){
	
	# do nothing if null
	if( is.null(object) ) return()
	
	# use RNG data from object if available
	rng <- getRNG(object, ...)
	if( !is.null(rng) && !identical(rng, object) ) return( setRNG(rng, ...) )
	
	# get/restore .Random.seed on.exit in case of errors
	orseed <- getRNG()
	on.exit({
		message("Restoring RNG settings probably due to an error in setRNG")
		RNGscope(orseed) 
	})

	# call S4 method on object
	.setRNG(object, ...) 
	
	# cancel RNG restoration
	on.exit()
	if( verbose ) .showRNG()			
	
	invisible(orseed)
}

#' \code{.setRNG} is an S4 generic that sets the current RNG settings, from a 
#' variety of format.
#' Its methods define the workhorse functions that are called by \code{setRNG}.
#' 
#' @inline
#' @rdname rng
#' @export 
setGeneric('.setRNG', function(object, ...) standardGeneric('.setRNG') )
#' Sets the RNG to kind \code{object}, assuming is a valid RNG kind:
#' it is equivalent to \code{RNGkind(object, ...}.
#' All arguments in \code{...} are passed to \code{\link{RNGkind}}.
#' 
#' @param verbose a logical that indicates if the new RNG settings should
#' be displayed.
#' 
#' @examples
#' # set RNG kind
#' old <- setRNG('Marsaglia')
#' # restore
#' setRNG(old)
setMethod('.setRNG', 'character',
	function(object, ...){
		RNGkind(kind=object, ...)
	}
)

#' Sets the RNG settings using \code{object} directly the new value for 
#' \code{.Random.seed} or to initialise it with \code{\link{set.seed}}.
#' 
#' @examples 
#' 
#' # directly set .Random.seed
#' rng <- getRNG()
#' r <- runif(10)
#' setRNG(rng)
#' rng.equal(rng)
#' 
#' # initialise from a single number (<=> set.seed)
#' setRNG(123)
#' rng <- getRNG()
#' runif(10)
#' set.seed(123)
#' rng.equal(rng)
#' 
setMethod('.setRNG', 'numeric',
	function(object, ...){
		
		seed <- as.integer(object)
		if( length(seed) == 1L ){
			set.seed(seed, ...)
		}else{
			assign('.Random.seed', seed, envir=.GlobalEnv)
			# check validity of the seed
			tryCatch(runif(1)
			, error=function(err){					
				stop("setRNG - Invalid value for .Random.seed ["
					, .collapse(seed, n=5), "]: ", err$message, '.'
					, call.=FALSE)
			})
			assign('.Random.seed', seed, envir=.GlobalEnv)			
		}
	}
)

#' \code{RNGdigest} computes a hash from the RNG settings associated with an 
#' object. 
#' 
#' @rdname rng
#' @export
RNGdigest <- function(x){
	
	object <- if( missing(x) )	getRNG() else getRNG(x)
	
	# exit if no RNG was extracted
	if( is.null(object) )
		return(digest(NULL)) # TODO: return NULL
		
	digest(object)
	
}

#' \code{rng.equal} and \code{rng1.equal} return \code{TRUE} the RNG settings 
#' associated with two objects are identical, and \code{FALSE} otherwise. 
#' The comparison is made between the hashes returned by \code{RNGdigest}.
#' 
#' @param x objects from which RNG settings are extracted
#' @param y object from which RNG settings are extracted
#' 
#' @return \code{rng.equal} and \code{rng.equal1} return a \code{TRUE} or 
#' \code{FALSE}.
#' 
#' @rdname rng
#' @export
rng.equal <- function(x, y){
	if( missing(y) )
		y <- getRNG()
	identical(RNGdigest(x), RNGdigest(y))
}

#' The function \code{rng1.equal} tests whether two objects have identical 
#' \strong{initial} RNG settings.
#' 
#' @rdname rng
#' @export
rng1.equal <- function(x, y){
	if( missing(y) )
		y <- getRNG()
	rng.equal(getRNG1(x), getRNG1(y))
}
