# TODO: Add comment
# 
# Author: Renaud Gaujoux
# Creation: 08 Nov 2011
###############################################################################

RNGrecovery <- function(){
	
	s <- as.integer(c(401,0,0))
	assign(".Random.seed", s, envir=.GlobalEnv)
	RNGkind("default")
	
}

baseRNGkind <- base::RNGkind

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
	
	kind <- baseRNGkind()
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
				assign('.Random.seed', seed, .GlobalEnv)
			}
	invisible(res)
}

###% Returns a single string that describes the given RNG state
###% 
RNGdesc <- function(seed){
	
	if( missing(seed) ){
		rp <- RNGprovider()
		if( rp == 'base' || length(.Random.seed) > 1L )
			seed <- .Random.seed
		else 
			return( "Unknown" )		
	}
	
	if( is.null(seed) )
		return('NULL')
	else if( is.numeric(seed) ){
		paste(paste(head(seed, 7), collapse=', ')
			, if( length(seed) > 7L ) paste(', ... [', digest(seed), ']', sep=''), sep='')
	}
	else
		paste(class(seed), ' [', digest(seed), ']', sep='')
}

RNGinfo <- function(object, prefix=''){
	if( missing(object) ){
		cat(prefix, "RNG kind: ", paste(baseRNGkind(), collapse=" / "), "\n")								
		cat(prefix, "RNG state: ", RNGdesc(), "\n")
	}else{
		object <- getRNG(object)
		kinds <- c("Wichmann-Hill", "Marsaglia-Multicarry", "Super-Duper", 
				"Mersenne-Twister", "Knuth-TAOCP", "user-supplied", "Knuth-TAOCP-2002", 
				"L'Ecuyer-CMRG", "default")
		K <- kinds[object[1] %% 100 + 1]
		
		n.kinds <- c("Buggy Kinderman-Ramage", "Ahrens-Dieter", "Box-Muller", 
				"user-supplied", "Inversion", "Kinderman-Ramage", "default")
		NK <- n.kinds[floor(object[1]/100)+1] 
		
		prov <- RNGprovider()
		prov <- 
		if( prov == 'base')	K
		else paste('pkg:', prov, sep='') 		
		cat(prefix, "RNG kind: ", paste(prov, NK, sep=" / "), "\n")
		cat(prefix, "RNG state:", RNGdesc(object), "\n")
	}
} 


# Get the RNG Settings 
setGeneric('getRNG', function(object, ...) standardGeneric('getRNG') ) 
setMethod('getRNG', 'missing',
		function(object){#}, packed=FALSE){
			
			# return current value of .Random.seed
			return( RNGseed() )
			
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
#setMethod('getRNG', 'rstream',
#		function(object){
#			object	
#		}
#)
setMethod('getRNG', 'integer',
		function(object){
			object	
		}
)


# Returns the current or a future random seed 
RNGseed <- function(seed, ...){
	
	if( missing(seed) ){
		if( !exists('.Random.seed', .GlobalEnv) )
			sample(NA)
		return( get('.Random.seed', .GlobalEnv) )
	}
	
	# only work for numeric seeds
	if( !is.numeric(seed) )
		stop("Invalid seed: expecting a numeric seed.")
	
	# get/restore .Random.seed on.exit
	orseed <- RNGscope()
	on.exit( RNGscope(orseed), add=TRUE)
	
	setRNG(seed)
	RNGseed()	
}


###% Generic function that sets the Random Number Generator
###% 
###% 
setGeneric('setRNG', function(object, ...) standardGeneric('setRNG') )
setMethod('setRNG', 'character',
		function(object, ..., verbose=FALSE){
			
			old <- RNGseed()
			baseRNGkind(kind=object, ...)
			old
		}
)

setMethod('setRNG', 'numeric',
		function(object, verbose=FALSE, ...){
			
			# get/restore .Random.seed on.exit in case of errors
			orseed <- RNGscope()
			on.exit({
				if( verbose ) message("Restoring RNG settings probably due to an error")
				RNGscope(orseed) 
			})
			
			# retrieve current seed
			old <- RNGseed()
			
			seed <- as.integer(object)
			if( length(seed) == 1L ){
				set.seed(seed, ...)
			}else{			
				assign('.Random.seed', seed, .GlobalEnv)
				# check validity of the seed
				tryCatch(runif(1)
				, error=function(err){
					message('')
					stop("setRNG - Invalid value for .Random.seed [", err$message, ']', call.=FALSE)
				})
				assign('.Random.seed', seed, .GlobalEnv)			
			}
			
			# cancel RNG restoration
			on.exit()
			if( verbose ) .showRNG()			
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

RNGdigest <- function(x){
	
	object <- if( missing(x) )	getRNG() else getRNG(x)
	
	# exit if no RNG was extracted
	if( is.null(object) )
		return(digest(NULL)) # TODO: return NULL
		
	digest(object)
	
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