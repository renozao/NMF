# Definitions used in the parallel computations of NMF
#
# - reproducible backend
# - reproducible %dopar% operator: %dorng%
# 
# Author: Renaud Gaujoux
# Creation: 08-Feb-2011
###############################################################################

###% Returns a single string that describes the given RNG state
###% 
RNGdesc <- function(seed){
	
	if( missing(seed) ){
		rp <- RNGprovider()
		if( rp == 'base' )
			seed <- .Random.seed
		else if( rp == 'rstream' )
			return( paste(RNGstate(which='state'), collapse=', ') )
		else if( rp == 'rlecuyer' ){
			stream <- .Call("r_get_current_stream",PACKAGE="rlecuyer" )
			return( paste(stream[[1]][1:6], collapse=', ') )
		}			
		else
			return(as.character(.Random.seed))
			
	}
	
	if( is.null(seed) )
		return('NULL')
	else if( is.numeric(seed) ){
		if( length(seed) == 1 )
			as.character(seed)
		else if( length(seed) == 6 )
			paste(seed, collapse=', ')
		else
			digest(seed)
	}
	else if( is(seed, 'rstream') ){
		if( is(seed, 'rstream.mrg32k3a') )
			paste(class(seed), ' - ', paste(RNGstate(seed, 'seed'), collapse=", "), ' [', RNGdigest(seed), ']', sep='')
		else
			paste(class(seed), ' [', RNGdigest(seed), ']', sep='')
	}else
		paste(class(seed), ' [', digest(seed), ']', sep='')
}


###% Returns the state of a random stream
###% 
RNGstate <- function(object, which){
	
	object <- if( missing(object) ) getRNG() else getRNG(object)
	
	if( !is(object, 'rstream.mrg32k3a') )
		stop("RNGstate - Invalid object: only class 'rstream.mrg32k3a' is supported")
	
	if( !rstream.packed(object) ){
		state <- .Call("R_RngStreams_GetData", object@stream, PACKAGE="rstream")
	}else
		state <- object@pack$state
	
	state <- list(seed=state[13:18], substream=state[7:12], state=state[1:6])
	if( missing(which) )
		state
	else{
		which <- match.arg(which, names(state))
		state[[which]]
	}
}

###% Returns the default seed used to generate the next random stream
.rstream.get.seed <- function(){
	get(".rstream.mrg32k3a.DefaultSeed", envir=rstream:::.rstream.envir)
}

.rstream.set.seed <- function(seed){
	
	# check sed validity
	seed <- rstream:::.rstream.mrg32k3a.CheckSeed(seed)
	# retrieve current value
	old <- get(".rstream.mrg32k3a.DefaultSeed", envir=rstream:::.rstream.envir)
	
	## save seed in rstream library
	.Call("R_RngStreams_SetPackageSeed", as.double(seed), PACKAGE="rstream")	
	## save seed as R variable
	assign(".rstream.mrg32k3a.DefaultSeed",	as.double(seed), envir=rstream:::.rstream.envir)
	assign(".rstream.mrg32k3a.HasSeed", TRUE, envir=rstream:::.rstream.envir)
	
	# return old seed
	invisible(old)
}

###% Set the RNG for use by the function nmf.
###% 
###% It returns the old RNG as an rstream object or the result of set.seed 
###% if the RNG is not changed due to one of the following reason:
###% - the settings are not compatible with rstream  
###%
.setupRNG <- function(seed, n=1, use.streams, verbose=FALSE, ...){

	if( use.streams || is(seed, 'rstream') || (is.numeric(seed) && length(seed) == 6) )
		.setupRNGrstream(seed=seed, n=n, verbose=verbose, ...)
	else{
		# immediately setup the RNG in the standard way (except if seed is NA)
		if( isNA(seed) ) seed <- 'random' 
		rng <- .setupRNGstd(seed=seed, verbose=verbose, ...)
		
		# build fake sequence: RNG for first run is the RNG that's just been setup 		
		seq <- if( n == 1 )	rng
				else c(list(rng), replicate(n-1, NULL))
		
		# reset seeding method only if object `seed` was meaningful for seeding
		if( !is.null(rng) )
			list(method='random', seq=seq)
		else
			list(method=seed, seq=seq)
	}
} 

.setupRNGrstream <- function(seed, n=1, unlist=TRUE, verbose=FALSE){
	
	
	if( verbose > 2 ) message("# RNG mode: reproducible [using RNGstream]")
	
	# init result list
	res <- list(method=seed, seq=NULL)
	
	# check and possibly generate a numeric seed for RNGstream
	use.random <-
	if( is.numeric(seed) ) TRUE
	else if( is(seed, 'rstream') ){ # use the current state of the given stream
		
		## For single run: directly set the rstream object and return the current RNG
		# or NULL if the RNG was not set
		if( n==1 ){
			if( verbose > 2 ) message("# Set current RNG ... ", appendLF=verbose>3) 
			.tryCatch_setRNG(seed, verbose=verbose>3)
			if( verbose > 2 ) message("OK")
			res$method <- 'random'
			res$seq <- seed
			return(res)
		}
		
		
		if( !is(seed, 'rstream.mrg32k3a') )
			stop("NMF::nmf - Invalid seed object for multiple runs: only rstream objects of class 'rstream.mrg32k3a' are supported")
		
		# generate a sequence using the stream as the first element 
		seed <- RNGstate(seed, 'seed')
		
		TRUE
	}else if( isNA(seed) ) TRUE # keep value NA => using the current next seed of the rstream
	else{
		seed <- NULL # will generate a random seed for rstream
		FALSE
	}

#		# show the new state of the current RNG in verbose mode
#		if( verbose > 2 ){
#			message("OK")
#			message("# ** New current RNG settings:")
#			.showRNG()	
#		}
		
	# change the actual seeding method to 'random' if relevant
	if( use.random ) res$method <- 'random'
	
	# show the used random seed in verbose mode
	if( verbose > 2 ) message("# Generate RNG stream(s) using seed: ", RNGdesc(seed), " ... ", appendLF=verbose>3)
	# generate a sequence of streams
	res$seq <- RNGseq(n, seed, packed=TRUE
					, prefix=paste(basename(tempfile('NMF')), 'run', sep='_')
					, unlist=unlist, verbose=verbose>3)	
	if( verbose > 2 ) message("OK")
	
	## For single run: directly set the rstream object and return the current RNG
	# or NULL if the RNG was not set
	if( n==1 ){
		if( verbose > 2 ) message("# Set current RNG ... ", appendLF=verbose>3) 
		.tryCatch_setRNG(res$seq, verbose=verbose>3)
		if( verbose > 2 ) message("OK")
	}
	
	# return the sequence of RNGs
	res
}

# Try-catch wrapped setRNG used for setting rstream object and get an informative 
# message in case of errors
.tryCatch_setRNG <- function(...){
	tryCatch(setRNG(...)
			, error = function(e){				
				stop(e$message, "\n\tLibrary rstream is probably in conflict with another RNG library."
						,"\n\tRe-running with .options='v4' may provide more debugging information."
						,"\n\tUsing .options='-R' with a single numeric seed should remove this error.", call.=FALSE)
			}
	)
}

.setupRNGstd <- function(seed, verbose=FALSE){
	
	if( verbose > 2 ) message("# RNG mode: standard [using set.seed]")
	
	rng <- 
	if( is.numeric(seed) ){
		
		if( length(seed) != 1 )
			stop("NMF::nmf - Invalid numeric seed: expects a single numeric value or a 6-length numeric vector")
		
		if( verbose > 2 )
			message("# Set the current RNG with set.seed(", seed,") ... ", appendLF=FALSE)
			
		# set the RNG with standard set.seed
		set.seed(seed)
		getRNG()
		
	}
	
	# in verbose > 2: ouptut the new RNG settings
	if( !is.null(rng) && verbose > 2 ){
		message("OK")
		message("# ** New RNG settings:")
		.showRNG()		
	}
	
	# return new RNG (it is NULL if no changes were made)
	rng
}

###% Create a given number of rstream objects to be used as random number generators
###% for each NMF run when performing multiple runs.
###% 
###% This ensures complete reproducibility of the set of run. 
###% The streams are created using the RNGstream C++ package (from P. L'Ecuyer), 
###% using the interface provided by the R package rstream.
###% 
###% If a seed is provided, the original rstream seed is restored on exit, so that
###% the consistency of the global sequence of streams generated is not jeopardised.
###% 
###% @param n Number of streams to be created
###% @param seed Numerical seed used to initialise the set of streams.
###% @param packed Logical. If TRUE the streams are returned packed (see \code{\link{rstream.packed}}).
###%  
RNGseq <- function(n, seed=NA, packed=TRUE, prefix=NULL, unlist=TRUE, verbose=FALSE){
	
	# check parameters
	if( n <= 0 )
		stop("NMF::createStream - invalid value for 'n' [positive value expected]")
	
	# force the initial seed if provided
	if( !isNA(seed) ){
		oldseed <- RNGseed(seed, verbose=verbose)
		on.exit({.rstream.set.seed(oldseed)}, add=TRUE)
	}
		
	# generate the sequence of streams
	res <- lapply(1:n, function(i){
				s <- new('rstream.mrg32k3a', name=if( !is.null(prefix) ) paste(prefix, i, sep="_") else NULL );
				rstream.packed(s) <- packed;
				s}
	)
	
	# return list or single RNG
	if( n==1 && unlist )
		res[[1]]
	else
		invisible(res)
	
}
#seq.rng <- RNGseq

##################################################################
## END
##################################################################
