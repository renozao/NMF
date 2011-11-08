# Definitions used in the parallel computations of NMF
#
# - reproducible backend
# - reproducible %dopar% operator: %dorng%
# 
# Author: Renaud Gaujoux
# Creation: 08-Feb-2011
###############################################################################


####% Returns the default seed used to generate the next random stream
#.rstream.get.seed <- function(){
#	get(".rstream.mrg32k3a.DefaultSeed", envir=rstream:::.rstream.envir)
#}
#
#.rstream.set.seed <- function(seed){
#	
#	# check sed validity
#	seed <- rstream:::.rstream.mrg32k3a.CheckSeed(seed)
#	# retrieve current value
#	old <- get(".rstream.mrg32k3a.DefaultSeed", envir=rstream:::.rstream.envir)
#	
#	## save seed in rstream library
#	.Call("R_RngStreams_SetPackageSeed", as.double(seed), PACKAGE="rstream")	
#	## save seed as R variable
#	assign(".rstream.mrg32k3a.DefaultSeed",	as.double(seed), envir=rstream:::.rstream.envir)
#	assign(".rstream.mrg32k3a.HasSeed", TRUE, envir=rstream:::.rstream.envir)
#	
#	# return old seed
#	invisible(old)
#}

###% Set the RNG for use by the function nmf.
###% 
###% It returns the old RNG as an rstream object or the result of set.seed 
###% if the RNG is not changed due to one of the following reason:
###% - the settings are not compatible with rstream  
###%
.setupRNG <- function(seed, n, verbose=FALSE){
	
	# for multiple runs one always uses RNGstreams
	if( n > 1 ){
		
		if( verbose > 2 ) message("# RNG setup: reproducible [using RNGstream]")
		# seeding with numeric values only
		if( is.numeric(seed) ){
			
			if( verbose > 2 ){
				message("# Generate RNGStream sequence using seed ("
						, RNGdesc(seed), ") ... "
						, appendLF=FALSE)
			}
			res <- RNGseq(n, seed)
			if( verbose > 2 ) message("OK")
			return(res)
			
		}else{ # create a sequence of RNGstream using a random seed
			if( verbose > 2 ){
				message("# Generate RNGStream sequence using a random seed ... "
						, appendLF=FALSE)
			}
			res <- RNGseq(n, NULL)
			if( verbose > 2 ) message("OK")
			return(res)
		}
	}else if( is.numeric(seed) ){ 
		# for single runs: 1-length seeds are used to set the current RNG
		# 6-length seeds are used to set RNGstream
		
		# convert to an integer vector
		seed <- as.integer(seed)
		# immediately setup the RNG in the standard way		
		if( length(seed) == 1L ){
			if( verbose > 2 ){
				message("# RNG setup: standard [seeding current RNG]")
				message("# Seeding current RNG with seed (", seed, ") ... "
						, appendLF=FALSE)
			}
			set.seed(seed)
			if( verbose > 2 ) message("OK")				
			return( RNGseed() )
		}else if( length(seed) == 6L ){
			if( verbose > 2 ){
				message("# RNG setup: reproducible [using RNGstream]")
				message("# Generate RNGStream sequence using seed ("
						, RNGdesc(seed), ") ... "
						, appendLF=FALSE)
			}
			res <- RNGseq(1, seed)
			setRNG(res)
			if( verbose > 2 ) message("OK")
			return( res )
		}else{
			if( verbose > 2 ){
				message("# RNG setup: directly setting RNG")
				message("# Setting RNG with .Random.seed= ("
						, RNGdesc(seed), ") ... "
						, appendLF=FALSE)
			}
			setRNG(seed, verbose > 2)
			if( verbose > 2 ) message("OK")
			return( RNGseed() )
		}
			stop("NMF::nmf - Invalid numeric seed: expects a single numeric value or a 6-length numeric vector")		
	}else{
		if( verbose > 2 ) message("# RNG setup: standard [using current RNG]")
		NULL
	}
} 

#.setupRNGrstream <- function(seed, n=1, unlist=TRUE, verbose=FALSE){
#	
#	
#	if( verbose > 2 ) message("# RNG mode: reproducible [using RNGstream]")
#	
#	# init result list
#	res <- list(method=seed, seq=NULL)
#	
#	# check and possibly generate a numeric seed for RNGstream
#	use.random <-
#	if( is.numeric(seed) ) TRUE
#	else if( is(seed, 'rstream') ){ # use the current state of the given stream
#		
#		## For single run: directly set the rstream object and return the current RNG
#		# or NULL if the RNG was not set
#		if( n==1 ){
#			if( verbose > 2 ) message("# Set current RNG ... ", appendLF=verbose>3) 
#			.tryCatch_setRNG(seed, verbose=verbose>3)
#			if( verbose > 2 ) message("OK")
#			res$method <- 'random'
#			res$seq <- seed
#			return(res)
#		}
#		
#		
#		if( !is(seed, 'rstream.mrg32k3a') )
#			stop("NMF::nmf - Invalid seed object for multiple runs: only rstream objects of class 'rstream.mrg32k3a' are supported")
#		
#		# generate a sequence using the stream as the first element 
#		seed <- RNGstate(seed, 'seed')
#		
#		TRUE
#	}else if( isNA(seed) ) TRUE # keep value NA => using the current next seed of the rstream
#	else{
#		seed <- NULL # will generate a random seed for rstream
#		FALSE
#	}
#
##		# show the new state of the current RNG in verbose mode
##		if( verbose > 2 ){
##			message("OK")
##			message("# ** New current RNG settings:")
##			.showRNG()	
##		}
#		
#	# change the actual seeding method to 'random' if relevant
#	if( use.random ) res$method <- 'random'
#	
#	# show the used random seed in verbose mode
#	if( verbose > 2 ) message("# Generate RNG stream(s) using seed: ", RNGdesc(seed), " ... ", appendLF=verbose>3)
#	# generate a sequence of streams
#	res$seq <- RNGseq(n, seed, packed=TRUE
#					, prefix=paste(basename(tempfile('NMF')), 'run', sep='_')
#					, unlist=unlist, verbose=verbose>3)	
#	if( verbose > 2 ) message("OK")
#	
#	## For single run: directly set the rstream object and return the current RNG
#	# or NULL if the RNG was not set
#	if( n==1 ){
#		if( verbose > 2 ) message("# Set current RNG ... ", appendLF=verbose>3) 
#		.tryCatch_setRNG(res$seq, verbose=verbose>3)
#		if( verbose > 2 ) message("OK")
#	}
#	
#	# return the sequence of RNGs
#	res
#}

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

#.setupRNGstd <- function(seed, verbose=FALSE){
#	
#	if( verbose > 2 ) message("# RNG mode: standard [using set.seed]")
#	
#	rng <- 
#	if( is.numeric(seed) ){
#		
#		if( length(seed) != 1 )
#			stop("NMF::nmf - Invalid numeric seed: expects a single numeric value or a 6-length numeric vector")
#		
#		if( verbose > 2 )
#			message("# Set the current RNG with set.seed(", seed,") ... ", appendLF=FALSE)
#			
#		# set the RNG with standard set.seed
#		set.seed(seed)
#		getRNG()
#		
#	}
#	
#	# in verbose > 2: ouptut the new RNG settings
#	if( !is.null(rng) && verbose > 2 ){
#		message("OK")
#		message("# ** New RNG settings:")
#		.showRNG()		
#	}
#	
#	# return new RNG (it is NULL if no changes were made)
#	rng
#}

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
RNGseq <- function(n, seed=NULL, unlist=TRUE){
	
	# check parameters
	if( n <= 0 )
		stop("RNGseq - Invalid value for 'n' [positive value expected]")
	if( is.null(seed) )
		seed <- integer(0)	
	if( !is.numeric(seed) )
		stop("RNGseq - Invalid seed: expected NULL or a numeric value [", class(seed) ,"]")	
	
	# get/restore .Random.seed on.exit
	orseed <- RNGscope()
	on.exit( RNGscope(orseed), add=TRUE)
	
	seed <- as.integer(seed)
	
	# detect and extract RNGstream seed
	if( length(seed) == 7L && seed[1L] %% 100L == 7L )
		seed <- seed[-1]
	
	seed <- 
	if( length(seed) == 0L ){
		# draw one value of restored RNG on.exit
		on.exit( runif(1), add=TRUE)
		RNGkind(kind = "L'Ecuyer-CMRG")
		RNGseed()
	}else if( length(seed) == 1L ){		
		set.seed(seed, kind = "L'Ecuyer-CMRG")
		RNGseed()
	}else if( length(seed) == 6L ){
		RNGkind(kind = "L'Ecuyer-CMRG")
		c(RNGseed()[1], seed)
	}else{ # use directly as .Random.seed (we know it is not a CMRG seed)
		setRNG(seed)
		RNGkind(kind = "L'Ecuyer-CMRG")
		RNGseed()
	}		

	# generate the sequence of RNG streams
	res <- as.list(rep(NA, n))
	res[[1]] <- seed
	if( n > 1 ){		
		for(i in 2:n){
			res[[i]] <- nextRNGStream(res[[i-1]])
		}
	}
	stopifnot( !any(is.na(res)) )
	
	# return list or single RNG
	if( n==1 && unlist )
		res[[1]]
	else
		res
	
}
#seq.rng <- RNGseq

##################################################################
## END
##################################################################
