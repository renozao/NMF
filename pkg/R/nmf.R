###% @include NMFSet-class.R
###% @include NMFStrategy-class.R
NA

isNMFResult <- function(x){
	is(x, 'NMFfit') || is(x, 'NMFfitX')
}

###% Generic interface to Nonnegative Matrix Factorization algorithms.
###%
###% Function \code{nmf} is the main entry point to perform NMF algorithms defined within framework
###% set up by package \code{NMF}.
###% It provides an interface to combine the different algorithms with the different seeding methods. It returns
###% the result as an object of class \code{NMF} that can be directly passed to visualization or benchmarking methods.
###% 
###% The default behaviour of \code{nmf} when \code{method} is missing is to use the algorithm
###% from Brunet et al., which is implemented as a predefined NMFStrategy.
###%
###% @seealso NMFStrategy-class
###% setGeneric('nmf', function(x, rank, method=nmf.getOption('default.algorithm'), ...) standardGeneric('nmf') )
setGeneric('nmf', function(x, rank, method, ...) standardGeneric('nmf') )

###% Performs NMF on a data.frame: the target matrix is converted data.frame \code{as.matrix(x)}.
setMethod('nmf', signature(x='data.frame', rank='ANY', method='ANY'), 
	function(x, rank, method, ...)
	{
		# replace missing values by NULL values for correct dispatch
		if( missing(method) ) method <- NULL
		if( missing(rank) ) rank <- NULL
		
		# apply NMF to the the data.frame converted into a matrix	
		nmf(as.matrix(x), rank, method, ...)
	}
)

###% Performs NMF using an already defined NMF model as starting point (no rank provided)
setMethod('nmf', signature(x='matrix', rank='ANY', method='ANY'), 
		function(x, rank, method, seed, ...)
		{
			# if rank is missing or NULL: seed must be there to specify the model 
			# AND algorithm 
			if( missing(rank) || is.null(rank) ){
				# argument seed must be supplied as a NMF object			
				if( missing(seed) || !inherits(seed, 'NMF') )
					stop("NMF::nmf : when argument 'rank' is not provided, argument 'seed' is required to inherit from class 'NMF'. See ?nmf."
						, call.=FALSE)
			
				rank <- nbasis(seed)
			}
			else if( !is.numeric(rank) )
				stop("NMF::nmf - argument 'rank' must be numeric.", call.=FALSE)
			
			# if method is missing or NULL try to find the correct algorithm 
			# from the seed's model
			if( missing(method) || is.null(method) ){
				
				# a priori the default method will be used
				method <- nmf.getOption('default.algorithm')
				
				# try to find the algorithm suitable for the seed's NMF model
				if( !missing(seed) && inherits(seed, 'NMF') ){
					
					method.potential <- nmfAlgorithm(model=modelname(seed))
					if( is.null(method.potential) )
						stop("NMF::nmf - No algorithm is defined for model '", modelname(seed), "'")
					
					if( length(method.potential) == 1 ) # only one to choose
						method <- method.potential
					else if( !is.element(method, method.potential) ) # several options, none is default
						stop("NMF::nmf - Could not infer the algorithm to use with model '", modelname(seed), "'. Argument 'method' should be one of: "
							, paste(paste("'", method.potential, "'", sep=''), collapse=', ')
							, call.=FALSE)
				}
					
			}
			else if( !hasMethod('nmf', signature=c('matrix', 'numeric', class(method)) ) )
				stop("NMF::nmf - no 'nmf' method for signature 'matrix', 'numeric', '", class(method), "'."
					, call.=FALSE)
			
			# use default seeding method if seed is missing
			if( missing(seed) )
				seed <- nmf.getOption('default.seed')
			
			nmf(x, rank, method, seed=seed, ...)
		}
)


###% Performs NMF on an object using a given list of algorithms.
setMethod('nmf', signature(x='matrix', rank='numeric', method='list'), 
	function(x, rank, method, ...)
	{
		# apply each NMF algorithm
		k <- 0
		t <- system.time({
			res <- lapply(method, 
				function(meth, ...){
					k <<- k+1
					message("Compute NMF method ", k, " ... ", appendLF=FALSE)					
					#o <- capture.output( 
							res <- try( nmf(x, rank, meth, ...) , silent=TRUE) 
					#)
					if( is(res, 'try-error') )
						message("ERROR")
					else 
						message("OK")
					return(res)
				}
				, ...)
		})
		
		# filter out bad results
		ok <- sapply(res, isNMFfit)
		if( any(!ok) ){ # throw warning if some methods raised an error
			err <- lapply(which(!ok), function(i){ paste("'", method[[i]],"': ", res[[i]], sep='')})
			warning("NMF::nmf - Incomplete results due to ", sum(!ok), " errors: \n- ", paste(err, collapse="- "), call.=FALSE)
		}
		res <- res[ok]
		# TODO error if ok is empty

		# add names to the result list
		names(res) <- sapply(res, algorithm)
				
		# wrap the result in a NMFSet object
		# DO NOT WRAP anymore here: NMFSet objects are used only for results of multiple runs (single method)
		# the user can still join on the result if he wants to
		#res <- join(res, runtime=t)
		res <- new('NMFList', res, runtime=t)
		
		# return result
		return(res)
	}
)

###% Performs NMF on a matrix using a predefined named strategy. 
###%
###% The available strategies are:
###% - brunet : Based from Brunet et al.
setMethod('nmf', signature(x='matrix', rank='numeric', method='character'),
function(x, rank, method, ...)
{	
	# if there is more than one methods then treat the vector as a list
	if( length(method) > 1 ) return( nmf(x, rank, as.list(method), ...) )
	
	# create the NMFStrategy from its name
	strategy <- nmfAlgorithm(method)		
	# apply nmf using the retrieved strategy		
	nmf(x, rank, method=strategy, ...)
}
)

###% Performs NMF on a matrix using a given function.
###%
setMethod('nmf', signature(x='matrix', rank='numeric', method='function'),
	function(x, rank, method, name, objective='euclidean', model='NMFstd', mixed=FALSE, ...){

		# build a NMFStrategyFunction object on the fly to wrap function 'method'
		model.name <- model
		model.parameters <- NULL
		
		# if model is a list then its element will be used to instanciate the NMF model
		# => it will be passed in argument 'model' to the main 'nmf' method
		if( is.list(model) ){
			
			# all elements are used, unless the first element is the name of a 
			# class that extends class 'NMF', then only the remaining elements
			# are used for the instantiation
			model.parameters <- model
			if( length(model) > 0 
				&& is.character(model[[1]]) && extends(model[[1]], 'NMF')){
				model.name <- model[[1]]
				# use the remaining elements to instanciate the NMF model
				model.parameters <- model[-1]
			}
			else model.name <- 'NMFstd'
			
		}
		
		# if name is missing: generate a temporary unique name
		if( missing(name) ) name <- basename(tempfile("NMF.algo."))
		# check that the name is not a registered name
		if( existsNMFAlgorithm(name) )
			stop("Invalid name for custom NMF algorithm: '",name,"' is already a registered NMF algorithm")
		
		# only use first element of mixed
		if( length(mixed) > 1 ){
			mixed <- mixed[1]
			warning("NMF::nmf : Only the first element of argument 'mixed' will be used [val=",mixed,"]")
		}
				
		strategy <- new('NMFStrategyFunction'
						, name=name, objective=objective, model=model.name
						, algorithm=method
						, mixed=mixed)
		# valid the strategy
		validObject(strategy, complete=TRUE)
		
		# call method 'nmf' with the new object
		nmf(x, rank, strategy, model=model.parameters, ...)
	}
)

.as.numeric <- function(x){
	suppressWarnings( as.numeric(x) )
}

.translate.string <- function(string, dict){
	
	res <- list()
	dict <- as.list(dict)
	if( nchar(string) == 0 ) return(res)
	opt.val <- TRUE
	last.key <- NULL
	buffer <- ''
	lapply(strsplit(string, '')[[1]], 
		function(c){
			if( c=='-' ) opt.val <<- FALSE
			else if( c=='+' ) opt.val <<- TRUE
			else if( opt.val && !is.na(.as.numeric(c)) )
				buffer <<- paste(buffer, c, sep='')
			else if( !is.null(dict[[c]]) ){
				# flush the buffer into the last key if necessary
				if( nchar(buffer) > 0 && !is.null(last.key) && !is.na(buffer <- .as.numeric(buffer)) ){
					res[[dict[[last.key]]]] <<- buffer
					buffer <<- ''
				}
					
				res[[dict[[c]]]] <<- opt.val
				last.key <<- c
			}
		}
	)

	# flush the buffer into the last key
	if( nchar(buffer) > 0 && !is.null(last.key) && !is.na(buffer <- .as.numeric(buffer)) )
		res[[dict[[last.key]]]] <- buffer

	# return result	
	return(res)
}

.showRNG <- function(...){	
	message(paste(capture.output( RNGinfo(..., prefix="#") ), collapse="\n"))	
} 

# locally define readRDS and saveRDS 
if( !existsFunction('readRDS') )
	readRDS <- .readRDS
if( !existsFunction('saveRDS') )
	saveRDS <- .saveRDS 

###% Performs NMF on a matrix using a given NMF method.
###%
###% This method is the entry point for NMF. It is eventually called by any definition of the \code{nmf} function.
setMethod('nmf', signature(x='matrix', rank='numeric', method='NMFStrategy'),
#function(x, rank, method, seed='random', nrun=1, keep.all=FALSE, optimized=TRUE, init='NMF', track, verbose, ...)
function(x, rank, method
		, seed=nmf.getOption('default.seed'), nrun=1, model=NULL, .options=list()
		, .pbackend=nmf.getOption('parallel.backend')
		, .callback=NULL #callback function called after a run  
		, ...)
{
	# if options are given as a character string, translate it into a list of booleans
	if( is.character(.options) ){
		.options <- .translate.string(.options, 
				c(t='track', v='verbose', d='debug'
				, p='parallel', P='parallel.required'
				, k='keep.all', r='restore.seed', f='dry.run'
				, g='garbage.collect'))
	}
	
	# setup verbosity options
	debug <- if( !is.null(.options$debug) ) .options$debug else nmf.getOption('debug')
	verbose <- if( debug ) Inf
				else if( !is.null(.options$verbose) ) .options$verbose
				else nmf.getOption('verbose')
	
	# nmf over a range of values: pass the call to nmfEstimateRank
	if( length(rank) > 1 ){
		if( verbose <= 1 )
			.options$verbose <- FALSE
		if( missing(nrun) )
			nrun <- 30
		return( nmfEstimateRank(x, range = rank, method = method, nrun = nrun
								, seed = seed, model = model
								, .pbackend = .pbackend, .callback = .callback
								, verbose=verbose, .options=.options, ...) )
	}
		
	# dry run
	dry.run <- if( is.null(.options$dry.run) ) FALSE else .options$dry.run 
	# call the garbage collector regularly
	opt.gc <- if( is.null(.options$garbage.collect) ) FALSE else .options$garbage.collect
	if( is.logical(opt.gc) && opt.gc )
		opt.gc <- ceiling(nrun / 3)
	
	# parallel runs options
	# backend
	# TODO: disable foreach if .pbackend=NULL
	if( is.null(.pbackend) ) .pbackend <- 'seq'
	else if( is.na(.pbackend) ) .pbackend <- 'registered'
	# option require-parallel: use value (potentially numeric) if not equivalent to FALSE
	opt.parallel.required <- if( !is.null(.options$parallel.required) && .options$parallel.required)
								.options$parallel.required
							else 
								FALSE
	# option require-parallel implies and takes precedence over option try-parallel
	if( opt.parallel.required ) .options$parallel <- opt.parallel.required
	
	# option try-parallel: .pbackend='seq' forces sequential
	opt.parallel <- if( !is.null(.options$parallel) ) .options$parallel # prioritary on anything else
					else .pbackend != '' # run in parallel only if the backend is defined
			
	keep.all <- if( !is.null(.options$keep.all) ) .options$keep.all else FALSE
	if( is.function(.callback) ){
		if( nrun==1 )
			warning("NMF::nmf - argument '.callback' is not used when performing a single NMF run [nrun=1].", call.=FALSE)
		else if( keep.all ) 
			warning("NMF::nmf - argument '.callback' is not used when option 'keep.all' is TRUE", call.=FALSE)
	}
	
	# Set debug/verbosity option just for the time of the run
	old.opt <- nmf.options(debug=debug, verbose=verbose);
	on.exit({nmf.options(old.opt)}, add=TRUE)
	
	# make sure rank is an integer
	rank <- as.integer(rank)
	if( length(rank) != 1 ) stop("NMF::nmf - invalid argument 'rank': must be a single numeric value")
	if( rank <= 1 ) stop("NMF::nmf - invalid argument 'rank': must be greater than 1")
	
	# option 'restore.seed' is deprecated
	if( !is.null(.options$restore.seed) )
		warning("NMF::nmf - Option 'restore.seed' is deprecated and discarded since version 0.5.99.")
	
	if( verbose ){
		if( dry.run ) message("*** dry-run ***")
		message("NMF algorithm: '", name(method), "'")
	}
	
	##START_MULTI_RUN
	# if the number of run is more than 1, then call itself recursively
	if( nrun > 1 )
	{
		if( verbose ) message("Multiple runs: ", nrun)
		
		# allow overriding some options passed to the call that performs each single run 
		pass.options <- .options
		
		if( opt.parallel ){
			if( verbose > 1 )
				message("# Setting up requested `foreach` environment: "
						, if( opt.parallel.required ) 'require-parallel' else 'try-parallel'
						,' [', .pbackend, ']')
			
			# check for 'foreach' package: required for parallel computation
			if( !require.quiet(foreach) ){
				if( opt.parallel.required )
					stop("NMF::nmf - the 'foreach' package is required to run NMF parallel computation"
							, call.=FALSE)
				else if( verbose > 1 ) message("# NOTE: NMF parallel computation disabled ['foreach' package is missing].")
			
				# disable parallel computation
				opt.parallel <- FALSE
			}
		}
		
		# if opt.parallel is TRUE: check and setup everything is there to run in parallel mode
		if( opt.parallel ){
			
			## 0. SETUP PARALLEL MODE
			worker.type <- ''
			single.machine <- TRUE
			ncores <- min(Inf, getOption('cores')) # by default use all cores, or rather use the 'cores' option if not NULL
			
			# get number of cores requested from options p or P and argument .pbackend
			if( .pbackend != 'seq'){
				if( !is.numeric(ncores) || ncores <= 0 )
					stop("NMF::nmf - invalid number of core(s) specified in option 'cores' [",ncores,"]", call.=FALSE)
				# get the number of workers/cores to use from the options or argument '.pbackend'
				if( is.numeric(opt.parallel) && length(opt.parallel) > 0 ){
					ncores <- opt.parallel[1]
					if( ncores <= 0 )
						stop("NMF::nmf - invalid number of core(s) specified in option 'p' (parallel) or 'P' (parellel.required) [",ncores,"]", call.=FALSE)
					else if( ncores == 1 )
						.pbackend <- 'seq'
					else 
						.pbackend <- 'mc'
				}
				else if( is.numeric(.pbackend) && length(.pbackend) > 0 ){
					ncores <- .pbackend[1]
					if( ncores <= 0 )
						stop("NMF::nmf - invalid number of core(s) specified in argument '.pbackend' [",ncores,"]", call.=FALSE)
					else if( ncores == 1 )
						.pbackend <- 'seq'
					else 
						.pbackend <- 'mc'
				}
			}
			
			## Check the host capability of using the required parallel backend 
			switch( .pbackend,
				mc = {
					if( verbose > 1 )
						message("# Setting up `doMC` ... ", appendLF=FALSE)
					# test the OS: multicore package does not work on Windows
					if( .Platform$OS.type == 'windows' ){
						# error only if the parallel computation was explicitly asked by the user
						if( opt.parallel.required ){
							if( verbose > 1 ) message("ERROR")
							stop('NMF::nmf - multicore computation impossible [not available under MS Windows]'
									, call.=FALSE)
						}else if( verbose > 1 ){
							message("NOTE")
							message("# NOTE: NMF parallel computation disabled [not availbale under MS Windows].")
						}
						
						# disable parallel computation
						opt.parallel <- FALSE
					}
					else if( is.Mac(check.gui=TRUE) ){ # check if we're not running on MAC from GUI
						# error only if the parallel computation was explicitly asked by the user
						if( opt.parallel.required ){
							if( verbose > 1 ) message("ERROR")
							stop("NMF::nmf - multicore computation stopped [not safe from R.app on Mac OS X]."
								, "\n\t-> Please use a terminal session, starting R from the command line."
								, call.=FALSE)
						}else if( verbose > 1 ){ 
							message("NOTE")
							message("# NOTE: NMF parallel computation disabled [not safe from R.app on Mac OS X]."
									, "\n\t-> To be able to use it, please use a terminal session, starting R from the command line.")
						}
						
						# disable parallel computation
						opt.parallel <- FALSE
					}
					else if( require.quiet(doMC) ){
						
						ncores.machine <- parallel::detectCores()
						if( ncores.machine == 1 ){
							if( opt.parallel.required ){
								if( verbose > 1 ) message("ERROR")
								stop("NMF::nmf - multicore computation aborted : single core detected"
									, call.=FALSE)
							}else if( verbose > 1 ){
								message("NOTE")
								message("# NOTE: NMF parallel computation disabled [single core detected]")
							}
							opt.parallel <- FALSE
						}else{
							use.ncores <- min(ncores.machine, ncores)
							# warn the user if a reduced number of cores is used
							if( !is.infinite(ncores) && use.ncores < ncores ){
								if( opt.parallel.required ){
									if( verbose > 1 ) message("ERROR")
									stop("NMF::nmf - multicore computation aborted : "
											, ncores," cores were required but only ", ncores.machine," cores were detected"
											, call.=FALSE)
								}else if( verbose > 1 ){
									message("NOTE")
									message("# NOTE: NMF parallel computation will only use "
											, use.ncores, "/",ncores.machine," cores [requested ", ncores, " cores]")
								}
							}
							else if( verbose > 1 ) 
								message("OK") 
							registerDoMC(use.ncores)
							
							worker.type <- paste('/', ncores.machine, ' core(s)', sep='')
						}
					}
					else if( opt.parallel.required ){
						if( verbose > 1 ) message("ERROR")
						stop("NMF::nmf - missing required package for multicore computation: 'doMC'"
								, call.=FALSE)
					}else{
						if( verbose > 1 ){
							message("NOTE")
							message("# NOTE: NMF multicore computation disabled [package 'doMC' is required]")
						}
						
						# disable parallel computation
						opt.parallel <- FALSE
					}
					
				}
				, seq = {
					if( verbose > 1 ) message("# Setting up `doSEQ` ... ", appendLF=FALSE)
					registerDoSEQ()
					if( verbose > 1 ) message("OK")
					worker.type <- ' core'					
				}
#				, registered = {
#					worker.type <- 'node(s)'
#					if( !getDoParRegistered() )
#						stop("NMF::nmf - no relse if( opt.parallel.required )egistered backend to run NMF parallel computation")
#				}
				, stop("NMF::nmf - invalid backend ['", .pbackend, "'] for NMF parallel computation. Argument '.pbackend' must be one of 'mc' (or number of cores), 'seq' (or NULL). See ?nmf"
						, call.=FALSE)
			)

			# From this point, the backend is registered
			# => one knows if we'll run a sequential or parallel foreach loop
			.MODE_SEQ <- getDoParName() == 'doSEQ'
			.MODE_PAR <- !.MODE_SEQ

			# if one wants to keep only the best result one needs the package 'bigmemory'
			# to store the best residual in shared memory
			# From version 4 of bigmemory, the mutexes are in package synchronicity
			if( opt.parallel && !keep.all ){
				if( verbose > 1 )
					message("# Setting up `bigmemory` ... ", appendLF=FALSE)
				if( require.quiet(bigmemory) ){
					if( verbose > 1 )
						message("OK")
					
					# check if the installed version of 'bigmemory' provides mutexes (prior to version 4.0)
					use.bigmemory4 <- !existsFunction('rw.mutex', where=asNamespace('bigmemory'))
					# otherwise one requires the 'synchronicity' package for non-sequential foreach loops
					if( .MODE_PAR && use.bigmemory4 ){
						
						if( verbose > 1 )
							message("# Setting up `synchronicity` ... ", appendLF=FALSE)

						if(  !require.quiet(synchronicity) ){
							if( opt.parallel.required ){
								if( verbose > 1 ) message("ERROR")
								stop("NMF::nmf - the 'synchronicity' package is required to run in parallel mode with option 'keep.all'=FALSE"
										, call.=FALSE)					
							}else{
								if( verbose > 1 ){
									message("NOTE")
									message("# NOTE: NMF multicore computation disabled [package 'synchronicity' is required when option 'keep.all'=FALSE]")
								}
						
								# disable parallel computation
								opt.parallel <- FALSE
							}
						}else if( verbose > 1 )
							message("OK")
					}
				}else if( opt.parallel.required ){
					if( verbose > 1 ) message("ERROR")
					stop("NMF::nmf - the 'bigmemory' package is required to run in parallel mode with option 'keep.all'=FALSE"
						, call.=FALSE)
				}else{
					if( verbose > 1 ){
						message("NOTE")
						message("# NOTE: NMF multicore computation disabled [package 'bigmemory' is required when option 'keep.all'=FALSE]")
					}

					# disable parallel computation
					opt.parallel <- FALSE
				}
			}
		}
		
		# check seed method: fixed values are not sensible -> warning
		if( is.nmf(seed) && !is.empty.nmf(seed) )
			warning("NMF::nmf - You are running multiple NMF runs with a fixed seed")
		# start_RNG_all
		# if the seed is numerical or a rstream object,	then use it to set the 
		# initial state of the random number generator:			
		# build a sequence of RNGstreams: if no suitable seed is provided
		# then the sequence use a random seed generated with a single draw 
		# of the current active RNG. If the seed is valid, then the 
		# 
		if( verbose > 2 ){
			message("# ** Original RNG settings:")
			.showRNG()
		}
		# setup the RNG sequence
		.RNG.seed <- .setupRNG(seed, n = nrun, verbose=verbose)
		stopifnot( length(.RNG.seed) == nrun )
		# store the current RNG state
		#NB: if no seeding occured then the RNG has still been drawn once in RNGseq
		mainRNG <- getRNG()
		# update the seeding method if necessary
		if( is.numeric(seed) )			
			seed <- 'random'
		# restore RNG settings on exit
		on.exit({
			if( verbose > 2 )
				message("# Restoring RNG settings ... ", appendLF=FALSE)							
			setRNG(mainRNG)					
			if( verbose > 2 ){
				message("OK")
				.showRNG()
			}
		}, add=TRUE)
		#end_RNG_all
		
		####FOREACH_NMF
		if( opt.parallel ){
			
			# in parallel mode: verbose message from each run are only shown in debug mode
			pass.options$verbose <- FALSE 
			
			if( verbose > 1 )
					message("# Using foreach backend: ", getDoParName()
							," [version ", getDoParVersion(),"]")
				
			run.all <- function(...){
								
				## 1. SETUP			
				# Specific thing only if one wants only the best result
				if( !keep.all ){ 
					# - Define the shared memory objects
					.SHARED.err <- if( use.bigmemory4 ) big.matrix(1, 1, type='double', init=NA) 
									else shared.big.matrix(1, 1, type='double', init=NA)					
					.SHARED.err.desc <- bigmemory::describe(.SHARED.err)				
					# the consensus matrix is computed only if not all the results are kept				
					.SHARED.consensus <- if( use.bigmemory4 ) big.matrix(ncol(x), ncol(x), type='double', init=0)
										else shared.big.matrix(ncol(x), ncol(x), type='double', init=0)
					.SHARED.consensus.desc <- bigmemory::describe(.SHARED.consensus)
			
					## In MODE_PAR: define mutex management functions to control 
					# access to the shared memory objects
					if( .MODE_PAR ){
						if( use.bigmemory4 ){
							if( verbose > 1 ) message("# Mutex provider: `synchronicity`")
							initMutex <- boost.mutex
							lockMutex <- function(mut.desc){
								mut <- attach.mutex(mut.desc)
								lock(mut)
								mut
							}
							unlockMutex <- unlock
							#boost.mutex()
						}else{
							if( verbose > 1 ) message("# Mutex provider: `bigmemory`")
							initMutex <- rw.mutex
							lockMutex <- function(mut.desc){
								mut <- attach.rw.mutex(mut.desc)	
								rwlock(mut)
								mut
							}
							unlockMutex <- unlock
							#rw.mutex()
						}
						
						# initialize mutex
						mut <- initMutex()
						mut.desc <- bigmemory::describe(mut)			
						##
					}else{
						mut.desc <- NULL
						lockMutex <- function(...){}
						unlockMutex <- function(...){}
					}
					#
				
					# - Define a temporary file to store the best fit				
					best.filename <- paste(tempfile('nmf.run.'), 'rds', sep='.')
				}
				
				## 2. RUN
				if( verbose ){
					if( debug || (.MODE_SEQ && verbose > 1) )
						pass.options$verbose <- verbose

					message("Mode: ",
							if( getDoParWorkers() == 1 ) "sequential [foreach]" 
									else paste("parallel (", getDoParWorkers(), worker.type,")", sep='')
					)
					
					if( !.MODE_SEQ && !debug || (.MODE_SEQ && verbose == 1) )
						cat("Runs:")
				}
				
				res.runs <- foreach(n=1:nrun
								, RNGobj = .RNG.seed
								, .verbose=debug, .errorhandling='pass'
								#, .options.RNG=.RNG.seed
								) %dopar% { #START_FOREACH_LOOP
				
					# in mode sequential or debug: show details for each run
					if( .MODE_SEQ && verbose > 1 )
						cat("\n## Run: ",n, "/", nrun, "\n", sep='')
					
					# set the RNG if necessary and restore after each run 
#					if( !is.null(RNGobj) ){
						if( .MODE_SEQ && verbose > 2 )
							message("# Setting up loop RNG ... ", appendLF=FALSE)
						
						# unpack and set the RNG stream  
#						rstream.packed(RNGobj) <- FALSE
#						oldRNG <- .tryCatch_setRNG(RNGobj, verbose=verbose>3 && .MODE_SEQ)
						setRNG(RNGobj, verbose=verbose>3 && .MODE_SEQ)
						
						if( .MODE_SEQ && verbose > 2 )
							message("OK")
						
#						# setup restoration
#						on.exit({
#							if( .MODE_SEQ && verbose > 2 )
#								message("# Restoring from loop RNG settings ... ", appendLF=FALSE)
#							setRNG(oldRNG)
#							if( .MODE_SEQ && verbose > 2 ){
#								message("OK")
#								.showRNG()
#							}
#						}, add=TRUE)
#					}
									
					# limited verbosity in simple mode
					if( verbose && !(.MODE_SEQ && verbose > 1)){
						mut <- lockMutex(mut.desc)
						cat('', n)		
						unlockMutex(mut)
					}
					
					# fit a single NMF model
					res <- nmf(x, rank, method, nrun=1, seed=seed, .options=pass.options, ...)
					
					# if only the best fit must be kept then update the shared objects
					if( !keep.all ){
						# load shared objects
						.SHARED.err <- attach.big.matrix(.SHARED.err.desc)
						.SHARED.consensus <- attach.big.matrix(.SHARED.consensus.desc)
						
						if( .MODE_PAR ){
							##LOCK_MUTEX					
							mut <- lockMutex(mut.desc)
#							# retrieve and lock the mutex
#							if( use.bigmemory4 ){
#								mut <- attach.mutex(mut.desc)
#								lock(mut) 
#							}else{
#								mut <- attach.rw.mutex(mut.desc)	
#								rwlock(mut)
#							}
						}
						
						# check if the run found a better fit
						.STATIC.err <- .SHARED.err[]
						
						# retrieve the residual error TODO: call deviance here
						err <- residuals(res)
						
						if( is.na(.STATIC.err) || err < .STATIC.err ){
							
							if( n>1 && verbose ){
								if( .MODE_SEQ && verbose > 1 ) cat("## Better fit found [err=", err, "]\n")
								else cat('*')
							}
							
							# update residuals
							.SHARED.err[] <- err							
							# update best fit
							saveRDS(res, file=best.filename, compress=FALSE)
															
						}
						
						# update the consensus matrix
						.SHARED.consensus[] <- .SHARED.consensus[] + connectivity(res)[]
						
						# call the callback function if necessary
						if( is.function(.callback) ){
							cb.res <- tryCatch(.callback(res), error=function(e) e)
							if( is(cb.res, 'error') ){
								class(cb.res) <- c(class(cb.res), 'errorCB')
							}
						}
						
						if( .MODE_PAR ){
							# unlock the mutex
							unlockMutex(mut)
							##END_LOCK_MUTEX
						}
									
						# reset the result to NULL
						res <- NULL
						
						# if a callback function was provided, actually return the result of the callback
						if( is.function(.callback) )
							res <- cb.res
					}
					
					# garbage collection if requested
					if( opt.gc && n %% opt.gc == 0 ){
						if( verbose > 2 ){
							if( .MODE_SEQ )
								message("# Call garbage collector")
							else{
								mut <- lockMutex(mut.desc)
								cat('%')
								unlockMutex(mut)
							}
						}
						
						gc(verbose= .MODE_SEQ && verbose > 3)
					}
					
					# return the result
					res
				}				
				## END_FOREACH_LOOP
				
				## 3. CHECK FOR ERRORS
				# - in the run
				errors <- sapply(res.runs, function(x){ is(x, 'error') && !is(x, 'errorCB') })				
				nerrors <- sum(errors)
				if( nerrors > 0 ){										
					stop("NMF::nmf - ", nerrors,"/", nrun, " fit(s) threw an error.\n"
						,"# Error(s):\n-- "
						, paste(sapply(unique(res.runs[which(errors)]), function(x) x$message), collapse="\n-- ")
						, call.=FALSE)
				}
				# - in the callback
				if( is.function(.callback) ){
					errors <- sapply(res.runs, function(x) is(x, 'errorCB'))				
					nerrors <- nerrors.cb <- sum(errors)
					if( nerrors > 0 ){										
						warning("NMF::nmf - all NMF fits were successful but ", nerrors,"/", nrun, " callback call(s) threw an error.\n"
								,"# ", if(nerrors>10) "First 10 c" else "C", "allback error(s):\n-- "
								, paste(paste("Run #", 1:min(nerrors,10),': ', sapply(res.runs[which(errors)[1:min(nerrors,10)]], function(x) x$message), sep=''), collapse="\n-- ")
								, call.=FALSE)
					}
				}	
				## 4. WRAP UP
				res <- list(fit=res.runs)
				
				if( !keep.all ){
					# check existence of the result file
					if( !file_test('-f', best.filename) )
						stop("NMF::nmf - error in parallel mode: the result file does not exist")
					res <- readRDS(best.filename)					
					# NB: the object 'res' is now the best NMFfit
					#remove the result file
					unlink(best.filename)
					
					# wrap the result in a list: fit + consensus
					res <- list(fit=res, consensus=.SHARED.consensus[])
					
					# add the result of the callback function if necessary
					if( is.function(.callback) )
						res$.callback <- if( nerrors.cb > 0 ) res.runs else sapply(res.runs, identity)					
					
				}
				##
				
				if( .MODE_SEQ && verbose>1 ) cat("## DONE\n")
				
				# return result
				res
			}			
		}####END_FOREACH_NMF
		else{####SAPPLY_NMF
			
			# by default force no verbosity from the runs
			pass.options$verbose=FALSE
			if( verbose ){
				message("Mode: sequential [sapply]")
				if( verbose > 1 ){					
					# pass verbosity options in this case
					pass.options$verbose <- verbose
				}
			}
			
			run.all <- function(...){
				
				## 1. SETUP				
				# define static variables for the case one only wants the best result
				if( !keep.all ){
					# statis list with best result: fit, residual, consensus
					best.static <- list(fit=NULL, residuals=NA, consensus=matrix(0, ncol(x), ncol(x)))					
				}
											
				## 2. RUN:
				# perform a single run `nrun` times								
				if( verbose && !debug ) cat('Runs:')
				res.runs <- lapply(1:nrun, function(n){
					
					#start_verbose
					if( verbose ){
						# in mode verbose > 1: show details for each run
						if( verbose > 1 ){
							cat("\n## Run: ",n, "/", nrun, "\n", sep='')							
						}else{
						# otherwise only some details for the first run
							if(n == 1 ){
								.showRNG()								
								cat("Runs:")
							}
							cat('', n)
						}
					}#end_verbose
					
					# set the RNG if necessary and restore after each run
					RNGobj <- .RNG.seed[[n]]
#					if( !is.null(RNGobj) ){
						if( verbose > 2 )
							message("# Setting up loop RNG ... ", appendLF=FALSE)
						
						# set the RNG stream  						
#						oldRNG <- .tryCatch_setRNG(RNGobj, verbose=verbose>3)
						setRNG(RNGobj, verbose=verbose>3)
						
						if( verbose > 2 )
							message("OK")
						
#						# TODO: ONLY RESTORE IF NEXT RNG IS PROVIDED!!!
#						# setup restoration
#						on.exit({
#							if( verbose > 2 )
#								message("# Restoring from loop RNG settings ... ", appendLF=FALSE)
#							setRNG(oldRNG)
#							if( verbose > 2 ){
#								message("OK")
#								.showRNG()
#							}
#						}, add=TRUE)
#					}
				
					# fit a single NMF model
					res <- nmf(x, rank, method, nrun=1, seed=seed, .options=pass.options, ...)
					
					if( !keep.all ){						
						# check if the run found a better fit
						err <- residuals(res)
						best <- best.static$residuals
						if( is.na(best) || err < best ){
							if( n>1 && verbose ){
								if( debug ) cat("## Better fit found [err=", err, "]\n")
								else cat('*')
							}
							
							# update best fit (only if necessary)
							best.static$fit <<- res

							best.static$residuals <<- err				
						}
							
						# update the static consensus matrix (only if necessary)
						best.static$consensus <<- best.static$consensus + connectivity(res)
						
						# call the callback function if necessary
						if( is.function(.callback) ){
							res <- tryCatch(.callback(res), error=function(e) e)
							if( is(res, 'error') ){								
								class(res) <- c(class(res), 'errorCB')
							}							
						}
						else # reset the result to NULL
							res <- NULL
						
					}
					
					# garbage collection if requested
					if( opt.gc && n %% opt.gc == 0 ){
						if( verbose > 1 )
							message("# Call garbage collection NOW")
						else if( verbose )
							cat('%')
						
						gc(verbose= .MODE_SEQ && verbose > 3)
					}
					
					if( verbose > 1 ) cat("## DONE\n")
					
					# return the result
					res
				})
				##
				
				## 3. CHECK FOR ERRORS
				# - in the run
				errors <- sapply(res.runs, function(x){ is(x, 'error') && !is(x, 'errorCB') })				
				nerrors <- sum(errors)
				if( nerrors > 0 ){										
					stop("NMF::nmf - ", nerrors,"/", nrun, " fit(s) threw an error.\n"
							,"# Error(s):\n-- "
							, paste(sapply(unique(res.runs[which(errors)]), function(x) x$message), collapse="\n-- ")
							, call.=FALSE)
				}
				# - in the callback
				if( is.function(.callback) ){
					errors <- sapply(res.runs, function(x) is(x, 'errorCB'))				
					nerrors <- nerrors.cb <- sum(errors)
					if( nerrors > 0 ){										
						warning("NMF::nmf - all NMF fits were successful but ", nerrors,"/", nrun, " callback call(s) threw an error.\n"
								,"# ", if(nerrors>10) "First 10 c" else "C", "allback error(s):\n-- "
								, paste(paste("Run #", 1:min(nerrors,10),': ', sapply(res.runs[which(errors)[1:min(nerrors,10)]], function(x) x$message), sep=''), collapse="\n-- ")
								, call.=FALSE)
					}
				}
				
				## 4. WRAP UP
				res <- list(fit=res.runs)
				
				if( !keep.all ){
					res$fit <- best.static$fit
					res$consensus <- best.static$consensus
					# add the result of the callback function if necessary
					if( is.function(.callback) )
						res$.callback <- if( nerrors.cb > 0 ) res.runs else sapply(res.runs, identity)
				}
				##
				
				# return the result
				res
			}
			
		}####END_SAPPLY_NMF
			
		####END_DEFINE_RUN		
		
		# perform all the NMF runs
		t <- system.time({res <- run.all(...)})
		if( verbose && !debug ){
			cat(" ... DONE\n")
			cat("System time:\n")
			print(t)
		}
		
		# ASSERT the presence of the result
		stopifnot( !is.null(res$fit) )
		
		# if one just want the best result only return the best 
		if( !keep.all ){
			# ASSERT the presence of the consensus matrix
			stopifnot( !is.null(res$consensus) )
			res.final <- join(res$fit, consensus=res$consensus/nrun, runtime.all=t, nrun=as.integer(nrun), rng1=.RNG.seed[[1]])
			
			# set the callback result if necessary
			if( is.function(.callback) )
				res.final$.callback <- res$.callback
			
			return(res.final)
		}
		
		# when keeping all the fits: join the results into an NMFfitXn object
		# TODO: improve memory management here
		return( join(res$fit, runtime.all=t) )
		
	}##END_MULTI_RUN
	
	# start_RNG
	# show original RNG settings in verbose > 2
	if( verbose > 2 ){
		message("# ** Current RNG settings:")
		.showRNG()
	}
	
	# do something if the RNG was actually changed
	mainRNG <- getRNG()
	.RNG.seed <- .setupRNG(seed, 1, verbose=verbose)
	if( verbose > 2 ) .showRNG()
	
	# update the seeding method
	if( !is.null(.RNG.seed) ){
		seed <- 'random'
		# restore RNG settings
		on.exit({
			if( verbose > 2 )
				message("# Restoring original RNG settings ... ", appendLF=FALSE)							
			setRNG(mainRNG)					
			if( verbose > 2 ){
				message("OK")
				.showRNG()
			}
		}, add=TRUE)
	}
	#end_RNG
	
	# CHECK PARAMETERS:	
	# test for negative values in x only if the method is not mixed
	if( !is.mixed(method) && min(x) < 0 ) stop('Input matrix ', substitute(x),' contains some negative entries.');
	# test if one row contains only zero entries
	if( min(rowSums(x)) == 0) stop('Input matrix ', substitute(x),' contains at least one null row.');	

	# a priori the parameters for the run are all the one in '...'
	parameters.method <- list(...)
	
	if( inherits(seed, 'NMF') ){
		# if the seed is a NMFfit object then only use the fit (i.e. the NMF model)
		# => we want a fresh and clean NMFfit object
		if( isNMFfit(seed) )
			seed <- fit(seed)
		
		# Wrap up the seed into a NMFfit object
		seed <- newNMFfit(fit=seed, seed='none')
	}
	else if( !inherits(seed, 'NMFfit') ){
		
		# retrieve the NMF model to use (from the provided method)
		init <- modelname(method)
		stopifnot( extends(init, 'NMF') )
		
		## MODEL INSTANTIATION :
	
		# some of the instantiation parameters are set internally
		# TODO: change target into x (=> impact on nmfModel
		parameters.model.internal <- list(rank=rank, target=0, model=init)
		parameters.model <- list()
		
		# if 'model' is NULL: initialization parameters are searched in '...' 
		if( is.null(model) ){
			
			# extract the parameters from '...' that correspond to slots in the given class
			parameters <- .extract.slots.parameters(init, ...)	
					
			# restrict parameters.method to the ones that won't be used to instantiate the model
			overriden <- is.element(names(parameters$slots), names(parameters.model.internal))
			parameters.method <- c(parameters$extra, parameters$slots[overriden])
			
			#- the model parameters come from the remaining elements
			parameters.model <- parameters$slots
			
		}else if( is.list(model) ){  # otherwise argument 'model' must be a list
			
			# if the list is not empty then check all elements are named and 
			# not conflicting with the internally set values
			if( length(model) > 0 ){
				# all the elements must be named
				if( is.null(names(model)) || any(names(model)=='') )  
					stop("NMF::nmf : Invalid argument 'model' [elements must all be named]. See ?nmf."
					, call.=FALSE)
				
				# warn the user if some elements are conflicting and won't be used
				overriden <- is.element(names(model), names(parameters.model.internal))
				if( any(overriden) )
					warning("NMF::nmf : Model parameter(s) " 
							, paste( paste("'", names(model)[overriden], "'", sep=''), collapse=', ')
							, " discarded. Internal values are used instead."
							, call.=FALSE)
			}
			
			# all the instanciation parameters come from argument 'model'
			parameters.model <- model
			
		}else{ 			
			stop("NMF::nmf : Invalid argument 'model' [expected NULL or a list to set slots in the NMF model class '",init,"']. See ?nmf."
					, call.=FALSE)
		}	
			
			
		#- force the value of the internally set arguments for the instantiation of the model
		init <- .merge.override(parameters.model, parameters.model.internal)		
				
	# at this point 'init' should be the list of the initialization parameters
	if( !is.list(init) ) stop("Invalid object: 'init' must be a list")
	if( !is.element('model', names(init)) ) stop("Invalid object: 'init' must contain an element named 'model'")	
	
	## SEEDING:
	# the seed must either be an instance of class 'NMF', the name of a seeding method as a character string
	# or a list of parameters to pass to the 'seed' function.
			parameters.seed <- list()
			seed.method <- NULL
			if( (is.character(seed) && length(seed) == 1) 
				|| is.numeric(seed) 
				|| is.null(seed) 
#				|| is(seed, 'rstream') 
				) seed.method <- seed
			else if( is.function(seed) ) seed.method <- seed
			else if( is.list(seed) ){ # seed is a list...
				
				if( !is.null(seed$method) ){ # 'seed' must contain an element giving the method...
					seed.method <- seed$method
					parameters.seed <- seed[-which(names(seed)=='method')]
				}
				else if ( is.null(names(seed)) || names(seed)[1] == '' ){ # ... or the first element must be a method
					seed.method <- seed[[1]]
					if( length(seed) > 1 ) parameters.seed <- seed[2:length(seed)]
				}
				else stop("Invalid parameter: list 'seed' must contain the seeding method through its first element or through an element named 'method'")
				
				# check validity of the method provided via the list
				if( !is.function(seed.method) && !(is.character(seed.method) && length(seed.method)==1) )
					stop("The seeding method provided by parameter 'seed' is invalid: a valid function or a character string is expected")
			}
			else stop("Invalid parameter 'seed'. Acceptable values are:\n\t- ",
						paste("an object that inherits from class 'NMF'"
							, "the name of a seeding method (see ?nmfSeed)"
							, "a valid seed method definition"
							, "a list containing the seeding method (i.e. a function or a character string) as its first element\n\tor as an element named 'method' [and optionnally extra arguments it will be called with]"
							, "a numerical value used to set the seed of the random generator"
							, "NULL to directly pass the model instanciated from arguments 'model' or '...'."
							, sep="\n\t- "))
						 			
			# call the 'seed' function passing the necessary parameters
			if( verbose )
				message("NMF seeding method: ", 
						if( is.character(seed.method) || is.numeric(seed.method) ) seed.method
						else if( is.null(seed.method) ) 'NULL'
						else if( !is.null(attr(seed.method, 'name')) ) attr(seed.method, 'name') 
						else if( is.function(seed.method) ) '<function>'
						else NA)
			
			#seed <- do.call(getGeneric('seed', package='NMF')
			seed <- do.call(getGeneric('seed')
					, c(list(x=x, model=init, method=seed.method), parameters.seed))
			
			# check the validity of the seed
			if( !inherits(seed, 'NMFfit') ) 
				stop("The seeding method function should return class 'NMF' ["
					, if( is.character(seed.method) ) paste('method "', seed.method, "' ", sep='') else NULL 
					, "returned class: '", class(seed), "']")
	}
	# -> at this point the 'seed' object is an instance of class 'NMFfit'
	nmf.debug('nmf', "Seed is of class: '", class(seed), "'")
	# ASSERT just to be sure
	if( !inherits(seed, 'NMFfit') )
		stop("NMF::nmf - Invalid class '", class(seed), "' for the computed seed: object that inherits from class 'NMFfit' expected.")
	
	# check the consistency of the NMF model expected by the algorithm and 
	# the one defined by the seed
	#if( none( sapply(model(method), function(c) extends(model(seed), c)) ) )
	if( all( !inherits(fit(seed), modelname(method)) ) )
		stop("NMF::nmf - Invalid NMF model '", modelname(seed),"': algorithm '", name(method), "' expects model(s) "
			, paste(paste("'", modelname(method),"'", sep=''), collapse=', ')
			, " or extension."
			, call.=FALSE)
	
	# get the complete seeding method's name 
	seed.method <- seeding(seed)
	
	## FINISH SETUP OF THE SEED OBJECT: store some data within the seed so
	# that strategy methods can access them directly
	algorithm(seed) <- name(method) # algorithm name
	seed@distance <- objective(method) # distance name
	seed@parameters <- parameters.method # extra parameters
	run.options(seed) <- nmf.options.runtime() # set default run options
	if( !is.null(.options$track) ) run.options(seed, 'error.track') <- .options$track
	run.options(seed, 'verbose') <- verbose

	## print options if in verbose > 3
	if( verbose >= 3 ){
		cat("## OPTIONS:\n")		
		sapply(seq_along(.options)
				, function(i){
					r <- i %% 4
					cat(if(r!=1) '\t| ' else "# ", names(.options)[i],': ', .options[[i]], sep='')
					if(r==0) cat("\n")
				})
		if( length(.options) %% 4 != 0 )cat("\n")
	}
	## RUN NMF METHOD:
	# call the strategy's run method [and time it] using the element of 'parameters.method' as parameters	
	parameters.run <- c(list(method=method, x=x, seed=seed), parameters.method)
	t <- system.time({				
		res <- if( !dry.run )
					do.call('run', parameters.run)
				else
					seed
	})
	
	## CHECK RESULT
	# check the result is of the right type
	if( !inherits(res, 'NMFfit') ) stop("NMF method should return an instance of class 'NMFfit' [returned class:", class(res), "]")			

	## ENSURE SOME SLOTS ARE STILL CORRECTLY SET
	# slot 'method'
	algorithm(res) <- name(method)	
	# slot 'distance'
	res@distance <- objective(method)	
	# slot 'seed'
	if( seed.method != '' ) seeding(res) <- seed.method
	# slot 'parameters'
	res@parameters <- parameters.method
	# set dimnames of the result only if necessary
	if( is.null(dimnames(res)) )
		dimnames(res) <- dimnames(seed)
	
	## CLEAN-UP + EXTRAS:
	# add extra information to the object	 
	if( length(residuals(res)) == 0 ) residuals(res) <- objective(method, x, res)
	if( is.na(residuals(res)) ) warning("NMF residuals: final objective value is NA")
	res@runtime <- t
	
	# return the result
	res
})

###% Common interface for seeding methods for Nonnegative Matrix Factorization (NMF) algorithms.
###% 
###% This function calls the different seeding methods that define a starting point for NMF methods.
###% These methods at least set the slots \code{W} and \code{H} of \code{object} to valid nonnegative matrices.
###% They will be used as a starting point by any NMF algorithm that accept initialization.
###%
###% @param x The target matrix one wants to approximate with NMF
###% @param rank The rank of the factorization to seed
###% @param method Name of the seeding method to call
###% @param object Either an instance (resp. the name) of a NMF (sub)class, that will be seeded (resp. instanciated and seeded).
###% @param ... Parameters to be passed to the seeding method
###%  
###% @return \code{object} initialized using method \code{name}.
###%  
setGeneric('seed', function(x, model, method, ...) standardGeneric('seed') )
setMethod('seed', signature(x='ANY', model='ANY', method='missing'),
	function(x, model, method, ...){
					
		seed(x, model, nmf.getOption('default.seed'), ...)
		
	}
)
setMethod('seed', signature(x='ANY', model='ANY', method='NULL'),
	function(x, model, method, ...){
		
		seed(x, model, 'none', ...)

	}
)
setMethod('seed', signature(x='ANY', model='ANY', method='numeric'),
		function(x, model, method, ...){
			
			# set the seed using the numerical value by argument 'method'
			orng <- setRNG(method)
			#TODO: restore the RNG state?
			
			# call seeding method 'random'
			res <- seed(x, model, 'random', ...)
			
			# return result
			return(res)
		}
)
#setMethod('seed', signature(x='ANY', model='ANY', method='rstream'),
#		function(x, model, method, ...){
#			
#			# set the seed using the numerical value by argument 'method'
#			orng <- setRNG(method)
#			#TODO: restore the RNG state? 
#			
#			# call seeding method 'random'
#			res <- seed(x, model, 'random', ...)
#			
#			# return result
#			return(res)
#		}
#)
setMethod('seed', signature(x='ANY', model='ANY', method='character'),
		function(x, model, method, ...){
			
			# get the seeding method from the registry
			seeding.fun <- nmfSeed(method)
				
			#Vc#Use seeding method: '${method}'
			# call 'seed' with the seeding.function			
			seed(x, model, method=seeding.fun, ...)
			
		}
)
setMethod('seed', signature(x='ANY', model='list', method='NMFSeed'), 
	function(x, model, method, ...){	
		
		## check validity of the list: there should be at least the NMF (sub)class name and the rank
		if( length(model) < 2 )
			stop("Invalid parameter: list 'model' must contain at least two elements giving the model's class name and the factorization rank")
					
		# 'model' must contain an element giving the class to instanciate
		if( is.null(model$model) ){
			
			err.msg <- "Invalid parameter: list 'model' must contain a valid NMF model classname in an element named 'model' or in its first un-named element"			
			unamed <- if( !is.null(names(model)) ) which(names(model) %in% c('', NA)) else 1			
			if ( length(unamed) > 0 ){ # if not the first unamed element is taken as the class name
				idx <- unamed[1]
				val <- unlist(model[idx], rec=FALSE)				
				if( is.character(val) && length(val)==1 && extends(val, 'NMF') )
					names(model)[idx] <- 'model'
				else stop(err.msg)
			}else stop(err.msg)
		}
		
		# 'model' must contain an element giving the factorization rank
		if( is.null(model$rank) ){
			err.msg <- "Invalid parameter: list 'model' must contain the factorization rank in an element named 'rank' or in its second un-named element"
			unamed <- if( !is.null(names(model)) ) which(names(model) %in% c('', NA)) else 1
			if ( length(unamed) > 0 ){ # if not the second element is taken as the factorization rank
				idx <- unamed[1]
				val <- unlist(model[idx], rec=FALSE)
				if( is.numeric(val) && length(val)==1 )
					names(model)[idx] <- 'rank'
				else stop(err.msg)
			}
			else stop(err.msg)
		}
					
		s.init <- capture.output(print(model))
		nmf.debug('seed', "using model parameters:\n", s.init)			
		# instantiate the object using the factory method		
		model <- do.call('nmfModel', model)
		nmf.debug('seed', "using NMF model '", class(model), "'")
		
		# check that model is from the right type, i.e. inherits from class NMF
		if( !inherits(model, 'NMF') ) stop("Invalid object returned by model: object must inherit from class 'NMF'")
		
		seed(x, model, method, ...)
	}
)

setMethod('seed', signature(x='ANY', model='numeric', method='NMFSeed'), 
	function(x, model, method, ...){	

		seed(x, nmfModel(model), method, ...)
	}
)

setMethod('seed', signature(x='ANY', model='ANY', method='function'),
	function(x, model, method, name, ...){
		
		# generate runtime name if necessary
		if( missing(name) ) name <- basename(tempfile("NMF.seed."))
		# check that the name is not a registered name		
		if( existsNMFSeed(name) )
			stop("Invalid name for custom seeding method: '",name,"' is already a registered seeding method")
		
		# wrap function method into a new NMFSeed object		 						
		seedObj <- new('NMFSeed', name=name, method=method)
		# call version with NMFSeed 
		seed(x, model, seedObj, ...)
	}
)

setMethod('seed', signature(x='matrix', model='NMF', method='NMFSeed'), 
	function(x, model, method, rng, ...){	
				
		# debug message
		nmf.debug('seed', "use seeding method: '", name(method), "'")
				
		# temporarly set the RNG if provided 
		if( !missing(rng) ){
			orng <- setRNG(rng)
			on.exit(setRNG(orng))
		}
		
		# save the current RNG numerical seed
		rng.s <- getRNG()
		# create the result NMFfit object, storing the RNG numerical seed
		res <- newNMFfit()
		# call the seeding function passing the extra parameters, and store the result into slot 'fit'
		fit(res) <- do.call(method(method), c(list(model, x), ...))				
		# ASSERT: check that the RNG seed is correctly set
		stopifnot( rng.equal(res,rng.s) )
		
		# if not already set: store the seeding method's name in the resulting object
		if( seeding(res) == '' ) seeding(res) <- name(method)
		# set the dimnames from the target matrix
		dimnames(res) <- dimnames(x)
		
		# return the seeded object
		res
	}
)

###% Extract from a list the elements that can be used to initialize the slot of a class.
###% 
###% This function only extract named elements.
###% 
###% @param class.name Name of the class from whose slots will be search into '...'
###% @param ... The parameters in which the slot names will be search for
###% 
###% @return a list with two elements:
###% - \code{slots}: is a list that contains the named parameters that can be used to instantiate an object of class \code{class.name} 
###% - \code{extra}: is a list of the remaining parameters from \code{parameters} (i.e. the ones that do not correspond to a slot).
###%  
.extract.slots.parameters <- function(class.name, ...){
		
	# check validity of class.name
	if( !isClass(class.name) ) stop("Invalid class name: class '", class.name, "' dose not exist")
	
	# transform '...' into a list
	parameters <- list(...)	
	# get the slots from the class name
	slots <- slotNames(class.name)
	# get the named parameters that correspond to a slot
	in.slots <- is.element(names(parameters), slots)
	# return the two lists	
	list( slots=parameters[in.slots], extra=parameters[!in.slots])
}

###% Merges two lists, but overriding with the values of the second list in the case
###% of duplicates.
.merge.override <- function(l1, l2, warning=FALSE){
	sapply(names(l2), function(name){
				if( warning && !is.null(l1[[name]]) )
					warning("overriding element '", name, "'")
				l1[[name]] <<- l2[[name]]
			})
	
	# return updated list
	return(l1)
}

###% Estimate the factorization rank
nmfEstimateRank <- function(x, range, method=nmf.getOption('default.algorithm')
					, nrun=30, verbose=FALSE, stop=FALSE, ...){
	
	# check if there is no duplicates in the range
	if( any(duplicated(range)) )
		stop("duplicated values in argument 'range' are not allowed")	
	
	# initiate the list of consensus matrices: start with single NA values
	c.matrices <- setNames(lapply(range, function(x) NA), as.character(range))
	fit <- setNames(lapply(range, function(x) NA), as.character(range))
	bootstrap.measures <- list()

	# combine function: take all the results at once and merge them into a big matrix
	comb <- function(...){
		measures <- list(...)
		
		err <- which( sapply(measures, is.character) )		
		if( length(err) == length(measures) ){ # all runs produced an error
		
			# build an warning using the error messages
			msg <- paste(paste("#", seq_along(range),' ', measures, sep=''), collapse="\n\t-")
			stop("All the runs produced an error:\n\t-", msg)
		
		}else if( length(err) > 0 ){ # some of the runs returned an error
			
			# simplify the results with no errors into a matrix
			measures.ok <- sapply(measures[-err], function(x) x)
			
			# build a NA matrix for all the results
			n <- nrow(measures.ok)
			tmp.res <- matrix(NA, n, length(range))
			rownames(tmp.res) <- rownames(measures.ok)
			
			# set the results that are ok
			tmp.res[,-err] <- measures.ok
			# set only the rank for the error results 
			tmp.res['rank', err] <- range[err]
			# build an warning using the error messages
			msg <- paste(paste("#", err, measures[err], ' ', sep=''), collapse="\n\t-")
			warning("NAs were produced due to errors in some of the runs:\n\t-", msg)
			
			# return full matrix
			tmp.res
		}
		else # all the runs are ok 
			sapply(measures, function(x) x)
	}
	
#	measures <- foreach(r = range, .combine=comb, .multicombine=TRUE, .errorhandling='stop') %do% {
	measures <- sapply(range, function(r, ...){
			if( verbose ) cat("Compute NMF rank=", r, " ... ")
			
			res <- tryCatch({ #START_TRY
				
				res <- nmf(x, r, method, nrun=nrun, ...)					
				# directly return the result if a valid NMF result
				if( !isNMFResult(res) )
					return(res)
				
				# store the consensus matrix
				c.matrices[[as.character(r)]] <<- consensus(res)
				# store the fit
				fit[[as.character(r)]] <<- res				
				
				# if confidence intervals must be computed then do it
	#			if( conf.interval ){
	#				# resample the tries
	#				samp <- sapply(seq(5*nrun), function(i){ sample(nrun, nrun, replace=TRUE) })
	#				
	#				bootstrap.measures[[as.character(r)]] <<- apply(samp, 2, function(s){
	#					res.sample <- join(res[s])
	#					summary(res.sample, target=x)
	#				})
	#			}
				
				# compute quality measures
				if( verbose ) cat('+ measures ... ')
				measures <- summary(res, target=x)
				
				if( verbose ) cat("OK\n")
				
				# return the measures
				measures
			} #END_TRY
	
			, error = function(e) {
					mess <- if( is.null(e$call) ) e$message else paste(e$message, " [in call to '", e$call[1],"']", sep='')
					mess <- paste('[r=', r, '] -> ', mess, sep='')
					if( stop ){ # throw the error
						if( verbose ) cat("\n")
						stop(mess, call.=FALSE)
					} # pass the error message
					if( verbose ) message("ERROR")					
					return(mess)
				}
			)
			
			# return the result
			res
		}
	, ..., simplify=FALSE)
	
	measures <- do.call(comb, measures)
	
	# reformat the result into a data.frame
	measures <- as.data.frame(t(measures))	
	
	# wrap-up result into a 'NMF.rank' S3 object
	res <- list(measures=measures, consensus=c.matrices, fit=fit)
	#if( conf.interval ) res$bootstrap.measure <- bootstrap.measures
	class(res) <- 'NMF.rank'
	return(res)
	
}

plot.NMF.rank <- function(x, what=c('all', 'cophenetic', 'rss', 'residuals'
									, 'dispersion', 'evar', 'sparseness'
									, 'sparseness.basis', 'sparseness.coef')
						, ref=NULL, na.rm=FALSE, ... ){

	
	# little trick not to display a title when it is not needed
	if( !exists('.NMF.rank.plot.notitle', parent.frame()) ){
		.NMF.rank.plot.notitle <- FALSE
	}
	else .NMF.rank.plot.notitle <- TRUE
	
	what <- match.arg(what)	
	if( what == 'all' ){
		opar <- par(mfrow=c(2,3), oma=c(0,0,3,0))
		on.exit( par(opar), add=TRUE)
		sapply(c('cophenetic', 'rss', 'residuals'
				, 'dispersion', 'evar', 'sparseness'),
				function(w, ...){
					.NMF.rank.plot.notitle <- TRUE
					plot(x, w, ref, na.rm, ...) 
				}
		)
		title("NMF rank estimation", outer=TRUE)
		return(invisible())
	}
	
	measures <- x$measures
	iwhat <- grep(paste('^',what,sep=''), colnames(measures))
	
	# remove NA values if required
	if( na.rm )
		measures <- measures[ apply(measures, 1, function(row) !any(is.na(row[iwhat]))), ]
	
	vals <- measures[,iwhat, drop=FALSE]
	x <- as.numeric(measures$rank)
	xlim <- range(x)
	
	vals.ref <- NULL
	xref <- NULL
	if( !missing(ref) && is(ref, 'NMF.rank') ){
		
		# remove NA values if required
		if( na.rm )
			ref$measures <- ref$measures[ apply(ref$measures, 1, function(row) !any(is.na(row[iwhat]))), ]
			
		xref <- as.numeric(ref$measures$rank)
		xlim <- range(xlim, xref)
		vals.ref <- ref$measures[,iwhat, drop=FALSE]
	}
	
	# compute the ylim from main and ref values
	ylim <- range(vals, na.rm=TRUE)
	if( !is.null(vals.ref) )
		ylim <- range(ylim, vals.ref, na.rm=TRUE)
	
	# detect if the values should be between 0 and 1
	#if( all(ylim >=0 & ylim <= 1) )
	#	ylim <- c(0,1)

	# retreive the graphical parameters and match them to the sub-sequent call to 'plot.default'
	graphical.params <- list(...)
	names(graphical.params) <- .match.call.args(names(graphical.params), 'plot.default', call='NMF::plot.NMF.rank')

	# set default graphical parameters for type 'consensus'
	graphical.params <- .set.list.defaults(graphical.params
			, main=paste(if(!.NMF.rank.plot.notitle) "NMF rank estimation\n", "-", what, '-')
			, ylab=paste('Quality measure:', what)
			, xlab='Factorization rank'
			, xlim=xlim, ylim = ylim
	)
	
	#init the plot
	do.call(plot, c(list(x=NULL, axes=FALSE, frame=TRUE), graphical.params))
	
	opar <- par(lwd=2)
	on.exit( par(opar), add=TRUE)
	
	lines(x, vals[,1], type = 'b', col='blue')
	if( ncol(vals) > 1 ){
		lines(x, vals[,2], type = 'b', col='green', pch=24)
		legend('bottomright', legend=colnames(vals), pch=c(1, 24))
	}
	if( !is.null(vals.ref) ){
		lines(xref, vals.ref[,1], type='b', col='red')
		if( ncol(vals.ref) > 1 )
			lines(xref, vals.ref[,2], type = 'b', col='pink', pch=24)
	}
	axis(1, at=unique(c(x,xref)), ...)
	axis(2, ...)

}

