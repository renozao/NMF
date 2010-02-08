#' @include NMFSet-class.R
#' @include NMFStrategy-class.R

#' Generic interface to Nonnegative Matrix Factorization algorithms.
#'
#' Function \code{nmf} is the main entry point to perform NMF algorithms defined within framework
#' set up by package \code{NMF}.
#' It provides an interface to combine the different algorithms with the different seeding methods. It returns
#' the result as an object of class \code{NMF} that can be directly passed to visualization or benchmarking methods.
#' 
#' The default behaviour of \code{nmf} when \code{method} is missing is to use the algorithm
#' from Brunet et al., which is implemented as a predefined NMFStrategy.
#'
#' @seealso NMFStrategy-class
#' if ( !isGeneric("nmf") ) setGeneric('nmf', function(x, rank, method=nmf.getOption('default.algorithm'), ...) standardGeneric('nmf') )
if ( !isGeneric("nmf") ) setGeneric('nmf', function(x, rank, method, ...) standardGeneric('nmf') )

#' Performs NMF on a data.frame: the target matrix is converted data.frame \code{as.matrix(x)}.
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

#' Performs NMF using an already defined NMF model as starting point (no rank provided)
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
					
					method.potential <- nmfAlgorithm(model=model(seed))
					if( is.null(method.potential) )
						stop("NMF::nmf - No algorithm is defined for model '", model(seed), "'")
					
					if( length(method.potential) == 1 ) # only one to choose
						method <- method.potential
					else if( !is.element(method, method.potential) ) # several options, none is default
						stop("NMF::nmf - Could not infer the algorithm to use with model '", model(seed), "'. Argument 'method' should be one of: "
							, paste(paste("'", method.potential, "'", sep=''), collapse=', ')
							, call.=FALSE)
				}
					
			}
			else if( !existsMethod('nmf', signature=c('matrix', 'numeric', class(method)) ) )
				stop("NMF::nmf - no 'nmf' method for signature 'matrix', 'numeric', '", class(method), "'."
					, call.=FALSE)
			
			# use default seeding method if seed is missing
			if( missing(seed) )
				seed <- nmf.getOption('default.seed')
			
			nmf(x, rank, method, seed=seed, ...)
		}
)


#' Performs NMF on an object using a given list of algorithms.
setMethod('nmf', signature(x='matrix', rank='numeric', method='list'), 
	function(x, rank, method, ...)
	{
		# apply each NMF algorithm
		k <- 1
		t <- system.time({
			res <- lapply(method, 
				function(meth, ...){
					message("Compute NMF method ", k, " ... ", appendLF=FALSE)
					k <<- k+1
					res <- nmf(x, rank, meth, ...)
					message("done")
					return(res)
				}
				, ...)
		})		
		names(res) <- sapply(res, algorithm)		
		
		# wrap the result in a NMFSet object
		res <- join(res, runtime=t)
		
		# return result
		return(res)
	}
)

#' Performs NMF on a matrix using a predefined named strategy. 
#'
#' The available strategies are:
#' - brunet : Based from Brunet et al.
setMethod('nmf', signature(x='matrix', rank='numeric', method='character'),
function(x, rank, method, ...)
{	
	# create the NMFStrategy from its name
	strategy <- nmfAlgorithm(method)		
	# apply nmf using the retrieved strategy		
	nmf(x, rank, method=strategy, ...)
}
)

#' Performs NMF on a matrix using a given function.
#'
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

.translate.string <- function(string, dict){
	
	res <- list()
	dict <- as.list(dict)
	if( nchar(string) == 0 ) return(res)
	opt.val <- TRUE
	lapply(strsplit(string, '')[[1]], 
		function(c){
			if( c=='-' ) opt.val <<- FALSE
			else if( c=='+' ) opt.val <<- TRUE
			else if( !is.null(dict[[c]]) ) res[[dict[[c]]]] <<- opt.val
		}
	)
	
	# return result
	return(res)
}

#' Performs NMF on a matrix using a given NMF method.
#'
#' This method is the entry point for NMF. It is eventually called by any definition of the \code{nmf} function.
setMethod('nmf', signature(x='matrix', rank='numeric', method='NMFStrategy'),
#function(x, rank, method, seed='random', nrun=1, keep.all=FALSE, optimized=TRUE, init='NMF', track, verbose, ...)
function(x, rank, method
		, seed=nmf.getOption('default.seed'), nrun=1, model=NULL, .options=list()
		, .pbackend=nmf.getOption('parallel.backend')
		, ...)
{
	# if options are given as a character string, translate it into a list of booleans
	if( is.character(.options) ){
		.options <- .translate.string(.options, 
				c(t='track', v='verbose', d='debug', p='parallel', P='parallel.required', k='keep.all'))
	}
	
	# setup verbosity options
	debug <- if( !is.null(.options$debug) ) .options$debug else nmf.getOption('debug')
	verbose <- if( !is.null(.options$verbose) ) .options$verbose else nmf.getOption('verbose') || debug
	# parallel runs options
	# backend
	if( is.null(.pbackend) ) .pbackend <- 'seq'
	else if( is.na(.pbackend) ) .pbackend <- 'registered'
	# parallel required: implies option simple parallel
	opt.parallel.required <- if( !is.null(.options$parallel.required) && .options$parallel.required)
								TRUE else FALSE
	if( opt.parallel.required ) .options$parallel <- TRUE
	
	opt.parallel <- if( !is.null(.options$parallel) ) .options$parallel # prioritary on anything else
					else .pbackend != '' # run in parallel only if the backend is defined
			
	keep.all <- if( !is.null(.options$keep.all) ) .options$keep.all else FALSE
	
	# Set debug/verbosity option just for the time of the run
	old.opt <- nmf.options(debug=debug, verbose=verbose); 
	on.exit({nmf.options(old.opt)}, add=TRUE)
	
	# make sure rank is an integer
	rank <- as.integer(rank)
	if( length(rank) != 1 ) stop("NMF::nmf - invalid argument 'rank': must be a single numeric value")
	if( rank <= 1 ) stop("NMF::nmf - invalid argument 'rank': must be greater than 1")
	
	##START_MULTI_RUN
	# if the number of run is more than 1, then call itself recursively
	if( nrun > 1 )
	{
		if( verbose ) message("# NMF Start - Multiple runs: ", nrun)
		
		# check seed method: fixed values are not sensible -> warning
		if( inherits(seed, 'NMF') && !is.empty.nmf(seed) )
			warning("nmf: it looks like you are running multiple NMF runs with a fixed seed")
		# if the seed is numerical use it to set the random number generator
		if( is.numeric(seed) ){
			if( verbose ) message("Set random seed: ", seed)
			set.seed(seed)
			seed <- 'random'
		}
		
		if( opt.parallel ){
			if( verbose ) message("# ", if( opt.parallel.required ) 'Setup required' else 'Try to setup'
						," parallel computation environment")
			
			# check for 'foreach' package: required for parallel computation
			if( !require(foreach) ){
				if( opt.parallel.required )
					stop("NMF::nmf - the 'foreach' package is required to run NMF parallel computation"
							, call.=FALSE)
				else if( verbose ) message("# NOTE: NMF parallel computation disabled ['foreach' package is missing].")
			
				# disable parallel computation
				opt.parallel <- FALSE
			}
		}
		
		# if opt.parallel is TRUE: check and setup everything is there to run in parallel mode
		if( opt.parallel ){
			
			## 0. SETUP PARALLEL MODE
			worker.type <- ''
			single.machine <- TRUE
			ncores <- Inf
			if( is.numeric(.pbackend) && length(.pbackend) > 0 ){
				.pbackend <- .pbackend[1]
				if( .pbackend <= 0 )
					stop("NMF::nmf - invalid number of core(s) specified in argument '.pbackend' [",.pbackend,"]")
				ncores <- .pbackend
				.pbackend <- 'mc'
			}
			switch( .pbackend,
				mc = {
					# test the OS: multicore package does not work on Windows
					if( .Platform$OS.type == 'windows' ){
						# error only if the parallel computation was explicitly asked by the user
						if( opt.parallel.required )
							stop('NMF::nmf - multicore computation impossible [not available under MS Windows]')
						else if( verbose ) 
							message("# NOTE: NMF parallel computation disabled [not availbale under MS Windows].")
						
						# disable parallel computation
						opt.parallel = FALSE
					}
					else if( require(doMC) ){
						
						ncores.machine <- multicore:::detectCores()
						if( ncores.machine == 1 ){
							if( opt.parallel.required )
								stop("NMF::nmf - multicore computation aborted : single core detected"
									, call.=FALSE)
							else if( verbose )
								message("# NOTE: NMF parallel computation disabled [single core detected]")
							opt.parallel = FALSE
						}else{
							
							ncores <- min(ncores.machine, ncores)
							registerDoMC(ncores)
							worker.type <- 'core(s)'
						}
					}
					else if( opt.parallel.required )
						stop("NMF::nmf - missing required package for multicore computation: 'doMC'"
								, call.=FALSE)
					else{
						if( verbose )
							message("# NOTE: NMF multicore computation disabled ['multicore' package not installed]")
						
						# disable parallel computation
						opt.parallel <- FALSE
					}
					
				}
				, seq = {
					registerDoSEQ()
					worker.type <- 'core'
					if( verbose )
						message("# NOTE: this is a SEQUENTIAL computation")
				}
#				, registered = {
#					worker.type <- 'node(s)'
#					if( !getDoParRegistered() )
#						stop("NMF::nmf - no registered backend to run NMF parallel computation")
#				}
				, stop("NMF::nmf - invalid backend ['", .pbackend, "'] for NMF parallel computation. Argument '.pbackend' must be one of 'mc' (or number of cores), 'seq' (or NULL). See ?nmf"
						, call.=FALSE)
			)
		}
		
		
		####PARALLEL_NMF
		if( opt.parallel){
			
			if( !keep.all && !require(bigmemory) )
				stop("NMF::nmf - the 'bigmemory' package is required to run in parallel mode with option 'keep.all'=FALSE"
					, call.=FALSE)
			
			if( verbose ) message("# Using foreach backend: ",getDoParName()
								," [version ", getDoParVersion(),"] / "
								, getDoParWorkers(), " ", worker.type)
	
			run.all <- function(...){
								
				## 1. SETUP				
				if( !keep.all ){ #Specific thing only if one wants only the best result
					# - Define a shared memory objects
					best.shared <- shared.big.matrix(1, 1, type='double', init=NA)			
					best.desc <- bigmemory::describe(best.shared)				
					# the consensus matrix is computed only if not all the results are kept				
					consensus.shared <- shared.big.matrix(ncol(x), ncol(x), type='double', init=0)
					consensus.desc <- bigmemory::describe(consensus.shared)
					
					# - Define a mutex to control the access to the shared memory objects
					mut <- rw.mutex()
					mut.desc <- bigmemory::describe(mut)			
					# - Define a temporary file to store the best fit				
					best.filename <- paste(tempfile('nmf.run.'), 'RData', sep='.')
					##
				}
	
				## 2. RUN
				if( verbose ) cat('Runs:')
				res.runs <- foreach(n=1:nrun, .verbose=debug) %dopar% {
					
					if( verbose ) cat('', n)
					if( debug ) cat("\n")
					res <- nmf(x, rank, method, nrun=1, seed=seed, .options='-v', ...)
					
					# if only the best fit must be kept then update the shared objects
					if( !keep.all ){
						# load shared objects
						best.shared <- attach.big.matrix(best.desc)
						consensus.shared <- attach.big.matrix(consensus.desc)
						
						##LOCK_MUTEX					
						# retrieve and lock the mutex
						mut <- attach.rw.mutex(mut.desc)	
						rwlock(mut)
						
						# check if the run found a better fit
						best <- best.shared[]
						err <- residuals(res)					
						if( is.na(best) || err < best ){
							
							# update residuals
							best.shared[] <- err
							
							# update best fit
							save(res, file=best.filename)
							
						}
						
						# update the consensus matrix
						consensus.shared[] <- consensus.shared[] + connectivity(res)[]
						
						# unlock the mutex
						bigmemory::unlock(mut)
						##END_LOCK_MUTEX
									
						# reset the result to NULL
						res <- NULL
					}
					
					# return the result
					res
				}
				##
				
				
				## 3. WRAP UP
				res <- list(fit=res.runs)
				
				if( !keep.all ){
					# check existence of the result file
					if( !file_test('-f', best.filename) )
						stop("NMF::nmf - error in parallel mode: the result file does not exist")
					load(best.filename)					
					#remove the result file
					unlink(best.filename)
					
					res$fit <- res
					res$consensus <- consensus.shared[]
				}
				##
				
				# return result
				res
			}			
		}####END_PARALLEL_NMF
		else{####SEQUENTIAL_NMF
			
			if( verbose ) message("# Using standard sequential computation")
			run.all <- function(...){
				
				## 1. SETUP
				
				# define static variables for the case one only wants the best result
				if( !keep.all ){
					# statis list with best result: fit, residual, consensus
					best.static <- list(fit=NULL, residuals=NA, consensus=matrix(0, ncol(x), ncol(x)))					
				}
								
				# define a function that performs a single NMF and update the static variables
				single.run <- function(n, ...){
					if( verbose ) cat('', n)
					if( debug ) cat("\n")
					res <- nmf(x, rank, method, nrun=1, seed=seed, .options='-v', ...)
					
					if( !keep.all ){						
						# check if the run found a better fit
						err <- residuals(res)
						best <- best.static$residuals
						if( is.na(best) || err < best ){
							if( n>1 && verbose ){
								if( debug ) cat(": better fit found [err=", err, "]")
								else cat('*')
							}
							
							# update best fit (only if necessary)
							best.static$fit <<- res
							
							best.static$residuals <<- err				 
						}
						
						# update the static consensus matrix (only if necessary)
						best.static$consensus <<- best.static$consensus + connectivity(res)
						
						# reset the result to NULL
						res <- NULL
					}
					
					if( debug ) cat("\n")
					
					# return the result
					res
				}
				##
				
				## 2. RUN:
				# run 'single.run' nrun times
				if( verbose ) cat('Runs:')
				res.runs <- lapply(seq(nrun), single.run, ...)
				##
				
				## 3. WRAP UP
				res <- list(fit=res.runs)
				
				if( !keep.all ){
					res$fit <- best.static$fit
					res$consensus <- best.static$consensus
				}
				##
				
				# return the result
				res
			}
			
		}####END_SEQUENTIAL_NMF
		
		# perform all the NMF runs
		t <- system.time({res <- run.all(...)})
		if( verbose ) cat(" ... DONE\n")
		
		# ASSERT the presence of the result
		stopifnot( !is.null(res$fit) )
		
		# if one just want the best result only return the best 
		#Note: in this case res is full of NULL any way (see function single.run above)
		if( !keep.all ){
			# ASSERT the presence of the consensus matrix
			stopifnot( !is.null(res$consensus) )
			return( join(list(res$fit), consensus=res$consensus/nrun, runtime=t, nrun=nrun) )
		}
		
		return( join(res$fit, runtime=t) )
		
	}##END_MULTI_RUN
	
	if( verbose ) message("NMF algorithm: '", name(method), "'")
	
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
		if( inherits(seed, 'NMFfit') )
			seed <- fit(seed)
		
		# Wrap up the seed into a NMFfit object
		seed <- new('NMFfit', fit=seed, seed='none')
	}
	else if( !inherits(seed, 'NMFfit') ){
		
		# retrieve the NMF model to use (from the provided method)
		init <- model(method)
		stopifnot( extends(init, 'NMF') )
		
		## MODEL INSTANTIATION :
	
		# some of the instantiation parameters are set internally
		# TODO: change target into x (=> impact on newNMF)
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
				|| is.null(seed) ) seed.method <- seed
			else if( is.function(seed) ) seed.method <- seed
			else if( is.list(seed) ){ # seed is a list...
				
				if( !is.null(seed$method) ){ # 'seed' must contain an element giving the method...
					seed.method <- seed$method
					parameters.seed <- seed[-which(names(seed)=='method')]
				}
				else if ( names(seed)[1] == '' ){ # ... or the first element must be a method
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
							, "the name of a seeding method"
							, "a valid function definition"
							, "a list containing the seeding method (i.e. a function or a character string) as its first element\n\tor as an element named 'method' [and optionnally extra arguments it will be called with]"
							, "a numerical value used internally to set the seed of the random generator [via 'set.seed']"
							, "NULL to directly pass the model instanciated from arguments 'model' or '...'."
							, sep="\n\t- "))
						 			
			# call the 'seed' function passing the necessary parameters
			if( verbose )
				message("NMF seeding method: ", 
						if( is.character(seed.method) ) seed.method
						else if( is.null(seed.method) ) 'NULL'
						else if( !is.null(attr(seed.method, 'name')) ) attr(seed.method, 'name') 
						else if( is.function(seed.method) ) '<function>'
						else NA)
						
			#seed <- do.call(getGeneric('seed', package='NMF')
			seed <- do.call(getGeneric('seed')
					, list(x=x, model=init, method=seed.method, parameters=parameters.seed))
			
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
	if( all( !inherits(fit(seed), model(method)) ) )
		stop("NMF::nmf - Invalid NMF model '", model(seed),"': algorithm '", name(method), "' expects model(s) "
			, paste(paste("'", model(method),"'", sep=''), collapse=', ')
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
	if( !missing(verbose) ) run.options(seed, 'verbose') <- verbose

	## RUN NMF METHOD:
	# call the strategy's run method [and time it] using the element of 'parameters.method' as parameters	
	parameters.run <- c(list(method=method, x=x, seed=seed), parameters.method)
	t <- system.time({res <- do.call('run', parameters.run)})
	
	## CHECK RESULT
	# check the result is of the right type
	if( !inherits(res, 'NMFfit') ) stop("NMF method should return an instance of class 'NMF' [returned class:", class(res), "]")			

	## ENSURE SOME SLOTS ARE STILL CORRECTLY SET
	# slot 'method'
	algorithm(res) <- name(method)	
	# slot 'distance'
	res@distance <- objective(method)	
	# slot 'seed'
	if( seed.method != '' ) seeding(res) <- seed.method
	# slot 'parameters'
	res@parameters <- parameters.method	
	
	## CLEAN-UP + EXTRAS:
	# add extra information to the object	 
	if( length(residuals(res)) == 0 ) residuals(res) <- objective(method, x, res)
	if( is.na(residuals(res)) ) stop("NMF residuals: final objective value is NA")
	res@runtime <- t
	
	# return the result
	res
})

#' Common interface for seeding methods for Nonnegative Matrix Factorization (NMF) algorithms.
#' 
#' This function calls the different seeding methods that define a starting point for NMF methods.
#' These methods at least set the slots \code{W} and \code{H} of \code{object} to valid nonnegative matrices.
#' They will be used as a starting point by any NMF algorithm that accept initialization.
#'
#' @param x The target matrix one wants to approximate with NMF
#' @param rank The rank of the factorization to seed
#' @param method Name of the seeding method to call
#' @param object Either an instance (resp. the name) of a NMF (sub)class, that will be seeded (resp. instanciated and seeded).
#' @param ... Parameters to be passed to the seeding method
#'  
#' @return \code{object} initialized using method \code{name}.
#'  
if ( is.null(getGeneric('seed')) ) setGeneric('seed', function(x, model, method, ...) standardGeneric('seed') )
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
			
			# check value
			if( length(method) == 0 )
				stop('NMF::seed - numeric seed is empty [single value expected]')
			if( length(method) > 1 ) 
				warning('NMF::seed - numeric seed has length > 1, only the first element will be used')
			# only use first element
			method <- method[1]
			
			# set the seed using the numerical value by argument 'method'
			set.seed(method)
			# call seeding method 'random'
			res <- seed(x, model, 'random', ...)
			# set seeding method to the numeric seed
			seeding(res) <- as.character(method)
			
			# return result
			return(res)
		}
)
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
		model <- do.call('newNMF', model)
		nmf.debug('seed', "using NMF model '", class(model), "'")
		
		# check that model is from the right type, i.e. inherits from class NMF
		if( !inherits(model, 'NMF') ) stop("Invalid object returned by model: object must inherit from class 'NMF'")
		
		seed(x, model, method, ...)
	}
)

setMethod('seed', signature(x='ANY', model='numeric', method='NMFSeed'), 
	function(x, model, method, ...){	

		seed(x, list(model='NMFstd', rank=model), method, ...)
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
	function(x, model, method, parameters=list(), ...){	
				
		# debug message
		nmf.debug('seed', "use seeding method: '", name(method), "'")
		
		# use argument '...' for the method's parameters if they are not supplied by argument 'parameters'
		if( !is.list(parameters) || length(parameters)==0 ) parameters <- list(...)						
		
		# call the seeding function passing the extra parameters
		res <- do.call(method(method), c(list(model, x), parameters))
		
		res <- new('NMFfit', fit=res)
		# if not already set: store the seeding method's name in the resulting object
		if( seeding(res) == '' ) seeding(res) <- name(method)
		
		# return the seeded object
		res
	}
)

#' Extract from a list the elements that can be used to initialize the slot of a class.
#' 
#' This function only extract named elements.
#' 
#' @param class.name Name of the class from whose slots will be search into '...'
#' @param ... The parameters in which the slot names will be search for
#' 
#' @returnType list
#' @return a list with two elements:
#' - \code{slots}: is a list that contains the named parameters that can be used to instantiate an object of class \code{class.name} 
#' - \code{extra}: is a list of the remaining parameters from \code{parameters} (i.e. the ones that do not correspond to a slot).
#'  
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

#' Merges two lists, but overriding with the values of the second list in the case
#' of duplicates.
.merge.override <- function(l1, l2, warning=FALSE){
	sapply(names(l2), function(name){
				if( warning && !is.null(l1[[name]]) )
					warning("overriding element '", name, "'")
				l1[[name]] <<- l2[[name]]
			})
	
	# return updated list
	return(l1)
}

#' Estimate the factorization rank
nmfEstimateRank <- function(x, range, method=nmf.getOption('default.algorithm')
					, nrun=30, verbose=FALSE, ...){
	
	concat.to <- function(base, r, value){
		base <- c(base, value)
		names(base)[length(base)] <- r
		return(base)
	}
	
	
	# initiate the list of consensus matrices
	c.matrices <- list()
	bootstrap.measures <- list()
	measures <- sapply(range, function(r, ...){
			if( verbose ) message("Compute NMF rank=", r, " ... ", appendLF=FALSE)
			res <- nmf(x, r, method, nrun=nrun
				, .options=list(keep.all=FALSE, verbose=verbose), ...)
			
			# sotre the consensus matrix
			c.matrices[[as.character(r)]] <<- consensus(res)
			
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
			if( verbose ) message('+ measures... ', appendLF=FALSE)
			measures <- summary(res, target=x)

			if( verbose ) message('done')
			return(measures)
		}
	, ...)
	
	# set column names for measures
	colnames(measures) <- as.character(range)
	measures <- as.data.frame(t(measures))
	
	# wrap-up result into a 'NMF.rank' S3 object
	res <- list(measures=measures, consensus=c.matrices)
	#if( conf.interval ) res$bootstrap.measure <- bootstrap.measures
	class(res) <- 'NMF.rank'
	return(res)
	
}

plot.NMF.rank <- function(x, what=c('all', 'cophenetic', 'rss', 'residuals'
									, 'dispersion', 'evar', 'sparseness'
									, 'sparseness.basis', 'sparseness.coef')
						, ref=NULL, ... ){

	measures <- x$measures
	what <- match.arg(what)
	
	if( !exists('.NMF.rank.plot.notitle', parent.frame()) ){
		.NMF.rank.plot.notitle <- FALSE
	}
	else .NMF.rank.plot.notitle <- TRUE
	
	if( what == 'all' ){
		opar <- par(mfrow=c(2,3), oma=c(0,0,3,0))
		on.exit( par(opar), add=TRUE)
		sapply(c('cophenetic', 'rss', 'residuals'
				, 'dispersion', 'evar', 'sparseness'),
				function(w, ...){
					.NMF.rank.plot.notitle <- TRUE
					plot(x, w, ref, ...) 
				}
		)
		title("NMF rank estimation", outer=TRUE)
		return(invisible())
	}
	
	iwhat <- grep(paste('^',what,sep=''), colnames(measures))
	vals <- measures[,iwhat, drop=FALSE]
	x <- as.numeric(rownames(measures))
	xlim <- range(x)
	
	vals.ref <- NULL
	if( !missing(ref) && is(ref, 'NMF.rank') ){
		xref <- as.numeric(rownames(ref$measures))
		xlim <- range(xlim, xref)
		vals.ref <- ref$measures[,iwhat, drop=FALSE]
	}
	
	# compute the ylim from main and ref values
	ylim <- range(vals)
	if( !is.null(vals.ref) )
		ylim <- range(ylim, vals.ref)
	
	# detect if the values should be between 0 and 1
	#if( all(ylim >=0 & ylim <= 1) )
	#	ylim <- c(0,1)
	
	#init the plot
	plot(x=NULL, axes=FALSE, frame=TRUE
		, xlim=xlim, ylim = ylim
		, ylab=paste('Quality measure:', what), xlab='Factorization rank'
		, main=paste(if(!.NMF.rank.plot.notitle) "NMF rank estimation\n"
					,"-", what, '-'), ...)
	
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
	axis(1, at=x, ...)
	axis(2, ...)

}

