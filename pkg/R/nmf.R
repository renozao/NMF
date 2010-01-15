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
if ( !isGeneric("nmf") ) setGeneric('nmf', function(x, rank, method=nmf.getOption('default.algorithm'), ...) standardGeneric('nmf') )
#' Performs a NMF using the default algorithm: brunet.
setMethod('nmf', signature(x='ANY', rank='ANY', method='missing'), 
	function(x, rank, method, ...)
	{
		# apply default algorithm (see default in the generic definition)
		nmf(x, rank, method, ...)
	}
)

#' Performs NMF on an object using a given list of algorithms.
setMethod('nmf', signature(x='ANY', rank='ANY', method='list'), 
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

#' Performs NMF on a data.frame: the target matrix is converted data.frame \code{as.matrix(x)}.
setMethod('nmf', signature(x='data.frame', rank='ANY', method='ANY'), 
function(x, rank, method, ...)
{
	# apply NMF to the the data.frame converted into a matrix
	nmf(as.matrix(x), rank, method, ...)
})

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
		model.parameters <- list()
		if( is.list(model) ){
			model.name <- model[[1]]
			if( !is.character(model.name) || !extends(model.name, 'NMF') )
				stop("First element of argument 'model' must be the name of class that extends class 'NMF'")
			if( length(model) > 1 ) model.parameters <- model[2:length(model)]
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
		validObject(strategy)
		
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
function(x, rank, method, seed=nmf.getOption('default.seed'), nrun=1, model=list(), .options=list(), ...)
{
	# if options are given as a character string, translate it into a list of booleans
	if( is.character(.options) ){
		.options <- .translate.string(.options, 
				c(t='track', v='verbose', d='debug', o='optimized', k='keep.all'))
	}
	
	# setup verbosity options
	debug <- if( !is.null(.options$debug) ) .options$debug else nmf.getOption('debug')
	verbose <- if( !is.null(.options$verbose) ) .options$verbose else nmf.getOption('verbose') || debug
	optimized <- if( !is.null(.options$optimized) ) .options$optimized else nmf.getOption('optimize.mc')
	keep.all <- if( !is.null(.options$keep.all) ) .options$keep.all else FALSE
	
	# make sure rank is an integer
	rank <- as.integer(rank)
	if( length(rank) != 1 ) stop("Invalid rank: must be a single numeric value")
	
	# if the number of run is more than 1, then call itself recursively
	if( nrun > 1 )
	{
		# check seed method: fixed values are not sensible -> warning
		if( inherits(seed, 'NMF') && !is.empty.nmf(seed) )
			warning("nmf: it looks like you are running multiple NMF runs with a fixed seed")
		# if the seed is numerical use it to set the random number generator
		if( is.numeric(seed) ){
			if( verbose ) cat('Set random seed: ', seed, "\n")
			set.seed(seed)
			seed <- 'random'
		}
		
		# define the single run function to use in the (mc)lapply call
		best <- list(res=NULL, residuals=Inf, consensus=NULL)
		if( verbose ) cat('Runs: ')
		single.run <- function(n, ...){
			if( verbose ) cat('', n)
			if( debug ) cat("\n")
			res <- nmf(x, rank, method, nrun=1, seed=seed, .options='-v', ...)
			
			# check if the run found a better fit
			err <- residuals(res)
			if( err < best$residuals ){
				if( n>1 && verbose ){
					if( debug ) cat(": better fit found [err=", err, "]")
					else cat('*')
				}
				best$res <<- res
				best$residuals <<- err				 
			}
			conn <- connectivity(res)
			best$consensus <<- if( is.null(best$consensus) ) conn else (best$consensus + conn)
			if( debug ) cat("\n")
			
			if( !keep.all ) res <- NULL
			
			return(invisible(res))
		}
		
		# test the OS: multicore package does not work on Windows
		# FORCE optimized = FALSE
		#TODO: implement multicore runs for nix machines (see package doMC)
		#if( .Platform$OS.type != 'unix' ){
		if( optimized ){
			#warning('NMF::nmf - turned off multicore optimized mode [only available on *nix machines]')
			warning('NMF::nmf - turned off multicore optimized mode [not yet implemented]')
			optimized = FALSE
		}
				
		# if optimized is TRUE: use mclapply [from multicore package], otherwise use lapply
		if( optimized ){
			if( verbose ) message("Start Multicore mode...")
			library(multicore)
			lapply.fun <- mclapply
		}
		else lapply.fun <- lapply		
		
		# run NMF nrun times
		t <- system.time({res <- lapply.fun(seq(nrun), single.run, ...)})		
		if( verbose ) cat(" ...done\n")
		
		# if one just want the best result only return the best 
		#Note: in this case res is full of NULL any way (see function single.run above)
		if( !keep.all ) 
			return( join(list(best$res), consensus=best$consensus/nrun, runtime=t, nrun=nrun) )
		
		return( join(res, runtime=t) )
	}
	
	if( verbose ) message("NMF algorithm: '", name(method), "'")
	
	# CHECK PARAMETERS:	
	# test for negative values in x only if the method is not mixed
	if( !is.mixed(method) && min(x) < 0 ) stop('Input matrix ', substitute(x),' contains some negative entries.');
	# test if one row contains only zero entries
	if( min(rowSums(x)) == 0) stop('Input matrix ', substitute(x),' contains at least one null row.');	

	# a priori the parameters for the run are all the one in '...'
	parameters.method <- list(...)
	
	if( inherits(seed, 'NMF') ){
		# Wrap up the seed into a NMFfit object
		seed <- new('NMFfit', fit=seed, seed='none')
	}
	else if( !inherits(seed, 'NMFfit') ){
		
		# retrieve the NMF model to use (from the provided method)
		init <- model(method)
		stopifnot( extends(init, 'NMF') )
		
	## MODEL INSTANTIATION :
	# Extra initialization parameters are searched into '...' and argument 'model'
		# extract the parameters from '...' that correspond to slots in the given class
		parameters <- .extract.slots.parameters(init, ...)	
		
		# a priori restrict parameters.method to the ones that won't be used to instantiate the model
		parameters.method <- parameters$extra
		
		# build the list of the model's parameters
		#- a priori all potential arguments is used
		parameters.model <- parameters$slots
		#- override with arguments passed in argument model
		if( !is.list(model) ) 
			stop("invalid argument 'model' [must be a list of values to initialize the NMF model]")
		parameters.model <- .merge.override(parameters.model, model)
		#- force the value of the required arguments for the model
		init <- .merge.override(parameters.model, list(rank=rank, target=0, model=init))		
		#- a posteriori restore overriden arguments: they'll be used by the algorithm
		overriden <- is.element(names(parameters$slots), c(names(model), 'rank', 'target', 'model'))
		if( any( overriden ) )
			parameters.method <- c(parameters.method, parameters$slots[overriden])
		
	# at this point 'init' should be the list of the initialization parameters
	if( !is.list(init) ) stop("Invalid object: 'init' must be a list")
	if( !is.element('model', names(init)) ) stop("Invalid object: 'init' must contain an element named 'model'")	
	
	## SEEDING:	
	# the seed must either be an instance of class 'NMF', the name of a seeding method as a character string
	# or a list of parameters to pass to the 'seed' function.
			parameters.seed <- list()
			seed.method <- NULL
			if( (is.character(seed) && length(seed) == 1) 
				|| is.numeric(seed) ) seed.method <- seed
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
							, "a numerical value used internally to set the seed of the random generator [via 'set.seed']."
							, sep="\n\t- "))
						 			
			# call the 'seed' function passing the necessary parameters
			if( verbose )
				message("NMF seeding method: '", 
						if( is.character(seed.method) ) seed.method
						else if( !is.null(attr(seed.method, 'name')) ) attr(seed.method, 'name') 
						else if( is.function(seed.method) ) 'function'
						else NA
						, "'")
						
			#seed <- do.call(getGeneric('seed', package='NMF')
			seed <- do.call(getGeneric('seed')
					, list(x=x, model=init, method=seed.method, parameters=parameters.seed))
			
			# check the validity of the seed
			if( !inherits(seed, 'NMFfit') ) 
				stop("The seeding method function should return class 'NMF' ["
					, if( is.character(seed.method) ) paste('method "', seed.method, "' ", sep='') else NULL 
					, "returned class: '", class(seed), "']")
	}
	# -> at this point the 'seed' object is an instance of class 'NMF'
	nmf.debug('nmf', "Seed is of class: '", class(seed), "'")
	
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
	parameters.run <- c(list(method, x, seed), parameters.method)
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
if ( is.null(getGeneric('seed')) ) setGeneric('seed', function(x, model, method=nmf.getOption('default.seed'), ...) standardGeneric('seed') )
setMethod('seed', signature(x='ANY', model='ANY', method='missing'),
		function(x, model, method, ...){
			seed(x, model, method, ...)
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
setMethod('seed', signature(x='ANY', model='list', method='ANY'), 
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

setMethod('seed', signature(x='ANY', model='numeric', method='ANY'), 
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
nmfEstimateRank <- function(x, range, method=nmf.getOption('default.algorithm'), nrun=30, conf.interval=FALSE, ...){
	
	concat.to <- function(base, r, value){
		base <- c(base, value)
		names(base)[length(base)] <- r
		return(base)
	}
	
	
	# initiate the list of consensus matrices
	c.matrices <- list()
	bootstrap.measures <- list()
	measures <- sapply(range, function(r, ...){
			message("Compute NMF rank=", r, " ... ", appendLF=FALSE)
			res <- nmf(x, r, method, nrun=nrun
				, .options=list(keep.all=conf.interval, verbose=TRUE), ...)
			
			# sotre the consensus matrix
			c.matrices[[as.character(r)]] <<- consensus(res)
			
			# if confidence intervals must be computed then do it
			if( conf.interval ){
				# resample the tries
				samp <- sapply(seq(5*nrun), function(i){ sample(nrun, nrun, replace=TRUE) })
				
				bootstrap.measures[[as.character(r)]] <<- apply(samp, 2, function(s){
					res.sample <- join(res[s])
					summary(res.sample, target=x)
				})
			}
			
			# compute quality measures
			message('+ measures... ', appendLF=FALSE)
			measures <- summary(res, target=x)

			message('done')
			return(measures)
		}
	, ...)
	
	# set column names for measures
	colnames(measures) <- as.character(range)
	
	# wrap-up result into a 'NMF.rank' S3 object
	res <- list(measures=measures, consensus=c.matrices)
	if( conf.interval ) res$bootstrap.measure <- bootstrap.measures
	class(res) <- 'NMF.rank'
	return(res)
	
}

plot.NMF.rank <- function(x, what=c('all', 'cophenetic', 'rss', 'residuals', 'dispersion'), ... ){

	measures <- x$measures
	what <- match.arg(what)
	
	if( what == 'all' ){
		opar <- par(mfrow=c(2,2))
		on.exit( par(opar), add=TRUE)
		sapply(c('cophenetic', 'rss', 'residuals', 'dispersion'),
				function(w, ...){ plot(x, w, ...) }
		)
		return(invisible())
	}
	
	what <- match.arg(what, rownames(measures))
	vals <- measures[what,]
	x <- colnames(measures)
	
	plot(x=x, y=vals, axes=FALSE, frame=TRUE
		, type='b', lwd=2, col='blue'
		, ylab=paste('Quality measure:', what), xlab='Factorization rank'
		, main=paste("NMF rank estimation\nmeasure:", what), ...)
	axis(1, at=x, ...)
	axis(2, ...)

}

