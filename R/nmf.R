#' @include NMFstd-class.R
#' @include NMFSet-class.R
#' @include registry-seed.R
#' @include registry-algorithms.R
#' @include parallel.R
NULL

#' Running NMF algorithms
#'
#' @description
#' The function \code{nmf} is a S4 generic defines the main interface to run NMF 
#' algorithms within the framework defined in package \code{NMF}.
#' It has many methods that facilitates applying, developing and testing NMF 
#' algorithms.
#' 
#' The package vignette \code{vignette('NMF')} contains an introduction to the 
#' interface, through a sample data analysis.
#' 
#' @details
#'
#' The \code{nmf} function has multiple methods that compose a very flexible 
#' interface allowing to:
#' \itemize{
#' \item combine NMF algorithms with seeding methods and/or stopping/convergence 
#' criterion at runtime;
#' 
#' \item perform multiple NMF runs, which are computed in parallel whenever the host 
#' machine allows it;
#' 
#' \item run multiple algorithms with a common set of parameters, ensuring a 
#' consistent environment (notably the RNG settings).
#' }
#' 
#' The workhorse method is \code{nmf,matrix,numeric,NMFStrategy}, which is eventually 
#' called by all other methods.
#' The other methods provides convenient ways of specifying the NMF algorithm(s),  
#' the factorization rank, or the seed to be used.
#' Some allow to directly run NMF algorithms on different types of objects, such 
#' as \code{data.frame} or \code{\link[Biobase]{ExpressionSet}} objects.
#' 
#' @section Optimized C++ vs. plain R: 
#' Lee and Seung's multiplicative updates are used by several NMF algorithms. To improve 
#' speed and memory usage, a C++ implementation of the specific matrix products is used 
#' whenever possible. It directly computes the updates for each entry in the updated matrix, 
#' instead of using multiple standard matrix multiplication.
#' 
#' The algorithms that benefit from this optimization are: 'brunet', 'lee', 'nsNMF' and 'offset'. % and 'lnmf'
#' However there still exists plain R versions for these methods, which implement the updates 
#' as standard matrix products. These are accessible by adding the prefix '.R#' to their name: 
#' '.R#brunet', '.R#lee', '.R#nsNMF' and '.R#offset'.
#' 
#' @param x target data to fit, i.e. a matrix-like object
#' @param rank specification of the factorization rank.
#' It is usually a single numeric value, but other type of values are possible 
#' (e.g. matrix), for which specific methods are implemented.
#' See for example methods \code{nmf,matrix,matrix,ANY}.
#' 
#' If \code{rank} is a numeric vector with more than one element, e.g. a range of ranks, 
#' then \code{\link{nmf}} performs the estimation procedure described in 
#' \code{\link{nmfEstimateRank}}. 
#' 
#' @param method specification of the NMF algorithm.
#' The most common way of specifying the algorithm is to pass the access key 
#' (i.e. a character string) of an algorithm stored in the package's dedicated registry, 
#' but methods exists that handle other types of values, such as \code{function} or \code{list} 
#' object. See their descriptions in section \emph{Methods}.
#' 
#' If \code{method} is missing the algorithm to use is obtained from the option 
#' \code{nmf.getOption('default.algorithm')}, unless it can be infer from the type of NMF model 
#' to fit, if this later is available from other arguments. 
#' Factory fresh default value is \sQuote{brunet}, which corresponds to the standard NMF 
#' algorithm from \cite{Brunet2004} (see section \emph{Algorithms}).
#' 
#' Cases where the algorithm is inferred from the call are when an NMF model is passed in arguments \code{rank} 
#' or \code{seed} (see description for \code{nmf,matrix,numeric,NULL} in section \emph{Methods}).
#'  
#' @param ... extra arguments to allow extension of the generic.
#' Arguments that are not used in the chain of internal calls to \code{nmf} methods 
#' are passed to the function that effectively implements the algorithm that fits 
#' an NMF model on \code{x}.
#' 
#' @export
#' @inline
#'
#' @examples
#' 
#' # Only basic calls are presented in this manpage.
#' # Many more examples are provided in the demo file nmf.R
#' \dontrun{
#' demo('nmf')
#' }
#' 
#' # random data
#' x <- rmatrix(20,10)
#' 
#' # run default algorithm with rank 2
#' res <- nmf(x, 2)
#' 
#' # specify the algorithm
#' res <- nmf(x, 2, 'lee')
#' 
#' # get verbose message on what is going on
#' res <- nmf(x, 2, .options='v') 
#' \dontrun{ 
#' # more messages
#' res <- nmf(x, 2, .options='v2')
#' # even more
#' res <- nmf(x, 2, .options='v3')
#' # and so on ... 
#' }
#'  
#' @demo Using the main function nmf()
#' 
#' # generate a synthetic dataset with known classes: 50 features, 23 samples (10+5+8)
#' n <- 20; counts <- c(5, 3, 2);
#' p <- sum(counts)
#' x <- syntheticNMF(n, counts)
#' dim(x)
#' 
#' # build the true cluster membership
#' groups <- unlist(mapply(rep, seq(counts), counts))
#' 
setGeneric('nmf', function(x, rank, method, ...) standardGeneric('nmf') )
#' Fits an NMF model on a \code{data.frame}.
#' 
#' The target \code{data.frame} is coerced into a matrix with \code{\link{as.matrix}}.
#' 
#' @demo
#' 
#' # run on a data.frame
#' res <- nmf(data.frame(x), 3)
#' 
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
#' Fits an NMF model using an appropriate algorithm when \code{method} is not supplied.
#' 
#' This method tries to select an appropriate algorithm amongst the NMF algorithms 
#' stored in the internal algorithm registry, which contains the type of NMF models 
#' each algorithm can fit.
#' This is possible when the type of NMF model to fit is available from argument \code{seed}, 
#' i.e. if it is an NMF model itself.
#' Otherwise the algorithm to use is obtained from \code{nmf.getOption('default.algorithm')}.
#' 
#' This method is provided for internal usage, when called from other \code{nmf} methods 
#' with argument \code{method} missing in the top call (e.g. \code{nmf,matrix,numeric,missing}).
#' 
#' @demo
#' 
#' # missing method: use algorithm suitable for seed
#' res <- nmf(x, 2, seed=rnmf(2, x))
#' algorithm(res)
#' res <- nmf(x, 2, seed=rnmf(2, x, model='NMFns'))
#' algorithm(res)
#' 
setMethod('nmf', signature(x='matrix', rank='numeric', method='NULL'), 
		function(x, rank, method, seed=NULL, model=NULL, ...)
		{
			
			# a priori the default method will be used
			method <- nmf.getOption('default.algorithm')
			
			# use default seeding method if seed is missing
			if( is.null(seed) ){
#				seed <- nmf.getOption('default.seed')
			}else{
				# get reference object from which to infer model type
				refobj <- if( is.nmf(seed) ) seed else if( is.nmf(model) ) model
				
				if( !is.null(refobj) ){
					mtype <- modelname(refobj)
					# try to find the algorithm suitable for the seed's NMF model
					method.potential <- selectNMFMethod(model=mtype, exact=TRUE, quiet=TRUE)
					if( is.null(method.potential) )
						stop("NMF::nmf - Found no algorithm defined for model '", mtype, "'")
					
					if( length(method.potential) == 1 ) # only one to choose
						method <- method.potential
					else if( !is.element(method, method.potential) ){# several options, none is default
						method <- method.potential[1]
						warning("NMF::nmf - Selected algorithm '", method, "' to fit model '", mtype, "'."
							, "\n  Alternatives are: "
							, str_out(method.potential[-1], Inf)
							, call.=FALSE, immediate.=TRUE)
					}
				}
			}
			
			nmf(x, rank, method, seed=seed, model=model, ...)
		}
)


#' Fits multiple NMF models on a common matrix using a list of algorithms.
#' 
#' The models are fitted sequentially with \code{nmf} using the same options 
#' and parameters for all algorithms.
#' In particular, irrespective of the way the computation is seeded, this method 
#' ensures that all fits are performed using the same initial RNG settings.
#' 
#' This method returns an object of class \code{\linkS4class{NMFList}}, that is  
#' essentially a list containing each fit.  
#' 
#' @param .parameters list of method-specific parameters.
#' Its elements must have names matching a single method listed in \code{method},
#' and be lists of named values that are passed to the corresponding method. 
#' 
#' @demo 
#' # compare some NMF algorithms (tracking the approximation error)
#' res <- nmf(x, 2, list('brunet', 'lee', 'nsNMF'), .options='t')
#' res
#' summary(res, class=groups)
#' 
#' # plot the track of the residual errors
#' plot(res)
#' 
setMethod('nmf', signature(x='matrix', rank='numeric', method='list'), 
	function(x, rank, method, ..., .parameters = list())
	{
		# apply each NMF algorithm
		k <- 0
		n <- length(method)
        
        # setup/check method specific parameters
        ARGS <- NULL
        .used.parameters <- character()
        if( !is.list(.parameters) )
            stop("NMF::nmf - Invalid value for argument `.parameters`: must be a named list.")
        if( length(.parameters) && (is.null(names(.parameters)) || any(names(.parameters) == '')) )
            stop("NMF::nmf - Invalid value for argument `.parameters`: all elements must be named.") 
        
        t <- system.time({
			res <- lapply(method, 
				function(meth, ...){
					k <<- k+1
					methname <- if( isString(meth) ) meth else name(meth)
					cat("Compute NMF method '", methname, "' [", k, "/", n, "] ... ", sep='')
					# restore RNG on exit (except after last method)
					# => this ensures the methods use the same stochastic environment
					orng <- RNGseed()
					if( k < n ) on.exit( RNGseed(orng), add = TRUE)
					
                    # look for method-specific arguments
                    i.param <- 0L
                    if( length(.parameters) ){
                        i.param <- charmatch(names(.parameters), methname)
                        if( !length(i.param <- seq_along(.parameters)[!is.na(i.param)]) )
                            i.param <- 0L
                        else if( length(i.param) > 1L ){
                            stop("Method name '", methname, "' matches multiple method-specific parameters "
                                    , "[", str_out(names(.parameters)[i.param], Inf), "]")
                        }
                    }
					#o <- capture.output( 
                        if( !i.param ){
                            res <- try( nmf(x, rank, meth, ...) , silent=TRUE)
                        }else{
                            if( is.null(ARGS) ) ARGS <<- list(x, rank, ...)
                            .used.parameters <<- c(.used.parameters, names(.parameters)[i.param])
                            res <- try( do.call(nmf, c(ARGS, method = meth, .parameters[[i.param]])) 
                                        , silent=TRUE)
                        } 
					#)
					if( is(res, 'try-error') )
						cat("ERROR\n")
					else 
						cat("OK\n")
					return(res)
				}
				, ...)
		})
		
		# filter out bad results
		ok <- sapply(res, function(x){
					if( is(x, 'NMF.rank') ) all(sapply(x$fit, isNMFfit))
					else isNMFfit(x)
			})
		if( any(!ok) ){ # throw warning if some methods raised an error
			err <- lapply(which(!ok), function(i){ paste("'", method[[i]],"': ", res[[i]], sep='')})
			warning("NMF::nmf - Incomplete results due to ", sum(!ok), " errors: \n- ", paste(err, collapse="- "), call.=FALSE)
		}
		res <- res[ok]
		# TODO error if ok is empty

        # not-used parameters
        if( length(.used.parameters) != length(.parameters) ){
            warning("NMF::nmf - Did not use methods-specific parameters ", str_out(setdiff(names(.parameters), .used.parameters), Inf))
        }

		# add names to the result list
		names(res) <- sapply(res, function(x){
					if( is(x, 'NMF.rank') ) x <- x$fit[[1]]
					algorithm(x)
				})
				
		# return list as is if surveying multiple ranks 
		if( length(rank) > 1 ) return(res)
		
		# wrap the result in a NMFList object
		# DO NOT WRAP anymore here: NMFfitX objects are used only for results of multiple runs (single method)
		# the user can still join on the result if he wants to
		#res <- join(res, runtime=t)
		res <- new('NMFList', res, runtime=t)
		
		# return result
		return(res)
	}
)
#' Fits an NMF model on \code{x} using an algorithm registered with access key 
#' \code{method}.
#' 
#' Argument \code{method} is partially match against the access keys of all 
#' registered algorithms (case insensitive).
#' Available algorithms are listed in section \emph{Algorithms} below or the 
#' introduction vignette. 
#' A vector of their names may be retrieved via \code{nmfAlgorithm()}.
#' 
#' @section Algorithms:
#' All algorithms are accessible by their respective access key as listed below.
#' The following algorithms are available:
#' \describe{
#' 
#' \item{\sQuote{brunet}}{ Standard NMF, based on the Kullback-Leibler divergence, 
#' from \cite{Brunet2004}.
#' It uses simple multiplicative updates from \cite{Lee2001}, enhanced to avoid 
#' numerical underflow.
#' 
#' Default stopping criterion: invariance of the connectivity matrix
#' (see \code{\link{nmf.stop.connectivity}}).
#' }
#' 
#' \item{\sQuote{lee}}{ Standard NMF based on the Euclidean distance from \cite{Lee2001}.
#' It uses simple multiplicative updates.
#' 
#' Default stopping criterion: invariance of the connectivity matrix
#' (see \code{\link{nmf.stop.connectivity}}).
#' }
#' 
#' \item{ls-nmf}{ Least-Square NMF from \cite{Wang2006}.
#' It uses modified versions of Lee and Seung's multiplicative updates for the 
#' Euclidean distance, which incorporates weights on each entry of the target 
#' matrix, e.g. to reflect measurement uncertainty.
#' 
#' Default stopping criterion: stationarity of the objective function
#' (see \code{\link{nmf.stop.stationary}}).
#' }
#' 
#' \item{\sQuote{nsNMF}}{ Nonsmooth NMF from \cite{Pascual-Montano2006}. 
#' It uses a modified version of Lee and Seung's multiplicative updates for the 
#' Kullback-Leibler divergence \cite{Lee2001}, to fit a extension of the standard 
#' NMF model, that includes an intermediate smoothing matrix, meant meant to produce 
#' sparser factors.
#' 
#' Default stopping criterion: invariance of the connectivity matrix
#' (see \code{\link{nmf.stop.connectivity}}).
#' }
#' 
#' \item{\sQuote{offset}}{ NMF with offset from \cite{Badea2008}.
#' It uses a modified version of Lee and Seung's multiplicative 
#' updates for Euclidean distance \cite{Lee2001}, to fit an NMF model that includes 
#' an intercept, meant to capture a common baseline and shared patterns, in 
#' order to produce cleaner basis components.
#' 
#' Default stopping criterion: invariance of the connectivity matrix
#' (see \code{\link{nmf.stop.connectivity}}).
#' }
#' 
#' \item{\sQuote{pe-nmf}}{ Pattern-Expression NMF from \emph{Zhang2008}. 
#' It uses multiplicative updates to minimize an objective function based on the 
#' Euclidean distance, that is regularized for effective expression of patterns 
#' with basis vectors.
#' 
#' Default stopping criterion: stationarity of the objective function
#' (see \code{\link{nmf.stop.stationary}}).
#' }
#' 
#' \item{\sQuote{snmf/r}, \sQuote{snmf/l}}{ Alternating Least Square (ALS) approach 
#' from \cite{KimH2007}.
#' It applies the nonnegative least-squares algorithm from \cite{VanBenthem2004} 
#' (i.e. fast combinatorial nonnegative least-squares for multiple right-hand), 
#' to estimate the basis and coefficient matrices alternatively 
#' (see \code{\link{fcnnls}}).
#' It minimises an Euclidean-based objective function, that is regularized to 
#' favour sparse basis matrices (for \sQuote{snmf/l}) or sparse coefficient matrices 
#' (for \sQuote{snmf/r}).
#' 
#' Stopping criterion: built-in within the internal workhorse function \code{nmf_snmf}, 
#' based on the KKT optimality conditions.
#' }
#' 
#' }
#' 
#' @section Seeding methods:
#' The purpose of seeding methods is to compute initial values for the factor 
#' matrices in a given NMF model. 
#' This initial guess will be used as a starting point by the chosen NMF algorithm.
#' 
#' The seeding method to use in combination with the algorithm can be passed 
#' to interface \code{nmf} through argument \code{seed}.
#' The seeding seeding methods available in registry are listed by the function 
#' \code{\link{nmfSeed}} (see list therein). 
#' 
#' Detailed examples of how to specify the seeding method and its parameters can 
#' be found in the \emph{Examples} section of this man page and in the package's 
#' vignette.		
#' 
#' @seealso \code{\link{nmfAlgorithm}}
#' 
#' @demo 
#' 
#' # specify algorithm by its name
#' res <- nmf(x, 3, 'nsNMF', seed=123) # nonsmooth NMF 
#' # names are partially matched so this also works  
#' identical(res, nmf(x, 3, 'ns', seed=123))
#' 
#' res <- nmf(x, 3, 'offset') # NMF with offset
#' 
#' 
setMethod('nmf', signature(x='matrix', rank='numeric', method='character'),
function(x, rank, method, ...)
{	
	# if there is more than one methods then treat the vector as a list
	if( length(method) > 1 ){
		return( nmf(x, rank, as.list(method), ...) )
	}
	
	# create the NMFStrategy from its name
	strategy <- nmfAlgorithm(method)		
	# apply nmf using the retrieved strategy		
	nmf(x, rank, method=strategy, ...)
}
)
#' Fits an NMF model on \code{x} using a custom algorithm defined the function 
#' \code{method}. 
#' 
#' The supplied function must have signature \code{(x=matrix, start=NMF, ...)} 
#' and return an object that inherits from class \code{\linkS4class{NMF}}. 
#' It will be called internally by the workhorse \code{nmf} method, with an NMF model 
#' to be used as a starting point passed in its argument \code{start}.
#' 
#' Extra arguments in \code{...} are passed to \code{method} from the top 
#' \code{nmf} call.
#' Extra arguments that have no default value in the definition of the function 
#' \code{method} are required to run the algorithm (e.g. see argument \code{alpha} 
#' of \code{myfun} in the examples).
#' 
#' If the algorithm requires a specific type of NMF model, this can be specified 
#' in argument \code{model} that is handled as in the workhorse \code{nmf} 
#' method (see description for this argument).
#' 
#' @param name name associated with the NMF algorithm implemented by the function
#' \code{method} [only used when \code{method} is a function].
#' @param objective specification of the objective function associated with the 
#' algorithm implemented by the function \code{method} 
#' [only used when \code{method} is a function].
#' 
#' It may be either \code{'euclidean'} or \code{'KL'} for specifying the euclidean 
#' distance (Frobenius norm) or the Kullback-Leibler divergence respectively, 
#' or a function with signature \code{(x="NMF", y="matrix", ...)} that computes
#' the objective value for an NMF model \code{x} on a target matrix \code{y}, 
#' i.e. the residuals between the target matrix and its NMF estimate.
#' Any extra argument may be specified, e.g. \code{function(x, y, alpha, beta=2, ...)}.
#' 
#' @param mixed a logical that indicates if the algorithm implemented by the function
#' \code{method} support mixed-sign target matrices, i.e. that may contain negative 
#' values [only used when \code{method} is a function].
#' 
#' @demo 
#' 
#' # run a custom algorithm defined as a standard function
#' myfun <- function(x, start, alpha){
#' 	# update starting point
#' 	# ...
#' 	basis(start) <- 3 * basis(start)
#' 	# return updated point
#' 	start 
#' }
#' 
#' res <- nmf(x, 2, myfun, alpha=3)
#' algorithm(res)
#' # error: alpha missing
#' try( nmf(x, 2, myfun) )
#' 
#' # possibly the algorithm fits a non-standard NMF model, e.g. NMFns model
#' res <- nmf(x, 2, myfun, alpha=3, model='NMFns')
#' modelname(res)
#' 
setMethod('nmf', signature(x='matrix', rank='numeric', method='function'),
	function(x, rank, method, seed, model='NMFstd', ..., name, objective='euclidean', mixed=FALSE){

		model_was_a_list <- is.list(model)
		if( is.character(model) )
			model <- list(model=model) 
		if( !is.list(model) ){
			stop("nmf - Invalid argument `model`: must be NULL or a named list of initial values for slots in an NMF model.")
		}
				
		
		# arguments passed to the call to NMFStrategyFunction
		strat <- list('NMFStrategyFunction'
					, algorithm = method
					, objective = objective
					, mixed = mixed[1]
					)
		
		## Determine type of NMF model associated with the NMFStrategy
		# All elements of `model` (except the model class) will be passed to 
		# argument `model` of the workhorse `nmf` method, which will use them  
		# to create the NMF model in a call to `nmfModel`
		if( length(model) > 0L ){
			if( !is.null(model$model) ){
				strat$model <- model$model
				model$model <- NULL
			}else if( isNMFclass(model[[1]]) ){
				strat$model <- model[[1]]
				# use the remaining elements to instanciate the NMF model
				model <- model[-1]
			}
			# all elements must be named
			if( !hasNames(model, all=TRUE) ){
				stop("NMF::nmf - Invalid argument `model`: all elements must be named, except the first one which must then be an NMF model class name")
			}
		}
		##
		
		# if name is missing: generate a temporary unique name
		if( missing(name) ) name <- basename(tempfile("nmf_"))
		# check that the name is not a registered name
		if( existsNMFMethod(name) )
			stop("Invalid name for custom NMF algorithm: '",name,"' is already a registered NMF algorithm")
		strat$name <- name

		# create NMFStrategy
		strategy <- do.call('new', strat)
		# full validation of the strategy
		validObject(strategy, complete=TRUE)
		
		if( missing(seed) ) seed <- NULL
		if( !model_was_a_list && length(model) == 0L ) model <- NULL
		# call method 'nmf' with the new object
		nmf(x, rank, strategy, seed=seed, model=model, ...)
	}
)

#' Fits an NMF model using the NMF model \code{rank} to seed the computation, 
#' i.e. as a starting point.
#' 
#' This method is provided for convenience as a shortcut for 
#' \code{nmf(x, nbasis(object), method, seed=object, ...)} 
#' It discards any value passed in argument \code{seed} and uses the NMF model passed 
#' in \code{rank} instead.
#' It throws a warning if argument \code{seed} not missing.
#' 
#' If \code{method} is missing, this method will call the method 
#' \code{nmf,matrix,numeric,NULL}, which will infer an algorithm suitable for fitting an 
#' NMF model of the class of \code{rank}.
#' 
#' @demo
#' 
#' # assume a known NMF model compatible with the matrix `x`
#' y <- rnmf(3, x)
#' # fits an NMF model (with default method) on some data using y as a starting point
#' res <- nmf(x, y)
#' # the fit can be reproduced using the same starting point
#' nmf.equal(nmf(x, y), res)
#'  
setMethod('nmf', signature(x='matrix', rank='NMF', method='ANY'),
	function(x, rank, method, seed, ...){
		
		if( !missing(seed) ){		
			if( isNumber(seed) ){
				set.seed(seed)
			}else if( !is.null(seed) ){
				warning("NMF::nmf - Discarding value of argument `seed`: directly using NMF model supplied in `rank` instead.\n"
						, "  If seeding is necessary, please use argument `model` pass initial model slots, which will be filled by the seeding method.")
			}
#			# pass the model via a one-off global variable
#			.nmf_InitModel(rank)
		}
		
		# replace missing method by NULL for correct dispatch
		if( missing(method) ) method <- NULL
		
		nmf(x, nbasis(rank), method, seed=rank, ...)
	}
)
.nmf_InitModel <- oneoffVariable()

#' Fits an NMF model using the NMF model supplied in \code{seed}, to seed the computation,
#' i.e. as a starting point.
#' 
#' This method is provided for completeness and is equivalent to 
#' \code{nmf(x, seed, method, ...)}.
#'   
setMethod('nmf', signature(x='matrix', rank='NULL', method='ANY'),
	function(x, rank, method, seed, ...){
		
		if( missing(seed) || !is.nmf(seed) )
			stop("NMF::nmf - Argument `seed` must be an NMF model when argument `rank` is missing.")
		
		# replace missing method by NULL for correct dispatch
		if( missing(method) ) method <- NULL
		
		nmf(x, nbasis(seed), method, seed=seed, ...)
	}
)
#' Method defined to ensure the correct dispatch to workhorse methods in case
#' of argument \code{rank} is missing.
setMethod('nmf', signature(x='matrix', rank='missing', method='ANY'),
	function(x, rank, method, ...){
		# replace missing method by NULL for correct dispatch
		if( missing(method) ) method <- NULL
		nmf(x, NULL, method, ...)
	}
)
#' Method defined to ensure the correct dispatch to workhorse methods in case
#' of argument \code{method} is missing.
#' 
#' @demo
#' # missing method: use default algorithm
#' res <- nmf(x, 3)
#' 
setMethod('nmf', signature(x='matrix', rank='numeric', method='missing'),
	function(x, rank, method, ...){
		nmf(x, rank, NULL, ...)
	}
)
#' Fits an NMF model partially seeding the computation with a given matrix passed 
#' in \code{rank}.
#' 
#' The matrix \code{rank} is used either as initial value for the basis or mixture 
#' coefficient matrix, depending on its dimension.
#' 
#' Currently, such partial NMF model is directly used as a seed, meaning that 
#' the remaining part is left uninitialised, which is not accepted by all NMF algorithm.
#' This should change in the future, where the missing part of the model will be 
#' drawn from some random distribution.
#' 
#' Amongst built-in algorithms, only \sQuote{snmf/l} and \sQuote{snmf/r} support 
#' partial seeds, with only the coefficient or basis matrix initialised 
#' respectively.
#' 
#' @demo
#' 
#' # Fit a 3-rank model providing an initial value for the basis matrix
#' nmf(x, rmatrix(nrow(x), 3), 'snmf/r')
#'  
#' # Fit a 3-rank model providing an initial value for the mixture coefficient matrix
#' nmf(x, rmatrix(3, ncol(x)), 'snmf/l')
#' 
setMethod('nmf', signature(x='matrix', rank='matrix', method='ANY'), 
	function(x, rank, method, seed, model=list(), ...)
	{
		if( is.character(model) )
			model <- list(model=model) 
		if( !is.list(model) )
			stop("nmf - Invalid argument `model`: must be NULL or a named list of initial values for slots in an NMF object.")
		if( !hasNames(model, all=TRUE) )
			stop("nmf - Invalid argument `model`: all elements must be named")
			
		# remove rank specification if necessary
		if( !is.null(model$rank) ){
			warning("nmf - Discarding rank specification in argument `model`: use value inferred from matrix supplied in argument `rank`")
			model$rank <- NULL
		}
		# check compatibility of dimensions
		newseed <- 
		if( nrow(rank) == nrow(x) ){
			# rank is the initial value for the basis vectors
			if( length(model)==0L ) nmfModel(W=rank)
			else{
				model$W <- rank
				do.call('nmfModel', model)
			}
		}else if( ncol(rank) == ncol(x) ){
            # rank is the initial value for the mixture coefficients
            if( length(model)==0L ) nmfModel(H=rank)
            else{
                model$H <- rank
                do.call('nmfModel', model)
            }
        }else
			stop("nmf - Invalid argument `rank`: matrix dimensions [",str_out(dim(x),sep=' x '),"]"
				, " are incompatible with the target matrix [", str_out(dim(x),sep=' x '),"].\n"
				, "  When `rank` is a matrix it must have the same number of rows or columns as the target matrix `x`.")
		
		# replace missing values by NULL values for correct dispatch
		if( missing(method) ) method <- NULL
		if( missing(seed) ) seed <- NULL
		#nmf(x, nbasis(newseed), method, seed=seed, model=newseed, ...)
		nmf(x, newseed, method, seed=seed, ...)
	}
)
#' Shortcut for \code{nmf(x, as.matrix(rank), method, ...)}.
setMethod('nmf', signature(x='matrix', rank='data.frame', method='ANY'),
	function(x, rank, method, ...){
		# replace missing values by NULL values for correct dispatch
		if( missing(method) ) method <- NULL
		
		nmf(x, as.matrix(rank), method, ...)
	}
)

#' This method implements the interface for fitting formula-based NMF models.
#' See \code{\link{nmfModel}}.
#' 
#' Argument \code{rank} target matrix or formula environment.
#' If not missing, \code{model} must be a \code{list}, a \code{data.frame} or 
#' an \code{environment} in which formula variables are searched for.    
#' 
setMethod('nmf', signature(x='formula', rank='ANY', method='ANY'),
	function(x, rank, method, ..., model=NULL){
		# replace missing values by NULL values for correct dispatch
		if( missing(method) ) method <- NULL
		if( missing(rank) ) rank <- NULL
		
		# if multiple numeric rank: use nmfRestimateRank
		if( is.vector(rank) && is.numeric(rank)  ){
			if( length(rank) > 1L ){
				return( nmfEstimateRank(x, rank, method, ..., model=model) )
			}
		}
		
		# build formula based model
		model <- nmfModel(x, rank, data=model)
		nmf(attr(model, 'target'), nbasis(model), method, ..., model=model)
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

#' Error Checks in NMF Runs
#' 
#' Auxiliary function for internal error checks in nmf results.
#' 
#' @param object a list of lists
#' @param element name of an element of the inner lists 
#' 
#' @keywords internal
checkErrors <- function(object, element=NULL){
	
	# extract error messages
	errors <- 
			if( is.null(element) ){
				lapply(seq_along(object), function(i){ 
							x <- object[[i]]
							if( is(x, 'error') ) c(i, x)
							else NA 
						})
			}else{
				lapply(seq_along(object), function(i){
							x <- object[[i]][[element, exact=TRUE]]
							if( is(x, 'error') ) c(i, x) 
							else NA
						})
			}
	errors <- errors[!is.na(errors)]
	nerrors <- length(errors)
	res <- list(n = nerrors)
	
	# format messages
	if( nerrors ){
		ierrors <- sapply(errors, '[[', 1L)
		msg <- sapply(errors, '[[', 2L)
		ierrors_unique <- ierrors[!duplicated(msg)]
		res$msg <- str_c("  - ", str_c("run #", ierrors_unique, ': ', msg[ierrors_unique], collapse="\n  - "))
	}
	
	# return error data
	res
}

###% Performs NMF on a matrix using a given NMF method.
###%
###% This method is the entry point for NMF. It is eventually called by any definition of the \code{nmf} function.

#' @param seed specification of the starting point or seeding method, which will 
#' compute a starting point, usually using data from the target matrix in order to 
#' provide a good guess. 
#' 
#' The seeding method may be specified in the following way:
#' 		
#' \describe{
#' 
#' \item{a \code{character} string:}{ giving the name of a \emph{registered} 
#' seeding method. The corresponding method will be called to compute 
#' the starting point. 
#' 
#' Available methods can be listed via \code{nmfSeed()}.
#' See its dedicated documentation for details on each available registered methods
#' (\code{\link{nmfSeed}}).
#' }
#' 		
#' \item{a \code{list}:}{ giving the name of a \emph{registered} 
#' seeding method and, optionally, extra parameters to pass to it.}
#' 
#' \item{a single \code{numeric}:}{ that is used to seed the random number 
#' generator, before generating a random starting point.
#' 
#' Note that when performing multiple runs, the L'Ecuyer's RNG is used in order to 
#' produce a sequence of random streams, that is used in way that ensures 
#' that parallel computation are fully reproducible.
#' }
#' 
#' \item{an object that inherits from \code{\linkS4class{NMF}}:}{ it should 
#' contain the	data of an initialised NMF model, i.e. it must contain valid 
#' basis and mixture coefficient matrices, directly usable by the algorithm's 
#' workhorse function.}
#' 
#' \item{a \code{function}:}{ that computes the starting point. It must have 
#' 	signature \code{(object="NMF", target="matrix", ...)} and return an object that 
#' inherits from class \code{NMF}. 
#' It is recommended to use argument \code{object} as a template for the returned object, 
#' by only updating the basis and coefficient matrices, using \code{\link{basis<-}} and 
#' \code{\link{coef<-}} respectively.
#' }
#' 
#' }   
#' 
#' @param rng rng specification for the run(s).
#' This argument should be used to set the the RNG seed, while still specifying the seeding 
#' method argument \var{seed}.
#' 
#' @param model specification of the type of NMF model to use.
#' 
#' It is used to instantiate the object that inherits from class \code{\linkS4class{NMF}}, 
#' that will be passed to the seeding method.
#' The following values are supported:
#' \itemize{
#'   	
#' \item \code{NULL}, the default model associated to the NMF algorithm is 
#'  instantiated and \code{...} is looked-up for arguments with names that 
#' 	correspond to slots in the model class, which are passed to the function 
#'  \code{\link{nmfModel}} to instantiate the model.
#' 	Arguments in \code{...} that do not correspond to slots are passed to the
#' 	algorithm.  
#' 
#' \item a single \code{character} string, that is the name of the NMF model 
#' 	class to be instantiate. 
#'  In this case, arguments in \code{...} are handled in the same way as 
#' 	when \code{model} is \code{NULL}.
#'   	 
#' \item a \code{list} that contains named values that are passed to the 
#'  function \code{\link{nmfModel}} to instantiate the model.
#'  In this case, \code{...} is not looked-up at all, and passed entirely to 
#'  the algorithm.
#'  This means that all necessary model parameters must be specified in 
#' 	\code{model}.
#' 
#'  }
#'  	
#'  \strong{Argument/slot conflicts:}
#'  In the case a parameter of the algorithm has the same name as a model slot, 
#'  then \code{model} MUST be a list -- possibly empty --, if one wants this 
#'  parameter to be effectively passed to the algorithm.
#'    	
#'  If a variable appears in both arguments \code{model} and \code{\dots}, 
#'  the former will be used to initialise the NMF model, the latter will be 
#'  passed to the NMF algorithm. 
#'  See code examples for an illustration of this situation.
#' 
#' @param nrun number of runs to perform. 
#' It specifies the number of runs to perform.
#' By default only one run is performed, except if \code{rank} is a numeric vector 
#' with more than one element, in which case a default of 30 runs per value of the 
#' rank are performed, allowing the computation of a consensus matrix that is used 
#' in selecting the appropriate rank (see \code{\link{consensus}}).
#' 
#' When using a random seeding method, multiple runs are generally required to 
#' achieve stability and avoid \emph{bad} local minima. 
#' 
#' @param .options this argument is used to set runtime options. 
#' 
#' It can be a \code{list} containing named options with their values, or, in 
#' the case only boolean/integer options need to be set, a character string 
#' that specifies which options are turned on/off or their value, in a unix-like 
#' command line argument way.
#'  
#' The string must be composed of characters that correspond to a given option 
#' (see mapping below), and modifiers '+' and '-' that toggle options on and off respectively. 
#' E.g. \code{.options='tv'} will toggle on options \code{track} and \code{verbose}, 
#' while \code{.options='t-v'} will toggle on option \code{track} and toggle off 
#' option \code{verbose}. 
#' 
#' Modifiers '+' and '-' apply to all option character found after them:
#' \code{t-vp+k} means \code{track=TRUE}, \code{verbose=parallel=FALSE}, 
#' and \code{keep.all=TRUE}.  
#' The default behaviour is to assume that \code{.options} starts with a '+'.
#' 
#' for options that accept integer values, the value may be appended to the 
#' option's character e.g. \code{'p4'} for asking for 4 processors or \code{'v3'} 
#' for showing verbosity message up to level 3.
#' 
#' The following options are available (the characters after \dQuote{-} are those 
#' to use to encode \code{.options} as a string):
#' \describe{
#' 
#' \item{debug - d}{ Toggle debug mode (default: \code{FALSE}). 
#' Like option \code{verbose} but with more information displayed.}
#' 
#' \item{keep.all - k}{ used when performing multiple runs (\code{nrun}>1): if
#' \code{TRUE}, all factorizations are saved and returned (default: \code{FALSE}).
#' Otherwise only the factorization achieving the minimum residuals is returned.}
#' 
#' \item{parallel - p}{ this option is useful on multicore *nix or Mac machine
#' only, when performing multiple runs (\code{nrun} > 1) (default: \code{TRUE}).  
#' If \code{TRUE}, the runs are performed using the parallel foreach backend 
#' defined in argument \code{.pbackend}. 
#' If this is set to \code{'mc'} or \code{'par'} then \code{nmf} tries to 
#' perform the runs using multiple cores with package 
#' \code{link[doParallel]{doParallel}} -- which therefore needs to be installed.
#' 
#' If equal to an integer, then \code{nmf} tries to perform the computation on 
#' the specified number of processors.
#' When passing options as a string the number is appended to the option's character 
#' e.g. \code{'p4'} for asking for 4 processors.
#' 
#' If \code{FALSE}, then the computation is performed sequentially using the base 
#' function \code{\link{sapply}}.
#' 
#' Unlike option 'P' (capital 'P'), if the computation cannot be performed in
#' parallel, then it will still be carried on sequentially.
#' 
#' \strong{IMPORTANT NOTE FOR MAC OS X USERS:} The parallel computation is
#' based on the \code{doMC} and \code{multicore} packages, so the same care
#' should be taken as stated in the vignette of \code{doMC}: \emph{\dQuote{it
#' is not safe to use doMC from R.app on Mac OS X. Instead, you should use doMC
#' from a terminal session, starting R from the command line.}} }
#' 
#' \item{parallel.required - P}{ Same as \code{p}, but an error is thrown if
#' the computation cannot be performed in parallel or with the specified number 
#' of processors.}
#' 
#' \item{shared.memory - m}{ toggle usage of shared memory (requires the 
#' \pkg{synchronicity} package).
#' Default is as defined by \code{nmf.getOption('shared.memory')}.}
#' 
#' \item{restore.seed - r}{ deprecated option since version 0.5.99.
#' Will throw a warning if used.}
#' 
#' \item{simplifyCB - S}{ toggle simplification of the callback results. 
#' Default is \code{TRUE}}
#' 
#' \item{track - t}{ enables error tracking (default: FALSE).
#' If \code{TRUE}, the returned object's slot \code{residuals} contains the 
#' trajectory of the objective values, which can be retrieved via 
#' \code{residuals(res, track=TRUE)}
#' This tracking functionality is available for all built-in algorithms.
#' }
#' 
#' \item{verbose - v}{ Toggle verbosity (default: \code{FALSE}). 
#' If \code{TRUE}, messages about the configuration and the state of the 
#' current run(s) are displayed.
#' The level of verbosity may be specified with an integer value, the greater 
#' the level the more messages are displayed.
#' Value \code{FALSE} means no messages are displayed, while value \code{TRUE} 
#' is equivalent to verbosity level 1.  
#' }
#' 
#' }
#' 
#' @param .pbackend specification of the \code{\link{foreach}} parallel backend 
#' to register and/or use when running in parallel mode. 
#' See options \code{p} and \code{P} in argument \code{.options} for how to 
#' enable this mode.
#' Note that any backend that is internally registered is cleaned-up on exit, 
#' so that the calling foreach environment should not be affected by a call to 
#' \code{nmf} -- except when \code{.pbackend=NULL}. 
#' 
#' Currently it accepts the following values: 
#' \describe{
#' 
#' \item{\sQuote{par}}{ use the backend(s) defined by the package 
#' \code{\link{doParallel}};}
#' \item{a numeric value}{ use the specified number of cores with \code{doParallel}
#' backend;}
#' \item{\sQuote{seq}}{ use the foreach sequential backend \code{doSEQ};}
#' \item{\code{NULL}}{ use currently registered backend;} 
#' \item{\code{NA}}{ do not compute using a foreach loop -- and therefore not in
#'  parallel --  but rather use a call to standard \code{\link{sapply}}.
#' This is useful for when developing/debugging NMF algorithms, as foreach loop
#' handling may sometime get in the way.
#'   
#' Note that this is equivalent to using \code{.options='-p'} or \code{.options='p0'}, 
#' but takes precedence over any option specified in \code{.options}: 
#' e.g. \code{nmf(..., .options='P10', .pbackend=NA)} performs all runs sequentially 
#' using \code{sapply}.
#' Use \code{nmf.options(pbackend=NA)} to completely disable foreach/parallel computations 
#' for all subsequent \code{nmf} calls.}
#' 
#' \item{\sQuote{mc}}{ identical to \sQuote{par} and defined to ensure backward 
#' compatibility.}
#' }
#' 
#' @param .callback Used when option \code{keep.all=FALSE} (default).  It
#' allows to pass a callback function that is called after each run when
#' performing multiple runs (i.e. with \code{nrun>1}).  
#' This is useful for example if one is also interested in saving summary 
#' measures or process the result of each NMF fit before it gets discarded. 
#' After each run, the callback function is called with two arguments, the
#' \code{\linkS4class{NMFfit}} object that as just been fitted and the run 
#' number: \code{.callback(res, i)}.
#' For convenience, a function that takes only one argument or has 
#' signature \code{(x, ...)} can still be passed in \code{.callback}.
#' It is wrapped internally into a dummy function with two arguments, 
#' only the first of which is passed to the actual callback function (see example 
#' with \code{summary}).
#' 
#' The call is wrapped into a tryCatch so that callback errors do not stop the 
#' whole computation (see below).
#' 
#' The results of the different calls to the callback function are stored in a
#' miscellaneous slot accessible using the method \code{$} for \code{NMFfit} 
#' objects: \code{res$.callback}.
#' By default \code{nmf} tries to simplify the list of callback result using 
#' \code{sapply}, unless option \code{'simplifyCB'} is \code{FASE}.
#' 
#' If no error occurs \code{res$.callback} contains the list of values that 
#' resulted from the calling the callback function --, ordered as the fits.
#' If any error occurs in one of the callback calls, then the whole computation is 
#' \strong{not} stopped, but the error message is stored in \code{res$.callback}, 
#' in place of the result.
#' 
#' See the examples for sample code.
#' 
#' @return The returned value depends on the run mode:
#' 
#' \item{Single run:}{An object of class \code{\linkS4class{NMFfit}}.} 
#' 
#' \item{Multiple runs, single method:}{When \code{nrun > 1} and \code{method} 
#' is not \code{list}, this method returns an object of class \code{\linkS4class{NMFfitX}}.} 
#' 
#' \item{Multiple runs, multiple methods:}{When \code{nrun > 1} and \code{method} 
#' is a \code{list}, this method returns an object of class \code{\linkS4class{NMFList}}.}
#' 
#' @demo
#' 
#' # default fit
#' res <- nmf(x, 2)
#' summary(res, class=groups)
#' 
#' # run default algorithm multiple times (only keep the best fit)
#' res <- nmf(x, 3, nrun=10)
#' res
#' summary(res, class=groups)
#' 
#' # run default algorithm multiple times keeping all the fits
#' res <- nmf(x, 3, nrun=10, .options='k')
#' res
#' summary(res, class=groups)
#' 
#' ## Note: one could have equivalently done
#' # res <- nmf(V, 3, nrun=10, .options=list(keep.all=TRUE))
#'  
#' # use a method that fit different model
#' res <- nmf(x, 2, 'nsNMF')
#' fit(res)
#' 
#' # pass parameter theta to the model via `...`
#' res <- nmf(x, 2, 'nsNMF', theta=0.2)
#' fit(res)
#' 
#' ## handling arguments in `...` and model parameters
#' myfun <- function(x, start, theta=100){ cat("theta in myfun=", theta, "\n\n"); start }
#' # no conflict: default theta
#' fit( nmf(x, 2, myfun) ) 
#' # no conlfict: theta is passed to the algorithm
#' fit( nmf(x, 2, myfun, theta=1) )  
#' # conflict: theta is used as model parameter
#' fit( nmf(x, 2, myfun, model='NMFns', theta=0.1) )
#' # conflict solved: can pass different theta to model and algorithm
#' fit( nmf(x, 2, myfun, model=list('NMFns', theta=0.1), theta=5) )
#' 
#' ## USING SEEDING METHODS
#' 
#' # run default algorithm with the Non-negative Double SVD seeding method ('nndsvd')
#' res <- nmf(x, 3, seed='nndsvd')
#' 
#' ## Note: partial match also works
#' identical(res, nmf(x, 3, seed='nn'))
#' 
#' # run nsNMF algorithm, fixing the seed of the random number generator 
#' res <- nmf(x, 3, 'nsNMF', seed=123456)
#' nmf.equal(nmf(x, 3, 'nsNMF', seed=123456), res)
#' 
#' # run default algorithm specifying the starting point following the NMF standard model
#' start.std <- nmfModel(W=matrix(0.5, n, 3), H=matrix(0.2, 3, p))   
#' nmf(x, start.std)
#' 
#' # to run nsNMF algorithm with an explicit starting point, this one
#' # needs to follow the 'NMFns' model:
#' start.ns <- nmfModel(model='NMFns', W=matrix(0.5, n, 3), H=matrix(0.2, 3, p))   
#' nmf(x, start.ns)
#' # Note: the method name does not need to be specified as it is infered from the 
#' # when there is only one algorithm defined for the model.
#' 
#' # if the model is not appropriate (as defined by the algorihtm) an error is thrown 
#' # [cf. the standard model doesn't include a smoothing parameter used in nsNMF] 
#' try( nmf(x, start.std, method='nsNMF') )
#' 
#' ## Callback functions
#' # Pass a callback function to only save summary measure of each run
#' res <- nmf(x, 3, nrun=3, .callback=summary)
#' # the callback results are simplified into a matrix
#' res$.callback
#' res <- nmf(x, 3, nrun=3, .callback=summary, .opt='-S')
#' # the callback results are simplified into a matrix
#' res$.callback
#' 
#' # Pass a custom callback function
#' cb <- function(obj, i){ if( i %% 2 ) sparseness(obj) >= 0.5 }
#' res <- nmf(x, 3, nrun=3, .callback=cb)
#' res$.callback
#' 
#' # Passs a callback function which throws an error
#' cb <- function(){ i<-0; function(object){ i <<- i+1; if( i == 1 ) stop('SOME BIG ERROR'); summary(object) }}
#' res <- nmf(x, 3, nrun=3, .callback=cb())
#' 
#' ## PARALLEL COMPUTATIONS
#' # try using 3 cores, but use sequential if not possible 
#' res <- nmf(x, 3, nrun=3, .options='p3')
#' 
#' # force using 3 cores, error if not possible
#' res <- nmf(x, 3, nrun=3, .options='P3')
#' 
#' # use externally defined cluster
#' library(parallel)
#' cl <- makeCluster(6)
#' res <- nmf(x, 3, nrun=3, .pbackend=cl)
#' 
#' # use externally registered backend
#' registerDoParallel(cl)
#' res <- nmf(x, 3, nrun=3, .pbackend=NULL)
#' 
setMethod('nmf', signature(x='matrix', rank='numeric', method='NMFStrategy'),
#function(x, rank, method, seed='random', nrun=1, keep.all=FALSE, optimized=TRUE, init='NMF', track, verbose, ...)
function(x, rank, method
		, seed=nmf.getOption('default.seed'), rng = NULL
		, nrun=if( length(rank) > 1L ) 30 else 1, model=NULL, .options=list()
		, .pbackend=nmf.getOption('pbackend')
		, .callback=NULL #callback function called after a run  
		, ...)
{
	fwarning <- function(...) nmf_warning('nmf', ...)
	fstop <- function(...) nmf_stop('nmf', ...)
	
	# if options are given as a character string, translate it into a list of booleans
	if( is.character(.options) ){
		.options <- .translate.string(.options, 
				c(t='track', v='verbose', d='debug'
				, p='parallel', P='parallel.required'
				, k='keep.all', r='restore.seed', f='dry.run'
				, g='garbage.collect'
				, c='cleanup', S='simplifyCB'
				, R='RNGstream', m='shared.memory'))
	}
	
	# get seeding method from the strategy's defaults if needed
	seed <- defaultArgument(seed, method, nmf.getOption('default.seed'), force=is.null(seed))
	.method_defaults <- method@defaults
	.method_defaults$seed <- NULL
	#
	# RNG specification
	if( isRNGseed(seed) ){
		if( !is.null(rng) )
			warning("Discarding RNG specification in argument `rng`: using those passed in argument `seed`.")
		rng <- seed
		seed <- 'random'
	}
	#

	# setup verbosity options
	debug <- if( !is.null(.options$debug) ) .options$debug else nmf.getOption('debug')
	verbose <- if( debug ) Inf
				else if( !is.null(.options$verbose) ) .options$verbose
				else nmf.getOption('verbose')
	
	# show call in debug mode
	if( debug ){
		.ca <- match.call()
		message('# NMF call: ', paste(capture.output(print(.ca)), collapse="\n  "))
	}
	# nmf over a range of values: pass the call to nmfEstimateRank
	if( length(rank) > 1 ){
		if( verbose <= 1 )
			.options$verbose <- FALSE
		return( nmfEstimateRank(x, range = rank, method = method, nrun = nrun
								, seed = seed, rng = rng, model = model
								, .pbackend = .pbackend, .callback = .callback
								, verbose=verbose, .options=.options, ...) )
	}
	
	.OPTIONS <- list()
	# cleanup on exit
	.CLEANUP <- .options$cleanup %||% TRUE
	
	# tracking of objective value
	.OPTIONS$track <- if( !is.null(.options$track) ) .options$track 
					else nmf.getOption('track')
	# dry run
	dry.run <- .options$dry.run %||% FALSE 
	# call the garbage collector regularly
	opt.gc <- if( !is.null(.options$garbage.collect) ) .options$garbage.collect
			  else nmf.getOption('gc')
	if( is.logical(opt.gc) && opt.gc )
		opt.gc <- ceiling(max(nrun,50) / 3)
	.options$garbage.collect <- opt.gc
	
	# keep results from all runs?
	keep.all <- .options$keep.all %||% FALSE
    # shared memory?
    shared.memory <- if( !is.null(.options$shared.memory) ) .options$shared.memory else nmf.getOption('shared.memory')
	# use RNG stream
	.options$RNGstream <- .options$RNGstream %||% TRUE
	
	# discard .callback when not used
	if( is.function(.callback) ){
		w <- if( nrun==1 ) "discarding argument `.callback`: not used when `nrun=1`."
			else if( keep.all )	
				"discarding argument `.callback`: not used when option `keep.all=TRUE`."
		if( !is.null(w) ){
			.callback <- NULL
			fwarning(w, immediate.=TRUE)
		}
		
		# wrap into another function if necessary
		if( is.function(.callback) ){
			# default is to simplify
			.options$simplifyCB <- .options$simplifyCB %||% TRUE 
			args <- formals(.callback)
			if( length(args) <= 2L ){
				if( length(args) < 2L || '...' %in% names(args) ){
					.CALLBACK <- .callback
					.callback <- function(object, i) .CALLBACK(object)
				}
			}
			
			# define post-processing function
			processCallback <- function(res){
				# check errors
				errors <- checkErrors(res, '.callback')
				if( errors$n > 0 ){
					fwarning("All NMF fits were successful but ", errors$n, "/", nrun, " callback call(s) threw an error.\n"
							,"# ", if(errors$n>10) "First 10 c" else "C", "allback error(s) thrown:\n"
							, errors$msg
					)
				}
				# add callback values to result list
				sapply(res, '[[', '.callback'
					, simplify=.options$simplifyCB && errors$n == 0L)
			}
		}
	}
	
	## ROLLBACK PROCEDURE
	exitSuccess <- exitCheck()
	on.exit({ 
		if( verbose > 1 ) message("# NMF computation exit status ... ", if( exitSuccess() ) 'OK' else 'ERROR')
		if( verbose > 2 ){
			if( exitSuccess() ){
				message('\n## Running normal exit clean up ... ')
			}else{ 
				message('\n## Running rollback clean up ... ')
			}
		}
	}, add=TRUE)
	# RNG restoration on error
	.RNG_ORIGIN <- getRNG()
	on.exit({
		if( !exitSuccess() ){
			if( verbose > 2 ) message("# Restoring RNG settings ... ", appendLF=verbose>3)
			setRNG(.RNG_ORIGIN)
			if( verbose > 3 ) showRNG(indent=' #')
			if( verbose > 2 ) message("OK")
		}
	}, add=TRUE)

	# Set debug/verbosity option just for the time of the run
	old.opt <- nmf.options(debug=debug, verbose=verbose, shared.memory = shared.memory);
	on.exit({
		if( verbose > 2 ) message("# Restoring NMF options ... ", appendLF=FALSE)
		nmf.options(old.opt)
		if( verbose > 2 ) message("OK")
	}, add=TRUE)
	
	# make sure rank is an integer
	rank <- as.integer(rank)
	if( length(rank) != 1 ) fstop("invalid argument 'rank': must be a single numeric value")
	if( rank < 1 ) fstop("invalid argument 'rank': must be greater than 0")
	
	# option 'restore.seed' is deprecated
	if( !is.null(.options$restore.seed) )
		fwarning("Option 'restore.seed' is deprecated and discarded since version 0.5.99.")
	
	if( verbose ){
		if( dry.run ) message("*** fake/dry-run ***")
		message("NMF algorithm: '", name(method), "'")
	}
	
	##START_MULTI_RUN
	# if the number of run is more than 1, then call itself recursively
	if( nrun > 1 )
	{
		if( verbose ) message("Multiple runs: ", nrun)
		
		if( verbose > 3 ){
			cat("## OPTIONS:\n")
			sapply(seq_along(.options)
					, function(i){
						r <- i %% 4
						cat(if(r!=1) '\t| ' else "# ", names(.options)[i],': ', .options[[i]], sep='')
						if(r==0) cat("\n# ")
					})
			if( length(.options) %% 4 != 0 )cat("\n")
		}
				
		## OPTIONS: parallel computations 
		# option require-parallel: parallel computation is required if TRUE or numeric != 0
		opt.parallel.required <- !is.null(.options$parallel.required) && .options$parallel.required
		# determine specification for parallel computations 
		opt.parallel.spec <- 
				if( opt.parallel.required ){ # priority over try-parallel 
					# option require-parallel implies and takes precedence over option try-parallel
					.options$parallel.required
				}else if( !is.null(.options$parallel) ) .options$parallel # priority over .pbackend
				else !is_NA(.pbackend) # required only if backend is not trivial

		# determine if one should run in parallel at all: TRUE or numeric != 0, .pbackend not NA
		opt.parallel <- !is_NA(.pbackend) && (isTRUE(opt.parallel.spec) || opt.parallel.spec)
		##
		if( opt.parallel ){
			if( verbose > 1 )
				message("# Setting up requested `foreach` environment: "
						, if( opt.parallel.required ) 'require-parallel' else 'try-parallel'
						, ' [', quick_str(.pbackend) , ']')
			
			
			# switch doMC backend to doParallel
			if( isString(.pbackend, 'MC', ignore.case=TRUE) ){ 
				.pbackend <- 'par'
			}
			# try setting up parallel foreach backend
			oldBackend <- setupBackend(opt.parallel.spec, .pbackend, !opt.parallel.required, verbose=verbose)
			opt.parallel <- !isFALSE(oldBackend)
			# setup backend restoration if using one different from the current one
			if( opt.parallel && !is_NA(oldBackend) ){
				on.exit({
						if( verbose > 2 ){
							message("# Restoring previous foreach backend '", getDoBackendName(oldBackend) ,"' ... ", appendLF=FALSE)
						}
						setDoBackend(oldBackend, cleanup=TRUE)
						if( verbose > 2 ) message('OK')
					}, add=TRUE)
			}#
			
			# From this point, the backend is registered
			# => one knows if we'll run a sequential or parallel foreach loop
			.MODE_SEQ <- is.doSEQ()
			MODE_PAR <- .MODE_PAR <- !.MODE_SEQ
			
		}
		
		# check seed method: fixed values are not sensible -> warning
		.checkRandomness <- FALSE
		if( is.nmf(seed) && !is.empty.nmf(seed) ){
			.checkRandomness <- TRUE
		}
		# start_RNG_all
		# if the seed is numerical or a rstream object,	then use it to set the 
		# initial state of the random number generator:			
		# build a sequence of RNGstreams: if no suitable seed is provided
		# then the sequence use a random seed generated with a single draw 
		# of the current active RNG. If the seed is valid, then the 
		# 
		# setup the RNG sequence

		# override with standard RNG if .options$RNGstream=FALSE
		resetRNG <- NULL
		if( !.options$RNGstream && (!opt.parallel || .MODE_SEQ) ){
			
			.RNG.seed <- rep(list(NULL), nrun)
			if( isNumber(rng) ){
				resetRNG <- getRNG()
				if( verbose > 2 ) message("# Force using current RNG settings seeded with: ", rng)
				set.seed(rng)
			}else if( verbose > 2 ) 
				message("# Force using current RNG settings")
			
		}else{
			.RNG.seed <- setupRNG(rng, n = nrun, verbose=verbose)
			# restore the RNG state on exit as after RNGseq:
			# - if no seeding occured then the RNG has still been drawn once in RNGseq
			# which must be reflected so that different unseeded calls use different RNG states
			# - one needs to restore the RNG because it switched to L'Ecuyer-CMRG. 
			resetRNG <- getRNG()
		}
		stopifnot( length(.RNG.seed) == nrun )
		# update RNG settings on exit if necessary
		# and only if no error occured
		if( !is.null(resetRNG) ){
			on.exit({
				if( exitSuccess() ){
					if( verbose > 2 ) message("# Updating RNG settings ... ", appendLF=FALSE)							
					setRNG(resetRNG)
					if( verbose > 2 ) message("OK")
					if( verbose > 3 ) showRNG()
				}
			}, add=TRUE)
		}
		#end_RNG_all
		
		####FOREACH_NMF
		if( opt.parallel ){
			
			if( verbose ){
				if( verbose > 1 )
						message("# Using foreach backend: ", getDoParName()
								," [version ", getDoParVersion(),"]")
				# show number of processes
				if( getDoParWorkers() == 1 ) message("Mode: sequential [foreach:",getDoParName(),"]")
				else message("Mode: parallel ", str_c("(", getDoParWorkers(), '/', parallel::detectCores()," core(s))"))
			}
			
			# check shared memory capability
			.MODE_SHARED <- !keep.all && setupSharedMemory(verbose)
			
			# setup temporary directory when not keeping all fits
			if( !keep.all || verbose ){
				NMF_TMPDIR <- setupTempDirectory(verbose)
				# delete on exit
				if( .CLEANUP ){
					on.exit({
						if( verbose > 2 ) message("# Deleting temporary directory '", NMF_TMPDIR, "' ... ", appendLF=FALSE)
						unlink(NMF_TMPDIR, recursive=TRUE)
						if( verbose > 2 ) message('OK')
					}, add=TRUE)
				}
			}
			
			run.all <- function(x, rank, method, seed, model, .options, ...){
								
				## 1. SETUP
				# load some variables from parent environment to ensure they 
				# are exported in the foreach loop
				MODE_SEQ <- .MODE_SEQ
				MODE_SHARED <- .MODE_SHARED
				verbose <- verbose
				keep.all <- keep.all
				opt.gc <- .options$garbage.collect
				CALLBACK <- .callback
				.checkRandomness <- .checkRandomness
				
				# check if single or multiple host(s)
				hosts <- unique(getDoParHosts())
				if( verbose > 2 ) message("# Running on ", length(hosts), " host(s): ", str_out(hosts))
				SINGLE_HOST <- length(hosts) <= 1L
				MODE_SHARED <- MODE_SHARED && SINGLE_HOST
				if( verbose > 2 ) message("# Using shared memory ... ", MODE_SHARED)
				
				# setup mutex evaluation function
				mutex_eval <- if( MODE_SHARED ) ts_eval(verbose = verbose > 4) else force
				
				# Specific thing only if one wants only the best result
				if( !keep.all ){ 
					NMF_TMPDIR <- NMF_TMPDIR
					# - Define the shared memory objects
					vOBJECTIVE <- gVariable(as.numeric(NA), MODE_SHARED) 
					# the consensus matrix is computed only if not all the results are kept				
					vCONSENSUS <- gVariable(matrix(0, ncol(x), ncol(x)), MODE_SHARED)			
				}
				
				## 2. RUN
				# ensure that the package NMF is in each worker's search path
				.packages <- setupLibPaths('NMF', verbose>3)
                
                # export all packages that contribute to NMF registries, 
                # e.g., algorithms or seeding methods.
                # This is important so that these can be found in worker nodes
                # for non-fork clusters.
                if( !is.null(contribs <- registryContributors(package = 'NMF')) ){
                    .packages <- c(.packages, contribs)
                }
                
				# export dev environment if in dev mode 
#				.export <- if( isDevNamespace('NMF') && !is.doSEQ() ) ls(asNamespace('NMF'))
				
				# in parallel mode: verbose message from each run are only shown in debug mode
				.options$verbose <- FALSE 
				if( verbose ){
					if( debug || (.MODE_SEQ && verbose > 1) )
						.options$verbose <- verbose
					
					if( (!.MODE_SEQ && !debug) || (.MODE_SEQ && verbose == 1) ){
						if( verbose == 1 ){
							# create progress bar
							pbar <- txtProgressBar(0, nrun+1, width=50, style=3, title='Runs:'
													, shared=NMF_TMPDIR)
						}else{
							cat("Runs: ")
						}
					}
				}
				
				# get options from master process to pass to workers
				nmf.opts <- nmf.options()
				
				# load extra required packages for shared mode 
				if( MODE_SHARED ) 
					.packages <- c(.packages, 'bigmemory', 'synchronicity')
				
				res.runs <- foreach(n=1:nrun
								, RNGobj = .RNG.seed
								, .verbose = debug
								, .errorhandling = 'pass'
								, .packages = .packages
#								, .export = .export
#								, .options.RNG=.RNG.seed
								) %dopar% { #START_FOREACH_LOOP
				
					# Pass options from master process
					nmf.options(nmf.opts)
					
					# in mode sequential or debug: show details for each run
					if( MODE_SEQ && verbose > 1 )
						cat("\n## Run: ",n, "/", nrun, "\n", sep='')
					
					# set the RNG if necessary and restore after each run 
					if( MODE_SEQ && verbose > 2 )
							message("# Setting up loop RNG ... ", appendLF=FALSE)
					setRNG(RNGobj, verbose=verbose>3 && MODE_SEQ)
					if( MODE_SEQ && verbose > 2 )
							message("OK")
						
					# limited verbosity in simple mode
					if( verbose && !(MODE_SEQ && verbose > 1)){
						if( verbose >= 2 ) mutex_eval( cat('', n) )		
						else{
							# update progress bar (in mutex)
							mutex_eval(setTxtProgressBar(pbar, n))
							#
						}
					}
					
					# check RNG changes
					if( n == 1 && .checkRandomness ){
						.RNGinit <- getRNG()
					}
					
					# fit a single NMF model
					res <- nmf(x, rank, method, nrun=1, seed=seed, model=model, .options=.options, ...)
					
					if( n==1 && .checkRandomness && rng.equal(.RNGinit) ){
						warning("NMF::nmf - You are running multiple non-random NMF runs with a fixed seed")
					}
					
					# if only the best fit must be kept then update the shared objects
					if( !keep.all ){
						
						# initialise result list
						resList <- list(filename=NA, residuals=NA, .callback=NULL)
						
						##LOCK_MUTEX
						mutex_eval({
												
							# check if the run found a better fit
							.STATIC.err <- vOBJECTIVE()
							
							# retrieve approximation error
							err <- deviance(res)
							
							if( is.na(.STATIC.err) || err < .STATIC.err ){
								
								if( n>1 && verbose ){
									if( MODE_SEQ && verbose > 1 ) cat("## Better fit found [err=", err, "]\n")
									else if( verbose >= 2 ) cat('*')
								}
								
								# update residuals
								vOBJECTIVE(err)
								
								# update best fit on disk: use pid if not using shared memory
								resfile <- hostfile("fit", tmpdir=NMF_TMPDIR, fileext='.rds', pid=!MODE_SHARED)
								if( MODE_SEQ && verbose > 2 )
									message("# Serializing fit object in '", resfile, "' ... ", appendLF=FALSE)
								saveRDS(res, file=resfile, compress=FALSE)
								if( MODE_SEQ && verbose > 2 ){
									message(if( file.exists(resfile) ) 'OK' else 'ERROR')
								}
								# store the filename and achieved objective value in the result list
								resList$filename <- resfile
								resList$residuals <- err
																
							}
							
							## CONSENSUS
							# update the consensus matrix
							if( MODE_SHARED && SINGLE_HOST ){								
								# on single host: shared memory already contains consensus
								vCONSENSUS(vCONSENSUS() + connectivity(res, no.attrib=TRUE))
							}else{
								# on multiple hosts: must return connectivity and aggregate at the end
								resList$connectivity <- connectivity(res, no.attrib=TRUE)
							}
							
							## CALLBACK
							# call the callback function if necessary (return error as well)
							if( is.function(CALLBACK) ){
								resList$.callback <- tryCatch(CALLBACK(res, n), error=function(e) e)
							}
						
						})
						##END_LOCK_MUTEX
									
						# discard result object
						res <- NULL
						# return description list 
						res <- resList
					}
										
					# garbage collection if requested
					if( opt.gc && n %% opt.gc == 0 ){
						if( verbose > 2 ){
							if( MODE_SEQ )
								message("# Call garbage collector")
							else{
								mutex_eval( cat('%') )
							}
						}
						
						gc(verbose= MODE_SEQ && verbose > 3)
					}
					
					# return the result
					res
				}				
				## END_FOREACH_LOOP
				
				if( verbose && !debug ){
					if( verbose >= 2 ) cat(" ... DONE\n")
					else{
						setTxtProgressBar(pbar, nrun+1)
						pbar$kill(.CLEANUP)
					}
				}
				
				## 3. CHECK FIT ERRORS
				errors <- checkErrors(res.runs)
				if( errors$n > 0 ){
					fstop(errors$n,"/", nrun, " fit(s) threw an error.\n"
							,"# Error(s) thrown:\n", errors$msg)
				}
				
				## 4. WRAP UP
				if( keep.all ){ # result is a list of fits
					# directly return the list of fits
					res <- res.runs
					
				}else{ # result is a list of lists: filename, .callback 
					# loop over the result files to find the best fit
					if( verbose > 2 ) message("# Processing partial results ... ", appendLF=FALSE)
					ffstop <- function(...){ message('ERROR'); fstop(...) }
					# get best fit index
					idx <- which.min(sapply(res.runs, '[[', 'residuals'))
					if( length(idx) == 0L )
						ffstop("Unexpected error: no partial result seem to have been saved.")
					resfile <- res.runs[[idx]]$filename
					# check existence of the result file
					if( !file_test('-f', resfile) )
						ffstop("could not find temporary result file '", resfile, "'")
						
					# update res with a better fit
					res <- readRDS(resfile)
					if( !isNMFfit(res) ) 
						ffstop("invalid object found in result file '", resfile, "'")
					if( verbose > 2 ) message('OK')
					# wrap the result in a list: fit + consensus
					res <- list(fit=res, consensus=NA)
					
					# CONSENSUS MATRIX
					if( !is.null(res.runs[[1]]$connectivity) ){ # not MODE_SHARED
						# aggregate connectivity matrices
						con <- matrix(0, ncol(x), ncol(x))
						sapply(res.runs, function(x){
							con <<- con + x$connectivity 
						})
						res$consensus <- con
						
					}else{ # in MODE_SHARED: get consensus from global shared variable
						res$consensus <- vCONSENSUS()
						cn <- colnames(x) 
						if( is.null(cn) ) dimnames(res$consensus) <- NULL
						else dimnames(res$consensus) <- list(cn, cn)
					}
					
					# CALLBACKS
					if( !is.null(.callback) ){
						res$.callback <- processCallback(res.runs)
					}
				}
				##
				
				if( MODE_SEQ && verbose>1 ) cat("## DONE\n")
				
				# return result
				res
			}			
		}####END_FOREACH_NMF
		else{####SAPPLY_NMF
			
			run.all <- function(x, rank, method, seed, model, .options, ...){
				
				# by default force no verbosity from the runs
				.options$verbose <- FALSE
				if( verbose ){
					message("Mode: sequential [sapply]")
					if( verbose > 1 ){					
						# pass verbosity options in this case
						.options$verbose <- verbose
					}
				}
				
				## 1. SETUP				
				# define static variables for the case one only wants the best result
				if( !keep.all ){
					# statis list with best result: fit, residual, consensus
					best.static <- list(fit=NULL, residuals=NA, consensus=matrix(0, ncol(x), ncol(x)))					
				}
											
				## 2. RUN:
				# perform a single run `nrun` times
				if( verbose == 2 ){
					showRNG()
				}
				if( verbose && !debug ) cat('Runs:')
				res.runs <- mapply(1:nrun, .RNG.seed, FUN=function(n, RNGobj){
					
					#start_verbose
					if( verbose ){
						# in mode verbose > 1: show details for each run
						if( verbose > 1 ){
							cat("\n## Run: ",n, "/", nrun, "\n", sep='')							
						}else{
						# otherwise only some details for the first run
							cat('', n)
						}
					}#end_verbose
					
					# set the RNG for each run
					if( verbose > 2 ) message("# Setting up loop RNG ... ", appendLF=FALSE)
					setRNG(RNGobj, verbose=verbose>3)
					if( verbose > 2 ) message("OK")
					
					# check RNG changes
					if( n == 1 && .checkRandomness ){
						.RNGinit <- getRNG()
					}
					
					# fit a single NMF model
					res <- nmf(x, rank, method, nrun=1, seed=seed, model=model, .options=.options, ...)
					
					if( n==1 && .checkRandomness && rng.equal(.RNGinit) ){
						warning("NMF::nmf - You are running multiple non-random NMF runs with a fixed seed"
								, immediate.=TRUE)
					}
					
					if( !keep.all ){
						
						# initialise result list
						resList <- list(residuals=NA, .callback=NULL)
						
						# check if the run found a better fit
						err <- residuals(res)
						best <- best.static$residuals
						if( is.na(best) || err < best ){
							if( verbose ){
								if( verbose > 1L ) cat("## Updating best fit [deviance =", err, "]\n", sep='')
								else cat('*')
							}
							
							# update best fit (only if necessary)
							best.static$fit <<- res
							best.static$residuals <<- err
							
							resList$residuals <- err
						}
							
						# update the static consensus matrix (only if necessary)
						best.static$consensus <<- best.static$consensus + connectivity(res, no.attrib=TRUE)
						
						# call the callback function if necessary
						if( !is.null(.callback) ){
							resList$.callback <- tryCatch(.callback(res, n), error=function(e) e)							
						}
						
						# reset the result to NULL
						res <- resList
						
					}
					
					# garbage collection if requested
					if( opt.gc && n %% opt.gc == 0 ){
						if( verbose > 1 )
							message("# Call garbage collection NOW")
						else if( verbose )
							cat('%')
						
						gc(verbose = verbose > 3)
					}
					
					if( verbose > 1 ) cat("## DONE\n")
					
					# return the result
					res
				}, SIMPLIFY=FALSE)
				##
				
				if( verbose && !debug ) cat(" ... DONE\n")
				
				## 3. ERROR CHECK / WRAP UP
				
				if( keep.all ){
					res <- res.runs
				}else{
					res <- list(fit=best.static$fit, consensus=best.static$consensus)
					
					# CALLBACKS
					if( !is.null(.callback) ){
						res$.callback <- processCallback(res.runs)
					}
				}
				
				res
			}
			
		}####END_SAPPLY_NMF
			
		####END_DEFINE_RUN		
		
		# perform all the NMF runs
		t <- system.time({res <- run.all(x=x, rank=rank, method=method, seed=seed, model=model, .options, ...)})
		if( verbose && !debug ){
			cat("System time:\n")
			print(t)
		}
		
		if( keep.all ){
			
			# when keeping all the fits: join the results into an NMFfitXn object
			# TODO: improve memory management here
			res <- NMFfitX(res, runtime.all=t)
			
			return( exitSuccess(res) )
			
		}else{# if one just want the best result only return the best
			# ASSERT the presence of the result
			stopifnot( !is.null(res$fit) )
			# ASSERT the presence of the consensus matrix
			stopifnot( !is.null(res$consensus) )
			
			res.final <- NMFfitX(res$fit, consensus=res$consensus/nrun
								, runtime.all=t, nrun=as.integer(nrun)
								, rng1=.RNG.seed[[1]])
			
			# ASSERT and add callback if necessary 
			if( !is.null(.callback) ){
				stopifnot( !is.null(res$.callback) )
				res.final$.callback <- res$.callback
			}
			
			return( exitSuccess(res.final) )
		}
		
	}##END_MULTI_RUN
	
	# start_RNG
	# show original RNG settings in verbose > 2
	if( verbose > 3 ){
		message("# ** Current RNG settings:")
		showRNG()
	}
	
	# do something if the RNG was actually changed
	newRNG <- getRNG()
	.RNG.seed <- setupRNG(rng, 1, verbose=verbose-1)
	# setup restoration
	if( isRNGseed(rng) ){
		if( verbose > 3 ) showRNG()
		
		# restore RNG settings
		on.exit({
			if( verbose > 2 ) message("# Restoring RNG settings ... ", appendLF=FALSE)							
			setRNG(newRNG)					
			if( verbose > 2 ) message("OK")
			if( verbose > 3 ) showRNG()
		}, add=TRUE)
	}
	#end_RNG
	
	# CHECK PARAMETERS:	
	# test for negative values in x only if the method is not mixed
	if( !is.mixed(method) && min(x, na.rm = TRUE) < 0 )
        fstop('Input matrix ', substitute(x),' contains some negative entries.');
	# test if one row contains only zero entries
    if( min(rowSums(x, na.rm = TRUE), na.rm = TRUE) == 0 )
        fstop('Input matrix ', substitute(x),' contains at least one null or NA-filled row.');	

	# a priori the parameters for the run are all the one in '...'
	# => expand with the strategy's defaults (e.g., maxIter)
	parameters.method <- expand_list(list(...), .method_defaults)
	#
	
	if( is.nmf(seed) ){
		
		if( !is.null(model) )
			fwarning("Discarding argument `model`: directly using NMF model supplied in argument `seed`")
		
		# if the seed is a NMFfit object then only use the fit (i.e. the NMF model)
		# => we want a fresh and clean NMFfit object
		if( isNMFfit(seed) )
			seed <- fit(seed)
		
		# Wrap up the seed into a NMFfit object
		seed <- NMFfit(fit=seed, seed='NMF')
	}
	else if( !inherits(seed, 'NMFfit') ){
		
		## MODEL INSTANTIATION :
	
		# default NMF model is retrieved from the NMF strategy
		.modelClass <- modelname(method)
		# if a character string then use this type of NMF model, but still look 
		# for slots in `...`
		if( is.character(model) ){
			.modelClass <- model
			model <- NULL
		}
		
		# some of the instantiation parameters are set internally
		# TODO: change target into x (=> impact on nmfModel ?
		parameters.model.internal <- list(rank=rank, target=0)
		parameters.model <- list()
		
		init <- 
		if( is.nmf(model) ){
			model
		}else{
			# if 'model' is NULL: initialization parameters are searched in '...' 
			if( is.null(model) ){
				
				# extract the parameters from '...' that correspond to slots in the given class
				stopifnot( isNMFclass(.modelClass) )
				parameters <- .extract.slots.parameters(.modelClass, parameters.method)	
				
				# restrict parameters.method to the ones that won't be used to instantiate the model
				overriden <- is.element(names(parameters$slots), names(parameters.model.internal))
				parameters.method <- c(parameters$extra, parameters$slots[overriden])
							
				#- the model parameters come from the remaining elements
				parameters.model <- c(model=.modelClass, parameters$slots)
				
			} else if( is.list(model) ){  # otherwise argument 'model' must be a list
				
				# if the list is not empty then check all elements are named and 
				# not conflicting with the internally set values			
				if( length(model) > 0 ){
					# all the elements must be named
					if( !hasNames(model, all=TRUE) )  
						fstop("Invalid argument `model` [elements must all be named]. See ?nmf.")
					
					# warn the user if some elements are conflicting and won't be used
					overriden <- is.element(names(model), names(parameters.model.internal))
					if( any(overriden) )
						warning("NMF::nmf - Model parameter(s) [" 
								, str_out(model[overriden], use.names=TRUE, max=Inf)
								, "] discarded. Used internally set value(s) ["
								, str_out(parameters.model.internal[names(model[overriden])], use.names=TRUE, max=Inf)
								, "]"
								, call.=FALSE)
				}
				
				# add default model class if necessary
				if( is.null(model$model) )
					model$model <- .modelClass
				# all the instantiation parameters come from argument 'model'
				parameters.model <- model
				
			}else{ 			
				fstop("Invalid argument 'model' [expected NULL, a character string, or a list to set slots in the NMF model class '",.modelClass,"']. See ?nmf.")
			}	
				
			
			#- force the value of the internally set arguments for the instantiation of the model
			parameters.model <- .merge.override(parameters.model, parameters.model.internal)		
			
			# at this point 'init' should be the list of the initialization parameters
			if( !is.list(parameters.model) ){
				fstop("Unexpected error: object 'parameters.model' must be a list")
			}
			if( !is.element('model', names(parameters.model)) ){
				fstop("Unexpected error: object 'parameters.model' must contain an element named 'model'")
			}
			
			parameters.model
		}
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
				else fstop("Invalid parameter: list 'seed' must contain the seeding method through its first element or through an element named 'method' [", str_desc(seed, 2L), "]")
				
				# check validity of the method provided via the list
				if( !is.function(seed.method) && !(is.character(seed.method) && length(seed.method)==1) )
					fstop("The seeding method provided by parameter 'seed' [", str_desc(seed.method), "] is invalid: a valid function or a character string is expected")
			}
			else fstop("Invalid parameter 'seed'. Acceptable values are:\n\t- ",
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
				fstop("The seeding method function should return class 'NMF' ["
					, if( is.character(seed.method) ) paste('method "', seed.method, "' ", sep='') else NULL 
					, "returned class: '", class(seed), "']")
	}
	# -> at this point the 'seed' object is an instance of class 'NMFfit'
	nmf.debug('nmf', "Seed is of class: '", class(seed), "'")
	# ASSERT just to be sure
	if( !inherits(seed, 'NMFfit') )
		fstop("Invalid class '", class(seed), "' for the computed seed: object that inherits from class 'NMFfit' expected.")
	
	# check the consistency of the NMF model expected by the algorithm and 
	# the one defined by the seed
	#if( none( sapply(model(method), function(c) extends(model(seed), c)) ) )
	if( all( !inherits(fit(seed), modelname(method)) ) )
		fstop("Invalid NMF model '", modelname(seed),"': algorithm '", name(method), "' expects model(s) "
			, paste(paste("'", modelname(method),"'", sep=''), collapse=', ')
			, " or extension.")
	
	# get the complete seeding method's name 
	seed.method <- seeding(seed)
	
	## FINISH SETUP OF THE SEED OBJECT: store some data within the seed so
	# that strategy methods can access them directly
	algorithm(seed) <- name(method) # algorithm name
	seed@distance <- objective(method) # distance name
	seed@parameters <- parameters.method # extra parameters
	run.options(seed) <- nmf.options() # set default run options
	run.options(seed, 'error.track') <- .OPTIONS$track
	if( is.numeric(.OPTIONS$track) )
		run.options(seed, 'track.interval') <- .OPTIONS$track
	run.options(seed, 'verbose') <- verbose
	# store ultimate nmf() call
	seed@call <- match.call()
	##

	## print options if in verbose > 3
	if( verbose > 3 ){
		cat("## OPTIONS:\n")		
		sapply(seq_along(.options)
				, function(i){
					r <- i %% 4
					cat(if(r!=1) '\t| ' else "# ", names(.options)[i],': ', .options[[i]], sep='')
					if(r==0) cat("\n")
				})
		if( length(.options) %% 4 != 0 )cat("\n")
	}
	
	
	## run parameters: 
	parameters.run <- c(list(object=method, y=x, x=seed), parameters.method)
	## Compute the initial residuals if tracking is enabled
	init.resid <- if( .OPTIONS$track && !is.partial.nmf(seed) ){
		do.call('deviance', parameters.run)
	}
	
	## RUN NMF METHOD:
	# call the strategy's run method [and time it]
	t <- system.time({				
		res <- if( !dry.run ){
				do.call('run', parameters.run)
				
			}else{ 
				seed
			}
	})

	## WRAP/CHECK RESULT
	res <- .wrapResult(x, res, seed, method=method, seed.method=seed.method, t)
	if( !isNMFfit(res) ){ # stop if error
		fstop(res)
	}
	##
	
	## CLEAN-UP + EXTRAS:
	# add extra information to the object
	# slot 'parameters'
	if( length(res@parameters) == 0L && length(parameters.method)>0L )
		res@parameters <- parameters.method
	# last residuals
	if( length(residuals(res)) == 0 && !is.partial.nmf(seed) ){
		parameters.run$x <- res
		residuals(res, niter=niter(res)) <- do.call('deviance', parameters.run)
	}
	# first residual if tracking is enabled
	if( .OPTIONS$track && !is.null(init.resid) ){
		if( !hasTrack(res, niter=0) )
			residuals(res, track=TRUE) <- c('0'=init.resid, residuals(res, track=TRUE))
	}
	
	if( length(residuals(res)) && is.na(residuals(res)) ) warning("NMF residuals: final objective value is NA")
	res@runtime <- t
	
	# return the result
	exitSuccess(res)
})

# wrap result
.wrapResult <- function(x, res, seed, method, seed.method, t){
	
	## wrap into an NMFfit object (update seed)
	if( !isNMFfit(res) ){
		# extract expression data if necessary
		if( is(res, 'ExpressionSet') ) res <- exprs(res)
		if( is(x, 'ExpressionSet') ) x <- exprs(x)
		# wrap
		if( is.matrix(res) ){
			if( ncol(res) == ncol(x) ){# partial fit: coef
				# force dimnames
				colnames(res) <- colnames(x)
				res <- nmfModel(H=res)
			}else if( nrow(res) == nrow(x) ){# partial fit: basis
				# force dimnames
				rownames(res) <- rownames(x)
				res <- nmfModel(W=res)
			}
		}else if( is.list(res) ){ # build NMF model from result list
			res <- do.call('nmfModel', res)
		}
		
		# substitute model in fit object
		if( is.nmf(res) ){
			tmp <- seed
			fit(tmp) <- res
			tmp@runtime <- t
			res <- tmp
		}
	}
	
	## check result
	if( !isTRUE(err <- .checkResult(res, seed)) ) return(err) 
	
	## Enforce some slot values
	# slot 'method'
	algorithm(res) <- name(method)	
	# slot 'distance'
	res@distance <- objective(method)	
	# slot 'seed'
	if( seed.method != '' ) seeding(res) <- seed.method
	# set dimnames of the result only if necessary
	if( is.null(dimnames(res)) )
		dimnames(res) <- dimnames(seed)
	
	res
}

# check result
.checkResult <- function(fit, seed){
	# check the result is of the right type
	if( !inherits(fit, 'NMFfit') ){ 
		return(str_c("NMF algorithms should return an instance of class 'NMFfit' [returned class:", class(fit), "]"))
	}
	
	# check that the model has been fully estimated
	if( is.partial.nmf(fit) ){
		warning("nmf - The NMF model was only partially estimated [dim = (", str_out(dim(fit), Inf),")].")
	}
	# check that the fit conserved all fixed terms (only warning)
	if( nterms(seed) ){
		if( length(i <- icterms(seed)) && !identical(coef(fit)[i,], coef(seed)[i,]) ){
			warning("nmf - Fixed coefficient terms were not all conserved in the fit: the method might not support them.")
		}
		if( length(i <- ibterms(seed)) && !identical(basis(fit)[,i], basis(seed)[,i]) ){
			warning("nmf - Fixed basis terms were not all conserved in the fit: the method might not support them.")
		}
	}
	TRUE
}

#' Interface for NMF Seeding Methods
#' 
#' @description
#' The function \code{seed} provides a single interface for calling all seeding 
#' methods used to initialise NMF computations.
#' These methods at least set the basis and coefficient matrices of the initial 
#' \code{object} to valid nonnegative matrices.
#' They will be used as a starting point by any NMF algorithm that accept 
#' initialisation.
#' 
#' IMPORTANT: this interface is still considered experimental and is subject 
#' to changes in future release. 
#' 
#' @param x target matrix one wants to approximate with NMF
#' @param model specification of the NMF model, e.g., the factorization rank.
#' @param method specification of a seeding method.
#' See each method for details on the supported formats. 
#' @param ... extra to allow extensions and passed down to the actual seeding method.
#'  
#' @return an \code{\linkS4class{NMFfit}} object.
#' 
#' @inline
#' @export 
setGeneric('seed', function(x, model, method, ...) standardGeneric('seed') )

#' This is the workhorse method that seeds an NMF model object using a given 
#' seeding strategy defined by an \code{NMFSeed} object, to fit a given 
#' target matrix.
#' 
#' @param rng rng setting to use. 
#' If not missing the RNG settings are set and restored on exit using 
#' \code{\link{setRNG}}. 
#' 
#' All arguments in \code{...} are passed to teh seeding strategy.
#' 
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
		res <- NMFfit()
		# ASSERT: check that the RNG seed is correctly set
		stopifnot( rng.equal(res,rng.s) )
		# call the seeding function passing the extra parameters
		f <- do.call(algorithm(method), c(list(model, x), ...))
		# set the dimnames from the target matrix
		dimnames(f) <- dimnames(x)
		# set the basis names from the model if any
		if( !is.null(basisnames(model)) )
			basisnames(f) <- basisnames(model)
		# store the result into the NMFfit object
		fit(res) <- f
		
		# if not already set: store the seeding method's name in the resulting object
		if( seeding(res) == '' ) seeding(res) <- name(method)
		
		# return the seeded object
		res
	}
)
#' Seeds an NMF model using a custom seeding strategy, defined by a function.
#' 
#' \code{method} must have signature \code{(x='NMFfit', y='matrix', ...)}, where 
#' \code{x} is the unseeded NMF model and \code{y} is the target matrix to fit.
#' It must return an \code{\linkS4class{NMF}} object, that contains the seeded 
#' NMF model.
#' 
#' @param name optional name of the seeding method for custom seeding strategies.   
#'  
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
#' Seeds the model with the default seeding method given by 
#' \code{nmf.getOption('default.seed')} 
setMethod('seed', signature(x='ANY', model='ANY', method='missing'),
	function(x, model, method, ...){
		seed(x, model, nmf.getOption('default.seed'), ...)
	}
)
#' Use NMF method \code{'none'}.
setMethod('seed', signature(x='ANY', model='ANY', method='NULL'),
	function(x, model, method, ...){
		seed(x, model, 'none', ...)
	}
)
#' Use \code{method} to set the RNG with \code{\link{setRNG}} and use method 
#' \dQuote{random} to seed the NMF model.
#' 
#' Note that in this case the RNG settings are not restored.
#' This is due to some internal technical reasons, and might change in future 
#' releases. 
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
#' Use the registered seeding method whose access key is \code{method}. 
setMethod('seed', signature(x='ANY', model='ANY', method='character'),
		function(x, model, method, ...){
			
			# get the seeding method from the registry
			seeding.fun <- nmfSeed(method)
				
			#Vc#Use seeding method: '${method}'
			# call 'seed' with the seeding.function			
			seed(x, model, method=seeding.fun, ...)
			
		}
)
#' Seed a model using the elements in \code{model} to instantiate it with 
#' \code{\link{nmfModel}}. 
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
				val <- unlist(model[idx], recursive=FALSE)				
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
				val <- unlist(model[idx], recursive=FALSE)
				if( is.numeric(val) && length(val)==1 )
					names(model)[idx] <- 'rank'
				else stop(err.msg)
			}
			else stop(err.msg)
		}
					
		nmf.debug('seed', "using model parameters:\n", capture.output(print(model)) )
		# instantiate the object using the factory method		
		model <- do.call('nmfModel', model)
		nmf.debug('seed', "using NMF model '", class(model), "'")
		
		# check that model is from the right type, i.e. inherits from class NMF
		if( !inherits(model, 'NMF') ) stop("Invalid object returned by model: object must inherit from class 'NMF'")
		
		seed(x, model, method, ...)
	}
)
#' Seeds a standard NMF model (i.e. of class \code{\linkS4class{NMFstd}}) of rank 
#' \code{model}. 
setMethod('seed', signature(x='ANY', model='numeric', method='NMFSeed'), 
	function(x, model, method, ...){	

		seed(x, nmfModel(model), method, ...)
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
	
	if( length(parameters) == 1L && is.null(names(parameters)) ){
		parameters <- parameters[[1L]]
	}

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

#' Estimate Rank for NMF Models
#' 
#' A critical parameter in NMF algorithms is the factorization rank \eqn{r}.
#' It defines the number of basis effects used to approximate the target
#' matrix. 
#' Function \code{nmfEstimateRank} helps in choosing an optimal rank by
#' implementing simple approaches proposed in the literature.
#' 
#' Note that from version \emph{0.7}, one can equivalently call the 
#' function \code{\link{nmf}} with a range of ranks.
#' 
#' @details 
#' Given a NMF algorithm and the target matrix, a common way of estimating
#' \eqn{r} is to try different values, compute some quality measures of the
#' results, and choose the best value according to this quality criteria. See
#' \cite{Brunet2004} and \cite{Hutchins2008}.
#' 
#' The function \code{nmfEstimateRank} allows to perform this estimation
#' procedure. 
#' It performs multiple NMF runs for a range of rank of
#' factorization and, for each, returns a set of quality measures together with
#' the associated consensus matrix.
#' 
#' In order to avoid overfitting, it is recommended to run the same procedure on 
#' randomized data.
#' The results on the original and the randomised data may be plotted on the 
#' same plots, using argument \code{y}.
#' 
#' @param x For \code{nmfEstimateRank} a target object to be estimated, in one
#' of the format accepted by interface \code{\link{nmf}}.
#' 
#' For \code{plot.NMF.rank} an object of class \code{NMF.rank} as returned by
#' function \code{nmfEstimateRank}.
#' @param range a \code{numeric} vector containing the ranks of factorization
#' to try.
#' Note that duplicates are removed and values are sorted in increasing order. 
#' The results are notably returned in this order. 
#' 
#' @param method A single NMF algorithm, in one of the format accepted by
#' the function \code{\link{nmf}}.
#' 
#' @param nrun a \code{numeric} giving the number of run to perform for each
#' value in \code{range}.
#' 
#' @param model model specification passed to each \code{nmf} call.
#' In particular, when \code{x} is a formula, it is passed to argument 
#' \code{data} of \code{\link{nmfModel}} to determine the target matrix -- and 
#' fixed terms.
#' 
#' @param verbose toggle verbosity.  This parameter only affects the verbosity
#' of the outer loop over the values in \code{range}. 
#' To print verbose (resp. debug) messages from each NMF run, one can use 
#' \code{.options='v'} (resp. \code{.options='d'}) 
#' that will be passed to the function \code{\link{nmf}}.
#' 
#' @param stop logical flag for running the estimation process with fault
#' tolerance.  When \code{TRUE}, the whole execution will stop if any error is
#' raised.  When \code{FALSE} (default), the runs that raise an error will be
#' skipped, and the execution will carry on. The summary measures for the runs
#' with errors are set to NA values, and a warning is thrown.
#' 
#' @param ... For \code{nmfEstimateRank}, these are extra parameters passed
#' to interface \code{nmf}. Note that the same parameters are used for each
#' value of the rank.  See \code{\link{nmf}}.
#' 
#' For \code{plot.NMF.rank}, these are extra graphical parameter passed to the
#' standard function \code{plot}. See \code{\link{plot}}.
#' 
#' @return 
#' \code{nmfEstimateRank} returns a S3 object (i.e. a list) of class 
#' \code{NMF.rank} with the following elements:
#'  
#' \item{measures }{a \code{data.frame} containing the quality
#' measures for each rank of factorizations in \code{range}. Each row
#' corresponds to a measure, each column to a rank. } 
#' \item{consensus }{ a
#' \code{list} of consensus matrices, indexed by the rank of factorization (as
#' a character string).}
#' \item{fit }{ a \code{list} of the fits, indexed by the rank of factorization
#' (as a character string).}
#' 
#' @export
#' @examples
#'
#' if( !isCHECK() ){
#'  
#' set.seed(123456)
#' n <- 50; r <- 3; m <- 20
#' V <- syntheticNMF(n, r, m)
#' 
#' # Use a seed that will be set before each first run
#' res <- nmfEstimateRank(V, seq(2,5), method='brunet', nrun=10, seed=123456)
#' # or equivalently
#' res <- nmf(V, seq(2,5), method='brunet', nrun=10, seed=123456)
#' 
#' # plot all the measures
#' plot(res)
#' # or only one: e.g. the cophenetic correlation coefficient
#' plot(res, 'cophenetic')
#' 
#' # run same estimation on randomized data
#' rV <- randomize(V)
#' rand <- nmfEstimateRank(rV, seq(2,5), method='brunet', nrun=10, seed=123456)
#' plot(res, rand)
#' }
#' 
nmfEstimateRank <- function(x, range, method=nmf.getOption('default.algorithm')
					, nrun=30, model=NULL, ..., verbose=FALSE, stop=FALSE){
	
	# fix method if passed NULL (e.g., from nmf('formula', 'numeric'))
	if( is.null(method) )
		method <- nmf.getOption('default.algorithm')
	
	# special handling of formula: get target data from the formula 
	if( is(x, 'formula') ){
		# dummy model to resolve formula
		dummy <- nmfModel(x, 0L, data=model)
		# retrieve target data
		V <- attr(dummy, 'target')
	}else{
		V <- x
	}
	
	# remove duplicates and sort
	range <- sort(unique(range))
	
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
			tmp.res <- matrix(as.numeric(NA), n, length(range))
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
	k.rank <- 0 
	measures <- sapply(range, function(r, ...){
			k.rank <<- k.rank + 1L
			if( verbose ) cat("Compute NMF rank=", r, " ... ")
			
			# restore RNG on exit (except after last rank)
			# => this ensures the methods use the same stochastic environment
			orng <- RNGseed()
			if( k.rank < length(range) ) on.exit( RNGseed(orng), add = TRUE)
			
			res <- tryCatch({ #START_TRY
				
				res <- nmf(x, r, method, nrun=nrun, model=model, ...)					
				# directly return the result if a valid NMF result
				if( !isNMFfit(res, recursive = FALSE) )
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
				measures <- summary(res, target=V)
				
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

#' @S3method summary NMF.rank
summary.NMF.rank <- function(object, ...){
	s <- summary(new('NMFList', object$fit), ...)
	# NB: sort measures in the same order as required in ...
	i <- which(!names(s) %in% names(object$measures))
	cbind(s[, i], object$measures[match(object$measures$rank, s$rank), ])
}


#' \code{plot.NMF.rank} plots the result of rank estimation survey.
#' 
#' In the plot generated by \code{plot.NMF.rank}, each curve represents a 
#' summary measure over the range of ranks in the survey.
#' The colours correspond to the type of data to which the measure is related:
#' coefficient matrix, basis component matrix, best fit, or consensus matrix. 
#' 
#' @param y reference object of class \code{NMF.rank}, as returned by
#' function \code{nmfEstimateRank}. 
#' The measures contained in \code{y} are used and plotted as a reference.
#' It is typically used to plot results obtained from randomized data.  
#' The associated curves are drawn in \emph{red} (and \emph{pink}), 
#' while those from \code{x} are drawn in \emph{blue} (and \emph{green}).
#' @param what a \code{character} vector whose elements partially match 
#' one of the following item, which correspond to the measures computed
#' by \code{\link{summary}} on each -- multi-run -- NMF result: 
#' \sQuote{all}, \sQuote{cophenetic}, \sQuote{rss},
#' \sQuote{residuals}, \sQuote{dispersion}, \sQuote{evar}, 
#' \sQuote{silhouette} (and more specific *.coef, *.basis, *.consensus), 
#' \sQuote{sparseness} (and more specific *.coef, *.basis).
#' It specifies which measure must be plotted (\code{what='all'} plots 
#' all the measures).
#' @param na.rm single logical that specifies if the rank for which the
#' measures are NA values should be removed from the graph or not (default to
#' \code{FALSE}).  This is useful when plotting results which include NAs due
#' to error during the estimation process. See argument \code{stop} for
#' \code{nmfEstimateRank}.
#' @param xname,yname legend labels for the curves corresponding to measures from 
#' \code{x} and \code{y} respectively 
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param main main title
#' 
#' @S3method plot NMF.rank
#' @rdname nmfEstimateRank
#' @import ggplot2
#' @import reshape2
plot.NMF.rank <- function(x, y=NULL, what=c('all', 'cophenetic', 'rss', 'residuals'
									, 'dispersion', 'evar', 'sparseness'
									, 'sparseness.basis', 'sparseness.coef'
                                    , 'silhouette'
                                    , 'silhouette.coef', 'silhouette.basis'
                                    , 'silhouette.consensus')
						, na.rm=FALSE
                        , xname = 'x'
                        , yname = 'y'
                        , xlab = 'Factorization rank'
                        , ylab = ''
                        , main = 'NMF rank survey'
                        , ... ){

	
    # trick for convenience 
	if( is.character(y) && missing(what) ){
		what <- y
		y <- NULL
	}
	
	what <- match.arg(what, several.ok=TRUE)
    if( 'all' %in% what ){
        what <- c('cophenetic', 'rss', 'residuals', 'dispersion', 'evar', 'sparseness', 'silhouette')
    }
    
    .getvals <- function(x, xname){
    	measures <- x$measures
    	iwhat <- unlist(lapply(paste('^',what,sep=''), grep, colnames(measures)))
    	
    	# remove NA values if required
    	if( na.rm )
    		measures <- measures[ apply(measures, 1, function(row) !any(is.na(row[iwhat]))), ]
    	
    	vals <- measures[,iwhat, drop=FALSE]
    	x <- as.numeric(measures$rank)
    	xlim <- range(x)
        
        # define measure type
        measure.type <- setNames(rep('Best fit', ncol(measures)), colnames(measures))
        cons.measures <- c('silhouette.consensus', 'cophenetic', 'cpu.all')
        measure.type[match(cons.measures, names(measure.type))] <- 'Consensus'
        measure.type[grep("\\.coef$", names(measure.type))] <- 'Coefficients'
        measure.type[grep("\\.basis$", names(measure.type))] <- 'Basis'
        measure.type <- factor(measure.type)
        
        pdata <- melt(cbind(rank = x, vals), id.vars = 'rank')
        # set measure type
        pdata$Type <- measure.type[as.character(pdata$variable)]
        # define measure groups
        pdata$Measure <- gsub("^([^.]+).*", "\\1", pdata$variable)
        pdata$Data <- xname
        pdata
    }
    
    pdata <- .getvals(x, xname)
    
    # add reference data
    if( is(y, 'NMF.rank') ){
        pdata.y <- .getvals(y, yname)
        pdata <- rbind(pdata, pdata.y)
    }
    
    p <- ggplot(pdata, aes_string(x = 'rank', y = 'value')) +
            geom_line( aes_string(linetype = 'Data', colour = 'Type') ) +
            geom_point(size = 2, aes_string(shape = 'Data', colour = 'Type') ) +
            theme_bw() +
            scale_x_continuous(xlab, breaks = unique(pdata$rank)) +
            scale_y_continuous(ylab) +
            ggtitle(main)
    # remove legend if not necessary
    if( !is(y, 'NMF.rank') ){
        p <- p + scale_shape(guide = 'none') + scale_linetype(guide = 'none')
    }
    
    # use fix set of colors
    myColors <- brewer.pal(5,"Set1")
    names(myColors) <- levels(pdata$Type)
    p <- p + scale_colour_manual(name = "Measure type", values = myColors)
    
    # add facet
    p <- p + facet_wrap( ~ Measure, scales = 'free')
    
    # return plot
    p
}

