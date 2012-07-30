#' @include NMFStrategy-class.R
#' @include NMFfit-class.R
NULL

#' Interface for Algorithms: Implementation for Iterative NMF Algorithms
#' 
#' @description
#' This class provides a specific implementation for the function \code{run} -- completing the 
#' the interface class \code{\linkS4class{NMFStrategy}}, 
#' for NMF algorithms that conform to the following iterative schema:
#'  
#' \itemize{
#' \item 1. Initialisation
#' \item 2. Update the model at each iteration
#' \item 3. Stop if some criterion is satisfied
#' \item 4. Wrap up
#' }
#' 
#' This schema could possibly apply to all NMF algorithms, since these are essentially optimisation algorithms, 
#' almost all of which use iterative methods to approximate a solution of the optimisation problem.
#' The main advantage is that it allows to implement updates and stopping criterion separately, and combine them
#' in different ways.
#' In particular, many NMF algorithms are based on multiplicative updates, following the approach from  
#' \cite{Lee2001}, which are specially suitable to be cast into this simple schema. 
#'  
#' @slot onInit function that performs some initialisation or pre-processing on 
#' the model, before starting the iteration loop.
#' @slot Update function that implement the update step, which computes new values for the model, based on its
#' previous value.
#' It is called at each iteration, until the stopping criterion is met or the maximum number of iteration is 
#' achieved.
#' @slot Stop function that implements the stopping criterion.
#' It is called \emph{before} each Update step.
#' @slot onReturn function that wraps up the result into an NMF object.
#' It is called just before returning the 
#'  
setClass('NMFStrategyIterative'
	, representation(
                onInit = '.functionSlot.null',
				Update = '.functionSlot', # update method	
				Stop = '.functionSlot.null', # method called just after the update
				onReturn = '.functionSlot.null' # method called just before returning the resulting NMF object
				)	
  , prototype=prototype(
          		onInit = NULL
				, Update = ''
				, Stop = NULL
				, onReturn = NULL
			)
	, contains = 'NMFStrategy'
	, validity = function(object){
		
		if( is.character(object@Update) && object@Update == '' )
			return("Slot 'Update' is required")
		
		# check the arguments of methods 'Update' and 'Stop'
		# (except for the 3 mandatory ones)
		n.update <- names(formals(object@Update))
		
		# at least 3 arguments for 'Update'
		if( length(n.update) < 3 )
			return("Invalid 'Update' method - must have at least 3 arguments: 'i' [iteration number], 'target' [target matrix], 'data' [current NMF model]")
		
		n.update <- n.update[-seq(3)]
		# argument '...' must be present in method 'Update'
		if( !is.element('...', n.update) )
			return("Invalid 'Update' method: must have argument '...' (even if not used)")

		# at least 3 arguments for 'Stop'
		if( !is.null(object@Stop) ){
			
			# retrieve the stopping criterion 
			.stop <- NMFStop(object@Stop)
			n.stop <- names(formals(.stop))
			
			if( length(n.stop) < 4 )
				return("Invalid 'Stop' method - must have at least 4 arguments: 'strategy', 'i', 'target' [target matrix], 'data' [current NMF model]")
			
			n.stop <- n.stop[-seq(4)]
			# argument '...' must be present in method 'Stop'
			if( !is.element('...', n.stop) )
				return("Invalid 'Stop' method: must have argument '...' (even if not used)")
		
			# Update and Stop methods cannot have overlapping arguments 
			overlap <- intersect(n.update, n.stop)
			overlap <- overlap[which(overlap!='...')]
			if( length(overlap) > 0 )
				return("Invalid 'Update' and 'Stop' methods - conflicting arguments")
		}
		
		TRUE
	}
)


#' Show method for objects of class \code{NMFStrategyIterative}
#' @export
setMethod('show', 'NMFStrategyIterative',
	function(object){
		
		#cat('<object of class: NMFStrategyIterative>')
		callNextMethod()
		cat(" <Iterative schema>\n")
		# go through the slots
		s.list <- names(getSlots('NMFStrategyIterative'))
		s.list <- setdiff(s.list, names(getSlots('NMFStrategy')))
		#s.list <- s.list[s.list=='ANY']
#		s.list <- c('Update', 'Stop', 'WrapNMF')
		out <-
		sapply(s.list, function(sname){
					svalue <- slot(object,sname)
					svalue <- 
					if( is.function(svalue) ) {
						str_args(svalue, exdent=12)
					} else if( is.null(svalue) ){
						'none'
					} else { 
						paste("'", svalue,"'", sep='')
					}
					str_c(sname, ": ", svalue)
				})
		cat(str_c('  ', out, collapse='\n'), "\n", sep='')
		return(invisible())
	}
)
###% This class is an auxiliary class that defines the strategy's methods by directly callable functions. 
setClass('NMFStrategyIterativeX'
	, contains = 'NMFStrategyIterative'
	, representation = representation(
				workspace = 'environment' # workspace to use persistent variables accross methods
				)
)


###% Creates a NMFStrategyIterativeX object from a NMFStrategyIterative object.
xifyStrategy <- function(strategy, workspace){	
	
	# first check the strategy's validity
	if( is.character(err <- validObject(strategy, test=TRUE)) ){
		stop("Invalid strategy definition:\n\t- ", err)
	}
	
	# intanciate the NMFStrategyIterativeX, creating the strategy's workspace
	strategyX <- new('NMFStrategyIterativeX', strategy, workspace=workspace)
	
	# define auxiliary function to preload the 'function' slots in class NMFStrategyIterativeX
	preload.slot <- function(strategy, sname, default){
		
		# get the content of the slot
		svalue <- slot(strategy,sname)
		
		# if the slot is valid (i.e. it's a non-empty character string), then process the name into a valid function
		fun <-
		if( is.null(svalue) && !missing(default) ) default
		else if( sname == 'Stop' ) NMFStop(svalue)
		else if( is.character(svalue) && nchar(svalue) > 0 ){
			# set the slot with the executable version of the function			
			getFunction(svalue)
		}else if( is.function(svalue) )	svalue		
		else
			stop("NMFStrategyIterativeX - could not pre-load slot '", sname, "'")		

		# return the loaded function
		fun
	}
	
	# preload the function slots
	slot(strategyX, 'Update') <- preload.slot(strategyX, 'Update')
	slot(strategyX, 'Stop') <- preload.slot(strategyX, 'Stop', function(strategy, i, target, data, ...){FALSE})
	slot(strategyX, 'onReturn') <- preload.slot(strategyX, 'onReturn', identity)
	
	# load the objective function
	objective(strategyX) <- nmfDistance(objective(strategy))

	# valid the preloaded object
	validObject(strategyX)
	
	# return the executable strategy 
	strategyX
}

#
#setGeneric('Update', function(object, v, ...) standardGeneric('Update') )
#setMethod('Update', signature(object='NMFStrategyIterative', v='matrix'), function(object, v, ...){ object@data <- object@Update(v, object@data, ...) })
#
#setGeneric('Stop', function(object, i) standardGeneric('Stop') )
#setMethod('Stop', signature(object='NMFStrategyIterative', i='integer'), function(object, i){ object@Stop(i, object@data) })
#
#setGeneric('WrapNMF', function(object) standardGeneric('WrapNMF') )
#setMethod('WrapNMF', signature(object='NMFStrategyIterative'), function(object){ object@WrapNMF(object@data) })

###% Hook to initialize built-in iterative methods when the package is loaded


###% Hook to initialize old R version built-in iterative methods

#' Get/Set a Static Variable in NMF Algorithms
#' 
#' @description
#' This function is used in iterative NMF algorithms to manage variables
#' stored in a local workspace, that are accessible to all functions that 
#' define the iterative schema described in \code{\linkS4class{NMFStrategyIterative}}.
#' 
#' It is specially useful for computing stopping criteria, which often require model data from   
#' different iterations.
#' 
#' @param name Name of the static variable (as a single character string)
#' @param value New value of the static variable
#' @param init a logical used when a \code{value} is provided, that specifies 
#' if the variable should be set to the new value only if it does not exist yet 
#' (\code{init=TRUE}). 
#' @return The value of the static variable
#' @export
staticVar <- local({
	
	.Workspace <- NULL
	function(name, value, init=FALSE){	
		
		# return last workspace
		if( missing(name) ) return(.Workspace)			
		else if( is.environment(name) ){ # setup up static environment			
			nmf.debug('Strategy Workspace', "initialize static workspace: ", capture.output(.Workspace), "=", capture.output(name))
			.Workspace <<- name
		}else if( missing(value) ){
			get(name, envir=.Workspace, inherits=FALSE)
		}else{
			if( !init || !exists(name, envir=.Workspace, inherits=FALSE) )
			{
				if( init ) nmf.debug('Strategy Workspace', "initialize variable '", name, "'")
				assign(name, value, envir=.Workspace)
			}
		}
		
	}
})

#' Runs an NMF iterative algorithm on a target matrix \code{y}.
#' 
#' @param .stop specification of a stopping criterion, that is used instead of the 
#' one associated to the NMF algorithm.
#' It may be specified as:
#' \itemize{
#' \item the access key of a registered stopping criterion;
#' \item a single integer that specifies the exact number of iterations to perform, which will 
#' be honoured unless a lower value is explicitly passed in argument \code{maxIter}.
#' \item a single numeric value that specifies the stationnarity threshold for the 
#' objective function, used in with \code{\link{nmf.stop.stationary}}; 
#' \item a function with signature \code{(object="NMFStrategy", i="integer", y="matrix", x="NMF", ...)}, 
#' where \code{object} is the \code{NMFStrategy} object that describes the algorithm being run, 
#' \code{i} is the current iteration, \code{y} is the target matrix and \code{x} is the current value of 
#' the NMF model.  
#' }
#' @param maxIter maximum number of iterations to perform.
#'   
#' @rdname NMFStrategy
setMethod('run', signature(object='NMFStrategyIterative', y='matrix', x='NMFfit'),
	function(object, y, x, .stop=NULL, maxIter=2000L, ...){
	
	method <- object
	# override the stop method on runtime
	if( !is.null(.stop) ){
		method@Stop <- NMFStop(.stop)
		if( is.integer(.stop) && missing(maxIter) )
			maxIter <- .stop[1]
	}
	
	# debug object in debug mode
	if( nmf.getOption('debug') ) show(method)		
	
	#Vc# Define local workspace for static variables
	# this function can be called in the methods to get/set/initialize 
	# variables that are persistent within the strategy's workspace
	.Workspace <- new.env()	
	staticVar(.Workspace)
		
	# runtime resolution of the strategy's functions by their names if necessary
	strategyX = xifyStrategy(method, .Workspace)
	run(strategyX, y, x, maxIter=maxIter, ...)
})

#' @rdname NMFStrategy
setMethod('run', signature(object='NMFStrategyIterativeX', y='matrix', x='NMFfit'),
	function(object, y, x, maxIter, ...){
				
	strategy <- object
	v <- y
	seed <- x
	#V!# NMFStrategyIterativeX::run
	
	#Vc# Define workspace accessor function
	# this function can be called in the methods to get/set/initialize 
	# variables that are persistent within the strategy's workspace
#	.Workspace <- strategy@workspace
#	assign('staticVar', function(name, value, init=FALSE){
#			if( missing(value) ){
#				get(name, envir=.Workspace, inherits=FALSE)
#			}else{
#				if( !init || !exists(name, envir=.Workspace, inherits=FALSE) )
#				{
#					if( init ) nmf.debug('Strategy Workspace', "initialize variable '", name, "'")
#					assign(name, value, envir=.Workspace)
#				}
#			}
#		}
#		, envir=.Workspace)
	
	#Vc# initialize the strategy
	# check validity of arguments if possible
	update.args <- formals(strategy@Update)
	stop.args <- formals(strategy@Stop)
	internal.args <- names(c(update.args[1:3], stop.args[1:4]))
	expected.args <- c(update.args[-(1:3)], stop.args[-(1:4)])
	passed.args <- names(list(...))
	forbidden.args <- is.element(passed.args, c(internal.args))
	if( any(forbidden.args) )
		stop("NMF::run - Update/Stop method : formal argument(s) "
			, paste( paste("'", passed.args[forbidden.args],"'", sep=''), collapse=', ')
			, " already set internally.", call.=FALSE)
	# !is.element('...', expected.args) && 
	if( any(t <- !pmatch(passed.args, names(expected.args), nomatch=FALSE)) )
		stop("NMF::run - Update/Stop method for algorithm '", name(strategy),"': unused argument(s) "
			, paste( paste("'", passed.args[t],"'", sep=''), collapse=', '), call.=FALSE)
	# check for required arguments
	required.args <- sapply(expected.args, function(x){ x <- as.character(x); length(x) == 1 && nchar(x) == 0 } )
	required.args <- names(expected.args[required.args])
	required.args <- required.args[required.args!='...']
	
	if( any(t <- !pmatch(required.args, passed.args, nomatch=FALSE)) )
		stop("NMF::run - Update/Stop method for algorithm '", name(strategy),"': missing required argument(s) "
			, paste( paste("'", required.args[t],"'", sep=''), collapse=', '), call.=FALSE)
	
	# set default value for missing argument TODO
	missing.args <- expected.args[!pmatch(names(expected.args), passed.args, nomatch=FALSE)]
	
	
	#Vc# Start iterations
	nmfData <- seed
	# cache verbose level
	verbose <- verbose(nmfData)
	
	# clone the object to allow the updates to work in place
	if( verbose > 1 ) 
		message("Duplicating the seed NMF object ", appendLF=FALSE)
	nmfFit <- clone(fit(nmfData))
	if( verbose > 1 )
		message("[", C.ptr(fit(nmfData)), " -> ", C.ptr(nmfFit), "]")		
	
	# pre-load slots
	updateFun <- strategy@Update
	stopFun <- strategy@Stop
	
	if( verbose ) cat('Iterations:')
	i <- 0L
	while( TRUE ){
		
		#Vc# Stopping criteria
		# check convergence (generally do not stop for i=0L, but only initialise static variables
		stop.signal <- stopFun(strategy, i, v, nmfFit, ...)
		
		# if the strategy ask for stopping, then stop the iteration
		if( stop.signal || i >= maxIter ) break;
		
		# increment i
		i <- i+1L
		
		if( verbose && (i==1L || i %% 50 == 0) ) cat('', i)
		
		#Vc# update the matrices
		nmfFit <- updateFun(i, v, nmfFit, ...)
		
		# every now and then track the error if required
		nmfData <- trackError(nmfData, deviance(strategy, nmfFit, v, ...), i)
				
	}
	if( verbose ) cat("\nDONE (stopped at ",i,'/', maxIter," iterations)\n", sep='')
	
	# force to compute last error if not already done
	nmfData <- trackError(nmfData, deviance(strategy, nmfFit, v, ...), i, force=TRUE)
	
	# store the fitted model
	fit(nmfData) <- nmfFit
	
	#Vc# wrap up
	# let the strategy build the result
	nmfData <- strategy@onReturn(nmfData)
	if( !inherits(nmfData, 'NMFfit') ){
		stop('NMFStrategyIterative[', name(strategy), ']::onReturn did not return a "NMF" instance [returned: "', class(nmfData), '"]')
	}
	
	# set the number of iterations performed
	niter(nmfData) <- i
	
	#return the result
	nmf.debug('NMFStrategyIterativeX::run', 'Done')
	invisible(nmfData)
})

################################################################################################
# INITIALIZATION METHODS
################################################################################################

################################################################################################
# UPDATE METHODS
################################################################################################

#' NMF Multiplicative Updates for Kullback-Leibler Divergence
#' 
#' Multiplicative updates from \cite{Lee2001} for standard Nonnegative Matrix Factorization 
#' models \eqn{V \approx W H}, where the distance between the target matrix and its NMF 
#' estimate is measured by the Kullback-Leibler divergence.
#' 
#' \code{nmf_update.KL.w} and \code{nmf_update.KL.h} compute the updated basis and coefficient 
#' matrices respectively.
#' They use a \emph{C++} implementation which is optimised for speed and memory usage. 
#' 
#' @details
#' The coefficient matrix (\code{H}) is updated as follows:
#' \deqn{
#' H_{kj} \leftarrow H_{kj}  \frac{\left( sum_i \frac{W_{ik} V_{ij}}{(WH)_{ij}} \right)}{ sum_i W_{ik} }.
#' }{
#' H_kj <- H_kj ( sum_i [ W_ik V_ij / (WH)_ij ] ) / ( sum_i W_ik )
#' }
#' 
#' @param v target matrix
#' @param w current basis matrix
#' @param h current coefficient matrix
#' @param nbterms number of fixed basis terms
#' @param ncterms number of fixed coefficient terms
#' @param copy logical that indicates if the update should be made on the original
#' matrix directly (\code{FALSE}) or on a copy (\code{TRUE} - default).
#' With \code{copy=FALSE} the memory footprint is very small, and some speed-up may be 
#' achieved in the case of big matrices.
#'
#' @return a matrix of the same dimension as the input matrix to update 
#' (i.e. \code{w} or \code{h}).
#' If \code{copy=FALSE}, the returned matrix uses the same memory as the input object.
#' 
#' @author 
#' Update definitions by \cite{Lee2001}. 
#' 
#' C++ optimised implementation by Renaud Gaujoux. 
#' 
#' @rdname nmf_update_KL
#' @aliases nmf_update.KL
#' @export
nmf_update.KL.h <- std.divergence.update.h <- function(v, w, h, nbterms=0L, ncterms=0L, copy=TRUE)
{	
	.Call("divergence_update_H", v, w, h, nbterms, ncterms, copy, PACKAGE='NMF')
}
#' \code{nmf_update.KL.w_R} and \code{nmf_update.KL.h_R} implement the same updates 
#' in \emph{plain R}.
#' 
#' @param wh already computed NMF estimate used to compute the denominator term.
#' 
#' @rdname nmf_update_KL
#' @export   
nmf_update.KL.h_R <- R_std.divergence.update.h <- function(v, w, h, wh=NULL)
{	
	# compute WH if necessary	
	if( is.null(wh) ) wh <- w %*% h
	
	# divergence-reducing NMF iterations
	# H_au = H_au ( sum_i [ W_ia V_iu / (WH)_iu ] ) / ( sum_k W_ka ) -> each row of H is divided by a the corresponding colSum of W
	h * crossprod(w, v / wh) / colSums(w)	
}

#' @details
#' The basis matrix (\code{W}) is updated as follows:
#' \deqn{
#' W_{ik} \leftarrow W_{ik} \frac{ sum_j [\frac{H_{kj} A_{ij}}{(WH)_{ij}} ] }{sum_j H_{kj} }
#' }{
#' W_ik <- W_ik ( sum_u [H_kl A_il / (WH)_il ] ) / ( sum_l H_kl )
#' }
#' @rdname nmf_update_KL
#' @export
nmf_update.KL.w <- std.divergence.update.w <- function(v, w, h, nbterms=0L, ncterms=0L, copy=TRUE)
{	
	.Call("divergence_update_W", v, w, h, nbterms, ncterms, copy, PACKAGE='NMF')
}
#' @rdname nmf_update_KL
#' @export
nmf_update.KL.w_R <- R_std.divergence.update.w <- function(v, w, h, wh=NULL)
{			
	# compute WH if necessary	
	if( is.null(wh) ) wh <- w %*% h
	
	# W_ia = W_ia ( sum_u [H_au A_iu / (WH)_iu ] ) / ( sum_v H_av ) -> each column of W is divided by a the corresponding rowSum of H
	#x2 <- matrix(rep(rowSums(h), nrow(w)), ncol=ncol(w), byrow=TRUE); 
	#w * tcrossprod(v / wh, h) / x2;
	sweep(w * tcrossprod(v / wh, h), 2L, rowSums(h), "/", check.margin = FALSE) # optimize version?
	
}



#' NMF Multiplicative Updates for Euclidean Distance
#' 
#' Multiplicative updates from \cite{Lee2001} for standard Nonnegative Matrix Factorization 
#' models \eqn{V \approx W H}, where the distance between the target matrix and its NMF 
#' estimate is measured by the -- euclidean -- Frobenius norm.
#' 
#' \code{nmf_update.euclidean.w} and \code{nmf_update.euclidean.h} compute the updated basis and coefficient 
#' matrices respectively.
#' They use a \emph{C++} implementation which is optimised for speed and memory usage. 
#' 
#' @details
#' The coefficient matrix (\code{H}) is updated as follows:
#' \deqn{
#' H_{kj} \leftarrow \frac{\max(H_{kj} W^T V)_{kj}, \varepsilon) }{(W^T W H)_{kj} + \varepsilon}
#' }{
#' H_kj <- max(H_kj (W^T V)_kj, eps) / ( (W^T W H)_kj + eps )
#' }
#' 
#' @param v target matrix
#' @param w current basis matrix
#' @param h current coefficient matrix
#' @param eps small numeric value used to ensure numeric stability
#' @param nbterms number of fixed basis terms (see \code{\link{bterms}}).
#' @param ncterms number of fixed coefficient terms (see \code{\link{cterms}})
#' @param copy logical that indicates if the update should be made on the original
#' matrix directly (\code{FALSE}) or on a copy (\code{TRUE} - default).
#' With \code{copy=FALSE} the memory footprint is very small, and some speed-up may be 
#' achieved in the case of big matrices.
#'
#' @return a matrix of the same dimension as the input matrix to update 
#' (i.e. \code{w} or \code{h}).
#' If \code{copy=FALSE}, the returned matrix uses the same memory as the input object.
#' 
#' @author 
#' Update definitions by \cite{Lee2001}. 
#' 
#' C++ optimised implementation by Renaud Gaujoux. 
#' 
#' @rdname nmf_update_euclidean
#' @aliases nmf_update.euclidean
#' @export
nmf_update.euclidean.h <- std.euclidean.update.h <- 
function(v, w, h, eps=10^-9, nbterms=0L, ncterms=0L, copy=TRUE){
	.Call("euclidean_update_H", v, w, h, eps, nbterms, ncterms, copy, PACKAGE='NMF')
}
#' \code{nmf_update.euclidean.w_R} and \code{nmf_update.euclidean.h_R} implement the same updates 
#' in \emph{plain R}.
#' 
#' @param wh already computed NMF estimate used to compute the denominator term. 
#' 
#' @rdname nmf_update_euclidean
#' @export
nmf_update.euclidean.h_R <- R_std.euclidean.update.h <- function(v, w, h, wh=NULL, eps=10^-9){
	# compute WH if necessary	
	den <- if( is.null(wh) ) crossprod(w) %*% h
			else{ t(w) %*% wh}
	
	# H_au = H_au (W^T V)_au / (W^T W H)_au
	pmax(h * crossprod(w,v),eps) / (den + eps);
}

#' @details
#' The basis matrix (\code{W}) is updated as follows:
#' \deqn{
#' W_ik \leftarrow \frac{\max(W_ik (V H^T)_ik, \varepsilon) }{ (W H H^T)_ik + \varepsilon}
#' }{
#' W_ik <- max(W_ik (V H^T)_ik, eps) / ( (W H H^T)_ik + eps )
#' }
#' 
#' @param weight numeric vector of sample weights, e.g., used to normalise samples 
#' coming from multiple datasets.
#' It must be of the same length as the number of samples/columns in \code{v} 
#' -- and \code{h}. 
#' 
#' @rdname nmf_update_euclidean
#' @export
nmf_update.euclidean.w <- std.euclidean.update.w <-
function(v, w, h, eps=10^-9, nbterms=0L, ncterms=0L, weight=NULL, copy=TRUE){
	.Call("euclidean_update_W", v, w, h, eps, weight, nbterms, ncterms, copy, PACKAGE='NMF')
}
#' @rdname nmf_update_euclidean
#' @export
nmf_update.euclidean.w_R <- R_std.euclidean.update.w <- function(v, w, h, wh=NULL, eps=10^-9){
	# compute WH if necessary	
	den <- if( is.null(wh) ) w %*% tcrossprod(h)
			else{ wh %*% t(h)}
	
	# W_ia = W_ia (V H^T)_ia / (W H H^T)_ia and columns are rescaled after each iteration	
	pmax(w * tcrossprod(v, h), eps) / (den + eps);
}


################################################################################################
# AFTER-UPDATE METHODS
################################################################################################

#' Stopping Criteria for NMF Iterative Strategies
#' 
#' The function documented here implement stopping/convergence criteria 
#' commonly used in NMF algorithms.
#' 
#' \code{NMFStop} acts as a factory method that create stopping criterion functions 
#' to be used with \code{\link{nmf}}, from different types of values.
#' 
#' @details
#' \code{NMFStop} returns functions unchanged, integer values are used to 
#' create a stopping criterion of that number of iterations via \code{nmf.stop.iteration}, 
#' numeric values are used to create a stopping criterion of that stationary threshold 
#' via \code{nmf.stop.threshold}.
#' Character strings are assumed to be access keys for registered criteria (currently 
#' available: \dQuote{connectivity} and \dQuote{stationary}), or a function name in the 
#' global environment or the namespace of the loading package.
#' 
#' @param val access key that can be a character string, a single integer or 
#' numeric, or a function.
#' 
#' @return a function that can be passed to argument \code{.stop} of function 
#' \code{\link{nmf}}.
#' 
#' @aliases stop-NMF
#' @rdname stop-NMF
#' @export
NMFStop <- function(val){
	
	key <- val
	if( is.integer(key) )	nmf.stop.iteration(key)
	else if( is.numeric(key) ) nmf.stop.threshold(key)
	else if( is.function(key) ) key
	else if( is.character(key) ){
		# update .stop for back compatibility:
		if( key == 'nmf.stop.consensus') key <- 'connectivity'
		
		# first lookup for a `nmf.stop.*` function
		key2 <- paste('nmf.stop.', key, sep='')
		e <- pkgmaker::packageEnv()
		sfun <- getFunction(key2, mustFind=FALSE, where = e)
		if( is.null(sfun) ) # lookup for the function as such
			sfun <- getFunction(key, mustFind = FALSE, where = e)			
		if( is.null(sfun) )
			stop("Invalid key ['", key,"']: could not find functions '",key2, "' or '", key, "'")
		sfun
	}else if( identical(val, FALSE) ) # create a function that does not stop 
		function(strategy, i, target, data, ...){FALSE}
	else
		stop("Invalid key: should be a function, a character string or a single integer/numeric value. See ?NMFStop.")	
}

#' \code{nmf.stop.iteration} generates a function that implements the stopping 
#' criterion that limits the number of iterations to a maximum of \code{n}),
#' i.e. that returns \code{TRUE} if \code{i>=n}, \code{FALSE} otherwise.
#' 
#' @param n maximum number of iteration to perform.
#'   
#' @return a function that can be used as a stopping criterion for NMF algorithms 
#' defined as \code{\linkS4class{NMFStrategyIterative}} objects. 
#' That is a function with arguments \code{(strategy, i, target, data, ...)} 
#' that returns \code{TRUE} if the stopping criterion is satisfied -- which in 
#' turn stops the iterative process, and \code{FALSE} otherwise.
#'   
#' @export
#' @family NMFStrategyIterative
#' @rdname stop-NMF
nmf.stop.iteration <- function(n){
	
	nmf.debug("Using stopping criterion - Fixed number of iterations: ", n)
	if( !is.numeric(n) )
		stop("Invalid argument `n`: must be an integer value")
	if( length(n) > 1 )
		warning("NMF::nmf - Argument `n` [", deparse(substitute(n)), "] has length > 1: only using the first element.")
	
	.max <- n[1]
	function(object, i, y, x, ...) i >= .max
}

#' \code{nmf.stop.threshold} generates a function that implements the stopping 
#' criterion that stops when a given stationarity threshold is achieved by 
#' successive iterations. 
#' The returned function is identical to \code{nmf.stop.stationary}, but with 
#' the default threshold set to \code{threshold}.
#' 
#' @param threshold default stationarity threshold  
#' 
#' @export
#' @rdname stop-NMF
nmf.stop.threshold <- function(threshold){	
	
	nmf.debug("Using stopping criterion - Stationarity threshold: ", threshold)
	if( !is.numeric(threshold) )
		stop("Invalid argument `threshold`: must be a numeric value")
	if( length(threshold) > 1 )
		warning("NMF::nmf - Argument `threshold` [", deparse(substitute(threshold)), "] has length > 1: only using the first element.")
	
	eval(parse(text=paste("function(strategy, i, target, data, stationary.th=", threshold, ", ...)
		nmf.stop.stationary(strategy, i, target, data, stationary.th=stationary.th, ...)")))
}


#' \code{nmf.stop.stationary} implements the stopping criterion of stationarity 
#' of the objective value, which stops when the objective value does not change 
#' anymore with further iterations.
#' Objective values are compared at given regular intervals and using a given 
#' threshold.
#' 
#' @param object an NMF strategy object
#' @param i the current iteration
#' @param y the target matrix
#' @param x the current NMF model 
#' @param stationary.th maximum difference for two objective values to be 
#' considered equal.
#' @param check.interval interval (in number of iterations) on which stationarity 
#' is computed. 
#' @param ... extra arguments passed to the function \code{\link{objective}}, 
#' which computes the objective value between \code{x} and \code{y}.
#' 
#' @export
#' @rdname stop-NMF
nmf.stop.stationary <- function(object, i, y, x, stationary.th=10^-6, check.interval=10, ...){
		
		if( i == 0L ){ # initialisation call: compute initial objective value
			current.value <- deviance(object, x, y, ...)
			# check for NaN, i.e. probably infinitely small value (cf. bug reported by Nadine POUKEN SIEWE)
			if( is.nan(current.value) ) return(TRUE)
			# store value in workspace for later calls
			staticVar('objective.value', current.value, init=TRUE)
			return(FALSE)
		}
		
		# test convergence only every 10 iterations
		if( i %% check.interval != 0 ) return( FALSE );
		
		# get last objective value from workspace		
		last.value <- staticVar('objective.value')
		current.value <- deviance(object, x, y, ...)
		# check for NaN, i.e. probably infinitely small value (cf. bug reported by Nadine POUKEN SIEWE)
		if( is.nan(current.value) ) return(TRUE)
		# if the relative decrease in the objective value is to small then stop
		if( abs( (last.value - current.value)/check.interval ) <= stationary.th ) return( TRUE )
		
		# update the objective value
		staticVar('objective.value', current.value)
		
		# do NOT stop
		FALSE
}
#' \code{nmf.stop.connectivity} implements the stopping criterion that is based 
#' on the stationarity of the connectivity matrix.
#' 
#' @param stopconv number of iterations intervals over which the connectivity 
#' matrix must not change for stationarity to be achieved.
#'   
#' @export
#' @rdname stop-NMF
nmf.stop.connectivity <- function(object, i, y, x, stopconv=40, check.interval=10, ...){

		if( i == 0L ){ # initialisation call
			# Initialize consensus variables 
			# => they are static variables within the strategy's workspace so that
			# they are persistent and available throughout across the calls
			h <- coef(x)
			staticVar('consold', matrix(0, ncol(h), ncol(h)), init=TRUE)
			staticVar('inc', 0, init=TRUE)
			return(FALSE)
		}
	
		# test convergence only every 10 iterations
		if( i %% check.interval != 0 ) return( FALSE );
		
		# retrieve metaprofiles
		h <- coef(x, all=FALSE)
		
		# retrieve the last values of the consensus variables
		consold <- staticVar('consold')
		inc <- staticVar('inc')
						
		# construct connectivity matrix
		index <- apply(h, 2, function(x) which.max(x) )
		cons <- outer(index, index, function(x,y) ifelse(x==y, 1,0));

		changes <- cons!=consold
		if( !any(changes) ) inc <- inc+1 # connectivity matrix has not changed: increment the count
		else{
			consold=cons;			
			inc <- 0;                         # else restart counting
		}
		
		# prints number of changing elements 		
		#if( verbose(x) ) cat( sprintf('%d ', sum(changes)) ) 
		#cat( sprintf('%d ', sum(changes)) )
										
		# assume convergence is connectivity stops changing 
		if(inc>stopconv) return( TRUE );
		
		# update the consensus variables in the workspace
		staticVar('consold', consold)
		staticVar('inc', inc)
		
		# do NOT stop
		FALSE
}


################################################################################################
# WRAP-UP METHODS
################################################################################################
