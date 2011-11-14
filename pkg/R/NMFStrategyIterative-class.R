###% @include NMFStrategy-class.R
###% @include NMF-class.R
###% @include NMFOffset-class.R
###% @include NMFns-class.R
###% @include registry.R
NA

###% NMFStrategyIterative class definition
###%
###% An NMFStrategyIterative is the implementation of strategy design-pattern for NMF algorithms.
###% It implements the following interface:
###% - Initialization of the NMF object
###% - Update of variables at each iteration
###% - Stop specific task
###% - WrapUp the NMF object
###%
###% @author Renaud Gaujoux

#Initialise verbosity with VComments 
#V1# threshold=0

###% Base class to define NMF algorithms.
###%
###% This class defines the common strategy interface used by most NMF algorithms.
###%
###% @slot Update the update method that compute the values of the factors and parameters at each iteration.
###%
###% @slot Stop the stop method that implement the algorithm's stopping criteria.
###%
###% @slot WrapNMF the method that wrap up the result as an NMFClass instance
###%
setClass('NMFStrategyIterative'
	, representation(
				Update = '.functionSlot', # update method	
				Stop = '.functionSlot.null', # method called just after the update
				WrapNMF = '.functionSlot.null' # method called just before returning the resulting NMF object
				)	
  , prototype=prototype(
				Update = ''
				, Stop = NULL
				, WrapNMF = NULL
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


setMethod('show', 'NMFStrategyIterative',
	function(object){
		
		#cat('<object of class: NMFStrategyIterative>')
		callNextMethod()
		cat("<Iterative schema:>\n")
		# go through the slots
		#s.list <- getSlots('NMFStrategyIterative')
		#s.list <- s.list[s.list=='ANY']
		s.list <- c('Update', 'Stop', 'WrapNMF')
		names(s.list) <- s.list
		sapply(names(s.list), function(sname){
					svalue <- slot(object,sname)
					cat(sname, ": "
						, if( is.function(svalue) ) capture.output(print(args(svalue)))[1] else paste("'", svalue,"'", sep='')
						, "\n")
				})
		
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
	slot(strategyX, 'WrapNMF') <- preload.slot(strategyX, 'WrapNMF', identity)
	
	# load the objective function
	objective(strategyX) <- distance(method=objective(strategy))

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
.nmf.plugin.core <- function(){
	
	list(
		# Brunet
		new('NMFStrategyIterative', name='brunet', objective='KL'
					, Update='nmf.update.brunet'
					, Stop='connectivity'
			)

		# Lee	
		, new('NMFStrategyIterative', name='lee', objective='euclidean'
					, Update='nmf.update.lee'
					, Stop='connectivity'
			)
		
		# NMF with offset
		, new('NMFStrategyIterative', name='offset', objective='euclidean'
					, model = 'NMFOffset'
					, Update='nmf.update.offset'
					, Stop='connectivity'
			)
			
		# nsNMF
		, new('NMFStrategyIterative', name='nsNMF', objective='KL'
					, model='NMFns'
					, Update='nmf.update.ns'
					, Stop='connectivity'
			)
	)
}

###% Hook to initialize old R version built-in iterative methods
.nmf.plugin.core_R <- function(){
	
	list(
			# Brunet
			new('NMFStrategyIterative', name='.R#brunet', objective='KL'
					, Update='R_nmf.update.brunet'
					, Stop='connectivity'
			)
			
			# Lee	
			, new('NMFStrategyIterative', name='.R#lee', objective='euclidean'
					, Update='R_nmf.update.lee'
					, Stop='connectivity'
			)			
	
			# NMF with offset
			, new('NMFStrategyIterative', name='.R#offset', objective='euclidean'
					, model = 'NMFOffset'
					, Update='R_nmf.update.offset'
					, Stop='connectivity'
			)
			
			# nsNMF
			, new('NMFStrategyIterative', name='.R#nsNMF', objective='KL'
					, model='NMFns'
					, Update='R_nmf.update.ns'
					, Stop='connectivity'
			)
	
	)
}

#' Get/Set a Static Variable in NMF Algorithms
#' 
#' @param name Name of the static variable (as a single character string)
#' @param value New value of the static variable
#' @param init a logical used when a \code{value} is provided, that specifies 
#' if the variable should be set to the new value only if it does not exist yet 
#' (\code{init=TRUE}). 
#' @return The static variable's value
#' @export
staticVar <- function(){
	
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
}
staticVar <- staticVar()

setMethod('run', signature(method='NMFStrategyIterative', x='matrix', seed='NMFfit'),
	function(method, x, seed, .stop=NULL, ...){
	
	no.stop <- identical(.stop, FALSE)
	# override the stop method on runtime
	if( !is.null(.stop) ){
		method@Stop <- NMFStop(.stop)
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
	
	# call the xified startegy's run method				
	if( no.stop ){
		# issue a message for the user to manually stop the algorithm
		message("#######\n# NMF - Using no stopping criterion: press Ctrl+C to stop\n#######")
		run(strategyX, x, seed, maxIter=Inf, ...)
	}else
		run(strategyX, x, seed, ...)
})

###% Generic algorithm for NMF, based on NMFStrategyIterativeX object.
setMethod('run', signature(method='NMFStrategyIterativeX', x='matrix', seed='NMFfit'),
	function(method, x, seed, maxIter=2000, ...){
				
	strategy <- method
	v <- x
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
	for( i in 1:maxIter ){
		
		if( verbose && (i %% 50 == 0) ) cat('', i)
		
		#Vc# update the matrices
		nmfFit <- updateFun(i, v, nmfFit, ...)
		
		#Vc# Stopping criteria
		# give the strategy the opportunity to perform stuff after the update: modify the data and/or stop iteration
		stop.signal <- stopFun(strategy, i, v, nmfFit, ...)
	
		# every now and then track the error if required
		nmfData <- trackError(nmfData, objective(strategy, v, nmfFit, ...), i)
		
		# if the strategy ask for stopping, then stop the iteration
		if( stop.signal ) break;
				
	}
	if( verbose ) cat("\nDONE (stopped at ",i,'/', maxIter," iterations)\n", sep='')
	
	# force to compute last error if not already done
	nmfData <- trackError(nmfData, objective(strategy, v, nmfFit, ...), i, force=TRUE)
	
	# store the fitted model
	fit(nmfData) <- nmfFit
	
	#Vc# wrap up
	# let the strategy build the result
	nmfData <- strategy@WrapNMF(nmfData)
	if( !inherits(nmfData, 'NMFfit') ) stop('NMFStrategyIterative[', name(strategy), ']::WrapNMF did not return a "NMF" instance [returned: "', class(nmfData), '"]')
	
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

###% Standard multiplicative update for matrix \code{H} (i.e. the second factor) in a divergence based NMF.
###% 
###% The matrix \code{H} is updated as follows:
###% \deqn{%
###% H_{ij} \leftarrow H_{ij}  \frac{\left( sum_k \frac{W_ki V_kj}{(WH)_kj} \right)}{ sum_k W_ka }.%
###% }
###% 
###% @refrences ﻿Lee, D., & Seung, H. (2001)
###% , Algorithms for non-negative matrix factorization
###% , Advances in neural information processing systems,
###% http://scholar.google.com/scholar?q=intitle:Algorithms+for+non-negative+matrix+factorization#0
###% 
R_std.divergence.update.h <- function(v, w, h, wh=NULL)
{	
	# compute WH if necessary	
	if( is.null(wh) ) wh <- w %*% h
	
	# divergence-reducing NMF iterations
	# H_au = H_au ( sum_i [ W_ia V_iu / (WH)_iu ] ) / ( sum_k W_ka ) -> each row of H is divided by a the corresponding colSum of W
	h * crossprod(w, v / wh) / colSums(w)	
}
std.divergence.update.h <- function(v, w, h, copy=TRUE)
{	
	.Call("divergence_update_H", v, w, h, copy)
}

###% Standard multiplicative update for matrix \code{W} (i.e. the second factor) in a divergence based NMF.
###% 
###% The matrix \code{W} is updated as follows:
###% \deqn{%
###% W_ij \leftarrow W_ij \frac{ sum_k [H_jk A_ik / (WH)_ik ] }{sum_k H_jk } %
###% }
###% 
###% @refrences
###% Lee, D..D., and Seung, H.S., (2001), 'Algorithms for Non-negative Matrix Factorization'
###% , Advances in neural information processing systems 13, 556-562.
###% , http://scholar.google.com/scholar?q=intitle:Algorithms+for+non-negative+matrix+factorization#0
###% 
R_std.divergence.update.w <- function(v, w, h, wh=NULL)
{			
	# compute WH if necessary	
	if( is.null(wh) ) wh <- w %*% h
	
	# W_ia = W_ia ( sum_u [H_au A_iu / (WH)_iu ] ) / ( sum_v H_av ) -> each column of W is divided by a the corresponding rowSum of H
	#x2 <- matrix(rep(rowSums(h), nrow(w)), ncol=ncol(w), byrow=TRUE); 
	#w * tcrossprod(v / wh, h) / x2;
	sweep(w * tcrossprod(v / wh, h), 2L, rowSums(h), "/", check.margin = FALSE) # optimize version?
	
}
std.divergence.update.w <- function(v, w, h, copy=TRUE)
{	
	.Call("divergence_update_W", v, w, h, copy)
}

###% Computes the Nonegative Matrix Factorization of a matrix.
###%
###% This software and its documentation are copyright 2004 by the
###% Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
###% This software is supplied without any warranty or guaranteed support whatsoever. 
###% Neither the Broad Institute nor MIT can not be responsible for its use, misuse, 
###% or functionality. 
###%
###% @param v N (genes) x M (samples) original matrix to be factorized.
###%           Numerical data only. 
###%           Must be non negative. 
###%           Not all entries in a row can be 0. If so, add a small constant to the 
###%           matrix, eg.v+0.01*min(min(v)),and restart.
###%           
###% @param r the number of desired factors (i.e. the rank of factorization)
###% @param verbose prints iteration count and changes in connectivity matrix elements unless verbose is 0 
###% @return The two factors of the factorization:
###% \item{w }{N x r NMF factor}
###% \item{h }{r x M NMF factor}
###%
###% @note NMF iterations stop when connectivity matrix has not changed 
###%        for 10*stopconv interations. This is experimental and can be
###%        adjusted.
###% @author Jean-Philippe Brunet \email{brunet@@broad.mit.edu}
###% @author Renaud Gaujoux \email{renaud@@cbio.uct.ac.za}
###% @references
###% Metagenes and molecular pattern discovery using matrix factorization
###% , Brunet, J.~P., Tamayo, P., Golub, T.~R., and Mesirov, J.~P. (2004)
###% , Proc Natl Acad Sci U S A
###% , 101(12)
###% , 4164--4169.
###% NMF divergence update equations :
###% Lee, D..D., and Seung, H.S., (2001), 'Algorithms for Non-negative Matrix 
###% Factorization', Adv. Neural Info. Proc. Syst. 13, 556-562.
R_nmf.update.brunet <- function(i, v, data, eps=.Machine$double.eps, ...)
{
	# retrieve each factor
	w <- basis(data); h <- coef(data);
	
	# standard divergence-reducing NMF update for H
	h <- R_std.divergence.update.h(v, w, h)
	
	# standard divergence-reducing NMF update for W
	w <- R_std.divergence.update.w(v, w, h)
	
	#every 10 iterations: adjust small values to avoid underflow 
	if( i %% 10 == 0 ){
		#precision threshold for numerical stability
		#eps <- .Machine$double.eps
		h[h<eps] <- eps;
		w[w<eps] <- eps;
	}
	
	#return the modified data
	basis(data) <- w; coef(data) <- h;	
	return(data)
	
}
nmf.update.brunet <- function(i, v, data, copy=FALSE, eps=.Machine$double.eps, ...)
{
	# retrieve each factor
	w <- basis(data); h <- coef(data);	
	
	# standard divergence-reducing NMF update for H	
	h <- std.divergence.update.h(v, w, h, copy=copy)
		
	# standard divergence-reducing NMF update for W
	w <- std.divergence.update.w(v, w, h, copy=copy)
	
	#every 10 iterations: adjust small values to avoid underflow
	# NB: one adjusts in place even when copy=TRUE, as 'h' and 'w' are local variables
	if( i %% 10 == 0 ){
		#eps <- .Machine$double.eps
		h <- pmin.inplace(h, eps)
		w <- pmin.inplace(w, eps)
	}
	
	# update object if the updates duplicated the data
	if( copy ){		
		#return the modified data	
		basis(data) <- w; 
		coef(data) <- h;
	}
	return(data)
	
}

###% Updates for Euclidean norm reduction
###% based on Lee and Seung algorithm
###% 
R_std.euclidean.update.h <- function(v, w, h, wh=NULL, eps=10^-9){
	# compute WH if necessary	
	den <- if( is.null(wh) ) crossprod(w) %*% h
			else{ t(w) %*% wh}
	
	# H_au = H_au (W^T V)_au / (W^T W H)_au
	pmax(h * crossprod(w,v),eps) / (den + eps);
}
std.euclidean.update.h <- function(v, w, h, eps=10^-9, copy=TRUE){
	.Call("euclidean_update_H", v, w, h, eps, copy)
}
# with offset
offset.std.euclidean.update.h <- function(v, w, h, offset, eps=10^-9, copy=TRUE){
	.Call("offset_euclidean_update_H", v, w, h, offset, eps, copy)
}

R_std.euclidean.update.w <- function(v, w, h, wh=NULL, eps=10^-9){
	# compute WH if necessary	
	den <- if( is.null(wh) ) w %*% tcrossprod(h)
			else{ wh %*% t(h)}
	
	# W_ia = W_ia (V H^T)_ia / (W H H^T)_ia and columns are rescaled after each iteration	
	pmax(w * tcrossprod(v, h), eps) / (den + eps);
}
std.euclidean.update.w <- function(v, w, h, eps=10^-9, copy=TRUE){
	.Call("euclidean_update_W", v, w, h, eps, copy)
}
# with offset
offset.std.euclidean.update.w <- function(v, w, h, offset, eps=10^-9, copy=TRUE){
	.Call("offset_euclidean_update_W", v, w, h, offset, eps, copy)
}

###% Multiplicative update for reducing the euclidean distance.
###%
###% 
R_nmf.update.lee <- function(i, v, data, rescale=TRUE, eps=10^-9, ...)
{
	# retrieve each factor
	w <- basis(data); h <- coef(data);	
	
	#precision threshold for numerical stability
	#eps <- 10^-9
	
	# compute the estimate WH
	#wh <- estimate(data)
	
	# euclidean-reducing NMF iterations	
	# H_au = H_au (W^T V)_au / (W^T W H)_au
	#h <- pmax(h * (t(w) %*% v),eps) / ((t(w) %*% w) %*% h + eps);
	h <- R_std.euclidean.update.h(v, w, h, eps=eps)
	
	# update H and recompute the estimate WH
	#metaprofiles(data) <- h
	#wh <- estimate(data)

	# W_ia = W_ia (V H^T)_ia / (W H H^T)_ia and columns are rescaled after each iteration	
	#w <- pmax(w * (v %*% t(h)), eps) / (w %*% (h %*% t(h)) + eps);
	w <- R_std.euclidean.update.w(v, w, h, eps=eps)
	#rescale columns TODO: effect of rescaling? the rescaling makes the update with offset fail
	if( rescale ) w <- sweep(w, 2L, colSums(w), "/", check.margin=FALSE)
	
	#return the modified data
	basis(data) <- w; coef(data) <- h;	
	return(data)
}

nmf.update.lee <- function(i, v, data, rescale=TRUE, copy=FALSE, eps=10^-9, ...)
{
	# retrieve each factor
	w <- basis(data); h <- coef(data);	
	
	#precision threshold for numerical stability
	#eps <- 10^-9
	
	# compute the estimate WH
	#wh <- estimate(data)
	
	# euclidean-reducing NMF iterations	
	# H_au = H_au (W^T V)_au / (W^T W H)_au
	h <- std.euclidean.update.h(v, w, h, eps=eps, copy=copy)
	# update original object if not modified in place
	if( copy ) coef(data) <- h
	
	# W_ia = W_ia (V H^T)_ia / (W H H^T)_ia and columns are rescaled after each iteration	
	w <- std.euclidean.update.w(v, w, h, eps=eps, copy=copy)
	#rescale columns TODO: effect of rescaling? the rescaling makes the update with offset fail
	if( rescale ) w <- sweep(w, 2L, colSums(w), "/", check.margin=FALSE)
	
	#return the modified data
	basis(data) <- w; 	
	return(data)
}

###% Multiplicative update for reducing the euclidean distance including on offset.
###%
###% The method is a modified version of Lee's method that also fits an offset vector which model a common expression baseline for each gene accross all samples.
R_nmf.update.offset <- function(i, v, data, eps=10^-9, ...)
{	
	# retrieve each factor
	w <- basis(data); h <- coef(data);
	# retrieve offset and fill it if necessary (with mean of rows)
	off <- offset(data)
	if( i == 1 && length(off) == 0 )
		off <- rowMeans(v)
	
	#precision threshold for numerical stability
	#eps <- 10^-9
	
	# compute standard lee update (it will take the offset into account) without rescaling W's columns
	
	h <- R_std.euclidean.update.h(v, w, h, wh=w%*%h + off, eps=eps)
	w <- R_std.euclidean.update.w(v, w, h, wh=w%*%h + off, eps=eps)
	#data <- nmf.update.lee(i, v, data, rescale=FALSE, ...)
	
	# update the offset	
	# V0_i = V0_i ( sum_j V_ij ) / ( sum_j (V.off + W H)_ij )
	data@offset <- off * pmax(rowSums(v), eps) / (rowSums(w%*%h + off) + eps)
		
	#return the modified data
	basis(data) <- w; coef(data) <- h;
	return(data)
}

nmf.update.offset <- function(i, v, data, copy=FALSE, eps=10^-9, ...)
{	
	# retrieve each factor
	w <- basis(data); h <- coef(data);
	# retrieve offset and fill it if necessary (with mean of rows)
	off <- offset(data)
	if( i == 1 && length(off) == 0 )
		off <- rowMeans(v)
	
	#precision threshold for numerical stability
	#eps <- 10^-9
	
	# compute standard offset updates
	h <- offset.std.euclidean.update.h(v, w, h, off, eps=eps, copy=copy)
	w <- offset.std.euclidean.update.w(v, w, h, off, eps=eps, copy=copy)
	
	# update the offset	
	# V0_i = V0_i ( sum_j V_ij ) / ( sum_j (V.off + W H)_ij )
	data@offset <- off * pmax(rowSums(v), eps) / (rowSums(w%*%h + off) + eps)	
	
	# update the original object if not modified in place
	if( copy ){ 
		basis(data) <- w; 
		coef(data) <- h;
	}
	return(data)
}

###% Multiplicative update for Nonsmooth Nonnegative Matrix Factorization (nsNMF).
###%
###% The update rules are essentialy the same as in Brunet, but WH is replaced by WSH 
###% whereas W (resp. H) is replaced by WS (resp. SH) in the update of H (resp. of W).
###%
###% @references 
###% Alberto Pascual-Montano et al. (2006), 'Nonsmooth Nonnegative Matrix Factorization (nsNMF)'
###% , IEEE Transactions On Pattern Analysis And Machine Intelligence, Vol. 28, No. 3, March 2006 403
###%
nmf.update.ns <- function(i, v, data, copy=FALSE, ...)
{
	# retrieve and alter the factors for updating H
	S <- smoothing(data)
	w <- basis(data)
	h <- coef(data);
	
	# standard divergence-reducing update for H with modified W
	h <- std.divergence.update.h(v, w %*% S, h, copy=copy)
	
	# update H if not modified in place
	if( copy ) coef(data) <- h
	
	# standard divergence-reducing update for W with modified H
	w <- std.divergence.update.w(v, w, S %*% h, copy=copy)
	
	# rescale columns of W
	w <- sweep(w, 2L, colSums(w), '/', check.margin=FALSE)
	
	#return the modified data
	basis(data) <- w;
	return(data)
}

R_nmf.update.ns <- function(i, v, data, ...)
{
	# retrieve and alter the factors for updating H
	S <- smoothing(data)
	w <- basis(data)
	#w <- metagenes(data) %*% smoothing(fit(data)); # W <- WS
	h <- coef(data);
	
	# compute the estimate WH
	#wh <- estimate(data, W=w.init, H=h, S=S)
	
	# standard divergence-reducing update for H with modified W
	h <- R_std.divergence.update.h(v, w %*% S, h)
	
	# update H and recompute the estimate WH
	coef(data) <- h
	# retrieve and alter the factors for updating W
	#w <- tmp;
	#h <- smoothing(fit(data)) %*% metaprofiles(data); # H <- SH
	#h <- S %*% h; # H <- SH
	
	# standard divergence-reducing update for W with modified H
	w <- R_std.divergence.update.w(v, w, S %*% h)
	
	# rescale columns of W
	w <- sweep(w, 2L, colSums(w), '/', check.margin=FALSE)
	
	#return the modified data
	basis(data) <- w; #metaprofiles(data) <- h;
	return(data)
}


################################################################################################
# AFTER-UPDATE METHODS
################################################################################################

#' Stopping Criteria for NMF Iterative Strategies
#' 
#' 
#' \code{nmf.stop.iteration}: the stopping criterium is that a given number of 
#' iterations (\code{n}) are performed. That is it returns 
#' \code{TRUE} if \code{i>=n}, \code{FALSE} otherwise.
#' 
#' @param n The exact number of iteration to perform.
#'   
#' @return a function that can be used as a stopping criterion for NMF algorithms 
#' defined as \code{\linkS4class{NMFStrategyIterative}} objects. That is a function 
#' with arguments \code{(strategy, i, target, data, ...)} that returns 
#' \code{TRUE} if the stopping criterium is satisfied -- which in turn stops the 
#' iterative process, and \code{FALSE} otherwise.
#'   
#' @export
#' @aliases stop-NMF
#' @family NMFStrategyIterative
#' @rdname stop-NMF
nmf.stop.iteration <- function(n){
	
	nmf.debug("Using stopping criterion - Fixed number of iterations: ", n)
	if( !is.numeric(n) )
		stop("Invalid argument `n`: must be an integer value")
	if( length(n) > 1 )
		warning("NMF::nmf - Argument `n` [", deparse(substitute(n)), "] has length > 1: only using the first element.")
	
	.max <- n[1]
	function(strategy, i, target, data, ...) i >= .max
}

nmf.stop.threshold <- function(threshold){	
	
	nmf.debug("Using stopping criterion - Stationatiry threshold: ", threshold)
	if( !is.numeric(threshold) )
		stop("Invalid argument `threshold`: must be a numeric value")
	if( length(threshold) > 1 )
		warning("NMF::nmf - Argument `threshold` [", deparse(substitute(threshold)), "] has length > 1: only using the first element.")
	
	eval(parse(text=paste("function(strategy, i, target, data, stationary.th=", threshold, ", ...)
		nmf.stop.stationary(strategy, i, target, data, stationary.th=stationary.th, ...)")))
}

nmf.stop.stationary <- function(strategy, i, target, data, stationary.th=10^-6, check.interval=10, ...){
		
		# first call compute the initial error
		if( i == 1 ) staticVar('objective.value', objective(strategy, target, data, ...), init=TRUE)
		
		# test convergence only every 10 iterations
		if( i %% check.interval != 0 ) return( FALSE );
		
		# initialize static variables to store the error across the calls		
		last.value <- staticVar('objective.value')
		current.value <- objective(strategy, target, data, ...)
		# if the relative decrease in the objective value is to small then stop
		if( abs( (last.value - current.value)/check.interval ) <= stationary.th ) return( TRUE )
		
		# update the objective value
		staticVar('objective.value', current.value)
		
		# do NOT stop
		FALSE
}

nmf.stop.connectivity <- function(strategy, i, target, data, stopconv=40, ...){
			
		# test convergence only every 10 iterations
		if( i %% 10 != 0 ) return( FALSE );
		
		# retrieve metaprofiles
		h <- coef(data)
		
		# Initialize consensus variables 
		# => they are static variables within the strategy's workspace so that
		# they are persistent and available throughout across the calls
		staticVar('consold', matrix(0, ncol(h), ncol(h)), init=TRUE)
		staticVar('inc', 0, init=TRUE)
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
		#if( verbose(data) ) cat( sprintf('%d ', sum(changes)) ) 
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
