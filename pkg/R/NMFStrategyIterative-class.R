#' @include NMFStrategy-class.R
#' @include NMF-class.R
#' @include NMFOffset-class.R
#' @include NMFns-class.R
#' @include registry.R
NA

#' NMFStrategyIterative class definition
#'
#' An NMFStrategyIterative is the implementation of strategy design-pattern for NMF algorithms.
#' It implements the following interface:
#' - Initialization of the NMF object
#' - Update of variables at each iteration
#' - Stop specific task
#' - WrapUp the NMF object
#'
#' @author Renaud Gaujoux

#Initialise verbosity with VComments 
#V1# threshold=0

#' Base class to define NMF algorithms.
#'
#' This class defines the common strategy interface used by most NMF algorithms.
#'
#' @slot Update the update method that compute the values of the factors and parameters at each iteration.
#'
#' @slot Stop the stop method that implement the algorithm's stopping criteria.
#'
#' @slot WrapNMF the method that wrap up the result as an NMFClass instance
#'
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
			n.stop <- names(formals(object@Stop))
			
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
						, if( is.function(svalue) ) '<function>' else paste("'", svalue,"'", sep='')
						, "\n")
				})
		
		return(invisible())
	}
)
#' This class is an auxiliary class that defines the strategy's methods by directly callable functions. 
setClass('NMFStrategyIterativeX'
	, contains = 'NMFStrategyIterative'
	, representation = representation(
				workspace = 'environment' # workspace to use persistent variables accross methods
				)
)


#' Creates a NMFStrategyIterativeX object from a NMFStrategyIterative object.
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
		if( is.character(svalue) && nchar(svalue) > 0 ){
			# set the slot with the executable version of the function 
			fun <- getFunction(svalue)
		}else if( is.function(svalue) )
			fun <- svalue
		else if( !missing(default) )
			fun <- default
		else
			stop("NMFStrategyIterativeX - could not preload slot '", sname, "'")
		
		# setup a dedicated evaluation environment: use the strategy's workspace
		environment(fun) <- strategy@workspace
		
		# return the loaded function
		fun
	}
	
	# preload the function slots
	slot(strategyX, 'Update') <- preload.slot(strategyX, 'Update')
	slot(strategyX, 'Stop') <- preload.slot(strategyX, 'Stop', function(strategy, i, x, data, ...){FALSE})
	slot(strategyX, 'WrapNMF') <- preload.slot(strategyX, 'WrapNMF', function(data){data})
	
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

#' Hook to initialize built-in iterative methods when the package is loaded
.nmf.plugin.core <- function(){
	
	list(
		# Brunet
		new('NMFStrategyIterative', name='brunet', objective='KL'
					, Update='nmf.update.brunet'
					, Stop='nmf.stop.consensus'
			)

		# Lee	
		, new('NMFStrategyIterative', name='lee', objective='euclidean'
					, Update='nmf.update.lee'
					, Stop='nmf.stop.consensus'
			)
		
		# NMF with offset
		, new('NMFStrategyIterative', name='offset', objective='euclidean'
					, model = 'NMFOffset'
					, Update='nmf.update.offset'
					, Stop='nmf.stop.consensus'
			)
			
		# nsNMF
		, new('NMFStrategyIterative', name='nsNMF', objective='KL'
					, model='NMFns'
					, Update='nmf.update.ns'
					, Stop='nmf.stop.consensus'
			)
	)
}

#' Hook to initialize old R version built-in iterative methods
.nmf.plugin.core_R <- function(){
	
	list(
			# Brunet
			new('NMFStrategyIterative', name='.R#brunet', objective='KL'
					, Update='R_nmf.update.brunet'
					, Stop='nmf.stop.consensus'
			)
			
			# Lee	
			, new('NMFStrategyIterative', name='.R#lee', objective='euclidean'
					, Update='R_nmf.update.lee'
					, Stop='nmf.stop.consensus'
			)			
	
			# NMF with offset
			, new('NMFStrategyIterative', name='.R#offset', objective='euclidean'
					, model = 'NMFOffset'
					, Update='R_nmf.update.offset'
					, Stop='nmf.stop.consensus'
			)
			
			# nsNMF
			, new('NMFStrategyIterative', name='.R#nsNMF', objective='KL'
					, model='NMFns'
					, Update='R_nmf.update.ns'
					, Stop='nmf.stop.consensus'
			)
	
	)
}

setMethod('run', signature(method='NMFStrategyIterative', x='matrix', seed='NMFfit'),
	function(method, x, seed, .stop=NULL, ...){
	
	# override the stop method on runtime
	if( !is.null(.stop) )
		method@Stop = .stop
	
	# debug object in debug mode
	if( nmf.getOption('debug') ) show(method)		
	
	#Vc# Define local workspace for static variables
	# this function can be called in the methods to get/set/initialize 
	# variables that are persistent within the strategy's workspace
	.Workspace <- new.env()
	staticVar <- function(name, value, init=FALSE){		
		if( missing(value) ){
			get(name, envir=.Workspace, inherits=FALSE)
		}else{
			if( !init || !exists(name, envir=.Workspace, inherits=FALSE) )
			{
				if( init ) nmf.debug('Strategy Workspace', "initialize variable '", name, "'")
				assign(name, value, envir=.Workspace)
			}
		}
	}
	
	# runtime resolution of the strategy's functions by their names if necessary
	strategyX = xifyStrategy(method, .Workspace)
	# call the xified startegy's run method			
	run(strategyX, x, seed, ...)
})

#' Generic algorithm for NMF, based on NMFStrategyIterativeX object.
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
	required.args <- names(expected.args[which(sapply(expected.args, function(x) x==''))])
	required.args <- required.args[required.args!='...']
	
	if( any(t <- !pmatch(required.args, passed.args, nomatch=FALSE)) )
		stop("NMF::run - Update/Stop method for algorithm '", name(strategy),"': missing required argument(s) "
			, paste( paste("'", required.args[t],"'", sep=''), collapse=', '), call.=FALSE)
	
	# set default value for missing argument TODO
	missing.args <- expected.args[!pmatch(names(expected.args), passed.args, nomatch=FALSE)]
	
	
	#Vc# Start iterations
	nmfData <- seed
	nmfFit <- fit(nmfData)
	for( i in 1:maxIter ){
		
		#Vc# update the matrices
		nmfFit <- strategy@Update(i, v, nmfFit, ...)
		
		#Vc# Stopping criteria
		# give the strategy the opportunity to perform stuff after the update: modify the data and/or stop iteration
		stop.signal <- strategy@Stop(strategy, i, v, nmfFit, ...)
	
		# every now and then track the error if required
		nmfData <- trackError(nmfData, objective(strategy, v, nmfFit, ...), i)
		
		# if the strategy ask for stopping, then stop the iteration
		if( stop.signal ) break;
				
	}
	if( verbose(nmfData) ) cat("\n")
	
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

#' Standard multiplicative update for matrix \code{H} (i.e. the second factor) in a divergence based NMF.
#' 
#' The matrix \code{H} is updated as follows:
#' \deqn{%
#' H_{ij} \leftarrow H_{ij}  \frac{\left( sum_k \frac{W_ki V_kj}{(WH)_kj} \right)}{ sum_k W_ka }.%
#' }
#' 
#' @refrences ﻿Lee, D., & Seung, H. (2001)
#' , Algorithms for non-negative matrix factorization
#' , Advances in neural information processing systems,
#' http://scholar.google.com/scholar?q=intitle:Algorithms+for+non-negative+matrix+factorization#0
#' 
R_std.divergence.update.h <- function(v, w, h, wh=NULL)
{	
	# compute WH if necessary	
	if( is.null(wh) ) wh <- w %*% h
	
	# divergence-reducing NMF iterations
	# H_au = H_au ( sum_i [ W_ia V_iu / (WH)_iu ] ) / ( sum_k W_ka ) -> each row of H is divided by a the corresponding colSum of W
	h * crossprod(w, v / wh) / colSums(w)	
}
std.divergence.update.h <- function(v, w, h)
{	
	.Call("divergence_update_H", v, w, h)
}

#' Standard multiplicative update for matrix \code{W} (i.e. the second factor) in a divergence based NMF.
#' 
#' The matrix \code{W} is updated as follows:
#' \deqn{%
#' W_ij \leftarrow W_ij \frac{ sum_k [H_jk A_ik / (WH)_ik ] }{sum_k H_jk } %
#' }
#' 
#' @refrences
#' Lee, D..D., and Seung, H.S., (2001), 'Algorithms for Non-negative Matrix Factorization'
#' , Advances in neural information processing systems 13, 556-562.
#' , http://scholar.google.com/scholar?q=intitle:Algorithms+for+non-negative+matrix+factorization#0
#' 
R_std.divergence.update.w <- function(v, w, h, wh=NULL)
{			
	# compute WH if necessary	
	if( is.null(wh) ) wh <- w %*% h
	
	# W_ia = W_ia ( sum_u [H_au A_iu / (WH)_iu ] ) / ( sum_v H_av ) -> each column of W is divided by a the corresponding rowSum of H
	#x2 <- matrix(rep(rowSums(h), nrow(w)), ncol=ncol(w), byrow=TRUE); 
	#w * tcrossprod(v / wh, h) / x2;
	sweep(w * tcrossprod(v / wh, h), 2L, rowSums(h), "/", check.margin = FALSE) # optimize version?
	
}
std.divergence.update.w <- function(v, w, h)
{
	.Call("divergence_update_W", v, w, h)
}

#' Computes the Nonegative Matrix Factorization of a matrix.
#'
#' This software and its documentation are copyright 2004 by the
#' Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
#' This software is supplied without any warranty or guaranteed support whatsoever. 
#' Neither the Broad Institute nor MIT can not be responsible for its use, misuse, 
#' or functionality. 
#'
#' @param v N (genes) x M (samples) original matrix to be factorized.
#'           Numerical data only. 
#'           Must be non negative. 
#'           Not all entries in a row can be 0. If so, add a small constant to the 
#'           matrix, eg.v+0.01*min(min(v)),and restart.
#'           
#' @param r the number of desired factors (i.e. the rank of factorization)
#' @param verbose prints iteration count and changes in connectivity matrix elements unless verbose is 0 
#' @return The two factors of the factorization:
#' \item{w }{N x r NMF factor}
#' \item{h }{r x M NMF factor}
#'
#' @note NMF iterations stop when connectivity matrix has not changed 
#'        for 10*stopconv interations. This is experimental and can be
#'        adjusted.
#' @author Jean-Philippe Brunet \email{brunet@@broad.mit.edu}
#' @author Renaud Gaujoux \email{renaud@@cbio.uct.ac.za}
#' @references
#' Metagenes and molecular pattern discovery using matrix factorization
#' , Brunet, J.~P., Tamayo, P., Golub, T.~R., and Mesirov, J.~P. (2004)
#' , Proc Natl Acad Sci U S A
#' , 101(12)
#' , 4164--4169.
#' NMF divergence update equations :
#' Lee, D..D., and Seung, H.S., (2001), 'Algorithms for Non-negative Matrix 
#' Factorization', Adv. Neural Info. Proc. Syst. 13, 556-562.
R_nmf.update.brunet <- function(i, v, data, ...)
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
		eps <- .Machine$double.eps
		h[h<eps] <- eps;
		w[w<eps] <- eps;
	}
	
	#return the modified data
	basis(data) <- w; coef(data) <- h;	
	return(data)
	
}
nmf.update.brunet <- function(i, v, data, ...)
{
	# retrieve each factor
	w <- basis(data); h <- coef(data);
	
	# standard divergence-reducing NMF update for H
	h <- std.divergence.update.h(v, w, h)
	
	# standard divergence-reducing NMF update for W
	w <- std.divergence.update.w(v, w, h)
	
	#every 10 iterations: adjust small values to avoid underflow 
	if( i %% 10 == 0 ){
		#precision threshold for numerical stability
		eps <- .Machine$double.eps
		h[h<eps] <- eps;
		w[w<eps] <- eps;
	}
		
	#return the modified data
	basis(data) <- w; coef(data) <- h;	
	return(data)
	
}

#' Updates for Euclidean norm reduction
#' based on Lee and Seung algorithm
#' 
R_std.euclidean.update.h <- function(v, w, h, wh=NULL, eps=10^-9){
	# compute WH if necessary	
	den <- if( is.null(wh) ) crossprod(w) %*% h
			else{ t(w) %*% wh}
	
	# H_au = H_au (W^T V)_au / (W^T W H)_au
	pmax(h * crossprod(w,v),eps) / (den + eps);
}
std.euclidean.update.h <- function(v, w, h, eps=10^-9){
	.Call("euclidean_update_H", v, w, h, eps)
}
# with offset
offset.std.euclidean.update.h <- function(v, w, h, offset, eps=10^-9){
	.Call("offset_euclidean_update_H", v, w, h, offset, eps)
}

R_std.euclidean.update.w <- function(v, w, h, wh=NULL, eps=10^-9){
	# compute WH if necessary	
	den <- if( is.null(wh) ) w %*% tcrossprod(h)
			else{ wh %*% t(h)}
	
	# W_ia = W_ia (V H^T)_ia / (W H H^T)_ia and columns are rescaled after each iteration	
	pmax(w * tcrossprod(v, h), eps) / (den + eps);
}
std.euclidean.update.w <- function(v, w, h, eps=10^-9){
	.Call("euclidean_update_W", v, w, h, eps)
}
# with offset
offset.std.euclidean.update.w <- function(v, w, h, offset, eps=10^-9){
	.Call("offset_euclidean_update_W", v, w, h, offset, eps)
}

#' Multiplicative update for reducing the euclidean distance.
#'
#' 
R_nmf.update.lee <- function(i, v, data, rescale=TRUE, ...)
{
	# retrieve each factor
	w <- basis(data); h <- coef(data);	
	
	#precision threshold for numerical stability
	eps <- 10^-9
	
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

nmf.update.lee <- function(i, v, data, rescale=TRUE, ...)
{
	# retrieve each factor
	w <- basis(data); h <- coef(data);	
	
	#precision threshold for numerical stability
	eps <- 10^-9
	
	# compute the estimate WH
	#wh <- estimate(data)
	
	# euclidean-reducing NMF iterations	
	# H_au = H_au (W^T V)_au / (W^T W H)_au
	h <- std.euclidean.update.h(v, w, h, eps=eps)
	
	# W_ia = W_ia (V H^T)_ia / (W H H^T)_ia and columns are rescaled after each iteration	
	w <- std.euclidean.update.w(v, w, h, eps=eps)
	#rescale columns TODO: effect of rescaling? the rescaling makes the update with offset fail
	if( rescale ) w <- sweep(w, 2L, colSums(w), "/", check.margin=FALSE)
	
	#return the modified data
	basis(data) <- w; coef(data) <- h;	
	return(data)
}

#' Multiplicative update for reducing the euclidean distance including on offset.
#'
#' The method is a modified version of Lee's method that also fits an offset vector which model a common expression baseline for each gene accross all samples.
R_nmf.update.offset <- function(i, v, data, ...)
{	
	# retrieve each factor
	w <- basis(data); h <- coef(data);
	# retrieve offset and fill it if necessary (with mean of rows)
	off <- offset(data)
	if( i == 1 && length(off) == 0 )
		off <- rowMeans(v)
	
	#precision threshold for numerical stability
	eps <- 10^-9
	
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

nmf.update.offset <- function(i, v, data, ...)
{	
	# retrieve each factor
	w <- basis(data); h <- coef(data);
	# retrieve offset and fill it if necessary (with mean of rows)
	off <- offset(data)
	if( i == 1 && length(off) == 0 )
		off <- rowMeans(v)
	
	#precision threshold for numerical stability
	eps <- 10^-9
	
	# compute standard offset updates
	h <- offset.std.euclidean.update.h(v, w, h, off, eps=eps)
	w <- offset.std.euclidean.update.w(v, w, h, off, eps=eps)
	
	# update the offset	
	# V0_i = V0_i ( sum_j V_ij ) / ( sum_j (V.off + W H)_ij )
	data@offset <- off * pmax(rowSums(v), eps) / (rowSums(w%*%h + off) + eps)	
	
	#return the modified data
	basis(data) <- w; coef(data) <- h;
	return(data)
}

#' Multiplicative update for Nonsmooth Nonnegative Matrix Factorization (nsNMF).
#'
#' The update rules are essentialy the same as in Brunet, but WH is replaced by WSH 
#' whereas W (resp. H) is replaced by WS (resp. SH) in the update of H (resp. of W).
#'
#' @references 
#' Alberto Pascual-Montano et al. (2006), 'Nonsmooth Nonnegative Matrix Factorization (nsNMF)'
#' , IEEE Transactions On Pattern Analysis And Machine Intelligence, Vol. 28, No. 3, March 2006 403
#'
nmf.update.ns <- function(i, v, data, ...)
{
	# retrieve and alter the factors for updating H
	S <- smoothing(data)
	w <- basis(data)
	#w <- metagenes(data) %*% smoothing(fit(data)); # W <- WS
	h <- coef(data);
	
	# compute the estimate WH
	#wh <- estimate(data, W=w.init, H=h, S=S)
		
	# standard divergence-reducing update for H with modified W
	h <- std.divergence.update.h(v, w %*% S, h)
	
	# update H and recompute the estimate WH
	coef(data) <- h
	# retrieve and alter the factors for updating W
	#w <- tmp;
	#h <- smoothing(fit(data)) %*% metaprofiles(data); # H <- SH
	#h <- S %*% h; # H <- SH
	
	# standard divergence-reducing update for W with modified H
	w <- std.divergence.update.w(v, w, S %*% h)
	
	# rescale columns of W
	w <- sweep(w, 2L, colSums(w), '/', check.margin=FALSE)
	
	#return the modified data
	basis(data) <- w; #metaprofiles(data) <- h;
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

nmf.stop.stationary <- function(strategy, i, target, data, stationary.th=10^-6, ...){
		
		# first call compute the initial error
		if( i == 1 ) staticVar('objective.value', objective(strategy, target, data, ...), init=TRUE)
		
		# test convergence only every 10 iterations
		interval <- 10
		if( i %% interval != 0 ) return( FALSE );
		
		# initialize static variables to store the error across the calls		
		last.value <- staticVar('objective.value')
		current.value <- objective(strategy, target, data, ...)
		# if the relative decrease in the objective value is to small then stop
		if( abs( (last.value - current.value)/interval ) <= stationary.th ) return( TRUE )
		
		# update the objective value
		staticVar('objective.value', current.value)
		
		# do NOT stop
		FALSE
}

nmf.stop.consensus <- function(strategy, i, target, data, ...){
			
		# test convergence only every 10 iterations
		if( i %% 10 != 0 ) return( FALSE );
				
		stopconv <- 40		
		
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
