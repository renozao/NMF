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
				Preprocess = '.functionSlot', # Pre-processing before entering the update loop
				Update = '.functionSlot', # update method	
				Stop = '.functionSlot', # method called just after the update
				WrapNMF = '.functionSlot' # method called just before returning the resulting NMF object
				)	
  , prototype=prototype(
				Preprocess = '',
				Update = '',
  				Stop = '',
				WrapNMF = ''
				)
	, contains = 'NMFStrategy'
#	, validity = function(object){
#		
#		# slots must be either character strings or functions
#		#s.list <- getSlots('NMFStrategyIterative')
#		#s.list <- s.list[s.list=='ANY']
#		s.list <- c('Preprocess', 'Update', 'Stop', 'WrapNMF')
#		names(s.list) <- s.list
#		check <- sapply(names(s.list), function(sname){
#					svalue <- slot(object,sname)
#					if( !is.character(svalue) && !is.function(svalue) )
#						return(paste("slot '"
#									, sname
#									, "' must be either a function or a character string (possibly empty)"
#									, sep=''))
#					return(NA)
#				})
#		
#		err <- which(!is.na(check))
#		if( length(err) > 0 ){
#			check <- check[err]
#			if( length(check) > 1 ) check <- c('', check)  
#			return(paste(check, collapse="\n\t- "))
#		}
#		
#		return(TRUE)
#	}
)


setMethod('show', 'NMFStrategyIterative',
	function(object){
		
		#cat('<object of class: NMFStrategyIterative>')
		callNextMethod()
		cat("<Iterative schema:>\n")
		# go through the slots
		#s.list <- getSlots('NMFStrategyIterative')
		#s.list <- s.list[s.list=='ANY']
		s.list <- c('Preprocess', 'Update', 'Stop', 'WrapNMF')
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
	, representation(
				definition = 'NMFStrategyIterative', 
				workspace = 'environment', # workspace to use persistent variables accross methods
				Preprocess = 'function', # Pre-processing before entering the update loop
				Update = 'function', # update method	
				Stop = 'function', # method called just after the update
				WrapNMF = 'function' # method called just before returning the resulting NMF object
				)
	, prototype=prototype(
				Preprocess = function(data, target) data,
				Update = function(...) stop('Update function is required but not defined'),
  				Stop = function(i, data){ FALSE }, # simply return the a FALSE stop signal
				WrapNMF = function(data){ data } # by default do not do anything special
				)
)


#' Creates a NMFStrategyIterativeX object from a NMFStrategyIterative object.
xifyStrategy <- function(strategy, workspace){	
	
	# first check the strategy's validity
	if( is.character(err <- validObject(strategy, test=TRUE)) ){
		stop("Invalid strategy definition:\n\t- ", err)
	}
	
	# intanciate the NMFStrategyIterativeX, creating the strategy's workspace
	strategyX <- new('NMFStrategyIterativeX', definition=strategy, workspace=workspace)
	
	# build the list of 'function' slots in class NMFStrategyIterativeX
	s.list <- getSlots('NMFStrategyIterativeX')
	s.list <- s.list[s.list=='function']
	s.base <- slotNames(strategy)
	for( sname in names(s.list) ){
		# check the the slot match a slot in class NMFStrategyIterative
		if( !is.element(sname, s.base) ) 
			stop("NMFStrategyIterativeX runtime check: slot '", sname, "' not defined in class NMFStrategyIterative")
		
		# get the content of the slot
		svalue <- slot(strategy,sname)
		# if the slot is valid (i.e. it's a non-empty character string), then process the name into a valid function
		if( is.character(svalue) && nchar(svalue) > 0 ){
			# set the slot with the executable version of the function 
			slot(strategyX,sname) <- getFunction(svalue)			 
		}else if( is.function(svalue) ) slot(strategyX,sname) <- svalue		
		# setup a dedicated evaluation environment: use the strategy's workspace
		environment(slot(strategyX,sname)) <- strategyX@workspace			 
	}
	
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
.load.algorithm.NMFStrategyIterative <- function(){
			
	# Brunet
	nmfRegisterAlgorithm(new('NMFStrategyIterative', name='brunet', objective='KL'
					, Update='nmf.update.brunet'
					, Stop='nmf.stop.consensus'
					, WrapNMF='')
			, overwrite=TRUE)
	
	# Lee	
	nmfRegisterAlgorithm(new('NMFStrategyIterative', name='lee', objective='euclidean'
					, Update='nmf.update.lee'
					, Stop='nmf.stop.consensus'
					, WrapNMF='')
			, overwrite=TRUE)

	# NMF with offset
	nmfRegisterAlgorithm(new('NMFStrategyIterative', name='offset', objective='euclidean'
					, model = 'NMFOffset'
					, Update='nmf.update.offset'
					, Stop='nmf.stop.consensus'
					, WrapNMF='')
			, overwrite=TRUE)

	# nsNMF
	nmfRegisterAlgorithm(new('NMFStrategyIterative', name='nsNMF', objective='KL'
					, model='NMFns'
					, Update='nmf.update.ns'
					, Stop='nmf.stop.consensus'
					, WrapNMF='')
			, overwrite=TRUE)	
}

setMethod('run', signature(object='NMFStrategyIterative', target='matrix', start='NMFfit'),
	function(object, target, start, ...){
	
	# debug object in debug mode
	if( nmf.getOption('debug') ) show(object)		
	
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
	strategyX = xifyStrategy(object, .Workspace)
	# call the xified startegy's run method			
	run(strategyX, target, start, ...)
})

#' Generic algorithm for NMF, based on NMFStrategyIterativeX object.
setMethod('run', signature(object='NMFStrategyIterativeX', target='matrix', start='NMFfit'),
	function(object, target, start, maxIter=2000, ...){
				
	strategy <- object
	v <- target
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
	nmfData <- strategy@Preprocess(start, v, ...)
	if( !inherits(nmfData, 'NMFfit') ) stop('NMFStrategyIterative[', name(strategy@definition), ']::Preprocess did not return a "NMF" instance [returned: "', class(nmfData), '"]')	
	#message('NMFStrategyIterative:: object class:', class(nmfData))
	
	#Vc# Start iterations
	nmfFit <- fit(nmfData)
	for( i in 1:maxIter ){
		
		#Vc# update the matrices
		nmfFit <- strategy@Update(i, v, nmfFit)
		
		#Vc# Stopping criteria
		# give the strategy the opportunity to perform stuff after the update: modify the data and/or stop iteration
		stop.signal <- strategy@Stop(i, v, nmfFit)
		
		# every now and then track the error if required
		nmfData <- trackError(nmfData, objective(strategy@definition, v, nmfFit), i)
		
		# if the strategy ask for stopping, then stop the iteration
		if( stop.signal ) break;
				
	}
	fit(nmfData) <- nmfFit
	if( verbose(nmfData) ) cat("\n")	
	
	#Vc# wrap up
	# let the strategy build the result
	nmfData = strategy@WrapNMF(nmfData)	
	if( !inherits(nmfData, 'NMFfit') ) stop('NMFStrategyIterative[', name(strategy@definition), ']::WrapNMF did not return a "NMF" instance [returned: "', class(nmfData), '"]')
	
	# force to compute last error if not already done	
	nmfData <- trackError(nmfData, objective(strategy@definition, v, nmfData), i, force=TRUE)
	
	#return the result
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
std.divergence.update.h <- function(v, w, h, wh=NULL)
{	
	# compute WH if necessary	
	if( is.null(wh) ) wh <- w %*% h
	
	# divergence-reducing NMF iterations
	# H_au = H_au ( sum_i [ W_ia V_iu / (WH)_iu ] ) / ( sum_k W_ka ) -> each row of H is divided by a the corresponding colSum of W
	x1 <- colSums(w); # division will recycle the elements (=> for each column c of H we'll have: c <- c/x1 )
	h * crossprod(w, v / wh) / x1;
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
std.divergence.update.w <- function(v, w, h, wh=NULL)
{			
	# compute WH if necessary	
	if( is.null(wh) ) wh <- w %*% h
	
	# W_ia = W_ia ( sum_u [H_au A_iu / (WH)_iu ] ) / ( sum_v H_av ) -> each column of W is divided by a the corresponding rowSum of H
	x2 <- matrix(rep(rowSums(h), nrow(w)), ncol=ncol(w), byrow=TRUE); 
	w * tcrossprod(v / wh, h) / x2;
	
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
nmf.update.brunet <- function(i, v, data, ...)
{
	# retrieve each factor
	w <- basis(data); h <- coef(data);
	
	# standard divergence-reducing NMF update for H
	h <- std.divergence.update.h(v, w, h, w %*% h)
	
	# standard divergence-reducing NMF update for W
	w <- std.divergence.update.w(v, w, h, w %*% h)
	
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


std.euclidean.update.h <- function(v, w, h, wh=NULL, eps){
	# compute WH if necessary	
	den <- if( is.null(wh) ) crossprod(w) %*% h
			else{ t(w) %*% wh}
	
	# H_au = H_au (W^T V)_au / (W^T W H)_au
	pmax(h * crossprod(w,v),eps) / (den + eps);
}

std.euclidean.update.w <- function(v, w, h, wh=NULL, eps){
	# compute WH if necessary	
	den <- if( is.null(wh) ) w %*% tcrossprod(h)
			else{ wh %*% t(h)}
	
	# W_ia = W_ia (V H^T)_ia / (W H H^T)_ia and columns are rescaled after each iteration	
	pmax(w * tcrossprod(v, h), eps) / (den + eps);
}
#' Multiplicative update for reducing the euclidean distance.
#'
#' 
nmf.update.lee <- function(i, v, data, rescale=TRUE)
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
	h <- std.euclidean.update.h(v, w, h, eps=eps)
	
	# update H and recompute the estimate WH
	#metaprofiles(data) <- h
	#wh <- estimate(data)

	# W_ia = W_ia (V H^T)_ia / (W H H^T)_ia and columns are rescaled after each iteration	
	#w <- pmax(w * (v %*% t(h)), eps) / (w %*% (h %*% t(h)) + eps);
	w <- std.euclidean.update.w(v, w, h, eps=eps)
	#rescale columns TODO: effect of rescaling? the rescaling makes the update with offset fail
	if( rescale ) w <- apply(w, 2, function(x) x/sum(x))
	
	#return the modified data
	basis(data) <- w; coef(data) <- h;	
	return(data)
}

#' Multiplicative update for reducing the euclidean distance including on offset.
#'
#' The method is a modified version of Lee's method that also fits an offset vector which model a common expression baseline for each gene accross all samples.
nmf.update.offset <- function(i, v, data, ...)
{	
	# retrieve each factor
	w <- basis(data); h <- coef(data);
	
	#precision threshold for numerical stability
	eps <- 10^-9
	
	# compute standard lee update (it will take the offset into account) without rescaling W's columns
	
	h <- std.euclidean.update.h(v, w, h, wh=w%*%h + offset(data), eps=eps)
	w <- std.euclidean.update.w(v, w, h, wh=w%*%h + offset(data), eps=eps)
	#data <- nmf.update.lee(i, v, data, rescale=FALSE, ...)
	
	# update the offset	
	# V0_i = V0_i ( sum_j V_ij ) / ( sum_j (V.off + W H)_ij )
	data@offset <- data@offset * pmax(rowSums(v), eps) / (rowSums(w%*%h + offset(data)) + eps)
		
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
	tmp <- w %*% S
	h <- std.divergence.update.h(v, tmp, h, wh= tmp %*% h)
	
	# update H and recompute the estimate WH
	coef(data) <- h
	# retrieve and alter the factors for updating W
	#w <- tmp;
	#h <- smoothing(fit(data)) %*% metaprofiles(data); # H <- SH
	#h <- S %*% h; # H <- SH
	
	# standard divergence-reducing update for W with modified H
	tmp <- S %*% h
	w <- std.divergence.update.w(v, w, tmp, wh=w %*% tmp)
	
	# rescale columns of W
	w <- apply(w, 2, function(x) x/sum(x))		
		
	#return the modified data
	basis(data) <- w; #metaprofiles(data) <- h;
	return(data)
}


################################################################################################
# AFTER-UPDATE METHODS
################################################################################################

nmf.stop.stationnary <- function(i, target, data){
		
		# first call compute the initial error
		if( i == 1 ) staticVar('objective.value', distance(target, data), init=TRUE)
		
		# test convergence only every 10 iterations
		interval <- 10		
		if( i %% interval != 0 ) return( FALSE );
		
		# initialize static variables to store the error across the calls		
		last.value <- staticVar('objective.value')
		current.value <- distance(target, data)
		
		# if the relative decrease in the objective value is to small then stop
		threshold <- if( !is.null(data@parameters$threshold) ) data@parameters$threshold else 10^-6
		if( abs( (last.value - current.value)/interval ) <= threshold ) return( TRUE )
		
		# update the objective value
		staticVar('objective.value', current.value)
		
		# do NOT stop
		FALSE
}

nmf.stop.consensus <- function(i, target, data){
			
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
