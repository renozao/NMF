# Standard NMF algorithms
# 
# Author: Renaud Gaujoux
# Creation: 30 Apr 2012
###############################################################################

#' @include NMFstd-class.R
#' @include NMFOffset-class.R
#' @include NMFns-class.R
#' @include registry-algorithms.R
NULL


################################################################################
# BRUNET (standard KL-based NMF)
################################################################################

#' NMF Algorithm/Updates for Kullback-Leibler Divergence
#' 
#' The built-in NMF algorithms described here minimise 
#' the Kullback-Leibler divergence (KL) between an NMF model and a target matrix. 
#' They use the updates for the basis and coefficient matrices (\eqn{W} and \eqn{H}) 
#' defined by \cite{Brunet2004}, which are essentially those from \cite{Lee2001}, 
#' with an stabilisation step that shift up all entries from zero every 10 iterations, 
#' to a very small positive value.
#' 
#' @param i current iteration number.
#' @param v target matrix.
#' @param x current NMF model, as an \code{\linkS4class{NMF}} object.
#' @param eps small numeric value used to ensure numeric stability, by shifting up
#' entries from zero to this fixed value.
#' @param ... extra arguments. These are generally not used and present
#' only to allow other arguments from the main call to be passed to the 
#' initialisation and stopping criterion functions (slots \code{onInit} and 
#' \code{Stop} respectively). 
#' @inheritParams nmf_update.KL.h
#' 
#' @author 
#' Original implementation in MATLAB: Jean-Philippe Brunet \email{brunet@@broad.mit.edu}
#' 
#' Port to R and optimisation in C++: Renaud Gaujoux
#' 
#' @source
#' 
#' Original MATLAB files and references can be found at:
#' 
#' \url{http://www.broadinstitute.org/mpr/publications/projects/NMF/nmf.m}
#' 
#' \url{http://www.broadinstitute.org/publications/broad872}
#' 
#' Original license terms:
#'  
#' This software and its documentation are copyright 2004 by the
#' Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
#' This software is supplied without any warranty or guaranteed support whatsoever. 
#' Neither the Broad Institute nor MIT can not be responsible for its use, misuse, 
#' or functionality. 
#' 
#' @details
#' \code{nmf_update.brunet_R} implements in pure R a single update step, i.e. it updates 
#' both matrices.
#' 
#' @export
#' @rdname KL-nmf
#' @aliases KL-nmf
nmf_update.brunet_R <- function(i, v, x, eps=.Machine$double.eps, ...)
{
	# retrieve each factor
	w <- .basis(x); h <- .coef(x);
	
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
	
	#return the modified model
	.basis(x) <- w; .coef(x) <- h;	
	return(x)
	
}

#' \code{nmf_update.brunet} implements in C++ an optimised version of the single update step.
#'  
#' @export
#' @rdname KL-nmf
nmf_update.brunet <- function(i, v, x, copy=FALSE, eps=.Machine$double.eps, ...)
{
	# retrieve each factor
	w <- .basis(x); h <- .coef(x);
	# fixed terms
	nb <- nbterms(x); nc <- ncterms(x)
	
	# standard divergence-reducing NMF update for H	
	h <- std.divergence.update.h(v, w, h, nbterms=nb, ncterms=nc, copy=copy)
	
	# standard divergence-reducing NMF update for W
	w <- std.divergence.update.w(v, w, h, nbterms=nb, ncterms=nc, copy=copy)
	
	#every 10 iterations: adjust small values to avoid underflow
	# NB: one adjusts in place even when copy=TRUE, as 'h' and 'w' are local variables
	if( i %% 10 == 0 ){
		#eps <- .Machine$double.eps
		h <- pmax.inplace(h, eps, icterms(x))
		w <- pmax.inplace(w, eps, ibterms(x))
	}
	
	# update object if the updates duplicated the model
	if( copy ){		
		#return the modified model	
		.basis(x) <- w; 
		.coef(x) <- h;
	}
	return(x)
	
}

#' Algorithms \sQuote{brunet} and \sQuote{.R#brunet} provide the complete NMF algorithm from \cite{Brunet2004}, 
#' using the C++-optimised and pure R updates \code{\link{nmf_update.brunet}} and \code{\link{nmf_update.brunet_R}} 
#' respectively.
#' 
#' @inheritParams run,NMFStrategyIterative,matrix,NMFfit-method
#' @inheritParams nmf.stop.connectivity
#' 
#' @rdname KL-nmf
#' @aliases brunet_R-nmf
nmfAlgorithm.brunet_R <- setNMFMethod('.R#brunet'
		, objective='KL' 
		, Update=nmf_update.brunet_R
		, Stop='connectivity')

# Optimised version
#' @rdname KL-nmf
#' @aliases brunet-nmf
nmfAlgorithm.brunet <- setNMFMethod('brunet', '.R#brunet', Update=nmf_update.brunet)

#' Algorithm \sQuote{KL} provides an NMF algorithm based on the C++-optimised version of 
#' the updates from \cite{Brunet2004}, which uses the stationarity of the objective value 
#' as a stopping criterion \code{\link{nmf.stop.stationary}}, instead of the  
#' stationarity of the connectivity matrix \code{\link{nmf.stop.connectivity}} as used by 
#' \sQuote{brunet}.
#' 
#' @inheritParams nmf.stop.stationary
#' 
#' @rdname KL-nmf
nmfAlgorithm.KL <- setNMFMethod('KL'
		, objective='KL' 
		, Update=nmf_update.brunet
		, Stop='stationary')

################################################################################
# LEE (standard Euclidean-based NMF)
################################################################################

#' NMF Algorithm/Updates for Frobenius Norm
#' 
#' The built-in NMF algorithms described here minimise 
#' the Frobenius norm (Euclidean distance) between an NMF model and a target matrix. 
#' They use the updates for the basis and coefficient matrices (\eqn{W} and \eqn{H}) 
#' defined by \cite{Lee2001}.
#' 
#' @inheritParams nmf_update.brunet
#' @inheritParams nmf_update.euclidean.h
#' @param rescale logical that indicates if the basis matrix \eqn{W} should be 
#' rescaled so that its columns sum up to one. 
#' 
#' @author 
#' Original update definition: D D Lee and HS Seung
#' 
#' Port to R and optimisation in C++: Renaud Gaujoux
#' 
#' @details
#' \code{nmf_update.lee_R} implements in pure R a single update step, i.e. it updates 
#' both matrices.
#' 
#' @export
#' @rdname Frobenius-nmf
#' @aliases Frobenius-nmf
nmf_update.lee_R <- function(i, v, x, rescale=TRUE, eps=10^-9, ...)
{
	# retrieve each factor
	w <- .basis(x); h <- .coef(x);	
	
	#precision threshold for numerical stability
	#eps <- 10^-9
	
	# compute the estimate WH
	#wh <- estimate(x)
	
	# euclidean-reducing NMF iterations	
	# H_au = H_au (W^T V)_au / (W^T W H)_au
	#h <- pmax(h * (t(w) %*% v),eps) / ((t(w) %*% w) %*% h + eps);
	h <- R_std.euclidean.update.h(v, w, h, eps=eps)
	
	# update H and recompute the estimate WH
	#metaprofiles(x) <- h
	#wh <- estimate(x)
	
	# W_ia = W_ia (V H^T)_ia / (W H H^T)_ia and columns are rescaled after each iteration	
	#w <- pmax(w * (v %*% t(h)), eps) / (w %*% (h %*% t(h)) + eps);
	w <- R_std.euclidean.update.w(v, w, h, eps=eps)
	#rescale columns TODO: effect of rescaling? the rescaling makes the update with offset fail
	if( rescale ) w <- sweep(w, 2L, colSums(w), "/", check.margin=FALSE)
	
	#return the modified model
	.basis(x) <- w; .coef(x) <- h;	
	return(x)
}

#' \code{nmf_update.lee} implements in C++ an optimised version of the single update step.
#'  
#' @export
#' @rdname Frobenius-nmf	
nmf_update.lee <- function(i, v, x, rescale=TRUE, copy=FALSE, eps=10^-9, weight=NULL, ...)
{
	# retrieve each factor
	w <- .basis(x); h <- .coef(x);
	# fixed terms
	nb <- nbterms(x); nc <- ncterms(x)
	
	#precision threshold for numerical stability
	#eps <- 10^-9
	
	# compute the estimate WH
	#wh <- estimate(x)
	
	# euclidean-reducing NMF iterations	
	# H_au = H_au (W^T V)_au / (W^T W H)_au
	h <- std.euclidean.update.h(v, w, h, eps=eps, nbterms=nb, ncterms=nc, copy=copy)
	# update original object if not modified in place
	if( copy ) .coef(x) <- h
	
	# W_ia = W_ia (V H^T)_ia / (W H H^T)_ia and columns are rescaled after each iteration	
	w <- std.euclidean.update.w(v, w, h, eps=eps, weight=weight, nbterms=nb, ncterms=nc, copy=copy)
	#rescale columns TODO: effect of rescaling? the rescaling makes the update with offset fail
	if( rescale ){
		w <- sweep(w, 2L, colSums(w), "/", check.margin=FALSE)
    }
	
	#return the modified model
	.basis(x) <- w; 	
	return(x)
}

#' Algorithms \sQuote{lee} and \sQuote{.R#lee} provide the complete NMF algorithm from \cite{Lee2001}, 
#' using the C++-optimised and pure R updates \code{\link{nmf_update.lee}} and \code{\link{nmf_update.lee_R}}
#' respectively.
#' 
#' @inheritParams run,NMFStrategyIterative,matrix,NMFfit-method
#' @inheritParams nmf.stop.connectivity
#' 
#' @rdname Frobenius-nmf
#' @aliases lee_R-nmf
nmfAlgorithm.lee_R <- setNMFMethod('.R#lee', objective='euclidean'
		, Update=nmf_update.lee_R
		, Stop='connectivity')	

# Optimised version
#' @rdname Frobenius-nmf
#' @aliases lee-nmf
nmfAlgorithm.lee <- setNMFMethod('lee', '.R#lee', Update=nmf_update.lee)

#' Algorithm \sQuote{Frobenius} provides an NMF algorithm based on the C++-optimised version of 
#' the updates from \cite{Lee2001}, which uses the stationarity of the objective value 
#' as a stopping criterion \code{\link{nmf.stop.stationary}}, instead of the  
#' stationarity of the connectivity matrix \code{\link{nmf.stop.connectivity}} as used by 
#' \sQuote{lee}.
#' 
#' @inheritParams nmf.stop.stationary
#' 
#' @rdname Frobenius-nmf
nmfAlgorithm.Frobenius <- setNMFMethod('Frobenius', objective='euclidean'
		, Update=nmf_update.lee
		, Stop='stationary')

################################################################################
# OFFSET (Euclidean-based NMF with offset) [Badea (2008)]
################################################################################


#' NMF Multiplicative Update for NMF with Offset Models
#' 
#' These update rules proposed by \cite{Badea2008} are modified version of 
#' the updates from \cite{Lee2001}, that include an offset/intercept vector, 
#' which models a common baseline for each feature accross all samples: 
#' \deqn{V \approx W H + I}
#' 
#' \code{nmf_update.euclidean_offset.h} and \code{nmf_update.euclidean_offset.w} 
#' compute the updated NMFOffset model, using the optimized \emph{C++} implementations.
#' 
#' @details 
#' The associated model is defined as an \code{\linkS4class{NMFOffset}} object. 
#' The details of the multiplicative updates can be found in \cite{Badea2008}.
#' Note that the updates are the ones defined for a single datasets, not the 
#' simultaneous NMF model, which is fit by algorithm \sQuote{siNMF} from 
#' formula-based NMF models.
#' 
#' @inheritParams nmf_update.brunet
#' @inheritParams nmf_update.euclidean.h
#' 
#' @param offset current value of the offset/intercept vector.
#' It must be of length equal to the number of rows in the target matrix.
#' 
#' @author 
#' Original update definition: Liviu Badea
#' 
#' Port to R and optimisation in C++: Renaud Gaujoux
#' 
#' @return an \code{\linkS4class{NMFOffset}} model object.
#' 
#' @export
#' @rdname offset-nmf
nmf_update.euclidean_offset.h <- function(v, w, h, offset, eps=10^-9, copy=TRUE){
	.Call("offset_euclidean_update_H", v, w, h, offset, eps, copy, PACKAGE='NMF')
}
#' @export 
#' @rdname offset-nmf
nmf_update.euclidean_offset.w <- function(v, w, h, offset, eps=10^-9, copy=TRUE){
	.Call("offset_euclidean_update_W", v, w, h, offset, eps, copy, PACKAGE='NMF')
}
#' \code{nmf_update.offset_R} implements a complete single update step, 
#' using plain R updates.
#' @export 
#' @rdname offset-nmf
nmf_update.offset_R <- function(i, v, x, eps=10^-9, ...)
{	
	# retrieve each factor
	w <- .basis(x); h <- .coef(x);
	# retrieve offset and fill it if necessary (with mean of rows)
	off <- offset(x)
	if( i == 1 && length(off) == 0 )
		off <- rowMeans(v)
	
	#precision threshold for numerical stability
	#eps <- 10^-9
	
	# compute standard lee update (it will take the offset into account) without rescaling W's columns
	
	h <- R_std.euclidean.update.h(v, w, h, wh=w%*%h + off, eps=eps)
	w <- R_std.euclidean.update.w(v, w, h, wh=w%*%h + off, eps=eps)
	#x <- nmf_update.lee(i, v, x, rescale=FALSE, ...)
	
	# update the offset	
	# V0_i = V0_i ( sum_j V_ij ) / ( sum_j (V.off + W H)_ij )
	x@offset <- off * pmax(rowSums(v), eps) / (rowSums(w%*%h + off) + eps)
	
	#return the modified model
	.basis(x) <- w; .coef(x) <- h;
	return(x)
}
#' \code{nmf_update.offset} implements a complete single update step, 
#' using C++-optimised updates.
#' @export 
#' @rdname offset-nmf
nmf_update.offset <- function(i, v, x, copy=FALSE, eps=10^-9, ...)
{	
	# retrieve each factor
	w <- .basis(x); h <- .coef(x);
	# retrieve offset and fill it if necessary (with mean of rows)
	off <- offset(x)
	if( i == 1 && length(off) == 0 )
		off <- rowMeans(v)
	
	#precision threshold for numerical stability
	#eps <- 10^-9
	
	# compute standard offset updates
	h <- nmf_update.euclidean_offset.h(v, w, h, off, eps=eps, copy=copy)
	w <- nmf_update.euclidean_offset.w(v, w, h, off, eps=eps, copy=copy)
	
	# update the offset	
	# V0_i = V0_i ( sum_j V_ij ) / ( sum_j (V.off + W H)_ij )
	x@offset <- off * pmax(rowSums(v), eps) / (rowSums(w%*%h + off) + eps)	
	
	# update the original object if not modified in place
	if( copy ){ 
		.basis(x) <- w; 
		.coef(x) <- h;
	}
	return(x)
}

#' Algorithms \sQuote{offset} and \sQuote{.R#offset} provide the complete NMF-with-offset algorithm 
#' from \cite{Badea2008}, using the C++-optimised and pure R updates \code{\link{nmf_update.offset}} 
#' and \code{\link{nmf_update.offset_R}} respectively.
#' 
#' @inheritParams run,NMFStrategyIterative,matrix,NMFfit-method
#' @inheritParams nmf.stop.connectivity
#' 
#' @rdname offset-nmf
#' @aliases offset_R-nmf
nmfAlgorithm.offset_R <- setNMFMethod('.R#offset', objective='euclidean'
		, model = 'NMFOffset'
		, Update=nmf_update.offset_R
		, Stop='connectivity')

# NMF with offset (optimised version)
#' @rdname offset-nmf
nmfAlgorithm.offset <- setNMFMethod('offset', '.R#offset', Update=nmf_update.offset)

################################################################################
# Non-smooth NMF (KL-based NMF) [Pascual-Montano (2006)]
################################################################################

#' NMF Multiplicative Update for Nonsmooth Nonnegative Matrix Factorization (nsNMF).
#' 
#' These update rules, defined for the \code{\linkS4class{NMFns}} model \eqn{V \approx W S H} from 
#' \cite{Pascual-Montano2006}, that introduces an intermediate smoothing matrix to enhance
#' sparsity of the factors.  
#' 
#' \code{nmf_update.ns} computes the updated nsNMF model.
#' It uses the optimized \emph{C++} implementations \code{\link{nmf_update.KL.w}} and 
#' \code{\link{nmf_update.KL.h}} to update \eqn{W} and \eqn{H} respectively.
#' 
#' @details
#' The multiplicative updates are based on the updates proposed by \cite{Brunet2004}, 
#' except that the NMF estimate \eqn{W H} is replaced by \eqn{W S H} and \eqn{W} 
#' (resp. \eqn{H}) is replaced by \eqn{W S} (resp. \eqn{S H}) in the update of 
#' \eqn{H} (resp. \eqn{W}).
#' 
#' See \code{\link{nmf_update.KL}} for more details on the update formula.
#' 
#' @inheritParams nmf_update.brunet
#' 
#' @return an \code{\linkS4class{NMFns}} model object.
#' 
#' @export
#' @rdname nsNMF-nmf
nmf_update.ns <- function(i, v, x, copy=FALSE, ...)
{
	# retrieve and alter the factors for updating H
	S <- smoothing(x)
	w <- .basis(x)
	h <- .coef(x);
	
	# standard divergence-reducing update for H with modified W
	h <- std.divergence.update.h(v, w %*% S, h, copy=copy)
	
	# update H if not modified in place
	if( copy ) .coef(x) <- h
	
	# standard divergence-reducing update for W with modified H
	w <- std.divergence.update.w(v, w, S %*% h, copy=copy)
	
	# rescale columns of W
	w <- sweep(w, 2L, colSums(w), '/', check.margin=FALSE)
	
	#return the modified model
	.basis(x) <- w;
	return(x)
}
#' \code{nmf_update.ns_R} implements the same updates in \emph{plain R}.
#' 
#' @export
#' @rdname nsNMF-nmf
nmf_update.ns_R <- function(i, v, x, ...)
{
	# retrieve and alter the factors for updating H
	S <- smoothing(x)
	w <- .basis(x)
	#w <- metagenes(x) %*% smoothing(fit(x)); # W <- WS
	h <- .coef(x);
	
	# compute the estimate WH
	#wh <- estimate(x, W=w.init, H=h, S=S)
	
	# standard divergence-reducing update for H with modified W
	h <- R_std.divergence.update.h(v, w %*% S, h)
	
	# update H and recompute the estimate WH
	.coef(x) <- h
	# retrieve and alter the factors for updating W
	#w <- tmp;
	#h <- smoothing(fit(x)) %*% metaprofiles(x); # H <- SH
	#h <- S %*% h; # H <- SH
	
	# standard divergence-reducing update for W with modified H
	w <- R_std.divergence.update.w(v, w, S %*% h)
	
	# rescale columns of W
	w <- sweep(w, 2L, colSums(w), '/', check.margin=FALSE)
	
	#return the modified model
	.basis(x) <- w; #metaprofiles(x) <- h;
	return(x)
}

## REGISTRATION
#' Algorithms \sQuote{nsNMF} and \sQuote{.R#nsNMF} provide the complete NMF algorithm from \cite{Pascual-Montano2006}, 
#' using the C++-optimised and plain R updates \code{\link{nmf_update.brunet}} and \code{\link{nmf_update.brunet_R}} 
#' respectively.
#' The stopping criterion is based on the stationarity of the connectivity matrix.
#' 
#' @inheritParams run,NMFStrategyIterative,matrix,NMFfit-method
#' @inheritParams nmf.stop.connectivity
#' 
#' @rdname nsNMF-nmf
#' @aliases nsNMF_R-nmf
nmfAlgorithm.nsNMF_R <- setNMFMethod('.R#nsNMF', objective='KL'
		, model='NMFns'
		, Update=nmf_update.ns_R
		, Stop='connectivity')

# Optmized version
#' @rdname nsNMF-nmf
nmfAlgorithm.nsNMF <- setNMFMethod('nsNMF', '.R#nsNMF', Update=nmf_update.ns)
