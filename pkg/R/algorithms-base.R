# Standard NMF algorithms
# 
# Author: Renaud Gaujoux
# Creation: 30 Apr 2012
###############################################################################

#' @include NMFStrategyIterative-class.R
#' @include NMFstd-class.R
#' @include NMFOffset-class.R
#' @include NMFns-class.R
NULL


################################################################################
# BRUNET (standard KL-based NMF)
################################################################################

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

# Brunet (R version)
setNMFMethod('.R#brunet'
			, objective='KL' 
			, Update=R_nmf.update.brunet
			, Stop='connectivity')

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

# Brunet (optimised version)
setNMFMethod('brunet', '.R#brunet', Update=nmf.update.brunet)

################################################################################
# LEE (standard Euclidean-based NMF)
################################################################################

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

# Lee (R version)	
setNMFMethod('.R#lee', objective='euclidean'
			, Update=R_nmf.update.lee
			, Stop='connectivity')	

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

# Lee (optimised version)	
setNMFMethod('lee', '.R#lee', Update=nmf.update.lee)

################################################################################
# OFFSET (Euclidean-based NMF with offset) [Badea (2008)]
################################################################################


# Updates for NMF with offset
# H
offset.std.euclidean.update.h <- function(v, w, h, offset, eps=10^-9, copy=TRUE){
	.Call("offset_euclidean_update_H", v, w, h, offset, eps, copy, PACKAGE='NMF')
}
# W
offset.std.euclidean.update.w <- function(v, w, h, offset, eps=10^-9, copy=TRUE){
	.Call("offset_euclidean_update_W", v, w, h, offset, eps, copy, PACKAGE='NMF')
}

###% Multiplicative update for reducing the euclidean distance including on offset.
###%
###% The method is a modified version of Lee's method that also fits an offset 
###% vector which model a common expression baseline for each gene accross all 
###% samples.
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

# NMF with offset (R version)
setNMFMethod('.R#offset', objective='euclidean'
		, model = 'NMFOffset'
		, Update=R_nmf.update.offset
		, Stop='connectivity')

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

# NMF with offset (optimised version)
setNMFMethod('offset', '.R#offset', Update=R_nmf.update.offset)

################################################################################
# Non-smooth NMF (KL-based NMF) [Pascual-Montano (2006)]
################################################################################

###% Multiplicative update for Nonsmooth Nonnegative Matrix Factorization (nsNMF).
###%
###% The update rules are essentialy the same as in Brunet, but WH is replaced by WSH 
###% whereas W (resp. H) is replaced by WS (resp. SH) in the update of H (resp. of W).
###%
###% @references 
###% Alberto Pascual-Montano et al. (2006), 'Nonsmooth Nonnegative Matrix Factorization (nsNMF)'
###% , IEEE Transactions On Pattern Analysis And Machine Intelligence, Vol. 28, No. 3, March 2006 403
###%
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


# nsNMF (R version)
setNMFMethod('.R#nsNMF', objective='KL'
		, model='NMFns'
		, Update=R_nmf.update.ns
		, Stop='connectivity')

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


# nsNMF (optimised version)
setNMFMethod('nsNMF', '.R#nsNMF', Update=R_nmf.update.ns)
