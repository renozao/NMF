#' Algorithm for Nonnegative Matrix Factorization: Local NMF (LNMF)
#'
#' @author Renaud Gaujoux
#' @created 21 Jul 2009


#' Algorithm for Nonnegative Matrix Factorization: Local NMF (LNMF).
#'
#' The local NMF algorithm is minimizes use the following Kullback-Leibler divergence based objective function:
#' $$ 
#' \sum_{i=1}^m\sum_{j=1}^n\left(X_{ij} \log\frac{X_{ij}}{(WH)_{ij}} - X_{ij} + (WH)_{ij} + \alpha U_{ij}\right) - \beta \sum_i V_{ij},
#' $$
#' where $\alpha, \beta > 0$ are some constants, $U = W^TW$ and $V = HH^T$.
#'
#' TODO: add explaination for each terms (see Wild 2002)
#'
#' @references Learning spatially localized, parts-based representation
#' , S.Z. Li, X.W. Hou, and H.J. Zhang.
#' , In Proceedings of IEEE International Conference on Computer Vision and Pattern Recognition
#' , December 2001
R_nmf.update.lnmf <- function(i, v, data, ...){
	
	# retrieve each factor
	w <- basis(data); h <- coef(data);
	
	# update H 
	h <- sqrt( h * crossprod(w, v / (w %*% h)) )
	
	# update W using the standard divergence based update
	w <- R_std.divergence.update.w(v, w, h, w %*% h)
	
	# scale columns of W
	w <- sweep(w, 2L, colSums(w), "/", check.margin=FALSE)	
	
	#every 10 iterations: adjust small values to avoid underflow 
	if( i %% 10 == 0 ){
		#precision threshold for numerical stability
		eps <- .Machine$double.eps
		h[h<eps] <- eps;
		w[w<eps] <- eps;
	}
		
	# return updated data	
	basis(data) <- w; coef(data) <- h
	return(data)
}

nmf.update.lnmf <- function(i, v, data, ...){
	
	# retrieve each factor
	w <- basis(data); h <- coef(data);
	
	# update H 
	h <- sqrt( h * crossprod(w, v / (w %*% h)) )
	
	# update W using the standard divergence based update
	w <- std.divergence.update.w(v, w, h)
	
	# scale columns of W
	w <- apply(w, 2, function(x) x/sum(x))
	
	#every 10 iterations: adjust small values to avoid underflow 
	if( i %% 10 == 0 ){
		#precision threshold for numerical stability
		eps <- .Machine$double.eps
		h[h<eps] <- eps;
		w[w<eps] <- eps;
	}
	
	# return updated data	
	basis(data) <- w; coef(data) <- h
	return(data)
}

# Hook to register the algorithm when the package is loaded
.nmf.plugin.lnmf <- function(){
	return(NULL)
	list(
			new('NMFStrategyIterative', name='lnmf', objective='KL'
						, Update=nmf.update.lnmf
						, Stop='nmf.stop.consensus')
			
			, new('NMFStrategyIterative', name='.R#lnmf', objective='KL'
				, Update=R_nmf.update.lnmf
				, Stop='nmf.stop.consensus')
	)
}
