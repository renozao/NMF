#' NMF Algorithm: Pattern Expression NMF
#'
#' Implements the PE-NMF algorithm from Zhang et al (2008).
#'
#' It is implemented using the iterative schema defined by the 
#' NMFStrategyIterative class.
#' The algorithm minimizes the Frobenius norm, with two regularization terms
#' (one for each matrix factor) parametrized by two parameters:
#' 
#' min_{W,H} 1/2 ||V - WH||² 
#' 			+ alpha \sum_{i<>j} W_i^T W_j 
#' 			+ beta \sum_{i,j} H_{ij}
#' 
#' So there is two parameters: alpha and beta.
#' The updates for the matrix factors are (in R notations):
#' 
#' H_{i+1} = H_i ( W_i^T %*% V ) / ( W_i^T %*% W_i %*% H_i + beta)
#' W_{i+1} = W_i ( V %*% H_i^T ) / ( W_i %*% H_i %*% H_i^T + alpha W_i %*% M )
#'
#' with matrix M is full of one with diagonal zero.
#' 
#' @author Renaud Gaujoux
#' @creation 17 Jan 2010
#' 

penmf.objective <- function(x, fit, alpha, beta, ...)
{
	w <- basis(fit)
	1/2 * sum( (x - fitted(fit))^2 )
		+ alpha * ( crossprod(w) - sum(w^2) )
		+ beta * sum(coef(fit))
}

nmf.update.penmf <- function(i, x, data, alpha, beta, ...){
	
	# retrieve each factor
	w <- basis(data); h <- coef(data);
	
	# At the first iteration initialise matrix M
	if( TRUE || i == 1 ){
		r <- ncol(w)
		M <- matrix(1, nrow=r, ncol=r) - diag(1, r)
		#staticVar('M', M, init=TRUE)
	}
	#else M <- staticVar('M')
	
	#precision threshold for numerical stability
	eps <- 10^-9
	
	# H_{i+1} = H_i ( W_i^T %*% V ) / ( W_i^T %*% W_i %*% H_i + beta)
	h <- h * crossprod(w, x) / ( crossprod(w) %*% h + beta)
	
	# W_{i+1} = W_i ( V %*% H_i^T ) / ( W_i %*% H_i %*% H_i^T + alpha W_i %*% M )
	w <- w * tcrossprod(x, h) / ( w %*% tcrossprod(h) + alpha * w %*% M )
	
	#return the modified data
	basis(data) <- w; coef(data) <- h;
	data
}

# Hook to register the algorithm
.nmf.plugin.penmf <- function(){
	
	# PE-NMF
	new('NMFStrategyIterative'
		, name='pe-nmf', objective = penmf.objective
		, model='NMFstd'
		, Update= nmf.update.penmf
		, Stop='nmf.stop.stationary')
}