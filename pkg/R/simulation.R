# Functions to simulate NMF data
# 
# Author: Renaud Gaujoux 
###############################################################################

#' @include utils.R
NULL

#' Simulating Datasets
#' 
#' The function \code{syntheticNMF} generates random target matrices that follow
#' some defined NMF model, and may be used to test NMF algorithms.
#' It is designed to designed to produce data with known or clear classes of 
#' samples.
#' 
#' @param n number of rows of the target matrix. 
#' @param r specification of the factorization rank. 
#' It may be a single \code{numeric}, in which case argument \code{p} is required
#' and \code{r} groups of samples are generated from a draw from a multinomial 
#' distribution with equal probabilities, that provides their sizes.
#' 
#' It may also be a numerical vector, which contains the number of samples in 
#' each class (i.e integers). In this case argument \code{p} is discarded
#' and forced to be the sum of \code{r}.
#' @param p number of columns of the synthetic target matrix. 
#' Not used if parameter \code{r} is a vector (see description of argument \code{r}).
#' @param offset specification of a common offset to be added to the synthetic target
#' matrix, before noisification.
#' Its may be a numeric vector of length \code{n}, or a single numeric value that
#' is used as the standard deviation of a centred normal distribution from which 
#' the actual offset values are drawn.
#' @param noise a logical that indicate if noise should be added to the 
#' matrix.
#' @param factors a logical that indicates if the NMF factors should be return 
#' together with the matrix.
#' @param seed a single numeric value used to seed the random number generator 
#' before generating the matrix.
#' The state of the RNG is restored on exit.
#' 
#' @return a matrix, or a list if argument \code{factors=TRUE}.
#' 
#' When \code{factors=FALSE}, the result is a matrix object, with the following attributes set:
#' \describe{
#' \item{coefficients}{the true underlying coefficient matrix (i.e. \code{H});}
#' \item{basis}{the true underlying coefficient matrix (i.e. \code{H});}
#' \item{offset}{the offset if any;}
#' \item{pData}{a \code{list} with one element \code{'Group'} that contains a factor 
#' that indicates the true groups of samples, i.e. the most contributing basis component for each sample;} 
#' \item{fData}{a \code{list} with one element \code{'Group'} that contains a factor 
#' that indicates the true groups of features, i.e. the basis component 
#' to which each feature contributes the most.}
#' } 
#' 
#' Moreover, the result object is an \code{\link{ExposeAttribute}} object, which means that 
#' relevant attributes are accessible via \code{$}, e.g., \code{res$coefficients}.
#' In particular, methods \code{\link{coef}} and \code{\link{basis}} will work as expected
#' and return the true underlying coefficient and basis matrices respectively.
#' 
#' @export
#' @examples
#' 
#' # generate a synthetic dataset with known classes: 50 features, 18 samples (5+5+8)
#' n <- 50
#' counts <- c(5, 5, 8)
#' 
#' # no noise
#' V <- syntheticNMF(n, counts, noise=FALSE)
#' \dontrun{aheatmap(V)}
#' 
#' # with noise
#' V <- syntheticNMF(n, counts)
#' \dontrun{aheatmap(V)}
#' 
syntheticNMF <- function(n, r, p, offset=NULL, noise=TRUE, factors=FALSE, seed=NULL){
	
	# set seed if necessary
	if( !is.null(seed) ){
		os <- RNGseed()
		on.exit( RNGseed(os) )
		set.seed(seed)
	}
	
	# internal parameters
	mu.W <- 1; sd.W <- 1
	if( isTRUE(noise) ){
		noise <- list(mean=0, sd=1)
	}else if( isNumber(noise) ){
		noise <- list(mean=0, sd=noise)
	}else if( is.list(noise) ){
		stopifnot( length(noise) == 2L )
		noise <- setNames(noise, c('mean', 'sd'))
	}else
		noise <- FALSE
	
	if( length(r) == 1 ){
		g <- rmultinom(1, p, rep(1, r))			
	}else{ # elements of r are the number of samples in each class 
		g <- r		
		p <- sum(r) # total number of samples
		r <- length(r) # number of class
	}
	
	# generate H
	H <- matrix(0, r, p)
	tmp <- 0
	for( i in 1:r ){
		H[i,(tmp+1):(tmp+g[i])] <- 1
		tmp <- tmp+g[i]
	} 	
	
	if( length(n) == 1 ){
		b <- rmultinom(1, n, rep(1, r))		
	}else{ # elements of n are the number of genes in each class 
		b <- n
		n <- sum(n)
	}
	
	# generate W
	W <- matrix(0, n, r)
	tmp <- 0
	for( i in 1:r ){		
		W[(tmp+1):(tmp+b[i]),i] <- abs(rnorm(b[i], mu.W, sd.W))
		tmp <- tmp + b[i]
	}	
	
	# build the composite matrix
	res <- W %*% H
	# add the offset if necessary
	if( !is.null(offset) ){
		if( length(offset) == 1L )
			offset <- rnorm(n, mean=0, sd=offset)
		
		stopifnot(length(offset)==n)
		res <- res +  offset
	}
	
	# add some noise if required
	if( !isFALSE(noise) )
		res <- pmax(res + rmatrix(res, dist=rnorm, mean=noise$mean, sd=noise$sd), 0)	
	
	# return the factors if required
	pData <- list(Group=factor(unlist(mapply(rep, 1:r, g, SIMPLIFY=FALSE))))
	fData <- list(Group=factor(unlist(mapply(rep, 1:r, b, SIMPLIFY=FALSE))))
	if( factors ) res <- list(res, W=W, H=H, offset=offset, pData=pData, fData=fData)
	
	# wrap results and expose relevant attributes
	ExposeAttribute(res, coefficients=H, basis=W, offset=offset
						, pData = pData, fData = fData
						, .VALUE=TRUE, .MODE='r')
	
}

