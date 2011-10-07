###% Class for \emph{Nonsmooth Nonnegative Matrix Factorization} (nsNMF).
###%
###% The Nonsmooth NMF algorithm is a modification of the standard divergence based NMF algorithm (cf. TODO: NMFstandard). 
###% Given a non-negative \eqn{n \times m} matrix \eqn{V} and a factorization rank \eqn{r}, it fits the following model:
###% \deqn{V \equiv WS(\theta)H,}
###% where:
###% \itemize{
###% \item{\eqn{W} and \eqn{H} are such as in the standard model, that is non-negative matrices of dimension 
###% \eqn{n \times r} and \eqn{r \times m} respectively;}
###% \item{\eqn{S} is a \eqn{r \times r} square matrix whose entries depends on an extra parameter \eqn{0\leq\theta\leq 1} 
###% in the following way: \deqn{S = (1-\theta)I + \frac{\theta}{r} 11^T },} where \eqn{I} is the identity matrix and \eqn{1}
###% is a vector of ones.
###% }
###% 
###% The interpretation of S as a smoothing matrix can be explained as follows: Let \eqn{X} be a positive, nonzero, vector.
###% Consider the transformed vector \eqn{Y = SX}. If \eqn{\theta = 0}, then \eqn{Y = X} and no smoothing on \eqn{X} has occurred. 
###% However, as \eqn{\theta \rightarrow 1}, the vector \eqn{Y} tends to the constant vector with all elements almost equal 
###% to the average of the elements of \eqn{X}. This is the smoothest possible vector in the sense of  nonsparseness 
###% because all entries are equal to the same nonzero value, instead of having some values close to zero and 
###% others clearly nonzero.
###% 
###% @include NMF-class.R
###% @references 
###% Alberto Pascual-Montano et al. (2006), 'Nonsmooth Nonnegative Matrix Factorization (nsNMF)'
###% , IEEE Transactions On Pattern Analysis And Machine Intelligence, Vol. 28, No. 3, March 2006 403
###% 
###% @slot theta the smoothing parameter (a single \code{numeric}), that controls the sparseness of the Nonsmooth NMF model. 
###% @author Renaud Gaujoux \email{renaud@@cbio.uct.ac.za}
###% @export
setClass('NMFns'
	, representation(
				theta = 'numeric' # smoothing matrix
				)
	, contains = 'NMFstd'
  	, prototype = prototype(
				theta = 0.5
				)
	, validity = function(object){
		if( object@theta < 0 || object@theta > 1 ) 
			return(paste("Invalid value for theta (",object@theta,"): must be between 0 and 1", sep=''))
		TRUE
	}
	
)

###% Show method for objects of class \code{NMFns}
setMethod('show', 'NMFns', 
		function(object)
		{			
			callNextMethod()
			cat("theta:", object@theta, "\n")
		}
)

###% Compute estimate for an NMFns object, according to the Nonsmooth NMF model (cf. \code{\link{NMFns-class}}).
###% 
###% @param x a \code{NMFns} object
###% @param ... extra parameter passed to method \code{smoothing}. Typically used to pass parameter \code{theta} to
###% compute the smoothing matrix with a value of \code{theta} different from the one stored in \code{x}.
###% 
###% @return the estimate matrix from \code{NMFns} object \code{x}
###% @seealso smoothing, NMFns-class
###% @author Renaud Gaujoux \email{renaud@@cbio.uct.ac.za} 
###% @export
setMethod('fitted', signature(object='NMFns'), 
	function(object, W, H, S, ...){
		if( missing(W) ) W <- object@W
		if( missing(H) ) H <- object@H
		if( missing(S) ) S <- smoothing(object, ...)
		W %*% (S %*% H)		
	}
)

###% Returns the smoothing matrix in the Nonsmooth NMF model.
###% 
###% For a \eqn{r}-rank NMF, the smoothing matrix of parameter \eqn{\theta} is computed as follows:
###% \deqn{S = (1-\theta)I + \frac{\theta}{r} 11^T },} where \eqn{I} is the identity matrix and \eqn{1}
###% is a vector of ones (cf \code{\link{NMFns-class}} for more details).
###% 
###% @param x a \code{NMFns} object
###% @param theta the smoothing parameter (numeric) between 0 and 1.   
###% 
###% @return if \code{x} estimates a \eqn{r}-rank NMF, then the result is a \eqn{r \times r} square matrix.  
###% @author Renaud Gaujoux \email{renaud@@cbio.uct.ac.za}
###% @export
setGeneric('smoothing', function(x, ...) standardGeneric('smoothing') )
setMethod('smoothing', signature(x='NMFns'), 
	function(x, theta=x@theta){	
		# check validity of theta
		if( theta < 0 || theta > 1 ) stop("Invalid smoothing parameter theta('",theta,"'): theta must be susch that 0 <= theta <=1")
		diag(1-theta, nbasis(x)) + theta / nbasis(x)		
	}
)

