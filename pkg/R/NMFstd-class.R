# Class that implements the standard NMF model
# 
# Author: Renaud Gaujoux \email{renaud@@cbio.uct.ac.za}
###############################################################################

#' Implementation of the standard NMF model 
setClass('NMFstd'
		, representation(
				W = 'matrix', # basis matrix
				H = 'matrix' # mixture coefficients matrix	
		)
		
		, prototype = prototype(
				W = matrix(NA, 0, 0),
				H = matrix(NA, 0, 0)
		)
		
		, validity = function(object){
			
			# dimension compatibility: W and H must be compatible for matrix multiplication
			if( ncol(object@W) != nrow(object@H) ) return(paste('Dimensions of W and H are not compatible [ncol(W)=', ncol(object@W) , '!= nrow(H)=', nrow(object@H), ']'))
			# give a warning if the dimensions look strange: rank greater than the number of samples
			if( !is.empty.nmf(object) && ncol(object@W) > ncol(object@H) ) warning(paste('Dimensions of W and H look strange [ncol(W)=', ncol(object@W) , '> ncol(H)=', ncol(object@H), ']'))
			
			# everything went fine: return TRUE
			return(TRUE)
		}
		, contains = 'NMF'
)


#' Get/Set the basis matrix
setMethod('basis', 'NMFstd',
	function(object){ 
		object@W 
	}
)
setReplaceMethod('basis', signature(object='NMFstd', value='matrix'), 
	function(object, value){ 
		object@W <- value		
		object # TODO: valid object before returning it (+param check=TRUE or FALSE)
	} 
)

#' Get/Set the mixture matrix
setMethod('coef', 'NMFstd',
	function(object){
		object@H
	}
)
setReplaceMethod('coef', signature(object='NMFstd', value='matrix'), 
		function(object, value){ 
			object@H <- value # TODO: valid object before returning it (+param check=TRUE or FALSE)			
			object
		} 
)

#' Compute the estimate target matrix for NMF \emph{standard model}.
#' 
#' The estimate matrix is computed as the porduct of the two matrix slots \code{W} and \code{H}:
#' \deqn{\hat{V} = WH.} 
#' 
#' @param x 
#' @returnType matrix
#' @return 
#' @author Renaud Gaujoux \email{renaud@@cbio.uct.ac.za}
#' @export
setMethod('fitted', signature(object='NMFstd'), 
	function(object, W, H, ...){
		if( missing(W) ) W <- object@W
		if( missing(H) ) H <- object@H
		return(W %*% H)
	}
)

