#' @include NMF-class.R
#' @include fixed-terms.R
NULL

#' TRANSFAC Model - Multidimensional Nonnegative Matrix Factorization
#' 
#' This class provides data structure for NMF array models (e.g., TRANSFAC models), 
#' where the basis and/or the coefficient data is stored in a multidimensional array.
#' 
#' @export
#' @family NMF-model
#' @examples
#'  
#' # create a completely empty NMF object
#' new('NMFarray')
#' 
#' # create a NMF object based on random (compatible) arrays
#' n <- 4; r <- 3; p <- 5; q <- 2
#' w <- array(seq(n*r*q), dim = c(n, r, q))
#' h <- rmatrix(r, p)
#' new('NMFarray', W=w, H=h)
#' 
setClass('NMFarray'
	, representation(
	    W = 'array' # basis array
        , H = 'array' # coef array
	)
	, contains = 'NMF'
  	, prototype=prototype(
        W = matrix(as.numeric(NA), 0, 0),
		H = matrix(as.numeric(NA), 0, 0)
	)
	, validity = function(object){
		
		# dimension compatibility: W and H must be compatible for matrix multiplication
		if( ncol(object@W) != nrow(object@H) ){
			return(paste('Dimensions of W and H are not compatible [ncol(W)=', ncol(object@W) , '!= nrow(H)=', nrow(object@H), ']'))
		}
		# give a warning if the dimensions look strange: rank greater than the number of samples
		if( !is.empty.nmf(object) && ncol(object@H) && nrow(object@W) && ncol(object@W) > ncol(object@H) ){
			warning(paste('Dimensions of W and H look strange [ncol(W)=', ncol(object@W) , '> ncol(H)=', ncol(object@H), ']'))
		}
		
		# everything went fine: return TRUE
		return(TRUE)
	}
)

#' Show method for objects of class \code{NMFarray}
#' @export
setMethod('show', 'NMFarray', 
	function(object)
	{		
		callNextMethod()
        l <- c(dim(basis(object))[3L], dim(coef(object))[3L])
        l[is.na(l)] <- 1L
        cat("levels: ", l[1L], "|", l[2L], "\n")
	}
)

setMethod('ibterms', 'NMFarray', function(object){})
setMethod('icterms', 'NMFarray', function(object){})

#' Get the basis matrix in NMF array models 
#' 
#' This function returns slot \code{W} of \code{object}.
#' 
setMethod('.basis', 'NMFarray',
    function(object, slice = NULL){
        if( is.null(slice) ) object@W
        else object@W[, , slice]
    }
)

#' Replaces a slice of the basis array.
setReplaceMethod('.basis', signature(object='NMFarray', value='matrix'), 
    function(object, ..., slice = 1L, value){
        # error if passed extra arguments
        if( length(xargs<- list(...)) ){
            stop(".basis<-,NMFarray - Unused arguments: ", str_out(xargs, Inf, use.names = TRUE))
        }
        if( length(dim(object@W)) > 2L ) object@W[,, slice] <- value
        else if( slice == 1L ) object@W <- value
        else stop("Invalid slice argument (>1): basis data is a matrix.")
        object
    }
)

#' Set the basis array in NMF array models 
#' 
#' This function sets slot \code{W} of \code{object}.
setReplaceMethod('.basis', signature(object='NMFarray', value='array'), 
    function(object, value){ 
        object@W <- value
        object
    }
)

#' Get the mixture coefficient array in standard NMF array models 
#' 
#' This function returns slot \code{H} of \code{object}.
setMethod('.coef', 'NMFarray',
    function(object, slice = NULL){
        if( is.null(slice) ) object@H
        else object@H[, , slice]
    }
)
#' Set the mixture coefficient matrix in standard NMF models 
#' 
#' This function sets slot \code{H} of \code{object}.
setReplaceMethod('.coef', signature(object='NMFarray', value='array'), 
    function(object, value){
        object@H <- value
        object
    }
)

#' Replaces a slice of the coefficent array.
setReplaceMethod('.coef', signature(object='NMFarray', value='matrix'), 
    function(object, ..., slice = 1L, value){
        # error if passed extra arguments
        if( length(xargs<- list(...)) ){
            stop(".coef<-,NMFarray - Unused arguments: ", str_out(xargs, Inf, use.names = TRUE))
        }
        if( length(dim(object@H)) > 2L ) object@H[,, slice] <- value
        else if( slice == 1L ) object@H <- value
        else stop("Invalid slice argument (>1): coefficient data is a matrix.")
        object
    }
)

#' Computes the target matrix estimate for an NMF array model.
#' @inline
setMethod('fitted', signature(object='NMFarray'), 
	function(object, W, H){
		if( missing(W) ) W <- object@W
		if( missing(H) ) H <- object@H
        
        # handle different dimension cases
        if( is.matrix(object@W) &&  is.matrix(object@H) ) return(object@W %*% object@H)
        else{
            nl <- dim(object)[4L]
            res <- if( is.matrix(object@H) ) apply(object@W, 3L, `%*%`, object@H)
                    else{
                        l <- setNames(seq(nl), dimnames(object@H)[3L])            
                        if( is.matrix(object@W) ) sapply(l, function(k) object@W %*% object@H[,,k])
                        else sapply(l, function(k) object@W[,,k] %*% object@H[,,k])
                    }
            # reshape into an array
            array(res, dim = c(nrow(W), ncol(H), nl))
        }
	}
)

##' Generates a random NMF model with offset, from class \code{NMFOffset}.
##' 
##' The offset values are drawn from a uniform distribution between 0 and 
##' the maximum entry of the basis and coefficient matrices, which are drawn 
##' by the next suitable \code{\link{rnmf}} method, which is the workhorse 
##' method \code{rnmf,NMF,numeric}.
##'   
##' @examples
##' 
##' # random NMF model with offset
##' x <- rnmf(2, 3, model='NMFOffset')
##' x
##' offset(x)
##' # from a matrix
##' x <- rnmf(2, rmatrix(5,3, max=10), model='NMFOffset')
##' offset(x)
##' 
#setMethod('rnmf', signature(x='NMFarray', target='numeric'), 
#function(x, target, ...){	
#	
#	# call the parent's 'rnmf' method to build a standard random NMF factorization
#	res <- callNextMethod()
#		
#	#Vc# Initialize a random offset of length the number of genes
#	res@offset <- runif(nrow(res), min=0, max=max(basis(res), coef(res)));
#	
#	# return the initialized NMFOffset object
#	res
#})
