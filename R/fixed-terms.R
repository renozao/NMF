# Interface for NMF models that contain fixed terms 
# 
# Author: Renaud Gaujoux
# Creation: 03 Jul 2012
###############################################################################

#' @include NMF-class.R
#' @include NMFstd-class.R
NULL

#' Concatenating NMF Models
#' 
#' Binds compatible matrices and NMF models together.
#' 
#' @param x an NMF model
#' @param ... other objects to concatenate. Currently only two objects at a time 
#' can be concatenated (i.e. \code{x} and \code{..1}).
#' @param margin integer that indicates the margin along which to concatenate 
#' (only used when \code{..1} is a matrix):
#' \describe{
#' \item{1L}{}
#' \item{2L}{}
#' \item{3L}{}
#' \item{4L}{}
#' }
#' If missing the margin is heuristically determined by looking at common 
#' dimensions between the objects. 
#' 
#' @keywords internal
setMethod('c', 'NMF',
	function(x, ..., margin=3L, recursive=FALSE){
		
		y <- ..1
		if( is.matrix(y) ){
			
			if( missing(margin) ){
				if( nrow(y) == nrow(x) ){
					if( ncol(y) == ncol(x) ){
						warning("NMF::`c` - Right argument match both target dimensions: concatenating basis columns."
								, " Use `margin=4L` to concatenate coefficient rows.")
					}
					margin <- 3L
				}else if( ncol(y) == ncol(x) ){
					margin <- 4L
				}else{
					stop("NMF::`c` - Incompatible argument dimensions: could not infer concatenation margin.")
				}
			}
			
			if( margin == 1L ){ # extend basis vectors
				
				if( nbterms(x) ){ # cannot extend models with fixed basis terms
					stop("NMF::`c` - Could not extend basis vectors:"
							, " NMF model has fixed basis terms [", nbterms(x), "]")
				}
				if( ncol(y) != nbasis(x) ){
					stop("NMF::`c` - Could not extend basis vectors:"
							, " incompatible number of columns [", nbasis(x), '!=', ncol(y), "].")
				}
				
				# extend basis vectors
				basis(x) <- rbind(basis(x), y)
				
			} else if( margin == 2L ){ # extend basis profiles
				
				if( ncterms(x) ){ # cannot extend models with fixed coef terms
					stop("NMF::`c` - Could not extend basis profiles:"
							, " NMF model has fixed coefficient terms [", ncterms(x), "]")
				}
				if( nrow(y) != nbasis(x) ){
					stop("NMF::`c` - Could not extend basis profiles:"
							, " incompatible number of rows [", nbasis(x), '!=', nrow(y), "].")
				}
				
				# extend basis profiles
				coef(x) <- cbind(coef(x), y)
				
			} else if( margin == 3L ){ # add basis vectors
				
				if( nrow(y) != nrow(x) ){
					stop("NMF::`c` - Could not concatenate basis vectors:"
						, " incompatible number of rows [", nrow(x), '!=', nrow(y), "].")
				}
				
				# bind basis terms
				.basis(x) <- cbind(basis(x), y)
				dn <- colnames(.basis(x))
				# bind dummy coef
				.coef(x) <- rbind(coef(x), matrix(NA, ncol(y), ncol(x)))
				basisnames(x) <- dn
				
			} else if( margin == 4L ){ # add basis profiles
				
				if( ncol(y) != ncol(x) ){
					stop("NMF::`c` - Could not concatenate basis profiles:"
						, " incompatible number of columns [", ncol(x), '!=', ncol(y), "].")
				}
				
				# bind coef terms
				.coef(x) <- rbind(coef(x), y)
				dn <- rownames(.coef(x))
				# bind dummy basis
				.basis(x) <- cbind(basis(x), matrix(NA, nrow(x), nrow(y)))
				basisnames(x) <- dn
				
			}else{
				stop("NMF::`c` - Invalid concatenation margin: should be either"
								, " 1L (basis rows), 2L (coef columns), 3L (basis vectors/columns) or 4L (basis profiles/coef rows).")
			}
		}else if( is.nmf(y) ){
			# check dimensions
			if( nrow(x) != nrow(y) )
				stop("NMF::`c` - Could not concatenate NMF objects:"
					, " incompatible number of rows [", nrow(x), '!=', nrow(y), "]")
			if( ncol(x) != ncol(y) )
				stop("NMF::`c` - Could not concatenate NMF objects:"
					, " incompatible number of columns [", ncol(x), '!=', ncol(y), "]")
			.basis(x) <- cbind(basis(x), basis(y))
			.coef(x) <- rbind(coef(x), coef(y))
		}else{
			stop("NMF::`c` - Concatenation of an NMF object with objects of class '", class(y), "' is not supported.")
		}
		
		# return augmented object
		x
	}
)


fterms <- function(value){

	res <- list(n=0L, terms=NULL, df=NULL, i=integer())
	if( !length(value) ) return(res)
	
	# convert into a data.frame
	if( is.factor(value) ) value <- data.frame(Group=value)
	else if( is.numeric(value) ) value <- data.frame(Var=value)
	else if( !is.data.frame(value) ) value <- as.data.frame(value)
	
	res$n <- length(value)
    res$df <- value
	# generate fixed term matrix
	terms <- model.matrix(~ -1 + ., data=value)
	res$terms <- terms
	# build indexes
	res$i <- 1:ncol(terms)
		
	res
}

## #' Annotations in NMF Models
## #' 
## #' NMF models may contain annotations for columns/rows and/or rows/features, in 
## #' a similar way gene expression data are annotated 
## #' \code{\linkS4class{ExpressionSet}} objects in Bioconductor.
## #'
## NULL

#' Fixed Terms in NMF Models
#' 
#' These functions are for internal use and should not be called by the end-user.
#' 
#' They use \code{\link{model.matrix}(~ -1 + ., data=value)} to generate suitable
#' term matrices.
#' 
#' @param object NMF object to be updated.
#' @param value specification of the replacement value for fixed-terms.
#' 
#' @rdname terms-internal
#' @keywords internal
#' @inline
setGeneric('bterms<-', function(object, value) standardGeneric('bterms<-'))
#' Default method tries to coerce \code{value} into a \code{data.frame} with 
#' \code{\link{as.data.frame}}.
setReplaceMethod('bterms', signature('NMFstd', 'ANY'),
	function(object, value){
		
		if( nterms(object) ){
			stop("Cannot set fixed basis terms on an object that already has fixed terms:",
					" these can be set only once and before setting any fixed coefficient term",
					" [coef=", ncterms(object), ", basis=", nbterms(object), "].")
		}
		# build terms
		t <- fterms(value)
		
		if( !t$n ) return(object)
		
		# check dimension
		if( nrow(t$terms) != nrow(object) ){
			stop("Invalid fixed basis terms: all terms should have length the number of target rows"
					, "[terms=", nrow(t$terms), " != ", nrow(object), "=target]")
		}
		
		# set data
		object@bterms <- t$df
		# set indexes
		i <- t$i
		nv <- nbasis(object)
		object@ibterms <- nv + i
		# set terms
		object <- c(object, t$terms, margin=3L)
		
		object
	}
)
#' \code{cterms<-} sets fixed coefficient terms or indexes and should only be 
#' called on a newly created NMF object, i.e. in the constructor/factory generic 
#' \code{\link{nmfModel}}.
#' 
#' @rdname terms-internal
#' @inline
setGeneric('cterms<-', function(object, value) standardGeneric('cterms<-'))
#' Default method tries to coerce \code{value} into a \code{data.frame} with 
#' \code{\link{as.data.frame}}.
setReplaceMethod('cterms', signature('NMFstd', 'ANY'),
	function(object, value){
		
		if( ncterms(object) ){
			stop("Cannot set fixed coef terms on an object that already has fixed coef terms:",
				" these can be set only once", 
				" [coef=", ncterms(object), ", basis=", nbterms(object), "].")
		}
		# build terms
		t <- fterms(value)
		
		if( !t$n ) return(object)
		
		# check dimension
		if( nrow(t$terms) != ncol(object) ){
			stop("Invalid fixed coefficient terms: all terms should have length the number of target columns"
					, "[terms=", nrow(t$terms), " != ", ncol(object), "=target]")
		}
		
		# transpose term matrix
		t$terms <- t(t$terms)
		# set data
		object@cterms <- t$df
		# set indexes
		i <- t$i
		nv <- nbasis(object)
		object@icterms <- nv + i
		# set terms
		object <- c(object, t$terms, margin=4L)
					
		object
	}
)

#' Fixed Terms in NMF Models
#' 
#' @description
#' Formula-based NMF models may contain fixed basis and/or coefficient terms.
#' The functions documented here provide access to these data, which are 
#' read-only and defined when the model object is instantiated 
#' (e.g., see \code{\link[=nmfModel,formula,ANY-method]{nmfModel,formula-method}}).
#' 
#' \code{ibterms}, \code{icterms} and \code{iterms} respectively return the 
#' indexes of the fixed basis terms, the fixed coefficient terms and all fixed 
#' terms, within the basis and/or coefficient matrix of an NMF model.
#' 
#' @param object NMF object
#' @param ... extra parameters to allow extension (currently not used)
#' 
#' @export
#' @rdname terms
setGeneric('ibterms', function(object, ...) standardGeneric('ibterms') )
#' Default pure virtual method that ensure a method is defined for concrete 
#' NMF model classes.
setMethod('ibterms', 'NMF', 
	function(object, ...){
		stop("NMF::ibterms is a pure virtual method of interface 'NMF'."
			," It should be overloaded in class '", class(object),"'.")
	}
)
#' Method for standard NMF models, which returns the integer vector that is 
#' stored in slot \code{ibterms} when a formula-based NMF model is instantiated. 
setMethod('ibterms', 'NMFstd', 
	function(object){
		object@ibterms
	}
)
#' @export
#' @rdname terms
setGeneric('icterms', function(object, ...) standardGeneric('icterms') )
#' Default pure virtual method that ensure a method is defined for concrete 
#' NMF model classes.
setMethod('icterms', 'NMF', 
	function(object, ...){
		stop("NMF::icterms is a pure virtual method of interface 'NMF'."
				," It should be overloaded in class '", class(object),"'.")
	}
)
#' Method for standard NMF models, which returns the integer vector that is 
#' stored in slot \code{icterms} when a formula-based NMF model is instantiated. 
setMethod('icterms', 'NMFstd', 
	function(object){
		object@icterms
	}
)
#' @export
#' @rdname terms
iterms <- function(object, ...){
	c(ibterms(object), icterms(object))
}

#' \code{nterms}, \code{nbterms}, and \code{ncterms} return, respectively, 
#' the number of all fixed terms, fixed basis terms and fixed coefficient terms 
#' in an NMF model.
#' In particular: i.e. \code{nterms(object) = nbterms(object) + ncterms(object)}.
#' @export
#' @rdname terms
nterms <- function(object){
	length(ibterms(object)) + length(icterms(object))
}
#' @export
#' @rdname terms
nbterms <- function(object){
	length(ibterms(object))
}
#' @export
#' @rdname terms
ncterms <- function(object){
	length(icterms(object))
}

#' \code{bterms} and \code{cterms} return, respectively, the primary data for 
#' fixed basis and coefficient terms in an NMF model -- as stored in slots 
#' \code{bterms} and \code{cterms} .
#' These are factors or numeric vectors which define fixed basis components, 
#' e.g., used for defining separate offsets for different \emph{a priori} groups 
#' of samples, or to incorporate/correct for some known covariate.
#' 
#' @export
#' @rdname terms
bterms <- function(object){
	object@bterms		
}
#' @export
#' @rdname terms
cterms <- function(object){
	object@cterms	
}

#' \code{ibasis} and \code{icoef} return, respectively, the 
#' indexes of all latent basis vectors and estimated coefficients within the 
#' basis or coefficient matrix of an NMF model.
#' @export
#' @rdname terms
ibasis <- function(object, ...){
	i <- 1:nbasis(object)
	if( length(idx <- ibterms(object, ...)) ) i[-idx]
	else i
}
#' @export
#' @rdname terms
icoef <- function(object, ...){
	i <- 1:nbasis(object)
	if( length(idx <- icterms(object, ...)) ) i[-idx]
	else i
}

#' @S3method t NMFstd
t.NMFstd <- function(x){
    # transpose and swap factors
    x <- t.NMF(x)
    # swap fixed terms
    bt <- bterms(x)
    ibt <- ibterms(x)
    x@bterms <- cterms(x)
    x@ibterms <- icterms(x)
    x@cterms <- bt
    x@icterms <- ibt
    
    # returns
    x
}
