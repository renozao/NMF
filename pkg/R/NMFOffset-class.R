#' @include NMF-class.R
NA

#' Class for NMF factorization with offset.
#'
#' @references Badea (2008)
#'	, 'Extracting Gene Expression Profiles Common to Colon and Pancreatic Adenocarcinoma Using Simultaneous Nonnegative Matrix Factorization
#'	, Pacific Symposium on Biocomputing 13:279-290 (2008)
#'
setClass('NMFOffset'
	, representation(
				offset = 'numeric' # offset vector
				)
	, contains = 'NMFstd'
  	, prototype=prototype(
  				offset = numeric()
				)
	
)

#' Show method for objects of class \code{NMFOffset}
setMethod('show', signature(object='NMFOffset'), 
	function(object)
	{
		callNextMethod()
		cat("offset: ")
		if( length(object@offset) > 0 ){
			cat('[', head(object@offset, 5)
				, if( length(object@offset) > 5 ) "..." else NULL
				, ']')
		}
		else cat('none')
		cat("\n")
	}
)

#' Initialize method for class \code{NMFOffset}.
#' 
#' It tries to correct slot \code{offset} with a value consistent with the parent \code{NMF} object.
#' It will complete the offset with zeros to get the length equal to the number of rows in slot \code{W}.
#'  
setMethod("initialize", 'NMFOffset', 
		function(.Object, ..., offset){			
			.Object <- callNextMethod()
			# correct the offset slot if possible
			if( missing(offset) ) offset <- numeric()
			if( !is.numeric(offset) ) stop("Unvalid value for parameter 'offset': a numeric vector is expected")			
			 
			# force length to be consistent with the factorization's dimension
			n <- nrow(.Object)
			if( n > 0 ) .Object@offset <- c( offset, rep(NA, max(0, n - length(offset))) )[1:n]
			
			# return the initialized valid object
			.Object
		}
)

#' Returns the offset vector of an \code{NMFOffset} instance.
#' 
#' @param object an instance of class \code{NMFOffset}
#' @returnType numeric
#' @return the offset vector of \code{object} (from slot \code{offset})
#' 
if ( !isGeneric('offset') ) setGeneric('offset', package='stats')
setMethod('offset', signature(object='NMFOffset'), 
	function(object){
		object@offset
	}
)

#' Compute estimate for an NMFOffset object
setMethod('fitted', signature(object='NMFOffset'), 
	function(object, offset=object@offset){ 
		object@W %*% object@H + offset
	}
)

#' Initialize a random factorization with offset.
setMethod('rnmf', signature(x='NMFOffset', target='numeric'), 
function(x, target, ...){	
	
	# call the parent's 'rnmf' method to build a standard random NMF factorization
	res <- callNextMethod(x, target, ...)
		
	#Vc# Initialize a random offset of length the number of genes
	res@offset <- runif(nrow(res), min=0, max=max(basis(res), coef(res)));
	
	# return the initialized NMFOffset object
	res
})
