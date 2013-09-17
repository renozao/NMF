# Factory/Constructor Methods for NMF models
# 
# Author: Renaud Gaujoux
# Creation: 03 Jul 2012
###############################################################################

#' @include NMFstd-class.R
#' @include NMFns-class.R
#' @include NMFOffset-class.R
NULL

#' Factory Methods NMF Models
#'
#' \code{nmfModel} is a S4 generic function which provides a convenient way to 
#' build NMF models. 
#' It implements a unified interface for creating \code{NMF} objects from any 
#' NMF models, which is designed to resolve potential dimensions inconsistencies.
#' 
#' All \code{nmfModel} methods return an object that inherits from class \code{NMF}, 
#' that is suitable for seeding NMF algorithms via arguments \code{rank} or 
#' \code{seed} of the \code{\link{nmf}} method, in which case the factorisation 
#' rank is implicitly set by the number of basis components in the seeding 
#' model (see \code{\link{nmf}}).
#' 
#' For convenience, shortcut methods and internal conversions for working on 
#' \code{data.frame} objects directly are implemented.
#' However, note that conversion of a \code{data.frame} into a \code{matrix} 
#' object may take some non-negligible time, for large datasets. 
#' If using this method or other NMF-related methods several times, consider 
#' converting your data \code{data.frame} object into a matrix once for good, 
#' when first loaded.
#' 
#' @param rank specification of the target factorization rank 
#' (i.e. the number of components). 
#' @param target an object that specifies the dimension of the estimated target matrix.
#' @param ... extra arguments to allow extension, that are passed down to the 
#' workhorse method \code{nmfModel,numeric.numeric}, where they are used to 
#' initialise slots specific to the instantiating NMF model class.
#'
#' @return an object that inherits from class \code{\linkS4class{NMF}}. 
#' @family NMF-interface
#' @export
#' @inline
setGeneric('nmfModel', function(rank, target=0L, ...) standardGeneric('nmfModel'))

#' Main factory method for NMF models
#' 
#' This method is the workhorse method that is eventually called by all other methods.
#' See section \emph{Main factory method} for more details.
#' 
#' @param ncol a numeric value that specifies the number 
#' of columns of the target matrix, fitted the NMF model.
#' It is used only if not missing and when argument \code{target} is a single 
#' numeric value.
#' @param model the class of the object to be created. 
#' It must be a valid class name that inherits from class \code{NMF}. 
#' Default is the standard NMF model \code{\linkS4class{NMFstd}}.
#' @param W value for the basis matrix. 
#' \code{data.frame} objects are converted into matrices with \code{\link{as.matrix}}.
#' @param H value for the mixture coefficient matrix
#' \code{data.frame} objects are converted into matrices with \code{\link{as.matrix}}.
#' @param force.dim logical that indicates whether the method should try 
#' lowering the rank or shrinking dimensions of the input matrices to 
#' make them compatible 
#' @param order.basis logical that indicates whether the basis components should
#' reorder the rows of the mixture coefficient matrix to match the order of the 
#' basis components, based on their respective names. It is only used if the 
#' basis and coefficient matrices have common unique column and row names 
#' respectively.
#' 
#' @section Main factory method:
#' The main factory engine of NMF models is implemented by the method with 
#' signature \code{numeric, numeric}.
#' Other factory methods provide convenient ways of creating NMF models from e.g. a 
#' given target matrix or known basis/coef matrices (see section \emph{Other Factory Methods}). 
#' 
#' This method creates an object of class \code{model}, using the extra 
#' arguments in \code{...} to initialise slots that are specific to the given model.
#' 
#' All NMF models implement get/set methods to access the matrix factors 
#' (see \code{\link{basis}}), which are called to initialise them from arguments 
#' \code{W} and \code{H}.
#' These argument names derive from the definition of all built-in models that 
#' inherit derive from class \code{\linkS4class{NMFstd}}, which has two slots, 
#' \var{W} and \var{H}, to hold the two factors -- following the notations used 
#' in \cite{Lee1999}.
#' 
#' If argument \code{target} is missing, the method creates a standard NMF 
#' model of dimension 0x\code{rank}x0.
#' That is that the basis and mixture coefficient matrices, \var{W} and \var{H}, 
#' have dimension 0x\code{rank} and \code{rank}x0 respectively.
#' 
#' If target dimensions are also provided in argument \code{target} as a 
#' 2-length vector, then the method creates an \code{NMF} object compatible to 
#' fit a target matrix of dimension \code{target[1]}x\code{target[2]}.
#' That is that the basis and mixture coefficient matrices, \var{W} and \var{H}, 
#' have dimension \code{target[1]}x\code{rank} and \code{rank}x\code{target[2]} 
#' respectively.
#' The target dimensions can also be specified using both arguments \code{target} 
#' and \code{ncol} to define the number of rows and the number of columns of the 
#' target matrix respectively.
#' If no other argument is provided, these matrices are filled with NAs.
#' 
#' If arguments \code{W} and/or \code{H} are provided, the method creates a NMF 
#' model where the basis and mixture coefficient matrices, \var{W} and \var{H}, 
#' are initialised using the values of \code{W} and/or \code{H}.
#' 
#' The dimensions given by \code{target}, \code{W} and \code{H}, must be compatible.
#' However if \code{force.dim=TRUE}, the method will reduce the dimensions to the achieve 
#' dimension compatibility whenever possible.
#' 
#' When \code{W} and \code{H} are both provided, the \code{NMF} object created is 
#' suitable to seed a NMF algorithm in a call to the \code{\link{nmf}} method.
#' Note that in this case the factorisation rank is implicitly set by the number 
#' of basis components in the seed.
#' 
#' @examples 
#'
#' # data
#' n <- 20; r <- 3; p <- 10  
#' V <- rmatrix(n, p) # some target matrix
#' 
#' # create a r-ranked NMF model with a given target dimensions n x p as a 2-length vector
#' nmfModel(r, c(n,p)) # directly
#' nmfModel(r, dim(V)) # or from an existing matrix <=> nmfModel(r, V)
#' # or alternatively passing each dimension separately
#' nmfModel(r, n, p)
#' 
#' # trying to create a NMF object based on incompatible matrices generates an error
#' w <- rmatrix(n, r) 
#' h <- rmatrix(r+1, p)
#' try( new('NMFstd', W=w, H=h) )
#' try( nmfModel(w, h) )
#' try( nmfModel(r+1, W=w, H=h) )
#' # The factory method can be force the model to match some target dimensions
#' # but warnings are thrown
#' nmfModel(r, W=w, H=h)
#' nmfModel(r, n-1, W=w, H=h)
#' 
setMethod('nmfModel', signature(rank='numeric', target='numeric'),
	function(rank, target, ncol=NULL, model='NMFstd', W, H, ..., force.dim=TRUE, order.basis=TRUE){
		
		if( is.null(model) ) model <- 'NMFstd'
		# check validity of the provided class
		if( !isClass(model) ) stop("nmfModel - Invalid model name: class '", model,"' is not defined.")
		if( !extends(model, 'NMF') ) stop("nmfModel - Invalid model name: class '", model,"' does not extend class 'NMF'.")
		
		# check the validity of the target
		if( length(target) == 0 ) stop('nmfModel - Invalid dimensions: `target` must be at least of length 1')
		if( length(target) > 2 ) stop('nmfModel - Invalid dimensions: `target` must be at most of length 2')
		if( !missing(ncol) && !is.null(ncol) && (!is.vector(ncol) || length(ncol) > 1 || !is.numeric(ncol) || ncol<0 ) )
			stop('nmfModel - Invalid dimensions: `ncol` must be a single nonnegative integer')
				
		# compute the target dimension
		target <- as.integer(target)
		n <- target[1]
		m <- if( length(target) == 2 ) target[2] 
			 else if( !missing(ncol) && !is.null(ncol) ) ncol
			 else if( !missing(H) ) ncol(H)
	 		 else n 
	 	if( n < 0 )
			stop("nmfModel - Invalid target number of rows: nonnegative value expected")
		if( m < 0 )
			stop("nmfModel - Invalid target number of columns: nonnegative value expected")
		# force rank to be an integer
		r <- as.integer(rank)
		
		# check the validity of the rank
		if( length(r) != 1 ) stop("Invalid argument 'rank': single numeric expected")
		if( r < 0 ) stop("nmfModel - Invalid argument 'rank': nonnegative value expected")
		
		# do not allow dimension incompatibility if required
		if( !force.dim && !missing(W) && !missing(H) && ncol(W) != nrow(H) ){
			stop('nmfModel - Invalid number of columns in the basis matrix [', ncol(W), ']: '
						, 'it should match the number of rows in the mixture coefficient matrix [', nrow(H), ']')
		}
				
		# build dummy compatible W and H if necessary
		W.was.missing <- FALSE
		if( missing(W) ){
			W <- matrix(as.numeric(NA), n, r)
			W.was.missing <- TRUE
		}
		else{
			if( is.vector(W) ) # convert numerical vectors into a matrix
				W <- matrix(W, n, r)
			else if( is.data.frame(W) ) # convert data.frame into matrix
				W <- as.matrix(W) 
			
			if( r == 0 ) r <- ncol(W)
			else if( r < ncol(W) ){
				if( !force.dim ){
					stop('nmfModel - Invalid number of columns in the basis matrix [', ncol(W), ']: ',
							'it should match the factorization rank [', r, ']')
				}
				warning("Objective rank is [",r,"] lower than the number of columns in W [",ncol(W),"]: "
						, "only the first ", r," columns of W will be used")
				W <- W[,1:r, drop=FALSE]				
			}
			else if( r > ncol(W) ){
				stop("nmfModel - Objective rank [",r,"] is greater than the number of columns in W [",ncol(W),"]")
			}
			
			# resolve consistency with target
			if( n == 0 ) n <- nrow(W)
			else if( n < nrow(W) ){
				if( !force.dim ){
					stop('nmfModel - Invalid number of rows in the basis matrix [', nrow(W), ']: '
							, 'it should match the target number of rows [', n, ']')
				}
				
				warning("nmfModel - Number of rows in target is lower than the number of rows in W [",nrow(W),"]: ",
							"only the first ", n," rows of W will be used")
				W <- W[1:n, , drop=FALSE]				
			}
			else if( n > nrow(W) ){
				stop("nmfModel - Number of rows in target [",n,"] is greater than the number of rows in W [",nrow(W),"]")
			}
		}
		
		if( missing(H) ) 
			H <- matrix(as.numeric(NA), ncol(W), m)
		else{
			# convert numerical vectors into a matrix
			if( is.vector(H) )
				H <- matrix(H, r, m)
			else if( is.data.frame(H) ) # convert data.frame into matrix
				H <- as.matrix(H)
			
			if( r == 0 ) r <- nrow(H)
			else if( r < nrow(H) ){
				if( !force.dim ){
					stop('nmfModel - Invalid number of rows in the mixture coefficient matrix [', nrow(H), ']: '
						, 'it should match the factorization rank [', r, ']')
				}
				warning("nmfModel - Objective rank [",r,"] is lower than the number of rows in H [",nrow(H),"]: "
								, "only the first ", r," rows of H  will be used")
				H <- H[1:r,, drop=FALSE]				
			}
			else if( r > nrow(H) ) stop("nmfModel - Objective rank [",r,"] is greater than the number of rows in H [",nrow(H),"]")
			# force dummy W to be at least compatible with H
			if( W.was.missing ) W <- matrix(as.numeric(NA), n, r)

			# resolve consistency with target
			if( m == 0 ) m <- ncol(H)
			else if( m < ncol(H) ){
				if( !force.dim ){
					stop('nmfModel - Invalid number of columns in the mixture coefficient matrix [', ncol(H), ']:'
						, ' it should match the target number of columns [', m, ']')
				}
				
				warning("nmfModel - Number of columns in target is lower than the number of columns in H [",ncol(H),"]:"
								, " only the first ", m," columns of H will be used")
				H <- H[, 1:m, drop=FALSE]				
			}
			else if( m > ncol(H) ){ 
				stop("nmfModel - Number of columns in target [",m,"]"
						," is greater than the number of columns in H [",ncol(H),"]")
			}
		}
		
		# check validity of matrices W and H (only if one of the target dimension is not null)
		if( n + m > 0 ){
			if( nrow(W) != n ) stop('nmfModel - Invalid number of rows for W: should match number of rows in target [', n, ']')
			if( ncol(W) != r ) stop('nmfModel - Invalid number of columns for W: should match factorization rank [', r, ']')
			if( nrow(H) != r ) stop('nmfModel - Invalid number of rows for H: should match factorization rank [', r, ']')
			if( ncol(H) != m ) stop('nmfModel - Invalid number of columns for H: should match number of columns in target [', m, ']')
		}
		
		# build and return a dummy NMF object
		nmf.debug('nmfModel', "Instantiate NMF model:", model)
		res <- new(model, ...)
		nmf.debug('nmfModel', "Set factors in model:", model)
		# set the dimnames if possible
		cW <- !is.null(colnames(W))
		rH <- !is.null(rownames(H))
		if( cW && !rH )# use colnames of W as basisnames
			rownames(H) <- colnames(W)
		else if( !cW && rH )# use rownames of H as basisnames
			colnames(W) <- rownames(H)
		else if( cW && rH ){# try to match names or use colnames of W (with a warning)
			
			# reorder as in the basis matrix if it makes sense, i.e. if the names are the same
			if( order.basis && !anyDuplicated(rownames(H)) && length(setdiff(rownames(H), colnames(W)))==0 ){ 
				H <- H[match(rownames(H), colnames(W)),]
			}
			else{
				rownames(H) <- colnames(W)
				warning("nmfModel - The rownames of the mixture matrix were set to match the colnames of the basis matrix")
			}
			
		}
		# set the basis and coef matrices
		.basis(res) <- W; .coef(res) <- H
		# check validity
		validObject(res)

		# return the model
		res
	}
)

#' Creates an empty NMF model of a given rank.
#' 
#' This call is equivalent to \code{nmfModel(rank, 0L, ...)}, which 
#' creates \emph{empty} \code{NMF} object with a basis and mixture coefficient matrix  
#' of dimension 0 x \code{rank} and \code{rank} x 0 respectively.
#' 
#' @seealso \code{\link{is.empty.nmf}}
#' @examples 
#' ## Empty model of given rank
#' nmfModel(3)
#'  
setMethod('nmfModel', signature(rank='numeric', target='missing'),
		function(rank, target, ...){
			nmfModel(rank, 0L, ...)
		}
)

#' Creates an empty NMF model of null rank and a given dimension. 
#' 
#' This call is equivalent to \code{nmfModel(0, target, ...)}.
#' 
#' @examples 
#' nmfModel(target=10) #square
#' nmfModel(target=c(10, 5))
#' 
setMethod('nmfModel', signature(rank='missing', target='ANY'),
		function(rank, target, ...){
			nmfModel(0L, target, ...)
		}
)

#' Creates an empty NMF model of null rank and given dimension. 
#' 
#' This call is equivalent to \code{nmfModel(0, target, ...)}, and is meant for 
#' internal usage only.
setMethod('nmfModel', signature(rank='NULL', target='ANY'),
		function(rank, target, ...){
			nmfModel(0L, target, ...)
		}
)

#' Creates an empty NMF model or from existing factors
#' 
#' This method is equivalent to \code{nmfModel(0, 0, ..., force.dim=FALSE)}.	
#' This means that the dimensions of the NMF model will be taken from the optional 
#' basis and mixture coefficient arguments \code{W} and \code{H}.
#' An error is thrown if their dimensions are not compatible.
#' 
#' Hence, this method may be used to generate an NMF model from existing factor 
#' matrices, by providing the named arguments \code{W} and/or \code{H}: 
#' 
#' \code{nmfModel(W=w)} or \code{nmfModel(H=h)} or \code{nmfModel(W=w, H=h)} 
#' 
#' Note that this may be achieved using the more convenient interface is 
#' provided by the method \code{nmfModel,matrix,matrix} (see its dedicated description).
#' 
#' See the description of the appropriate method below.
#' 
#' @examples 
#' 
#' # Build an empty NMF model 
#' nmfModel()
#' 
#' # create a NMF object based on one random matrix: the missing matrix is deduced
#' # Note this only works when using factory method NMF 
#' n <- 50; r <- 3; 
#' w <- rmatrix(n, r) 
#' nmfModel(W=w)
#' 
#' # create a NMF object based on random (compatible) matrices
#' p <- 20
#' h <- rmatrix(r, p)
#' nmfModel(H=h)
#' 
#' # specifies two compatible matrices
#' nmfModel(W=w, H=h)
#' # error if not compatible
#' try( nmfModel(W=w, H=h[-1,]) ) 
#' 
setMethod('nmfModel', signature(rank='missing', target='missing'),
		function(rank, target, ...){
			# build an a priori empty model (extra args may provide the true dimension)
			# NB: do not allow dimension incompatibilities
			nmfModel(0L, 0L, ..., force.dim=FALSE)
		}
)

#' Creates an NMF model compatible with a target matrix.
#' 
#' This call is equivalent to \code{nmfModel(rank, dim(target), ...)}.
#' That is that the returned NMF object fits a target matrix of the same 
#' dimension as \code{target}.
#' 
#' Only the dimensions of \code{target} are used to construct the \code{NMF} object. 
#' The matrix slots are filled with \code{NA} values if these are not specified 
#' in arguments \code{W} and/or \code{H}.
#' However, dimension names are set on the return NMF model if present in 
#' \code{target} and argument \code{use.names=TRUE}.
#' 
#' @param  use.names a logical that indicates whether the dimension names of the 
#' target matrix should be set on the returned NMF model. 
#' 
#' @inline
#' @examples
#'  
#' # create a r-ranked NMF model compatible with a given target matrix
#' obj <- nmfModel(r, V)
#' all(is.na(basis(obj)))
#' 
setMethod('nmfModel', signature(rank='numeric', target='matrix'),
		function(rank, target, ..., use.names=TRUE){
			# build an object compatible with the target's dimensions
			res <- nmfModel(rank, dim(target), ...)
			
			# try to set dimnames if it makes sense: 
			# set on target and not somehow already set on the result
			if( use.names && !is.null(dimnames(target)) ){
				dn <- dimnames(res)
				if( is.null(dn) )
					dn <- list(NULL, NULL, NULL)
				if( is.null(rownames(res)) && !is.null(rownames(target)) )
					dimnames(res) <- c(dimnames(target)[1], dn[2:3])
				if( is.null(colnames(res)) && !is.null(colnames(target)) )					
					dimnames(res) <- c(dimnames(res)[1], dimnames(target)[2], dimnames(res)[3])				
			}
			res
		}	
)

#' Creates an NMF model based on two existing factors.
#' 
#' This method is equivalent to \code{nmfModel(0, 0, W=rank, H=target..., force.dim=FALSE)}.
#' This allows for a natural shortcut for wrapping existing \strong{compatible} 
#' matrices into NMF models:
#' \samp{nmfModel(w, h)}
#' 
#' Note that an error is thrown if their dimensions are not compatible.
#' 
#' @examples 
#' ## From two existing factors
#' 
#' # allows a convenient call without argument names
#' w <- rmatrix(n, 3); h <- rmatrix(3, p)
#' nmfModel(w, h)
#' 
#' # Specify the type of NMF model (e.g. 'NMFns' for non-smooth NMF)
#' mod <- nmfModel(w, h, model='NMFns')
#' mod
#' 
#' # One can use such an NMF model as a seed when fitting a target matrix with nmf()
#' V <- rmatrix(mod)
#' res <- nmf(V, mod)
#' nmf.equal(res, nmf(V, mod))
#' 
#' # NB: when called only with such a seed, the rank and the NMF algorithm 
#' # are selected based on the input NMF model. 
#' # e.g. here rank was 3 and the algorithm "nsNMF" is used, because it is the default 
#' # algorithm to fit "NMFns" models (See ?nmf). 
#' 
setMethod('nmfModel', signature(rank='matrix', target='matrix'),
		function(rank, target, ...){
			# use rank and target as W and H respectively
			# NB: do not allow dimension incompatibilities
			nmfModel(0L, 0L, W=rank, H=target, ..., force.dim=FALSE)
			
		}	
)

#' Same as \code{nmfModel('matrix', 'matrix')} but for \code{data.frame} objects,
#' which are generally produced by \code{\link{read.delim}}-like functions.
#' 
#' The input \code{data.frame} objects are converted into matrices with 
#' \code{\link{as.matrix}}.
setMethod('nmfModel', signature(rank='data.frame', target='data.frame'),
	function(rank, target, ...){
		nmfModel(as.matrix(rank), as.matrix(target), ...)
	}
)

#' Creates an NMF model with arguments \code{rank} and \code{target} swapped.
#' 
#' This call is equivalent to \code{nmfModel(rank=target, target=rank, ...)}.
#' This allows to call the \code{nmfModel} function with arguments \code{rank} 
#' and \code{target} swapped. 
#' It exists for convenience: 
#' \itemize{
#'  \item allows typing \code{nmfModel(V)} instead of \code{nmfModel(target=V)} to create 
#' a model compatible with a given matrix \code{V} (i.e. of dimension \code{nrow(V), 0, ncol(V)})
#' \item one can pass the arguments in any order (the one that comes to the user's mind first) 
#' 	and it still works as expected.
#' }
#' 
#' @examples
#' ## swapped arguments `rank` and `target`
#' V <- rmatrix(20, 10)
#' nmfModel(V) # equivalent to nmfModel(target=V)
#' nmfModel(V, 3) # equivalent to nmfModel(3, V) 
#' 
setMethod('nmfModel', signature(rank='matrix', target='ANY'),
		function(rank, target, ...){
			if( missing(target) ) target <- NULL
			# call nmfModel with swapping the arguments
			nmfModel(target, rank, ...)
			
		}	
)

#' Simple Parsing of Formula
#' 
#' Formula parser for formula-based NMF models.
#' 
#' @param x formula to parse
#' @return a list with the following elements:
#' \item{response}{ logical that indicates if the formula has a response term.}
#' \item{y}{ name of the response variable.}
#' \item{x}{ list of regressor variable names.}
#' \item{n}{ number of regressor variables.}
#' 
#' @keywords internal
parse_formula <- function(x){
	
	res <- list()
	# parse formula
	f <- as.character(x)
	hasResponse <- length(f) == 3L
	# response
	res$response <- hasResponse
	res$y <- if( hasResponse ) f[2L]
	# regressors
	reg <- if( hasResponse ) f[3L] else f[2L]
	res$x <- strsplit(reg, ' ')[[1]]
	res$n <- length(res$reg)
	# as a tring
	res$string <- paste(res$y, '~', reg, collapse='')
	
	res	
}
#' Build a formula-based NMF model, that can incorporate fixed basis or 
#' coefficient terms.
#' 
#' @param data Optional argument where to look for the variables used in the 
#' formula.
#' @param no.attrib logical that indicate if attributes containing data related
#' to the formula should be attached as attributes. 
#' If \code{FALSE} attributes \code{'target'} and \code{'formula'} contain the 
#' target matrix, and a list describing each formula part (response, regressors, 
#' etc.). 
#' 
#' @inline 
#' 
#' @examples
#' 
#' # empty 3-rank model
#' nmfModel(~ 3)
#' 
#' # 3-rank model that fits a given data matrix
#' x <- rmatrix(20,10)
#' nmfModel(x ~ 3)
#' 
#' # add fixed coefficient term defined by a factor
#' gr <- gl(2, 5)
#' nmfModel(x ~ 3 + gr)
#' 
#' # add fixed coefficient term defined by a numeric covariate
#' nmfModel(x ~ 3 + gr + b, data=list(b=runif(10)))
#' 
#' # 3-rank model that fits a given ExpressionSet (with fixed coef terms)
#' e <- ExpressionSet(x)
#' pData(e) <- data.frame(a=runif(10))
#' nmfModel(e ~ 3 + gr + a) # `a` is looked up in the phenotypic data of x pData(x)
#' 
setMethod('nmfModel', signature(rank='formula', target='ANY'),
	function(rank, target, ..., data=NULL, no.attrib=FALSE){
		
		# missing target is NULL
		if( missing(target) ) target <- NULL
		
		# data is a model class name (passed from nmf)
		if( is.character(data) ){
			model <- data
			data <- NULL
		}else model <- NULL
		
		# parse formula
		f <- parse_formula(rank)
		enclos <- environment(rank)
		
		rank <- 0L
		if( is.vector(target) && is.numeric(target) ){
			rank <- target
			target <- NULL
		}
		
		# utility function to merge data and pData
		merge_pdata <- function(x, data){
			pd <- pData(x)
			if( length(pd) ){
				if( is.null(data) ) pd
				else{
					cbind(data, pd)
				}
			}else data
		}
		
		# determine formula data
		if( is.null(data) ){
			# target data.frame taken as data if a response variable if defined
			if( is.data.frame(target) && f$response ){
				data <- target
				target <- NULL
			}else if( is.environment(target) ){ # use target as enclosure
				enclos <- target
				target <- NULL 
			}
		}
		
		# determine target matrix:
		X <- 0L
		# if a response term is present, lookup target data in other arguments
		if( f$response ){
			X <- eval(parse(text=f$y), enclos)
			if( is.eset(target) && !identical(X, target) ){
				warning("Conflicting response term and target: the ExpressionSet in `target` will only be used for covariates.")
				data <- merge_pdata(target, data)
			}
		} 
		else if( is.null(target) ){
			# no response, no target: try ExpressionSet in data 
			if( is.eset(data) ){
				X <- exprs(data)
			}
		}else{
			X <- target
		}
		
		# merge data and pData from ExpressionSet target
		if( is.eset(X) ){
			data <- merge_pdata(X, data)
			X <- exprs(X)
		}
		
		r <- rank
		cterms <- bterms <- list()
		
		# dimensions are also inferred from the formula
		n <- if( identical(X, 0L) ) 0L else nrow(X)
		p <- if( identical(X, 0L) ) 0L else ncol(X)
		
		for( v in f$x ){
			if( grepl("^[0-9]+$", v) ){
				if( rank == 0L ){ # rank not specified in target 
					r <- as.numeric(v)
				}else{
					warning("NMF::nmfModel - Discarding rank specified in the formula [", v,"]:"
							, " using value specified in target rank instead [", rank, "].")
				}
			}else if( grepl("^[+-]$", v) ) next
			else {
				val <- eval(parse(text=v), data, enclos)
                .add_term <- function(v, val, type = NULL){
    				if( p==0L || length(val) ==  p || identical(type, 'coef') ){
    					cterms[[v]] <<- val
    					if( p==0L ) p <<- length(val)
    				}else if( n==0L || length(val) ==  n  || identical(type, 'basis') ){
    					bterms[[v]] <<- val
    					if( n==0L ) n <<- length(val)
    				}else
    					stop("Invalid", type," term '", v, "' length [", length(val), "]:"
                                , " length must either be the number of target columns [", p, "]"
    							, " or rows [", n, "]")
                }
                
                if( is.null(dim(val)) ) .add_term(v, val)
                else if( n == 0L || nrow(val) == n ){
                    lapply(1:ncol(val), function(i){
                        if( !is.null(cname <- colnames(val)[i]) && nzchar(cname) ) vname <- cname
                        else vname <- paste0(v, i)
                        .add_term(vname, val[, i], type = 'basis')   
                    })
                }else{
                    # special handling of data.frames: 
                    # -> coef terms are passed as column variables
                    if( is.data.frame(val) && (p == 0L || nrow(val) == p)){
                        val <- t(val)
                    } 
                    if( p == 0L || ncol(val) == p ){
                    lapply(1:nrow(val), function(i){
                        if( !is.null(cname <- rownames(val)[i]) && nzchar(cname) ) vname <- cname
                        else vname <- paste0(v, i)
                        .add_term(vname, val[i, ], type = 'coef')   
                    })
                    }else{
                        stop("Incompatible matrix-like term '", v, "' dimensions [", str_dim(val), "]:"
                                , " number of rows or columns must match the ones of the target matrix [", str_dim(X, dims = c(n, p)) ,"]")
                    }                
                }
			}
		}
		# try to fixup X if possible
		if( identical(X, 0L) ) X <- c(n, p)
		
		# call nmfModel with cterms
		if( hasArg(model) || is.null(model) ) object <- nmfModel(r, X, ...)
		else object <- nmfModel(r, X, ..., model=model)
		# set fixed basis terms
		if( length(bterms) ){
			bterms(object) <- as.data.frame(bterms)	
		}
		# set fixed coef terms
		if( length(cterms) ){
			cterms(object) <- as.data.frame(cterms)
		}
		
		# valid object
		validObject(object)
		# attach formula data
		if( !no.attrib ){
			attr(object, 'target') <- X
			attr(object, 'formula') <- f
		}
		# return object
		object
	}
)

#' Listing NMF Models
#' 
#' \code{nmfModels} lists all available NMF models currently defined that can be 
#' used to create NMF objects, i.e. -- more or less -- all S4 classes that 
#' inherit from class \code{\linkS4class{NMF}}.
#' 
#' @param builtin.only logical that indicates whether only built-in NMF models, 
#' i.e. defined within the NMF package, should be listed.
#' 
#' @return a list
#' 
#' @export
#' @family NMF-interface
#' @rdname nmfModel
#' @examples
#' 
#' # show all the NMF models available (i.e. the classes that inherit from class NMF)
#' nmfModels()
#' # show all the built-in NMF models available
#' nmfModels(builtin.only=TRUE)
#' 
nmfModels <- function(builtin.only=FALSE){
	
	if( builtin.only ) return( .nmf.Models.Builtin )
	
	# return all subclasses of class 'NMF' (minus class 'NMFfit' and its subclasses)
	models <- names(methods::getClass('NMF')@subclasses)
	models.wraps <- c('NMFfit', names(methods::getClass('NMFfit')@subclasses))
	return( models[!is.element(models, models.wraps)] )
	
}

###% Initialization function for NMF models
.nmf.Models.Builtin <- NULL
.init.nmf.models <- function(){	
	.nmf.Models.Builtin <<- nmfModels()
}

