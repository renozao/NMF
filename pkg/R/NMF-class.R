#library(R.utils)

# declare old S3 class 'proc_time' to use it as a slot for class NMF 
setOldClass('proc_time', prototype=numeric())

################################
# Class: NMF
################################

###% Base virutal interface to store \strong{Non-negative Matrix Factorization} models.
setClass('NMF'
		, representation(
			misc = 'list' # misceleneaous data used during fitting
		)
		, contains = 'VIRTUAL')

###% Compute the fitted target matrix based on a \code{NMF} model.
###% 
###% The estimation depends on the underlying NMF model. For example in the standard model, the target
###% matrix is estimated by the matrix product \eqn{WH}. In other models, the estimate may depend on extra 
###% parameters/matrix (cf. Non-smooth NMF in \code{\link{NMFns-class}}).
###% 
###% @param x an object that inherits from class \code{NMF} 
###% @param ... extra parameters
###%
###% @return the fitted target matrix from object \code{x} 
###% @author Renaud Gaujoux \email{renaud@@cbio.uct.ac.za}
###% @export
setGeneric('fitted', package='stats')
setMethod('fitted', signature(object='NMF'),
		function(object){
			stop("NMF::fitted is a pure virtual method of interface 'NMF'. It should be overloaded in class '", class(object),"'.")
		}
)

###% Get/Set the matrix factors in a NMF model
###% 
###% \code{basis} and \code{basis<-} are S4 generic functions which respectively
###% extract and set the matrix of basis vectors (i.e. the first matrix factor)
###% of a NMF model.  For example, in the case of the standard NMF model \eqn{V
###% \equiv WH}, method \code{basis} will return matrix \eqn{W}.
###% 
###% @name basis/coef-methods: Accessing NMF matrix factors
###% @rdname basis-methods
###% @family NMF-interface
setGeneric('basis', function(object, ...) standardGeneric('basis') )
###% @autoRd
setMethod('basis', signature(object='NMF'),
		function(object){
			stop("NMF::basis is a pure virtual method of interface 'NMF'. It should be overloaded in class '", class(object),"'.")
		}
)
###% @autoRd
setGeneric('basis<-', function(object, value) standardGeneric('basis<-') )
###% @autoRd
setReplaceMethod('basis', signature(object='NMF', value='matrix'), 
		function(object, value){ 
			stop("NMF::basis<- is a pure virtual method of interface 'NMF'. It should be overloaded in class '", class(object),"'.")
		} 
)

###% @description \code{coef} and \code{coef<-} are S4 methods defined for the associated
###% generic functions from package \code{stats} (See \link[stats]{coef}).  They
###% respectively extract and set the matrix of mixture coefficients (i.e. the
###% second matrix factor) of a NMF model.  For example, in the case of the
###% standard NMF model \eqn{V \equiv WH}, method \code{coef} will return matrix
###% \eqn{H}.
###% 
###% Methods \code{coefficients} and \code{coefficients<-} are simple aliases for
###% methods \code{coef} and \code{coef<-} respectively.
###% 
###% @rdname basis-methods
setGeneric('coef', package='stats')
###% @autoRd
setMethod('coef', signature(object='NMF'),
		function(object){
			stop("NMF::coef is a pure virtual method of interface 'NMF'. It should be overloaded in class '", class(object),"'.")
		}
)
###% @autoRd
setGeneric('coef<-', function(object, value) standardGeneric('coef<-') )
###% @autoRd
setReplaceMethod('coef', signature(object='NMF', value='matrix'), 
		function(object, value){ 
			stop("NMF::coef<- is a pure virtual method of interface 'NMF'. It should be overloaded in class '", class(object),"'.")
		} 
)
###% @autoRd
setGeneric('coefficients', package='stats')
###% @autoRd
setMethod('coefficients', signature(object='NMF'),
		function(object){
			coef(object)
		}
)

###% @description \code{scoef} returns the mixture coefficient matrix of a NMF 
###% model with the columns scaled so that they sum up to a given value (1 by default).
###% 
###% @rdname basis-methods
setGeneric('scoef', function(object, ...) standardGeneric('scoef') )
###% @autoRd
setMethod('scoef', signature(object='NMF'),
	function(object, scale=1){
		sweep(coef(object), 2, scale * colSums(coef(object)), '/')
	}
)

###% Rescales an NMF Model.
###% 
###% Standard NMF models are identifiable modulo a scaling factor, meaning that the
###% basis vectors and basis profiles can be rescaled without changing the fitted 
###% values.
###% 
rescale <- function(x){
	wc <- colSums(basis(x))
	basis(x) <- sweep(basis(x), 2, wc, '/')
	coef(x) <- sweep(coef(x), 1, wc, '*')
	x
}


###% Returns a random sample of the implemented NMF model
setGeneric('rnmf', function(x, target, ...) standardGeneric('rnmf') )

###% Returns all the NMF model available
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

# Define the loading namespace
.PKG.NAMESPACE <- getLoadingNamespace(TRUE)
if( is.null(.PKG.NAMESPACE) )
	.PKG.NAMESPACE <- .GlobalEnv

###% Advanced usage of package NMF
###% 
###% The functions documented here provide advanced functionalities useful when
###% developing within the framework implemented in the NMF package.
###%
###% @details \emph{is.nmf} tests if an object is an NMF model or a class that extends
###% the class NMF.
###% 
###% NB for is.nmf: one has to work with the namespace as this function needs to return  
###% correct results when called in \code{.onLoad}.
###% See discussion on r-devel: \url{https://stat.ethz.ch/pipermail/r-devel/2011-June/061357.html}
###%
###% @name Advanced usage of package NMF 
###% @rdname advanced
is.nmf <- function(object){
	
	# load definition for base class NMF
	clref <- getClass('NMF', .Force=TRUE, where=.PKG.NAMESPACE)
	
	if( is.character(object) ){ # test object is a class that extends NMF		
		cl <- getClass(object, .Force=TRUE, where=.PKG.NAMESPACE)
		if( is.null(cl) )
			cl <- getClass(object, .Force=TRUE)
		extends(cl, clref)
	}else
		is(object, clref)
	
}

################################

# deprecated version of newNMF
newNMF <- function(...){
	.Deprecated('nmfModel')
	nmfModel(...)
}
###% Factory method for objects that inherit from class \code{NMF}.
###% 
###% @param rank the factorization rank (i.e. the number of metagenes). It is coerced into an integer.
###% @param target an object that defines the dimension of the target matrix to estimate.
###% @param ... extra parameters passed to method \code{new}.
###% 
###% @return an instance of class \code{class}. 
###% @seealso NMF-class 
###% @author Renaud Gaujoux \email{renaud@@cbio.uct.ac.za}
###% @export
setGeneric('nmfModel', function(rank, target=0, ...) standardGeneric('nmfModel'))
###% Main constructor method that is eventually called by all \code{NMF} methods. 
###% 
###% @param rank the factorization rank (i.e. the number of metagenes). It is coerced into an integer.
###% @param target a numeric vector (at most 2-length) giving the dimension of the target matrix to estimate. If \code{target} is of 1-length then a square target matrix is assumed.
###% @param class the class of the object to be created. It must be a valid class name that inherits from class \code{NMF}. Default is \code{'NMF'}.
###% @param ... extra parameters passed to class \code{class} default constructor.
###% @examples 
###% # create a NMF object of factorization rank 3 adapted to fit a 1000 x 30 matrix 
###% obj <- nmfModel3, c(1000,30))
###% obj
###% 
###% @export
setMethod('nmfModel', signature(rank='numeric', target='numeric'),
	function(rank, target, ncol=NULL, model='NMFstd', W, H, ..., force.dim=TRUE, order.basis=TRUE){
		
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
		if( !force.dim && !missing(W) && !missing(H) && ncol(W) != nrow(H) )
			stop('nmfModel - Invalid number of columns in the basis matrix (', ncol(W), '): it should match the number of rows in the mixture coefficient matrix [', nrow(H), ']')
				
		# build dummy compatible W and H if necessary
		W.was.missing <- FALSE
		if( missing(W) ){
			W <- matrix(NA, n, r)
			W.was.missing <- TRUE
		}
		else{
			# convert numerical vectors into a matrix
			if( is.vector(W) )
				W <- matrix(W, n, r)
			
			if( r == 0 ) r <- ncol(W)
			else if( r < ncol(W) ){
				if( !force.dim )
					stop('nmfModel - Invalid number of columns in the basis matrix (', ncol(W), '): it should match the factorization rank [', r, ']')
				
				warning("Objective rank is (",r,") lower than the number of columns in W (",ncol(W),"): only the first ", r," columns of W will be used")
				W <- W[,1:r, drop=FALSE]				
			}
			else if( r > ncol(W) ) stop("nmfModel - Objective rank (",r,") is greater than the number of columns in W (",ncol(W),")")
			
			# resolve consistency with target
			if( n == 0 ) n <- nrow(W)
			else if( n < nrow(W) ){
				if( !force.dim )
					stop('nmfModel - Invalid number of rows in the basis matrix (', nrow(W), '): it should match the target number of rows [', n, ']')
				
				warning("nmfModel - Number of rows in target is lower than the number of rows in W (",nrow(W),"): only the first ", n," rows of W will be used")
				W <- W[1:n, , drop=FALSE]				
			}
			else if( n > nrow(W) ) stop("nmfModel - Number of rows in target (",n,") is greater than the number of rows in W (",nrow(W),")")
		}
		
		if( missing(H) ) 
			H <- matrix(NA, ncol(W), m)
		else{
			# convert numerical vectors into a matrix
			if( is.vector(H) )
				H <- matrix(H, r, m)
			
			if( r == 0 ) r <- nrow(H)
			else if( r < nrow(H) ){
				if( !force.dim )
					stop('nmfModel - Invalid number of rows in the mixture coefficient matrix (', nrow(H), '): it should match the factorization rank [', r, ']')
				
				warning("nmfModel - Objective rank (",r,") is lower than the number of rows in H (",nrow(H),"): only the first ", r," rows of H  will be used")
				H <- H[1:r,, drop=FALSE]				
			}
			else if( r > nrow(H) ) stop("nmfModel - Objective rank (",r,") is greater than the number of rows in H (",nrow(H),")")
			# force dummy W to be at least compatible with H
			if( W.was.missing ) W <- matrix(NA, n, r)

			# resolve consistency with target
			if( m == 0 ) m <- ncol(H)
			else if( m < ncol(H) ){
				if( !force.dim )
					stop('nmfModel - Invalid number of columns in the mixture coefficient matrix (', ncol(H), '): it should match the target number of columns [', m, ']')
				
				warning("nmfModel - Number of columns in target is lower than the number of columns in H (",ncol(H),"): only the first ", m," columns of H will be used")
				H <- H[, 1:m, drop=FALSE]				
			}
			else if( m > ncol(H) ) stop("nmfModel - Number of columns in target (",m,") is greater than the number of columns in H (",ncol(H),")")
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
		basis(res) <- W; coef(res) <- H
		
		# return the model
		res
	}
)
###% Constructor method with \code{rank} as a single argument. 
###% 
###% Because the target dimension it creates a new \em{empty} \code{NMF} object where slots \code{W} 
###% and \code{H} have dimension 0 x \code{rank} and \code{rank} x 0 respectively.
###% 
###% @return an empty \code{NMF} object.
###% 
###% @examples
###% # create an empty NMF object of factorization rank 3
###% obj <- nmfModel3)
###% obj
###% 
###% @seealso is.empty.nmf 
###% @export
setMethod('nmfModel', signature(rank='numeric', target='missing'),
		function(rank, target, ...){
			nmfModel(rank, 0, ...)
		}
)

###% Creates an empty NMF object
setMethod('nmfModel', signature(rank='missing', target='ANY'),
		function(rank, target, ...){
			nmfModel(0, target, ...)
		}
)

###% Creates an empty NMF object
setMethod('nmfModel', signature(rank='NULL', target='ANY'),
		function(rank, target, ...){
			nmfModel(0, target, ...)
		}
)

###% List the currently defined NMF models
setMethod('nmfModel', signature(rank='missing', target='missing'),
		function(rank, target, ...){
			# build an a priori empty model (extra args may provide the true dimension)
			# NB: do not allow dimension incompatibilities
			nmfModel(0, 0, ..., force.dim=FALSE)
		}
)

###% Constructor method to create a new \code{NMF} object to fit a given target matrix.
###% 
###% Only the dimensions of \code{target} are used to construct the \code{NMF} object. The matrix slots 
###% \code{W} and \code{H} are filled  with \code{NA}. 
###% 
###% @return a \code{NMF} object with dimensions compatible with matrix \code{target}.
###% 
###% @examples
###% # create a NMF object of factorization rank 3 adapted to fit a given matrix
###% n <- 1000; m <-30 # the target matrix dimensions 
###% V <- matrix(runif(n*m), n, m) # create a random matrix 
###% obj <- nmfModel3, V)
###% obj
###% all(is.na(obj@@W))
###% 
###% @export
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

###% Swapped arguments rank/target when rank is a matrix 
setMethod('nmfModel', signature(rank='matrix', target='matrix'),
		function(rank, target, ...){
			# use rank and target as W and H respectively
			# NB: do not allow dimension incompatibilities
			nmfModel(0, 0, W=rank, H=target, ..., force.dim=FALSE)
			
		}	
)

###% Swapped arguments rank/target when rank is a matrix 
setMethod('nmfModel', signature(rank='matrix', target='ANY'),
		function(rank, target, ...){
			# call nmfModel with swapping the arguments
			nmfModel(target, rank, ...)
			
		}	
)

###% Show method for an object of class \code{NMF}.
setMethod('show', 'NMF', 
		function(object)
		{
			cat("<Object of class:", class(object), ">\n")
			cat("features:", nrow(object), "\n")
			cat("basis/rank:", nbasis(object), "\n")
			cat("samples:", ncol(object), "\n")
			# show the miscellaneous model parameters
			if( length(object@misc) > 0 ){
				cat("miscellaneous:\n")
				print(object@misc)
			}
		}
)


###% Dims method for an NMF object. It returns the dimension of the NMF model.
###% 
###% @note This method calls internally \code{\link{dims}}
###% 
###% @param x a object that inherits from class \code{NMF}
###% 
###% @return a 3-length numeric vector giving the dimension of the estimated target matrix, and 
###% the factorization rank (i.e. the number of metagenes).
###% @author Renaud Gaujoux \email{renaud@@cbio.uct.ac.za}
###% @seealso dims
###% @export
setMethod('dim', signature(x='NMF'), 
function(x){
	c(nrow(basis(x)), ncol(coef(x)), nbasis(x))	
})

###% Dimnames method for an NMF object. It returns the dimension names of the target matrix.
###% 
###% @note This method calls internally \code{\link{dims}}
###% 
###% @param x a object that inherits from class \code{NMF}
###% 
###% @return a 2-length numeric vector giving the dimension of the estimated target matrix. An attribute 
###% \code{'rank'} is set to the factorization rank (i.e. the number of metagenes).
###% @author Renaud Gaujoux \email{renaud@@cbio.uct.ac.za}
###% @seealso dims
###% @export
setMethod('dimnames', 'NMF', 
	function(x){
		b <- dimnames(basis(x))
		if( is.null(b) )
			b <- list(NULL, NULL)
		c <- dimnames(coef(x))
		if( is.null(c) )
			c <- list(NULL, NULL)
		l <- c(b[1],c[2],b[2])
		if( all(sapply(l, is.null)) ) NULL else l
	}
)

setReplaceMethod('dimnames', 'NMF', 
	function(x, value){
		if( is.list(value) ){
			if( length(value) == 0 )
				value <- NULL
			else if( length(value) == 1 )
				value <- c(value, list(NULL, NULL))			
			else if( length(value) == 2 ) # if only the two first dimensions reset the third one
				value <- c(value, list(NULL))
			else if( length(value)!=3 ) # check length of value
				stop("NMF::dimnames - invalid argument 'value' [a 2 or 3-length list is expected]")
		}
		
		dimnames(basis(x)) <- value[c(1,3)]		
		dimnames(coef(x)) <- value[c(3,2)]		
		x	
	}
)

###% Returns the factorization rank of a \code{NMF} object.
###% 
###% For a factorization \deqn{V \equiv WH,} the factorization rank is the number of columns of matrix \eqn{W}.
###% In the context of microarray data, the factorization rank is the number of metagenes. 
###% See \code{\link{NMF-class}} for more details.
###% 
###% @param x an object that inherits from class \code{NMF} 
###% @param ... extra parameters
###% 
###% @return a single numeric value
###% @author Renaud Gaujoux \email{renaud@@cbio.uct.ac.za}
###% @seealso NMF-class, dims
###% @export
setGeneric('nbasis', function(x, ...) standardGeneric('nbasis') )
setMethod('nbasis', signature(x='matrix'), 
	function(x)
	{
		attr(x, 'nbasis')
	}
)
setMethod('nbasis', signature(x='NMF'), 
	function(x)
	{
		ncol(basis(x))
	}
)

###% return the basis names
setGeneric('basisnames', function(x, ...) standardGeneric('basisnames') )
setMethod('basisnames', signature(x='NMF'), 
		function(x)
		{
			rownames(coef(x))
		}
)
###% set the basis names
setGeneric('basisnames<-', function(x, ..., value) standardGeneric('basisnames<-') )
setReplaceMethod('basisnames', signature(x='NMF'), 
		function(x, ..., value)
		{
			rownames(coef(x)) <- value
			colnames(basis(x)) <- value
			x
		}
)

###% Subsetting method for NMF objects
###% i = features, j = samples, k = basis
setMethod('[', 'NMF', 
	function (x, i, j, ..., drop = FALSE)
	{
		k <- NULL
		mdrop <- missing(drop)
		
		# compute number of arguments: x and drop are always passed
		Nargs <- nargs() - !mdrop
		single.arg <- FALSE
		k.notmissing <- FALSE
		if( !missing(i) && Nargs < 3L ){
			k <- i
			single.arg <- TRUE
		}
		else if( Nargs > 3L ){
			dots <- list(...)
			if( length(dots) != 1 )
				warning("NMF::[ - using only the first extra subset index, the remaining ", length(dots)-1," are discarded.")
			k <- dots[[1]]
			k.notmissing <- TRUE
		}
		
		# no indice was provided => return the object unchanged
		if ( missing(i) && missing(j) && !k.notmissing ) {
			# check if there is other arguments
			if (length(list(...)) != 0)
				stop("NMF::[] method - please specify which features, samples or basis to subset. See class?NMF.")
			# otherwise return the untouched object
			return(x)
		}
		
		# subset the rows of the basis matrix		
		if ( !missing(i) && !single.arg )			
			basis(x) <- basis(x)[i, , drop = FALSE]
		
		# subset the columns of mixture coefficient matrix		
		if (!missing(j)) 
			coef(x) <- coef(x)[, j, drop = FALSE]
		
		# subset the basis: columns of basis matrix and row of mixture coefficient matrix		
		if( single.arg )# return the selected metagenes (as a matrix)
			return(basis(x)[, k, drop = if( mdrop ) TRUE else drop])			
		else if( k.notmissing ){
			basis(x) <- basis(x)[, k, drop = FALSE]
			coef(x) <- coef(x)[k, , drop = FALSE]
		}
		
		# if drop is TRUE and only one dimension is missing then return affected matrix
		if( !single.arg && drop ){
			if( missing(i) && !missing(j) )
				return( drop(coef(x)) )
			else if( missing(j) && !missing(i) )
				return( drop(basis(x)) )
		}
		
		# return subset object
		return(x)
	}
)

###% Get/Set methods for slot 'misc'
setMethod('$', 'NMF', 
		function(x, name){ 
			x@misc[[name, exact=TRUE]]; 
		} 
)

setReplaceMethod('$', 'NMF',
		function(x, name, value) {
			x@misc[[name]] <- value
			x
		}
)

###% Checks whether an \code{NMF} object describes an empty NMF.
###%
###% An empty \code{NMF} object has slots \code{W} and \code{H} with 0 rows and 0 columns respectively.
###% Such objects are usually created to save memory space and often require a specific treatment. 
###% 
###% @note While an empty \code{NMF} estimates a target matrix of dimension \eqn{0 \times 0}, 
###% its factorization rank is usually not null.  
###%  
###% @param object An instance of NMF-class
###% 
###% @return \code{TRUE} if the NMF is empty, \code{FALSE} otherwise
###% 
setGeneric('is.empty.nmf', function(object) standardGeneric('is.empty.nmf') )
setMethod('is.empty.nmf', signature(object='NMF'), 
		function(object){
			nrow(object) == 0 && ncol(object) == 0
		}
)


###% Generates a random matrix of the same dimension as the target matrix of an NMF 
###% model
setMethod('rmatrix', 'NMF', 
	function(x, ...){
		rmatrix(nrow(x), ncol(x), ...)
	}
)

###% Initialize a random Nonnegative Matrix Factorization (NMF).
###%
###% The methods sets the slots \code{W} and \code{H} of \code{x} to matrices whose entries are drawn from the uniform distribution.
###%
###% @param x an instance of class \code{NMF} that will be used as template
###% @param target a pair of integer value representing the dimension of the target matrix OR the target matrix itself.
###% When a \code{matrix} is used as the target, the method uses the maximum of its entries as an upperbound.
###% @param r the rank of the factorization (i.e. the number of columns of the first factor)
###% @param ... extra parameters passed to \code{runif}
###% 
###% @return a modified an instance of \code{x}. This implies that the result is an instance of the same class as \code{x}.
###%
###% @seealso runif, NMF-class
###%
setMethod('rnmf', signature(x='NMF', target='numeric'), 
	function(x, target, ncol=NULL, keep.names=TRUE, ...){
		
		# store original dimnames
		if( keep.names ) dn <- dimnames(x)
		
		# valid parameter 'target'
		if( length(target) != 1 && length(target) != 2 )
			stop('NMF::rnmf - invalid target dimensions [length must be 1 or 2. Here length = ', length(target) ,']')
		if( any(is.na(target)) ) 
			stop('NMF::rnmf - invalid target dimensions [NA values in element(s): ', paste(which(is.na(target)), collapse=' and '), ']')		
		# shortcut for symetric case: provide only one dimension
		if( length(target) == 1 ){
			ncol <- if( !is.null(ncol) ){
				if( !is.numeric(ncol) || length(ncol) != 1 || is.na(ncol) )
					stop("NMF::rnmf - invalid argument `ncol`: must be a single numeric value")
				ncol
			}else target
			target <- c(target, ncol)
		}
		
		# retrieve dimension of the target matrix
		n <- target[1]; m <- target[2];
		# retrieve the factorization rank					
		r <- nbasis(x)
		
		#Vc# Initialize random matrix: W
		basis(x) <- matrix(runif(n*r, min=0, ...), n, r);
		#Vc# Initialize random matrix: H
		coef(x) <- matrix(runif(r*m, min=0, ...), r, m);
		
		# if one needs to keep the names (possibly or reducing/increasing) 
		if( keep.names && !is.null(dn) )
			dimnames(x) <- list(dn[[1]][1:n], dn[[2]][1:m], dn[[3]][1:r])
		
		# return the modified object
		x
	}
)

###% Method 'rnmf' for 'matrix' target objects: 
###% - use the dimension of 'target' as target dimensions
###% - its maximum entry as a limsup for drawing the random entries 
###%   
setMethod('rnmf', signature(x='ANY', target='matrix'), 
	function(x, target, use.dimnames=TRUE, ...){	
				
		# compute the upper-bound of the random entries (if not provided)
		no.na <- abs(target[!is.na(target)])
		max <- if( length(no.na) == 0 ) 1 else max(no.na)
		# build a random NMF with the dimensions of the target matrix upper-bounded by the target's maximum entry.
		res <- rnmf(x, dim(target), max=max, ...)
		# set the dimnames from the target matrix if necessary
		if( use.dimnames )
			dimnames(res) <- dimnames(target)
		
		# return result
		res
	}
)

setMethod('rnmf', signature(x='NMF', target='missing'), 
	function(x, target, ...){
		rnmf(x, c(nrow(x),ncol(x)), ...)
	}
)

###% Generate a random NMF model of given dimensions
setMethod('rnmf', signature(x='numeric', target='numeric'), 
		function(x, target, ncol=NULL, ...){		
			rnmf(nmfModel(x, target, ncol), ...)
		}
)

###% Returns the NMF model's name: i.e the class name
setGeneric('modelname', function(object, ...) standardGeneric('modelname'))
setMethod('modelname', signature(object='NMF'), 
		function(object)
		{
			as.character(class(object))
		}
)

###% Produces different kind of plots.
#setGeneric('plot', package='graphics')
#setMethod('plot', signature( x='NMF', y='missing'), 
#	function(x, y, type=c('hist', 'heatmap'), ...)
#	{
#		# retrieve what to plot	
#		type = match.arg(type)
#		
#		# save graphical parameters
#		oldpar = par(no.readonly=TRUE)
#		on.exit( {par(oldpar)} ) # reset the graphical parameters on exit
#		
#		if( what == 'heatmap' ){
#			#basicHM(metaprofiles(x), ...)
#			heatmap.2(metaprofiles(x), trace='none', ...)
#		}
#		else if( what == 'hist' ) hist(x, ...)
#						
#	}
#)
setGeneric('metaHeatmap', function(object, ...) standardGeneric('metaHeatmap') )
setMethod('metaHeatmap', signature(object='matrix'),
		function(object, type=c('plain', 'consensus'), class
				, unit.scaling=c('none', 'row', 'column'), palette="YlOrRd"
				, rev.palette=FALSE, show.prediction=TRUE, ...){
			
			.Deprecated('metaHeatmap', 'NMF', "The S4-Method 'metaHeatmap' is deprecated, use 'coefmap', 'basismap', 'consensusmap' or 'aheatmap' instead.")
			
			# load libary RColorBrewer
			library(RColorBrewer)
			
			# retreive the graphical parameters and match them to the sub-sequent call to 'heatmap.plus.2'
			graphical.params <- list(...)
			names(graphical.params) <- .match.call.args(names(graphical.params), 'heatmap.plus.2', in.fun='metaHeatmap', call='NMF::metaHeatmap')
			
			type <- match.arg(type)
			if( type == 'consensus' ){
				# set default graphical parameters for type 'consensus'
				graphical.params <- .set.list.defaults(graphical.params
						, distfun = function(x){ as.dist(1-x) }
						, main='Consensus matrix'
						, symm=TRUE
						, Rowv=TRUE
						, revC=TRUE
				)
				
				if( missing(palette) ) palette <- 'RdYlBu'
				if( missing(rev.palette) ) rev.palette <- TRUE
				if( missing(unit.scaling) ) unit.scaling <- 'none'
				show.prediction <- FALSE # not used for consensus matrices
			}

			# apply unit scaling if necessary
			unit.scaling <- match.arg(unit.scaling)
			if( unit.scaling == 'column' )
				object <- apply(object, 2, function(x) x/sum(x))
			else if ( unit.scaling == 'row' )
				object <- t(apply(object, 1, function(x) x/sum(x)))

			# check validity of palette
			col.palette <- brewer.pal(brewer.pal.info[palette,'maxcolors'],palette)
			if( rev.palette ) col.palette <- rev(col.palette) 
			
			# set default graphical parameters (if those are not already set)
			graphical.params <- .set.list.defaults(graphical.params
					, cexRow=0.8, cexCol=0.8
					, hclustfun = function(m) hclust(m,method="average")
					, dendrogram='none'
					, col=col.palette
					, scale='none', trace="none"
					, keysize=1, margins=c(5,10)
			)
			
			# if a known class is provided, add a side color over the top row
			if( !missing(class) ){
				if( !is.factor(class) ) class <- as.factor(class)
				class.num <- as.numeric(class)
				legend.pal <- palette(rainbow(max(2,nlevels(class))))[1:nlevels(class)]
				col.matrix <- matrix(legend.pal[class.num], ncol(object), 1)
				
				# show association with metagenes
				if( show.prediction ){
					# only if there is less than 9 metagenes 
					# cf. limitation of brewer color palette
					if( nrow(object) <= 9 ){
						prediction <- .predict.nmf(object)
						prediction.num <- as.numeric(prediction)
						pal.pred <- brewer.pal(max(3,nrow(object)),'Set2')[1:nrow(object)]
						col.matrix <- cbind(pal.pred[prediction.num], col.matrix)
						graphical.params <- .set.list.defaults(graphical.params
								, RowSideColors=pal.pred
						)
				}
					else warning("NMF::metaHeatmap - cannot not show prediction for more than 9 metagenes.")
				}
				# do that otherwise heatmap.plus complains
				if( ncol(col.matrix) < 2 )
					col.matrix <- cbind(col.matrix, col.matrix)
				
				# add the ColSideColors
				graphical.params <- .set.list.defaults(graphical.params
								, ColSideColors=col.matrix
						)
			}
			
			
			res.heatmap <- do.call('heatmap.plus.2', c(list(object), graphical.params))
			
			if( !missing(class) ){
				# order properly the legend boxes
				class.num <- as.numeric(class[res.heatmap$colInd])
			
				occ <- NA # will store the current number of occurences
				class.max.occ <- rep(0, nlevels(class)) # will store the current maximum number of occurences per class
				class.start <- rep(NA, nlevels(class)) # will store the current start of the longer stretch per class
				last.l <- ''
				sapply( seq(length(class.num), 1, -1), 
						function(i){
							l <- class.num[i]
							if(l==last.l){
								occ <<- occ + 1
							}else{
								occ <<- 1
							}
							if(occ > class.max.occ[l]){
								class.max.occ[l] <<- occ
								class.start[l] <<- i
							}
							last.l <<- l
						}
				)
				
				class.ord <- order(class.start)
				l.names <- levels(class)[class.ord]
				l.color <- legend.pal[class.ord]
				legend('top', title='Classes'
						, legend=l.names, fill=l.color
						, horiz=TRUE, bty='n')
			}
			
			# return invisible
			invisible(res.heatmap)
		}
)

setMethod('metaHeatmap', signature(object='NMF'),
	function(object, what=c('samples', 'features'), filter=FALSE, ...){

		what <- match.arg(what)
		if( what == 'samples' ){
			# send deprecated warning
			.Deprecated('coefmap', 'NMF', "Direct use of the S4-Method 'metaHeatmap' for 'NMF' objects is deprecated, use 'coefmap' instead.")
			
			# call the new function 'coefmap'
			return( coefmap(object, ...) )			
			
		}else if( what == 'features' ){
			# send deprecated warning
			.Deprecated('basismap', 'NMF', "Direct use of the S4-Method 'metaHeatmap' for 'NMF' objects is deprecated, use 'basismap' instead.")
			
			# call the new function 'basismap'
			return( basismap(object, filter=filter, ...) )
			
		}
	}
)

###% Heatmap of the basis vector matrix
setGeneric('basismap', function(object, ...) standardGeneric('basismap') )
setMethod('basismap', signature(object='NMF'),
	function(object, subsetRow=FALSE, tracks = 'basis', info = FALSE, ...){
	
		# retreive the graphical parameters and match them to the sub-sequent call to 'heatmap.plus.2'
		graphical.params <- list(...)
		names(graphical.params) <- .match.call.args(names(graphical.params), 'aheatmap', call='NMF::basismap')
			
		# for backward compatibility: look for argument filter
		if( !is.null(graphical.params$filter) ){
			if( missing(subsetRow) ){
				warning("NMF::basismap - Argument `filter` is deprecated and was renamed to `subsetRow`.")
				subsetRow <- graphical.params$filter
			}else
				warning("NMF::basismap - Discarding deprecated argument `filter`.")
			graphical.params$filter <- NULL
		}
						
		# resolve subsetRow if its a single value
		if( length(subsetRow) == 1 ){
			subsetRow <- 
			if( identical(subsetRow, FALSE) )
				NULL
			else if( isTRUE(subsetRow) ) # use Kim and Park scoring scheme for filtering 			
				extractFeatures(object, format='combine')
			else if( is.character(subsetRow) ) # use subsetRow as a filtering method
				extractFeatures(object, method=subsetRow, format='combine')
			else if( is.numeric(subsetRow) ){ # numerical value for filter
			
				# only keep the specified number of feature for each basis vector
				ord <- apply( basis(object), 2, 
						function(v, limit){
							order(v, decreasing=TRUE)[1:limit]
						}
						, as.integer(subsetRow))
				
				unique(as.integer(ord))
			}
			else stop("NMF::basismap - invalid single value for argument 'subsetRow' [logical, numeric or character expected]")
			
		}
		
		# extract the basis vector matrix
		x <- basis(object)
			
		# add side information if requested
		info <- if( isTRUE(info) && isNMFfit(object) ) paste("Method: ", algorithm(object))
				else if( isFALSE(info) ) NULL
				else info
		
		# set default graphical parameters for type 'basis'
		graphical.params <- .set.list.defaults(graphical.params
				, main="Basis components"				
				, Colv=NA
				, color = 'YlOrRd'
				# custom parameters
				, scale = 'r1'
				, info = info
		)
		
		# add annotation tracks		
		if( length(tracks) >0 && !isNA(tracks) ){
			
			# create extra annotation tracks			
			tr <- sapply( tracks, function(t){
				switch(t
						, basis = predict(object, 'features')					
						, stop("NMF::basismap - Invalid annotation track: ", t)
				)
			}, simplify=FALSE)
			
			graphical.params[['annRow']] <- atrack(tr, graphical.params[['annRow']])			
		}
		
		# call metaHeatmap on matrix
		do.call('aheatmap', c(list(x, subsetRow=subsetRow), graphical.params))	
	}
)

###% Heatmap of the mixture coefficient matrix
setGeneric('coefmap', function(object, ...) standardGeneric('coefmap') )
setMethod('coefmap', signature(object='NMF'),
	function(object, tracks = 'basis', info = FALSE, ...){
	
		# retreive the graphical parameters and match them to the sub-sequent call to 'heatmap.plus.2'
		graphical.params <- list(...)	
		names(graphical.params) <- .match.call.args(names(graphical.params), 'aheatmap', call='NMF::coefmap')
				
		# use the mixture coefficient matrix
		x <- coef(object)
		
		
		# add side information if requested
		info <- if( isTRUE(info) && isNMFfit(object) ) paste("Method: ", algorithm(object))
				else if( isFALSE(info) ) NULL
				else info
		
		# set default graphical parameters for type 'coefficients'
		graphical.params <- .set.list.defaults(graphical.params
				, main="Mixture coefficients"
				, Rowv = NA
				# default ordering: order the columns by their dominant basis
				, Colv=order(as.numeric(predict(object))) 
				, color = 'YlOrRd'				
				, scale = 'c1'
				, info = info
		)
		
		# add annotation tracks
		if( length(tracks) >0 && !isNA(tracks) ){
			
			# create extra annotation tracks
			tr <- sapply( tracks, function(t){
				switch(t
					, basis = predict(object)					
					, stop("NMF::coefmap - Invalid annotation track: ", t)
			)
			}, simplify=FALSE)
			
			graphical.params[['annCol']] <- atrack(tr, graphical.params[['annCol']])
		}
		
		
		# call metaHeatmap on matrix
		do.call('aheatmap', c(list(x), graphical.params))
	}
)

#setGeneric('hist', package='graphics')
#setMethod('hist', signature(x='NMF'), 
#	function(x, ref=1, alpha=20, ...)
#	{
#		stopifnot( ref >= 1 && ref <= ncol(metagenes(x)) )
#		alpha = sprintf("%02d", alpha) #add leading zero to alpha if nessecary
#		
#		# save graphical parameters
#		oldpar = par(no.readonly=TRUE)
#		on.exit( {par(oldpar)} ) # reset the graphical parameters on exit
#
#		# order genes by decreasing contribution to the reference factor
#		M = metagenes(x)[order(metagenes(x)[,ref], decreasing=T), ] 
#		
#		#plot the contributions to the reference factor
#		par(lwd = 0.5)
#		x = seq(nrow(M))
#		html.colors = apply( col2rgb( seq(ncol(M))+1 ), 2, function(x) paste("#", paste(intToHex(x), collapse=''), alpha, sep='') )
#		plot(x=x, y=M[,ref], type='h'
#			, col=html.colors[ref], ylim=c(min(M), max(M))
#			, main='Contribution to metagenes', xlab=paste('Genes ordered based on factor', ref), ylab='Contribution')
#			
#		# plot the remaining metagenes		
#		remaining.factor = seq(ncol(M))[seq(ncol(M)) != ref]						
#		sapply(remaining.factor, 
#			function(f){				
#				lines(x=x, M[,f], type='h', col=html.colors[f])
#			}
#		)
#		
#		#put the legend
#		legend('top', legend=paste('Factor', seq(ncol(M))), fill=sub("^(#[a-f0-9]{6}).*", "\\1", html.colors, ignore.case=TRUE) )
#		
#		invisible()
#	}
#)

###% Utility function used to sets default elements in a list if they are
###% not already set
###% The default values are given in argument ...
.set.list.defaults <- function(input.list, ...){
	defaults <- list(...)
	n <- names(defaults)[!names(defaults) %in% names(input.list)]
	lapply(n, function(name){
			input.list[[name]] <<- defaults[[name]]
		}
	)
	input.list
}

###% Partially match arguments for a given function 
.match.call.args <- function(x, fun, in.fun=NULL, call=NULL){
	stopifnot( is.character(fun) && length(fun) == 1 )
	if( length(x) == 0 ) return(x)
	x.ind <- charmatch(x, args <- formalArgs(getFunction(fun)))
	
	sapply(seq(length(x)), function(i){
			ind <- x.ind[i]
			# the argument is not part of the call: keep it unchanged
			if( is.na(ind) ) return(x[i])
			# multiple matches: error
			if( ind == 0 ){
				alt <- paste(grep(paste('^', x[i], sep=''), args, value=TRUE), collapse=', ')
				stop(if( !is.null(call) ) c(call, ' - '), "Multiple match for argument '", x[i], "' of function '"
					, if( is.null(in.fun) ) fun else in.fun, "' [use one of: ", alt, "]"
					, call.=FALSE)
			}
			# return the matched full names
			args[ind]
	})
}


###% Computes a set of measures usefull to assess the factorization's quality.
###% 
###% 
###% @param object a \code{NMF} object
###% @return a numeric vector of the measures. 
###% 
setGeneric('summary', package='base')
setMethod('summary', signature(object='NMF'), 
		function(object, class, target){
			
			res <- numeric()
			
			## IMPORTANT: if adding a summary measure also add it in the sorting 
			## schema of method NMFSet::compare to allow ordering on it
			
			# rank
			res <- c(res, rank=nbasis(object))
			# compute sparseness
			res <- c(res, sparseness=sparseness(object))
			
			# if class is provided: also computes entropy and purity
			if( !missing(class) ){
				# compute purity
				res <- c(res, purity=purity(object, class=class))
				# compute entropy
				res <- c(res, entropy=entropy(object, class=class))
			}
			
			# if the target is provided compute the RSS
			if( !missing(target) ){
				RSS <- rss(object, target)
				res <- c(res, rss=RSS)
				# explained variance
				res <- c(res, evar=evar(object, target))
			}
			
			# return result
			return(res)
		}
)

###% Computes the sparseness of a vector as defined by Hoyer (2004).
###%
###% This sparseness measure quantifies how much energy of a vector is packed into only few components.
###% It is defined by:
###% \deqn{Sparseness(x) = \frac{\sqrt{n} - \frac{\sum |x_i|}{\sqrt{\sum x_i^2}}}{\sqrt{n}-1}
###%, where \eqn{n} is the length of \code{x}.
###%
###% The sparseness is a real number in \eqn{[0,1]}. It is equal to 1 if and only if \code{x} contains 
###% a single nonzero component, and is equal to 0 if and only if all components of \code{x} are equal.
###% It interpolates smoothly between these two extreme values. The closer to 1 is the spraseness the sparser is the vector.
###%
###% @param x a vector.
###% @param ... extra parameters
###% @return the sparseness of \code{x} when \code{x} is a vector
###% , the mean sparseness of the columns of \code{x} when \code{x} is a matrix,
###% , the mean sparseness of the basis (resp. encoding vectors (i.e. rows of the mixture coeffcients matrix))
###% when \code{x} is a \code{NMF} object and \code{what} is \code{'genes'} (resp. \code{what} is \code{'samples'}).
###%
###% @references Hoyer, P. O. (2004)
###% , 'Non-negative Matrix Factorization with Sparseness Constraints'
###% , Journal of Machine Learning Research 5 (2004) 1457–1469
###%
setGeneric('sparseness', function(x, ...) standardGeneric('sparseness') )
setMethod('sparseness', signature(x='numeric'), 
	function(x){
		# get length of x
		n <- length(x)
		# compute and return the sparseness
		( sqrt(n) - sum(abs(x)) / sqrt(sum(x^2)) ) / (sqrt(n)-1)
	}
)
setMethod('sparseness', signature(x='matrix'), 
	function(x){
		# compute the sparseness of each column
		s <- apply(x, 2, sparseness)
		
		# return the mean sparseness
		mean(s)
	}
)
setMethod('sparseness', signature(x='NMF'), 
	function(x){		
		# return the sparseness of the basis and coef matrix
		c(basis=sparseness(basis(x)), coef=sparseness(coef(x)))
	}
)


###% Computes the purity of a clustering given a known factor.
###%
###% The purity is a measure of performance of a clustering method, in recovering the classes defined by a factor a priori known (i.e. one knows the true class labels).
###% Suppose we are given \eqn{l} categories, while the clustering method generates \eqn{k} clusters. Purity is given by:
###% \deqn{Purity = \frac{1}{n} \sum_{q=1}^k \max_{1 \leq j \leq l} n_q^j}
###%, where:
###% - \eqn{n} is the total number of samples;
###% - \eqn{n_q^j} is the number of samples in cluster \eqn{q} that belongs to original class \eqn{j} (\eqn{1 \leq j \leq l}).
###%
###% The purity is therefore a real number in \eqn{[0,1]}. The larger the purity, the better the clustering performance.
###%
###% @param x a factor or integer vector that gives the cluster membership for each sample.
###% @param class a factor or integer vector that gives the true class labels for each sample.
###% @return the purity (i.e. a numerical value)
###%
###% @references Kim, H. & Park, H. 
###%	Sparse non-negative matrix factorizations via alternating non-negativity-constrained least squares for microarray data analysis.
###%	Bioinformatics (2007). 
###%	\url{http://dx.doi.org/10.1093/bioinformatics/btm134}.
###%	
setGeneric('purity', function(x, class, ...) standardGeneric('purity') )
setMethod('purity', signature(x='table', class='missing'), 
	function(x, class, ...){
		#for each cluster: compute maximum number of samples common to a class
		t <- apply(x, 1, max)
		# average and return the result
		sum(t) / sum(x)
	}
)
setMethod('purity', signature(x='factor', class='factor'), 
	function(x, class, ...){		
		#compute the purity on the contingency table between clusters and true classes (clusters are in rows)
		purity(table(x, class))
	}
)

setMethod('purity', signature(x='NMF', class='factor'), 
	function(x, class, ...){
		# compute the purity for the samples clusters defined by the profiles
		purity(predict(x, what='samples'), class)
	}
)
setMethod('purity', signature(x='NMF', class='ANY'), 
		function(x, class, ...){
			# convert class into a factor
			class <- as.factor(class)
			# compute the purity
			purity(x, class)
		}
)


###% Computes the entropy of a clustering given a known factor.
###%
###% The entropy is a measure of performance of a clustering method, in recovering classes defined by factor a priori known (i.e. one knows the true class labels).
###% Suppose we are given \eqn{l} categories, while the clustering method generates \eqn{k} clusters. Entropy is given by:
###% \deqn{Entropy = - \frac{1}{n \log_2 l} \sum_{q=1}^k \sum_{j=1}^l n_q^j \log_2 \frac{n_q^j}{n_q}
###%, where:
###% - \eqn{n} is the total number of samples;
###% - \eqn{n} is the total number of samples in cluster \eqn{q};
###% - \eqn{n_q^j} is the number of samples in cluster \eqn{q} that belongs to original class \eqn{j} (\eqn{1 \leq j \leq l}).
###%
###% The smaller the entropy, the better the clustering performance.
###%
###% @param x a factor or integer vector that gives the cluster membership for each sample.
###% @param class a factor or integer vector that gives the true class labels for each sample.
###% @return the entropy (i.e. a numerical value)
###%
###% @references Kim, H. & Park, H. 
###%	Sparse non-negative matrix factorizations via alternating non-negativity-constrained least squares for microarray data analysis.
###%	Bioinformatics (2007). 
###%	\url{http://dx.doi.org/10.1093/bioinformatics/btm134}.
###%	
setGeneric('entropy', function(x, class, ...) standardGeneric('entropy') )
setMethod('entropy', signature(x='table', class='missing'), 
	function(x, class, ...){
		#for each cluster: compute the inner sum
		t <- apply(x, 1, function(n){ c.size <- sum(n); n %*% ifelse( n!=0, log2(n/c.size), 0)} )
		
		# weight and return the result
		- sum(t) / ( sum(x) * log2(ncol(x)) )
	}
)

setMethod('entropy', signature(x='factor', class='factor'), 
	function(x, class, ...){
		#copmute entropy on contingency table between clusters and true classes (clusters are in rows)
		entropy(table(x, class))
	}
)

setMethod('entropy', signature(x='NMF', class='factor'), 
	function(x, class, ...){
		# compute the entropy for the samples clusters defined by the metagenes expression matrix
		entropy(predict(x, what='samples'), class)
	}
)

setMethod('entropy', signature(x='NMF', class='ANY'), 
		function(x, class, ...){
			# convert class into a factor
			class <- as.factor(class)
			# compute the entropy
			entropy(x, class)
		}
)

###% Extract the genes that characterize each factor.
###%
###% For each factor the genes are first sorted by decreasing contribution. The first successive ones whose contribution to the factor
###% is greater than their contribution to all other metagenes are selected.
###%
###% @param x the matrix of metagenes. That is a matrix with metagenes in column, genes in row, contain the genes' contribution to each factor
###% @return a list with number of metagenes elements, each being a vector containing the indexes of the characterizing genes
#setGeneric('computeContrib', function(x, ...) standardGeneric('computeContrib') )
#setMethod('computeContrib', signature(x='matrix'), 
#	function(x, ...){
#		# determine the specific genes for each factor
#		lapply(1:ncol(x), 
#		function(i){		
#			g <- x[,i]		
#			#order by decreasing contribution to factor i
#			index.sort <- order(g, decreasing=TRUE)		
#			
#			for( k in seq_along(index.sort) )
#			{
#				index <- index.sort[k]
#				#if the gene contributes more to any other factor then return the genes above it
#				if( any(x[index,-i] >= g[index]) )
#				{
#					if( k == 1 ) return(NULL)				
#					else return(rownames(x)[index.sort[1:(k-1)]])
#				}
#			}
#			
#			# all genes are meeting the criteria
#			rownames(x)
#		})			
#	}
#)
#

###% Computes the genes scores as suggested in Kim, H. & Park, H (Bioinformatics 2007).
###%
###% The score for gene \eqn{i} is defined as:
###% \deqn{S_i = 1 + \frac{1}{\log_2 k} \sum_{q=1}^k p(i,q) \log_2 p(i,q),}
###% where \eqn{p(i,q)} is the probability that the \eqn{i}-th gene contributes to cluster \eqn{q}:
###% \deqn{p(i,q) = \frac{W(i,q)}{\sum_{r=1}^k W(i,r)} }
###%
###% The gene score is a real value within the range [0,1]. The higher the gene score the more cluster-specific the corresponding gene.
###%
###% @references Kim, H. & Park, H. 
###%	Sparse non-negative matrix factorizations via alternating non-negativity-constrained least squares for microarray data analysis.
###%	Bioinformatics (2007). 
###%	\url{http://dx.doi.org/10.1093/bioinformatics/btm134}.
###%	
setGeneric('featureScore', function(object, ...) standardGeneric('featureScore') )
setMethod('featureScore', signature(object='matrix'), 
	function(object, method=c('kim', 'max')){
		
		method <- match.arg(method)
		score <- switch(method,
		
		kim = {
			#for each row compute the score
			s <- apply(object, 1, function(g){
					g <- abs(g)
					p_i <- g/sum(g)
					crossprod(p_i, log2(p_i))
			})
			# scale, translate and return the result
			1 + s / log2(ncol(object))		
			}
		, max = {
			apply(object, 1, function(x) max(abs(x)))
			}
		)
		
		# return the computed score
		return(score)
	}
)

setMethod('featureScore', signature(object='NMF'), 
	function(object, method=c('kim', 'max')){
		featureScore(basis(object), method)
	}
)


###% Identify the most metagenes-specific genes as suggested in Kim, H. & Park, H (Bioinformatics 2007).
###%
###% The genes are first scored using the function \code{featureScore}. Then only the genes that fullfil both following criteria are retained:
###% - score greater than \eqn{\hat{\mu} + 3 \hat{\sigma}}, where \eqn{\hat{\mu}} and \eqn{\hat{\sigma}} are the median and the median absolute deviation (MAD) of the scores respectively;
###% - the maximum contribution to a metagene is greater than the median of all contributions (i.e. of all elements of W)
###%
###% @references Kim, H. & Park, H. 
###%	Sparse non-negative matrix factorizations via alternating non-negativity-constrained least squares for microarray data analysis.
###%	Bioinformatics (2007). 
###%	\url{http://dx.doi.org/10.1093/bioinformatics/btm134}.
###%	
setGeneric('extractFeatures', function(object, ...) standardGeneric('extractFeatures') )
setMethod('extractFeatures', signature(object='NMF'), 
	function(object, method=c('kim', 'max'), format=c('list', 'combine', 'subset')){
		
		res <- method <- match.arg(method)
		res <- switch(method,
			kim = { # KIM & PARK method
				
				# first score the genes
				s <- featureScore(object, method='kim')
				
				# filter for the genes whose score is greater than \mu + 3 \sigma
				th <- median(s) + 3 * mad(s)
				sel <- s >= th
				#print( s[sel] )
				#print(sum(sel))
				
				# build a matrix with:
				#-> row#1=max column index, row#2=max value in row, row#3=row index
				temp <- 0;
				g.mx <- nmfApply(object, 1L, 
						function(x){
							temp <<- temp +1
							i <- which.max(abs(x));
							#i <- sample(c(1,2), 1)
							c(i, x[i], temp)
						}
				)
				
				# test the second criteria
				med <- median(abs(basis(object)))
				sel2 <- g.mx[2,] >= med
				#print(sum(sel2))
				
				# subset the indices
				g.mx <- g.mx[, sel & sel2, drop=FALSE] 
				
				# return the indexes of the features that fullfil both criteria
				#Note: make sure there is an element per basis
				cl <- factor(g.mx[1,], levels=seq(nbasis(object))) 
				res <- split(g.mx[3,], cl)
				res <- lapply(res, function(ind){ if(length(ind)==0) ind<-NA; as.integer(ind)} )
				
				# add the threshold used
				attr(res, 'threshold') <- th
				
				# return result
				res
				
			},
			max = { # MAX method from bioNMF
				
				# determine the specific genes for each basis vector
				res <- lapply(1:nbasis(object), 
					function(i){
						mat <- basis(object)
						vect <- mat[,i]
						#order by decreasing contribution to factor i
						index.sort <- order(vect, decreasing=TRUE)		
						
						for( k in seq_along(index.sort) )
						{
							index <- index.sort[k]
							#if the feature contributes more to any other factor then return the features above it
							if( any(mat[index,-i] >= vect[index]) )
							{
								if( k == 1 ) return(as.integer(NA))
								else return( index.sort[1:(k-1)] )
							}
						}
						
						# all features meet the criteria
						seq_along(vect)
					}
				)
				
				# return res
				res
			}
		)
		
		# apply the desired format
		format <- match.arg(format)
		res <- switch(format
				#combine: return all the indices in a single vector
				, combine = res <- unique(unlist(res, use.names=FALSE)) 
				#subset: return the object subset with the selected indices
				, subset = {ind <- unique(unlist(res, use.names=FALSE));  object[ind,] }
				#else: leave as a list
				, res
				)
		
		# return result
		return( res )
	}
)

###% 'apply' method for NMF objects
###% When MARGIN = 1, the apply is performed on the rows of the basis matrix
###% When MARGIN = 2, the apply is performed on the columns of the mixture coefficient matrix
###% When MARGIN = 3, the apply is performed on the columns of the basis coefficient matrix (i.e. the basis vector)
###% When MARGIN = 4, the apply is performed on the rows of the mixture coefficient matrix (i.e. the basis profiles)
setGeneric('nmfApply', function(object, ...) standardGeneric('nmfApply') )
setMethod('nmfApply', signature(object='NMF'),
	function(object, MARGIN, FUN, ...){
		if( MARGIN == 1 )
			apply(basis(object), 1, FUN, ...)
		else if( MARGIN == 3 )
			apply(basis(object), 2, FUN, ...)
		else if( MARGIN == 2 )
			apply(coef(object), 2, FUN, ...)
		else if( MARGIN == 4 )
			apply(coef(object), 1, FUN, ...)
		else stop("NMF::nmfApply : invalid argument 'MARGIN' (expected value is: 1-basis rows, 2-coef columns, 3-basis columns, or 4-coef rows)")
	}
)

###% Utility function to compute the dominant column for each row for a matrix.
.predict.nmf <- function(x, prob=FALSE){
	
	if( !is.matrix(x) ) stop('NMF:::.predict.nmf : only works on matrices')
	if( !prob ){
		#for each column return the (row) index of the maximum
		return( as.factor(apply(x, 2, function(v) which.max(abs(v)))) )
	}
	else{
		#for each column return the (row) index of the maximum AND the associated probaility
		res <- apply(x, 2,
				function(p){
					p <- abs(p)
					i <- which.max(p)
					c(i, p[i]/sum(p))
				}
		)
		# return the result as a list of two elements
		return( list(predict=as.factor(res[1,]), prob=res[2,]) )
	}
}

###% \code{predict} returns the predicted cluster for each row/column
###%
###% @param object an NMF object
###% @param what a character string specifying which clusters should be returned:
###% \code{'rows'} and \code{'features'} return row clusters, \code{'columns'} 
###% and \code{'samples'} columns clusters.
###% The clusters defined by the index of the dominant column (resp. row) of the 
###% basis matrix (resp. mixture coefficient matrix). 
###% In the case of \code{\link{NMFfitX}} objects (resulting from multiple NMF runs), 
###% then \code{what} can also be \code{'consensus'} to return the clusters as defined 
###% by the hierarchical clustering of the consensus matrix, using the euclidean 
###% distance and average linkage. These correspond to column clusters 
###% (i.e. groups of samples) and might be different from the column clusters 
###% returned by the best fit of the multiple runs.  
###% 
###% @return a factor of length the number of columns (resp. rows)
###% 
###% @export
###% @rdname stats
###% 
setGeneric('predict', package='stats')
setMethod('predict', signature(object='NMF'),
		function(object, what=c('columns', 'rows', 'samples', 'features'), prob=FALSE){
			# determine which matrix to use for the prediction
			what <- match.arg(what)
			x <- if( what %in% c('features', 'rows') ) t(basis(object)) else coef(object)
			
			# compute the indice of the dominant row for each column
			return( .predict.nmf(x, prob) )
		}
)

####% Compute the dominant column for each row.
####%
####% @param x a matrix containing the mixture coefficients (basis vector in rows, samples in columns)
####% @return a factor of length the number of columns, giving the dominant column for each row
####% @note This function is now deprecated
#setGeneric('clusters', function(object, newdata, ...) standardGeneric('clusters') )
####% Compute the dominant metagene for each sample.
####%
####% @param x a NMF object
####% @return a factor of length the number of samples, giving the dominant metagene for each sample
####% @note This function is now deprecated
#setMethod('clusters', signature(object='NMF', newdata='missing'), 
#	function(object, newdata, ...){
#		predict(object, ...)
#	}
#)

###% Computes the Correlation Between NMF Basis Matrices
setGeneric('basiscor', function(x, y, ...) standardGeneric('basiscor') )
setMethod('basiscor', signature(x='NMF', y='matrix'),
	function(x, y, ...){
		cor(basis(x), y, ...)
	}
)
setMethod('basiscor', signature(x='NMF', y='NMF'),
		function(x, y, ...){
			basiscor(x, basis(y), ...)
		}
)
setMethod('basiscor', signature(x='NMF', y='missing'),
		function(x, y, ...){
			basiscor(x, x, ...)
		}
)
###% Method 'basiscor' for swapped arguments 
setMethod('basiscor', signature(x='matrix', y='NMF'),
		function(x, y, ...){
			cor(x, basis(y), ...)
		}
)

###% Computes the Correlation Between NMF Mixture Coefficient Matrices
setGeneric('profcor', function(x, y, ...) standardGeneric('profcor') )
setMethod('profcor', signature(x='NMF', y='matrix'),
		function(x, y, ...){
			cor(t(coef(x)), t(y), ...)
		}
)
setMethod('profcor', signature(x='NMF', y='NMF'),
		function(x, y, ...){
			profcor(x, coef(y), ...)
		}
)
setMethod('profcor', signature(x='NMF', y='missing'),
		function(x, y, ...){
			profcor(x, x, ...)
		}
)
###% Method 'profcor' for swapped arguments
setMethod('profcor', signature(x='matrix', y='NMF'),
		function(x, y, ...){
			cor(t(x), t(coef(y)), ...)
		}
)


###% Computes the connectivity matrix for the samples based on their metagenes expression.
###%
###% The connectivity matrix of a clustering is a matrix \eqn{C} containing only 0 or 1 entries such that:
###% \deqn{C_{ij} = \left\{\begin{array}{l}1\mbox{ if sample }i\mbox{ belongs to the same cluster as sample }j\\0\mobx{ otherwise}\end{array}\right..}
###%
setGeneric('connectivity', function(x, ...) standardGeneric('connectivity') )
setMethod('connectivity', signature(x='NMF'), 
	function(x, ...){
		c <- predict(x, what='samples');
		C <- outer(c, c, function(x,y) ifelse(x==y, 1,0));
		class(C) <- c(class(C), 'NMF.consensus')
		attr(C, 'nrun') <- 1
		attr(C, 'nbasis') <- nbasis(x)
		C
	}
)

setGeneric('rss', function(object, ...) standardGeneric('rss'))
setMethod('rss', 'NMF', 
	function(object, target){
		
		# make sure the target is provided
		if( missing(target) ) stop("NMF::rss - Argument 'target' is missing and required to compute the residual sum of squares.")
		
		# use the expression matrix if necessary
		if( inherits(target, 'ExpressionSet') ){
			# requires Biobase
			if( !require.quiet(Biobase) ) 
				stop("NMF::rss - The 'Biobase' package is required to extract expression data from 'ExpressionSet' objects [see ?'nmf-bioc']")
			
			target <- Biobase::exprs(target)
		}else if( is.data.frame(target) )
			target <- as.matrix(target)
		
		# return rss using the optimized C function
		.rss(fitted(object),target)
	}
)

setGeneric('evar', function(object, ...) standardGeneric('evar'))
setMethod('evar', 'NMF', 
		function(object, target){
			
			# make sure the target is provided
			if( missing(target) ) stop("NMF::evar - Argument 'target' is missing and required to compute the explained variance.")
			
			# use the expression matrix if necessary
			if( inherits(target, 'ExpressionSet') ){
				# requires Biobase
				if( !require.quiet(Biobase) ) 
					stop("NMF::evar - The 'Biobase' package is required to extract expression data from 'ExpressionSet' objects [see ?'nmf-bioc']")
				
				target <- Biobase::exprs(target)
			}
			
			t <- as.numeric(target)
			1 - rss(object, target) / sum(t^2)
		}
)


###% Common interface to compute matrix distances registered in the NMF registry.
###% 
###% This method provides a single interface to compute the error between a target matrix and its estimate. The 
###% distance method can either be a registered NMF-distance (cf TODO: registry distance) or any defined \code{function} (cf. TODO: details for distance function). 
###% 
###% @author Renaud Gaujoux \email{renaud@@cbio.uct.ac.za}
###% @export
setGeneric('distance', function(target, x,...) standardGeneric('distance') )
###% Computes the distance between a target matrix and its estimate given by a \code{NMF} object.
###% 
###% ###% The target matrix \code{target} and its NMF estimate \code{x} must have the same dimension.
###% 
###% @param target the target \code{matrix} with the same dimension as \code{x}.  
###% @param x a \code{NMF} object with the same dimension as \code{target}.
###% @param method if missing, slot \code{distance} of \code{x} is used.
###% @param ... extra 
###% 
###% @return the distance (a nonnegative numerical value) between matrices \code{target} and \code{fitted(x)}. 
###% @author Renaud Gaujoux \email{renaud@@cbio.uct.ac.za}
###% @export
###% 

setMethod('distance', signature(target='matrix', x='NMF'), 
	function(target, x, method=c('', 'KL', 'euclidean'), ...){

	fun <- distance(method=method)
	
	if( is.null(fun) ){
		warning('Undefined distance method: distance cannot be computed [returned NA]')
		return(as.numeric(NA))
	}
	
	# apply the function and return the result
	fun(target, x, ...)

	}
)

setMethod('distance', signature(target='missing', x='missing'), 
	function(target, x, method=c('', 'KL', 'euclidean')){
	
		#message('compute distance')
		# determinate the distance measure to use
		if( is.null(method) ) return(NULL)
		
		if( is.character(method) ){
			errMeth <- try(method <- match.arg(method), silent=TRUE)
			# if the method is not predefined, try to find a function with the given name
			if( inherits(errMeth, 'try-error') ){			
				#TODO: this is not working with local functions
				if( is.character(method) ){
					errFun <- try(fun <- getFunction(method), silent=TRUE)
					if( inherits(errFun, 'try-error') ) stop("Could not find distance measure '", method, "':\n\t- not a predefined measures -> ", errMeth,"\t- not a function -> ", errFun)
				}
				else fun <- method
				
				if( !is.function(fun) )
					stop('Invalid distance measure: should be a character string or a valid function definition')
			}
			else{
				# compute and return the distance measure		
				fun <- switch(method,
						euclidean = function(target, x, ...){
							# call optimized C function
							.rss(target, fitted(x))/2							
						},
						KL = function(target, x, ...){							
							# call optimized C function
							.KL(target, fitted(x))					
						}
				)
			}
		}
		else if( is.function(method) )
			fun <- method
		else
			stop('Invalid distance measure: should be a character string or a valid function definition')
		
		# return the distance function
		fun
		
	}
		
)


###% Compare two NMF models
setGeneric('nmf.equal', function(x, y, ...) standardGeneric('nmf.equal') )
setMethod('nmf.equal', signature(x='NMF', y='NMF'), 
		function(x, y, identical=TRUE, ...){
			if( identical )
				identical(x, y)
			else
				all.equal(x, y, ...)
		}
)

# Match and Order Basis Components
#
# 
match.basis <- function(object, return.table=FALSE){
	
	# compute the contingency table
	#pcmap <- predict(object, 'cmap')
	# build the tree from consensus matrix
	h <- hclust(as.dist(1-consensus(object)), method='average')
	# extract membership from the tree
	cl <- cutree(h, k=nbasis(object))
	# change the class indexed to match the order of the consensus clusters 
	cl <- match(cl, unique(cl[h$order]))
	pcmap <- as.factor(cl)
	occ <- table(consensus=pcmap, fit=predict(object))		
	
	# add names if present
#	if( !is.null(basisnames(object)) ){
#		rownames(occ) <- colnames(occ) <- basisnames(object)
#	}
		
	# for each estimated component look for the maximum agreement		
	T.tmp <- occ				
	res <- rep(0, ncol(T.tmp))	
	for( i in 1:ncol(T.tmp) ){
		# get the row and column index of the maximum over the remaining entries 
		xm <- which.max(T.tmp)-1
		jm <- xm %/% nrow(T.tmp) + 1
		im <- xm - (jm-1) * nrow(T.tmp) + 1
		
		# assign the estimate row to the inferred reference column
		stopifnot( res[im]==0 )
		res[im] <- jm
		
		# erase the assigned estimate row
		T.tmp[im,] <- NA
		# erase the assigned reference column
		T.tmp[,jm] <- NA		
	}
	
	# return the mapping as an integer vector
	res <- as.integer(res)
	if( return.table )
		res <- list(match=res, table=occ)
	
	# return result
	res
}