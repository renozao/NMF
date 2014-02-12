#library(R.utils)

#' @include utils.R
#' @include versions.R
#' @include algorithmic.R
#' @include aheatmap.R
NULL

#' Advanced Usage of the Package NMF
#' 
#' The functions documented here provide advanced functionalities useful when
#' developing within the framework implemented in the NMF package.
#'
#' @rdname advanced
#' @name advanced-NMF
NULL 

# declare old S3 class 'proc_time' to use it as a slot for class NMF 
setOldClass('proc_time', prototype=numeric())

################################
# Class: NMF
################################

#' Generic Interface for Nonnegative Matrix Factorisation Models
#' 
#' The class \code{NMF} is a \emph{virtual class} that defines a common 
#' interface to handle Nonnegative Matrix Factorization models (NMF models) 
#' in a generic way.
#' Provided a minimum set of generic methods is implemented by concrete 
#' model classes, these benefit from a whole set of functions and utilities
#' to perform common computations and tasks in the context of Nonnegative Matrix
#' Factorization.
#' 
#' Class \code{NMF} makes it easy to develop new models that integrate well
#' into the general framework implemented by the \emph{NMF} package.
#' 
#' Following a few simple guidelines, new types of NMF models benefit from all the
#' functionalities available for the built-in NMF models -- that derive themselves
#' from class \code{NMF}.
#' See section \emph{Implementing NMF models} below.
#' 
#' See \code{\linkS4class{NMFstd}}, and references and links therein for
#' details on the built-in implementations of the standard NMF model and its 
#' extensions.
#' 
#' @slot misc A list that is used internally to temporarily store algorithm 
#' parameters during the computation.
#' 
#' @export
#' @family NMF-interface 
#' 
#' @section Implementing NMF models:
#' 
#' The class \code{NMF} only defines a basic data/low-level interface for NMF models, as 
#' a collection of generic methods, responsible with data handling, upon which 
#' relies a comprehensive set of functions, composing a rich higher-level interface.
#' 
#' Actual NMF models are defined as sub-classes that inherits from class 
#' \code{NMF}, and implement the management of data storage, providing 
#' definitions for the interface's pure virtual methods.
#' 
#' The minimum requirement to define a new NMF model that integrates into 
#' the framework of the \emph{NMF} package are the followings:
#' 	
#' \itemize{
#' 
#' \item Define a class that inherits from class \code{NMF} and implements the 
#' new model, say class \code{myNMF}.
#' 
#' \item Implement the following S4 methods for the new class \code{myNMF}:
#'     \describe{
#'     \item{fitted}{\code{signature(object = "myNMF", value = "matrix")}: 
#'     Must return the estimated target matrix as fitted by the NMF model 
#'     \code{object}.
#'     }
#'     \item{basis}{\code{signature(object = "myNMF")}: 
#'     Must return the basis matrix(e.g. the first matrix factor in 
#'     the standard NMF model).
#'     }
#'     \item{basis<-}{\code{signature(object = "myNMF", value = "matrix")}: 
#'     Must return \code{object} with the basis matrix set to
#'     \code{value}.
#'     } 
#'     \item{coef}{\code{signature(object = "myNMF")}: 
#'     Must return the matrix of mixture coefficients (e.g. the second matrix 
#'     factor in the standard NMF model).
#'     }
#'     \item{coef<-}{\code{signature(object = "myNMF", value = "matrix")}:
#'     Must return \code{object} with the matrix of mixture coefficients set to 
#'     \code{value}.
#'     } 
#'     }
#'     
#' 	The \emph{NMF} package provides "pure virtual" definitions of these
#'  methods for class \code{NMF} (i.e. with signatures \code{(object='NMF', ...)} 
#'  and \code{(object='NMF', value='matrix')}) that throw an error if called, so 
#'  as to force their definition for model classes. 
#' 
#' \item Optionally, implement method \code{rnmf}(signature(x="myNMF", target="ANY")). 
#' This method should call \code{callNextMethod(x=x, target=target, ...)} and 
#' fill the returned NMF model with its specific data suitable random values. 
#' }
#' 
#' For concrete examples of NMF models implementations, see class 
#' \code{\linkS4class{NMFstd}} and its extensions (e.g. classes 
#' \code{\linkS4class{NMFOffset}} or \code{\linkS4class{NMFns}}).
#' 
#' @section Creating NMF objects:
#' Strictly speaking, because class \code{NMF} is virtual, no object of class 
#' \code{NMF} can be instantiated, only objects from its sub-classes. 
#' However, those objects are sometimes shortly referred in the documentation and 
#' vignettes as "\code{NMF} objects" instead of "objects that inherits from 
#' class \code{NMF}".
#' 
#' For built-in models or for models that inherit from the standard model class 
#' \code{\linkS4class{NMFstd}}, the factory method \code{nmfModel} enables to easily create 
#' valid \code{NMF} objects in a variety of common situations. 
#' See documentation for the the factory method \code{\link{nmfModel}} for 
#' more details.
#' 
#' @references
#' Definition of Nonnegative Matrix Factorization in its modern formulation: \cite{Lee1999}
#' 
#' Historical first definition and algorithms: \cite{Paatero1994}
#' 
#' @family NMF-model Implementations of NMF models
#' @seealso
#' Main interface to perform NMF in \code{\link{nmf-methods}}.
#'  
#' Built-in NMF models and factory method in \code{\link{nmfModel}}.
#' 
#' Method \code{\link{seed}} to set NMF objects with values suitable to start 
#' algorithms with.
#' 
#' @examples
#' 
#' # show all the NMF models available (i.e. the classes that inherit from class NMF)
#' nmfModels()
#' # show all the built-in NMF models available
#' nmfModels(builtin.only=TRUE)
#' 
#' # class NMF is a virtual class so cannot be instantiated: 
#' try( new('NMF') )
#' 
#' # To instantiate an NMF model, use the factory method nmfModel. see ?nmfModel
#' nmfModel()
#' nmfModel(3)
#' nmfModel(3, model='NMFns')
#' 
setClass('NMF'
		, representation(
			misc = 'list' # misceleneaous data used during fitting
		)
		, contains = 'VIRTUAL')

#' Fitted Matrix in NMF Models
#' 
#' Computes the estimated target matrix based on a given \emph{NMF} model.
#' The estimation depends on the underlying NMF model. 
#' For example in the standard model \eqn{V \equiv W H}{V ~ W H}, the target matrix is 
#' estimated by the matrix product \eqn{W H}.
#' In other models, the estimate may depend on extra parameters/matrix 
#' (cf. Non-smooth NMF in \code{\link{NMFns-class}}).
#' 
#' This function is a S4 generic function imported from \link[stats]{fitted} in 
#' the package \emph{stats}.
#' It is implemented as a pure virtual method for objects of class 
#' \code{NMF}, meaning that concrete NMF models must provide a 
#' definition for their corresponding class (i.e. sub-classes of 
#' class \code{NMF}). 
#' See \code{\linkS4class{NMF}} for more details.
#' 
#' @param object an object that inherit from class \code{NMF}
#' @param ... extra arguments to allow extension
#' 
#' @return the target matrix estimate as fitted by the model \code{object} 
#' @export
setGeneric('fitted', package='stats')
#' @template VirtualNMF  
setMethod('fitted', signature(object='NMF'),
		function(object, ...){
			stop("NMF::fitted is a pure virtual method of interface 'NMF'. It should be overloaded in class '", class(object),"'.")
		}
)

#' Accessing NMF Factors
#' 
#' \code{basis} and \code{basis<-} are S4 generic functions which respectively
#' extract and set the matrix of basis components of an NMF model 
#' (i.e. the first matrix factor).
#' 
#' For example, in the case of the standard NMF model \eqn{V \equiv W H}{V ~ W H}, 
#' the method \code{basis} will return the matrix \eqn{W}.
#' 
#' \code{basis} and \code{basis<-} are defined for the top 
#' virtual class \code{\linkS4class{NMF}} only, and rely internally on the low-level 
#' S4 generics \code{.basis} and \code{.basis<-} respectively that effectively 
#' extract/set the coefficient data.
#' These data are post/pre-processed, e.g., to extract/set only their 
#' non-fixed terms or check dimension compatibility.
#' 
#' @param object an object from which to extract the factor matrices, typically an 
#' object of class \code{\linkS4class{NMF}}.
#' @param ... extra arguments to allow extension and passed to the low-level 
#' access functions \code{.coef} and \code{.basis}.
#' 
#' Note that these throw an error if used in replacement functions \code{}.
#'   
#' @rdname basis-coef-methods
#' @family NMF-interface
#' @export
#'  
setGeneric('basis', function(object, ...) standardGeneric('basis') )
#' Default method returns the value of S3 slot or attribute \code{'basis'}.
#' It returns \code{NULL} if none of these are set. 
#' 
#' Arguments \code{...} are not used by this method.
setMethod('basis', signature(object='ANY'),
	function(object, ...){
		if( is.list(object) && 'basis' %in% names(object) ) object[['basis']]
		else attr(object, 'basis')
	}
)
#' @param all a logical that indicates whether the complete matrix factor 
#' should be returned (\code{TRUE}) or only the non-fixed part.
#' This is relevant only for formula-based NMF models that include fixed basis or 
#' coefficient terms.
#' 
#' @inline
setMethod('basis', signature(object='NMF'),
	function(object, all=TRUE, ...){
		if( all || !length(i <- ibterms(object)) ){
			# return all coefficients
			.basis(object, ...)
		} else {
			# remove fixed basis
			.basis(object, ...)[, -i]
		}
	}
)

#' \code{.basis} and \code{.basis<-} are the low-level S4 generics that simply 
#' return/set basis component data in an object.
#' They are defined so that some common processing may be implemented in 
#' \code{basis} and \code{basis<-}.
#' 
#' The methods \code{.basis}, \code{.coef} and their replacement versions
#' are implemented as pure virtual methods for the interface class 
#' \code{NMF}, meaning that concrete NMF models must provide a 
#' definition for their corresponding class (i.e. sub-classes of 
#' class \code{NMF}).
#' See \code{\linkS4class{NMF}} for more details.
#' 
#' @rdname basis-coef-methods
#' @export
setGeneric('.basis', function(object, ...) standardGeneric('.basis') )
#' @template VirtualNMF
setMethod('.basis', signature(object='NMF'),
	function(object, ...){
		stop("NMF::.basis is a pure virtual method of interface 'NMF'. It should be overloaded in class '", class(object),"'.")
	}
)

#' @export
#' @rdname basis-coef-methods
#' @inline
setGeneric('basis<-', function(object, ..., value) standardGeneric('basis<-') )
#' Default methods that calls \code{.basis<-} and check the validity of the 
#' updated object. 
#' @param use.dimnames logical that indicates if the object's dim names should be 
#' set using those from the new value, or left unchanged -- after truncating 
#' them to fit new dimensions if necessary.
#' This is useful to only set the entries of a factor.
#' 
setReplaceMethod('basis', signature(object='NMF', value='ANY'), 
	function(object, use.dimnames = TRUE, ..., value){
		
        # error if passed extra arguments
        if( length(xargs<- list(...)) ){
            stop("basis<-,NMF - Unused arguments: ", str_out(xargs, Inf, use.names = TRUE))
        }
        
        # backup old dimnames to reapply them on exit
        if( !use.dimnames ) odn <- dimnames(object)
        nb_old <- nbasis(object)
        
		# only set non-fixed terms
		if( !nbterms(object) ) .basis(object) <- value
		else{
			i <- ibasis(object)
			.basis(object)[,i] <- value[, i]
		}
        # adapt coef if empty
        if( !hasCoef(object) ){
            x <- basis(object)
            .coef(object) <- rbind(coef(object)[1:min(nb_old, ncol(x)), , drop = FALSE], matrix(NA, max(ncol(x)-nb_old, 0), 0))
#            .coef(object) <- coef(object)[1:ncol(x), , drop = FALSE] 
        }
		# check object validity
		validObject(object)
        
        # update other factor if necessary
        if( use.dimnames ) basisnames(object) <- colnames(basis(object))
        else if( !length(odn) ) dimnames(object) <- NULL
        else dimnames(object) <- mapply(head, odn, dim(object), SIMPLIFY = FALSE)
        
		object
	}	
)
#' @param value replacement value 
#' @rdname basis-coef-methods
#' @export
setGeneric('.basis<-', function(object, value) standardGeneric('.basis<-') )
#' @template VirtualNMF
setReplaceMethod('.basis', signature(object='NMF', value='matrix'), 
	function(object, value){ 
		stop("NMF::.basis<- is a pure virtual method of interface 'NMF'. It should be overloaded in class '", class(object),"'.")
	} 
)

#' @export
setGeneric('loadings', package='stats')
#' Method loadings for NMF Models
#' 
#' The method \code{loadings} is identical to \code{basis}, but do 
#' not accept any extra argument. 
#' 
#' The method \code{loadings} is provided to standardise the NMF interface 
#' against the one defined in the \code{\link{stats}} package, 
#' and emphasises the similarities between NMF and PCA or factorial analysis 
#' (see \code{\link{loadings}}).
#' 
#' @rdname basis-coef-methods
setMethod('loadings', 'NMF', function(x) basis(x) )

#' Get/Set the Coefficient Matrix in NMF Models
#' 
#' \code{coef} and \code{coef<-} respectively extract and set the 
#' coefficient matrix of an NMF model (i.e. the second matrix factor).  
#' For example, in the case of the standard NMF model \eqn{V \equiv WH}{V ~ W H}, 
#' the method \code{coef} will return the matrix \eqn{H}.
#' 
#' \code{coef} and \code{coef<-} are S4 methods defined for the corresponding 
#' generic functions from package \code{stats} (See \link[stats]{coef}).
#' Similarly to \code{basis} and \code{basis<-}, they are defined for the top 
#' virtual class \code{\linkS4class{NMF}} only, and rely internally on the S4 
#' generics \code{.coef} and \code{.coef<-} respectively that effectively 
#' extract/set the coefficient data.
#' These data are post/pre-processed, e.g., to extract/set only their 
#' non-fixed terms or check dimension compatibility.
#' 
#' @rdname basis-coef-methods
#' @export
setGeneric('coef', package='stats')
#' @inline
setMethod('coef', 'NMF',
	function(object, all=TRUE, ...){
		
		if( all || !length(i <- icterms(object)) ){
			# return all coefficients
			.coef(object, ...)
		} else {
			# remove fixed coefficients
			.coef(object, ...)[-i, ]
		}
		
	}
)

#' \code{.coef} and \code{.coef<-} are low-level S4 generics that simply 
#' return/set coefficient data in an object, leaving some common processing 
#' to be performed in \code{coef} and \code{coef<-}. 
#'    
#' @rdname basis-coef-methods
#' @export
setGeneric('.coef', function(object, ...) standardGeneric('.coef'))
#' @template VirtualNMF
setMethod('.coef', signature(object='NMF'),
	function(object, ...){
		stop("NMF::.coef is a pure virtual method of interface 'NMF'. It should be overloaded in class '", class(object),"'.")
	}
)

#' @export
#' @rdname basis-coef-methods
#' @inline
setGeneric('coef<-', function(object, ..., value) standardGeneric('coef<-') )
#' Default methods that calls \code{.coef<-} and check the validity of the 
#' updated object. 
setReplaceMethod('coef', signature(object='NMF', value='ANY'), 
	function(object, use.dimnames = TRUE, ..., value){
		
        # error if passed extra arguments
        if( length(xargs<- list(...)) ){
            stop("coef<-,NMF - Unused arguments: ", str_out(xargs, Inf, use.names = TRUE))
        }
        # backup old dimnames to reapply them on exit
        if( !use.dimnames ) odn <- dimnames(object)
        nb_old <- nbasis(object)
        
		# only set non-fixed terms
		if( !ncterms(object) ) .coef(object) <- value
		else{
			i <- icoef(object)
			.coef(object)[i, ] <- value[i, ]
		}
        # adapt basis if empty before validation
        if( !hasBasis(object) ){
            x <- coef(object)
            .basis(object) <- cbind(basis(object)[, 1:min(nb_old, nrow(x)), drop = FALSE], matrix(NA, 0, max(nrow(x)-nb_old, 0))) 
        }
		# check object validity
		validObject(object)
        
        # update other factor if necessary
        if( use.dimnames ) basisnames(object) <- rownames(coef(object))
        else if( !length(odn) ) dimnames(object) <- NULL
        else dimnames(object) <- mapply(head, odn, dim(object), SIMPLIFY = FALSE)
        
            
		object
	}	
)

#' @export
#' @rdname basis-coef-methods
setGeneric('.coef<-', function(object, value) standardGeneric('.coef<-') )
#' @template VirtualNMF
setReplaceMethod('.coef', signature(object='NMF', value='matrix'), 
	function(object, value){ 
		stop("NMF::.coef<- is a pure virtual method of interface 'NMF'. It should be overloaded in class '", class(object),"'.")
	} 
)

#' @description Methods \code{coefficients} and \code{coefficients<-} are 
#' simple aliases for methods \code{coef} and \code{coef<-} respectively.
#' 
#' @export
#' @rdname basis-coef-methods
setGeneric('coefficients', package='stats')
#' Alias to \code{coef,NMF}, therefore also pure virtual.
#' @inline
setMethod('coefficients', signature(object='NMF'), selectMethod('coef', 'NMF'))

#' @description \code{scoef} is similar to \code{coef}, but returns the mixture 
#' coefficient matrix of an NMF model, with the columns scaled so that they 
#' sum up to a given value (1 by default).
#' 
#' @param scale scaling factor, which indicates to the value the columns of the 
#' coefficient matrix should sum up to.
#' 
#' @rdname basis-coef-methods
#' @export
#' 
#' @examples 
#' 
#' # Scaled coefficient matrix
#' x <- rnmf(3, 10, 5)
#' scoef(x)
#' scoef(x, 100)
#' 
setGeneric('scoef', function(object, ...) standardGeneric('scoef') )
#' @inline
setMethod('scoef', 'NMF',
	function(object, scale=1){
		sweep(coef(object), 2L, colSums(coef(object)) / scale, '/')
	}
)
#' @inline
setMethod('scoef', 'matrix',
	function(object, scale=1){
		sweep(object, 2L, colSums(object) / scale, '/')
	}
)

unit.test(scoef, {
	x <- rnmf(3, 10, 5)
	checkIdentical(colSums(scoef(x)), rep(1, nbasis(x))
		, "Default call: columns are scaled to sum-up to one")
	checkIdentical(colSums(scoef(x, 100)), rep(1, nbasis(x))
		, "Scale=10: columns are scaled to sum-up to 10")
})

#' Rescaling NMF Models
#' 
#' Rescales an NMF model keeping the fitted target matrix identical. 
#' 
#' Standard NMF models are identifiable modulo a scaling factor, meaning that the
#' basis components and basis profiles can be rescaled without changing the fitted 
#' values:
#' 
#' \deqn{X = W_1 H_1 = (W_1 D) (D^{-1} H_1) = W_2 H_2}{X = W H = (W D) (D^-1 H)}
#' with \eqn{D= \alpha diag(1/\delta_1, \ldots, 1\delta_r)}{D= alpha * diag(1/delta_1, ..., 1/delta_r)}
#' 
#' The default call \code{scale(object)} rescales the basis NMF object so that each   
#' column of the basis matrix sums up to one. 
#' 
#' @param x an NMF object
#' @param center either a numeric normalising vector \eqn{\delta}{delta}, or either 
#' \code{'basis'} or \code{'coef'}, which respectively correspond to using the 
#' column sums of the basis matrix or the inverse of the row sums of the 
#' coefficient matrix as a normalising vector.
#' If numeric, \code{center} should be a single value or a vector of length the 
#' rank of the NMF model, i.e. the number of columns in the basis matrix.
#' @param scale scaling coefficient applied to \eqn{D}, i.e. the value of \eqn{\alpha}{alpha}, 
#' or, if \code{center='coef'}, the value of \eqn{1/\alpha}{1/alpha} (see section \emph{Details}).
#' 
#' @return an NMF object
#' 
#' @S3method scale NMF
#' @examples
#' 
#' # random 3-rank 10x5 NMF model
#' x <- rnmf(3, 10, 5)
#' 
#' # rescale based on basis
#' colSums(basis(x))
#' colSums(basis(scale(x)))
#' 
#' rx <- scale(x, 'basis', 10)
#' colSums(basis(rx))
#' rowSums(coef(rx))
#' 
#' # rescale based on coef
#' rowSums(coef(x))
#' rowSums(coef(scale(x, 'coef')))
#' rx <- scale(x, 'coef', 10)
#' rowSums(coef(rx))
#' colSums(basis(rx))
#' 
#' # fitted target matrix is identical but the factors have been rescaled
#' rx <- scale(x, 'basis')
#' all.equal(fitted(x), fitted(rx))
#' all.equal(basis(x), basis(rx))
#' 
scale.NMF <- function(x, center=c('basis', 'coef'), scale=1){
	
	# determine base value
	if( missing(center) ) center <- match.arg(center)
	base <- center
	delta <-
	if( is.character(base) ){
		base <- match.arg(center)
		if( base == 'basis' ) colSums(basis(x))
		else{
			scale <- 1/scale 
			1 / rowSums(coef(x))
		}	
	}else if( is.numeric(base) ) base
	else stop("Invalid base value: should be a numeric or one of "
			, str_out(c('none', 'basis', 'coef')))

	# scale
	D <- scale/delta 
	# W <- W * D
	basis(x) <- sweep(basis(x), 2L, D, '*')
	# H <- D^-1 * H
	coef(x) <- sweep(coef(x), 1L, D, '/')
	x
	
}

unit.test("scale", {
			
	r <- 3
	x <- rnmf(r, 10, 5)
	
	.lcheck <- function(msg, rx, ref, target){
		.msg <- function(...) paste(msg, ':', ...)
		checkTrue(!identical(basis(x), basis(rx)), .msg("changes basis matrix"))
		checkTrue(!identical(coef(x), coef(rx)), .msg("changes coef matrix"))
		checkEqualsNumeric(fitted(x), fitted(rx), .msg("fitted target is identical"))
		
		brx <- colSums(basis(rx))
		crx <- rowSums(coef(rx)) 
		if( target == 1 ){
			checkEquals(brx, ref, .msg("correctly scales basis components"))
			checkTrue(!all(crx==ref), .msg("does not force scale on coefficient matrix"))
		}else{
			checkTrue(!all(brx==ref), .msg("does not force scale on basis matrix"))
			checkEquals(crx, ref
					, .msg("correctly scales rows of coef matrix"))
		}
	}
	
	.check <- function(msg, ref, ...){
		.lcheck(str_c(msg, " + argument center='basis'")
				, scale(x, center='basis', ...), ref, 1)
		.lcheck(str_c(msg, " + argument center='coef'")
				, scale(x, center='coef', ...), ref, 2)
	}
	
	.lcheck("Default call", scale(x), rep(1, r), 1)
	.check("Missing argument scale", rep(1, r))
	.check("Argument scale=10", rep(10, r), scale=10)
	s <- runif(r)
	.check("Argument scale=numeric", s, scale=s)
 
})

#' Generating Random NMF Models
#' 
#' Generates NMF models with random values drawn from a uniform distribution.
#' It returns an NMF model with basis and mixture coefficient matrices filled
#' with random values. 
#' The main purpose of the function \code{rnmf} is to provide a common 
#' interface to generate random seeds used by the \code{\link{nmf}} function.
#' 
#' If necessary, extensions of the standard NMF model or custom models must
#' define a method "rnmf,<NMF.MODEL.CLASS>,numeric" for initialising their
#' specific slots other than the basis and mixture coefficient matrices. 
#' In order to benefit from the complete built-in interface, the overloading
#' methods should call the generic version using function
#' \code{\link{callNextMethod}}, prior to set the values of the specific slots.
#' See for example the method \code{\link[=rnmf,NMFOffset,numeric-method]{rnmf}} 
#' defined for \code{\linkS4class{NMFOffset}} models:
#' \code{showMethods(rnmf, class='NMFOffset', include=TRUE))}.
#' 
#' For convenience, shortcut methods for working on \code{data.frame} objects 
#' directly are implemented.
#' However, note that conversion of a \code{data.frame} into a \code{matrix} 
#' object may take some non-negligible time, for large datasets. 
#' If using this method or other NMF-related methods several times, consider 
#' converting your data \code{data.frame} object into a matrix once for good, 
#' when first loaded.
#' 
#' @param x an object that determines the rank, dimension and/or class of the 
#' generated NMF model, e.g. a numeric value or an object that inherits from class 
#' \code{\linkS4class{NMF}}.
#' See the description of the specific methods for more details on the supported 
#' types.
#' @param target optional specification of target dimensions. 
#' See section \emph{Methods} for how this parameter is used by the different 
#' methods.
#' @param ... extra arguments to allow extensions and passed to the next method 
#' eventually down to \code{\link{nmfModel}}, where they are used to initialise 
#' slots that are specific to the instantiating NMF model.
#' 
#' @return An NMF model, i.e. an object that inherits from class 
#' \code{\linkS4class{NMF}}.
#' 
#' @inline
#' @export
#' @seealso \code{\link{rmatrix}}
#' @family NMF-interface
setGeneric('rnmf', function(x, target, ...) standardGeneric('rnmf') )

# Define the loading namespace
.PKG.NAMESPACE <- packageEnv()

#' Testing NMF Objects
#' 
#' @description
#' The functions documented here tests different characteristics of NMF objects.
#'
#' \code{is.nmf} tests if an object is an NMF model or a class that extends
#' the class NMF.
#' 
#' @details
#' 
#' \code{is.nmf} tests if \code{object} is the name of a class (if a \code{character} 
#' string), or inherits from a class, that extends \code{\linkS4class{NMF}}.
#' 
#' @note The function \code{is.nmf} does some extra work with the namespace as 
#' this function needs to return correct results even when called in \code{.onLoad}.
#' See discussion on r-devel: \url{https://stat.ethz.ch/pipermail/r-devel/2011-June/061357.html}
#' 
#' @param x an R object. See section \emph{Details}, for how each function
#' uses this argument.
#'  
#' @rdname types
#' @export
#' 
#' @examples
#' 
#' # test if an object is an NMF model, i.e. that it implements the NMF interface
#' is.nmf(1:4)
#' is.nmf('NMFstd')
#' is.nmf('NMFblah')
#' is.nmf( nmfModel(3) )
#' is.nmf( nmf(rmatrix(20,10), 3) )
#' 
is.nmf <- function(x){
	
	# load definition for base class NMF
	clref <- getClass('NMF', .Force=TRUE, where=.PKG.NAMESPACE)
	is(x, clref)
	
}

unit.test(is.nmf,{
	checkTrue(!is.nmf(1:4), "on vector: FALSE")
	checkTrue(!is.nmf(list(1:4)), "on list: FALSE")
	checkTrue(is.nmf('NMF'), "on 'NMF': TRUE")
	checkTrue(is.nmf('NMFstd'), "on 'NMFstd': TRUE")
	checkTrue( is.nmf( nmfModel(3) ), "on empty model: TRUE")
	checkTrue( is.nmf( rnmf(3, 20, 10) ), "on random model: TRUE")
	checkTrue( is.nmf( nmf(rmatrix(20,10), 3) ), "on NMFfit object: TRUE") 
})

isNMFclass <- function(x){
	
	if( is.character(x) ){ # test object is a class that extends NMF
		# load definition for base class NMF
		clref <- getClass('NMF', .Force=TRUE, where=.PKG.NAMESPACE)
		cl <- getClass(x, .Force=TRUE, where=.PKG.NAMESPACE)
		if( is.null(cl) )
			cl <- getClass(x, .Force=TRUE)
		extends(cl, clref)
	}else 
		FALSE
	
}

################################

# Taken from Biobase
selectSome <- function (obj, maxToShow = 5) 
{
    len <- length(obj)
    if (maxToShow < 3) 
        maxToShow <- 3
    if (len > maxToShow) {
        maxToShow <- maxToShow - 1
        bot <- ceiling(maxToShow/2)
        top <- len - (maxToShow - bot - 1)
        nms <- obj[c(1:bot, top:len)]
        c(as.character(nms[1:bot]), "...", as.character(nms[-c(1:bot)]))
    }
    else if (is.factor(obj)) 
        as.character(obj)
    else obj
}


.showFixedTerms <- function(x, ...){
	
	s <- 
	sapply(x, function(t){
		s <- 
		if( is.factor(t) ) selectSome(levels(t), ...)
		else selectSome(t, ...)
		s <- str_out(s, Inf, quote=FALSE)
		if( is.factor(t) )
    		s <- str_c('<', s, ">")
		s
	})
	paste(names(s), '=', s)
}

#' Show method for objects of class \code{NMF}
#' @export
setMethod('show', 'NMF', 
		function(object)
		{
			cat("<Object of class:", class(object), ">\n", sep='')
			cat("features:", nrow(object), "\n")
			cat("basis/rank:", nbasis(object), "\n")
			cat("samples:", ncol(object), "\n")
			# show fixed terms
			if( (n <- ncterms(object)) ){
				cat("fixed coef [", n, "]:\n" 
					, str_c('  ', .showFixedTerms(cterms(object), 4), collapse="\n")
					, "\n", sep='')
			}
			if( (n <- nbterms(object)) ){
				cat("fixed basis [", n, "]:\n" 
						, str_c('  ', .showFixedTerms(bterms(object), 4), collapse="\n")
						, "\n", sep='')
			}
			# show the miscellaneous model parameters
			if( length(object@misc) > 0L ){
				cat("miscellaneous:", str_desc(object@misc, exdent=12L), ". (use 'misc(object)')\n")
			}
		}
)


#' Dimension of NMF Objects
#' 
#' @description
#' The methods \code{dim}, \code{nrow}, \code{ncol} and \code{nbasis} return
#' the different dimensions associated with an NMF model.
#' 
#' \code{dim} returns all dimensions in a length-3 integer vector:
#' the number of row and columns of the estimated target matrix, 
#' as well as the factorization rank (i.e. the number of basis components).
#'  
#' \code{nrow}, \code{ncol} and \code{nbasis} provide separate access to each 
#' of these dimensions respectively.
#' 
#' @details
#' The NMF package does not implement specific functions \code{nrow} and \code{ncol},
#' but rather the S4 method \code{dim} for objects of class \code{\linkS4class{NMF}}.
#' This allows the base methods \code{\link{nrow}} and \code{\link{ncol}} to 
#' directly work with such objects, to get the number of rows and columns of 
#' the target matrix estimated by an NMF model.
#' 
#' The function \code{nbasis} is a new S4 generic defined in the package NMF, that
#' returns the number of basis components of an object.
#' Its default method should work for any object, that has a suitable 
#' \code{basis} method defined for its class.
#' 
#' @param x an object with suitable \code{basis} and \code{coef} methods, such 
#' as an object that inherit from \code{\linkS4class{NMF}}.
#' @param ... extra arguments to allow extension.
#' 
#' @return a single integer value or, for \code{dim}, a length-3 integer vector, 
#' e.g. \code{c(2000, 30, 3)} for an \code{NMF} model that fits a 2000 x 30 
#' matrix using 3 basis components.
#' 
#' @export
#' @rdname dims
#' @aliases dim-NMF
setGeneric('nbasis', function(x, ...) standardGeneric('nbasis') )
#' Default method which returns the number of columns of the basis matrix extracted 
#' from \code{x} using a suitable method \code{basis}, or, if the latter is \code{NULL}, 
#' the value of attributes \code{'nbasis'}.
#' 
#' For NMF models, this also corresponds to the number of rows in the coefficient
#' matrix.
#'    
setMethod('nbasis', signature(x='ANY'), 
	function(x, ...)
	{
		if( !is.null(n <- ncol(basis(x, ...))) ) n
		else if( is.list(x) && 'nbasis' %in% names(x) ) x[['nbasis']]
		else attr(x, 'nbasis')
	}
)
#' method for NMF objects for the base generic \code{\link{dim}}.
#' It returns all dimensions in a length-3 integer vector:
#' the number of row and columns of the estimated target matrix, 
#' as well as the factorization rank (i.e. the number of basis components).
#' 
#' @rdname dims
#' @export
setMethod('dim', signature(x='NMF'), 
	function(x){
		c(nrow(basis(x)), ncol(coef(x)), nbasis(x))	
	}
)


#' Dimension names for NMF objects
#' 
#' @description
#' The methods \code{dimnames}, \code{rownames}, \code{colnames} and
#' \code{basisnames} and their respective replacement form allow to get and set
#' the dimension names of the matrix factors in a NMF model.
#' 
#' \code{dimnames} returns all the dimension names in a single list.
#' Its replacement form \code{dimnames<-} allows to set all dimension names at once.
#' 
#' \code{rownames}, \code{colnames} and \code{basisnames} provide separate access 
#' to each of these dimension names respectively.
#' Their respective replacement form allow to set each dimension names separately.
#' 
#' @details
#' 
#' The function \code{basisnames} is a new S4 generic defined in the package NMF, 
#' that returns the names of the basis components of an object.
#' Its default method should work for any object, that has a suitable \code{basis} 
#' method defined for its class.
#' 
#' The method \code{dimnames} is implemented for the base generic \code{\link{dimnames}}, 
#' which make the base function \code{\link{rownames}} and \code{\link{colnames}} 
#' work directly.
#' 
#' Overall, these methods behave as their equivalent on \code{matrix} objects.
#' The function \code{basisnames<-} ensures that the dimension names are handled 
#' in a consistent way on both factors, enforcing the names on both matrix factors
#' simultaneously.
#' 
#' @param x an object with suitable \code{basis} and \code{coef} methods, such 
#' as an object that inherit from \code{\linkS4class{NMF}}.
#' @param ...  extra argument to allow extension.
#' 
#' @export
#' @rdname dimnames
#' @aliases dimnames-NMF
#' 
#' @examples
#' # create a random NMF object
#' a <- rnmf(2, 5, 3)
#' 
#' # set dimensions
#' dims <- list( features=paste('f', 1:nrow(a), sep='')
#' 				, samples=paste('s', 1:ncol(a), sep='')
#' 				, basis=paste('b', 1:nbasis(a), sep='') )
#' dimnames(a) <- dims
#' dimnames(a)
#' basis(a)
#' coef(a)
#' 
#' # access the dimensions separately
#' rownames(a)
#' colnames(a)
#' basisnames(a)
#' 
#' # set only the first dimension (rows of basis): the other two dimnames are set to NULL
#' dimnames(a) <- dims[1]
#' dimnames(a)
#' basis(a)
#' coef(a)
#' 
#' # set only the two first dimensions (rows and columns of basis and coef respectively):
#' # the basisnames are set to NULL 
#' dimnames(a) <- dims[1:2]
#' dimnames(a)
#' basis(a)
#' 
#' # reset the dimensions
#' dimnames(a) <- NULL
#' dimnames(a)
#' basis(a)
#' coef(a)
#' 
#' # set each dimensions separately
#' rownames(a) <- paste('X', 1:nrow(a), sep='') # only affect rows of basis
#' basis(a)
#' 
#' colnames(a) <- paste('Y', 1:ncol(a), sep='') # only affect columns of coef
#' coef(a)
#' 
#' basisnames(a) <- paste('Z', 1:nbasis(a), sep='') # affect both basis and coef matrices
#' basis(a)
#' coef(a)
#' 
setGeneric('basisnames', function(x, ...) standardGeneric('basisnames') )
#' Default method which returns the column names of the basis matrix extracted from 
#' \code{x}, using the \code{basis} method.
#' 
#' For NMF objects these also correspond to the row names of the coefficient matrix.
setMethod('basisnames', signature(x='ANY'), 
	function(x)
	{
		colnames(basis(x))
	}
)
#' The generic \code{basisnames<-} simultaneously sets the names of the basis 
#' components and coefficients of an object, for which suitable \code{basis} 
#' and \code{coef} methods are defined.
#' 
#' @details
#' The function \code{basisnames<-} is a new S4 generic defined in the package NMF, 
#' that sets the names of the basis components of an object.
#' Its default method should work for any object, that has suitable \code{basis<-} 
#' and \code{coef<-} methods method defined for its class.
#' 
#' @export
#' @inline
#' @rdname dimnames 
setGeneric('basisnames<-', function(x, ..., value) standardGeneric('basisnames<-') )
#' Default method which sets, respectively, the row and the column names of the basis 
#' matrix and coefficient matrix of \code{x} to \code{value}.
setReplaceMethod('basisnames', 'ANY', 
	function(x, ..., value)
	{
		rownames(.coef(x)) <- value
		colnames(.basis(x)) <- value
		x
	}
)

#' Returns the dimension names of the NMF model \code{x}.
#' 
#' It returns either NULL if no dimnames are set on the object, 
#' or a 3-length list containing the row names of the basis matrix,
#' the column names of the mixture coefficient matrix, and the column names of
#' the basis matrix (i.e. the names of the basis components).
#' 
#' @param value a character vector, or \code{NULL} or, in the case of
#' \code{dimnames<-}, a list 2 or 3-length list of character vectors.
#'
#' @rdname dimnames
#' @export
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
#' sets the dimension names of the NMF model \code{x}.  
#' 
#' \code{value} can be \code{NULL} which resets all dimension names, or a 
#' 1, 2 or 3-length list providing names at least for the rows of the basis 
#' matrix.
#'   
#' The optional second element of \code{value} (NULL if absent) is used to set 
#' the column names of the coefficient matrix.  
#' The optional third element of \code{value} (NULL if absent) is used to set 
#' both the column names of the basis matrix and the row names of the 
#' coefficient matrix.
#' 
#' @rdname dimnames
#' @export
setReplaceMethod('dimnames', 'NMF', 
	function(x, value){
		if( !is.list(value) && !is.null(value) )
			stop("NMF::dimnames - Invalid value: must be a list or NULL.")
		
		if( length(value) == 0 )
			value <- NULL
		else if( length(value) == 1 )
			value <- c(value, list(NULL, NULL))			
		else if( length(value) == 2 ) # if only the two first dimensions reset the third one
			value <- c(value, list(NULL))
		else if( length(value)!=3 ) # check length of value
			stop("NMF::dimnames - invalid argument 'value' [a 2 or 3-length list is expected]")
		
		# only set relevant dimensions
		if( length(w <- which(dim(x) == 0)) ){
			value[w] <- sapply(value[w], function(x) NULL, simplify=FALSE)
		}
		# set dimnames 
		dimnames(.basis(x)) <- value[c(1,3)]
		dimnames(.coef(x)) <- value[c(3,2)]
		# return updated model
		x	
	}
)

#' Sub-setting NMF Objects
#' 
#' This method provides a convenient way of sub-setting objects of class \code{NMF}, 
#' using a matrix-like syntax.
#' 
#' It allows to consistently subset one or both matrix factors in the NMF model, as well 
#' as retrieving part of the basis components or part of the mixture coefficients with 
#' a reduced amount of code.
#' 
#' @details
#' The returned value depends on the number of subset index passed and the
#' value of argument \code{drop}:
#' 
#' \itemize{ \item No index as in \code{x[]} or \code{x[,]}: the value is the
#' object \code{x} unchanged.
#' 
#' \item One single index as in \code{x[i]}: the value is the complete NMF 
#' model composed of the selected basis components, subset by \code{i}, 
#' except if argument \code{drop=TRUE}, or if it is missing and \code{i} is of length 1.
#' Then only the basis matrix is returned with dropped dimensions: 
#' \code{x[i, drop=TRUE]} <=> \code{drop(basis(x)[, i])}.
#' 
#' This means for example that \code{x[1L]} is the first basis vector, 
#' and \code{x[1:3, drop = TRUE]} is the matrix composed of the 3 first basis vectors -- in columns.
#' 
#' Note that in version <= 0.18.3, the call \code{x[i, drop = TRUE.or.FALSE]} was equivalent to
#' \code{basis(x)[, i, drop=TRUE.or.FALSE]}.
#' 
#' \item More than one index with \code{drop=FALSE} (default) as in
#' \code{x[i,j]}, \code{x[i,]}, \code{x[,j]}, \code{x[i,j,k]}, \code{x[i,,k]},
#' etc...: the value is a \code{NMF} object whose basis and/or mixture
#' coefficient matrices have been subset accordingly. The third index \code{k}
#' affects simultaneously the columns of the basis matrix AND the rows of the
#' mixture coefficient matrix. In this case argument \code{drop} is not used.
#' 
#' \item More than one index with \code{drop=TRUE} and \code{i} xor \code{j}
#' missing: the value returned is the matrix that is the more affected by the
#' subset index. That is that \code{x[i, , drop=TRUE]} and \code{x[i, , k,
#' drop=TRUE]} return the basis matrix subset by \code{[i,]} and \code{[i,k]}
#' respectively, while \code{x[, j, drop=TRUE]} and \code{x[, j, k, drop=TRUE]}
#' return the mixture coefficient matrix subset by \code{[,j]} and \code{[k,j]}
#' respectively.
#' 
#' }
#' 
#' @param i index used to subset on the \strong{rows} of the basis matrix (i.e.
#' the features).
#' It can be a \code{numeric}, \code{logical}, or \code{character} vector 
#' (whose elements must match the row names of \code{x}). 
#' In the case of a \code{logical} vector the entries are recycled if necessary.
#' @param j index used to subset on the \strong{columns} of the mixture
#' coefficient matrix (i.e. the samples).  
#' It can be a \code{numeric}, \code{logical}, or \code{character} vector 
#' (whose elements must match the column names of \code{x}).
#' In the case of a \code{logical} vector the entries are recycled if necessary.
#' @param ...  used to specify a third index to subset on the basis components,
#' i.e. on both the columns and rows of the basis matrix and mixture
#' coefficient respectively.  
#' It can be a \code{numeric}, \code{logical}, or \code{character} vector 
#' (whose elements must match the basis names of \code{x}).
#' In the case of a \code{logical} vector the entries are recycled if necessary.
#' 
#' Note that only the first extra subset index is used.
#' A warning is thrown if more than one extra argument is passed in \code{...}.
#' @param drop single \code{logical} value used to drop the \code{NMF-class}
#' wrapping and only return subsets of one of the factor matrices (see \emph{Details})
#' 
#' @rdname subset-NMF
#' @export
#' @examples
#' # create a dummy NMF object that highlight the different way of subsetting
#' a <- nmfModel(W=outer(seq(1,5),10^(0:2)), H=outer(10^(0:2),seq(-1,-10)))
#' basisnames(a) <- paste('b', 1:nbasis(a), sep='')
#' rownames(a) <- paste('f', 1:nrow(a), sep='')
#' colnames(a) <- paste('s', 1:ncol(a), sep='')
#' 
#' # or alternatively:
#' # dimnames(a) <- list( features=paste('f', 1:nrow(a), sep='')
#' #					, samples=paste('s', 1:ncol(a), sep='')
#' #					, basis=paste('b', 1:nbasis(a)) )
#' 
#' # look at the resulting NMF object 
#' a
#' basis(a)
#' coef(a)
#' 
#' # extract basis components
#' a[1]
#' a[1, drop=FALSE] # not dropping matrix dimension
#' a[2:3]
#' 
#' # subset on the features
#' a[1,]
#' a[2:4,]
#' # dropping the NMF-class wrapping => return subset basis matrix
#' a[2:4,, drop=TRUE]
#' 
#' # subset on the samples
#' a[,1]
#' a[,2:4]
#' # dropping the NMF-class wrapping => return subset coef matrix
#' a[,2:4, drop=TRUE]
#' 
#' # subset on the basis => subsets simultaneously basis and coef matrix
#' a[,,1]
#' a[,,2:3]
#' a[4:5,,2:3]
#' a[4:5,,2:3, drop=TRUE] # return subset basis matrix
#' a[,4:5,2:3, drop=TRUE] # return subset coef matrix
#' 
#' # 'drop' has no effect here
#' a[,,2:3, drop=TRUE]
#' 
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
			.basis(x) <- basis(x)[i, , drop = FALSE]
		
		# subset the columns of mixture coefficient matrix		
		if (!missing(j)) 
			.coef(x) <- coef(x)[, j, drop = FALSE]
		
		# subset the basis: columns of basis matrix and row of mixture coefficient matrix		
		if( single.arg || k.notmissing ){
			.basis(x) <- basis(x)[, k, drop = FALSE]
			# return basis only single arg and drop=TRUE 
			if( single.arg && ((mdrop && length(k) == 1L) || drop) ) return( drop(basis(x)) )
			.coef(x) <- coef(x)[k, , drop = FALSE]
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

#' The function \code{misc} provides access to miscellaneous data members stored 
#' in slot \code{misc} (as a \code{list}), which allow extensions of NMF models  
#' to be implemented, without defining a new S4 class.
#' 
#' @param object an object that inherit from class \code{NMF}
#' @param ... extra arguments (not used)
#' 
#' @rdname NMF-class
#' @export
misc <- function(object, ...){
	if( !isS4(object) && is.list(object) ) object[['misc']]
	else attr(object, 'misc')
}
#' shortcut for \code{x@@misc[[name, exact=TRUE]]} respectively.
#' @rdname NMF-class
#' @export 
setMethod('$', 'NMF', 
		function(x, name){ 
			x@misc[[name, exact=TRUE]]; 
		} 
)
#' shortcut for \code{x@@misc[[name]] <- value}
#' @rdname NMF-class
#' @export
setReplaceMethod('$', 'NMF',
	function(x, name, value) {
		x@misc[[name]] <- value
		x
	}
)

#' @importFrom utils .DollarNames
setGeneric('.DollarNames', package='utils')

#' @S3method .DollarNames NMF
.DollarNames.NMF <- function(x, pattern = "") grep(pattern, names(misc(x)), value=TRUE)

#' Auto-completion for \code{\linkS4class{NMF}} objects
#' @rdname NMF-class
#' @export
setMethod('.DollarNames', 'NMF', .DollarNames.NMF)

#' \code{is.empty.nmf} tests whether an \code{NMF} object describes an empty NMF model, 
#' i.e. it contains no data.
#' 
#' @details
#' \code{is.empty.nmf} returns \code{TRUE} if the basis and coefficient matrices of 
#' \code{x} have respectively zero rows and zero columns.
#' It returns \code{FALSE} otherwise.
#' 
#' In particular, this means that an empty model can still have a non-zero number 
#' of basis components, i.e. a factorization rank that is not null. 
#' This happens, for example, in the case of NMF models created calling the factory method 
#' \code{\link{nmfModel}} with a value only for the factorization rank.
#' 
#' @param ... extra parameters to allow extension or passed to subsequent calls
#' 
#' @rdname types
#' @export
#' 
#' @examples
#' 
#' # empty model
#' is.empty.nmf( nmfModel(3) )
#' # non empty models
#' is.empty.nmf( nmfModel(3, 10, 0) )
#' is.empty.nmf( rnmf(3, 10, 5) )
#'   
is.empty.nmf <- function(x, ...){
	nrow(x) == 0 && ncol(x) == 0
}


#' \code{hasBasis} tests whether an objects contains a basis matrix -- returned by 
#' a suitable method \code{basis} -- with at least one row.
#' 
#' @rdname types
#' @export
hasBasis <- function(x) nbasis(x) && nrow(basis(x)) != 0L

#' \code{hasBasis} tests whether an objects contains a coefficient matrix 
#' -- returned by a suitable method \code{coef} -- with at least one column.
#' 
#' @rdname types
#' @export
hasCoef <- function(x) nbasis(x) && ncol(coef(x)) != 0L

#' \code{is.partial.nmf} tests whether an NMF model object contains either an empty 
#' basis or coefficient matrix.
#' It is a shorcut for \code{!hasCoef(x) || !hasBasis(x)}.
#' 
#' @rdname types
#' @export
is.partial.nmf <- function(x) !hasCoef(x) || !hasBasis(x)

#' Returns the target matrix estimate of the NMF model \code{x}, perturbated by  
#' adding a random matrix generated using the default method of \code{rmatrix}: 
#' it is a equivalent to \code{fitted(x) + rmatrix(fitted(x), ...)}.
#' 
#' This method can be used to generate random target matrices that depart from 
#' a known NMF model to a controlled extend.
#' This is useful to test the robustness of NMF algorithms to the presence of 
#' certain types of noise in the data.
#' 
#' @examples
#' # generate noisy fitted target from an NMF model (the true model)
#' gr <- as.numeric(mapply(rep, 1:3, 3))
#' h <- outer(1:3, gr, '==') + 0 
#' x <- rnmf(10, H=h)
#' y <- rmatrix(x)
#' \dontrun{
#' # show heatmap of the noisy target matrix: block patterns should be clear
#' aheatmap(y) 
#' }
#' \dontshow{ stopifnot( identical(dim(y), dim(x)[1:2]) ) }
#' 
#' # test NMF algorithm on noisy data
#' # add some noise to the true model (drawn from uniform [0,1])
#' res <- nmf(rmatrix(x), 3)
#' summary(res)
#' 
#' # add more noise to the true model (drawn from uniform [0,10])
#' res <- nmf(rmatrix(x, max=10), 3)
#' summary(res)
#' 
setMethod('rmatrix', 'NMF', 
	function(x, ...){
		a <- fitted(x)
		a + rmatrix(a, ...)
	}
)

unit.test('rmatrix,NMF',{
	
	x <- nmfModel(3, 20, 5)
	checTrue(is.matrix(y <- rmatrix(x)), "default call: no error")
	checkIdentical(dim(y), dim(x)[1:2], "default call: correct dimension")
	checkTrue( !any(is.na(basis(y))), 'default call: no NAs in basis anymore')
	checkTrue( !any(is.na(coef(y))), 'default call: no NAs in coef anymore')
	checkTrue( max( max(abs(basis(y)-basis(x))), max(abs(coef(y)-coef(x))) ) <= 1
			, "default call: max difference is <= 1")
	
	set.seed(123)
	y <- rmatrix(x)
	set.seed(123)
	ref <- matrix(runif(nrow(x)*ncol(x)), nrow(x))
	checkIdentical(ref, y - fitted(x), "default call: add uniform random noise to fitted matrix")
	
	set.seed(123)
	ref <- matrix(rnorm(nrow(x)*ncol(x)), nrow(x))
	set.seed(123)
	y <- rmatrix(x, rnorm)	
	checkIdentical(ref, y - fitted(x), "dist is taken into account: add normal random noise to fitted matrix")
	set.seed(123)
	y <- rmatrix(x, dist=rnorm)
	checkIdentical(ref, y - fitted(x), "dist is taken into account: add normal random noise to fitted matrix")
	
	set.seed(123)
	checTrue(is.matrix(y <- rmatrix(x, max=10)), "call with arg max=10: no error")
	checkTrue( max( max(abs(basis(y)-basis(x))), max(abs(coef(y)-coef(x))) ) <= 10
			, "call with arg max=10: max difference is 10")
	checkTrue( max( max(abs(basis(y)-basis(x))), max(abs(coef(y)-coef(x))) ) >= 5
			, "call with arg max=10: max difference is >= 5")
	
})

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
	expand_list(input.list, ..., .exact=FALSE)
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
#' Assessing and Comparing NMF Models
#' 
#' @description
#' The NMF package defines \code{summary} methods for different classes of objects, 
#' which helps assessing and comparing the quality of NMF models by computing a set 
#' of quantitative measures, e.g. with respect to their ability to recover known 
#' classes and/or the original target matrix.
#' 
#' The most useful methods are for classes \code{\linkS4class{NMF}}, \code{\linkS4class{NMFfit}},
#' \code{\linkS4class{NMFfitX}} and \code{\linkS4class{NMFList}}, which compute summary measures  
#' for, respectively, a single NMF model, a single fit, a multiple-run fit and a list of heterogenous 
#' fits performed with the function \code{\link{nmf}}.
#' 
#' @details
#' Due to the somehow hierarchical structure of the classes mentionned in \emph{Description}, 
#' their respective \code{summary} methods call each other in chain, each super-class adding some 
#' extra measures, only relevant for objects of a specific class. 
#' 
#' @param object an NMF object. See available methods in section \emph{Methods}.
#' @param ... extra arguments passed to the next \code{summary} method.   
#' 
#' @export
#' @rdname assess
#' @aliases summary-NMF
#' 
#' @family assess Assessment measures for NMF models
#' 
setGeneric('summary', package='base')

#' Computes summary measures for a single NMF model.
#' 
#' The following measures are computed:
#' 
#' \describe{
#' \item{sparseness}{Sparseness of the factorization computed by the 
#' function \code{\link{sparseness}}.}
#' \item{entropy}{Purity of the clustering, with respect to known classes, 
#' computed by the function \code{\link{purity}}.}
#' \item{entropy}{Entropy of the clustering, with respect to known classes, 
#' computed by the function \code{\link{entropy}}.}
#' \item{RSS}{Residual Sum of Squares computed by the function \code{\link{rss}}.}
#' \item{evar}{Explained variance computed by the function \code{\link{evar}}.}
#' }
#' 
#' @param class known classes/cluster of samples specified in one of the formats
#' that is supported by the functions \code{\link{entropy}} and \code{\link{purity}}.
#' @param target target matrix specified in one of the formats supported by the 
#' functions \code{\link{rss}} and \code{\link{evar}}  
#' 
#' @rdname assess
#' 
#' @examples 
#' 
#' # random NMF model
#' x <- rnmf(3, 20, 12)
#' summary(x)
#' summary(x, gl(3, 4))
#' summary(x, target=rmatrix(x))
#' summary(x, gl(3,4), target=rmatrix(x))
#' 
setMethod('summary', signature(object='NMF'), 
		function(object, class, target){
			
			res <- numeric()
			
			## IMPORTANT: if adding a summary measure also add it in the sorting 
			## schema of method NMFList::summary to allow ordering on it
			
			# rank
			res <- c(res, rank=nbasis(object))
			# compute sparseness
			res <- c(res, sparseness=sparseness(object))
			
			# if class is provided: also computes entropy and purity
			if( !missing(class) ){
				# compute purity
				res <- c(res, purity=purity(object, class))
				# compute entropy
				res <- c(res, entropy=entropy(object, class))
			}
			
			# if the target is provided compute the RSS
			if( !missing(target) ){
				RSS <- rss(object, target)
				res <- c(res, rss=RSS)
				# explained variance
				res <- c(res, evar=evar(object, target))
			}
            
            # compute mean silhouette width
            siS <- silhouette(object, what = 'samples')
            siF <- silhouette(object, what = 'features')
            res <- c(res, silhouette.coef = if( !is_NA(siS) ) summary(siS)$avg.width else NA
                    , silhouette.basis = if( !is_NA(siF) ) summary(siF)$avg.width else NA)
			
			# return result
			return(res)
		}
)


#' Sparseness 
#' 
#' Generic function that computes the \emph{sparseness} of an object, as defined 
#' by \cite{Hoyer2004}. 
#' The sparseness quantifies how much energy of a vector is packed into only few components.
#'   
#' In \cite{Hoyer2004}, the sparseness is defined for a real vector \eqn{x} as: 
#' \deqn{Sparseness(x) = \frac{\sqrt{n} - \frac{\sum |x_i|}{\sqrt{\sum x_i^2}}}{\sqrt{n}-1}}{
#' (srqt(n) - ||x||_1 / ||x||_2) / (sqrt(n) - 1)}
#' 
#' , where \eqn{n} is the length of \eqn{x}.
#' 
#' The sparseness is a real number in \eqn{[0,1]}. 
#' It is equal to 1 if and only if \code{x} contains a single nonzero component, 
#' and is equal to 0 if and only if all components of \code{x} are equal.
#' It interpolates smoothly between these two extreme values.  
#' The closer to 1 is the sparseness the sparser is the vector.
#' 
#' The basic definition is for a \code{numeric} vector, and is extended for matrices as the 
#' mean sparseness of its column vectors.
#' 
#' @param x an object whose sparseness is computed.
#' @param ... extra arguments to allow extension
#' 
#' @return usually a single numeric value -- in [0,1], or a numeric vector. 
#' See each method for more details.
#' 
#' @export
#' @family assess  
setGeneric('sparseness', function(x, ...) standardGeneric('sparseness') )
#' Base method that computes the sparseness of a numeric vector.
#' 
#' It returns a single numeric value, computed following the definition 
#' given in section \emph{Description}. 
setMethod('sparseness', signature(x='numeric'), 
	function(x){
		# get length of x
		n <- length(x)
		# compute and return the sparseness
		( sqrt(n) - sum(abs(x)) / sqrt(sum(x^2)) ) / (sqrt(n)-1)
	}
)
#' Computes the sparseness of a matrix as the mean sparseness of its column vectors.
#' It returns a single numeric value.
setMethod('sparseness', signature(x='matrix'), 
	function(x){
		# compute the sparseness of each column
		s <- apply(x, 2, sparseness)
		
		# return the mean sparseness
		mean(s)
	}
)
#' Compute the sparseness of an object of class \code{NMF}, as the sparseness of 
#' the basis and coefficient matrices computed separately.
#' 
#' It returns the two values in a numeric vector with names \sQuote{basis} and \sQuote{coef}. 
setMethod('sparseness', signature(x='NMF'), 
	function(x){		
		# return the sparseness of the basis and coef matrix
		c(basis=sparseness(basis(x)), coef=sparseness(coef(x)))
	}
)

#' Purity and Entropy of a Clustering
#' 
#' The functions \code{purity} and \code{entropy} respectively compute the purity and the entropy 
#' of a clustering given \emph{a priori} known classes.
#' 
#' The purity and entropy measure the ability of a clustering method, to recover
#' known classes (e.g. one knows the true class labels of each sample), that are
#' applicable even when the number of cluster is different from the number of known classes.
#' \cite{KimH2007} used these measures to evaluate the performance of their alternate least-squares 
#' NMF algorithm.
#' 
#' @details
#' Suppose we are given \eqn{l} categories, while the clustering method generates 
#' \eqn{k} clusters.
#' 
#' The purity of the clustering with respect to the known categories is given by:
#' \deqn{Purity = \frac{1}{n} \sum_{q=1}^k \max_{1 \leq j \leq l} n_q^j} ,
#' 
#' where:
#' \itemize{
#' \item \eqn{n} is the total number of samples;
#' \item \eqn{n_q^j} is the number of samples in cluster \eqn{q} that belongs to
#' original class \eqn{j} (\eqn{1 \leq j \leq l}).
#' }
#' 
#' The purity is therefore a real number in \eqn{[0,1]}.  
#' The larger the purity, the better the clustering performance.
#'   
#' @param x an object that can be interpreted as a factor or can generate such an object, e.g. via
#' a suitable method \code{\link{predict}}, which gives the cluster membership for each sample.
#' @param y a factor or an object coerced into a factor that gives the true class labels for each sample.
#' It may be missing if \code{x} is a contingency table. 
#' @param ... extra arguments to allow extension, and usually passed to the next method. 
#' 
#' @return a single numeric value
#' @family assess
#' @export
#' 
#' @examples
#' # generate a synthetic dataset with known classes: 50 features, 18 samples (5+5+8)
#' n <- 50; counts <- c(5, 5, 8);
#' V <- syntheticNMF(n, counts)
#' cl <- unlist(mapply(rep, 1:3, counts))
#' 
#' # perform default NMF with rank=2
#' x2 <- nmf(V, 2)
#' purity(x2, cl)
#' entropy(x2, cl)
#' # perform default NMF with rank=2
#' x3 <- nmf(V, 3)
#' purity(x3, cl)
#' entropy(x3, cl)
#' 
setGeneric('purity', function(x, y, ...) standardGeneric('purity') )
#' Computes the purity directly from the contingency table \code{x}
setMethod('purity', signature(x='table', y='missing'), 
	function(x, y){
		#for each cluster: compute maximum number of samples common to a class
		t <- apply(x, 1, max)
		# average and return the result
		sum(t) / sum(x)
	}
)
#' Computes the purity on the contingency table of \code{x} and \code{y}, that is  
#' coerced into a factor if necessary.
setMethod('purity', 'factor', 
	function(x, y, ...){
		
		# coerce `y` into a factor if necessary
		if( !is.factor(y) ) y <- as.factor(y)
		#compute the purity on the contingency table between clusters and true classes (clusters are in rows)
		purity(table(x, y), ...)
	}
)
#' Default method that should work for results of clustering algorithms, that have a 
#' suitable \code{predict} method that returns the cluster membership vector:
#' the purity is computed between \code{x} and \code{predict{y}}
setMethod('purity', 'ANY', 
	function(x, y, ...){
		# compute the purity for the samples clusters defined by the profiles
		purity(predict(x), y, ...)
	}
)

#' Entropy of a Clustering 
#' 
#' @details
#' The entropy of the clustering with respect to the known categories is given by:
#' \deqn{Entropy = - \frac{1}{n \log_2 l} \sum_{q=1}^k \sum_{j=1}^l n_q^j
#' \log_2 \frac{n_q^j}{n_q}}{
#' - 1/(n log2(l) ) sum_q sum_j n(q,j) log2( n(q,j) / n_q )},
#' 
#' where:
#' \itemize{
#' \item \eqn{n} is the total number of samples;
#' \item \eqn{n}{n_q} is the total number of samples in cluster \eqn{q} (\eqn{1 \leq q \leq k});
#' \item \eqn{n_q^j}{n(q,j)} is the number of samples in cluster \eqn{q} that belongs to
#' original class \eqn{j} (\eqn{1 \leq j \leq l}).
#' }
#' 
#' The smaller the entropy, the better the clustering performance.
#' @inheritParams purity
#' 
#' @return the entropy (i.e. a single numeric value)
#' @family assess
#' @rdname purity 
#' @export 
#' 
setGeneric('entropy', function(x, y, ...) standardGeneric('entropy') )
#' Computes the purity directly from the contingency table \code{x}.
#' 
#' This is the workhorse method that is eventually called by all other methods.
setMethod('entropy', signature(x='table', y='missing'), 
	function(x, y, ...){
		#for each cluster: compute the inner sum
		t <- apply(x, 1, function(n){ c.size <- sum(n); n %*% ifelse( n!=0, log2(n/c.size), 0)} )
		
		# weight and return the result
		- sum(t) / ( sum(x) * log2(ncol(x)) )
	}
)
#' Computes the purity on the contingency table of \code{x} and \code{y}, that is 
#' coerced into a factor if necessary.
setMethod('entropy', 'factor', 
	function(x, y, ...){
		# coerce `y` into a factor if necessary
		if( !is.factor(y) ) y <- as.factor(y)
		#copmute entropy on contingency table between clusters and true classes (clusters are in rows)
		entropy(table(x, y))
	}
)
#' Default method that should work for results of clustering algorithms, that have a 
#' suitable \code{predict} method that returns the cluster membership vector:
#' the purity is computed between \code{x} and \code{predict{y}}
setMethod('entropy', 'ANY', 
	function(x, y, ...){
		# compute the entropy for the samples clusters defined by the metagenes expression matrix
		entropy(predict(x), y)
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


#' Apply Function for NMF Objects
#' 
#' The function \code{nmfApply} provides exteneded \code{apply}-like 
#' functionality for objects of class \code{NMF}.
#' It enables to easily apply a function over different margins of
#' NMF models.
#' 
#' The function \code{FUN} is applied via a call to \code{\link{apply}} 
#' or \code{\link{sapply}} according to the value of argument \code{MARGIN} 
#' as follows:
#' 
#' \describe{
#' \item{MARGIN=1}{ apply \code{FUN} to each \emph{row} of the basis matrix:
#' \code{apply(basis(X), 1L, FUN, ...)}.}
#'  
#' \item{MARGIN=2}{ apply \code{FUN} to each \emph{column} of the coefficient matrix:
#' \code{apply(coef(X), 2L, FUN, ...)}.}
#' 
#' \item{MARGIN=3}{ apply \code{FUN} to each \emph{pair} of associated basis component 
#' and basis profile:
#' more or less \code{sapply(seq(nbasis(X)), function(i, ...) FUN(basis(X)[,i], coef(X)[i, ], ...), ...)}.
#' 
#' In this case \code{FUN} must be have at least two arguments, to which are passed 
#' each basis components and basis profiles respectively -- as numeric vectors.}
#' 
#' \item{MARGIN=4}{ apply \code{FUN} to each \emph{column} of the basis matrix, i.e. to each 
#' basis component:
#' \code{apply(basis(X), 2L, FUN, ...)}.}
#' 
#' \item{MARGIN=5}{ apply \code{FUN} to each \emph{row} of the coefficient matrix:
#' \code{apply(coef(X), 1L, FUN, ...)}.}
#' 
#' }
#' 
#' 
#' @param X an object that has suitable \code{\link{basis}} and \code{coef} methods, 
#' e.g. an NMF model.
#' @param MARGIN a single numeric (integer) value that specifies over which margin(s)
#' the function \code{FUN} is applied.
#' See section \emph{Details} for a list of possible values.
#' @param FUN a function to apply over the specified margins.
#' @param ... extra arguments passed to \code{FUN}
#' @param simplify a logical only used when \code{MARGIN=3}, that indicates if \code{sapply} 
#' should try to simplify result if possible.
#' Since this argument follows \sQuote{...} its name cannot be abbreviated.
#' @param USE.NAMES a logical only used when \code{MARGIN=3}, that indicates if \code{sapply} 
#' should use the names of the basis components to name the results if present.
#' Since this argument follows \sQuote{...} its name cannot be abbreviated.
#'  
#' @return a vector or a list. 
#' See \code{\link[base]{apply}} and \code{\link[base]{sapply}} for more details on 
#' the output format.
#' 
#' @export  
#setGeneric('nmfApply', function(object, ...) standardGeneric('nmfApply') )
nmfApply <- function(X, MARGIN, FUN, ..., simplify = TRUE, USE.NAMES = TRUE){
	if( MARGIN == 1L )
		apply(basis(X), 1L, FUN, ...)
	else if( MARGIN == 4L )
		apply(basis(X), 2L, FUN, ...)
	else if( MARGIN == 2L )
		apply(coef(X), 2L, FUN, ...)
	else if( MARGIN == 5L )
		apply(coef(X), 1L, FUN, ...)
	else if( MARGIN == 3L ){
		b <- basis(X)
		p <- coef(X)
		sapply(setNames(seq(nbasis(X), basisnames(X)))
				, function(i, ...) FUN(b[,i], p[i,], ...)
				, simplify = simplify, USE.NAMES = USE.NAMES)
	}else stop("invalid argument 'MARGIN' (expected values are: 1-basis rows, 2-coef columns, 3-(basis columns, coef rows), or 4-basis columns or 5-coef rows)")
}

###% Utility function to compute the dominant column for each row for a matrix.
.predict.nmf <- function(x, prob=FALSE){
	
	if( !is.matrix(x) ) stop('NMF:::.predict.nmf : only works on matrices')
	if( !prob ){
		#for each column return the (row) index of the maximum
		return( as.factor(apply(x, 1L, function(v) which.max(abs(v)))) )
	}
	else{
		#for each column return the (row) index of the maximum AND the associated probaility
		res <- apply(x, 1L,
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

#' Clustering and Prediction
#' 
#' The methods \code{predict} for NMF models return the cluster membership 
#' of each sample or each feature.
#' Currently the classification/prediction of new data is not implemented. 
#' 
#' The cluster membership is computed as the index of the dominant basis 
#' component for each sample (\code{what='samples' or 'columns'}) or each feature
#' (\code{what='features' or 'rows'}), based on their corresponding
#' entries in the coefficient matrix or basis matrix respectively.
#' 
#' For example, if \code{what='samples'}, then the dominant basis component 
#' is computed for each column of the coefficient matrix as the row index
#' of the maximum within the column.
#' 
#' If argument \code{prob=FALSE} (default), the result is a \code{factor}.
#' Otherwise a list with two elements is returned: element \code{predict}
#' contains the cluster membership index (as a \code{factor}) and element 
#' \code{prob} contains the relative contribution of
#' the dominant component to each sample (resp. the relative contribution of 
#' each feature to the dominant basis component):
#' 
#' \itemize{
#' \item Samples: \deqn{p_j = x_{k_0} / \sum_k x_k}{p(j) = x(k0) / sum_k x(k)}, 
#' for each sample \eqn{1\leq j \leq p}, where \eqn{x_k}{x(k)} is the contribution 
#' of the \eqn{k}-th basis component to \eqn{j}-th sample (i.e. \code{H[k ,j]}), and 
#' \eqn{x_{k_0}}{x(k0)} is the maximum of these contributions.
#' 
#' \item Features: \deqn{p_i = y_{k_0} / \sum_k y_k}{p(i) = y(k0) / sum_k y(k)}, 
#' for each feature \eqn{1\leq i \leq p}, where \eqn{y_k}{y(k)} is the contribution 
#' of the \eqn{k}-th basis component to \eqn{i}-th feature (i.e. \code{W[i, k]}), and 
#' \eqn{y_{k_0}}{y(k0)} is the maximum of these contributions.
#' 
#' }
#' 
#' @param object an NMF model
#'  
#' @family stats Methods for the Interface Defined in Package stats
#' 
#' @cite Brunet2004,Pascual-Montano2006
#' @inline
#' @export
setGeneric('predict', package='stats')

#' Default method for NMF models
#' 
#' @param what a character string that indicates the type of cluster membership should
#' be returned: \sQuote{columns} or \sQuote{rows} for clustering the colmuns or the 
#' rows of the target matrix respectively.
#' The values \sQuote{samples} and \sQuote{features} are aliases for \sQuote{colmuns} 
#' and \sQuote{rows} respectively.
#' @param prob logical that indicates if the relative contributions of/to the dominant 
#' basis component should be computed and returned. See \emph{Details}.
#' @param dmatrix logical that indicates if a dissimiliarity matrix should be 
#' attached to the result.
#' This is notably used internally when computing NMF clustering silhouettes.
#' 
#' @examples
#'
#' # random target matrix
#' v <- rmatrix(20, 10)
#' # fit an NMF model
#' x <- nmf(v, 5)
#'  
#' # predicted column and row clusters
#' predict(x)
#' predict(x, 'rows')
#' 
#' # with relative contributions of each basis component
#' predict(x, prob=TRUE)
#' predict(x, 'rows', prob=TRUE)
#' 
setMethod('predict', 'NMF',
		function(object, what=c('columns', 'rows', 'samples', 'features'), prob=FALSE, dmatrix = FALSE){
			# determine which matrix to use for the prediction
			what <- match.arg(what)
			x <- if( what %in% c('features', 'rows') ) basis(object, all=FALSE) else t(coef(object, all=FALSE))
			
			# compute the indice of the dominant row for each column
            res <- .predict.nmf(x, prob)
            # attach dissimilarity matrix if requested
            if( dmatrix ){
                attr(res, 'dmatrix') <- 1 - cor(t(x))
            }
			return( res )
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

#' Correlations in NMF Models
#' 
#' \code{basiscor} computes the correlation matrix between basis vectors, i.e. 
#' the \emph{columns} of its basis matrix -- which is the model's first matrix factor.
#' 
#' @details
#' Each generic has methods defined for computing correlations between NMF models 
#' and/or compatible matrices.
#' The computation is performed by the base function \code{\link{cor}}.
#' 
#' @param x a matrix or an object with suitable methods \code{\link{basis}} 
#' or \code{\link{coef}}.
#' @param y a matrix or an object with suitable methods \code{\link{basis}} 
#' or \code{\link{coef}}, and dimensions compatible with \code{x}.
#' If missing the correlations are computed between \code{x} and \code{y=x}. 
#' @param ... extra arguments passed to \code{\link{cor}}.
#' 
#' @export
#' @family NMFplots Plotting functions for NMF objects 
#' 
#' @examples 
#' 
#' # generate two random NMF models
#' a <- rnmf(3, 100, 20)
#' b <- rnmf(3, 100, 20)
#' 
#' # Compute auto-correlations
#' basiscor(a)
#' profcor(a)
#' # Compute correlations with b
#' basiscor(a, b)
#' profcor(a, b)
#' 
#' # try to recover the underlying NMF model 'a' from noisy data 
#' res <- nmf(fitted(a) + rmatrix(a), 3)
#' 
#' # Compute correlations with the true model
#' basiscor(a, res)
#' profcor(a, res)
#' 
#' # Compute correlations with a random compatible matrix
#' W <- rmatrix(basis(a))
#' basiscor(a, W)
#' identical(basiscor(a, W), basiscor(W, a))
#' 
#' H <- rmatrix(coef(a))
#' profcor(a, H)
#' identical(profcor(a, H), profcor(H, a))
#' 
setGeneric('basiscor', function(x, y, ...) standardGeneric('basiscor') )
#' Computes the correlations between the basis vectors of \code{x} and 
#' the columns of \code{y}. 
setMethod('basiscor', signature(x='NMF', y='matrix'),
	function(x, y, ...){
		cor(basis(x), y, ...)
	}
)
#' Computes the correlations between the columns of \code{x} 
#' and the the basis vectors of \code{y}.  
setMethod('basiscor', signature(x='matrix', y='NMF'),
	function(x, y, ...){
		cor(x, basis(y), ...)
	}
)
#' Computes the correlations between the basis vectors of \code{x} and \code{y}.
setMethod('basiscor', signature(x='NMF', y='NMF'),
		function(x, y, ...){
			basiscor(x, basis(y), ...)
		}
)
#' Computes the correlations between the basis vectors of \code{x}.
setMethod('basiscor', signature(x='NMF', y='missing'),
		function(x, y, ...){
			basiscor(x, x, ...)
		}
)

#' Correlations of Basis Profiles
#' 
#' \code{profcor} computes the correlation matrix between basis profiles, 
#' i.e. the \emph{rows} of the coefficient matrix -- which is the model's second 
#' matrix factor.
#' 
#' @rdname basiscor
#' @export
#' 
setGeneric('profcor', function(x, y, ...) standardGeneric('profcor') )
#' Computes the correlations between the basis profiles of \code{x} and 
#' the rows of \code{y}.
setMethod('profcor', signature(x='NMF', y='matrix'),
		function(x, y, ...){
			cor(t(coef(x)), t(y), ...)
		}
)
#' Computes the correlations between the rows of \code{x} and the basis 
#' profiles of \code{y}.
setMethod('profcor', signature(x='matrix', y='NMF'),
		function(x, y, ...){
			cor(t(x), t(coef(y)), ...)
		}
)
#' Computes the correlations between the basis profiles of \code{x} and \code{y}.
setMethod('profcor', signature(x='NMF', y='NMF'),
		function(x, y, ...){
			profcor(x, coef(y), ...)
		}
)
#' Computes the correlations between the basis profiles of \code{x}.
setMethod('profcor', signature(x='NMF', y='missing'),
		function(x, y, ...){
			profcor(x, x, ...)
		}
)

#' Clustering Connectivity and Consensus Matrices
#' 
#' \code{connectivity} is an S4 generic that computes the connectivity matrix 
#' based on the clustering of samples obtained from a model's \code{\link{predict}} 
#' method.
#' 
#' The connectivity matrix of a given partition of a set of samples (e.g. given 
#' as  a cluster membership index) is the matrix \eqn{C} containing only 0 or 1 
#' entries such that: 
#' \deqn{C_{ij} = \left\{\begin{array}{l}
#' 1\mbox{ if sample }i\mbox{ belongs to the same cluster as sample }j\\
#' 0\mbox{ otherwise}
#' \end{array}\right..}{
#' C_{ij} = 1 if sample i belongs to the same cluster as sample j, 0 otherwise}
#' 
#' @param object an object with a suitable \code{\link{predict}} method.
#' @param ... extra arguments to allow extension.
#' They are passed to \code{\link{predict}}, except for the \code{vector} and 
#' \code{factor} methods.
#'  
#' @return a square matrix of dimension the number of samples in the model, full 
#' of 0s or 1s.   
#' 
#' @seealso \code{\link{predict}}
#' 
#' @export
#' 
setGeneric('connectivity', function(object, ...) standardGeneric('connectivity') )
#' Default method which computes the connectivity matrix 
#' using the result of \code{predict(x, ...)} as cluster membership index.
#' 
#' @examples
#' 
#' # clustering of random data
#' h <- hclust(dist(rmatrix(10,20)))
#' connectivity(cutree(h, 2))
#' 
setMethod('connectivity', 'ANY', 
	function(object, ...){
		c <- predict(object, ...);
		outer(c, c, function(x,y) ifelse(x==y, 1,0));
	}
)

#' Computes the connectivity matrix using \code{x} as cluster membership index.
#' 
#' @examples
#' connectivity(gl(2, 4))
#' 
setMethod('connectivity', 'factor', 
	function(object, ...){
		outer(object, object, function(x,y) ifelse(x==y, 1,0));
	}
)
#' Equivalent to \code{connectivity(as.factor(x))}. 
setMethod('connectivity', 'numeric', 
	function(object, ...){
		connectivity(as.factor(object), ...)
	}
)
#' Computes the connectivity matrix for an NMF model, for which cluster 
#' membership is given by the most contributing basis component in each sample.
#' See \code{\link{predict,NMF-method}}.
#'
#' @inline 
#' @param no.attrib a logical that indicates if attributes containing information 
#' about the NMF model should be attached to the result (\code{TRUE}) or not
#' (\code{FALSE}). 
#' 
setMethod('connectivity', 'NMF', 
	function(object, no.attrib=FALSE){
		C <- callNextMethod(object=object, what='samples');
		if( !no.attrib ){
			class(C) <- c(class(C), 'NMF.consensus')
			attr(C, 'model') <- object
			attr(C, 'nrun') <- 1
			attr(C, 'nbasis') <- nbasis(object)
		}
		C
	}
)

# Unit test
unit.test(connectivity,{
	
	# build reference matrix
	n <- 10
	ref <- matrix(0, 2*n, 2*n)
	ref[1:n,1:n] <- 1
	ref[(n+1):(2*n),(n+1):(2*n)] <- 1
	
	checkIdentical(connectivity(gl(2, n)), ref, 'Factor')
	checkIdentical(connectivity(as.numeric(gl(2, n))), ref, 'Vector')
	# test with NMF model
	i <- gl(2, n)
	x <- nmfModel(H=matrix(c(rev(i), i), 2, byrow=TRUE))
	checkEquals(connectivity(x), ref, 'NMF model', check.attributes = FALSE)
	s <- sample.int(2*n)
	checkEquals(connectivity(x[,s]), ref[s,s], 'NMF model (shuffled)', check.attributes = FALSE)
	
})

#' Residual Sum of Squares and Explained Variance
#' 
#' \code{rss} and \code{evar} are S4 generic functions that respectively computes 
#' the Residual Sum of Squares (RSS) and explained variance achieved by a model.
#' 
#' @param object an R object with a suitable \code{\link{fitted}}, \code{rss} or
#' \code{evar} method.
#' @param ... extra arguments to allow extension, e.g. passed to \code{rss} 
#' in \code{evar} calls.
#' 
#' @return a single numeric value 
#' @inline
#' @export
#' 
setGeneric('rss', function(object, ...) standardGeneric('rss'))
#' Computes the RSS between a target matrix and its estimate \code{object}, 
#' which must be a matrix of the same dimensions as \code{target}.
#' 
#' The RSS between a target matrix \eqn{V} and its estimate \eqn{v} is computed as:
#' \deqn{RSS = \sum_{i,j} (v_{ij} - V_{ij})^2}
#' 
#' Internally, the computation is performed using an optimised C++ implementation, 
#' that is light in memory usage.
#' 
#' @param target target matrix
#' 
#' @examples
#' # RSS bewteeen random matrices
#' x <- rmatrix(20,10, max=50)
#' y <- rmatrix(20,10, max=50)
#' rss(x, y)
#' rss(x, x + rmatrix(x, max=0.1))
#' 
setMethod('rss', 'matrix', 
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
		.rss(object,target)
	}
)
#' Residual sum of square between a given target matrix and a model that has a 
#' suitable \code{\link{fitted}} method.
#' It is equivalent to \code{rss(fitted(object), ...)}
#' 
#' In the context of NMF, \cite{Hutchins2008} used the variation of the RSS 
#' in combination with the algorithm from \cite{Lee1999} to estimate the 
#' correct number of basis vectors. 
#' The optimal rank is chosen where the graph of the RSS first shows an inflexion
#' point, i.e. using a screeplot-type criterium.
#' See section \emph{Rank estimation} in \code{\link{nmf}}. 
#' 
#' Note that this way of estimation may not be suitable for all models. 
#' Indeed, if the NMF optimisation problem is not based on the Frobenius norm, 
#' the RSS is not directly linked to the quality of approximation of the NMF model.
#' However, it is often the case that it still decreases with the rank. 
#' 
#' @examples 
#' # RSS between an NMF model and a target matrix
#' x <- rmatrix(20, 10)
#' y <- rnmf(3, x) # random compatible model
#' rss(y, x)
#' 
#' # fit a model with nmf(): one should do better
#' y2 <- nmf(x, 3) # default minimizes the KL-divergence
#' rss(y2, x)
#' y2 <- nmf(x, 3, 'lee') # 'lee' minimizes the RSS
#' rss(y2, x)
#' 
setMethod('rss', 'ANY', 
	function(object, ...){
		rss(fitted(object), ...)
	}
)


unit.test(rss, {
	
	x <- rmatrix(20,10, max=50)
	y <- rmatrix(20,10, max=50)
	checkIdentical(rss(x, y), sum((x-y)^2), "Random matrices")
	
	y <- rnmf(3, x) # random compatible model
	r1 <- rss(y, x)
	checkIdentical(r, sum((x-fitted(y))^2), 'NMF model')
	checkIdentical(rss(y, ExpressionSet(x)), sum((x-fitted(y))^2), 'NMF model (ExpressionSet)')
	y <- nmf(x, 3)
	r2 <- rss(y, x)
	checkIdentical(r2, sum((x-fitted(y))^2), 'Fitted NMF model')
	checkTrue(r2 < r1, 'Fitted NMF model has better RSS')
	y <- nmf(x, 3, 'lee')
	checkTrue(rss(y, x) < r2, "Fitted NMF model with 'lee' has better RSS than 'brunet'")
	
})

#' Explained Variance
#' 
#' The explained variance for a target \eqn{V} is computed as:
#' \deqn{evar = 1 - \frac{RSS}{\sum_{i,j} v_{ij}^2} }{evar = 1 - RSS/sum v_{ij}^2},
#' 
#' where RSS is the residual sum of squares.
#' 
#' The explained variance is usefull to compare the performance of different 
#' models and their ability to accurately reproduce the original target matrix. 
#' Note, however, that a possible caveat is that some models explicitly aim at 
#' minimizing the RSS (i.e. maximizing the explained variance), while others do not.
#' 
#' @rdname rss
#' @inline
#' @export
#' 
setGeneric('evar', function(object, ...) standardGeneric('evar'))
#' Default method for \code{evar}.
#' 
#' It requires a suitable \code{rss} method to be defined 
#' for \code{object}, as it internally calls \code{rss(object, target, ...)}. 
setMethod('evar', 'ANY', 
	function(object, target, ...){
		
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
		1 - rss(object, target, ...) / sum(t^2)
	}
)

#' Distances and Objective Functions 
#' 
#' The NMF package defines methods for the generic \code{deviance} from the package \code{stats}, 
#' to compute approximation errors between NMF models and matrices, using a variety of 
#' objective functions.
#' 
#' @return \code{deviance} returns a nonnegative numerical value 
#' @family stats
#' 
#' @export 
setGeneric('deviance', package='stats')
#' Computes the distance between a matrix and the estimate of an \code{NMF} model.
#'  
#' @param y a matrix compatible with the NMF model \code{object}, i.e. \code{y} 
#' must have the same dimension as \code{fitted(object)}.
#' @param method a character string or a function with signature 
#' \code{(x="NMF", y="matrix", ...)} that implements a distance measure between 
#' an NMF model \code{x} and a target matrix \code{y}, i.e. an objective function
#' to use to compute the deviance. 
#' In \code{deviance}, it is passed to \code{nmfDistance} to get the function 
#' that effectively computes the deviance.
#' @param ... extra parameters passed to the objective function.
#' 
#' @inline 
#' @family stats
#' 
setMethod('deviance', 'NMF', 
	function(object, y, method=c('', 'KL', 'euclidean'), ...){

	fun <- nmfDistance(method)
	
	if( is.null(fun) ){
		warning('Undefined distance method: distance cannot be computed [returned NA]')
		return(as.numeric(NA))
	}
	
	# extract expression data from ExpressionSet objects
	if( is(y, 'ExpressionSet') )
		y <- Biobase::exprs(y)
	
	# apply the function and return the result
	fun(object, y, ...)

	}
)

#' \code{nmfDistance} returns a function that computes the distance between an NMF model and a 
#' compatible matrix.
#' 
#' @return \code{nmfDistance} returns a function with least two arguments: 
#' an NMF model and a matrix.  
#' 
#' @export 
#' @rdname deviance 
nmfDistance <- function(method=c('', 'KL', 'euclidean')){
	
		#message('compute distance')
		# determinate the distance measure to use
		if( is.null(method) ) return(NULL)
		
		if( is.character(method) ){
			errMeth <- try(method <- match.arg(method), silent=TRUE)
			# if the method is not predefined, try to find a function with the given name
			if( inherits(errMeth, 'try-error') ){			
				#TODO: this is not working with local functions
				if( is.character(method) ){
					errFun <- try(fun <- match.fun(method), silent=TRUE)
					if( inherits(errFun, 'try-error') ) 
						stop("Could not find distance measure '", method, "':\n\t- not a predefined measures -> ", errMeth,"\t- not a function -> ", errFun)
				}
				else fun <- method
				
				if( !is.function(fun) )
					stop('Invalid distance measure: should be a character string or a valid function definition')
			}
			else{
				# compute and return the distance measure		
				fun <- switch(method,
						euclidean = function(x, y, ...){
							# call optimized C function
							.rss(y, fitted(x))/2							
						},
						KL = function(x, y, ...){							
							# call optimized C function
							.KL(y, fitted(x))					
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
		


#' Testing Equality of NMF Models
#'  
#' The function \code{nmf.equal} tests if two NMF models are the same, i.e. they
#' contain -- almost -- identical data: same basis and coefficient matrices, as 
#' well as same extra parameters.
#' 
#' @details
#' \code{nmf.equal} compares two NMF models, and return \code{TRUE} iff they are 
#' identical acording to the function \code{\link{identical}} when \code{identical=TRUE}, 
#' or equal up to some tolerance acording to the function \code{\link{all.equal}}.
#' This means that all data contained in the objects are compared, which includes 
#' at least the basis and coefficient matrices, as well as the extra parameters
#' stored in slot \sQuote{misc}. 
#' 
#' If extra arguments are specified in \code{...}, then the comparison is performed
#' using \code{\link{all.equal}}, irrespective of the value of argument \code{identical}.    
#' 
#' @param x an NMF model or an object that is associated with an NMF model, e.g. 
#' the result from a fit with \code{\link{nmf}}.
#' @param y an NMF model or an object that is associated with an NMF model, e.g. 
#' the result from a fit with \code{\link{nmf}}.
#' @param identical a logical that indicates if the comparison should be made
#' using the function \code{\link{identical}} (\code{TRUE}) or \code{\link{all.equal}}
#' (\code{FALSE}). See description for method \code{nmf.equal,NMF,NMF}.
#' @param ... extra arguments to allow extension, and passed to subsequent calls
#'   
#' @export
#' 
setGeneric('nmf.equal', function(x, y, ...) standardGeneric('nmf.equal') )
#' Compares two NMF models.
#' 
#' Arguments in \code{...} are used only when \code{identical=FALSE} and are 
#' passed to \code{all.equal}.
#' @inline
setMethod('nmf.equal', signature(x='NMF', y='NMF'), 
		function(x, y, identical=TRUE, ...){
			
			dots <- list(...)
			if( identical && length(dots) == 0 )
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
