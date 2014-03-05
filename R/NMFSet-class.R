#' @include NMFfit-class.R
#' @include heatmaps.R
NULL

#' \code{isNMFfit} tells if an object results from an NMF fit. 
#' 
#' @details  \emph{isNMFfit} checks if \code{object} inherits from class 
#' \code{\linkS4class{NMFfit}} or \code{\linkS4class{NMFfitX}}, which are 
#' the two types of objects returned by the function \code{\link{nmf}}. 
#' If \code{object} is a plain \code{list} and \code{recursive=TRUE}, then 
#' the test is performed on each element of the list, and the return value 
#' is a logical vector (or a list if \code{object} is a list of list) of 
#' the same length as \code{object}.
#' 
#' @export
#' @rdname types
#' @param object any R object.
#' @param recursive if \code{TRUE} and \code{object} is a plain list then 
#' \code{isNMFfit} tests each element of the list. 
#' Note that the recursive test only applies in the case of lists that are 
#' not themselves NMFfit objects, like \code{NMFfitXn} objects for which 
#' the result of \code{isNMFfit} will always be \code{TRUE}, although they are 
#' list objects (a single logical value).
#' 
#' @return \code{isNMFfit} returns a \code{logical} vector (or a list if
#' \code{object} is a list of list) of the same length as \code{object}.
#' 
#' @seealso \code{\linkS4class{NMFfit}}, \code{\linkS4class{NMFfitX}}, 
#' \code{\linkS4class{NMFfitXn}}
#' @examples
#' 
#' ## Testing results of fits
#' # generate a random
#' V <- rmatrix(20, 10)
#' 
#' # single run -- using very low value for maxIter to speed up the example 
#' res <- nmf(V, 3, maxIter=3L)
#' isNMFfit(res)
#' 
#' # multiple runs - keeping single fit
#' resm <- nmf(V, 3, nrun=3, maxIter=3L)
#' isNMFfit(resm)
#' 
#' # multiple runs - keeping all fits
#' resM <- nmf(V, 3, nrun=3, .opt='k', maxIter=3L) 
#' isNMFfit(resM)
#' 
#' # with a list of results
#' isNMFfit(list(res, resm, resM, 'not a result'))
#' isNMFfit(list(res, list(resm, resM), 'not a result')) # list of list
#' isNMFfit(list(res, resm, resM, 'not a result'), recursive=FALSE)
#' 	 
isNMFfit <- function(object, recursive=TRUE){
	res <- is(object, 'NMFfit') || is(object, 'NMFfitX')
	# if the object is not a NMF result: apply to each element if a list (only in recursive mode)
	if( !res  && recursive && is.list(object) )
		sapply(object, isNMFfit)
	else
		res
}

#' Class for Storing Heterogeneous NMF fits
#' 
#' @description
#' This class wraps a list of NMF fit objects, which may come from different 
#' runs of the function \code{\link{nmf}}, using different parameters, methods, etc..
#' These can be either from a single run (NMFfit) or multiple runs (NMFfitX).
#' 
#' Note that its definition/interface is very likely to change in the future.
#' @export
#'   
setClass('NMFList'
		, representation(
			runtime='proc_time'
		)
		, contains='namedList'
		, validity=function(object){
			
			# the list must only contains NMFfit objects of the same dimensions
			ok <- isNMFfit(object)
			if( !is.logical(ok) )
				return("Could not validate elements in list: input is probably a complex structure of lists.")
			pb <- which(!ok)
			if( length(pb) ){
				return(paste("invalid class for element(s)"
							, str_out(i)
							, "of input list [all elements must be fitted NMF models]"))	
			}			
		}
)

#' Show method for objects of class \code{NMFList}
#' @export
setMethod('show', 'NMFList', 
	function(object)
	{
		cat("<Object of class:", class(object), ">\n")
		cat("Length:", length(object), "\n")
		if( length(object) > 0 ) cat("Method(s):", algorithm(object, string=TRUE), "\n")
		# show totaltime if present
		tt <- runtime(object)
		if( length(tt) > 0 ){
			cat("Total timing:\n"); show(tt);
		}
	}
)

#' Returns the method names used to compute the NMF fits in the list. 
#' It returns \code{NULL} if the list is empty.
#' 
#' @param string a logical that indicate whether the names should be collapsed
#' into a comma-separated string. 
#' @param unique a logical that indicates whether the result should contain the 
#' set of method names, removing duplicated names.
#' This argument is forced to \code{TRUE} when \code{string=TRUE}.
#' 
setMethod('algorithm', 'NMFList', 
	function(object, string=FALSE, unique=TRUE){
		l <- length(object)
		
		if( string ) unique <- TRUE
		
		if( l == 0 ) NULL
		else if( l == 1 ) algorithm(object[[1]])
		else{	
			# build the vector of the algorithm names (with no repeat)
			m <- sapply(object, algorithm)
			if( unique ) m <- unique(m)
			if( string ) m <- paste(m, collapse=', ')
			m
		}
	} 
)

.seqtime <- function(object){
	
	if( length(object) == 0 ) return(NULL)
	# sum up the time across the runs
	t.mat <- sapply(object, function(x){
		if( is(x, 'NMFfitXn') ) runtime.all(x)
		else runtime(x)
	})
	res <- rowSums(t.mat)
	class(res) <- 'proc_time'
	res
}

#' Returns the CPU time that would be required to sequentially compute all NMF 
#' fits stored in \code{object}.
#' 
#' This method calls the function \code{runtime} on each fit and sum up the 
#' results.
#' It returns \code{NULL} on an empty object.
setMethod('seqtime', 'NMFList', 
	function(object){	
		if( length(object) == 0 ) return(NULL)
		# sum up the time across the runs
		.seqtime(object)
	}
)

#' Returns the CPU time required to compute all NMF fits in the list.
#' It returns \code{NULL} if the list is empty.
#' If no timing data are available, the sequential time is returned. 
#' 
#' @param all logical that indicates if the CPU time of each fit should be 
#' returned (\code{TRUE}) or only the total CPU time used to compute all 
#' the fits in \code{object}.
setMethod('runtime', 'NMFList', 
	function(object, all=FALSE){
		if( !all ){
			t <- slot(object, 'runtime')
			if( length(t)==0 ) seqtime(object) else t
		}else
			sapply(object, runtime)
	}
)

as.NMFList <- function(..., unlist=FALSE){
	arg.l <- list(...)
	if( length(arg.l) == 1L && is.list(arg.l[[1]]) && !is(arg.l[[1]], 'NMFfitX') )
		arg.l <- arg.l[[1]]
	
	# unlist if required
	if( unlist )
		arg.l <- unlist(arg.l)
	
	# create a NMFList object from the input list 
	new('NMFList', arg.l)
}
			

#' Virtual Class to Handle Results from Multiple Runs of NMF Algorithms
#' 
#' This class defines a common interface to handle the results from multiple
#' runs of a single NMF algorithm, performed with the \code{\link{nmf}} method.
#' 
#' Currently, this interface is implemented by two classes,
#' \code{\linkS4class{NMFfitX1}} and \code{\linkS4class{NMFfitXn}}, which
#' respectively handle the case where only the best fit is kept, and the case
#' where the list of all the fits is returned.
#' 
#' See \code{\link{nmf}} for more details on the method arguments.
#' 
#' @slot runtime.all Object of class \code{\link[=proc.time]{proc_time}} that 
#' contains CPU times required to perform all the runs.
#' 
#' @export
#' @family multipleNMF
#' @examples
#' 
#' # generate a synthetic dataset with known classes
#' n <- 20; counts <- c(5, 2, 3);
#' V <- syntheticNMF(n, counts)
#' 
#' # perform multiple runs of one algorithm (default is to keep only best fit)
#' res <- nmf(V, 3, nrun=3)
#' res
#' 
#' # plot a heatmap of the consensus matrix
#' \dontrun{ consensusmap(res) }
#' 
#' # perform multiple runs of one algorithm (keep all the fits)
#' res <- nmf(V, 3, nrun=3, .options='k')
#' res
#'  
setClass('NMFfitX'
		, representation(
				runtime.all = 'proc_time' # running time to perform all the NMF runs
		)
		, contains='VIRTUAL'		
)

#' Returns the CPU time required to compute all the NMF runs.
#' It returns \code{NULL} if no CPU data is available. 
setMethod('runtime.all', 'NMFfitX', 
		function(object){	
			t <- slot(object, 'runtime.all')
			if( length(t) > 0 ) t else NULL
		}
)

#' Returns the number of NMF runs performed to create \code{object}.
#' 
#' It is a pure virtual method defined to ensure \code{nrun} is defined 
#' for sub-classes of \code{NMFfitX}, which throws an error if called.
#' 
#' Note that because the \code{\link{nmf}} function allows to run the NMF 
#' computation keeping only the best fit, \code{nrun} may return a value 
#' greater than one, while only the result of the best run is stored in 
#' the object (cf. option \code{'k'} in method \code{\link{nmf}}).
setMethod('nrun', 'NMFfitX', 
	function(object){	
		stop("NMF::NMFfitX - missing definition for pure virtual method 'nrun' in class '", class(object), "'")
	}
)
#' This method always returns 1, since an \code{NMFfit} object is obtained 
#' from a single NMF run.  
setMethod('nrun', 'NMFfit', 
	function(object){
		1L
	}
)

#' \code{consensus} is an S4 generic that computes/returns the consensus matrix 
#' from a model object, which is the mean connectivity matrix of all the runs.
#' 
#' The consensus matrix has been proposed by \cite{Brunet2004} to help 
#' visualising and measuring the stability of the clusters obtained by 
#' NMF approaches.
#' For objects of class \code{NMF} (e.g. results of a single NMF run, or NMF
#' models), the consensus matrix reduces to the connectivity matrix.
#' 
#' @rdname connectivity
#' @export
setGeneric('consensus', function(object, ...) standardGeneric('consensus') )
#' Pure virtual method defined to ensure \code{consensus} is defined for sub-classes of \code{NMFfitX}.
#' It throws an error if called.
setMethod('consensus', 'NMFfitX', 
	function(object, ...){	
		stop("NMF::NMFfitX - missing definition for pure virtual method 'consensus' in class '", class(object), "'")
	}
)

#' This method is provided for completeness and is identical to 
#' \code{\link{connectivity}}, and returns the connectivity matrix, 
#' which, in the case of a single NMF model, is also the consensus matrix.
setMethod('consensus', 'NMF', 
	function(object, ...){
		connectivity(object, ...)
	}
)

#' Hierarchical Clustering of a Consensus Matrix
#' 
#' The function \code{consensushc} computes the hierarchical clustering of 
#' a consensus matrix, using the matrix itself as a similarity matrix and 
#' average linkage.
#' It is 
#' 
#' @param object a matrix or an \code{NMFfitX} object, as returned by multiple 
#' NMF runs. 
#' @param ... extra arguments passed to next method calls
#' 
#' @return an object of class \code{dendrogram} or \code{hclust} depending on the  
#' value of argument \code{dendrogram}.   
#' 
#' @inline
#' @export
setGeneric('consensushc', function(object, ...) standardGeneric('consensushc'))
#' Workhorse method for matrices.
#' 
#' @param method linkage method passed to \code{\link{hclust}}.
#' @param dendrogram a logical that specifies if the result of the hierarchical 
#' clustering (en \code{hclust} object) should be converted into a dendrogram. 
#' Default value is \code{TRUE}.
setMethod('consensushc', 'matrix', 
	function(object, method='average', dendrogram=TRUE){
		
		# hierachical clustering based on the connectivity matrix
		hc <- hclust(as.dist(1-object), method=method)
		
		# convert into a dendrogram if requested
		if( dendrogram ) as.dendrogram(hc)
        else hc
	}
)
#' Compute the hierarchical clustering on the connectivity matrix of \code{object}.
setMethod('consensushc', 'NMF', 
	function(object, ...){		
		# hierachical clustering based on the connectivity matrix
		consensushc(connectivity(object), ...)		
	}
)
#' Compute the hierarchical clustering on the consensus matrix of \code{object}, 
#' or on the connectivity matrix of the best fit in \code{object}.
#' 
#' @param what character string that indicates which matrix to use in the 
#' computation. 
#' 
setMethod('consensushc', 'NMFfitX', 
	function(object, what=c('consensus', 'fit'), ...){
		
		what <- match.arg(what)
		if( what == 'consensus' ){
			# hierachical clustering on the consensus matrix
			consensushc(consensus(object), ...)
						
		}else if( what == 'fit' )
			consensushc(fit(object), ...)
		
	}
)

#' Returns the cluster membership index from an NMF model fitted with multiple 
#' runs.
#' 
#' Besides the type of clustering available for any NMF models 
#' (\code{'columns', 'rows', 'samples', 'features'}), this method can return 
#' the cluster membership index based on the consensus matrix, computed from 
#' the multiple NMF runs.
#' 
#' Argument \code{what} accepts the following extra types:
#' \describe{
#' \item{\code{'chc'}}{ returns the cluster membership based on the 
#' hierarchical clustering of the consensus matrix, as performed by 
#' \code{\link{consensushc}}.}
#' \item{\code{'consensus'}}{ same as \code{'chc'} but the levels of the membership 
#' index are re-labeled to match the order of the clusters as they would be displayed on the 
#' associated dendrogram, as re-ordered on the default annotation track in consensus 
#' heatmap produced by \code{\link{consensusmap}}.}
#' }
#' 
setMethod('predict', signature(object='NMFfitX'),
	function(object, what=c('columns', 'rows', 'samples', 'features', 'consensus', 'chc'), dmatrix = FALSE, ...){
		# determine which prediction to do
		what <- match.arg(what)
		res <- if( what %in% c('consensus', 'chc') ){
			# build the tree from consensus matrix
			h <- consensushc(object, what='consensus', dendrogram=FALSE)
			# extract membership from the tree
			cl <- cutree(h, k=nbasis(object))
			
			# rename the cluster ids in the case of a consensus map
			if( what != 'chc' ){
                dr <- as.dendrogram(h)
                o <- order.dendrogram(reorder(dr, rowMeans(consensus(object), na.rm=TRUE)))
				cl <- setNames(match(cl, unique(cl[o])), names(cl))
            }
            
			res <- as.factor(cl)
            # add dissimilarity matrix if requested
            if( dmatrix ){
                attr(res, 'dmatrix') <- 1 - consensus(object) 
            }
            if( what != 'chc' ) attr(res, 'iOrd') <- o
            
            # return
            res
		}
		else predict(fit(object), what=what, ..., dmatrix = dmatrix)
        attr(res, 'what') <- what
        res
	}
)

#' Returns the model object that achieves the lowest residual approximation 
#' error across all the runs.
#' 
#' It is a pure virtual method defined to ensure \code{fit} is defined 
#' for sub-classes of \code{NMFfitX}, which throws an error if called.
setMethod('fit', 'NMFfitX', 
	function(object){	
		stop("NMF::NMFfitX - missing definition for pure virtual method 'fit' in class '", class(object), "'")
	}
)

#' Returns the fit object that achieves the lowest residual approximation 
#' error across all the runs.
#' 
#' It is a pure virtual method defined to ensure \code{minfit} is defined 
#' for sub-classes of \code{NMFfitX}, which throws an error if called.
setMethod('minfit', 'NMFfitX',
	function(object){
		stop("NMF::NMFfitX - missing definition for pure virtual method 'minfit' in class '", class(object), "'")
	}
)

#' Show method for objects of class \code{NMFfitX}
#' @export
setMethod('show', 'NMFfitX', 
		function(object){
			cat("<Object of class:", class(object), ">\n")
			# name of the algorithm
			cat("  Method:", algorithm(object), "\n")
			# number of runs
			cat("  Runs: ", nrun(object),"\n");
			# initial state
			cat("  RNG:\n  ", RNGstr(getRNG1(object)),"\n");
			if( nrun(object) > 0 ){
				# show total timing			
				cat("  Total timing:\n"); show(runtime.all(object));
			}
		}
)

#' Extracting RNG Data from NMF Objects
#' 
#' The \code{\link{nmf}} function returns objects that contain embedded RNG data, 
#' that can be used to exactly reproduce any computation.
#' These data can be extracted using dedicated methods for the S4 generics 
#' \code{\link[rngtools]{getRNG}} and \code{\link[rngtools]{getRNG1}}.
#' 
#' @inheritParams rngtools::getRNG
#' @inheritParams rngtools::getRNG1
#' 
#' @inline
#' @rdname RNG
#' @export
setGeneric('getRNG1', package='rngtools')

#' Returns the RNG settings used for the first NMF run of multiple NMF runs. 
#' 
#' @examples
#' # For multiple NMF runs, the RNG settings used for the first run is also stored
#' V <- rmatrix(20,10)
#' res <- nmf(V, 3, nrun=4)
#' # RNG used for the best fit
#' getRNG(res)
#' # RNG used for the first of all fits
#' getRNG1(res)
#' # they may differ if the best fit is not the first one
#' rng.equal(res, getRNG1(res))
#' 
setMethod('getRNG1', signature(object='NMFfitX'),
	function(object){
		stop("NMF::getRNG1(", class(object), ") - Unimplemented pure virtual method: could not extract initial RNG settings.")
	}
)

#' Compares two NMF models when at least one comes from multiple NMF runs. 
setMethod('nmf.equal', signature(x='NMFfitX', y='NMF'), 
	function(x, y, ...){
		nmf.equal(fit(x), y, ...)
	}
)
#' Compares two NMF models when at least one comes from multiple NMF runs.
setMethod('nmf.equal', signature(x='NMF', y='NMFfitX'), 
		function(x, y, ...){
			nmf.equal(x, fit(y), ...)
		}
)

#' Returns the residuals achieved by the best fit object, i.e. the lowest 
#' residual approximation error achieved across all NMF runs.
setMethod('residuals', signature(object='NMFfitX'), 
		function(object, ...){
			residuals(minfit(object), ...)
		}
)
#' Returns the deviance achieved by the best fit object, i.e. the lowest 
#' deviance achieved across all NMF runs.
setMethod('deviance', signature(object='NMFfitX'), 
		function(object, ...){
			deviance(minfit(object), ...)
		}
)

#########################################################
# END_NMFfitX
#########################################################


#' Structure for Storing the Best Fit Amongst Multiple NMF Runs
#' 
#' This class is used to return the result from a multiple run of a single NMF
#' algorithm performed with function \code{nmf} with the -- default -- option
#' \code{keep.all=FALSE} (cf. \code{\link{nmf}}).
#' 
#' It extends both classes \code{\linkS4class{NMFfitX}} and
#' \code{\linkS4class{NMFfit}}, and stores a the result of the best fit in its
#' \code{NMFfit} structure.
#' 
#' Beside the best fit, this class allows to hold data about the computation of
#' the multiple runs, such as the number of runs, the CPU time used to perform
#' all the runs, as well as the consensus matrix.
#' 
#' Due to the inheritance from class \code{NMFfit}, objects of class
#' \code{NMFfitX1} can be handled exactly as the results of single NMF run --
#' as if only the best run had been performed.
#' 
#' 
#' @slot consensus object of class \code{matrix} used to store the
#' consensus matrix based on all the runs.
#' 
#' @slot nrun an \code{integer} that contains the number of runs
#' performed to compute the object.
#' 
#' @slot rng1 an object that contains RNG settings used for the first
#' run. See \code{\link{getRNG1}}.
#'
#' @export
#' @family multipleNMF 
#' @examples
#' 
#' # generate a synthetic dataset with known classes
#' n <- 20; counts <- c(5, 2, 3);
#' V <- syntheticNMF(n, counts)
#' 
#' # get the class factor
#' groups <- V$pData$Group
#' 
#' # perform multiple runs of one algorithm, keeping only the best fit (default)
#' #i.e.: the implicit nmf options are .options=list(keep.all=FALSE) or .options='-k'
#' res <- nmf(V, 3, nrun=3) 
#' res
#' 
#' # compute summary measures
#' summary(res)
#' # get more info
#' summary(res, target=V, class=groups)
#' 
#' # show computational time
#' runtime.all(res)
#' 
#' # plot the consensus matrix, as stored (pre-computed) in the object
#' \dontrun{ consensusmap(res, annCol=groups) }
#' 
setClass('NMFfitX1'
	, representation(
			#fit = 'NMFfit' # holds the best fit from all the runs
			consensus = 'matrix' # average connectivity matrix of all the NMF runs
			, nrun = 'integer'
			, rng1 = 'ANY'
	)
	, contains=c('NMFfitX', 'NMFfit')
	, prototype=prototype(
			consensus =	matrix(as.numeric(NA),0,0)
			, nrun = as.integer(0)
	)
)



#' Show method for objects of class \code{NMFfitX1}
#' @export
setMethod('show', 'NMFfitX1', 
	function(object){
		callNextMethod(object)
		
		# show details of the best fit
		#cat(" # Best fit:\n  ")
		#s <- capture.output(show(fit(object)))
		#cat(s, sep="\n  |")
	}
)

#' Returns the number of NMF runs performed, amongst which \code{object} was 
#' selected as the best fit.
setMethod('nrun', 'NMFfitX1', 
	function(object){
		slot(object,'nrun')
	}
)

#' Returns the consensus matrix computed while performing all NMF runs, 
#' amongst which \code{object} was selected as the best fit.
#' 
#' The result is the matrix stored in slot \sQuote{consensus}.
#' This method returns \code{NULL} if the consensus matrix is empty.
setMethod('consensus', signature(object='NMFfitX1'), 
	function(object, no.attrib = FALSE){
		
		C <- slot(object, 'consensus')
		if( length(C) > 0 ){
			if( !no.attrib ){
				class(C) <- c(class(C), 'NMF.consensus')
				attr(C, 'nrun') <- nrun(object)
				attr(C, 'nbasis') <- nbasis(object)
			}
			C
		}else NULL
		
	}
)

#' Returns the fit object associated with the best fit, amongst all the 
#' runs performed when fitting \code{object}.
#' 
#' Since \code{NMFfitX1} objects only hold the best fit, this method simply 
#' returns \code{object} coerced into an \code{NMFfit} object.
setMethod('minfit', 'NMFfitX1',
	function(object){	
		# coerce the object into a NMFfit object
		as(object, 'NMFfit')
	}
)

#' Returns the model object associated with the best fit, amongst all the 
#' runs performed when fitting \code{object}.
#' 
#' Since \code{NMFfitX1} objects only hold the best fit, this method simply 
#' returns the NMF model fitted by \code{object} -- that is stored in slot 
#' \sQuote{fit}.
setMethod('fit', signature(object='NMFfitX1'),
	function(object){
		slot(object, 'fit')
	}
)

#' Returns the RNG settings used to compute the first of all NMF runs, amongst
#' which \code{object} was selected as the best fit.
setMethod('getRNG1', signature(object='NMFfitX1'),
	function(object){
		object@rng1
	}
)

#' Compares the NMF models fitted by multiple runs, that only kept the best fits.
setMethod('nmf.equal', signature(x='NMFfitX1', y='NMFfitX1'), 
	function(x, y, ...){
		nmf.equal(fit(x), fit(y), ...)
	}
)

#########################################################
# END_NMFfitX1
#########################################################


#' Structure for Storing All Fits from Multiple NMF Runs
#' 
#' This class is used to return the result from a multiple run of a single NMF
#' algorithm performed with function \code{nmf} with option
#' \code{keep.all=TRUE} (cf. \code{\link{nmf}}).
#' 
#' It extends both classes \code{\linkS4class{NMFfitX}} and \code{list}, and
#' stores the result of each run (i.e. a \code{NMFfit} object) in its
#' \code{list} structure.
#' 
#' IMPORTANT NOTE: This class is designed to be \strong{read-only}, even though
#' all the \code{list}-methods can be used on its instances. Adding or removing
#' elements would most probably lead to incorrect results in subsequent calls.
#' Capability for concatenating and merging NMF results is for the moment only
#' used internally, and should be included and supported in the next release of
#' the package.
#' 
#' 
#' @slot .Data standard slot that contains the S3 \code{list} object data.
#' See R documentation on S3/S4 classes for more details (e.g., \code{\link{setOldClass}}).
#' 
#' @export
#' @family multipleNMF
#' @examples
#' 
#' # generate a synthetic dataset with known classes
#' n <- 20; counts <- c(5, 2, 3);
#' V <- syntheticNMF(n, counts)
#' 
#' # get the class factor
#' groups <- V$pData$Group
#' 
#' # perform multiple runs of one algorithm, keeping all the fits
#' res <- nmf(V, 3, nrun=3, .options='k') # .options=list(keep.all=TRUE) also works
#' res
#'  
#' summary(res)
#' # get more info
#' summary(res, target=V, class=groups)
#' 
#' # compute/show computational times
#' runtime.all(res)
#' seqtime(res)
#' 
#' # plot the consensus matrix, computed on the fly
#' \dontrun{ consensusmap(res, annCol=groups) }
#' 
setClass('NMFfitXn'
	, contains=c('NMFfitX', 'list')
	, validity=function(object){
				
		# the list must only contains NMFfit objects of the same dimensions
		ref.dim <- NULL
		ref.algo <- NULL
		for(i in seq_along(object)){
			# check class of the element
			item <- object[[i]]
			if( !(is(item, 'NMFfit') && !is(item, 'NMFfitX')) )
				return(paste("invalid class for element", i, "of input list [all elements must be a NMFfit object]"))
						
			# check dimensions
			if( is.null(ref.dim) ) ref.dim <- dim(item)	
			if( !identical(ref.dim, dim(item)) )
				return(paste("invalid dimension for element", i, "of input list [all elements must have the same dimensions]"))
			
			# check algorithm names
			if( is.null(ref.algo) ) ref.algo <- algorithm(item)	
			if( !identical(ref.algo, algorithm(item)) )
				return(paste("invalid algorithm for element", i, "of input list [all elements must result from the same algorithm]"))
			
		}
		
	}
)

# Updater for slot .Data
#objectUpdater('NMFfitXn', '0.5.06'
#	, vfun=function(object){ !.hasSlot(object, 'rng1') }
#	, function(x, y){
#		y@.Data <- lapply(x@.Data, nmfObject)
#	}
#)
 
#' Show method for objects of class \code{NMFfitXn}
#' @export
setMethod('show', 'NMFfitXn', 
	function(object){
		callNextMethod(object)
		
		# if the object is not empty and slot runtime.all is not null then show
		# the sequential time, as it might be different from runtime.all
		if( length(object) > 0 && !is.null(runtime.all(object, null=TRUE)) ){
			# show total sequential timing
			cat("  Sequential timing:\n"); show(seqtime(object));
		}
	}
)

#' Returns the number of basis components common to all fits.
#' 
#' Since all fits have been computed using the same rank, it returns the 
#' factorization rank of the first fit.
#' This method returns \code{NULL} if the object is empty.  
setMethod('nbasis', signature(x='NMFfitXn'), 
	function(x, ...){
		if( length(x) == 0 ) return(NULL)
		return( nbasis(x[[1]]) )
	}
)
#' Returns the dimension common to all fits.
#' 
#' Since all fits have the same dimensions, it returns the dimension of the 
#' first fit.
#' This method returns \code{NULL} if the object is empty.
#' 
#' @rdname dims
setMethod('dim', signature(x='NMFfitXn'), 
	function(x){
		if( length(x) == 0 ) return(NULL)
		return( dim(x[[1L]]) )
	}
)

#' Returns the coefficient matrix of the best fit amongst all the fits stored in 
#' \code{object}.
#' It is a shortcut for \code{coef(fit(object))}.  
setMethod('coef', signature(object='NMFfitXn'), 
	function(object, ...){
		coef(fit(object), ...)
	}
)

#' Returns the basis matrix of the best fit amongst all the fits stored in 
#' \code{object}.
#' It is a shortcut for \code{basis(fit(object))}.
setMethod('basis', signature(object='NMFfitXn'), 
	function(object, ...){
		basis(fit(object), ...)
	}
)

#' Method for multiple NMF fit objects, which returns the indexes of fixed basis 
#' terms from the best fitted model.
setMethod('ibterms', 'NMFfitX', 
	function(object){
		ibterms(fit(object))
	}
)
#' Method for multiple NMF fit objects, which returns the indexes of fixed 
#' coefficient terms from the best fitted model.
setMethod('icterms', 'NMFfit', 
	function(object){
		icterms(fit(object))
	}
)


#' Returns the number of runs performed to compute the fits stored in the list 
#' (i.e. the length of the list itself).
setMethod('nrun', 'NMFfitXn', 
	function(object){
		length(object)
	}
)

#' Returns the name of the common NMF algorithm used to compute all fits 
#' stored in \code{object}
#' 
#' Since all fits are computed with the same algorithm, this method returns the 
#' name of algorithm that computed the first fit.
#' It returns \code{NULL} if the object is empty.
setMethod('algorithm', 'NMFfitXn', 
	function(object){
		if( length(object) == 0 ) return(NULL)
		return( algorithm(object[[1]]) )
	} 
)
#' Returns the name of the common seeding method used the computation of all fits 
#' stored in \code{object}
#' 
#' Since all fits are seeded using the same method, this method returns the 
#' name of the seeding method used for the first fit.
#' It returns \code{NULL} if the object is empty.
setMethod('seeding', 'NMFfitXn',
	function(object){
		if( length(object) == 0 ) return(NULL)
		return( seeding(object[[1]]) )
	}
)
#' Returns the common type NMF model of all fits stored in \code{object}
#' 
#' Since all fits are from the same NMF model, this method returns the 
#' model type of the first fit.
#' It returns \code{NULL} if the object is empty.
setMethod('modelname', signature(object='NMFfitXn'), 
	function(object){
		if( length(object) == 0 ) return(NULL)
		return( modelname(object[[1]]) )
	}
)

#' Returns the CPU time that would be required to sequentially compute all NMF 
#' fits stored in \code{object}.
#' 
#' This method calls the function \code{runtime} on each fit and sum up the 
#' results.
#' It returns \code{NULL} on an empty object.
setMethod('seqtime', 'NMFfitXn', 
		function(object){	
			if( length(object) == 0 ) return(NULL)
			# sum up the time across the runs
			.seqtime(object)
		}
)

#' Returns the CPU time used to perform all the NMF fits stored in \code{object}.
#' 
#' If no time data is available from in slot \sQuote{runtime.all} and argument 
#' \code{null=TRUE}, then the sequential time as computed by 
#' \code{\link{seqtime}} is returned, and a warning is thrown unless \code{warning=FALSE}.
#' 
#' @param null a logical that indicates if the sequential time should be returned
#' if no time data is available in slot \sQuote{runtime.all}. 
#' @param warning a logical that indicates if a warning should be thrown if the 
#' sequential time is returned instead of the real CPU time.    
#' 
setMethod('runtime.all', 'NMFfitXn', 
	function(object, null=FALSE, warning=TRUE){
		
		if( length(object) == 0 ) return(NULL)
		stored.time <- slot(object, 'runtime.all')
		# if there is some time stored, return it
		if( length(stored.time) > 0 ) stored.time
		else if( null ) NULL
		else{
			if( warning )
				warning("NMFfitXn::runtime.all - computation time data not available [sequential time was used instead]")
			seqtime(object) # otherwise total sequential time
		}
		
	}
)

#' Returns the best NMF model in the list, i.e. the run that achieved the lower 
#' estimation residuals.
#' 
#' The model is selected based on its \code{deviance} value.
#'  
setMethod('minfit', 'NMFfitXn',
	function(object){
		
		b <- which.best(object, deviance)
		# test for length 0
		if( length(b) == 0 ) return(NULL)
				
		# return the run with the lower
		object[[ b ]]
	}
)

#' \code{which.best} returns the index of the best fit in a list of NMF fit, 
#' according to some quantitative measure.
#' The index of the fit with the lowest measure is returned.
#' 
#' @param object an NMF model fitted by multiple runs. 
#' @param FUN the function that computes the quantitative measure.
#' @param ... extra arguments passed to \code{FUN}. 
#' 
#' @export
#' @rdname advanced
which.best <- function(object, FUN=deviance, ...){
	
	# test for length 0
	if( length(object) == 0 ) 
		return(integer())
	
	# retrieve the measure for each run
	e <- sapply(object, FUN, ...)
	
	# return the run with the lower
	which.min(e)
}

#' Returns the RNG settings used for the first run.
#' 
#' This method throws an error if the object is empty.
setMethod('getRNG1', signature(object='NMFfitXn'),
	function(object){
		if( length(object) == 0 )
			stop("NMF::getRNG1 - Could not extract RNG data from empty object [class:", class(object), "]")
		
		getRNG(object[[1]])
	}
)

#' @inline
#' @rdname RNG
#' @export
setGeneric('.getRNG', package='rngtools')

#' Returns the RNG settings used for the best fit.
#' 
#' This method throws an error if the object is empty.
setMethod('.getRNG', signature(object='NMFfitXn'),
	function(object, ...){
		if( length(object) == 0 )
			stop("NMF::getRNG - Could not extract RNG data from empty object [class:", class(object), "]")
		
		getRNG(minfit(object), ...)
	}
)


#' Returns the best NMF fit object amongst all the fits stored in \code{object},
#' i.e. the fit that achieves the lowest estimation residuals.
setMethod('fit', signature(object='NMFfitXn'),
	function(object){
		fit( minfit(object) )
	}
)

#' Compares the results of multiple NMF runs.
#' 
#' This method either compare the two best fit, or all fits separately.
#' All extra arguments in \code{...} are passed to each internal call to 
#' \code{nmf.equal}. 
#' 
#' @param all a logical that indicates if all fits should be compared separately
#' or only the best fits
#' @param vector a logical, only used when \code{all=TRUE}, that indicates if 
#' all fits must be equal for \code{x} and \code{y} to be declared equal, or 
#' if one wants to return the result of each comparison in a vector.   
#' 
#' @inline
setMethod('nmf.equal', signature(x='list', y='list'), 
	function(x, y, ..., all=FALSE, vector=FALSE){
		if( !all )
			nmf.equal(x[[ which.best(x) ]], y[[ which.best(y) ]], ...)
		else{
			if( length(x) != length(y) )
				FALSE
			else
				res <- mapply(function(a,b,...) isTRUE(nmf.equal(a,b,...)), x, y, MoreArgs=list(...))
				if( !vector )
					res <- all( res )
				res
		}
	}
)

#' Compare all elements in \code{x} to \code{x[[1]]}.
setMethod('nmf.equal', signature(x='list', y='missing'), 
	function(x, y, ...){
		
		if( length(x) == 0L ){
			warning("Empty list argument `x`: returning NA")
			return(NA)
		}
		if( length(x) == 1L ){
			warning("Only one element in list argument `x`: returning TRUE")
			return(TRUE)
		}
		for( a in x ){
			if( !nmf.equal(x[[1]], a, ...) ) return(FALSE)
		}
		return(TRUE)
	}
)

#' Computes the consensus matrix of the set of fits stored in \code{object}, as
#' the mean connectivity matrix across runs.
#' 
#' This method returns \code{NULL} on an empty object.
#' The result is a matrix with several attributes attached, that are used by 
#' plotting functions such as \code{\link{consensusmap}} to annotate the plots.
#' 
#' @aliases plot.NMF.consensus
setMethod('consensus', signature(object='NMFfitXn'), 
	function(object, ..., no.attrib = FALSE){
		if( length(object) == 0 ) return(NULL)
		
		# init empty consensus matrix
		con <- matrix(0, ncol(object), ncol(object))
		# name the rows and columns appropriately: use the sample names of the first fit
		dimnames(con) <- list(colnames(object[[1]]), colnames(object[[1]]))
		
		# compute mean connectivity matrix
		sapply(object 
				, function(x, ...){
					con <<- con + connectivity(x, ..., no.attrib = TRUE)
					NULL
				}
				, ...
		)
		con <- con / nrun(object)
				
		# return result
		if( !no.attrib ){
			class(con) <- c(class(con), 'NMF.consensus')
			attr(con, 'nrun') <- nrun(object)
			attr(con, 'nbasis') <- nbasis(object)
		}
		con
	}
)

#' @S3method plot NMF.consensus 
plot.NMF.consensus <- function(x, ...){
	consensusmap(x, ...)
}

#' Dispersion of a Matrix
#' 
#' Computes the dispersion coefficient of a -- consensus -- matrix
#' \code{object}, generally obtained from multiple NMF runs.
#' 
#' The dispersion coefficient is based on the consensus matrix (i.e. the
#' average of connectivity matrices) and was proposed by \cite{KimH2007} to 
#' measure the reproducibility of the clusters obtained from NMF.
#' 
#' It is defined as: 
#' \deqn{\rho = \sum_{i,j=1}^n 4 (C_{ij} - \frac{1}{2})^2 , }
#' where \eqn{n} is the total number of samples.
#' 
#' By construction, \eqn{0 \leq \rho \leq 1} and \eqn{\rho = 1} only for a perfect
#' consensus matrix, where all entries 0 or 1. 
#' A perfect consensus matrix is obtained only when all the connectivity matrices 
#' are the same, meaning that the algorithm gave the same clusters at each run.  
#' See \cite{KimH2007}.
#' 
#' @param object an object from which the dispersion is computed
#' @param ... extra arguments to allow extension  
#'
#' @export 
setGeneric('dispersion', function(object, ...) standardGeneric('dispersion') )
#' Workhorse method that computes the dispersion on a given matrix.
setMethod('dispersion', 'matrix', 
	function(object, ...){
		stopifnot( nrow(object) == ncol(object) )
		sum( 4 * (object-1/2)^2 ) / nrow(object)^2
	}
)
#' Computes the dispersion on the consensus matrix obtained from multiple NMF
#' runs. 
setMethod('dispersion', 'NMFfitX', 
	function(object, ...){
		dispersion(consensus(object), ...)
	}
)

#' Factory Method for Multiple NMF Run Objects
#' 
#' @param object an object from which is created an \code{NMFfitX} object
#' @param ... extra arguments used to pass values for slots
#' 
#' @inline
#' @keywords internal
setGeneric('NMFfitX', function(object, ...) standardGeneric('NMFfitX') )
#' Create an \code{NMFfitX} object from a list of fits.
#' 
#' @param .merge a logical that indicates if the fits should be aggregated, only 
#' keeping the best fit, and return an \code{NMFfitX1} object.
#' If \code{FALSE}, an \code{NMFfitXn} object containing the data of all the fits
#' is returned.
#' 
setMethod('NMFfitX', 'list',
	function(object, ..., .merge=FALSE){
		
		if( length(object) == 0 )
			return(new('NMFfitXn'))
		else if( is(object, 'NMFfitXn') && !.merge)
			return(object)
		
		# retrieve the extra arguments
		extra <- list(...)
				
		# if runtime.all is provided: be sure it's of the right class
		tt <- extra$runtime.all
		compute.tt <- TRUE
		if( !is.null(tt) ){
			if( !is(tt, 'proc_time') ){
				if( !is.numeric(tt) || length(tt) != 5 )
					stop("NMF::NMFfitX - invalid value for 'runtime.all' [5-length numeric expected]")
				class(extra$runtime.all) <- 'proc_time'
			}
			compute.tt <- FALSE
		}else{
			extra$runtime.all <- rep(0,5)
			class(extra$runtime.all) <- 'proc_time'
		}
		
		# check validity and aggregate if required
		ref.algo <- NULL
		ref.class <- NULL
		nrun <- 0
		lapply( seq_along(object)
			, function(i){
				item <- object[[i]]
				
				# check the type of each element				
				if( !(is(item, 'NMFfitX') || is(item, 'NMFfit')) )
					stop("NMF::NMFfitX - invalid class for element ", i, " of input list [all elements must be NMFfit or NMFfitX objects]")
				
				# check that all elements result from the same algorithm
				if( is.null(ref.algo) ) ref.algo <<- algorithm(item)
				if( !identical(algorithm(item), ref.algo) )
					stop("NMF::NMFfitX - invalid algorithm for element ", i, " of input list [cannot join results from different algorithms]")
				
				# check if simple join is possible: only Ok if all elements are from the same class (NMFfit or NMFfitXn)
				if( length(ref.class) <= 1 ) ref.class <<- unique(c(ref.class, class(item)))
				
				# sum up the number of runs
				nrun <<- nrun + nrun(item)
				
				# compute total running time if necessary
				if( compute.tt )
					extra$runtime.all <<- extra$runtime.all + runtime.all(item)
				
			}
		)
		
		# force merging if the input list is hetergeneous or if it only contains NMFfitX1 objects
		if( length(ref.class) > 1 || ref.class == 'NMFfitX1' ){
			nmf.debug('NMFfitX', ".merge is forced to TRUE")
			.merge <- TRUE
		}
		
		# unpack all the NMFfit objects
		object.list <- unlist(object)
		nmf.debug('NMFfitX', "Number of fits to join = ", length(object.list))
					
		# one wants to keep only the best result
		if( .merge ){
			
			warning("NMF::NMFfitX - The method for merging lists is still in development")
			
			# set the total number of runs
			extra$nrun <- as.integer(nrun)			
									
			# consensus matrix
			if( !is.null(extra$consensus) )
				warning("NMF::NMFfitX - the value of 'consensus' was discarded as slot 'consensus' is computed internally")
			extra$consensus <- NULL
									
			consensus <- matrix(as.numeric(NA), 0, 0)
			best.res <- Inf		
			best.fit <- NULL
			sapply(object.list, function(x){
				if( !is(x, 'NMFfit') )
					stop("NMF::NMFfitX - all inner-elements of '",substitute(object),"' must inherit from class 'NMFfit'")
				
				# merge consensus matrices
				consensus <<- if( sum(dim(consensus)) == 0 ) nrun(x) * consensus(x)
							  else consensus + nrun(x) * consensus(x)
					  
				temp.res <- residuals(x)
				if( temp.res < best.res ){
					# keep best result
					best.fit <<- minfit(x)					
					best.res <<- temp.res
				}
			})
			# finalize consensus matrix
			consensus <- consensus/extra$nrun
			extra$consensus <- consensus 
									
			# return merged result
			return( do.call(NMFfitX, c(list(best.fit), extra)) )
		}
		else{
			# create a NMFfitXn object that holds the whole list			
			do.call('new', c(list('NMFfitXn', object.list), extra))
		}
	}
)
#' Creates an \code{NMFfitX1} object from a single fit.
#' This is used in \code{\link{nmf}} when only the best fit is kept in memory or 
#' on disk.
#'  
setMethod('NMFfitX', 'NMFfit',
		function(object, ...){
					
			extra <- list(...)
						
			# default value for nrun is 1 
			if( is.null(extra$nrun) ) extra$nrun = as.integer(1)
			
			# a consensus matrix is required (unless nrun is 1)
			if( is.null(extra$consensus) ){
				if( extra$nrun == 1 ) 
					extra$consensus <- connectivity(object)
				else
					stop("Slot 'consensus' is required to create a 'NMFfitX1' object where nrun > 1")				
			}
			
			# slot runtime.all is inferred if missing and nrun is 1
			if( is.null(extra$runtime.all) && extra$nrun == 1 )
				extra$runtime.all <- runtime(object)
			
			# create the NMFfitX1 object
			do.call('new', c(list('NMFfitX1', object), extra))
		}
)
#' Provides a way to aggregate \code{NMFfitXn} objects into an \code{NMFfitX1} 
#' object.
setMethod('NMFfitX', 'NMFfitX',
		function(object, ...){

			# nothing to do in the case of NMFfitX1 objects
			if( is(object, 'NMFfitX1') ) return(object)
			
			# retrieve extra arguments
			extra <- list(...)
			
			# take runtime.all from the object itself
			if( !is.null(extra$runtime.all) )
				warning("NMF::NMFfitX - argument 'runtime.all' was discarded as it is computed from argument 'object'")			
			extra$runtime.all <- runtime.all(object)						
			
			# create the NMFfitX1 object
			f <- selectMethod(NMFfitX, 'list')
			do.call(f, c(list(object), extra))
		}
)

#' Computes the best or mean purity across all NMF fits stored in \code{x}.
#' 
#' @param method a character string that specifies how the value is computed.
#' It may be either \code{'best'} or \code{'mean'} to compute the best or mean 
#' purity respectively.
#'  
#' @inline
setMethod('purity', signature(x='NMFfitXn', y='ANY'), 
	function(x, y, method='best', ...){
		c <- sapply(x, purity, y=y, ...)
		
		# aggregate the results if a method is provided
		if( is.null(method) ) c
		else aggregate.measure(c, method, decreasing=TRUE)		
	}
)

#' Computes the best or mean entropy across all NMF fits stored in \code{x}.
#' 
#' @inline
setMethod('entropy', signature(x='NMFfitXn', y='ANY'), 
	function(x, y, method='best', ...){
		c <- sapply(x, entropy, y=y, ...)		
		
		# aggregate the results if a method is provided
		if( is.null(method) ) c
		else aggregate.measure(c, method)
	}
)

###% Utility function to aggregate numerical quality measures from \code{NMFfitXn} objects.
###% 
###% Given a numerical vector, this function computes an aggregated value using one of the following methods:
###% - mean: the mean of the measures
###% - best: the best measure according to the specified sorting order (decreasing or not)
###%  
aggregate.measure <- function(measure, method=c('best', 'mean'), decreasing=FALSE){
	# aggregate the results
	method <- match.arg(method)
	res <- switch(method
			, mean = mean(measure)
			, best = if( decreasing ) max(measure) else min(measure)
	)
	
	# set the name to 
	names(res) <- method
	
	# return result
	res
} 

#' Computes a set of measures to help evaluate the quality of the \emph{best
#' fit} of the set. 
#' The result is similar to the result from the \code{summary} method of 
#' \code{NMFfit} objects. 
#' See \code{\linkS4class{NMF}} for details on the computed measures.  
#' In addition, the cophenetic correlation (\code{\link{cophcor}}) and 
#' \code{\link{dispersion}} coefficients of the consensus matrix are returned, 
#' as well as the total CPU time (\code{\link{runtime.all}}).
#'   
setMethod('summary', signature(object='NMFfitX'),
	function(object, ...){
		
		# compute summary measures for the best fit
		best.fit <- minfit(object)
		s <- summary(best.fit, ...)
		# get totaltime
		t <- runtime.all(object)		
		
		# replace cpu.all and nrun in the result (as these are set by the summary method of class NMFfit)
		s[c('cpu.all', 'nrun')] <- c(as.numeric(t['user.self']+t['user.child']), nrun(object))
		
		# compute cophenetic correlation coeff and dispersion
		C <- consensus(object)
		s <- c(s, cophenetic=cophcor(C), dispersion=dispersion(C))
		
        # compute mean consensus silhouette width
		si <- silhouette(object, what = 'consensus')
		s <- c(s, silhouette.consensus = if( !is_NA(si) ) summary(si)$avg.width else NA)
        
		# return result
		s
	}
)


#' Comparing Results from Different NMF Runs
#' 
#' The functions documented here allow to compare the fits computed in 
#' different NMF runs.
#' The fits do not need to be from the same algorithm, nor have the same
#' dimension.
#' 
#' The methods \code{compare} enables to compare multiple NMF fits either 
#' passed as arguments or as a list of fits.
#' These methods eventually call the method \code{summary,NMFList}, so that
#' all its arguments can be passed \strong{named} in \code{...}. 
#' 
#' @param ... extra arguments passed by \code{compare} to \code{summary,NMFList}
#' or to the \code{summary} method of each fit.
#' 
#' @name compare-NMF
#' @rdname nmf-compare
NULL

.compare_NMF <- function(...){
	args <- list(...)
	
	iargs <-
	if( is.null(names(args)) ){
		names(args) <- rep("", length(args)) 
		seq(args)
	}else{
		iargs <- which(names(args)=='')
		if( length(iargs) != length(args) )
			iargs <- iargs[ iargs < which(names(args)!='')[1L] ]
		iargs
	}

	lfit <- args[iargs]
	lfit <- unlist(lfit, recursive=FALSE)
	
	# wrap up into an NMFList object
	object <- as.NMFList(lfit)
	do.call('summary', c(list(object), args[-iargs]))
}

#' Compare multiple NMF fits passed as arguments.
#' 
#' @rdname nmf-compare
#' 
#' @examples
#' 
#' x <- rmatrix(20,10)
#' res <- nmf(x, 3)
#' res2 <- nmf(x, 2, 'lee')
#' 
#' # compare arguments
#' compare(res, res2, target=x)
#' 
setMethod('compare', signature(object='NMFfit'),
	function(object, ...){
		.compare_NMF(object, ...)
	}
)
#' Compares the fits obtained by separate runs of NMF, in a single 
#' call to \code{\link{nmf}}.
#' 
#' @rdname nmf-compare
#' 
#' # compare each fits in a multiple runs
#' res3 <- nmf(x, 2, nrun=3, .opt='k')
#' compare(res3)
#' compare(res3, res, res2)
#' compare(list(res3), res, res2, target=x)
#' 
setMethod('compare', signature(object='NMFfitXn'),
	function(object, ...){
		do.call(.compare_NMF, c(unlist(object), list(...)))
	}
)
#' Compares multiple NMF fits passed as a standard list.
#' 
#' @rdname nmf-compare
#'  
#' @examples
#' # compare elements of a list
#' compare(list(res, res2), target=x)
setMethod('compare', signature(object='list'),
	function(object, ...){
		do.call(.compare_NMF, c(list(object), list(...)))
	}
)

#' @details
#' \code{summary,NMFList} computes summary measures for each NMF result in the list 
#' and return them in rows in a \code{data.frame}. 
#' By default all the measures are included in the result, and \code{NA} values 
#' are used where no data is available or the measure does not apply to the 
#' result object (e.g. the dispersion for single' NMF runs is not meaningful). 
#' This method is very useful to compare and evaluate the performance of 
#' different algorithms.
#' 
#' @param select the columns to be output in the result \code{data.frame}.  The
#' column are given by their names (partially matched).  The column names are
#' the names of the summary measures returned by the \code{summary} methods of
#' the corresponding NMF results.
#' @param sort.by the sorting criteria, i.e. a partial match of a column name,
#' by which the result \code{data.frame} is sorted.  The sorting direction
#' (increasing or decreasing) is computed internally depending on the chosen
#' criteria (e.g. decreasing for the cophenetic coefficient, increasing for the
#' residuals). 
#' 
#' @rdname nmf-compare
setMethod('summary', signature(object='NMFList'),
	function(object, sort.by=NULL, select=NULL, ...){
		
		if( length(object) == 0L ) return()
		
		# define the sorting schema for each criteria (TRUE for decreasing, FALSE for increasing)
		sorting.schema <- list(method=FALSE, seed=FALSE, rng=FALSE, metric=FALSE
							, residuals=FALSE, cpu=FALSE, purity=TRUE, nrun=FALSE, cpu.all=FALSE
							, cophenetic=TRUE, dispersion=TRUE #NMFfitX only
							, entropy=FALSE, sparseness.basis=TRUE, sparseness.coef=TRUE, rank=FALSE, rss=FALSE
							, niter=FALSE, evar=TRUE
                            , silhouette.coef = TRUE, silhouette.basis = TRUE
                            , silhouette.consensus = TRUE)
				
		# for each result compute the summary measures
		measure.matrix <- sapply(object, summary, ...)		
		
		# the results from 'summary' might not have the same length => generate NA where necessary
		if( is.list(measure.matrix) ){
			name.all <- unique(unlist(sapply(measure.matrix, names)))
			measure.matrix <- sapply(seq_along(measure.matrix),
				function(i){
					m <- measure.matrix[[i]][name.all]
					names(m) <- name.all
					m
				}
			)
		}
		
		# transpose the results so that methods are in lines, measures are in columns
		measure.matrix <- t(measure.matrix)	
		
		# set up the resulting data.frame
		methods <- sapply(object, function(x, ...){
					x <- minfit(x)
					m <- algorithm(x)
					s <- seeding(x) 
					svalue <- objective(x)
					svalue <- if( is.function(svalue) ) '<function>' else svalue
					c(method=m, seed=s, rng=RNGdigest(x), metric=svalue)
				}
		)
		methods <- t(methods)	
		res <- as.data.frame(methods, stringsAsFactors=FALSE)	
		
		# add the measures to the result		
		res <- cbind(res, measure.matrix)
		res$rng <- as.numeric(factor(res$rng))
				
		# sort according to the user's preference
		# ASSERT FOR DEV: all columns measure must have a defined sorting schema 
		#if( !all( no.schema <- is.element(colnames(res), names(sorting.schema))) ) 
		#	warning("ASSERT: missing sorting schema for criteria(e): ", paste(paste("'", colnames(res)[!no.schema], "'", sep=''), collapse=', '))

		if( !is.null(sort.by) ){			
			sorting.criteria <- intersect(colnames(res), names(sorting.schema))
			sort.by.ind <- pmatch(sort.by, sorting.criteria)
			if( is.na(sort.by.ind) )
				stop("NMF::summary[NMFList] : argument 'sort.by' must be NULL or partially match one of "
					, paste( paste("'", names(sorting.schema), "'", sep=''), collapse=', ')
					, call.=FALSE)
			sort.by <- sorting.criteria[sort.by.ind]
			res <- res[order(res[[sort.by]], decreasing=sorting.schema[[sort.by]]) , ]
			
			# add an attribute to the result to show the sorting criteria that was used
			attr(res, 'sort.by') <- sort.by
		}
		
		# limit the output to the required measures
		if( !is.null(select) || !missing(select) ){
			select.full <- match.arg(select, colnames(res), several.ok=TRUE)
			if( length(select.full) <  length(select) )
				stop("NMF::summary[NMFList] - the elements of argument 'select' must partially match one of "
					, paste(paste("'", colnames(res),"'", sep=''), collapse=', ')
					, call.=FALSE)
			res <- subset(res, select=select.full)
		}
				
		# return result
		res
	}
)

#' @details
#' \code{plot} plot on a single graph the residuals tracks for each fit in \code{x}. 
#' See function \code{\link{nmf}} for details on how to enable the tracking of residuals.
#' 
#' @param x an \code{NMFList} object that contains fits from separate NMF runs.
#' @param y missing
#' @inheritParams plot,NMFfit,missing-method
#' 
#' @rdname nmf-compare
setMethod('plot', signature(x='NMFList', y='missing'), 
	function(x, y, skip=-1L, ...){
		
		# retrieve normalized residuals tracks
		max.iter <- 0
		tracks <- lapply( x, 
				function(res){
					res <- minfit(res)
					t <- residuals(res, track=TRUE)
					# skip some residuals(s) if requested
					if( skip == -1L && !is.null(names(t)) ) t <- t[names(t)!='0'] # remove initial residual
					else if( skip > 0 ) t <- t[-(1:skip)]
					#print(t)
					# update max iteration
					max.iter <<- max(max.iter, as.numeric(names(t)))
					# return normalized track
					t/t[1]
				}
		)
		minT <- min(sapply(tracks, min))
		maxT <- max(sapply(tracks, max))
		
		#print(tracks)
		# create an empty plot
		# set default graphical parameters (those can be overriden by the user)
		params <- .set.list.defaults(list(...)
				, xlab='Iterations', ylab='Normalised objective values'
				, main='NMF Residuals')
		
		# setup the plot
		do.call('plot', 
				c(list(0, xlim=c(0,max.iter+100), ylim=c(minT, maxT)), col='#00000000'
				, params)
				)
		
		# add legend
		cols <- seq_along(tracks)
		legend('topright', legend=names(tracks), fill=cols
				, title='Algorithm')
		
		# plot each tracks		
		lapply( seq_along(tracks),
				function(i){
					t <- tracks[[i]]
					points(names(t), t, col=cols[i], type='p', cex=0.5)
					points(names(t), t, col=cols[i], type='l', lwd=1.4)
				})
		
		
		# return invisible
		return(invisible())
	}		
)

#' Deprecated method subsituted by \code{\link{consensusmap}}.
setMethod('metaHeatmap', signature(object='NMFfitX'),
		function(object, ...){
			# send deprecated warning
			.Deprecated('metaHeatmap', 'NMF', "Direct use of the S4-Method 'metaHeatmap' for 'NMFfitX' objects is deprecated, use 'consensusmap' instead.")

			# call the new function 'consmap'
			return( consensusmap(object, ...) )
			
		}
)

#' \code{consensusmap} plots heatmaps of consensus matrices.
#' 
#' @details
#' \code{consensusmap} redefines default values for the following arguments of
#' \code{\link{aheatmap}}:
#' \itemize{
#' \item the colour palette;
#' \item the column ordering which is set equal to the row ordering, since
#' a consensus matrix is symmetric;
#' \item the distance and linkage methods used to order the rows (and columns).
#' The default is to use 1 minus the consensus matrix itself as distance, and 
#' average linkage. 
#' \item the addition of two special named annotation tracks, \code{'basis:'} and 
#' \code{'consensus:'}, that show, for each column (i.e. each sample), 
#' the dominant basis component in the best fit and the hierarchical clustering
#' of the consensus matrix respectively (using 1-consensus as distance and average 
#' linkage). 
#' 
#' These tracks are specified in argument \code{tracks}, which behaves as in 
#' \code{\link{basismap}}.
#' 
#' \item a suitable title and extra information like the type of NMF model or the 
#' fitting algorithm, when \code{object} is a fitted NMF model. 
#' }  
#' 
#' @rdname heatmaps
#' 
#' @examples
#' 
#' \dontrun{
#' res <- nmf(x, 3, nrun=3)
#' consensusmap(res)
#' }
#' 
#' @inline
#' @export
setGeneric('consensusmap', function(object, ...) standardGeneric('consensusmap') )
#' Plots a heatmap of the consensus matrix obtained when fitting an NMF model with multiple runs. 
setMethod('consensusmap', 'NMFfitX', 
	function(object, annRow=NA, annCol=NA
			, tracks=c('basis:', 'consensus:', 'silhouette:')
			, main = 'Consensus matrix', info = FALSE
			, ...){
			
		# add side information if requested
		info <- if( isTRUE(info) ){
					paste("NMF model: '", modelname(object)
					, "'\nAlgorithm: '", algorithm(object)
					, "'\nbasis: ", nbasis(object)
					,"\nnrun: ", nrun(object), sep='')
				}else if( isFALSE(info) ) NULL
				else info
		
		x <- consensus(object)
		
		# process annotation tracks
		ptracks <- process_tracks(x, tracks, annRow, annCol)
		annRow <- ptracks$row 
		annCol <- ptracks$col
		# set special annotation handler
		ahandlers <- list(
			basis = function() predict(object)
			, consensus = function() predict(object, what='consensus')
			, silhouette = function(){
				si <- silhouette(object, what='consensus', order = NA)
				if( is_NA(si) ) NA
				else si[, 'sil_width']
			}
		)
		specialAnnotation(1L, ahandlers)
		specialAnnotation(2L, ahandlers)
		#
		
		consensusmap(x, ..., annRow=annRow, annCol=annCol, main = main, info = info)	
	}
)
#' Plots a heatmap of the connectivity matrix of an NMF model.
setMethod('consensusmap', 'NMF', 
	function(object, ...){
		consensusmap(connectivity(object), ...)		
	}
)
#' Main method that redefines default values for arguments of \code{\link{aheatmap}}.
setMethod('consensusmap', 'matrix', 
	function(object, color='-RdYlBu'
			, distfun = function(x) as.dist(1-x), hclustfun = 'average'
			, Rowv = TRUE, Colv = "Rowv"
			, main = if( is.null(nr) || nr > 1 ) 'Consensus matrix' else 'Connectiviy matrix'
			, info = FALSE
			, ...){
				
		nr <- nrun(object)
		nb <- nbasis(object)
		info <- if( isTRUE(info) ){
					info <- NULL
					if( !is.null(nr) ) info <- c(info, paste("nrun:", nr))
					if( !is.null(nb) ) info <- c(info, paste("nbasis:", nb))
					info <- c(info, paste("cophcor:", round(cophcor(object), 3)))
				}else if( isFALSE(info) ) NULL
				else info
			
		aheatmap(object, color = color, ...
				, distfun = distfun, hclustfun = hclustfun
				, Rowv = Rowv, Colv = Colv
				, main = main
				, info = info)
	}
)

setOldClass('NMF.rank')
#' Draw a single plot with a heatmap of the consensus matrix obtained for each value of the rank, 
#' in the range tested with \code{\link{nmfEstimateRank}}.
#' 
#' @rdname nmf-compare
setMethod('consensusmap', 'NMF.rank', 
	function(object, ...){

		# plot the list of consensus matrix (set names to be used as default main titles)
		consensusmap(setNames(object$fit, paste("rank = ", lapply(object$fit, nbasis))), ...)
	}
)
#' Draw a single plot with a heatmap of the consensus matrix of each element in the list \code{object}.
#' 
#' @param layout specification of the layout.
#' It may be a single numeric or a numeric couple, to indicate a square or rectangular layout 
#' respectively, that is filled row by row.
#' It may also be a matrix that is directly passed to the function \code{\link[graphics]{layout}}
#' from the package \code{graphics}.
#' 
#' @rdname nmf-compare
setMethod('consensusmap', 'list', 
	function(object, layout
			, Rowv = FALSE, main = names(object)
			, ...){
				
		opar <- par(no.readonly=TRUE)
		on.exit(par(opar))
		
		# define default layout
		if (missing(layout) ){
			n <- length(object)
			nr <- nc <- floor(sqrt(n))
			if( nr^2 != n ){
				nc <- nr + 1
				if( nr == 1 && nr*nc < n )
					nr <- nr + 1
			}
			
			layout <- c(nr, nc)		
		}
		if( !is.matrix(layout) ){
			if( !is.numeric(layout) )
				stop("invalid layout specification: must be a matrix or a numeric")
			if( length(layout) == 1 )
				layout <- c(layout, layout)
			layout <- matrix(1:(layout[1]*layout[2]), layout[1], byrow=TRUE)
		}
		
		graphics::layout(layout)
		res <- sapply(seq_along(object), function(i, ...){
			x <- object[[i]]
			
			# set main title
			main <- if( !is.null(main) && length(main) > 1 ){
				if( length(main) != length(object) )
					stop("consensusmap - Invalid length for argument `main`: should be either a single character string, or a list or vector of same length as ", deparse(substitute(object)))
				main[[i]]
			}			
			
			# call method for the fit
			consensusmap(x, ..., Rowv=Rowv, main=main)
		}, ...)
		invisible(res)
	}
)

#' Plots a heatmap of the basis matrix of the best fit in \code{object}.
setMethod('basismap', signature(object='NMFfitX'),
	function(object, ...){
		# call the method on the best fit
		basismap(minfit(object), ...)	
	}
)

#' Plots a heatmap of the coefficient matrix of the best fit in \code{object}.
#' 
#' This method adds:
#' \itemize{
#' \item an extra special column annotation track for multi-run NMF fits,
#' \code{'consensus:'}, that shows the consensus cluster associated to each sample.
#' \item a column sorting schema \code{'consensus'} that can be passed
#' to argument \code{Colv} and orders the columns using the hierarchical clustering of the 
#' consensus matrix with average linkage, as returned by \code{\link{consensushc}(object)}.
#' This is also the ordering that is used by default for the heatmap of the consensus matrix 
#' as ploted by \code{\link{consensusmap}}. 
#' } 
setMethod('coefmap', signature(object='NMFfitX'),
	function(object
			, Colv=TRUE
			, annRow=NA, annCol=NA
			, tracks=c('basis', 'consensus:')
			, ...){
		
		x <- minfit(object)
		
		# process annotation tracks
		ptracks <- process_tracks(x, tracks, annRow, annCol)
		annRow <- ptracks$row 
		annCol <- ptracks$col
		# set special annotation handler
		specialAnnotation(2L, 'consensus', function() predict(object, what='consensus'))
		# row track handler is added in coefmap,NMF
		#
		
		## process ordering
		if( isString(Colv) ){
			if( Colv %in% c('consensus', 'cmap') )
				Colv <- consensushc(object, 'consensus')
		}
		##
		# call the method on the best fit
		coefmap(x, ..., Colv=Colv, annRow=annRow, annCol=annCol, tracks=NA)	
	}
)

#' Cophenetic Correlation Coefficient 
#'
#' The function \code{cophcor} computes the cophenetic correlation coefficient 
#' from consensus matrix \code{object}, e.g. as obtained from multiple NMF runs.
#' 
#' The cophenetic correlation coeffificient is based on the consensus matrix
#' (i.e. the average of connectivity matrices) and was proposed by 
#' \cite{Brunet2004} to measure the stability of the clusters obtained from NMF.
#' 
#' It is defined as the Pearson correlation between the samples' distances
#' induced by the consensus matrix (seen as a similarity matrix) and their
#' cophenetic distances from a hierachical clustering based on these very
#' distances (by default an average linkage is used).  
#' See \cite{Brunet2004}.
#' 
#' @param object an object from which is extracted a consensus matrix.
#' @param ... extra arguments to allow extension and passed to subsequent calls. 
#' 
#' @inline
#' @seealso \code{\link{cophenetic}}
#' @export
setGeneric('cophcor', function(object, ...) standardGeneric('cophcor') )
#' Workhorse method for matrices.
#' 
#' @param linkage linkage method used in the hierarchical clustering.
#' It is passed to \code{\link{hclust}}.
#' 
setMethod('cophcor', signature(object='matrix'),
	function(object, linkage='average'){
		
		# check for empty matrix
		if( nrow(object)==0  || ncol(object)==0 )
		{
			warning("NMF::cophcor - NA produced [input matrix is of dimension ", nrow(object), "x", ncol(object), "]"
					, call.=FALSE)
			return(NA)
		}
		
		# safe-guard for diagonal matrix: to prevent error in 'cor'
		if( all(object[upper.tri(object)]==0) && all(diag(object)==object[1,1]) )
			return(1)
		# convert consensus matrix into dissimilarities
		d.consensus <- as.dist(1 - object)
		# compute cophenetic distance based on these dissimilarities
		hc <- hclust(d.consensus, method=linkage)
		d.coph <- cophenetic(hc)
		
		# return correlation between the two distances
		res <- cor(d.consensus, d.coph, method='pearson')
		return(res)
	}
)
#' Computes the cophenetic correlation coefficient on the consensus matrix
#' of \code{object}.
#' All arguments in \code{...} are passed to the method \code{cophcor,matrix}.
setMethod('cophcor', signature(object='NMFfitX'),
	function(object, ...){
		# compute the consensus matrix
		C <- consensus(object)
		
		return( cophcor(C, ...))
	}
)

# TODO: uncomment this and make it compute the mean or best rss
#setMethod('rss', 'NMFfitXn', 
#	function(object, target, ...){
#		rss(fit(object, ...), target)
#	}
#)

