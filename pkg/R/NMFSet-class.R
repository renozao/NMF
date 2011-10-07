###% @include NMF-class.R
NA

###% Old NMFset class definition
###% 
###% The class wraps a list of NMFfit objects that results from a multiple NMF runs (of a single method)
setClass('NMFSet'
		, representation(
				consensus = 'matrix' # average connectivity matrix used when storing the best result only
				, runtime = 'proc_time' # running time to perform all the NMF runs
				, nrun = 'integer'
		)
		, contains='list'
		, prototype=prototype(
				consensus =	matrix(NA,0,0)
				, nrun = as.integer(0)
		)
)

###% @details  \emph{isNMFfit} tells if an object results from an NMF fit. That is it
###% checks if \code{object} inherits from class \code{\linkS4class{NMFfit}} or
###% form class \code{\linkS4class{NMFfitX}}, which are returned by the function
###% \code{\link{nmf}}.  If \code{object} is a \code{list} and
###% \code{recursive=TRUE}, then the check is performed on each element of the
###% list, and the return value is a vector (or a list if \code{object} is a list
###% of list) of the same length as \code{object}.
###% 
###% @rdname advanced
###% @param object any R object.
###% @param recursive if \code{TRUE} and \code{object} is a list then the check
###% is performed on each element of the list. Note that the recursivity only
###% applies in the case of lists that are not themselves NMFfit objects, unlike
###% \code{NMFfitXn} objects for which the result of \code{isNMFfit} will always
###% be \code{TRUE} (a single logical value).
###% 
###% @return \code{isNMFfit} returns a \code{logical} vector (or a list if
###% \code{object} is a list of list) of the same length as \code{object}.
###% 
###% @seealso \code{\linkS4class{NMFfit}}, \code{\linkS4class{NMFfitX}}
###% @examples
###% 
###% 	# generate a random 50 x 10 matrix
###% 	V <- rmatrix(50, 10)
###% 	
###% 	# single run
###% 	res <- nmf(V, 3)
###% 	isNMFfit(res)
###% 	
###% 	# multiple runs - keeping single fit
###% 	resm <- nmf(V, 3, nrun=3)
###% 	isNMFfit(resm)
###% 	
###% 	# multiple runs - keeping all fits
###% 	resM <- nmf(V, 3, nrun=3, .opt='k') 
###% 	isNMFfit(resM)
###% 	
###% 	# with a list of results
###% 	isNMFfit(list(res, resm, resM, 'not a result'))
###% 	isNMFfit(list(res, list(resm, resM), 'not a result')) # list of list
###% 	isNMFfit(list(res, resm, resM, 'not a result'), recursive=FALSE)
###% 	
###% 
isNMFfit <- function(object, recursive=TRUE){
	res <- is(object, 'NMFfitX') || is(object, 'NMFfit')
	# if the object is not a NMF result: apply to each element if a list (only in recursive mode)
	if( !res  && recursive && is.list(object) )
		sapply(object, isNMFfit)
	else
		res
}

###% NMFlist class definition
###% 
###% The class wraps a list of results of NMF runs.
###% These can be either from a single run (NMFfit) or multiple runs (NMFfitX).
###% It original aim is to hold NMF results from different methods  
setClass('NMFList'
		, representation(
			runtime='proc_time'
		)
		, contains='namedList'
		, validity=function(object){
			
			# the list must only contains NMFfit objects of the same dimensions
			for(i in seq_along(object)){
				# check class of the element
				item <- object[[i]]
				if( !(is(item, 'NMFfitX') || is(item, 'NMFfit')) )
					return(paste("invalid class for element", i, "of input list [all elements must be a NMFfit or NMFfitX object]"))
			}
			
		}
)

###% Show method for an NMFout object
setMethod('show', 'NMFList', 
	function(object)
	{
		cat("<Object of class:", class(object), ">\n")
		cat("Length:", length(object), "\n")
		if( length(object) > 0 ) cat("Method(s):", algorithm(object, collapse=TRUE), "\n")
		# show totaltime if present
		tt <- runtime(object)
		if( length(tt) > 0 ){
			cat("Total timing:\n"); show(tt);
		}
	}
)

###% Returns the name of the method(s) used to compute the element(s) in the list.
setMethod('algorithm', 'NMFList', 
	function(object, collapse=FALSE){
		l <- length(object)
		
		if( l == 0 ) NULL
		else if( l == 1 ) algorithm(object[[1]])
		else{	
			# build the vector of the algorithm names (with no repeat)
			m <- unique(sapply(object, algorithm))
			if( collapse )
				paste(m, collapse=', ')
			else m
		}
	} 
)

setMethod('runtime', 'NMFList', 
	function(object){
		t <- slot(object, 'runtime')
		if( length(t)==0 ) NULL else t
	}
)

setGeneric('as.NMFList', function(...) standardGeneric('as.NMFList') )
setMethod('as.NMFList', 'ANY' 
		, function(..., unlist=FALSE){
			arg.l <- list(...)
			if( length(arg.l) == 1 && is.list(arg.l[[1]]) && !is(arg.l[[1]], 'NMFfitX') )
				arg.l <- arg.l[[1]]
			
			# unlist if required
			if( unlist )
				arg.l <- unlist(arg.l)
			
			# create a NMFList object from the input list 
			new('NMFList', arg.l)
		}
)
			

###% NMFfitX class definition
###% 
###% Virtual class for the result from a multiple NMF runs (of a single method)
###% Its objective is to provide a common interface for the result of multiple runs 
###% wether it holds all the fits or only the best one 
setClass('NMFfitX'
		, representation(
				runtime.all = 'proc_time' # running time to perform all the NMF runs
		)
		, contains='VIRTUAL'		
)

setMethod('runtime.all', 'NMFfitX', 
		function(object){	
			t <- slot(object, 'runtime.all')
			if( length(t) > 0 ) t else NULL
		}
)

setGeneric('nrun', function(object, ...) standardGeneric('nrun') )
setMethod('nrun', 'matrix', 
	function(object){
		attr(object, 'nrun')
	}
)
setMethod('nrun', 'NMFfitX', 
		function(object){	
			stop("NMF::NMFfitX - missing definition for pure virtual method 'nrun' in class '", class(object), "'")
		}
)
###% Dummy 'nrun' method for NMFfit objects: returns 1 (i.e. the number of runs used to compute the fit)
setMethod('nrun', 'NMFfit', 
	function(object){
		1
	}
)

setGeneric('consensus', function(object, ...) standardGeneric('consensus') )
setMethod('consensus', 'NMFfitX', 
	function(object, ...){	
		stop("NMF::NMFfitX - missing definition for pure virtual method 'consensus' in class '", class(object), "'")
	}
)

#dummy 'concensus' method for class 'NMFfit'
setMethod('consensus', 'NMF', 
	function(object){	
		connectivity(object)
	}
)

#' Hierarchical Clustering of a Consensus Matrix
#' 
#' The function \code{consensustree} computes the hierarchical clustering of 
#' a consensus matrix, using the matrix itself as a similarity matrix and 
#' average linkage.
#' 
#' @param object a matrix or an \code{NMFfitX} object, as returned by multiple 
#' NMF runs. 
#' @param dendrogram a logical that specifies if the result of the hierarchical 
#' clustering (en \code{hclust} object) should be converted into a dendrogram. 
#' Default value is \code{TRUE}.
#' 
#' @return an object of class \code{dendrogram} or \code{hclust} depending on the  
#' value of argument \code{dendrogram}.   
#' 
#' @export
setGeneric('consensustree', function(object, ...) standardGeneric('consensustree'))
#' @autoRd
setMethod('consensustree', 'matrix', 
	function(object, dendrogram=TRUE, ...){
		
		# hierachical clustering based on the connectivity matrix
		hc <- hclust(as.dist(1-object), method='average')
		
		# convert into a dendrogram if requested
		if( dendrogram ) as.dendrogram(hc)
		else hc
	}
)
#' @autoRd
setMethod('consensustree', 'NMF', 
	function(object, ...){		
		# hierachical clustering based on the connectivity matrix
		coltree(connectivity(object), ...)		
	}
)
#' @autoRd
setMethod('consensustree', 'NMFfitX', 
	function(object, which=c('fit', 'consensus'), ...){
		
		which <- match.arg(which)
		if( which == 'consensus' ){
			# hierachical clustering on the consensus matrix
			coltree(consensus(object), ...)
						
		}else if( which == 'fit' )
			coltree(fit(object), ...)
		
	}
)

setMethod('predict', signature(object='NMFfitX'),
	function(object, what=c('columns', 'rows', 'samples', 'features', 'consensus', 'cmap'), ...){
		# determine which prediction to do
		what <- match.arg(what)
		if( what=='consensus' || what=='cmap' ){
			# build the tree from consensus matrix
			h <- hclust(as.dist(1-consensus(object)), method='average')
			# extract membership from the tree
			cl <- cutree(h, k=nbasis(object))
			
			# reorder the levels in the case of consensus map
			if( what=='cmap' )
				cl <- match(cl, unique(cl[h$order]))
			as.factor(cl)
		}
		else predict(fit(object), what=what, ...)
	}
)

setMethod('fit', 'NMFfitX', 
	function(object){	
		stop("NMF::NMFfitX - missing definition for pure virtual method 'fit' in class '", class(object), "'")
	}
)

setMethod('minfit', 'NMFfitX',
	function(object){
		stop("NMF::NMFfitX - missing definition for pure virtual method 'minfit' in class '", class(object), "'")
	}
)

setMethod('show', 'NMFfitX', 
		function(object){
			cat("<Object of class:", class(object), ">\n")
			# name of the algorithm
			cat("  Method:", algorithm(object), "\n")
			# number of runs
			cat("  Runs: ", nrun(object),"\n");
			# initial state
			cat("  RNG:\n  ", RNGdesc(getRNG1(object)),"\n");
			if( nrun(object) > 0 ){
				# show total timing			
				cat("  Total timing:\n"); show(runtime.all(object));
			}
		}
)

setMethod('getRNG1', signature(object='NMFfitX'),
	function(object){
		stop("NMF::getRNG1(", class(object), ") - Unimplemented pure virtual method: could not extract initial RNG settings.")
	}
)

###% compare two NMF models when at least one comes from a NMFfitX object 
setMethod('nmf.equal', signature(x='NMFfitX', y='NMF'), 
	function(x, y, ...){
		nmf.equal(fit(x), y, ...)
	}
)
setMethod('nmf.equal', signature(x='NMF', y='NMFfitX'), 
		function(x, y, ...){
			nmf.equal(x, fit(y), ...)
		}
)

#########################################################
# END_NMFfitX
#########################################################


###% NMFfitX1 class definition
###% 
###% The class holds a single NMFfit object that is the best result from multiple NMF runs (of a single method)
setClass('NMFfitX1'
	, representation(
			#fit = 'NMFfit' # holds the best fit from all the runs
			consensus = 'matrix' # average connectivity matrix of all the NMF runs
			, nrun = 'integer'
			, rng1 = 'ANY'
	)
	, contains=c('NMFfitX', 'NMFfit')
	, prototype=prototype(
			consensus =	matrix(NA,0,0)
			, nrun = as.integer(0)
	)
)



###% Show method for an NMFfitX1 object
setMethod('show', 'NMFfitX1', 
	function(object){
		callNextMethod(object)
		
		# show details of the best fit
		#cat(" # Best fit:\n  ")
		#s <- capture.output(show(fit(object)))
		#cat(s, sep="\n  |")
	}
)

setMethod('nrun', 'NMFfitX1', 
	function(object){
		slot(object,'nrun')
	}
)

setMethod('consensus', signature(object='NMFfitX1'), 
	function(object){
		
		C <- slot(object, 'consensus')
		if( length(C) > 0 ){
			class(C) <- c(class(C), 'NMF.consensus')
			attr(C, 'nrun') <- nrun(object)
			attr(C, 'nbasis') <- nbasis(object)
			C
		}else NULL
		
	}
)

###% Extract the best NMF fit from the object.
###%
###% @param x a \code{NMFfitX1} object from which to extract the best result
###% 
setMethod('minfit', 'NMFfitX1',
	function(object){	
		# coerce the object into a NMFfit object
		as(object, 'NMFfit')
	}
)

###% Extract the best NMF model from the object.
###%
###% @param x a \code{NMFfitX1} object from which to extract the best result
###%
.function.gen <- function(){	
	ncall <- 0
	function(object){
		if( ncall < 1 ){
			warning("NMF - Change in version 0.5.3: Method 'fit' for 'NMFfitX1' objects returns the best 'NMF' model object. Use 'minfit' to get the best 'NMFfit' object as before."
					, call.=FALSE)
			ncall <<- ncall + 1
		}
		slot(object, 'fit')
	}
}
setMethod('fit', signature(object='NMFfitX1'), .function.gen())
rm('.function.gen')


###% Get RNG Settings for first run 
setMethod('getRNG1', signature(object='NMFfitX1'),
	function(object){
		object@rng1
	}
)

###% nmf.equal for NMFfitX1: compare the best fitted models
setMethod('nmf.equal', signature(x='NMFfitX1', y='NMFfitX1'), 
		function(x, y, ...){
			nmf.equal(fit(x), fit(y), ...)
		}
)


#########################################################
# END_NMFfitX1
#########################################################


###% NMFfitXn class definition
###% 
###% The class holds a list of NMFfit objects from multiple NMF runs (of a single method)
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

###% Method 'nbasis' for 'NMFfitXn' objects: 
###% 
setMethod('nbasis', signature(x='NMFfitXn'), 
	function(x, ...){
		if( length(x) == 0 ) return(NULL)
		return( nbasis(x[[1]]) )
	}
)
###% Method 'dim' for 'NMFfitXn' objects: 
###% 
setMethod('dim', signature(x='NMFfitXn'), 
	function(x){
		if( length(x) == 0 ) return(NULL)
		return( dim(x[[1]]) )
	}
)

###% Method 'coef' for 'NMFfitXn' objects: 
###% 
setMethod('coef', signature(object='NMFfitXn'), 
		function(object){
			coef(fit(object))
		}
)

###% Method 'basis' for 'NMFfitXn' objects: 
###% 
setMethod('basis', signature(object='NMFfitXn'), 
	function(object){
			basis(fit(object))
	}
)


###% Returns the number of runs stored in the list (i.e. the list length)
setMethod('nrun', 'NMFfitXn', 
	function(object){
		length(object)
	}
)

setMethod('algorithm', 'NMFfitXn', 
	function(object){
		if( length(object) == 0 ) return(NULL)
		return( algorithm(object[[1]]) )
	} 
)

setMethod('seeding', 'NMFfitXn',
	function(object){
		if( length(object) == 0 ) return(NULL)
		return( seeding(object[[1]]) )
	}
)

setMethod('modelname', signature(object='NMFfitXn'), 
	function(object){
		if( length(object) == 0 ) return(NULL)
		return( modelname(object[[1]]) )
	}
)


setGeneric('seqtime', function(object, ...) standardGeneric('seqtime') )
setMethod('seqtime', 'NMFfitXn', 
	function(object){	
		if( length(object) == 0 ) return(NULL)
		# sum up the time accross the runs
		t.mat <- sapply(object, runtime)
		res <- rowSums(t.mat)
		class(res) <- 'proc_time'
		res
	}
)

###% Returns the total time used to perform all the runs.
###% If no time data has been set in slot 'runtime.all' then the total sequential time is returned
###% 
###% @seealso seqtime
###% 
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


###% Returns the best NMF model in the list, i.e. the run that have the lower estimation residuals.
###%
###% @param x a \code{NMFfitXn} object from which to extract the best result
###% 
setMethod('minfit', 'NMFfitXn',
	function(object){
		
		b <- which.best(object, deviance)
		# test for length 0
		if( length(b) == 0 ) return(NULL)
				
		# return the run with the lower
		object[[ b ]]
	}
)

# Returns the index of the best fit
which.best <- function(object, FUN=deviance, ...){
	
	# test for length 0
	if( length(object) == 0 ) 
		return(integer())
	
	# retrieve the measure for each run
	e <- sapply(object, FUN, ...)
	
	# return the run with the lower
	which.min(e)
}

###% Get RNG Settings 
setMethod('getRNG1', signature(object='NMFfitXn'),
	function(object){
		if( length(object) == 0 )
			stop("NMF::getRNG1 - Could not extract RNG data from empty object [class:", class(object), "]")
		
		getRNG(object[[1]])
	}
)
setMethod('getRNG', signature(object='NMFfitXn'),
	function(object){
		getRNG(minfit(object))
	}
)


###% Returns the best NMF run in the list, i.e. the run that have the lower estimation residuals.
###%
###% @param x a \code{NMFfitXn} object from which to extract the best result
###%
.function.gen <- function(){	
	ncall <- 0
	function(object){
		if( ncall < 1 ){
			warning("NMF - Change in version 0.6: Method 'fit' for 'NMFfitXn' objects returns the best 'NMF' model object. Use 'minfit' to get the best 'NMFfit' object as before."
					, call.=FALSE)
			ncall <<- ncall + 1
		}
		fit( minfit(object) )
	}
}
setMethod('fit', signature(object='NMFfitXn'), .function.gen())
#rm(.function.gen)

setMethod('dim', 'NMFfitXn', 
		function(x){
			if( length(x) == 0 ) return(NULL)
			dim(x[[1]])
		} 
)

###% Compare the Results of Multiple NMF Runs
###% 
###% Either compare the two best fit, or the result of of each run. 
###%  
setMethod('nmf.equal', signature(x='list', y='list'), 
	function(x, y, all=FALSE, vector=FALSE, ...){
		if( !all )
			nmf.equal(x[[ which.best(x) ]], y[[ which.best(y) ]], ...)
		else{
			if( length(x) != length(y) )
				FALSE
			else
				res <- mapply(function(a,b) isTRUE(nmf.equal(a,b)), x, y, MoreArgs=list(...))
				if( !vector )
					res <- all( res )
				res
		}
	}
)

###% Computes the consensus matrix, i.e. the average connectivity matrix of a set of NMF runs
###%
###% The consensus matrix is defined as:
###% \deqn{C_0 = \frac{1}{N} \sum_{r=1}^N C_r,}
###% where \eqn{C_r} is the connectivity matrix of run \eqn{r}, and \eqn{N} is the total number of runs.
###%
###% A perfect consensus matrix (all entries = 0 or 1) means that the same clusters are found by each run (i.e. independently of the initialization). 
###% The entries of the consensus matrix reflect the probability for each pair of samples to belong to the same cluster.
###%
setMethod('consensus', signature(object='NMFfitXn'), 
	function(object){
		if( length(object) == 0 ) return(NULL)
		
		# init empty consensus matrix
		con <- matrix(0, ncol(object), ncol(object))
		# name the rows and columns appropriately: use the sample names of the first fit
		dimnames(con) <- list(colnames(object[[1]]), colnames(object[[1]]))
		
		# compute mean connectivity matrix
		sapply(object 
				, function(x){
					con <<- con + connectivity(x)
					NULL
				}
		)
		con <- con / nrun(object)
				
		# return result
		class(con) <- c(class(con), 'NMF.consensus')
		attr(con, 'nrun') <- nrun(object)
		attr(con, 'nbasis') <- nbasis(object)
		con
	}
)

###% plot function for consensus matrices
plot.NMF.consensus <- function(x, ...){
	consensusmap(x, ...)
}

###% Computes the dispersion of the consensus matrix associated to a set of NMF run.
###%
###% The dispersion coeffificient of a consensus matrix (i.e. the average of connectivity matrices) is a measure of reproducibility of the clusters.
###% The dispersion coeffificient is given by:
###% \deqn{\rho = \sum_{i,j=1}^n 4 (C_{ij} - \frac{1}{2})^2 .}
###%, where \eqn{n} is the total number of samples.
###%
###% We have \eqn{0 \leq \rho \leq 1} and \eqn{\rho = 1} only for a perfect consensus matrix, where all entries 0 or 1. A perfect consensus matrix is obtained only when all the
###% the connectivity matrices are the same, meaning that the algorithm gave the same clusters at each run.
###%
###% @param x a consensus matrix.
###% @return the dispersion coefficient -- as a numeric value.
###%
###% @references Kim, H. & Park, H. 
###%	Sparse non-negative matrix factorizations via alternating non-negativity-constrained least squares for microarray data analysis.
###%	Bioinformatics (2007). 
###%	\url{http://dx.doi.org/10.1093/bioinformatics/btm134}.
###%	
setGeneric('dispersion', function(object, ...) standardGeneric('dispersion') )
setMethod('dispersion', signature(object='matrix'), 
	function(object, ...){
		stopifnot( nrow(object) == ncol(object) )
		sum( 4 * (object-1/2)^2 ) / nrow(object)^2
	}
)

setMethod('dispersion', signature(object='NMFfitX'), 
	function(object, ...){
		dispersion(consensus(object), ...)
	}
)

setGeneric('join', function(object, ...) standardGeneric('join') )
setMethod('join', 'list',
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
					stop("NMF::join - invalid value for 'runtime.all' [5-length numeric expected]")
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
					stop("NMF::join - invalid class for element ", i, " of input list [all elements must be NMFfit or NMFfitX objects]")
				
				# check that all elements result from the same algorithm
				if( is.null(ref.algo) ) ref.algo <<- algorithm(item)
				if( !identical(algorithm(item), ref.algo) )
					stop("NMF::join - invalid algorithm for element ", i, " of input list [cannot join results from different algorithms]")
				
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
			nmf.debug('join', ".merge is forced to TRUE")
			.merge <- TRUE
		}
		
		# unpack all the NMFfit objects
		object.list <- unlist(object)
		nmf.debug('join', "Number of fits to join = ", length(object.list))
					
		# one wants to keep only the best result
		if( .merge ){
			
			warning("NMF::join - The method for merging lists is still in development")
			
			# set the total number of runs
			extra$nrun <- as.integer(nrun)			
									
			# consensus matrix
			if( !is.null(extra$consensus) )
				warning("NMF::join - the value of 'consensus' was discarded as slot 'consensus' is computed internally")
			extra$consensus <- NULL
									
			consensus <- matrix(NA, 0, 0)
			best.res <- Inf		
			best.fit <- NULL
			sapply(object.list, function(x){
				if( !is(x, 'NMFfit') )
					stop("NMF::join - all inner-elements of '",substitute(object),"' must inherit from class 'NMFfit'")
				
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
			return( do.call(join, c(list(best.fit), extra)) )
		}
		else{
			# create a NMFfitXn object that holds the whole list			
			do.call('new', c(list('NMFfitXn', object.list), extra))
		}
	}
)
		
setMethod('join', 'NMFfit',
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

setMethod('join', 'NMFfitX',
		function(object, ...){

			# nothing to do in the case of NMFfitX1 objects
			if( is(object, 'NMFfitX1') ) return(object)
			
			# retrieve extra arguments
			extra <- list(...)
			
			# take runtime.all from the object itself
			if( !is.null(extra$runtime.all) )
				warning("NMF::join - argument 'runtime.all' was discarded as it is computed from argument 'object'")			
			extra$runtime.all <- runtime.all(object)						
			
			# create the NMFfitX1 object
			f <- selectMethod(join, 'list')
			do.call(f, c(list(object), extra))
		}
)

###% Computes the average purity of a set of NMF runs.
setMethod('purity', signature(x='NMFfitXn', class='ANY'), 
	function(x, class, method='best', ...){
		c <- sapply(x, purity, class=class, ...)
		
		# aggregate the results if a method is provided
		if( is.null(method) ) c
		else aggregate.measure(c, method, decreasing=TRUE)		
	}
)

###% Computes the average entropy of a set of NMF runs.
setMethod('entropy', signature(x='NMFfitXn', class='ANY'), 
	function(x, class, method='best', ...){
		c <- sapply(x, entropy, class=class, ...)		
		
		# aggregate the results if a method is provided
		if( is.null(method) ) c
		else aggregate.measure(c, method)
	}
)

###% Computes the average final residuals of a set of NMF runs.
setMethod('residuals', signature(object='NMFfitXn'), 
	function(object, method='best', ...){
		e <- sapply(object, residuals, ...)
		
		# aggregate the results if a method is provided
		if( is.null(method) ) e
		else aggregate.measure(e, method)
	}
)
###% Returns the deviance of a fitted NMF model
setMethod('deviance', 'NMFfitXn',
	function(object){
		setNames(residuals(object), NULL)
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


###% Returns the cluster prediction defined by the best fit
setMethod('predict', signature(object='NMFfitXn'),
	function(object, ...){
		predict(fit(object), ...)
	}
)

###% Summary method for class NMFfitX
###% Computes the summary for the best fit and add some extra measures specific to multiple runs
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
		
		# return result
		s
	}
)
###% Compare different runs of NMF.
###% 
###% This function compares the factorizations obtained by different runs of NMF. This is typically usefull
###% to evaluate and compare how different algorithms perform on the same data.
###% 
###% @return a \code{data.frame} with the different comparison criteriae in rows and the methods in column
###%  
setGeneric('compare', function(object, ...) standardGeneric('compare') )
setMethod('compare', signature(object='list'),
	function(object, ..., unlist=FALSE){
		# wrap up x into a NMFList object if necessary
		if( !is(object, 'NMFList') ) object <- as.NMFList(object, unlist=unlist)
		summary(object, ...)
	}
)

setMethod('summary', signature(object='NMFList'),
	function(object, sort.by=NULL, select=NULL, ...){
		
		# define the sorting schema for each criteria (TRUE for decreasing, FALSE for increasing)
		sorting.schema <- list(method=FALSE, seed=FALSE, rng=FALSE, metric=FALSE
							, residuals=FALSE, cpu=FALSE, purity=TRUE, nrun=FALSE, cpu.all=FALSE
							, cophenetic=TRUE, dispersion=TRUE #NMFfitX only
							, entropy=FALSE, sparseness.basis=TRUE, sparseness.coef=TRUE, rank=FALSE, rss=FALSE
							, niter=FALSE, evar=TRUE)
				
		# for each result compute the summary measures
		measure.matrix <- sapply(object, summary, ...)		
		
		# the results from 'summary' might not have the same length => generate NA where necessary
		if( is.list(measure.matrix) ){
			name.all <- unique(unlist(sapply(measure.matrix, names)))			
			measure.matrix <- sapply(names(measure.matrix),
				function(name){
					m <- measure.matrix[[name]][name.all]
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
					svalue <- if( is.function(svalue) ) '<function>' else paste("'", svalue,"'", sep='')
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

###% Compare objective value trajectories from different NMF algorithms.
setMethod('plot', signature(x='NMFList', y='missing'), 
	function(x, y, ...){
		
		# retrieve normalized residuals tracks
		max.iter <- 0
		tracks <- lapply( x, 
				function(res){
					res <- minfit(res)
					t <- residuals(res, track=TRUE)
					#t <- t[-1]
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
				, xlab='Iterations', ylab='Normalized objective values'
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

setMethod('metaHeatmap', signature(object='NMFfitX'),
		function(object, ...){
			# send deprecated warning
			.Deprecated('metaHeatmap', 'NMF', "Direct use of the S4-Method 'metaHeatmap' for 'NMFfitX' objects is deprecated, use 'consensusmap' instead.")

			# call the new function 'consmap'
			return( consensusmap(object, ...) )
			
		}
)

setGeneric('consensusmap', function(object, ...) standardGeneric('consensusmap') )
setMethod('consensusmap', 'NMFfitX', 
	function(object, tracks = c('basis', 'consensus'), info = FALSE, ...){
	
		# retreive the graphical parameters and match them to the sub-sequent call to 'aheatmap'
		graphical.params <- list(...)
		names(graphical.params) <- .match.call.args(names(graphical.params), 'aheatmap', call='NMF::consensusmap')
		
		# add side information if requested
		info <- if( isTRUE(info) ){
					paste("NMF model: '", modelname(object)
					, "'\nAlgorithm: '", algorithm(object)
					, "'\nbasis: ", nbasis(object)
					,"\nnrun: ", nrun(object), sep='')
				}else if( isFALSE(info) ) NULL
				else info
		
		graphical.params <- .set.list.defaults(graphical.params
			, info = info
			, main = 'Consensus matrix'
			)
			
		# add annotation tracks
		if( length(tracks) > 0 && !isNA(tracks) ){
			
			# create extra annotation tracks
			tr <- sapply( tracks, function(t){
				switch(t
						, consensus = predict(object, 'cmap')
						, basis = predict(object)
						, NA
				)
			}, simplify=FALSE)
			
			tracks <- tracks[is.na(tr)]
			# add to the other annotations 
			graphical.params[['annCol']] <- atrack(tr, graphical.params[['annCol']])			
		}		
			
		x <- consensus(object)
		do.call('consensusmap', c(list(x), graphical.params))	
	}
)

setMethod('consensusmap', 'NMF', 
	function(object, ...){
		consensusmap(connectivity(object), ...)		
	}
)

setMethod('consensusmap', 'matrix', 
	function(object, info = FALSE, ...){
		
		# retreive the graphical parameters and match them to the sub-sequent call to 'aheatmap'
		graphical.params <- list(...)
		names(graphical.params) <- .match.call.args(names(graphical.params), 'aheatmap', call='NMF::consensusmap')
		
		nr <- nrun(object)
		nb <- nbasis(object)
		info <- if( isTRUE(info) ){
					info <- NULL
					if( !is.null(nr) ) info <- c(info, paste("nrun:", nr))
					if( !is.null(nb) ) info <- c(info, paste("nbasis:", nb))
					info <- c(info, paste("cophcor:", round(cophcor(object), 3)))
				}else if( isFALSE(info) ) NULL
				else info
			
		# set default graphical parameters for type 'consensus'
		graphical.params <- .set.list.defaults(graphical.params
				, distfun = function(x) as.dist(1-x)
				, hclustfun = 'average'
				, Rowv = TRUE, Colv = "Rowv"						
				, color='-RdYlBu'
				, scale='none'
				, main = if( is.null(nr) || nr > 1 ) 'Consensus matrix' else 'Connectiviy matrix'
				, info = info
		)
		
		do.call('aheatmap', c(list(object), graphical.params))	
	}
)

setOldClass('NMF.rank')
setMethod('consensusmap', 'NMF.rank', 
	function(object, ...){

		# plot the list of consensus matrix (set names to be used as default main titles)
		consensusmap(setNames(object$consensus, paste("rank = ", lapply(object$fit, nbasis))), ...)
	}
)

setMethod('consensusmap', 'list', 
	function(object, layout, ...){
		
		# retreive the graphical parameters and match them to the sub-sequent call to 'aheatmap'
		graphical.params <- list(...)
		names(graphical.params) <- .match.call.args(names(graphical.params), 'aheatmap', call='NMF::consensusmap')
		
		# set default graphical parameters for type 'consensus'
		graphical.params <- .set.list.defaults(graphical.params
				, Rowv = FALSE # reorder but do not show the dendrograms
				, main = names(object) # default is to use the list's names as main titles
		)
		
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
		if( is.numeric(layout) ){
			if( length(layout) == 1 )
				layout <- c(layout, layout)
			layout <- matrix(1:(layout[1]*layout[2]), layout[1], byrow=TRUE)
		}
		
		graphics::layout(layout)
		res <- sapply(seq_along(object), function(i){
			x <- object[[i]]
			
			# set main title
			main <- if( !is.null(graphical.params$main) && length(graphical.params$main) > 1 ){
				if( length(graphical.params$main) != length(object) )
					stop("consensusmap - Invalid length for argument `main`: should be either a single character string, or a list or vector of same length as ", deparse(substitute(object)))
				graphical.params$main <- graphical.params$main[[i]]
			}			
			
			# call method for the fit
			do.call('consensusmap', c(list(x), graphical.params))
		})
		invisible(res)
	}
)

setMethod('basismap', signature(object='NMFfitX'),
	function(object, ...){
		# call the method on the best fit
		basismap(minfit(object), ...)	
	}
)

setMethod('coefmap', signature(object='NMFfitX'),
	function(object, tracks=c('basis', 'consensus'), ...){
		
		# retreive the graphical parameters and match them to the sub-sequent call to 'heatmap.plus.2'
		graphical.params <- list(...)	
		names(graphical.params) <- .match.call.args(names(graphical.params), 'aheatmap', call='NMF::coefmap')
		
		# add annotation tracks
		if( length(tracks) > 0 && !isNA(tracks) ){
			# create extra annotation tracks
			tr <- sapply( tracks, function(t){
				switch(t
					, consensus = predict(object, 'cmap')			
					, NA
				)
			}, simplify=FALSE)
			
			nt <- names(tr)
			names(tr)[nt==''] <- tracks[nt=='']
			tracks <- tracks[is.na(tr)]
			# add to the other annotations 
			graphical.params[['annCol']] <- atrack(tr, graphical.params[['annCol']])			
		}		
		# call the method on the best fit
		do.call('coefmap', c(list(minfit(object), tracks=tracks), graphical.params))	
	}
)

###% Computes cophenetic correlation coefficient
setGeneric('cophcor', function(object, ...) standardGeneric('cophcor') )
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

