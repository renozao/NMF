#' @include NMF-class.R
NA

#' Old NMFset class definition
#' 
#' The class wraps a list of NMFfit objects that results from a multiple NMF runs (of a single method)
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


isNMFfit <- function(object, recursive=TRUE){
	res <- is(object, 'NMFfitX') || is(object, 'NMFfit')
	# if the object is not a NMF result: apply to each element if a list (only in recursive mode)
	if( !res  && recursive && is.list(object) )
		sapply(object, isNMFfit)
	else
		res
}

#' NMFlist class definition
#' 
#' The class wraps a list of results of NMF runs.
#' These can be either from a single run (NMFfit) or multiple runs (NMFfitX).
#' It original aim is to hold NMF results from different methods  
setClass('NMFList'
		, representation(
			runtime='proc_time'
		)
		, contains='list'
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

#' Show method for an NMFout object
setMethod('show', signature(object='NMFList'), 
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

#' Returns the name of the method(s) used to compute the element(s) in the list.
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

if ( !isGeneric("as.NMFList") ) setGeneric('as.NMFList', function(...) standardGeneric('as.NMFList') )
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
			

#' NMFfitX class definition
#' 
#' Virtual class for the result from a multiple NMF runs (of a single method)
#' Its objective is to provide a common interface for the result of multiple runs 
#' wether it holds all the fits or only the best one 
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

if ( !isGeneric("nrun") ) setGeneric('nrun', function(object, ...) standardGeneric('nrun') )
setMethod('nrun', 'NMFfitX', 
		function(object){	
			stop("NMF::NMFfitX - missing definition for pure virtual method 'nrun' in class '", class(object), "'")
		}
)
#' Dummy 'nrun' method for NMFfit objects: returns 1 (i.e. the number of runs used to compute the fit)
setMethod('nrun', 'NMFfit', 
	function(object){
		1
	}
)

if ( !isGeneric('consensus') ) setGeneric('consensus', function(object, ...) standardGeneric('consensus') )
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

if ( !isGeneric("fit") ) setGeneric('fit', function(object, ...) standardGeneric('fit') )
setMethod('fit', 'NMFfitX', 
	function(object){	
		stop("NMF::NMFfitX - missing definition for pure virtual method 'fit' in class '", class(object), "'")
	}
)

setMethod('show', 'NMFfitX', 
		function(object){
			cat("<Object of class:", class(object), ">\n")
			# name of the algorithm
			cat("  Method:", algorithm(object), "\n")
			# number of runs
			cat("  Runs: ", nrun(object),"\n");
			if( nrun(object) > 0 ){
				# show total timing			
				cat("  Total timing:\n"); show(runtime.all(object));
			}
		}
)

#########################################################
# END_NMFfitX
#########################################################


#' NMFfitX1 class definition
#' 
#' The class holds a single NMFfit object that is the best result from multiple NMF runs (of a single method)
setClass('NMFfitX1'
	, representation(
			#fit = 'NMFfit' # holds the best fit from all the runs
			consensus = 'matrix' # average connectivity matrix of all the NMF runs
			, nrun = 'integer'
	)
	, contains=c('NMFfitX', 'NMFfit')
	, prototype=prototype(
			consensus =	matrix(NA,0,0)
			, nrun = as.integer(0)
	)
)



#' Show method for an NMFfitX1 object
setMethod('show', 'NMFfitX1', 
	function(object){
		cat("<Object of class:", class(object), ">\n")
		# name of the algorithm
		cat("  Method:", algorithm(object), "\n")
		# number of runs
		cat("  Runs: ", nrun(object),"\n");
		# show total timing
		cat("  Total timing:\n"); show(runtime.all(object));
		
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
		if( length(C) > 0 ) C
		else NULL
		
	}
)

#' Extract the best NMF fit from the object.
#'
#' @param x a \code{NMFfitX1} object from which to extract the best result
#' 
setMethod('fit', signature(object='NMFfitX1'),
	function(object){	
		# coerce the object into a NMFfit object
		as(object, 'NMFfit')
	}
)

#########################################################
# END_NMFfitX1
#########################################################


#' NMFfitXn class definition
#' 
#' The class holds a list of NMFfit objects from multiple NMF runs (of a single method)
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

#' Returns the number of runs stored in the list (i.e. the list length)
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

if ( !isGeneric("seqtime") ) setGeneric('seqtime', function(object, ...) standardGeneric('seqtime') )
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

#' Returns the total time used to perform all the runs.
#' If no time data has been set in slot 'runtime.all' then the total sequential time is returned
#' 
#' @seealso seqtime
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

#' Returns the best NMF run in the list, i.e. the run that have the lower estimation residuals.
#'
#' @param x a \code{NMFfitXn} object from which to extract the best result
#' 
setMethod('fit', signature(object='NMFfitXn'),
	function(object){
		
		# test for length 0
		if( length(object) == 0 ) return(NULL)
		
		# retirieve the estimation residuals for each run
		e <- sapply(object, residuals)
		
		# return the run with the lower
		object[[ which.min(e) ]]
	}
)

setMethod('dim', 'NMFfitXn', 
		function(x){
			if( length(x) == 0 ) return(NULL)
			return( dim(x[[1]]) )
		} 
)

#' Computes the consensus matrix, i.e. the average connectivity matrix of a set of NMF runs
#'
#' The consensus matrix is defined as:
#' \deqn{C_0 = \frac{1}{N} \sum_{r=1}^N C_r,}
#' where \eqn{C_r} is the connectivity matrix of run \eqn{r}, and \eqn{N} is the total number of runs.
#'
#' A perfect consensus matrix (all entries = 0 or 1) means that the same clusters are found by each run (i.e. independently of the initialization). 
#' The entries of the consensus matrix reflect the probability for each pair of samples to belong to the same cluster.
#'
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
		con
	}
)


#' Computes the dispersion of the consensus matrix associated to a set of NMF run.
#'
#' The dispersion coeffificient of a consensus matrix (i.e. the average of connectivity matrices) is a measure of reproducibility of the clusters.
#' The dispersion coeffificient is given by:
#' \deqn{\rho = \sum_{i,j=1}^n 4 (C_{ij} - \frac{1}{2})^2 .}
#', where \eqn{n} is the total number of samples.
#'
#' We have \eqn{0 \leq \rho \leq 1} and \eqn{\rho = 1} only for a perfect consensus matrix, where all entries 0 or 1. A perfect consensus matrix is obtained only when all the
#' the connectivity matrices are the same, meaning that the algorithm gave the same clusters at each run.
#'
#' @param x a consensus matrix.
#' @return the dispersion coefficient -- as a numeric value.
#'
#' @references Kim, H. & Park, H. 
#'	Sparse non-negative matrix factorizations via alternating non-negativity-constrained least squares for microarray data analysis.
#'	Bioinformatics (2007). 
#'	\url{http://dx.doi.org/10.1093/bioinformatics/btm134}.
#'	
if ( is.null(getGeneric("dispersion")) ) setGeneric('dispersion', function(object, ...) standardGeneric('dispersion') )
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

if ( is.null(getGeneric("join")) ) setGeneric('join', function(object, ...) standardGeneric('join') )
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
					if( is(x, 'NMFfitX1') ) best.fit <<- fit(x) # deal with the case of NMFfitX1 objects
					else best.fit <<- x
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

#' Computes the average purity of a set of NMF runs.
if ( is.null(getGeneric('purity')) ) setGeneric('purity', function(x, class, ...) standardGeneric('purity') )
setMethod('purity', signature(x='NMFfitXn', class='ANY'), 
	function(x, class, method=NULL, ...){
		c <- sapply(x, purity, class=class, ...)
		
		# aggregate the results
		aggregate.measure(c, method, decreasing=TRUE)		
	}
)

#' Computes the average entropy of a set of NMF runs.
if ( is.null(getGeneric('entropy')) ) setGeneric('entropy', function(x, class, ...) standardGeneric('entropy') )
setMethod('entropy', signature(x='NMFfitXn', class='ANY'), 
	function(x, class, method=NULL, ...){
		c <- sapply(x, entropy, class=class, ...)		
		
		# aggregate the results
		aggregate.measure(c, method)
	}
)

#' Computes the average final residuals of a set of NMF runs.
if( !isGeneric('residuals') ) setGeneric('residuals', package='stats')
setMethod('residuals', signature(object='NMFfitXn'), 
	function(object, method=NULL, ...){
		e <- sapply(object, residuals, ...)
		
		# aggregate the results
		aggregate.measure(e, method)		
	}
)

#' Utility function to aggregate numerical quality measures from \code{NMFfitXn} objects.
#' 
#' Given a numerical vector, this function computes an aggregated value using one of the following methods:
#' - mean: the mean of the measures
#' - best: the best measure according to the specified sorting order (decreasing or not)
#'  
aggregate.measure <- function(measure, method=c('mean', 'best'), decreasing=FALSE){
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


#' Returns the cluster prediction defined by the best fit
setMethod('predict', signature(object='NMFfitXn'),
	function(object, ...){
		predict(fit(object), ...)
	}
)

#' Summary method for class NMFfitX
#' Computes the summary for the best fit and add some extra measures specific to multiple runs
setMethod('summary', signature(object='NMFfitX'),
	function(object, ...){
		
		# compute summary measures for the best fit
		best.fit <- fit(object)
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
#' Compare different runs of NMF.
#' 
#' This function compares the factorizations obtained by different runs of NMF. This is typically usefull
#' to evaluate and compare how different algorithms perform on the same data.
#' 
#' @return a \code{data.frame} with the different comparison criteriae in rows and the methods in column
#'  
if ( is.null(getGeneric('compare')) ) setGeneric('compare', function(object, ...) standardGeneric('compare') )
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
		sorting.schema <- list(method=FALSE, seed=FALSE, metric=FALSE
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
					if( inherits(x, 'NMFfitX') ) x <- fit(x)
					m <- algorithm(x)
					s <- seeding(x) 
					svalue <- objective(x)
					svalue <- if( is.function(svalue) ) '<function>' else paste("'", svalue,"'", sep='')
					c(method=m, seed=s, metric=svalue)
				}
		)
		methods <- t(methods)	
		res <- as.data.frame(methods, stringsAsFactors=FALSE)	
		
		# add the measures to the result		
		res <- cbind(res, measure.matrix)
				
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

#' Compare objective value trajectories from different NMF algorithms.
setMethod('plot', signature(x='NMFList', y='missing'), 
	function(x, y, ...){
		
		# retrieve normalized residuals tracks
		max.iter <- 0
		tracks <- lapply( x, 
				function(res){
					if( is(res, 'NMFfitX') ) res <- fit(res)
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
			
			x <- consensus(object)
			do.call('metaHeatmap', c(list(x, type='consensus'), list(...)))
		}
)

#' Computes cophenetic correlation coefficient
if ( !isGeneric('cophcor') ) setGeneric('cophcor', function(object, ...) standardGeneric('cophcor') )
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
		C <- consensus(object, 'samples')
		
		return( cophcor(C, ...))
	}
)

# TODO: uncomment this and make it compute the mean or best rss
#setMethod('rss', 'NMFfitXn', 
#	function(object, target, ...){
#		rss(fit(object, ...), target)
#	}
#)

