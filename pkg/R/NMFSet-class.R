#' @include NMF-class.R
NA

#' NMFset class definition
#' 
#' The class wraps a list of NMF objects and defines some aggregating methods.
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




#' Show method for an NMFout object
setMethod('show', signature(object='NMFSet'), 
	function(object)
	{
		cat("<Object of class:", class(object), ">\n")
		if( length(object) > 0 ) cat("method:", algorithm(object), "\n")
		cat("runs:", nrun(object), "\n")
		if( length(object) != nrun(object) ) cat("fits:", length(object), "\n")
		if( length(object) > 1 )
			cat("average residuals:", residuals(object, method='mean'), "\n")	
		if( length(runtime(object)) > 0 ){ 
			cat("Timing:\n"); show(runtime(object));
			cat("Avg. timing:\n"); show(runtime(object)/nrun(object));
		}
	}
)

setMethod('runtime', 'NMFSet', 
	function(object, ...){
		stored.time <- slot(object, 'runtime')
		# if there is some time stored, return it
		if( length(stored.time) > 0 ) return(stored.time)
		
		# otherwise sum the time accross the runs
		t.mat <- sapply(object, runtime)
		res <- rowSums(t.mat)
		class(res) <- 'proc_time'
		return(res)
	}
)

if ( is.null(getGeneric("nrun")) ) setGeneric('nrun', function(object, ...) standardGeneric('nrun') )
setMethod('nrun', 'NMFSet', 
		function(object, ...){
			l <- length(object)
			# use slot nrun if the length is one 
			# -> maybe multiple runs were performed and only one was kept
			if( l == 1 && length(slot(object,'nrun')) > 0 )
				l <- object@nrun
			
			return(l)
		}
)

setMethod('algorithm', 'NMFSet', 
	function(object){
		if( length(object) == 0 ) return(invisible(NULL))
		
		m <- sapply(object, algorithm)
		if( any(m != m[1]) ) return(paste(m, collapse=', '))
		return(m[1]) 
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

setMethod('dispersion', signature(object='NMFSet'), 
	function(object, ...){
		dispersion(consensus(object), ...)
	}
)

if ( is.null(getGeneric("join")) ) setGeneric('join', function(object, ...) standardGeneric('join') )
setMethod('join', 'list',
	function(object, ...){
		# check validity
		lapply( seq_along(object)
			, function(i){
				if( !any(inherits(object[[i]], c('NMFfit', 'NMFSet'))) )
					stop("invalid class for element ", i, " of input list [expect objects of class NMFfit or NMFSet]")
			}
		)
		
		extra <- list(...)
		if( !is.null(extra$nrun) ){
			if( length(object) != 1 && extra$nrun != 0 ){
				warning("slot nrun force to 0 : setting slot 'nrun' is not allowed if the input list has more than one element")
				extra$nrun <- 0
			}
			extra$nrun <- as.integer(extra$nrun)
		}
		do.call('new', c(list('NMFSet', object), extra))
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
#' @param x 
#' @param ... 
#' @returnType matrix 
#' @return the average connectivity matrix computed on the set of NMF results.
#' @author Renaud Gaujoux \email{renaud@@cbio.uct.ac.za}
#' @export
if ( is.null(getGeneric('consensus')) ) setGeneric('consensus', function(object, ...) standardGeneric('consensus') )
setMethod('consensus', signature(object='NMFSet'), 
	function(object, ...){
		if( length(object) == 0 ) return(NULL)
		# when one stores only the best result, the consensus matrix is computed 
		# along the runs and stored in slot 'consensus' 
		if( length(object) == 1 && length(slot(object, 'consensus')) > 0 ){
			return(slot(object, 'consensus'))
		}
		
		# compute mean connectivity matrix
		c <- sapply(object, function(obj) connectivity(obj, ...)) 
		#-> c should be a matrix whose each column is the connectivity matrices (stored as a vector) for a single run
		if( is.list(c) ) stop('NMFSet::consensus : consensus matrix can only be computed for NMF results of the same dimension.')
		# the connectivity matrices should be square
		n <- trunc(sqrt(nrow(c)))
		if( n^2 != nrow(c) ) stop('NMFSet::consensus : connectivity matrices are not squared')
		con <- matrix( apply(c, 1, mean), n, n )
		
		# name the rows and columns appropriately: use the sample names of the first fit
		rownames(con) <- colnames(con) <- sampleNames(object[[1]])
		
		# return result
		con
	}
)

#' Computes the average purity of a set of NMF runs.
if ( is.null(getGeneric('purity')) ) setGeneric('purity', function(x, class, ...) standardGeneric('purity') )
setMethod('purity', signature(x='NMFSet', class='ANY'), 
	function(x, class, method=NULL, ...){
		c <- sapply(x, purity, class=class, ...)
		
		# aggregate the results
		aggregate.measure(c, method, decreasing=TRUE)		
	}
)

#' Computes the average entropy of a set of NMF runs.
if ( is.null(getGeneric('entropy')) ) setGeneric('entropy', function(x, class, ...) standardGeneric('entropy') )
setMethod('entropy', signature(x='NMFSet', class='ANY'), 
	function(x, class, method=NULL, ...){
		c <- sapply(x, entropy, class=class, ...)		
		
		# aggregate the results
		aggregate.measure(c, method)
	}
)

#' Computes the average final residuals of a set of NMF runs.
if( !isGeneric('residuals') ) setGeneric('residuals', package='stats')
setMethod('residuals', signature(object='NMFSet'), 
	function(object, method='best', ...){
		e <- sapply(object, residuals)
		
		# aggregate the results
		aggregate.measure(e, method)		
	}
)

#' Utility function to aggregate numerical quality measures from \code{NMFSet} objects.
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

#' Returns the best NMF run amongst the set, i.e. the run that have the lower estimation residuals.
#'
#' @param x a \code{NMFSet} object from which to extract the best result
#' @param criteria a character string giving the criteria to use to rank the results
#' @param ... extra parameters passed to the criteria method 
#' 
setMethod('fit', signature(object='NMFSet'),
	function(object){
		
		# shortcut in the case of a single element
		if( length(object) == 1 ) return(object[[1]])
		# retirieve the estimation residuals for each run
		e <- sapply(object, residuals, method='best')
		
		# return the run with the lower
		object[[ which.min(e) ]]
	}
)

#' Returns the cluster prediction defined by the best fit
setMethod('predict', signature(object='NMFSet'),
	function(object, ...){
		predict(fit(object), ...)
	}
)

#' Summary method for class NMFSet
#' Computes the summary for the best fit
setMethod('summary', signature(object='NMFSet'),
	function(object, ...){
		best.fit <- fit(object)
		c(summary(best.fit, ...), cophenetic=cophcor(object), dispersion=dispersion(object))
	}
)
#' Compare different runs of NMF.
#' 
#' This function compares the factorizations obtained by different runs of NMF. This is typically usefull
#' to evaluate and compare how different algorithms perform on the same data.
#' 
#' @param x a list of \code{NMF} objects to compare.
#' @param order a character string giving the ordering criteria
#' the purity and entropy
#' @param values a boolean specifying if the result \code{data.frame} should contain the values
#' of each criteria for all results (\code{values=TRUE}) or for each criteria a character '*' in the cell
#' coresponding to the best result.  
#' 
#' @return a \code{data.frame} with the different comparison criteriae in rows and the methods in column
#'  
if ( is.null(getGeneric('compare')) ) setGeneric('compare', function(x, ...) standardGeneric('compare') )
setMethod('compare', signature(x='list'),
	function(x, ...){
		# wrap up x into a NMFSet object 
		compare(join(x), ...)
	}
)
setMethod('compare', signature(x='NMFSet'),
	function(x, sort.by='residuals', values=FALSE, ...){
		
		# define the sorting schema for each criteria (TRUE for decreasing, FALSE for increasing)
		sorting.schema <- list(method=FALSE, seed=FALSE, metric=FALSE
							, residuals=FALSE, time=FALSE, purity=TRUE, cophenetic=TRUE
							, entropy=FALSE, sparseness=TRUE, rank=FALSE, rss=FALSE)
				
		# for each result compute the different criteriae
		measure.matrix <- sapply(x, 
				function(object, ...){
					#TODO: allow to choose the criteria for 'best'
					# take the best result of the set if necessary
					if( is(object, 'NMFSet') )
						object <- fit(object)
					
					summary(object, ...)
				}
				, ...
		)	
		measure.matrix <- t(measure.matrix)	
		
		# set up the resulting data.frame		
		methods <- sapply(x, function(object, ...){
					if( is(object, 'NMFSet') ) object <- fit(object)
					m <- algorithm(object)
					s <- seeding(object) 
					svalue <- objective(object)
					svalue <- if( is.function(svalue) ) '<function>' else paste("'", svalue,"'", sep='')
					c(method=m, seed=s, metric=svalue)
				}
		)
		methods <- t(methods)	
		res <- as.data.frame(methods, stringsAsFactors=FALSE)	
		
		# add the measures to the result		
		res <- cbind(res, measure.matrix)
				
		# sort according to the user's preference
		sorting.criteria <- colnames(res)
		# ASSERT: all columns measure must have a defined sorting schema 
		if( !all( no.schema <- is.element(sorting.criteria, names(sorting.schema))) ) 
			stop("ASSERT: missing sorting schema for criteria(e): ", paste(paste("'", sorting.criteria[!no.schema], "'", sep=''), collapse=', '))
		sort.by <- match.arg(sort.by, sorting.criteria)	
		res <- res[order(res[[sort.by]], decreasing=sorting.schema[[sort.by]]),]
		
		#TODO: handle argument 'values'
		
		# return result
		res
	}
)

#' Compare objective value trajectories from different NMF algorithms.
setMethod('errorPlot', signature(x='NMFSet'), 
	function(x, ...){
		
		# retrieve normalized residuals tracks
		max.iter <- 0
		tracks <- lapply( x, 
				function(res){
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
				, main='NMF Residuals plots')
		
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

setMethod('metaHeatmap', signature(object='NMFSet'),
		function(object, ...){
			
			x <- consensus(object)
			do.call('metaHeatmap', c(list(x, type='consensus'), list(...)))
		}
)

#' Computes cophenetic correlation coefficient
if ( !isGeneric('cophcor') ) setGeneric('cophcor', function(object, ...) standardGeneric('cophcor') )
setMethod('cophcor', signature(object='matrix'),
	function(object, linkage='average'){
		
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

setMethod('cophcor', signature(object='NMFSet'),
	function(object, ...){
		# compute the consensus matrix
		C <- consensus(object, 'samples')
		
		return( cophcor(C, ...))
	}
)

setMethod('rss', 'NMFSet', 
	function(object, target, ...){
		rss(fit(object, ...), target)
	}
)

