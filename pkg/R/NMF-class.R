#library(R.utils)

# declare old S3 class 'proc_time' to use it as a slot for class NMF 
setOldClass('proc_time', prototype=numeric())

################################
# Class: NMF
################################

#' Base virutal interface to store \strong{Non-negative Matrix Factorization} models.
setClass('NMF'
		, contains = 'VIRTUAL')

#' Compute the fitted target matrix based on a \code{NMF} model.
#' 
#' The estimation depends on the underlying NMF model. For example in the standard model, the target
#' matrix is estimated by the matrix product \eqn{WH}. In other models, the estimate may depend on extra 
#' parameters/matrix (cf. Non-smooth NMF in \code{\link{NMFns-class}}).
#' 
#' @param x an object that inherits from class \code{NMF} 
#' @param ... extra parameters
#' @returnType matrix 
#' @return the fitted target matrix from object \code{x} 
#' @author Renaud Gaujoux \email{renaud@@cbio.uct.ac.za}
#' @export
if ( !isGeneric('fitted') ) setGeneric('fitted', package='stats')
setMethod('fitted', signature(object='NMF'),
		function(object){
			stop("NMF::fitted is a pure virtual method of interface 'NMF'. It should be overloaded in class '", class(object),"'.")
		}
)

#' Get/Set the first factor matrix
if ( !isGeneric("basis")) setGeneric('basis', function(object, ...) standardGeneric('basis') )
setMethod('basis', signature(object='NMF'),
		function(object){
			stop("NMF::basis is a pure virtual method of interface 'NMF'. It should be overloaded in class '", class(object),"'.")
		}
)
if ( !isGeneric("basis<-") ) setGeneric('basis<-', function(object, value) standardGeneric('basis<-') )
setReplaceMethod('basis', signature(object='NMF', value='matrix'), 
		function(object, value){ 
			stop("NMF::basis<- is a pure virtual method of interface 'NMF'. It should be overloaded in class '", class(object),"'.")
		} 
)
#' Get/Set the second factor matrix
if ( !isGeneric("coef")) setGeneric('coef', package='stats')
setMethod('coef', signature(object='NMF'),
		function(object){
			stop("NMF::coef is a pure virtual method of interface 'NMF'. It should be overloaded in class '", class(object),"'.")
		}
)
if ( !isGeneric("coef<-") ) setGeneric('coef<-', function(object, value) standardGeneric('coef<-') )
setReplaceMethod('coef', signature(object='NMF', value='matrix'), 
		function(object, value){ 
			stop("NMF::coef<- is a pure virtual method of interface 'NMF'. It should be overloaded in class '", class(object),"'.")
		} 
)
if ( !isGeneric("coefficients")) setGeneric('coefficients', package='stats')
setMethod('coefficients', signature(object='NMF'),
		function(object){
			coef(object)
		}
)


#' Returns a random sample of the implemented NMF model
if ( !isGeneric('random') ) setGeneric('random', function(x, target, ...) standardGeneric('random') )

#' Returns all the NMF model available
nmf.models <- function(builtin.only=FALSE){
	
	if( builtin.only ) return( .nmf.Models.Builtin )
	
	# return all subclasses of class 'NMF' (minus class 'NMFfit' and its subclasses)
	models <- names(methods::getClass('NMF')@subclasses)
	models.wraps <- c('NMFfit', names(methods::getClass('NMFfit')@subclasses))
	return( models[!is.element(models, models.wraps)] )
	
}

#' Initialization function for NMF models
.nmf.Models.Builtin <- NULL
.init.nmf.models <- function(){	
	.nmf.Models.Builtin <<- nmf.models()
}
################################


#' Factory method for objects that inherit from class \code{NMF}.
#' 
#' @param rank the factorization rank (i.e. the number of metagenes). It is coerced into an integer.
#' @param target an object that defines the dimension of the target matrix to estimate.
#' @param ... extra parameters passed to method \code{new}.
#' @returnType \code{class}
#' @return an instance of class \code{class}. 
#' @seealso NMF-class 
#' @author Renaud Gaujoux \email{renaud@@cbio.uct.ac.za}
#' @export
if ( is.null(getGeneric('newNMF')) ) setGeneric('newNMF', function(rank, target=0, ...) standardGeneric('newNMF'))
#' Main constructor method that is eventually called by all \code{NMF} methods. 
#' 
#' @param rank the factorization rank (i.e. the number of metagenes). It is coerced into an integer.
#' @param target a numeric vector (at most 2-length) giving the dimension of the target matrix to estimate. If \code{target} is of 1-length then a square target matrix is assumed.
#' @param class the class of the object to be created. It must be a valid class name that inherits from class \code{NMF}. Default is \code{'NMF'}.
#' @param ... extra parameters passed to class \code{class} default constructor.
#' @examples 
#' # create a NMF object of factorization rank 3 adapted to fit a 1000 x 30 matrix 
#' obj <- newNMF(3, c(1000,30))
#' obj
#' 
#' @export
setMethod('newNMF', signature(rank='numeric', target='numeric'),
	function(rank, target, model='NMFstd', W, H, ...){
		
		# check validity of the provided class
		if( !isClass(model) ) stop("Invalid model name: class '", model,"' is not defined.")
		if( !extends(model, 'NMF') ) stop("Invalid model name: class '", model,"' does not extend class 'NMF'.")
		if( length(target) == 0 ) stop('Invalid dimensions: target must be at least of length 1')
		if( length(target) > 2 ) stop('Invalid dimensions: target must be at most of length 2')
		
		# compute the target dimension
		target <- as.integer(target)
		n <- target[1]
		m <- if( length(target) == 2 ) target[2] else n 
		# force rank to be an integer
		r <- as.integer(rank)
		
		# build dummy compatible W and H if necessary
		W.was.missing <- FALSE
		if( missing(W) ){
			W <- matrix(NA, n, r)
			W.was.missing <- TRUE
		}
		else{
			if( r == 0 ) r <- ncol(W)
			else if( r < ncol(W) ){
				warning("Objective rank is (",r,") lower than the number of columns in W (",ncol(W),"): only the first of W ", r," columns will be used")
				W <- W[,1:r, drop=FALSE]				
			}
			else if( r > ncol(W) ) stop("Objective rank (",r,") is greater than the number of columns in W (",ncol(W),")")
			
			# resolve consistency with target
			if( n == 0 ) n <- nrow(W)
			else if( n < nrow(W) ){
				warning("Number of rows in target is lower than the number of rows in W (",nrow(W),"): only the first ", n," rows of W will be used")
				W <- W[1:n, , drop=FALSE]				
			}
			else if( n > nrow(W) ) stop("Number of rows in target (",n,") is greater than the number of rows in W (",nrow(W),")")
		}
		
		if( missing(H) ) H <- matrix(NA, ncol(W), m)
		else{
			if( r == 0 ) r <- nrow(H)
			else if( r < nrow(H) ){
				warning("Objective rank (",r,") is lower than the number of rows in H (",nrow(H),"): only the first of H ", r," rows will be used")
				H <- H[1:r,, drop=FALSE]				
			}
			else if( r > nrow(H) ) stop("Objective rank (",r,") is greater than the number of rows in H (",nrow(H),")")
			# force dummy W to be at least compatible with H
			if( W.was.missing ) W <- matrix(NA, n, r)

			# resolve consistency with target
			if( m == 0 ) m <- ncol(H)
			else if( m < ncol(H) ){
				warning("Number of columns in target is lower than the number of columns in H (",ncol(H),"): only the first ", m," columns of H will be used")
				H <- H[, 1:m, drop=FALSE]				
			}
			else if( m > ncol(H) ) stop("Number of columns in target (",m,") is greater than the number of columns in H (",ncol(H),")")
		}
		
		# check validity of matrices W and H (only if one of the target dimension is not null)
		if( n + m > 0 ){
			if( nrow(W) != n ) stop('Invalid number of rows for W: should match number of rows in target [', n, ']')
			if( ncol(W) != r ) stop('Invalid number of columns for W: should match factorization rank [', r, ']')
			if( nrow(H) != r ) stop('Invalid number of rows for H: should match factorization rank [', r, ']')
			if( ncol(H) != m ) stop('Invalid number of columns for H: should match number of columns in target [', m, ']')
		}
		
		# build and return a dummy NMF object
		nmf.debug('newNMF', "Instantiate NMF model:", model)
		new(model, W=W, H=H, ...)		
	}
)
#' Constructor method with \code{rank} as a single argument. 
#' 
#' Because the target dimension it creates a new \em{empty} \code{NMF} object where slots \code{W} 
#' and \code{H} have dimension 0 x \code{rank} and \code{rank} x 0 respectively.
#' 
#' @return an empty \code{NMF} object.
#' 
#' @examples
#' # create an empty NMF object of factorization rank 3
#' obj <- newNMF(3)
#' obj
#' 
#' @seealso is.empty.nmf 
#' @export
setMethod('newNMF', signature(rank='numeric', target='missing'),
		function(rank, target, ...){
			newNMF(rank, 0, ...)
		}
)

#' Creates an empty NMF object
setMethod('newNMF', signature(rank='missing', target='ANY'),
		function(rank, target, ...){
			newNMF(0, target, ...)
		}
)

#' Constructor method to create a new \code{NMF} object to fit a given target matrix.
#' 
#' Only the dimensions of \code{target} are used to construct the \code{NMF} object. The matrix slots 
#' \code{W} and \code{H} are filled  with \code{NA}. 
#' 
#' @return a \code{NMF} object with dimensions compatible with matrix \code{target}.
#' 
#' @examples
#' # create a NMF object of factorization rank 3 adapted to fit a given matrix
#' n <- 1000; m <-30 # the target matrix dimensions 
#' V <- matrix(runif(n*m), n, m) # create a random matrix 
#' obj <- newNMF(3, V)
#' obj
#' all(is.na(obj@@W))
#' 
#' @export
setMethod('newNMF', signature(rank='numeric', target='matrix'),
		function(rank, target, seed=NULL, ...){
			#TODO: what is this seed argument???
			# build an object compatible with the target's dimensions
			newNMF(rank, dim(target), ...)
			
		}	
)

#' Show method for an object of class \code{NMF}.
if ( is.null(getGeneric('show')) ) setGeneric('show', function(object) standardGeneric('show')) 
setMethod('show', signature(object='NMF'), 
		function(object)
		{
			cat("<Object of class:", class(object), ">\n")
			cat("features:", nrow(object), "\n")
			cat("basis/rank:", nbasis(object), "\n")
			cat("samples:", ncol(object), "\n")
		}
)


#' Dims method for an NMF object. It returns the dimension of the target matrix.
#' 
#' @note This method calls internally \code{\link{dims}}
#' 
#' @param x a object that inherits from class \code{NMF}
#' @returnType numeric 
#' @return a 2-length numeric vector giving the dimension of the estimated target matrix. An attribute 
#' \code{'rank'} is set to the factorization rank (i.e. the number of metagenes).
#' @author Renaud Gaujoux \email{renaud@@cbio.uct.ac.za}
#' @seealso dims
#' @export

setMethod('dim', signature(x='NMF'), 
function(x){
	c(nrow(basis(x)), ncol(coef(x)), nbasis(x))	
})

#' Returns the factorization rank of a \code{NMF} object.
#' 
#' For a factorization \deqn{V \equiv WH,} the factorization rank is the number of columns of matrix \eqn{W}.
#' In the context of microarray data, the factorization rank is the number of metagenes. 
#' See \code{\link{NMF-class}} for more details.
#' 
#' @param x an object that inherits from class \code{NMF} 
#' @param ... extra parameters
#' @returnType numeric
#' @return a single numeric value
#' @author Renaud Gaujoux \email{renaud@@cbio.uct.ac.za}
#' @seealso NMF-class, dims
#' @export
setGeneric('nbasis', function(x, ...) standardGeneric('nbasis') )
setMethod('nbasis', signature(x='NMF'), 
	function(x)
	{
		ncol(basis(x))
	}
)


#' Subsetting method for NMF objects
setMethod('[', 'NMF', 
	function (x, i, j, ..., drop = FALSE)
	{
		if( missing(drop) )
			drop <- FALSE
		
		if (missing(i) && missing(j)) {
			# check if there is other arguments
			if (length(list(...)) != 0) 
				stop("NMF::[] method: specify which features or samples to subset")
			# otherwise return the untouched object
			return(x)
		}
		
		# subset the columns of mixture coefficient matrix
		if (!missing(j)) 
			coef(x) <- coef(x)[, j, ..., drop = drop]
		
		# subset the rows of the basis matrix
		if (!missing(i)) 
			basis(x) <- basis(x)[i, , ..., drop = drop]
		
		# return subset object
		return(x)
	}
)

#' Checks whether an \code{NMF} object describes an empty NMF.
#'
#' An empty \code{NMF} object has slots \code{W} and \code{H} with 0 rows and 0 columns respectively.
#' Such objects are usually created to save memory space and often require a specific treatment. 
#' 
#' @note While an empty \code{NMF} estimates a target matrix of dimension \eqn{0 \times 0}, 
#' its factorization rank is usually not null.  
#'  
#' @param object An instance of NMF-class
#' @returnType boolean
#' @return \code{TRUE} if the NMF is empty, \code{FALSE} otherwise
#' 
if ( is.null(getGeneric('is.empty.nmf')) ) setGeneric('is.empty.nmf', function(object) standardGeneric('is.empty.nmf') )
setMethod('is.empty.nmf', signature(object='NMF'), 
		function(object){
			nrow(object) == 0 && ncol(object) == 0
		}
)

#' Initialize a random Nonnegative Matrix Factorization (NMF).
#'
#' The methods sets the slots \code{W} and \code{H} of \code{x} to matrices whose entries are drawn from the uniform distribution.
#'
#' @param x an instance of class \code{NMF} that will be used as template
#' @param target a pair of integer value representing the dimension of the target matrix OR the target matrix itself.
#' When a \code{matrix} is used as the target, the method uses the maximum of its entries as an upperbound.
#' @param r the rank of the factorization (i.e. the number of columns of the first factor)
#' @param ... extra parameters passed to \code{runif}
#' 
#' @return a modified an instance of \code{x}. This implies that the result is an instance of the same class as \code{x}.
#'
#' @seealso runif, NMF-class
#'
setMethod('random', signature(x='NMF', target='numeric'), 
	function(x, target, ...){
					
		# valid parameter 'target'
		if( length(target) != 1 && length(target) != 2 )
			stop('NMF::random : invalid target dimensions [length must be 1 or 2. Here length = ', length(target) ,']')
		if( any(is.na(target)) ) 
			stop('NMF::random : invalid target dimensions [NA values in element(s): ', paste(which(is.na(target)), collapse=' and '), ']')		
		# shortcut for symetric case: provide only one dimension
		if( length(target) == 1 ) target <- c(target, target)
		
		# retrieve dimension of the target matrix
		n <- target[1]; m <- target[2];
		# retrieve the factorization rank					
		r <- nbasis(x)
		
		#Vc# Initialize random matrix: W
		basis(x) <- matrix(runif(n*r, min=0, ...), n, r);
		#Vc# Initialize random matrix: H
		coef(x) <- matrix(runif(r*m, min=0, ...), r, m);
		
		# return the modified object
		x
	}
)

setMethod('random', signature(x='NMF', target='matrix'), 
	function(x, target, max=NULL, ...){	
				
		# compute the upper-bound of the random entries (if not provided)
		if( is.null(max) ){
			no.na <- abs(target[!is.na(target)])
			max <- if( length(no.na) == 0 ) 1 else max(no.na)
		}		
		# build a random NMF with the dimensions of the target matrix upper-bounded by the target's maximum entry.
		random(x, dim(target), max=max, ...)
	}
)

setMethod('random', signature(x='NMF', target='missing'), 
	function(x, target, ...){
		random(x, c(nrow(x),ncol(x)), ...)
	}
)

#' Produces different kind of plots.
#if ( !isGeneric('plot') ) setGeneric('plot', package='graphics')
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
if ( !isGeneric('metaHeatmap') ) setGeneric('metaHeatmap', function(object, ...) standardGeneric('metaHeatmap') )
setMethod('metaHeatmap', signature(object='matrix'),
		function(object, type=c('plain', 'consensus'), class
				, unit.scaling=c('none', 'row', 'column'), palette="YlOrRd"
				, rev.palette=FALSE, show.prediction=TRUE, ...){
			# load libary RColorBrewer
			library(RColorBrewer)
			
			graphical.params <- list(...)
			
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

		graphical.params <- list(...)
		what <- match.arg(what)
		if( what == 'samples' ){
			# use the mixture coefficient matrix
			x <- coef(object)
			
			# impose the order induced by the predicted clusters
			d.fun <- function(){
				con <- connectivity(object)
				function(x){
					as.dist(1- con)
				}
			}
			d.fun <- d.fun()
			# set default graphical parameters for type 'coefficients'
			graphical.params <- .set.list.defaults(graphical.params
					, distfun = d.fun
					#, Colv= order(predict(object))
					, main="Sample view\n[mixture coefficients]"
					, Rowv = FALSE, Colv=TRUE
					# custom parameters
					, type = 'plain'
					, unit.scaling = 'column'
			)
			
		}else if( what == 'features' ){
			if( is.logical(filter) ){ # logical value for filter
				
				if( length(filter) == nrow(object) ) # use only the features specified in filter
					x <- basis(object)[filter,]
				else if( length(filter) == 1 ){ # use Kim and Park scoring scheme for filtering
					x <- if( filter ) basis(extractFeatures(object, format='subset')) else basis(object)
				}
				else stop("NMF::metaHeatmap - invalid logical value for argument 'filter' [must be either single or be of length the number of features in the target matrix]")
			}
			else if( is.numeric(filter) ){ # numerical value for filter
				
				if( length(filter) > 1 )# use only the features specified in filter
					x <- basis(object)[filter,]
				else{ # only keep the spcified number of feature for each basis vector
					ord <- apply( basis(object), 2, 
							function(v, limit){
								order(v, decreasing=TRUE)[1:limit]
							}
					, filter)
					ord <- unique(as.numeric(ord))
					x <- basis(object)[ord,]
				}
			}
			else stop("NMF::metaHeatmap - invalid value for argument 'filter' [logical or numeric expected]")
			
			# set default graphical parameters for type 'basis'
			graphical.params <- .set.list.defaults(graphical.params
					, main="Feature view\n[Basis components]"
					, Rowv = TRUE, Colv=FALSE
					# custom parameters
					, type = 'plain'
					, unit.scaling = 'row'
			)
		}
		else stop('invalid type of heatmap [', what,']')
		
		do.call('metaHeatmap', c(list(x), graphical.params))
	}
)

#if ( !isGeneric('hist') ) setGeneric('hist', package='graphics')
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
#		html.colors = apply( col2rgb( seq(ncol(M))+1 ), 2, function(x) paste('#', paste(intToHex(x), collapse=''), alpha, sep='') )
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

#' Utility function used to sets default elements in a list if they are
#' not already set
#' The default values are given in argument ...
.set.list.defaults <- function(input.list, ...){
	defaults <- list(...)
	lapply(names(defaults), 
		function(name){
			if( is.null(input.list[[name]]) ) input.list[[name]] <<- defaults[[name]]
		}
	)
	input.list
}


#' Computes a set of measures usefull to assess the factorization's quality.
#' 
#' 
#' @param object a \code{NMF} object
#' @return a numeric vector of the measures. 
#' 
if ( !isGeneric('summary') ) setGeneric('summary', package='base')
setMethod('summary', signature(object='NMF'), 
		function(object, class, target){
			
			res <- numeric()
			
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
				res <- c(res, rss=rss(object, target))
			}
			
			# return result
			return(res)
		}
)

#' Computes the sparseness of a vector as defined by Hoyer (2004).
#'
#' This sparseness measure quantifies how much energy of a vector is packed into only few components.
#' It is defined by:
#' \deqn{Sparseness(x) = \frac{\sqrt{n} - \frac{\sum |x_i|}{\sqrt{\sum x_i^2}}}{\sqrt{n}-1}
#', where \eqn{n} is the length of \code{x}.
#'
#' The sparseness is a real number in \eqn{[0,1]}. It is equal to 1 if and only if \code{x} contains 
#' a single nonzero component, and is equal to 0 if and only if all components of \code{x} are equal.
#' It interpolates smoothly between these two extreme values. The closer to 1 is the spraseness the sparser is the vector.
#'
#' @param x a vector.
#' @param ... extra parameters
#' @return the sparseness of \code{x} when \code{x} is a vector
#' , the mean sparseness of the columns of \code{x} when \code{x} is a matrix,
#' , the mean sparseness of the basis (resp. encoding vectors (i.e. rows of the mixture coeffcients matrix))
#' when \code{x} is a \code{NMF} object and \code{what} is \code{'genes'} (resp. \code{what} is \code{'samples'}).
#'
#' @references Hoyer, P. O. (2004)
#' , 'Non-negative Matrix Factorization with Sparseness Constraints'
#' , Journal of Machine Learning Research 5 (2004) 1457–1469
#'
if ( is.null(getGeneric('sparseness')) ) setGeneric('sparseness', function(x, ...) standardGeneric('sparseness') )
setMethod('sparseness', signature(x='numeric'), 
	function(x, ...){
		# get length of x
		n <- length(x)
		# compute and return the sparseness
		( sqrt(n) - sum(abs(x)) / sqrt(sum(x^2)) ) / (sqrt(n)-1)
	}
)
setMethod('sparseness', signature(x='matrix'), 
	function(x, ...){
		# compute the sparseness of each column
		s <- apply(x, 2, sparseness)
		
		# return the mean sparseness
		mean(s)
	}
)
setMethod('sparseness', signature(x='NMF'), 
	function(x, what=c('features', 'samples'), ...){
		# retireve full argument
		what <- match.arg(what)
		
		# return the required sparseness
		sparseness(if(what=='features') basis(x) else coef(x), ...)
	}
)


#' Computes the purity of a clustering given a known factor.
#'
#' The purity is a measure of performance of a clustering method, in recovering the classes defined by a factor a priori known (i.e. one knows the true class labels).
#' Suppose we are given \eqn{l} categories, while the clustering method generates \eqn{k} clusters. Purity is given by:
#' \deqn{Purity = \frac{1}{n} \sum_{q=1}^k \max_{1 \leq j \leq l} n_q^j}
#', where:
#' - \eqn{n} is the total number of samples;
#' - \eqn{n_q^j} is the number of samples in cluster \eqn{q} that belongs to original class \eqn{j} (\eqn{1 \leq j \leq l}).
#'
#' The purity is therefore a real number in \eqn{[0,1]}. The larger the purity, the better the clustering performance.
#'
#' @param x a factor or integer vector that gives the cluster membership for each sample.
#' @param class a factor or integer vector that gives the true class labels for each sample.
#' @return the purity (i.e. a numerical value)
#'
#' @references Kim, H. & Park, H. 
#'	Sparse non-negative matrix factorizations via alternating non-negativity-constrained least squares for microarray data analysis.
#'	Bioinformatics (2007). 
#'	\url{http://dx.doi.org/10.1093/bioinformatics/btm134}.
#'	
if ( is.null(getGeneric("purity")) ) setGeneric('purity', function(x, class, ...) standardGeneric('purity') )
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

#' Computes the entropy of a clustering given a known factor.
#'
#' The entropy is a measure of performance of a clustering method, in recovering classes defined by factor a priori known (i.e. one knows the true class labels).
#' Suppose we are given \eqn{l} categories, while the clustering method generates \eqn{k} clusters. Entropy is given by:
#' \deqn{Entropy = - \frac{1}{n \log_2 l} \sum_{q=1}^k \sum_{j=1}^l n_q^j \log_2 \frac{n_q^j}{n_q}
#', where:
#' - \eqn{n} is the total number of samples;
#' - \eqn{n} is the total number of samples in cluster \eqn{q};
#' - \eqn{n_q^j} is the number of samples in cluster \eqn{q} that belongs to original class \eqn{j} (\eqn{1 \leq j \leq l}).
#'
#' The smaller the entropy, the better the clustering performance.
#'
#' @param x a factor or integer vector that gives the cluster membership for each sample.
#' @param class a factor or integer vector that gives the true class labels for each sample.
#' @return the entropy (i.e. a numerical value)
#'
#' @references Kim, H. & Park, H. 
#'	Sparse non-negative matrix factorizations via alternating non-negativity-constrained least squares for microarray data analysis.
#'	Bioinformatics (2007). 
#'	\url{http://dx.doi.org/10.1093/bioinformatics/btm134}.
#'	
if ( is.null(getGeneric("entropy")) ) setGeneric('entropy', function(x, class, ...) standardGeneric('entropy') )
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

#' Extract the genes that characterize each factor.
#'
#' For each factor the genes are first sorted by decreasing contribution. The first successive ones whose contribution to the factor
#' is greater than their contribution to all other metagenes are selected.
#'
#' @param x the matrix of metagenes. That is a matrix with metagenes in column, genes in row, contain the genes' contribution to each factor
#' @return a list with number of metagenes elements, each being a vector containing the indexes of the characterizing genes
#if ( is.null(getGeneric("computeContrib")) ) setGeneric('computeContrib', function(x, ...) standardGeneric('computeContrib') )
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

#' Computes the genes scores as suggested in Kim, H. & Park, H (Bioinformatics 2007).
#'
#' The score for gene \eqn{i} is defined as:
#' \deqn{S_i = 1 + \frac{1}{\log_2 k} \sum_{q=1}^k p(i,q) \log_2 p(i,q),}
#' where \eqn{p(i,q)} is the probability that the \eqn{i}-th gene contributes to cluster \eqn{q}:
#' \deqn{p(i,q) = \frac{W(i,q)}{\sum_{r=1}^k W(i,r)} }
#'
#' The gene score is a real value within the range [0,1]. The higher the gene score the more cluster-specific the corresponding gene.
#'
#' @references Kim, H. & Park, H. 
#'	Sparse non-negative matrix factorizations via alternating non-negativity-constrained least squares for microarray data analysis.
#'	Bioinformatics (2007). 
#'	\url{http://dx.doi.org/10.1093/bioinformatics/btm134}.
#'	
if ( is.null(getGeneric("featureScore")) ) setGeneric('featureScore', function(object, ...) standardGeneric('featureScore') )
setMethod('featureScore', signature(object='matrix'), 
	function(object, method=c('kim', 'max')){
		
		method <- match.arg(method)
		score <- switch(method,
		
		kim = {
			#for each row compute the score
			s <- apply(object, 1, function(g){
				p_i <- g/sum(g)
				crossprod(p_i, log2(p_i))
			})
			# scale, translate and return the result
			1 + s / log2(ncol(object))		
			}
		, max = {
			apply(object, 1, max)
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


#' Identify the most metagenes-specific genes as suggested in Kim, H. & Park, H (Bioinformatics 2007).
#'
#' The genes are first scored using the function \code{featureScore}. Then only the genes that fullfil both following criteria are retained:
#' - score greater than \eqn{\hat{\mu} + 3 \hat{\sigma}}, where \eqn{\hat{\mu}} and \eqn{\hat{\sigma}} are the median and the median absolute deviation (MAD) of the scores respectively;
#' - the maximum contribution to a metagene is greater than the median of all contributions (i.e. of all elements of W)
#'
#' @references Kim, H. & Park, H. 
#'	Sparse non-negative matrix factorizations via alternating non-negativity-constrained least squares for microarray data analysis.
#'	Bioinformatics (2007). 
#'	\url{http://dx.doi.org/10.1093/bioinformatics/btm134}.
#'	
if ( is.null(getGeneric("extractFeatures")) ) setGeneric('extractFeatures', function(object, ...) standardGeneric('extractFeatures') )
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
				#print(sum(sel))
				
				# build a matrix with:
				#-> row#1=max column index, row#2=max value in row, row#3=row index
				temp <- 0;
				g.mx <- nmfApply(object, 1, 
						function(x){
							temp <<- temp +1
							i <- which.max(x);
							#i <- sample(c(1,2), 1)
							c(i, x[i], temp)
						}
				)
				
				# test the second criteria
				med <- median(basis(object))
				sel2 <- g.mx[2,] >= med
				#print(sum(sel2))
				
				# subset the indices
				g.mx <- g.mx[, sel & sel2] 
				
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

#' 'apply' method for NMF objects
#' When MARGIN = 1, the apply is performed on the lines of the basis matrix
#' When MARGIN = 2, the apply is performed on the columns of the mixture coefficient matrix
if ( is.null(getGeneric("nmfApply")) ) setGeneric('nmfApply', function(object, ...) standardGeneric('nmfApply') )
setMethod('nmfApply', signature(object='NMF'),
	function(object, MARGIN, FUN, ...){
		if( MARGIN == 1 )
			apply(basis(object), 1, FUN, ...)
		else if( MARGIN == 2 )
			apply(coef(object), 2, FUN, ...)
		else stop('NMF::nmfApply : invalid argument MARGIN (expected value is 1 or 2)')
	}
)

#' Utility function to compute the dominant column for each row for a matrix.
.predict.nmf <- function(x, prob=FALSE){
	
	if( !is.matrix(x) ) stop('NMF:::.predict.nmf : only works on matrices')
	if( !prob ){
		#for each column return the (row) index of the maximum
		return( as.factor(apply(x, 2, which.max)) )
	}
	else{
		#for each column return the (row) index of the maximum AND the associated probaility
		res <- apply(x, 2,
				function(p){
					i <- which.max(p)
					c(i, p[i]/sum(p))
				}
		)
		# return the result as a list of two elements
		return( list(predict=as.factor(res[1,]), prob=res[2,]) )
	}
}

#' Compute the dominant column for each row for an NMF object.
#'
#' @param x an NMF object
#' @return a factor of length the number of columns, giving the dominant column for each row
setGeneric('predict', package='stats')
setMethod('predict', signature(object='NMF'),
		function(object, what=c('samples', 'features'), prob=FALSE){
			# determine which matrix to use for the prediction
			what <- match.arg(what)
			x <- if( what=='features' ) t(basis(object)) else coef(object)
			
			# compute the indice of the dominant row for each column
			return( .predict.nmf(x, prob) )
		}
)

#' Compute the dominant column for each row.
#'
#' @param x a matrix containing the mixture coefficients (basis vector in rows, samples in columns)
#' @return a factor of length the number of columns, giving the dominant column for each row
#' @note This function is now deprecated
if ( is.null(getGeneric("clusters")) ) setGeneric('clusters', function(x, ...) standardGeneric('clusters') )
setMethod('clusters', signature(x='matrix'), 
	function(x){
		.Deprecated('.predict.nmf', 'NMF')
		#for each column return the (row) index of the maximum
		as.factor(apply(x, 2, function(p) which.max(p)))
	}
)

#' Compute the dominant metagene for each sample.
#'
#' @param x a NMF object
#' @return a factor of length the number of samples, giving the dominant metagene for each sample
#' @note This function is now deprecated
setMethod('clusters', signature(x='NMF'), 
	function(x, what=c('samples', 'features')){
		.Deprecated('predict', 'NMF')
		what <- match.arg(what)		
		clusters(if( what=='features' ) t(basis(x)) else coef(x))
	}
)

#' Computes the connectivity matrix for the samples based on their metagenes expression.
#'
#' The connectivity matrix of a clustering is a matrix \eqn{C} containing only 0 or 1 entries such that:
#' \deqn{C_{ij} = \left\{\begin{array}{l}1\mbox{ if sample }i\mbox{ belongs to the same cluster as sample }j\\0\mobx{ otherwise}\end{array}\right..}
#'
if ( is.null(getGeneric("connectivity")) ) setGeneric('connectivity', function(x, ...) standardGeneric('connectivity') )
setMethod('connectivity', signature(x='NMF'), 
	function(x, ...){
		c <- predict(x, what='samples');
		outer(c, c, function(x,y) ifelse(x==y, 1,0));
	}
)

if( !isGeneric('rss') ) setGeneric('rss', function(object, ...) standardGeneric('rss'))
setMethod('rss', 'NMF', 
	function(object, target){
		
		# make sure the target is provided
		if( missing(target) ) stop("NMF::rss - Argument 'target' is missing and required to compute the residual sum of squares.")
		
		# use the expression matrix if necessary
		if( inherits(target, 'ExpressionSet') ){
			# requires Biobase
			if( !suppressWarnings(require(Biobase, quietly=TRUE)) ) 
				stop("NMF::rss - The 'Biobase' package is required to extract expression data from 'ExpressionSet' objects [see ?'nmf-bioc']")
			
			target <- Biobase::exprs(target)
		}
		
		return( sum( (fitted(object) - target)^2 ) )
	}
)

#' Common interface to compute matrix distances registered in the NMF registry.
#' 
#' This method provides a single interface to compute the error between a target matrix and its estimate. The 
#' distance method can either be a registered NMF-distance (cf TODO: registry distance) or any defined \code{function} (cf. TODO: details for distance function). 
#' 
#' @author Renaud Gaujoux \email{renaud@@cbio.uct.ac.za}
#' @export
if ( is.null(getGeneric("distance")) ) setGeneric('distance', function(target, x,...) standardGeneric('distance') )
#' Computes the distance between a target matrix and its estimate given by a \code{NMF} object.
#' 
#' #' The target matrix \code{target} and its NMF estimate \code{x} must have the same dimension.
#' 
#' @param target the target \code{matrix} with the same dimension as \code{x}.  
#' @param x a \code{NMF} object with the same dimension as \code{target}.
#' @param method if missing, slot \code{distance} of \code{x} is used.
#' @param ... extra 
#' @returnType numeric
#' @return the distance (a nonnegative numerical value) between matrices \code{target} and \code{fitted(x)}. 
#' @author Renaud Gaujoux \email{renaud@@cbio.uct.ac.za}
#' @export
#' 

setMethod('distance', signature(target='matrix', x='NMF'), 
	function(target, x, method=c('', 'KL', 'euclidean'), ...){
		
		#message('compute distance')
		# determinate the distance measure to use
		if( is.null(method) ){
			warning('Undefined distance method: distance cannot be computed [returned NA]')
			return(as.numeric(NA))
		}
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
			# apply the function and return the result
			return(fun(target, x, ...))
		}
		
		# compute and return the distance measure		
		switch(method,
				euclidean = {
					estim <- fitted(x)
					sum((target - estim)^2)/2
				},
				KL = {
					estim <- fitted(x)
					# NB: treat zero entries separately
					sum( ifelse(target==0, estim, target * log(target/estim) - target + estim) );					
				}
		)					
	}
		
)
