# Feature selection functions 
# 
# Author: Renaud Gaujoux
# Created: Mar 18, 2013
###############################################################################

#' @include NMF-class.R
NULL

#' Feature Selection in NMF Models
#' 
#' The function \code{featureScore} implements different methods to computes 
#' basis-specificity scores for each feature in the data.
#' 
#' One of the properties of Nonnegative Matrix Factorization is that is tend to 
#' produce sparse representation of the observed data, leading to a natural 
#' application to bi-clustering, that characterises groups of samples by 
#' a small number of features.
#' 
#' In NMF models, samples are grouped according to the basis 
#' components that contributes the most to each sample, i.e. the basis 
#' components that have the greatest coefficient in each column of the coefficient 
#' matrix (see \code{\link{predict,NMF-method}}).
#' Each group of samples is then characterised by a set of features selected
#' based on basis-specifity scores that are computed on the basis matrix. 
#' 
#' @section Feature scores: 
#' The function \code{featureScore} can compute basis-specificity scores using 
#' the following methods:
#' 
#' \describe{
#' 
#' \item{\sQuote{kim}}{ Method defined by \cite{KimH2007}.
#' 
#' The score for feature \eqn{i} is defined as: 
#' \deqn{S_i = 1 + \frac{1}{\log_2 k} \sum_{q=1}^k p(i,q) \log_2 p(i,q)}{
#' S_i = 1 + 1/log2(k) sum_q [ p(i,q) log2( p(i,q) ) ] },
#' 
#' where \eqn{p(i,q)} is the probability that the \eqn{i}-th feature contributes 
#' to basis \eqn{q}: \deqn{p(i,q) = \frac{W(i,q)}{\sum_{r=1}^k W(i,r)} }{
#' p(i,q) = W(i,q) / (sum_r W(i,r)) }
#' 
#' The feature scores are real values within the range [0,1].
#' The higher the feature score the more basis-specific the corresponding feature.
#' }
#' 
#' \item{\sQuote{max}}{Method defined by \cite{Carmona-Saez2006}.
#' 
#' The feature scores are defined as the row maximums.   
#' }
#' 
#' }
#' 
#' @param object an object from which scores/features are computed/extracted
#' @param ... extra arguments to allow extension 
#' 
#' @return \code{featureScore} returns a numeric vector of the length the number 
#' of rows in \code{object} (i.e. one score per feature). 
#'  
#' @export
#' @rdname scores
#' @inline
#' 
setGeneric('featureScore', function(object, ...) standardGeneric('featureScore') )
#' Computes feature scores on a given matrix, that contains the basis component in columns. 
setMethod('featureScore', 'matrix', 
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
					apply(object, 1L, function(x) max(abs(x)))
				}
		)
		
		# return the computed score
		return(score)
	}
)
#' Computes feature scores on the basis matrix of an NMF model.
setMethod('featureScore', 'NMF', 
	function(object, ...){
		featureScore(basis(object), ...)
	}
)


#' The function \code{extractFeatures} implements different methods to select the 
#' most basis-specific features of each basis component.
#' 
#' @section Feature selection:
#' The function \code{extractFeatures} can select features using the following 
#' methods:
#' \describe{
#' \item{\sQuote{kim}}{ uses \cite{KimH2007} scoring schema and
#' feature selection method.
#' 
#' The features are first scored using the function
#' \code{featureScore} with method \sQuote{kim}.
#' Then only the features that fulfil both following criteria are retained:
#' 
#' \itemize{
#' \item score greater than \eqn{\hat{\mu} + 3 \hat{\sigma}}, where \eqn{\hat{\mu}}
#' and \eqn{\hat{\sigma}} are the median and the median absolute deviation
#' (MAD) of the scores respectively;
#' 
#' \item the maximum contribution to a basis component is greater than the median
#' of all contributions (i.e. of all elements of W).
#' }
#' 
#' }
#' 
#' \item{\sQuote{max}}{ uses the selection method used in the \code{bioNMF} 
#' software package and described in \cite{Carmona-Saez2006}.
#' 
#' For each basis component, the features are first sorted by decreasing 
#' contribution.
#' Then, one selects only the first consecutive features whose highest 
#' contribution in the basis matrix is effectively on the considered basis. 
#' }
#' 
#' }
#'  
#' @return \code{extractFeatures} returns the selected features as a list of indexes,
#' a single integer vector or an object of the same class as \code{object} 
#' that only contains the selected features. 
#' 
#' @rdname scores
#' @inline
#' @export
#' 
setGeneric('extractFeatures', function(object, ...) standardGeneric('extractFeatures') )

# internal functio to trick extractFeatures when format='subset'
.extractFeaturesObject <- local({
	.object <- NULL
	function(object){
		# first call resets .object
		if( missing(object) ){
			res <- .object
			.object <<- NULL
			res
		}else # set .object for next call
			.object <<- object
	}
})

#' Select features on a given matrix, that contains the basis component in columns.
#' 
#' @param method scoring or selection method.
#' It specifies the name of one of the method described in sections \emph{Feature scores} 
#' and \emph{Feature selection}. 
#' 
#' Additionally for \code{extractFeatures}, it may be an integer vector that 
#' indicates the number of top most contributing features to 
#' extract from each column of \code{object}, when ordered in decreasing order, 
#' or a numeric value between 0 and 1 that indicates the minimum relative basis 
#' contribution above which a feature is selected (i.e. basis contribution threshold).
#' In the case of a single numeric value (integer or percentage), it is used for all columns.
#' 
#' Note that \code{extractFeatures(x, 1)} means relative contribution threshold of
#' 100\%, to select the top contributing features one must explicitly specify 
#' an integer value as in \code{extractFeatures(x, 1L)}.
#' However, if all elements in methods are > 1, they are automatically treated as 
#' if they were integers: \code{extractFeatures(x, 2)} means the top-2 most 
#' contributing features in each component.
#' @param format output format. 
#' The following values are accepted:
#' \describe{
#' \item{\sQuote{list}}{(default) returns a list with one element per column in 
#' \code{object}, each containing the indexes of the selected features, as an 
#' integer vector.
#' If \code{object} has row names, these are used to name each index vector. 
#' Components for which no feature were selected are assigned a \code{NA} value.}
#' 
#' \item{\sQuote{combine}}{ returns all indexes in a single vector.
#' Duplicated indexes are made unique if \code{nodups=TRUE} (default).}
#' 
#' \item{\sQuote{subset}}{ returns an object of the same class as \code{object}, 
#' but subset with the selected indexes, so that it contains data only from 
#' basis-specific features.}
#' }
#' 
#' @param nodups logical that indicates if duplicated indexes, 
#' i.e. features selected on multiple basis components (which should in 
#' theory not happen), should be only appear once in the result.
#' Only used when \code{format='combine'}.
#' 
#' @examples
#' 
#' # random NMF model
#' x <- rnmf(3, 50,20)
#' 
#' # probably no feature is selected
#' extractFeatures(x)
#' # extract top 5 for each basis
#' extractFeatures(x, 5L)
#' # extract features that have a relative basis contribution above a threshold
#' extractFeatures(x, 0.5)
#' # ambiguity?
#' extractFeatures(x, 1) # means relative contribution above 100%
#' extractFeatures(x, 1L) # means top contributing feature in each component
#' 
setMethod('extractFeatures', 'matrix', 
	function(object, method=c('kim', 'max')
			, format=c('list', 'combine', 'subset'), nodups=TRUE){
		
		res <-
				if( is.numeric(method) ){
					# repeat single values
					if( length(method) == 1L ) method <- rep(method, ncol(object))
					
					# float means percentage, integer means count
					# => convert into an integer if values > 1						
					if( all(method > 1L) ) method <- as.integer(method)
					
					if( is.integer(method) ){ # extract top features
						
						# only keep the specified number of feature for each column
						mapply(function(i, l)	head(order(object[,i], decreasing=TRUE), l)
								, seq(ncol(object)), method, SIMPLIFY=FALSE)
						
					}else{ # extract features with contribution > threshold
						
						# compute relative contribution
						so <- sweep(object, 1L, rowSums(object), '/')
						# only keep features above threshold for each column
						mapply(function(i, l)	which(so[,i] >= l)
								, seq(ncol(object)), method, SIMPLIFY=FALSE)
						
					}
				}else{
					method <- match.arg(method)
					switch(method,
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
								g.mx <- apply(object, 1L, 
										function(x){
											temp <<- temp +1
											i <- which.max(abs(x));
											#i <- sample(c(1,2), 1)
											c(i, x[i], temp)
										}
								)
								
								# test the second criteria
								med <- median(abs(object))
								sel2 <- g.mx[2,] >= med
								#print(sum(sel2))
								
								# subset the indices
								g.mx <- g.mx[, sel & sel2, drop=FALSE]				
								# order by decreasing score
								g.mx <- g.mx[,order(s[sel & sel2], decreasing=TRUE)]
								
								# return the indexes of the features that fullfil both criteria
								cl <- factor(g.mx[1,], levels=seq(ncol(object))) 
								res <- split(g.mx[3,], cl)
								
								# add the threshold used
								attr(res, 'threshold') <- th
								
								# return result
								res
								
							},
							max = { # MAX method from bioNMF
								
								# determine the specific genes for each basis vector
								res <- lapply(1:ncol(object), 
										function(i){
											mat <- object
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
				}
		
		#Note: make sure there is an element per basis (possibly NA)
		res <- lapply(res, function(ind){ if(length(ind)==0) ind<-NA; as.integer(ind)} )
		
		# add names if possible
		if( !is.null(rownames(object)) ){
			noNA <- sapply(res, is_NA)
			res[noNA] <- lapply(res[noNA], function(x){
						setNames(x, rownames(object)[x])
					})
		}
		
		# apply the desired output format
		format <- match.arg(format)
		res <- switch(format
				#combine: return all the indices in a single vector
				, combine = { 
					# ensure that there is no names: for unlist no to mess up feature names
					names(res) <- NULL
					ind <- na.omit(unlist(res))
					if( nodups ) unique(ind) 
					else ind
				} 
				#subset: return the object subset with the selected indices
				, subset = {
					ind <- na.omit(unique(unlist(res)))
					sobject <- .extractFeaturesObject()
					{if( is.null(sobject) ) object else sobject}[ind, , drop=FALSE]
				}
				#else: leave as a list
				,{
					# add component names if any
					names(res) <- colnames(object)
					res
				}
		)
		
		# add attribute method to track the method used
		attr(res, 'method') <- method
		# return result
		return( res )
	}
)
#' Select basis-specific features from an NMF model, by applying the method 
#' \code{extractFeatures,matrix} to its basis matrix.
#' 
#'  
setMethod('extractFeatures', 'NMF', 
	function(object, ...){
		# extract features from the basis matrix, but subset the NMF model itself
		.extractFeaturesObject(object)
		extractFeatures(basis(object), ...)
	}
)

unit.test(extractFeatures, {
			
	.check <- function(x){
		msg <- function(...) paste(class(x), ':', ...)
		checkTrue( is.list(extractFeatures(x)), msg("default returns list"))
		checkTrue( is.list(extractFeatures(x, format='list')), msg("format='list' returns list"))
		checkTrue( is.integer(extractFeatures(x, format='combine')), msg("format='combine' returns an integer vector"))
		checkTrue( is(extractFeatures(x, format='subset'), class(x)), msg("format='subset' returns same class as object"))
	}
	
	.check(rmatrix(50, 5))
	.check(rnmf(3, 50, 5))
	
})
