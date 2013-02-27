# Heatmap functions
# 
# Author: Renaud Gaujoux
###############################################################################

#' @include NMF-class.R
#' @include aheatmap.R
NULL

#' @param object an R object 
#' @param ... other arguments 
#' 
#' @export
#' @inline
#' @rdname NMF-defunct
setGeneric('metaHeatmap', function(object, ...) standardGeneric('metaHeatmap') )
#' Defunct method substituted by \code{\link{aheatmap}}.
setMethod('metaHeatmap', signature(object='matrix'),
		function(object, ...){
			local <- function(object, type=c('plain', 'consensus'), class
				, unit.scaling=c('none', 'row', 'column'), palette="YlOrRd"
				, rev.palette=FALSE, show.prediction=TRUE, ...){
			
			.Defunct('metaHeatmap', 'NMF', "The S4 method 'metaHeatmap,matrix' is defunct, use 'aheatmap' instead.")
			
#				# load libary RColorBrewer
#				library(RColorBrewer)
#				
#				# retreive the graphical parameters and match them to the sub-sequent call to 'heatmap.plus.2'
#				graphical.params <- list(...)
#				names(graphical.params) <- .match.call.args(names(graphical.params), 'heatmap.plus.2', in.fun='metaHeatmap', call='NMF::metaHeatmap')
#				
#				type <- match.arg(type)
#				if( type == 'consensus' ){
#					# set default graphical parameters for type 'consensus'
#					graphical.params <- .set.list.defaults(graphical.params
#							, distfun = function(x){ as.dist(1-x) }
#							, main='Consensus matrix'
#							, symm=TRUE
#							, Rowv=TRUE
#							, revC=TRUE
#					)
#					
#					if( missing(palette) ) palette <- 'RdYlBu'
#					if( missing(rev.palette) ) rev.palette <- TRUE
#					if( missing(unit.scaling) ) unit.scaling <- 'none'
#					show.prediction <- FALSE # not used for consensus matrices
#				}
#				
#				# apply unit scaling if necessary
#				unit.scaling <- match.arg(unit.scaling)
#				if( unit.scaling == 'column' )
#					object <- apply(object, 2, function(x) x/sum(x))
#				else if ( unit.scaling == 'row' )
#					object <- t(apply(object, 1, function(x) x/sum(x)))
#				
#				# check validity of palette
#				col.palette <- brewer.pal(brewer.pal.info[palette,'maxcolors'],palette)
#				if( rev.palette ) col.palette <- rev(col.palette) 
#				
#				# set default graphical parameters (if those are not already set)
#				graphical.params <- .set.list.defaults(graphical.params
#						, cexRow=0.8, cexCol=0.8
#						, hclustfun = function(m) hclust(m,method="average")
#						, dendrogram='none'
#						, col=col.palette
#						, scale='none', trace="none"
#						, keysize=1, margins=c(5,10)
#				)
#				
#				# if a known class is provided, add a side color over the top row
#				if( !missing(class) ){
#					if( !is.factor(class) ) class <- as.factor(class)
#					class.num <- as.numeric(class)
#					legend.pal <- palette(rainbow(max(2,nlevels(class))))[1:nlevels(class)]
#					col.matrix <- matrix(legend.pal[class.num], ncol(object), 1)
#					
#					# show association with metagenes
#					if( show.prediction ){
#						# only if there is less than 9 metagenes 
#						# cf. limitation of brewer color palette
#						if( nrow(object) <= 9 ){
#							prediction <- .predict.nmf(object)
#							prediction.num <- as.numeric(prediction)
#							pal.pred <- brewer.pal(max(3,nrow(object)),'Set2')[1:nrow(object)]
#							col.matrix <- cbind(pal.pred[prediction.num], col.matrix)
#							graphical.params <- .set.list.defaults(graphical.params
#									, RowSideColors=pal.pred
#							)
#						}
#						else warning("NMF::metaHeatmap - cannot not show prediction for more than 9 metagenes.")
#					}
#					# do that otherwise heatmap.plus complains
#					if( ncol(col.matrix) < 2 )
#						col.matrix <- cbind(col.matrix, col.matrix)
#					
#					# add the ColSideColors
#					graphical.params <- .set.list.defaults(graphical.params
#							, ColSideColors=col.matrix
#					)
#				}
#				
#				
#				res.heatmap <- do.call('heatmap.plus.2', c(list(object), graphical.params))
#				
#				if( !missing(class) ){
#					# order properly the legend boxes
#					class.num <- as.numeric(class[res.heatmap$colInd])
#					
#					occ <- NA # will store the current number of occurences
#					class.max.occ <- rep(0, nlevels(class)) # will store the current maximum number of occurences per class
#					class.start <- rep(NA, nlevels(class)) # will store the current start of the longer stretch per class
#					last.l <- ''
#					sapply( seq(length(class.num), 1, -1), 
#							function(i){
#								l <- class.num[i]
#								if(l==last.l){
#									occ <<- occ + 1
#								}else{
#									occ <<- 1
#								}
#								if(occ > class.max.occ[l]){
#									class.max.occ[l] <<- occ
#									class.start[l] <<- i
#								}
#								last.l <<- l
#							}
#					)
#					
#					class.ord <- order(class.start)
#					l.names <- levels(class)[class.ord]
#					l.color <- legend.pal[class.ord]
#					legend('top', title='Classes'
#							, legend=l.names, fill=l.color
#							, horiz=TRUE, bty='n')
#				}
#				
#				# return invisible
#				invisible(res.heatmap)
			}
			local(object, ...)
		}
)
#' Deprecated method that is substituted by \code{\link{coefmap}} and \code{\link{basismap}}.
setMethod('metaHeatmap', signature(object='NMF'),
		function(object, ...){
			local <- function(object, what=c('samples', 'features'), filter=FALSE, ...){
			
				what <- match.arg(what)
				if( what == 'samples' ){
					# send deprecated warning
					.Defunct('coefmap', 'NMF', "Direct use of the S4-Method 'metaHeatmap' for 'NMF' objects is defunct, use 'coefmap' instead.")
					
					# call the new function 'coefmap'
					return( coefmap(object, ...) )			
					
				}else if( what == 'features' ){
					# send deprecated warning
					.Defunct('basismap', 'NMF', "Direct use of the S4-Method 'metaHeatmap' for 'NMF' objects is defunct, use 'basismap' instead.")
					
					# call the new function 'basismap'
					return( basismap(object, subsetRow=filter, ...) )
					
				}
			}
			local(object, ...)
	}
)

# match an annotation track against list of supported tracks
match_named_track <- function(annotation, tracks, msg, optional=FALSE){
	
	idx <- 
	if( is.character(annotation) ){
		i <- match(annotation, tracks, nomatch=if(optional) 0L else NA )
		if( any(!is.na(i)) ){
			if( !optional && any(is.na(i)) ){
				stop(msg, "invalid track(s) [", str_out(annotation[is.na(i)])
						, "]: should be one of ", str_out(tracks))
			}
		}
		i
	}else if( is.list(annotation) ){ 
		sapply(annotation, function(x){
					if( isString(x) ) match(x, tracks, nomatch=if(optional) 0L else NA )
					else NA
				})
	}
	
	if( is.null(idx) ) return()
	ok <- !is.na(idx)
	# result
	# remaining annotations
	ann <- annotation[!ok]
	if( length(ann) == 0L ) ann <- NULL
	# track annotations
	tr <- unlist(annotation[which(ok)])
	idx <- idx[which(ok)] 
	if( is.null(names(annotation)) ) names(tr) <- tr
	else{
		mn <- names(tr) == ''
		names(tr)[mn] <- tr[mn]
	}
	others <- tr[idx==0L]
	#
#	list(ann=ann, tracks=tr[idx>0L], others=if(length(others)) others else NULL)
	list(ann=as.list(ann), tracks=tr)
}

#' Heatmaps of NMF Factors
#'
#' The NMF package ships an advanced heatmap engine implemented by the function
#' \code{\link{aheatmap}}.
#' Some convenience heatmap functions have been implemented for NMF models, 
#' which redefine default values for some of the arguments of \code{\link{aheatmap}}, 
#' hence tuning the output specifically for NMF models.
#' 
#' \strong{IMPORTANT:} although they essentially have the same set of arguments, 
#' their order sometimes differ between them, as well as from \code{\link{aheatmap}}.
#' We therefore strongly recommend to use fully named arguments when calling these functions.
#' 
#' @rdname heatmaps
#' @name heatmap-NMF
#' 
#' @examples
#' 
#' ## More examples are provided in demo `heatmaps`
#' \dontrun{
#' demo(heatmaps)
#' }
#' ##
#' 
#' # random data with underlying NMF model
#' v <- syntheticNMF(20, 3, 10)
#' # estimate a model
#' x <- nmf(v, 3)
#' 
#' @demo Heatmaps of NMF objects
#' 
#' #' # random data with underlying NMF model
#' v <- syntheticNMF(20, 3, 10)
#' # estimate a model
#' x <- nmf(v, 3)
#'  
NULL
 
#' \code{basimap} draws an annotated heatmap of the basis matrix.
#' 
#' @details
#' \code{basimap} default values for the following arguments of \code{\link{aheatmap}}:
#' \itemize{
#' \item the color palette;
#' \item the scaling specification, which by default scales each 
#' row separately so that they sum up to one (\code{scale='r1'});
#' \item the column ordering which is disabled;
#' \item allowing for passing feature extraction methods in argument 
#' \code{subsetRow}, that are passed to \code{\link{extractFeatures}}.
#' See argument description here and therein.
#' \item the addition of a default named annotation track, that shows 
#' the dominant basis component for each row (i.e. each feature).
#' 
#' This track is specified in argument \code{tracks} (see its argument description).
#' By default, a matching column annotation track is also displayed, but may be 
#' disabled using \code{tracks=':basis'}.
#' 
#' \item a suitable title and extra information like the fitting algorithm, 
#' when \code{object} is a fitted NMF model. 
#' }
#' 
#' @param object an object from which is extracted NMF factors or a consensus 
#' matrix
#' @param ... extra arguments passed to \code{\link{aheatmap}}. 
#' 
#' @rdname heatmaps
#' @inline
#' @export
#' 
#' @examples 
#' 
#' # show basis matrix
#' basismap(x)
#' \dontrun{
#' # without the default annotation tracks
#' basismap(x, tracks=NA)
#' }
#' 
#' @demo
#'
#' # highligh row only (using custom colors)
#' basismap(x, tracks=':basis', annColor=list(basis=1:3))
#'   
#' ## character annotation vector: ok if it does not contain 'basis'
#' # annotate first and second row + automatic special track
#' basismap(x, annRow=c('alpha', 'beta'))
#' # no special track here
#' basismap(x, annRow=c('alpha', 'beta', ':basis'), tracks=NA)
#' # with special track `basis`
#' basismap(x, annRow=list(c('alpha', 'beta'), ':basis'), tracks=NA)
#' # highligh columns only (using custom colors)
#' basismap(x, tracks='basis:')
#' 
#' # changing the name of the basis annotation track
#' basismap(x, annRow=list(new_name=':basis'))
#' 
setGeneric('basismap', function(object, ...) standardGeneric('basismap') )
#' Plots a heatmap of the basis matrix of the NMF model \code{object}.
#' This method also works for fitted NMF models (i.e. \code{NMFfit} objects).
#' 
#' @inheritParams aheatmap
#' @param subsetRow Argument that specifies how to filter the rows that
#' will appear in the heatmap.  
#' When \code{FALSE} (default), all rows are used. 
#' Besides the values supported by argument \code{subsetRow} of 
#' \code{\link{aheatmap}}, other possible values are:
#' 
#' \itemize{ 
#' \item \code{TRUE}: only the rows that are basis-specific are used.  
#' The default selection method is from \cite{KimH2007}.
#' This is equivalent to \code{subsetRow='kim'}.
#' 
#' \item a single \code{character} string or numeric value that specifies 
#' the method to use to select the basis-specific rows, that should appear in the
#' heatmap (cf. argument \code{method} for function \code{\link{extractFeatures}}).
#' 
#' Note \code{\link{extractFeatures}} is called with argument \code{nodups=TRUE}, 
#' so that features that are selected for multiple components only appear once.
#' }
#' @param tracks Special additional annotation tracks to highlight associations between 
#' basis components and sample clusters:
#' \describe{
#' \item{basis}{matches each row (resp. column) to the most contributing basis component
#' in \code{basismap} (resp. \code{coefmap}).
#' In \code{basismap} (resp. \code{coefmap}), adding a track \code{':basis'} to 
#' \code{annCol} (resp. \code{annRow}) makes the column (resp. row) corresponding to 
#' the component being also highlited using the mathcing colours.}
#' }
#' @param info if \code{TRUE} then the name of the algorithm that fitted the NMF 
#' model is displayed at the bottom of the plot, if available.
#' Other wise it is passed as is to \code{aheatmap}.
#'  
#' 
setMethod('basismap', signature(object='NMF'),
	function(object, color = 'YlOrRd:50'
			, scale = 'r1' 
			, Rowv=TRUE, Colv=NA, subsetRow=FALSE
			, annRow=NA, annCol=NA, tracks = 'basis'
			, main="Basis components", info = FALSE
			, ...){
		
		# resolve subsetRow if its a single value
		if( is.atomic(subsetRow) && length(subsetRow) == 1 ){
			subsetRow <- 
				if( isFALSE(subsetRow) )
					NULL
				else if( isTRUE(subsetRow) ) # use Kim and Park scoring scheme for filtering 			
					extractFeatures(object, format='combine')
				else if( is.character(subsetRow) || is.numeric(subsetRow) ) # use subsetRow as a filtering method
					extractFeatures(object, method=subsetRow, format='combine')
				else stop("NMF::basismap - invalid single value for argument 'subsetRow' [logical, numeric or character expected]")
		}
		
		# extract the basis vector matrix
		x <- basis(object)
		
		# add side information if requested
		info <- if( isTRUE(info) && isNMFfit(object) ) 
					paste("Method:", algorithm(object))
				else if( isFALSE(info) ) NULL
				else info
		
		# process annotation tracks
		ptracks <- process_tracks(x, tracks, annRow, annCol)
		annRow <- ptracks$row 
		annCol <- ptracks$col
		# set special annotation handler
		specialAnnotation(1L, 'basis', function() predict(object, what='features'))
		specialAnnotation(2L, 'basis', function() as.factor(1:nbasis(object)))
		#
			
		# call aheatmap on matrix
		aheatmap(x, color = color, ...
				, scale = scale, Rowv=Rowv, Colv = Colv, subsetRow = subsetRow
				, annRow = annRow, annCol = annCol
				, main = main, info = info)	
	}
)

# check if an object contains some value
anyValue <- function(x){
	length(x) > 0L && !is_NA(x) 
}

grep_track <- function(x){
	list(
		both = grepl("^[^:].*[^:]$", x) | grepl("^:.*:$", x)
		, row = grepl("^:.*[^:]$", x)
		, col = grepl("^[^:].*:$", x)
	)
}

# process extra annotation tracks
process_tracks <- function(data, tracks, annRow=NA, annCol=NA){
	
	if( anyValue(tracks) ){
		
		# extract choices from caller function
		formal.args <- formals(sys.function(sys.parent()))
		choices <- eval(formal.args[[deparse(substitute(tracks))]])
		if( isTRUE(tracks) ) tracks <- choices
		else{
			if( !is.character(tracks) ) 
				stop("Special annotation tracks must be specified either as NA, TRUE or a character vector [", class(tracks), "].")
			
			# check validity
			pattern <- "^(:)?([^:]*)(:)?$"
			basech <- str_match(choices, pattern)
			basetr <- str_match(tracks, pattern)
			tr <- basetr[, 3L]
	#		print(basetr)
	#		print(basech)
			# extend base track name
			i <- charmatch(tr, basech[,3L])
			tr[!is.na(i)] <- basech[i[!is.na(i)],3L]
			tracks_long <- str_c(basetr[,2L], tr, basetr[,4L])
			# extend choices
			tty_choice <- grep_track(choices)
			if( any(tty_choice$both) )
				choices <- c(choices, str_c(':', choices[tty_choice$both]), str_c(choices[tty_choice$both], ':'))
			# look for exact match
			itr <- charmatch(tracks_long, choices)
			if( length(err <- which(is.na(itr))) ){
				stop("Invalid special annotation track name [", str_out(tracks[err], Inf)
					,"]. Should partially match one of ", str_out(choices, Inf), '.')
			}
			tracks[!is.na(itr)] <- choices[itr]
		}
#		print(tracks)
	}
	#
	tty <- grep_track(tracks)
	# create result object
	build <- function(x, ann, data, margin){
		t <- 
		if( anyValue(x) ) as.list(setNames(str_c(':', sub("(^:)|(:$)","",x)), names(x)))
		else NA
		# build annotations
		atrack(ann, t, .DATA=amargin(data,margin))
	}
	
	res <- list()
	res$row <- build(tracks[tty$both | tty$row], annRow, data, 1L)
	res$col <- build(tracks[tty$both | tty$col], annCol, data, 2L)
	#str(res)
	res
}

#' \code{coefmap} draws an annotated heatmap of the coefficient matrix.
#' 
#' @details
#' \code{coefmap} redefines default values for the following arguments of
#' \code{\link{aheatmap}}:
#' \itemize{
#' \item the color palette;
#' \item the scaling specification, which by default scales each 
#' column separately so that they sum up to one (\code{scale='c1'});
#' \item the row ordering which is disabled;
#' \item the addition of a default annotation track, that shows the most
#' contributing basis component for each column (i.e. each sample).
#' 
#' This track is specified in argument \code{tracks} (see its argument description).
#' By default, a matching row annotation track is also displayed, but can be disabled 
#' using \code{tracks='basis:'}.
#' \item a suitable title and extra information like the fitting algorithm, 
#' when \code{object} is a fitted NMF model. 
#' }  
#' 
#' @rdname heatmaps
#' @inline
#' @export
#' 
#' @examples 
#' 
#' # coefficient matrix
#' coefmap(x)
#' \dontrun{
#' # without the default annotation tracks
#' coefmap(x, tracks=NA)
#' }
#' 
#' @demo
#' 
#' # coefficient matrix
#' coefmap(x, annCol=c('alpha', 'beta')) # annotate first and second sample
#' coefmap(x, annCol=list('basis', Greek=c('alpha', 'beta'))) # annotate first and second sample + basis annotation
#' coefmap(x, annCol=c(new_name='basis'))
#' 
setGeneric('coefmap', function(object, ...) standardGeneric('coefmap') )
#' The default method for NMF objects has special default values for 
#' some arguments of \code{\link{aheatmap}} (see argument description).
setMethod('coefmap', signature(object='NMF'),
		function(object, color = 'YlOrRd:50'
				, scale = 'c1'
				, Rowv = NA, Colv = TRUE
				, annRow = NA, annCol = NA, tracks='basis'
				, main="Mixture coefficients", info = FALSE
				, ...){
						
			# use the mixture coefficient matrix
			x <- coef(object)
			
			# add side information if requested
			info <- if( isTRUE(info) && isNMFfit(object) ) 
						paste("Method: ", algorithm(object))
					else if( isFALSE(info) ) NULL
					else info
			
			# process annotation tracks
			ptracks <- process_tracks(x, tracks, annRow, annCol)
			annRow <- ptracks$row
			annCol <- ptracks$col
			# set special annotation handler
			specialAnnotation(1L, 'basis', function() as.factor(1:nbasis(object)))
			specialAnnotation(2L, 'basis', function() predict(object))
			#
			
			## process ordering
			if( isString(Colv) ){
				if( Colv == 'basis' ) Colv <- 'samples'
				if( Colv == 'samples' )
					Colv <- order(as.numeric(predict(object, Colv)))
			}
			##
			
			# call aheatmap on matrix
			aheatmap(x, ..., color = color
					, scale = scale, Rowv = Rowv, Colv=Colv
					, annRow=annRow, annCol = annCol
					, main=main, info = info)
		}
)

