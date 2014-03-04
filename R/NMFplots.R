# Plotting functions for NMF objects
# 
# Author: Renaud Gaujoux
# Creation: 16 Aug 2011
###############################################################################

#' @include NMFSet-class.R
NULL

# Scales a matrix so that its columns sum up to one. 
sum2one <- function(x){
	sweep(x, 2L, colSums(x), '/')
}

#' @import grDevices
corplot <- function(x, y, legend=TRUE, confint=TRUE, scales = 'fixed', ..., add=FALSE){
	
	cols <- rainbow(ncol(x))
	
	# set default arguments
	gpar <- .set.list.defaults(list(...)			
			, ylab=quote(substitute(y))
			, xlab=quote(substitute(x))
			, main="Correlation plot"
			, type='p'
			, pch=19
			, cex=0.8
			, col=alphacol(cols, alpha=90))
			
	if( is.null(colnames(x)) )
		colnames(x) <- paste("column", 1:ncol(x), sep='_')

	# draw plot using matplot
	pfun <- if( add ) matpoints else matplot
	#do.call(pfun, c(list(x, y), gpar))
	# add perfect match line
	#abline(a=0, b=1)	
	
	# initialise result
	res <- list(global=list())
	gco <- lm(as.numeric(y) ~ as.numeric(x))
	res$global$lm <- gco
	grsq <- CI.Rsqlm(gco)
	res$global$cortest <- cor.test( as.numeric(x), as.numeric(y) )
	grsq$rho <- res$global$cortest$estimate
	grsq$alpha <- res$global$lm$coef[2L]
	
	# add legend if requested
    x <- provideDimnames(x, base = list(as.character(1:max(dim(x)))))
    y <- provideDimnames(y, base = list(as.character(1:max(dim(y)))))
    ct.labs <- colnames(x) 
	if( legend ){
		# separate correlations
		res$local <- list(lm=list(), cortest=list())
		lco <- t(sapply(1:ncol(x), function(i){
				co <- lm(y[,i] ~ x[,i])
				res$local$lm[[i]] <<- co
				cotest <- cor.test( as.numeric(x[, i]), as.numeric(y[, i]) )
				res$local$cortest[[i]] <<- cotest
				rsq <- CI.Rsqlm(co)
				return(round(c(Rsq=rsq$Rsq
							, confint=rsq$UCL - rsq$Rsq
							, rho=cotest$estimate
							, alpha=co$coef[2L]), 2))
#				z <- as.numeric(cor.test(x[,i], y[,i])[c('estimate', 'p.value')])
#				z[1] <- round.pretty(z[1], 2)
#				z[2] <- round.pretty(z[2], 3)
#				z
			}
			))
		#
        ct.labs <- sapply(seq_along(ct.labs), function(i){
                    ci <- if( confint ) str_c(' +/- ', lco[i,2]) else ''
                    bquote(.(sprintf('%s (', colnames(y)[i])) 
                            ~ alpha == .(sprintf(' %0.2f | ', lco[i,4])) 
                            ~ rho == .(sprintf(' %.02f | ', lco[i,3]))
                             ~ R^2 == .(sprintf(' %0.2f %s)', lco[i,1], ci)))
                })
	}
    
    df <- data.frame(x = melt(x), y = melt(y))
    df[[5L]] <- factor(df[[5L]], levels = colnames(y))
    ct <- colnames(df)[5L]
    ct.title <- gsub('y.', '', ct, fixed = TRUE)
    p <- ggplot(df, aes_string(x='x.value', y='y.value'
                , color = ct)) + 
            geom_point() +
            xlab(gpar$xlab) + ylab(gpar$ylab) +
            scale_color_discrete(labels = ct.labs) + 
            stat_smooth(method = lm) +
            geom_abline(slope = 1, linetype = 3) +
            facet_grid(paste0('~ ', ct), scales = scales) + 
            labs(color = ct.title)
    if( legend ){
        p <- p + theme(legend.position = 'bottom') + 
                guides(color = guide_legend(ncol = 1))
    }else{
        p <- p + theme(legend.position = 'none')
    }
            
    p$correlations <- res
    p
}

#setMethod('corplot', signature(x='NMFfitXn', y='NMF')
#		, function(x, y, pch=19, ...){
#			
#			i <- 1
#			i0 <- which.best(x)
#			i2 <- which.best(x, maxAD, y)
#			.local <- function(f, skip, order, ...){
#				
#				# reorder if necessary
#				if( !missing(order) && !is.null(order) )
#					f <- match.nmf(f, order)
#				
#				# skip if needed
#				if( i == skip )
#					return()
#				
#				# compute correlations between profiles
#				co <- diag(cor(t(scoef(f)), t(scoef(y))))
#				if( i == 1 ){
#					mp <- plot(co, ylim=c(-1,1), xaxt='n', ...)
#					mtext(side = 1, basisnames(y), at= 1:nbasis(y), line = 1)				
#				}
#				else
#					lines(co, ...)
#				i <<- i+1
#				
#			}
#			lapply(x, .local, skip=i0, col="#00000010", type='l', ...)
#			.local(x[[i0]], 0, col="red", type='o', pch=19, ...)
#			.local(x[[i2]], 0, col="red", type='o', pch=19, lty='dashed', ...)
#			invisible()
#			
#		}
#)

#' Plotting Expression Profiles
#' 
#' @export
profplot <- function(x, ...){
	UseMethod('profplot')
}

#' 
#' The function \code{profplot} draws plots of the basis profiles, i.e. the rows
#' of the coefficient matrix of NMF models. 
#' A given profile is composed of the contribution of the corresponding 
#' basis to each sample.
#' 
#' When using NMF for clustering in particular, one looks for strong
#' associations between the basis and a priori known groups of samples.
#' Plotting the profiles may highlight such patterns.
#' 
#' The function can also be used to compare the profiles from two NMF models or
#' mixture coefficient matrices. In this case, it draws a scatter plot of the
#' paired profiles.
#' 
#' @param x a matrix or an NMF object from which is extracted the mixture
#' coefficient matrix. It is extracted from the best fit if \code{x} is the
#' results from multiple NMF runs.
#' @param y a matrix or an NMF object from which is extracted the mixture
#' coefficient matrix.
#' It is extracted from the best fit if \code{y} is the results from multiple NMF runs.
#' @param scale specifies how the data should be scaled before plotting.
#' If \code{'none'} or \code{NA}, then no scaling is applied and the "raw" data is plotted.
#' If \code{TRUE} or \code{'max'} then each row of both matrices
#' are normalised with their respective maximum values.
#' If \code{'c1'}, then each column of both matrix is scaled into proportions (i.e. to sum up to one).
#' Default is \code{'none'}.
#' @param match.names a logical that indicates if the profiles in \code{y} 
#' should be subset and/or re-ordered to match the profile names in \code{x} 
#' (i.e. the rownames). This is attempted only when both \code{x} and \code{y}
#' have names.
#' @param legend a logical that specifies whether drawing the legend or not, or
#' coordinates specifications passed to argument \code{x} of
#' \code{\link{legend}}, that specifies the position of the legend.
#' @param confint logical that indicates if confidence intervals for the 
#' R-squared should be shown in legend.
#' @param Colv specifies the way the columns of \code{x} are ordered before
#' plotting. It is used only when \code{y} is missing.  It can be: \itemize{
#' \item a single numeric value, specifying the index of a row of \code{x},
#' that is used to order the columns by \code{x[, order(x[abs(Colv),])]}.
#' Decreasing order is specified with a negative index.  \item an integer
#' vector directly specifying the order itself, in which case the columns are
#' ordered by \code{x[, Colv]} \item a factor used to order the columns by
#' \code{x[, order(Colv)]} and as argument \code{annotation} if this latter is
#' missing or not \code{NA}.  \item any other object with a suitable
#' \code{order} method. The columns are by \code{x[, order(Colv)]} }
#' @param labels a character vector containing labels for each sample (i.e.
#' each column of \code{x}). These are used for labelling the x-axis.
#' @param annotation a factor annotating each sample (i.e. each column of
#' \code{x}). If not missing, a coloured raw is plotted under the x-axis and
#' annotates each sample accordingly. If argument \code{Colv} is a factor, then
#' it is used to annotate the plot, unless \code{annotation=NA}.
#' @param ...  graphical parameters passed to \code{\link{matplot}} or \code{\link{matpoints}}.
#' @param add logical that indicates if the plot should be added as points to a previous plot 
#' 
#' @seealso \code{\link{profcor}}
#' @keywords aplot
#' @rdname profplot
#' @export
#' @S3method profplot default
#' @examples
#' 
#' # create a random target matrix
#' v <- rmatrix(50, 10)
#' 
#' # fit a single NMF model
#' res <- nmf(v, 3)
#' profplot(res)
#' 
#' # ordering according to first profile
#' profplot(res, Colv=1) # increasing
#' profplot(res, Colv=-1) # decreasing
#' 
#' # fit a multi-run NMF model
#' res2 <- nmf(v, 3, nrun=3)
#' profplot(res2)
#' 
#' # draw a profile correlation plot: this show how the basis components are 
#' # returned in an unpredictable order 
#' profplot(res, res2)
#' 
#' # looking at all the correlations allow to order the components in a "common" order
#' profcor(res, res2)
#' 
profplot.default <- function(x, y, scale=c('none', 'max', 'c1'), match.names=TRUE
							, legend=TRUE, confint=TRUE
							, Colv, labels, annotation, ..., add = FALSE){
	
	# initialise result list
	res <- list()
	# get extra graphical parameters
	gpar <- list(...)	
	
	# plot a correlation plot of y is not missing
	if( !missing(y) ){
		xvar <- deparse(substitute(x))
		# extract mixture coefficient from x 
		if( isNMFfit(x) ){
			gpar <- .set.list.defaults(gpar
					, xlab=paste("NMF model", xvar, "- Method:", algorithm(x)))
			x <- fit(x)
		}
		if( is.nmf(x) ){
			gpar <- .set.list.defaults(gpar
					, main="Mixture coefficient profile correlations"
					, xlab=paste("NMF model", xvar))
			x <- coef(x)
			
			if( is.null(rownames(x)) )
				rownames(x) <- paste("basis", 1:nrow(x), sep='_')
		}else if( is(x, 'ExpressionSet') ){
			x <- Biobase::exprs(x)
			gpar <- .set.list.defaults(gpar
					, main="Expression profile correlations"
					, xlab=paste("ExpressionSet", xvar))
		}else{
			gpar <- .set.list.defaults(gpar			
					, xlab=paste("Matrix ", xvar))
		}
		# at this stage x must be a matrix
		if( !is.matrix(x) )
			stop("NMF::profplot - Invalid argument `x`: could not extract mixture coefficient matrix")
		
		# extract mixture coefficient from y 
		yvar <- deparse(substitute(y))
		if( isNMFfit(y) ){
			gpar <- .set.list.defaults(gpar
					, ylab=paste("NMF model", yvar, "- Method:", algorithm(y)))
			y <- fit(y)
		}
		if( is.nmf(y) ){
			gpar <- .set.list.defaults(gpar
					, main="Mixture coefficient profile correlations"
					, ylab=paste("NMF model", yvar))			
			y <- coef(y)
		}else if( is(y, 'ExpressionSet') ){
			y <- Biobase::exprs(y)
			gpar <- .set.list.defaults(gpar
					, main="Expression profile correlations"
					, ylab=paste("ExpressionSet", yvar))
		}else{
			gpar <- .set.list.defaults(gpar			
					, ylab=paste("Matrix ", yvar))
		}
		# at this stage y must be a matrix
		if( !is.matrix(y) )
			stop("NMF::profplot - Invalid argument `y`: could not extract profile matrix")
		
		# match names if requested
		if( match.names && !is.null(rownames(x)) && !is.null(rownames(y)) ){
			# match the row in x to the rows in y 
			y.idx <- match(rownames(x), rownames(y), nomatch=0L)
			x.idx <- which(y.idx!=0L)
			# subset and reorder if possible
			if( length(x.idx) > 0L ){
				res$y.idx <- y.idx[x.idx]
				y <- y[y.idx, , drop = FALSE]
				res$x.idx <- x.idx				
				x <- x[x.idx, , drop = FALSE]
			}
		}
		
		# scale to proportions if requested
        if( missing(scale) ) scale <- NULL
        else if( isTRUE(scale) ) scale <- 'max'
        else if( isFALSE(scale) ) scale <- 'none'
        scale <- match.arg(scale)
        scales <- 'free'
		if( scale == 'max' ){
			gpar <- .set.list.defaults(gpar
					, xlim=c(0,1), ylim=c(0,1))
            # scale x
            iscale <- (xm <- apply(abs(x), 1L, max)) > 0 
		    x[iscale, ] <- sweep(x[iscale, , drop = FALSE], 1L, xm[iscale], '/')
            # scale y
            iscale <- (ym <- apply(abs(y), 1L, max)) > 0
            y[iscale, ] <- sweep(y[iscale, , drop = FALSE], 1L, ym[iscale], '/')
            scales <- 'fixed'
		} else if( scale == 'c1' ){
			gpar <- .set.list.defaults(gpar
					, xlim=c(0,1), ylim=c(0,1))
			x <- sum2one(x)
            y <- sum2one(y)
		}else{
			Mx <- max(x, y); mx <- min(x, y)
			# extend default limits by a 0.25 factor
			Mx <- Mx * 1.25
			mx <- mx * 0.75
			gpar <- .set.list.defaults(gpar
					, xlim=c(mx,Mx), ylim=c(mx,Mx))
		}
			
		
		gpar <- .set.list.defaults(gpar			
				, main="Profile correlations")
		# plot the correlation plot		
		p <- do.call(corplot, c(list(x=t(x), y=t(y), scales = scales, legend=legend, confint=confint, add=add), gpar))
        p <- expand_list(p, list(idx.map = res))
		
		# return result list
		return( p )
	}
		
	# extract mixture coefficient
	xvar <- deparse(substitute(x))
	if( isNMFfit(x) ){
		gpar <- .set.list.defaults(gpar, main=paste("Mixture coefficient profiles\nNMF method:", algorithm(x), "- runs:", nrun(x)))
		x <- fit(x)
	}
	if( is.nmf(x) ){
		gpar <- .set.list.defaults(gpar, main="Mixture coefficient profiles")
		x <- coef(x)
	}else if( is(x, 'ExpressionSet') ){
		x <- Biobase::exprs(x)
		gpar <- .set.list.defaults(gpar, main="Expression profiles")
	}
	
	# at this stage x must be a matrix
	if( !is.matrix(x) )
		stop("NMF::profplot - Invalid argument `x`: could not extract profile matrix")
	
	# scale to proportions if requested
    if( missing(scale) || !isTRUE(scale) ) scale <- FALSE
	if( scale ){
		gpar <- .set.list.defaults(gpar, ylim=c(0,1))
		x <- sum2one(x)
	}
	
	# reorder the samples if requested	
	if( missing(labels) ){
		labels <- 
		if( !is.null(colnames(x)) ) colnames(x)
		else 1:ncol(x)			
	} else if( length(labels) != ncol(x) ){
		labels <- rep(labels, length.out=ncol(x))
#	stop("NMF::profplot - Invalid argument `labels`: length should be equal to the number of columns in ", xvar, " [=", ncol(x),"]")
	}
	
	# check annotation
	if( !missing(annotation) && length(annotation) != ncol(x) )
		stop("NMF::profplot - Invalid argument `annotation`:: length should be equal to the number of columns in ", xvar, " [=", ncol(x),"]")
	
	# reorder the columns if requested
	if( !missing(Colv) && !is_NA(Colv) ){
		
		ord <- if( length(Colv) == 1 ){
			if( !is.numeric(Colv) || abs(Colv) > nrow(x) )
				stop("NMF::profplot - Invalid singel argument `Colv`: should be an integer between -nrow(x) and nrow(", xvar,") (i.e. [[-", nrow(x),",", nrow(x),"]])")			
			order(x[abs(Colv),], decreasing=Colv<0)
		}else{
			if( length(Colv) != ncol(x) )
				stop("NMF::profplot - Invalid length for argument `Colv`: should be of length ncol(", xvar, ") [=", nrow(x),"]")
		
			if( is.integer(Colv) && length(setdiff(Colv, 1:ncol(x)))==0 ) Colv
			else order(Colv)
		}
		
		# use Colv as annotation if not requested otherwise
		if( missing(annotation) && is.factor(Colv) )
			annotation <- Colv

		# reorder all relevant quantities
		x <- x[,ord]
		labels <- labels[ord]		
		if( !missing(annotation) && !is_NA(annotation) )
			annotation <- annotation[ord]
	}
	
	# set default arguments
	cols <- rainbow(nrow(x))
	gpar <- .set.list.defaults(gpar
			, xlab="Samples"
			, ylab="Mixture coefficient value"
			, main="Profile plot"
			, type='o'
			, lty=1
			, pch=19
			, cex=0.8
			, col=cols)
		
	# plot using matplot
	do.call(matplot, c(list(x=t(x)), gpar, xaxt='n'))
		
	# add legend if requested
	if( !isFALSE(legend) ){
		if( isTRUE(legend) )
			legend <- 'topleft'
		
		# use the rownames for the legend
		leg <- rownames(x)
		if( is.null(leg) )
			leg <- paste('basis', 1:nrow(x), sep='_')		
		legend(legend, legend=leg, col=gpar$col, lwd=1, pch=gpar$pch)
	}
	
	# axis ticks
	px <- 1:ncol(x)
	axis(1, at = px, labels = FALSE)
	
	# setup grid-base mixed graphic
	vps <- baseViewports()
	pushViewport(vps$inner, vps$figure, vps$plot)
	# clean up on exit
	on.exit(popViewport(3), add=TRUE)
	
	voffset <- 1
	# add sample annotation
	if( !missing(annotation) && !is_NA(annotation) && is.factor(annotation) ){
		
		grid.rect(x = unit(px, "native"), unit(-voffset, "lines")
			, width = unit(1, 'native'), height = unit(1, "lines")
			, gp = gpar(fill=alphacol(rainbow(nlevels(annotation))[annotation], 50), col = 'gray'))	
		voffset <- voffset+1		
	}
	
	# add labels
	if( !is_NA(labels) ){
		# setup grid-base mixed graphic
		#library(gridBase)
		#vps <- baseViewports()
		#pushViewport(vps$inner, vps$figure, vps$plot)
		
		# add axis
		adj <- if( is.character(labels) && max(nchar(labels)) >= 7 ) list(just='right', rot=45)
				else list(just='center', rot=0)
		grid.text(labels
				, x = unit(px, "native"), y = unit(-voffset,"lines")
				, just = adj$just, rot = adj$rot)
		voffset <- voffset+1
		# clean up on exit
		#popViewport(3)
	}
	
	invisible(nrow(x))
	# add xlab
	#if( nchar(xlab) > 0 )
	#	grid.text(xlab, x = unit(length(px)/2, "native"), y = unit(-voffset,"lines"), just = 'center')
		
}

#setGeneric('profplot', function(x, y, ...) standardGeneric('profplot'))
#setMethod('profplot', signature(x='matrix', y='missing')
#		, function(x, y, ...){
#			
#			gpar <- .set.list.defaults(list(...)
#					, xlim=c(0,1), ylim=c(0,1)
#					, main="Profile plot"
#					, type='b'
#					, pch=19)
#			
#			do.call(matplot, c(gpar, x=t(sum2one(x)), y=t(sum2one(y))))
#		}
#)
#setMethod('profplot', signature(x='matrix', y='matrix')
#		, function(x, y, scale=FALSE, ...){
#						
#			# x is the reference, y the estimation
#			if( scale ){
#				gpar <- .set.list.defaults(list(...)
#						, xlim=c(0,1), ylim=c(0,1)
#						, main="Profile correlation plot")
#				do.call(corplot, c(gpar, x=t(sum2one(x)), y=t(sum2one(y))))
#			}else
#				corplot(t(x), t(y), ...)
#	
#		}
#)
#setMethod('profplot', signature(x='matrix', y='NMF')
#		, function(x, y, ...){
#			profplot(x, coef(y), ...)
#		}
#)
#setMethod('profplot', signature(x='NMF', y='ANY')
#		, function(x, y, ...){
#			profplot(coef(x), y, ...)	
#		}
#)
#setMethod('profplot', signature(x='matrix', y='NMFfit')
#		, function(x, y, ...){
#			
#			if( !missing(y) ){ # x is the reference, y the estimation			
#				
#				# map components to the references
#				title <- paste("Profile correlation plot - Method:", algorithm(y))				
#				gpar <- .set.list.defaults(list(...),
#						list(main=title))				
#				do.call(profplot, c(gpar, x=x, y=fit(y)))
#			}
#		}
#)
#setMethod('profplot', signature(x='matrix', y='NMFfitXn')
#		, function(x, y, ...){
#			profplot(x, minfit(y), ...)
#		}
#)
#setMethod('profplot', signature(x='NMFfitXn', y='ANY')
#		, function(x, y, ...){
#			profplot(minfit(x), y, ...)
#		}
#)

#' Silhouette of NMF Clustering
#' 
#' @param x an NMF object, as returned by \code{\link{nmf}}.
#' @param what defines the type of clustering the computed silhouettes are 
#' meant to assess: \code{'samples'} for the clustering of samples 
#' (i.e. the columns of the target matrix),
#' \code{'features'} for the clustering of features (i.e. the rows of the
#' target matrix), and \code{'chc'} for the consensus clustering of samples as
#' defined by hierarchical clustering dendrogram, \code{'consensus'} for the 
#' consensus clustering of samples, with clustered ordered as in the 
#' \strong{default} hierarchical clustering used by 
#' \code{\link{consensusmap}} when plotting the heatmap of the consensus matrix 
#' (for multi-run NMF fits). 
#' That is \code{dist = 1 - consensus(x)}, average linkage and reordering based
#' on row means.  
#' @param order integer indexing vector that can be used to force the silhouette 
#' order.
#' @param ... extra arguments not used.  
#' 
#' @seealso \code{\link[NMF]{predict}}
#' @S3method silhouette NMF
#' @import cluster
#' @examples 
#' 
#' x <- rmatrix(100, 20, dimnames = list(paste0('a', 1:100), letters[1:20]))
#' # NB: using low value for maxIter for the example purpose only
#' res <- nmf(x, 4, nrun = 5, maxIter = 50)
#' 
#' # sample clustering from best fit
#' plot(silhouette(res))
#' 
#' # from consensus
#' plot(silhouette(res, what = 'consensus'))
#' 
#' # feature clustering
#' plot(silhouette(res, what = 'features')) 
#' 
#' # average silhouette are computed in summary measures
#' summary(res)
#' 
#' # consensus silhouettes are ordered as on default consensusmap heatmap
#' \dontrun{ op <- par(mfrow = c(1,2)) }
#' consensusmap(res)
#' si <- silhouette(res, what = 'consensus')
#' plot(si)
#' \dontrun{ par(op) }
#' 
#' # if the order is based on some custom numeric weights
#' \dontrun{ op <- par(mfrow = c(1,2)) }
#' cm <- consensusmap(res, Rowv = runif(ncol(res)))
#' # NB: use reverse order because silhouettes are plotted top-down
#' si <- silhouette(res, what = 'consensus', order = rev(cm$rowInd))
#' plot(si)
#' \dontrun{ par(op) }
#' 
#' # do the reverse: order the heatmap as a set of silhouettes
#' si <- silhouette(res, what = 'features')
#' \dontrun{ op <- par(mfrow = c(1,2)) } 
#' basismap(res, Rowv = si)
#' plot(si)
#' \dontrun{ par(op) }
#' 
silhouette.NMF <- function(x, what = NULL, order = NULL, ...){
    
    # compute prediction
    p <- predict(x, what = what, dmatrix = TRUE)
    # compute silhouette
    si <- silhouette(as.numeric(p), dmatrix = attr(p, 'dmatrix'))
    attr(si, 'call') <- match.call(call = sys.call(-1))
	if( is_NA(si) ) return(NA)
    # fix rownames if necessary
    if( is.null(rownames(si)) ){
        rownames(si) <- names(p)
		if( is.null(rownames(si)) )
			rownames(si) <- 1:nrow(si)
	}
    
    if( is.null(order) && !is.null(attr(p, 'iOrd')) ){
        # reorder as defined in prediction
        order <- attr(p, 'iOrd')
    }
    
    # order the silhouette
    if( !is.null(order) && !is_NA(order) ){
        si[1:nrow(si), ] <- si[order, , drop = FALSE]
        rownames(si) <- rownames(si)[order]
        attr(si, 'iOrd') <- order
        attr(si, 'Ordered') <- TRUE
    }
    
    si
}

#' @S3method silhouette NMFfitX
silhouette.NMFfitX <- function(x, ...){
    si <- silhouette.NMF(x, ...)
    attr(si, 'call') <- match.call(call = sys.call(-1))
    si
}
