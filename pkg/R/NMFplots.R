# Plotting functions for NMF objects
# 
# Author: Renaud Gaujoux
# Creation: 16 Aug 2011
###############################################################################


# Scales a matrix so that its columns sum up to one. 
sum2one <- function(x){
	sweep(x, 2L, colSums(x), '/')
}

# Add an alpha value to colours.
# 
alphacol <- function(x, alpha=FALSE){	
	apply(as.character(as.hexmode(col2rgb(x))), 2, function(x) paste("#", paste(x, collapse=''), alpha, sep=''))
}

corplot <- function(x, y, legend=TRUE, ...){
	
	library(grDevices)
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
	do.call(matplot, c(list(x, y), gpar))
	# add perfect match line
	abline(a=0, b=1)	
	
	gco <- cor( as.numeric(x), as.numeric(y) )
	
	# add legend if requested
	if( legend ){
		lco <- round(diag(cor(x, y)), 2)
		lpar <- .extract.args(gpar, graphics::legend)
		lpar$lty <- -1		
		lpar$pt.cex <- lpar$cex
		lpar$cex <- 1
		do.call('legend', c(list(x='topleft', legend=paste(colnames(x), ' (', lco, ')', sep='')), lpar))
		legend("bottomright", legend=bquote(r == .(round(gco, 2))) )
	}
	invisible(gco)
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

profplot <- function(x, y, scale=FALSE, legend=TRUE, Colv, labels, annotation, ...){
	
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
					, main="NMF profile correlation plot"
					, xlab=paste("NMF model", xvar))
			x <- coef(x)
			
			if( is.null(rownames(x)) )
				rownames(x) <- paste("basis", 1:nrow(x), sep='_')
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
					, main="NMF profile correlation plot"
					, ylab=paste("NMF model", yvar))			
			y <- coef(y)
		}else{
			gpar <- .set.list.defaults(gpar			
					, ylab=paste("Matrix ", yvar))
		}
		# at this stage y must be a matrix
		if( !is.matrix(y) )
			stop("NMF::profplot - Invalid argument `y`: could not extract mixture coefficient matrix")
		
		# scale to proportions if requested
		if( scale ){
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
				, main="Profile correlation plot")
		# plot the correlation plot
		return( do.call(corplot, c(list(x=t(x), y=t(y), legend=legend), gpar)) )
	}
		
	# extract mixture coefficient
	xvar <- deparse(substitute(x))
	if( isNMFfit(x) ){
		gpar <- .set.list.defaults(gpar, main=paste("NMF profile plot\nMethod:", algorithm(x), "- runs:", nrun(x)))
		x <- fit(x)
	}
	if( is.nmf(x) ){
		gpar <- .set.list.defaults(gpar, main="NMF profile plot")
		x <- coef(x)
	}
	
	# at this stage x must be a matrix
	if( !is.matrix(x) )
		stop("NMF::profplot - Invalid argument `x`: could not extract mixture coefficient matrix")
	
	# scale to proportions if requested
	if( scale ){
		gpar <- .set.list.defaults(gpar, ylim=c(0,1))
		x <- sum2one(x)
	}
	
	# reorder the samples if requested	
	labels <- if( missing(labels) ){
				if( !is.null(colnames(x)) ) colnames(x)
				else 1:ncol(x)			
			}else if( isNA(labels) ) NA
			else if( length(labels) != ncol(x) )
				stop("NMF::profplot - Invalid argument `labels`: length should be equal to the number of columns in ", xvar, " [=", ncol(x),"]")
			else
				labels
	
	# check annotation
	if( !missing(annotation) && length(annotation) != ncol(x) )
		stop("NMF::profplot - Invalid argument `annotation`:: length should be equal to the number of columns in ", xvar, " [=", ncol(x),"]")
	
	# reorder the columns if requested
	if( !missing(Colv) && !isNA(Colv) ){
		
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
		if( !missing(annotation) && !isNA(annotation) )
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
	if( !identical(legend, FALSE) ){
		if( isTRUE(legend) )
			legend <- 'topleft'
		
		# use the rownames for the legend
		leg <- rownames(x)
		if( is.null(leg) )
			leg <- paste('basis', 1:nrow(x), sep='_')		
		legend(legend, legend=leg, col=cols, lwd=1, pch=gpar$pch)
	}
	
	# axis ticks
	px <- 1:ncol(x)
	axis(1, at = px, labels = FALSE)
	
	# setup grid-base mixed graphic
	library(gridBase)
	vps <- baseViewports()
	pushViewport(vps$inner, vps$figure, vps$plot)
	# clean up on exit
	on.exit(popViewport(3), add=TRUE)
	
	voffset <- 1
	# add sample annotation
	if( !missing(annotation) && !isNA(annotation) && is.factor(annotation) ){
		
		grid.rect(x = unit(px, "native"), unit(-voffset, "lines")
			, width = unit(1, 'native'), height = unit(1, "lines")
			, gp = gpar(fill=alphacol(rainbow(nlevels(annotation))[annotation], 50), col = 'gray'))	
		voffset <- voffset+1		
	}
	
	# add labels
	if( !isNA(labels) ){
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