# Functions to define/extract compact colour specifications 
# 
# Author: "Renaud Gaujoux"
# Creation: 19 Sep 2011
###############################################################################

#' @import RColorBrewer
#' @import colorspace
#' @import grDevices
NULL

col2hex <- function (cname) 
{
	colMat <- col2rgb(cname)
	rgb(red = colMat[1, ]/255, green = colMat[2, ]/255, blue = colMat[3,]/255)
}

#cc <- function(x, cval=80, lval=30){
#	
#	sapply(x, function(co){
#		if( is.integer(co) ) col2hex(if( co <= 8 ) co else colors()[co])		
#		else if( is.numeric(co) ) hcl(co, c=cval, l=lval)
#		else if( !grepl("^#") )
#		else co
#	})
#	
#}


#' Flags a Color Palette Specification for Reversion
#' @keywords internal
revPalette <- function(x){
	attr(x, 'revPalette') <- TRUE
	x
}

#' Builds a Color Palette from Compact Color Specification
#' @keywords internal
ccPalette <- function(x, n=NA, verbose=FALSE){
		
	if( length(x)==1 ){
	
		# shortcut for single colors
		if( (is_NA(n) || n==1) && length(x) > 1L && all(grepl("^#", x)) ) return(x)
		
		sp <- ccSpec(x)
		x <- sp$palette
		if( is_NA(n) )
			n <- sp$n
		
		a <- attributes(x)
				
		if( is.integer(x) ) # integer code between 1 and 8: R basic colour
			x <- c("#F1F1F1", col2hex(x))
		else if( is.numeric(x) ) # numeric: hcl colour
			x <- rev(sequential_hcl(2, h = x, l = c(50, 95)))
		else if( is.character(x) ){ # Palette name: 
		
			if( require.quiet('RColorBrewer') && x %in% rownames(brewer.pal.info) ){
				if( verbose ) message("Load and generate ramp from RColorBrewer colour palette '", x, "'")			
				x <- brewer.pal(brewer.pal.info[x, 'maxcolors'], x)
			}else{
				cpal <- c('RdYlBu2', 'rainbow', 'heat', 'topo', 'terrain', 'cm', 'gray', 'grey')
				i <- pmatch(x, cpal)				
				if( is.na(i) && (x %in% colours() || grepl("^#[0-9a-fA-F]+$", x)) ){
					x <- c("#F1F1F1", x)
				}else{
					
					if( is.na(i) ){
												
						stop("Invalid palette name '", x, "': should be an RColorBrewer palette or one of "
							, paste("'", cpal ,"'", sep='', collapse=', ')
                            , ".\n  Available RColorBrewer palettes: ", str_out(sort(rownames(brewer.pal.info)), Inf), '.')
					}
					x <- cpal[i]
					
					# use default value of 10 for n if not specified
					np <- if( is_NA(n) ) 10 else n
					
					x <- switch(x
							, RdYlBu2 = c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF", "#E0F3F8", "#91BFDB", "#4575B4")
							, rainbow = rainbow(np)
							, gray = rev(gray.colors(np))
							, grey = rev(grey.colors(np))
							, heat = heat.colors(np)
							, topo = topo.colors(np)
							, terrain = terrain.colors(np)
							, cm = cm.colors(np)
							, stop("Unknown colour palette name: '", x, "'")
					)
				}
			}
		}		
		else
			stop("Invalid colour palette specification :", x)
		
		attributes(x) <- a
	}
	
	# revert the palette if requested
	if( !is.null(attr(x, 'revPalette')) ){
		x <- rev(x)
		attr(x, 'revPalette') <- NULL
	}
	
	# limit to the requested length
	if( !is_NA(n) )
		x <- x[1:n]
		
	# return converted palette
	x
}

#' Generate Break Intervals from Numeric Variables 
#' 
#' Implementation is borrowed from the R core function \code{\link{cut.default}}.
#' 
#' @keywords internal
ccBreaks <- function(x, breaks, center = NULL){
	
	if (!is.numeric(x)) 
		stop("'x' must be numeric")
	
	if (length(breaks) == 1L) {
        
        if (is.na(breaks) | breaks < 2L) 
    		stop("Invalid number of intervals: should be >= 2")
    	nb <- as.integer(breaks + 1)
    	dx <- diff(rx <- range(x, na.rm = TRUE))
    	if (dx == 0) 
    		dx <- abs(rx[1L])
    	
        if( is.null(center) ){
    		breaks <- seq.int(rx[1L] - dx/1000, rx[2L] + dx/1000, length.out = nb)
        }else{ # center the breaks on the requested value
            M <- max(abs(center - rx[1L]), abs(center - rx[2L]))
    		lb <- seq(center-M-dx/1000, center, length.out = floor(nb/2))
            n <- ceiling(nb/2)
    		rb <- seq(center, center+M+dx/1000, length.out = n+1)
    		breaks <- c(lb, rb[-1])
    	}
	}
	
	if (anyDuplicated(breaks)) 
		stop("'breaks' are not unique")
	
	breaks
}

#' Extract Colour Palette Specification
#' 
#' @param a character string that specify a colour palette.
#' @return a list with elements: palette, n and rev
#' 
#' @keywords internal
ccSpec <- function(x){

	n <- NA
	rv <- FALSE
	
	if( length(x) == 1 ){
		if( is.character(x) ){
			
			# flag as reversed if it starts with a '-'
			if( grepl('^-', x) ){		
				x <- substr(x, 2, nchar(x))
				rv <- TRUE
			}
			
			# extract length for string
			sm <- str_match(x, "([^:]+):([0-9]+).*")[1,]
			if( !is_NA(sm[1]) ){
				n <- as.integer(sm[3])
				x <- sm[2]
			}
			
			# convert to a colour code if possible
			# use maximum colour number of brewer sets 
			if( is_NA(n) && isString(x) && require.quiet('RColorBrewer') && x %in% rownames(brewer.pal.info) ){
				n <- brewer.pal.info[x,'maxcolors']
			}else if( grepl("^[1-8]$", x) ){# integer code between 1 and 8: R basic colour
				x <- palette()[as.integer(x)]
			}else if( grepl("^[0-9.]+$", x) ) # numeric: hcl colour
				x <- as.numeric(x)
		
		}else if( is.numeric(x) ){
			if( x < 0 ){
				x <- -x
				rv<- TRUE
			}
			# single integer specification: use R default colours  
			if( isInteger(x) ){
				if( x <= 8 ) x <- palette()[x]
				else x <- colours()[x]
			}
		}
		
	}
	
	if( rv )
		x <- revPalette(x)
	
	res <- list(palette=x, n=n, rev=rv)
#	print(res)
	res
}

#' Builds a Color Ramp from Compact Color Specification
#' 
#' @keywords internal 
ccRamp <- function(x, n = NA, breaks = NULL, data = NULL, ...){
	
    # generate random color specification if necessary
	if( missing(x) )
		x <- round(runif(1) * 360)
	
    # list specification
    if( is.list(x) && length(x) == 2L && isNumber(x[[2L]]) ){
        n <- x[[2L]]
        x <- x[[1L]]
    }
    
	# extract specifications
	sp <- ccSpec(x)
	x <- sp$palette
	# create a palette from specification x
	x <- ccPalette(x, ...)
    
    if( missing(n) ){
        
        # x is a complete colour scale 
        if( is.numeric(x) ) n <- length(x)
        else if( isInteger(breaks) ) n <- breaks
        else n <- sp$n
        
		if( is_NA(n) ){
            if( is.null(breaks) ) n <- 50
            else if( length(breaks) > 1L ) n <- length(breaks) - 1L
		}        
	}
    
	if( is_NA(n) ) n <- length(x)
    
	# create ramp from palette
    if( is.numeric(x) ){
        if( is.null(names(x)) )
            stop("Invalid colour specification: numeric scales must have names.")
        if( !is.null(breaks) ){
            if( isInteger(breaks) ){
                if( breaks < length(x) )
                    stop("Invalid number of breaks: must at least equal to the colour scale length.")
                if( breaks > length(x) ){
                    # generate breaks over full range
                    breaks <- ccBreaks(c(x, data), breaks)
                    y <- x[order(x)]
                    ib <- cut(y, breaks, include.lowest = TRUE, labels = FALSE)
                    cols <- lapply(seq(length(ib)-1), function(i){
                        n <- max(ib[i+1] - ib[i], 2)
                        col <- colorRampPalette(names(y)[c(i, i+1)])(n)
                        res <- setNames(rep(NA, length(col)), col)
                        res[1L] <- y[i]
                        if( n == 2L ) return(res[1L])
                        res[-1L] <- breaks[seq(ib[i]+1, ib[i+1]-1)]
                        res[-length(res)]
                    })
                    # combine result with breaks
                    x <- unlist(cols)
                }
            }else warning("Discarding non integer value of argument `break`: directly using complete colour scale specification.")
        }
        color <- x
    }else{
        color <- colorRampPalette(x)(n)
    
        if( !is.null(data) && (is.null(breaks) || isNumber(breaks)) ){
    		# if a single real number: center the breaks on this value
    		cbreaks <- if( isReal(breaks) ) breaks else NULL
    		breaks <- ccBreaks(data, length(color), center=cbreaks)
    	}
        
        if( !is.null(breaks) ) color <- setNames(breaks, color)
    }
    
    # return mapping
    color
}
