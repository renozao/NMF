#' @include atracks.R
#' @include grid.R
#' @include colorcode.R
NULL

library(grid)
library(gridBase)

# extends gpar objects
c_gpar <- function(gp, ...){
    x <- list(...)
    do.call(gpar, c(gp, x[!names(x) %in% names(gp)]))
}

lo <- function (rown, coln, nrow, ncol, cellheight, cellwidth
, treeheight_col, treeheight_row, legend, main, sub, info
, annTracks, annotation_legend, cexAnn, layout
, fontsize, fontsize_row, fontsize_col, gp){

    # pre-process layout to determine each component position/presence
    gl <- .aheatmap_layout(layout)
    
    # annotation data
	annotation_colors <- annTracks$colors
	row_annotation <- annTracks$annRow
	annotation <- annTracks$annCol
	
    gp0 <- gp
	coln_height <- unit(0, "bigpts")
	if(!is.null(coln)){
		longest_coln = which.max(nchar(coln))
		coln_height <- unit(10, "bigpts") +  unit(1.1, "grobheight", textGrob(coln[longest_coln], rot = 90, gp = c_gpar(gp, fontsize = fontsize_col)))
	}

	rown_width <- rown_width_min <- unit(10, "bigpts")
	if(!is.null(rown)){
		longest_rown = which.max(nchar(rown))
		rown_width <- rown_width_min + unit(1.2, "grobwidth", textGrob(rown[longest_rown], gp = c_gpar(gp, fontsize = fontsize_row)))
	}
	
	gp = c_gpar(gp, fontsize = fontsize)
	# Legend position
	if( !is_NA(legend) ) legend_width <- draw_legend(legend = legend, dims.only = TRUE)
	else legend_width <- unit(0, "bigpts")
    col_legend_height <- row_legend_width <- unit(0, 'bigpts')
    if( !isTRUE(gl$options$legend$horizontal) ) row_legend_width <- legend_width 
    else col_legend_height <- legend_width
    #
    
	
	.annLegend.dim <- function(annotation, fontsize){
		# Width of the corresponding legend
		longest_ann <- unlist(lapply(annotation, names))
		longest_ann <- longest_ann[which.max(nchar(longest_ann))]
		annot_legend_width = unit(1, "grobwidth", textGrob(longest_ann, gp = gp)) + unit(10, "bigpts")
		
		# width of the legend title
		annot_legend_title <- names(annotation)[which.max(nchar(names(annotation)))]
		annot_legend_title_width = unit(1, "grobwidth", textGrob(annot_legend_title, gp = c_gpar(gp, fontface = "bold")))
		
		# total width 
		max(annot_legend_width, annot_legend_title_width) + unit(5, "bigpts")
	}
	
	# Column annotations
	if( !is_NA(annotation) ){
		# Column annotation height		
		annot_height = unit(ncol(annotation) * (cexAnn[2L] * 8 + 2) + 2, "bigpts")
	}
	else{
		annot_height = unit(0, "bigpts")
	}
	
	# add a viewport for the row annotations
	if ( !is_NA(row_annotation) ) {
		# Row annotation width		
		row_annot_width = unit(ncol(row_annotation) * (cexAnn[1L] * 8 + 2) + 2, "bigpts")
	}
	else {
		row_annot_width = unit(0, "bigpts")
	}
	
	# Width of the annotation legend
	annot_legend_width <- 
		if( annotation_legend && !is_NA(annotation_colors) ){ 
			.annLegend.dim(annotation_colors, fontsize)
		}else unit(0, "bigpts")

	# Tree height
	treeheight_col = unit(treeheight_col, "bigpts") + unit(5, "bigpts")
	treeheight_row = unit(treeheight_row, "bigpts") + unit(5, "bigpts") 
	
	# main title
	main_height <- if(!is.null(main)) unit(1, "grobheight", main) + unit(20, "bigpts") else unit(0, "bigpts")
	# sub title
	sub_height <- if(!is.null(sub)) unit(1, "grobheight", sub) + unit(10, "bigpts")	else unit(0, "bigpts")
	# info panel
	if( !is.null(info) ){
		info_height <- unit(1, "grobheight", info) + unit(20, "bigpts")
		info_width <- unit(1, "grobwidth", info) + unit(10, "bigpts")
	}else{
		info_height <- unit(0, "bigpts")
		info_width <- unit(0, "bigpts")		
	}
	
	# Set cell sizes
	if(is.na(cellwidth)){
		matwidth = unit(1, "npc") - rown_width - row_legend_width - row_annot_width  - treeheight_row - annot_legend_width - gl$padding$h
	}
	else{
		matwidth = unit(cellwidth * ncol, "bigpts")
	}

	if(is.na(cellheight)){
		matheight = unit(1, "npc") - treeheight_col - annot_height - main_height - coln_height - col_legend_height - sub_height - info_height - gl$padding$v
	
		# recompute the cell width depending on the automatic fontsize
		if( is.na(cellwidth) && !is.null(rown) ){
			cellheight <- convertHeight(unit(1, "grobheight", rectGrob(0,0, matwidth, matheight)), "bigpts", valueOnly = T) / nrow
			fontsize_row <- convertUnit(min(unit(fontsize_row, 'points'), unit(0.6*cellheight, 'bigpts')), 'points')
			
			rown_width <- rown_width_min + unit(1.2, "grobwidth", textGrob(rown[longest_rown], gp = c_gpar(gp0, fontsize = fontsize_row)))
			matwidth <- unit(1, "npc") - rown_width - row_legend_width - row_annot_width  - treeheight_row - annot_legend_width - gl$padding$h
		}
	}
	else{
		matheight = unit(cellheight * nrow, "bigpts")
	}	
		
	# HACK: 
	# - use 6 instead of 5 column for the row_annotation
	# - take into account the associated legend's width
	# Produce layout()
    layout_size <- list(
                    list(rtree = treeheight_row, rann = row_annot_width, mat = matwidth, rnam = rown_width, leg = legend_width, aleg = annot_legend_width)
                    , list(main = main_height, ctree = treeheight_col, cann = annot_height, mat = matheight, cnam = coln_height, leg = legend_width
                            , sub = sub_height, info = info_height)
                )
    
#    aheatmap_layout(layout, size = layout_size); stop()
    glayout <- vplayout(NULL, layout = layout, size = layout_size)
    # reoder width/height according to layout
    unique.name <- glayout$name
    # push layout
	lo <- glayout$grid.layout
    # drop grid layout spec from result object
    glayout$grid.layout <- NULL
	hvp <- viewport( name=paste('aheatmap', unique.name, sep='-'), layout = lo)
	pushViewport(hvp)
	
	#grid.show.layout(lo); stop('sas')
	# Get cell dimensions
    cellwidth <- cellheight <- 0 
	if( vplayout('mat') ){
    	cellwidth = convertWidth(unit(1, "npc"), "bigpts", valueOnly = T) / ncol
    	cellheight = convertHeight(unit(1, "npc"), "bigpts", valueOnly = T) / nrow
    	upViewport()
    }
		
	height <- as.numeric(convertHeight(sum(lo$height), "inches"))
	width <- as.numeric(convertWidth(sum(lo$width), "inches"))
	# Return minimal cell dimension in bigpts to decide if borders are drawn
	mindim = min(cellwidth, cellheight) 
	return( list(width=width, height=height, vp=hvp
                , mindim=mindim, cellwidth=cellwidth, cellheight=cellheight
                , layout = glayout) )
}

.grid_dendrogram <- function(hc, horiz = FALSE, flip = FALSE){
	
	# convert into an hclust if necessary
	if( is(hc, 'dendrogram') ){
		hca <- attr(hc, 'hclust')
		hc <- if( !is.null(hca) ) hca else as.hclust(hc)
	}
	
	h = hc$height / max(hc$height) / 1.05
	m = hc$merge
	o = hc$order
	n = length(o)
	
	m[m > 0] = n + m[m > 0] 
	m[m < 0] = abs(m[m < 0])
	
	dist = matrix(0, nrow = 2 * n - 1, ncol = 2, dimnames = list(NULL, c("x", "y"))) 
	dist[1:n, 1] = 1 / n / 2 + (1 / n) * (match(1:n, o) - 1)
	
	for(i in 1:nrow(m)){
		dist[n + i, 1] = (dist[m[i, 1], 1] + dist[m[i, 2], 1]) / 2
		dist[n + i, 2] = h[i]
	}
    
    # flip around y-axis if requested
    if( flip ){
        dist[, 2] <- 1 - dist[, 2]
        h <- 1 - h
    }
	
	draw_connection = function(x1, x2, y1, y2, y){
		grid.lines(x = c(x1, x1), y = c(y1, y))
		grid.lines(x = c(x2, x2), y = c(y2, y))
		grid.lines(x = c(x1, x2), y = c(y, y))
	}
	
	# create a rotating viewport for row dendrogram 
	if( horiz ){
		gr = rectGrob()
        pushViewport(viewport(height = unit(1, "grobwidth", gr), width = unit(1, "grobheight", gr), angle = 90))
		on.exit(upViewport())
	}
	
	for(i in 1:nrow(m)){
		draw_connection(dist[m[i, 1], 1], dist[m[i, 2], 1], dist[m[i, 1], 2], dist[m[i, 2], 2], h[i])
	}		
	
}

.base_dendrogram <- function(hc, horiz = FALSE, flip = FALSE, ...){
    if( flip ) # not supported 
        stop("Could not draw ", if( horiz ) 'row' else 'column', " dendrogram: base function plot.dendrogram cannot draw flipped dendrogram")
    
#	suppressWarnings( opar <- par(plt = gridPLT(), new = TRUE) )
	( opar <- par(plt = gridPLT(), new = TRUE) )
	on.exit(par(opar))
	if( getOption('verbose') ) grid.rect(gp = gpar(col = "blue", lwd = 2))
	if( !is(hc, 'dendrogram') )
		hc <- as.dendrogram(hc)
	res <- plot(hc, horiz = horiz, xaxs="i", yaxs="i", axes=FALSE, leaflab="none", ...)
    
    if( !is.null(cluster_spec <- attr(hc, 'cluster.spec')) ){
        sapply(seq(cluster_spec$k)[cluster_spec$which], function(i){
            cluster_spec$which <- i
            d <- NULL
            if( !is.null(cluster_spec$col) ) cluster_spec$col <- alphacol(cluster_spec$col[i], .2)
            cluster_spec$text <- cluster_spec$text[[i]]
            cluster_spec <- expand_list(cluster_spec, horiz = horiz, density = d, lty = 5, lwd = 1.5, lower_rect = 0)
            .d <- function(...){
                rect.dendrogram(hc, ...)
            }
            do.call(.d, cluster_spec)
        })
    }
    invisible(res)
}


draw_dendrogram = function(hc, horizontal = FALSE, flip = FALSE){
    
    .draw.dendrodram <- if( flip ) .grid_dendrogram else .base_dendrogram 
	# create a margin viewport
	if( horizontal ){
        x <- if( flip ) 0.1 else 0
        vp <- viewport(x=x, y=0, width=0.9, height=1,just=c("left", "bottom"))
	}else{
        y <- if( !flip ) 0.1 else 0
        vp <- viewport(x=0, y=y, width=1, height=0.9,just=c("left", "bottom"))
    }
    pushViewport( vp )
	on.exit(upViewport())
	
    .draw.dendrodram(hc, horiz = horizontal, flip = flip)

}

# draw a matrix first row at bottom, last at top
draw_matrix = function(matrix, border_color, txt = NULL, gp = gpar()){
	n = nrow(matrix)
	m = ncol(matrix)
	x = (1:m)/m - 1/2/m
	y = (1:n)/n - 1/2/n
    
    # substitute NA values with empty strings
    if( !is.null(txt) ) txt[is.na(txt)] <- ''
     
    for(i in 1:m){
		grid.rect(x = x[i], y = y, width = 1/m, height = 1/n, gp = gpar(fill = matrix[,i], col = border_color))
        if( !is.null(txt) ){
            grid.text(label=txt[, i],
                                x=x[i],
                                y=y,
#                                just=just,
#                                hjust=hjust,
#                                vjust=vjust,
                                rot=0,
                                check.overlap= FALSE, #check.overlap,
                                default.units= 'npc', #default.units,
#                                name=name,
                                gp=gp,
#                                draw=draw,
#                                vp=vp
                  )
        }
	}
}

draw_colnames = function(coln, gp = gpar()){
	
	m = length(coln)
	
	# decide on the label orientation
	width <- m * unit(1, "grobwidth", textGrob(coln[i <- which.max(nchar(coln))], gp = gp))
	width <- as.numeric(convertWidth(width, "inches"))
	gwidth <- as.numeric(convertWidth(unit(1, 'npc'), "inches"))
	y <- NULL
	if( gwidth < width ){
		rot <- 270
		vjust <- 0.5
		hjust <- 0
		y <- unit(1, 'npc') - unit(5, 'bigpts')
	}else{
		rot <- 0
		vjust <- 0.5
		hjust <- 0.5
	}
	if( is.null(y) ){
		height <- unit(1, "grobheight", textGrob(coln[i], vjust = vjust, hjust = hjust, rot=rot, gp = gp))
		y <- unit(1, 'npc') - height
	}
	
	x = (1:m)/m - 1/2/m
	grid.text(coln, x = x, y = y, vjust = vjust, hjust = hjust, rot=rot, gp = gp)
}

# draw rownames first row at bottom, last on top
draw_rownames = function(rown, gp = gpar()){
	n = length(rown)
	y = (1:n)/n - 1/2/n
	grid.text(rown, x = unit(5, "bigpts"), y = y, vjust = 0.5, hjust = 0, gp = gp)	
}


draw_legend = function(color, breaks, legend, gp = gpar(), opts = NULL, dims.only = FALSE){
    
    # sizes
    padding <- unit(4, 'bigpts')
    thickness <- unit(10, 'bigpts')
    space <- unit(2, 'bigpts')
    
    if( dims.only ){
    	longest_break = which.max(nchar(as.character(legend)))
    	longest_break = unit(1.1, "grobwidth", textGrob(as.character(legend)[longest_break], gp = gp))
    	# minimum fixed width: plan for 2 decimals and a sign 
    	min_lw = unit(1.1, "grobwidth", textGrob("-00.00", gp = gp))
    	longest_break = min(longest_break, min_lw)
    	#title_length = unit(1.1, "grobwidth", textGrob("Scale", gp = c_gpar(gp0, fontface = "bold")))
    	legend_width <- padding + thickness + space + longest_break * 1
        return(legend_width)
    }
    
    if( !isTRUE(opts$expand) ){
        size <- min(unit(1, "npc"), unit(150, "bigpts"))
        shift <- unit(1, "npc") - size
        if( opts$pos == 'middle' ) shift <- shift * .5
        else if( opts$pos == 'bottom' && !isTRUE(opts$horizontal) ) shift <- unit(0, "npc")
        else if( opts$pos == 'top' && isTRUE(opts$horizontal) ) shift <- unit(0, "npc")
        
        if( !isTRUE(opts$horizontal) ){
            pushViewport(viewport(x = 0, y = shift, just = c(0, 0), height = size))
        }else{
            pushViewport(viewport(y = 0, x = shift, just = c(0, 0), width = size))
        }
        
        on.exit( upViewport() )
    }
    
    # compute raltive position for breaks and "ticks"
	tick_pos = (legend - min(breaks)) / (max(breaks) - min(breaks))
    breaks = (breaks - min(breaks)) / (max(breaks) - min(breaks))
    h <- diff(breaks)
    
    txt_shift <- thickness + space + padding
    
    flip_coord <- function(x, flip, max = unit(1, 'npc')){
        if( flip ) max - x 
        else x
    }
    
    if( !isTRUE(opts$horizontal) ){
    	grid.rect(x = flip_coord(padding, opts$flip$h, unit(1, 'npc') - thickness), y = breaks[-length(breaks)], width = thickness, height = h, hjust = 0, vjust = 0
                , gp = gpar(fill = color, col = "#FFFFFF00"))
        grid.text(legend, x = flip_coord(txt_shift, opts$flip$h, txt_shift), y = tick_pos, hjust = 0, gp = gp)
    }else{
        grid.rect(y = flip_coord(padding, !opts$flip$v), x = breaks[-length(breaks)], height = thickness, width = h, hjust = 0, vjust = 0
                , gp = gpar(fill = color, col = "#FFFFFF00"))
        grid.text(legend, y = flip_coord(txt_shift, !opts$flip$v), x = tick_pos, vjust = 0, gp = gp)
    }
}

convert_annotations = function(annotation, annotation_colors){
	
	#new = annotation
	x <- sapply(seq_along(annotation), function(i){
	#for(i in 1:length(annotation)){
		a = annotation[[i]]
		b <- attr(a, 'color')
		if( is.null(b) )
			b = annotation_colors[[names(annotation)[i]]]
		if(class(a) %in% c("character", "factor")){
			a = as.character(a)
      #print(names(b))
      #print(unique(a))
			if ( FALSE && length(setdiff(names(b), a)) > 0){
				stop(sprintf("Factor levels on variable %s do not match with annotation_colors", names(annotation)[i]))
			}
			#new[, i] = b[a]
			b[match(a, names(b))]
		}
		else{
			a = cut(a, breaks = 100)
			#new[, i] = colorRampPalette(b)(100)[a]
			ccRamp(b, 100)[a]
		}
	})

    # shape as a matrix
	if( !is.matrix(x) ) x <- as.matrix(x) 
                
	colnames(x) <- names(annotation)
	return(x)
	#return(as.matrix(new))
}

draw_annotations = function(converted_annotations, border_color, horizontal=TRUE, cex = 1){
	n = ncol(converted_annotations)
	m = nrow(converted_annotations)
    base_size <- 8
    size <- cex * base_size
    psize <- unit(size, "bigpts")
	if( horizontal ){
		x = (1:m)/m - 1/2/m
		y = cumsum(rep(size + 2, n)) - cex * base_size / 2
		for(i in 1:m){
			grid.rect(x = x[i], unit(y[n:1], "bigpts"), width = 1/m, height = psize, gp = gpar(fill = converted_annotations[i, ], col = border_color))
		}
	}else{
		x = cumsum(rep(size + 2, n)) - cex * base_size / 2
		y = (1:m)/m - 1/2/m
		for (i in 1:m) {
			grid.rect(x = unit(x[1:n], "bigpts"), y=y[i], width = psize, 
					height = 1/m, gp = gpar(fill = converted_annotations[i,]
					, col = border_color))
		}
	}
}

draw_annotation_legend = function(annotation_colors, border_color, gp = gpar()){
	
	y = unit(1, "npc")
	
	text_height = convertHeight(unit(1, "grobheight", textGrob("FGH", gp = gp)), "bigpts")	
	for(i in names(annotation_colors)){
		grid.text(i, x = 0, y = y, vjust = 1, hjust = 0, gp = c_gpar(gp, fontface = "bold"))
		y = y - 1.5 * text_height
		#if(class(annotation[[i]]) %in% c("character", "factor")){
		acol <- annotation_colors[[i]]
		if( attr(acol, 'afactor') ){
			sapply(seq_along(acol), function(j){
				grid.rect(x = unit(0, "npc"), y = y, hjust = 0, vjust = 1, height = text_height, width = text_height, gp = gpar(col = border_color, fill = acol[j]))
				grid.text(names(acol)[j], x = text_height * 1.3, y = y, hjust = 0, vjust = 1, gp = gp)
				y <<- y - 1.5 * text_height
			})
		}
		else{
			yy = y - 4 * text_height + seq(0, 1, 0.01) * 4 * text_height
			h = 4 * text_height * 0.02
			grid.rect(x = unit(0, "npc"), y = yy, hjust = 0, vjust = 1, height = h, width = text_height, gp = gpar(col = "#FFFFFF00", fill = ccRamp(acol, 100)))
			txt = c(tail(names(acol),1), head(names(acol))[1])
			yy = y - c(0, 3) * text_height
			grid.text(txt, x = text_height * 1.3, y = yy, hjust = 0, vjust = 1, gp = gp)
			y = y - 4.5 * text_height
		}
		y = y - 1.5 * text_height
	}
}

#' Annotated Heatmap Layout Preview
#' 
#' Shows a diagram of an annotated heatmap layout for given specification.
#' 
#' @param layout layout specification that indicates the relative position 
#' of the heatmap's components.
#' Two layouts can be defined: one horizontal, which relates to components associated to rows, 
#' and one vertical, which relates to components associated with columns.
#' Each layout is specified as a character strings, composed of characters 
#' that encode the order of each component: dendrogram (d), annotation tracks (a), 
#' data matrix (m), labels (l) and legend (L).
#' See section \emph{Layout syntax} for a complete specification  
#'
#' @param size list defining the size of each component (mainly for internal use).
#' 
#' 
#' @details Layout syntax:
#' 
#' Layouts are specified as character strings that can contain the following characters, 
#' each associated with a given component or behaviour:
#' 
#' \strong{Components}
#' \describe{
#' \item{\sQuote{d}}{ dendrogram component}
#' \item{\sQuote{a}}{ annotation tracks}
#' \item{\sQuote{m}}{ data matrix}
#' \item{\sQuote{l}}{ labels}
#' \item{\sQuote{L}}{ scale legend}
#' }
#' 
#' \strong{Behaviours}
#' \describe{
#' \item{\sQuote{^}}{ align top (resp. left) for horizontal (resp. vertical) layout.}
#' \item{\sQuote{-}}{ align middle (resp. center) for horizontal (resp. vertical) layout.}
#' \item{\sQuote{_}}{ align bottom (resp. right) for horizontal (resp. vertical) layout.
#' If used alone (i.e. \code{layout = "_"}), then this is equivalent to \code{"|.L_"}, 
#' which places the legend horizontally on the bottom-right corner.}  
#' \item{\sQuote{*}}{ used either alone or after after \sQuote{L} to specifiy that the 
#' legend should expand to full height/width.}
#' }
#' The specification must contain one instance of each of these character.
#' 
#' The default horizontal/vertical layout is \code{"daml"}, and can also be specified
#' as \code{"."}.
#' 
#' Separate layouts can be passed as a character vector with 2 element (e.g., \code{c("daml", "mald")}),
#' or as a single string, with layouts separated by \code{"|"} (e.g., \code{"almd | L."}).
#' When using this separator, a layout specification may be omitted, indicating
#' that the default layout shoud be used: \code{"almd|"} is equivalent to \code{"almd | ."}.
#' If only one layout specification is passed (i.e. a string without \code{"|"}), 
#' then it is used for both horizontal and vertical layouts.
#' 
#' \strong{Shortcuts}
#' \itemize{
#' \item \code{layout = "*"} is a shortcut for \code{layout = ".L*"}, which expands the legend to 
#' take up full height;
#' \item \code{layout = "_"} is a shortcut for \code{layout = "|.L_"}, which puts the legend at 
#' bottom-right corner;
#' \item \code{layout = "_^"} is a shortcut for \code{layout = "|.L^"}, which puts the legend at 
#' bottom-left corner;
#' \item \code{layout = "_*"} is a shortcut for \code{layout = "|.L*"}, which puts the legend 
#' at bottom, expanded to take up full width;
#' \item \code{layout = "^"} is a shortcut for \code{layout = "L^.|"}, which puts the legend 
#' on the top-left corner.
#' }
#' 
#' \strong{Examples:} 
#' \itemize{
#' \item \code{layout = "dlma"} puts labels at the leaves of the dendrograms and 
#' annotation track below or at the right of the data matrix
#' \item \code{layout = ". | amld"} use the default layout for rows, put
#' column annotation track on top of the data matrix, followed by column labels and 
#' dendrogram.
#' }
#' 
#' @export
#' @examples 
#' 
#' # default layout
#' aheatmap_layout()
#' 
#' # Common row/column layout: annotations > data > labels > dendrogram 
#' aheatmap_layout('amld')
#' 
#' # Separate row/column layout: row as above / column as default
#' aheatmap_layout('amld | .')
#' 
#' ## Legend
#' # horizontal bottom-right
#' aheatmap_layout('_')
#' # hotizontal top-left (equivalent to "|L^.")
#' aheatmap_layout('^')
#' 
aheatmap_layout <- function(layout = 'daml', size = NULL){
    
    # default sizes
    l <- unit(2, 'line')
    size0 <- list(list(rtree = 2 * l, rann = l, mat = unit(1, 'null'), rnam = l, leg = l, aleg = l)
            , list(main=l, ctree = 2 * l, cann = l, mat = unit(1, 'null'), cnam = l, leg = l, sub = l, info = 1.5 * l))
    
    if( is_NA(size) ){
        if( nargs() == 1L ) return( size0 )
        size <- NULL
    }
    
    # define dummy sizes
    if( is.null(size) ) size <- size0
    
    # compute layout
    gl <- vplayout(NULL, layout = layout, size = size)
    lo <- gl$grid.layout

    # plot layout diagram
    grid.newpage()
    grid.show.layout(lo, unit.col = NA, cell.label = FALSE)
    
    # label components
    pushViewport(viewport(0.5, 0.5, 0.8, 0.8, layout = lo))
    labels <- c(tree = 'dendrogram', ann = 'annotation tracks', nam = 'labels')
    ilabels <- function(x, y) setNames(paste(x, labels), paste0(y, names(labels)))
    labels <- c(ilabels('Row', 'r'), ilabels('Column', 'c'), leg = 'Scale legend', aleg = 'Annotation legend'
            , mat = 'Data'
            , main = 'Main title', sub = 'Subtitle', info = 'Extra info pane')
    
    # taken from internal loop in grid.show.layout
    label.vp <- function(x, rot = 0){
        vplayout(x)
        grid.text(labels[x], rot = rot)
        popViewport()
    }
    
    sapply(gl$h[!gl$h == 'mat'], label.vp, 90)
    sapply(gl$v, label.vp)
    popViewport()

}

.aheatmap_layout <- function(layout = 'daml', size = NULL){
    
    layout <- as.character(layout)
    
    # defaults
    default <- 'daml'
    defaultL <- paste0(default, 'L')
    v_default <- strsplit(default, '')[[1L]]
    cex.pad <- 1
    
    layout <- gsub(' ', '', layout, fixed = TRUE)
    x <- layout
    if( length(x) == 1L ){
        # special legend specification
        if( x == "*" ) x <- paste0(default, "L*")
        else if( x == "^" ) x <- c("L^.", default)
        else if( x == "_" ) x <- c(default, ".L_")
        else if( x == "_*" ) x <- c(default, ".L*")
        else if( x == "_^" ) x <- c(default, ".L^")
        else{
            x <- gsub("(\\|)?L?([-_*^])", "\\1L\\2", x)
            x <- strsplit(x, '|', fixed = TRUE)[[1L]]
            # deal with ending "|"
            if( grepl("\\|$", layout) ) x <- c(x, '')
        }
    }
    # resolve shortcuts
    if( length(x) == 1L ) x <- c(x, x)
    x <- gsub('.', defaultL, x, fixed = TRUE)
    
    # extract padding
    cex.pad <- sapply(x, function(x){
        m <- str_match(x, "^([0-9.]+)[0-9]")
        if( !is.na(m[, 1]) ) m[, 2]
        else cex.pad
    })
    cex.pad <- as.numeric(cex.pad)
    
    # split into letters
    x_v <- x # keep vector version for later use
    x <- strsplit(x, '')
    x <- lapply(x, unique)
    
    # process each layout
    x <- lapply(x, function(x){    
        xp <- if( !length(x) ) v_default # empty => default
        else if( identical(x, '!') ) 'm'
        else if( identical(x, 'L') ) c(v_default, 'L')
        else if( length(i <- grep('/', x, fixed = TRUE)) ){ # detect component skip
            skip <- tail(x, length(x) - i)
            if( any(skip == '*') ) skip <- setdiff(union(skip, v_default), 'm')
            if( i == 1L ) setdiff(v_default, skip)
            else setdiff(head(x, i-1), skip)
        } else x
        # only keep component characters
        xp <- intersect(xp, c(v_default, 'L'))
        if( !length(xp) ) v_default
        else xp
    })
    
#    if( !any(unlist(x) == 'L') ){
#        x[[1L]] <- c(x[[1L]], 'L')
#    }else 
    if( sum(unlist(x) == 'L') > 1L ){ # both
        x[[2L]] <- setdiff(x[[2L]], 'L')
    }
    
    e_order <- list(v=NULL, h=NULL)
    lexique <- c(tree = 'd', ann = 'a', nam = 'l')
    # vertical layout
    elements <- c(setNames(lexique, paste0('c', names(lexique))), mat = 'm', leg = 'L')
    ie <- match(elements, x[[2L]])
    i <- cbind(ie+1, match('m', x[[1L]]))
    rownames(i) <- names(elements)
    res <- i
    e_order$v <- c('main', names(elements)[setdiff(order(ie), which(is.na(ie)))], 'sub', 'info')
    
    # data matrix position
    xm <- res['mat', 1]
    ym <- res['mat', 2]
    
    # horizontal layout
    elements <- c(setNames(lexique, paste0('r', names(lexique))), mat = 'm', leg = 'L')
    ie <- match(elements, x[[1L]])
    i <- cbind(xm, ie)
    rownames(i) <- names(elements)
    res <- rbind(res, i)
    e_order$h <- names(elements)[setdiff(order(ie), which(is.na(ie)))]
    # add annotation legend only if necessary
    if( with_aleg <- any(unlist(x) == 'a') ){
        e_order$h <- c(e_order$h, 'aleg')
    }
    
    # fixed elements
    nc <- length(e_order$h)
    nr <- length(e_order$v)
    res <- rbind(res, main = c(1, ym)
			, sub = c(nr-1, ym)
			, info = c(nr, ym)
            , aleg = if( with_aleg ) c(xm, nc)
    )
    
    colnames(res) <- c('x', 'y')
    res <- res[!is.na(res[, 1]) & !is.na(res[, 2]), ]
    res <- res[!duplicated(rownames(res)), ]
    
    # build result list
    res <- c(list(layout = res), e_order, options = list())
    
    ## Component-specific options
    # dendrogram orientation
    flip <- sapply(e_order, function(x) grep('tree', x) > which(x=='mat'), simplify = FALSE)
    res$options$dendrogram <- list(flip = flip)
    
    # legend 
    hleg <- !'leg' %in% res$h
    res$options$legend <- list(horizontal = hleg)
    # detect legend positionning specifications
    leg_spec <- str_match(x_v, "L((\\^)|(-)|(_)|(\\*))")[1+hleg, -(1:2)]
    leg_spec <- !is.na(leg_spec) & nzchar(leg_spec)
    if( !length(ipos <- which(leg_spec[1:3])) ) ipos <- 2 * hleg + 1
    res$options$legend$pos <- c("top", "middle", "bottom")[ipos]
    res$options$legend$expand <- leg_spec[4L]
    # text orientation
    res$options$legend$flip <- sapply(e_order, function(x) grep('^leg', x) < which(x=='mat'), simplify = FALSE)
    ##
    
    ## grid.layout: generate and order component sizes according to layout specification
    padding <- unit(4, 'bigpts')
    res$padding <- list(h = (length(res$h) + 1) * cex.pad[1] * padding
                        , v = (length(res$v) + 1) * cex.pad[2] * padding)
    if( !is.null(size) ){
        wunits <- size[[1L]][res$h]
        hunits <- size[[2L]][res$v]
        padd <- function(x, size){
            tmp <- list(size)
            sapply(seq_along(x), function(i){ tmp[[2*i]] <<- x[[i]]; tmp[[2*i+1]] <<- tmp[[1L]]})
            tmp
        }
        
        # add padding
        res$cex.pad <- cex.pad	
        if( cex.pad[1] ) wunits <- padd(wunits, cex.pad[1] * padding)
        if( cex.pad[1] ) hunits <- padd(hunits, cex.pad[2] * padding)
        
        lo <- grid.layout(nrow = length(hunits), ncol = length(wunits)
	            , widths = do.call('unit.c', wunits)
	            , heights = do.call('unit.c', hunits))
        
        res$grid.layout <- lo
    }
    #
    
    invisible(res)
}

vplayout <- local(
{
	graphic.name <- NULL
	.index <- 0L
    .layout <- NULL
	function(x, y, verbose = getOption('verbose'), layout = NULL, ...){
		# initialize the graph name
		if( is.null(x) ){
            .index <<- .index + 1L
			graphic.name <<- paste0("AHEATMAP.VP.", .index) #grid:::vpAutoName()
            # determine layout
            if( is.null(layout) || identical(layout, 'default') ){ #default
                layout <- 'damlLA | daml'
            }
            .layout <<- .aheatmap_layout(layout, ...)
			return(c(list(name = graphic.name), .layout))
		}
		name <- NULL
		if( !is.numeric(x) ){
					
			name <- paste(graphic.name, x, sep='-')
			
			if( !missing(y) && is(y, 'viewport') ){
				y$name <- name
				return(pushViewport(y))			
			}
			if( !is.null(tryViewport(name, verbose=verbose)) ) return(TRUE)
			
            # lookup viewport 
            mlayout <- .layout$layout
            if( !x %in% rownames(mlayout) ) return(FALSE)
#                stop("aheatmap - invalid component name [", x, ']')
            xy <- mlayout[x, ]
            x <- xy[1L]
            if( .layout$cex.pad[1] ) x <- x * 2
            y <- xy[2L]
            if( .layout$cex.pad[2] ) y <- y * 2
		}
		if( verbose ) message("vp - create ", name)
		pushViewport(viewport(layout.pos.row = x, layout.pos.col = y, name=name))
        TRUE
	}	
})

#' Open a File Graphic Device
#'
#'  Opens a graphic device depending on the file extension
#' 
#' @keywords internal
gfile <- function(filename, width, height, ...){ 
	# Get file type
	r = regexpr("\\.[a-zA-Z]*$", filename)
	if(r == -1) stop("Improper filename")
	ending = substr(filename, r + 1, r + attr(r, "match.length"))
	
	f = switch(ending,
			pdf = function(x, ...) pdf(x, ...),
			svg = function(x, ...) svg(x, ...),			
			png = function(x, ...) png(x, ...),
			jpeg = function(x, ...) jpeg(x, ...),
			jpg = function(x, ...) jpeg(x, ...),
			tiff = function(x, ...) tiff(x, compression = "lzw", ...),
			bmp = function(x, ...) bmp(x, ...),
			stop("File type should be: pdf, svg, png, bmp, jpg, tiff")
	)
	
	args <- c(list(filename), list(...))	
	if( !missing(width) ){
		args$width <- as.numeric(width)
		args$height <- as.numeric(height)
		if( !ending %in% c('pdf','svg') && is.null(args[['res']]) ){
			args$units <- "in"
			args$res <- 300
		}
	}
	do.call('f', args)	
}

#gt <- function(){
#	
#	x <- rmatrix(20, 10)
#	z <- unit(0.1, "npc")
#	w <- unit(0.4, "npc")
#	h <- unit(0.3, "npc")
#	lo <- grid.layout(nrow = 7, ncol = 6
#			, widths = unit.c(z, z, w, z, z, z)
#			, heights = unit.c(z, z,  z, h, z, z, z))
#	
#	nvp <- 0
#	on.exit( upViewport(nvp) )
#	
#	u <- vplayout(NULL)
#	vname <- function(x) basename(tempfile(x))
#
#	hvp <- viewport( name=u, layout = lo)
#	pushViewport(hvp)
#	nvp <- nvp + 1
#	
#	pushViewport(viewport(layout.pos.row = 4, layout.pos.col = 3, name='test'))
#	#vplayout('mat')
#	nvp <- nvp + 1
#	
#	grid.rect()
#	NULL
#}

#gt2 <- function(){
#	
#	x <- rmatrix(10, 5)
#	lo(NULL, NULL, nrow(x), ncol(x), cellheight = NA, cellwidth = NA
#			, treeheight_col=0, treeheight_row=0, legend=FALSE, main = NULL, sub = NULL, info = NULL
#			, annTracks=list(colors=NA, annRow=NA, annCol=NA), annotation_legend=FALSE
#			, fontsize=NULL, fontsize_row=NULL, fontsize_col=NULL)
#	
#	#vplayout('mat')
#	vname <- function(x) basename(tempfile(x))
#	pushViewport(viewport(layout.pos.row = 4, layout.pos.col = 3, name=vname('test')))
#	print(current.vpPath())
#	grid.rect()
#	upViewport(2)
#	NULL
#}

d <- function(x){
	
	if( is.character(x) ) x <- rmatrix(dim(x))
	
	nvp <- 0
	on.exit(upViewport(nvp), add=TRUE)
	lo <- grid.layout(nrow = 4, ncol = 3)
	hvp <- viewport( name=basename(tempfile()), layout = lo)
	pushViewport(hvp)
	nvp <- nvp + 1
	
	pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
	nvp <- nvp + 1
	w = convertWidth(unit(1, "npc"), "bigpts", valueOnly = T) / 10
	h = convertHeight(unit(1, "npc"), "bigpts", valueOnly = T) / 10
	grid.rect()
	upViewport()
	nvp <- nvp - 1
	
	pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
	nvp <- nvp + 1
	# add inner padding viewport
	pushViewport( viewport(x=0,y=0,width=0.9,height=0.9,just=c("left", "bottom")) )
	nvp <- nvp + 1
	
	( opar <- par(plt = gridPLT(), new = TRUE) )
	on.exit(par(opar), add=TRUE)
		
	hc <- hclust(dist(x))
	plot(as.dendrogram(hc), xaxs="i", yaxs="i", axes=FALSE, leaflab="none")
	
	invisible(basename(tempfile()))
}

heatmap_motor = function(matrix, border_color, cellwidth, cellheight
	, tree_col, tree_row, treeheight_col, treeheight_row
	, filename=NA, width=NA, height=NA
	, breaks, color, legend, txt = NULL
	, annTracks, annotation_legend=TRUE, cexAnn = NA
	, new=TRUE, fontsize, fontsize_row, fontsize_col
	, main=NULL, sub=NULL, info=NULL
    , layout = NULL
	, verbose=getOption('verbose')
	, gp = gpar()){

	annotation_colors <- annTracks$colors
	row_annotation <- annTracks$annRow
	annotation <- annTracks$annCol
	writeToFile <- !is.na(filename)
    # annotation track size
    if( length(cexAnn) == 1L ) cexAnn <- c(cexAnn, cexAnn)
    cexAnn[is.na(cexAnn)] <- 1 
    
	# open graphic device (dimensions will be changed after computation of the correct height)
	if( writeToFile ){
		gfile(filename)
		on.exit(dev.off())
	}
	 
	# identify the plotting context: base or grid
	#NB: use custom function current.vpPath2 instead of official 
	# grid::current.vpPath as this one creates a new page when called 
	# on a fresh graphic device	
	vpp <- current.vpPath_patched()
	if( is.null(vpp) ){ # we are at the root viewport
		if( verbose ) message("Detected path: [ROOT]")
		mf <- par('mfrow')		
		#print(mf)
		# if in in mfrow/layout context: setup fake-ROOT viewports with gridBase
		# and do not call plot.new as it is called in grid.base.mix. 
		new <- if( !identical(mf, c(1L,1L)) ){ 
			if( verbose ) message("Detected mfrow: ", mf[1], " - ", mf[2], ' ... MIXED')
			opar <- grid.base.mix(trace=verbose>1)
			on.exit( grid.base.mix(opar) )
			FALSE
		}		
		else{
			if( verbose ){
				message("Detected mfrow: ", mf[1], " - ", mf[2])
				message("Honouring ", if( missing(new) ) "default "
						,"argument `new=", new, '` ... '
						, if( new ) "NEW" else "OVERLAY")				
			}
			new
		}
	}else{
		if( verbose ) message("Detected path: ", vpp)
		# if new is not specified: change the default behaviour by not calling 
		# plot.new so that drawing occurs in the current viewport
		if( missing(new) ){
			if( verbose ) message("Missing argument `new` ... OVERLAY")
			new <- FALSE
		}else if( verbose ) message("Honouring argument `new=", new, '` ... '
									, if( new ) "NEW" else "OVERLAY")
	}
	# reset device if necessary or requested
	if( new ){
		if( verbose ) message("Call: plot.new")
		#grid.newpage()
		plot.new()
	}	
	
	# define grob for main 
	mainGrob <- if( !is.null(main) && !is.grob(main) ) textGrob(main, gp = c_gpar(gp, fontsize = 1.2 * fontsize, fontface="bold"))
	subGrob <- if( !is.null(sub) && !is.grob(sub) ) textGrob(sub, gp = c_gpar(gp, fontsize = 0.8 * fontsize))
	infoGrob <- if( !is.null(info) && !is.grob(info) ){
#		infotxt <- paste(strwrap(paste(info, collapse=" | "), width=20), collapse="\n")
		grobTree(gList(rectGrob(gp = gpar(fill = "grey80"))
			,textGrob(paste(info, collapse=" | "), x=unit(5, 'bigpts'), y=0.5, just='left', gp = c_gpar(gp, fontsize = 0.8 * fontsize))))
	}
	
	# Set layout
	glo = lo(coln = colnames(matrix), rown = rownames(matrix), nrow = nrow(matrix), ncol = ncol(matrix)
	, cellwidth = cellwidth, cellheight = cellheight
	, treeheight_col = treeheight_col, treeheight_row = treeheight_row
	, legend = legend
	, annTracks = annTracks, annotation_legend = annotation_legend, cexAnn = cexAnn
	, fontsize = fontsize, fontsize_row = fontsize_row, fontsize_col = fontsize_col
    , layout = layout
	, main = mainGrob, sub = subGrob, info = infoGrob, gp = gp)
    
    # extract options
    loptions <- glo$layout$options 
	
	# resize the graphic file device if necessary
	if( writeToFile ){		
		if( verbose ) message("Compute size for file graphic device")
		m <- par('mar')
		if(is.na(height))
			height <- glo$height
		if(is.na(width))
			width <- glo$width
		
		dev.off()
		if( verbose ) message("Resize file graphic device to: ", width, " - ", height)
		gfile(filename, width=width, height=height)
		# re-call plot.new if it was called before
		if( new ){
			if( verbose ) message("Call again plot.new")
			op <- par(mar=c(0,0,0,0))
			plot.new()
			par(op)
		}		
		if( verbose ) message("Push again top viewport")
		# repush the layout
		pushViewport(glo$vp)
		if( verbose ) grid.rect(width=unit(glo$width, 'inches'), height=unit(glo$height, 'inches'), gp = gpar(col='blue'))
	}

	#grid.show.layout(glo$layout); return()
	mindim <- glo$mindim
	# Omit border color if cell size is too small 
	if(mindim < 3) border_color = NA

	# Draw tree for the columns
	if (!is_NA(tree_col) &&  treeheight_col != 0 && vplayout('ctree') ){
		draw_dendrogram(tree_col, horizontal = FALSE, flip = loptions$dendrogram$flip$v)
		upViewport()
	}

	# Draw tree for the rows
	if(!is_NA(tree_row) && treeheight_row !=0 && vplayout('rtree') ){
		draw_dendrogram(tree_row, horizontal = TRUE, flip = loptions$dendrogram$flip$h)
		upViewport()
	}

    # recompute margin fontsizes
    fontsize_row <- convertUnit(min(unit(fontsize_row, 'points'), unit(0.6*glo$cellheight, 'bigpts')), 'points')
    fontsize_col <- convertUnit(min(unit(fontsize_col, 'points'), unit(0.6*glo$cellwidth, 'bigpts')), 'points')
    
	# Draw matrix
	if( vplayout('mat') ){
    	draw_matrix(matrix, border_color, txt = txt, gp = gpar(fontsize = fontsize_row))
    	#d(matrix)
    	#grid.rect()
    	upViewport()
    }

	# Draw colnames
	if(length(colnames(matrix)) != 0 && vplayout('cnam') ){
		draw_colnames(colnames(matrix), gp = c_gpar(gp, fontsize = fontsize_col))
		upViewport()
	}
	
	# Draw rownames
	if(length(rownames(matrix)) != 0 && vplayout('rnam') ){
		draw_rownames(rownames(matrix), gp = c_gpar(gp, fontsize = fontsize_row))
		upViewport()
	}

	# Draw annotation tracks
	if( !is_NA(annotation) && vplayout('cann') ){
		draw_annotations(annotation, border_color, cex = cexAnn[2L])
		upViewport()
	}	
	
	# add row annotations if necessary	
	if ( !is_NA(row_annotation) && vplayout('rann') ) {
		draw_annotations(row_annotation, border_color, horizontal=FALSE, cex = cexAnn[1L])
		upViewport()
	}
	
	# Draw annotation legend
	if( annotation_legend && !is_NA(annotation_colors) && vplayout('aleg') ){
		draw_annotation_legend(annotation_colors, border_color, gp = c_gpar(gp, fontsize = fontsize))
		upViewport()
	}

	# Draw legend
	if(!is_NA(legend) && vplayout('leg') ){
		draw_legend(color, breaks, legend, gp = c_gpar(gp, fontsize = fontsize), opts = loptions$legend)
		upViewport()
	}

	# Draw main title
	if(!is.null(mainGrob) && vplayout('main') ){
		grid.draw(mainGrob)
		upViewport()
	}
	
	# Draw subtitle
	if(!is.null(subGrob) && vplayout('sub') ){
		grid.draw(subGrob)
		upViewport()
	}
	
	# Draw info
	if(!is.null(infoGrob) && vplayout('info') ){
		grid.draw(infoGrob)
		upViewport()
	}
		
	# return current vp tree
	#ct <- current.vpTree()
	#print(current.vpPath())
	upViewport()
	#popViewport()
	
	# grab current grob and return
#	gr <- grid.grab()
#	grid.draw(gr)
	#ct
	NULL
}

generate_breaks = function(x, n, center=NA){
	if( missing(center) || is_NA(center) )
		seq(min(x, na.rm = T), max(x, na.rm = T), length.out = n + 1)
	else{ # center the breaks on the requested value
		n2 <- ceiling((n+0.5)/2)
		M <- max(abs(center - min(x, na.rm = TRUE)), abs(center - max(x, na.rm = TRUE)))
		lb <- seq(center-M, center, length.out = n2)
		rb <- seq(center, center+M, length.out = n2)
		c(lb, rb[-1])
	}
}

scale_vec_colours = function(x, col = rainbow(10), breaks = NA){
    return(col[cut(x, breaks = breaks, include.lowest = T, labels = FALSE)])	
}

scale_colours = function(mat, col = rainbow(10), breaks = NA){
	mat = as.matrix(mat)
	return(matrix(scale_vec_colours(as.vector(mat), col = col, breaks = breaks), nrow(mat), ncol(mat), dimnames = list(rownames(mat), colnames(mat))))
}

cutheight <- function(x, n){
	
	# exit early if n <=1: nothing to do
	if( n <=1 ) return( attr(x, 'height') )
	
	res <- NULL
	.heights <- function(subtree, n){		
		if( is.leaf(subtree) ) return()
		if (!(K <- length(subtree))) 
			stop("non-leaf subtree of length 0")		
		# extract heights from each subtree 
		for( k in 1:K){
			res <<- c(res, attr(subtree[[k]], 'height'))			
		}
		
		# continue only if there is not yet enough subtrees
		if( length(res) < n ){
			for( k in 1:K){
				.heights(subtree[[k]], n)
			}
		}
	}
	
	# extract at least the top h heights
	.heights(x, n)
	
	# sort by decreasing order
	res <- sort(res, decreasing=TRUE)
	
	res[n-1]
}

.dendextend <- function(){
    qlibrary('colorspace')
    qlibrary('dendextend')
}

#' Fade Out the Upper Branches from a Dendrogram
#' 
#' @param x a dendrogram
#' @param n the number of groups
#' 
#' @import digest
#' @keywords internal
#' @import dendextend
cutdendro <- function(x, n){
	
    use_dendextend <- FALSE
    if( isString(n) || is.list(n) ){
        .dendextend()
        use_dendextend <- TRUE
        spec <- str_match(n[[1L]], "^#(-?[0-9]+)(@([0-9, ]+))?([|])?([!])?")
#        print(spec)
        n0 <- n
        n <- as.integer(spec[, 2L])
        n_cl <- abs(n)
            
        # border
        bd <- if( nzchar(spec[, 5L]) ) 8 else NA
        col <- if( nzchar(spec[, 6L]) ) NA
        # subset
        w <- if( nzchar(spec[, 4L]) ) eval(parse(text = sprintf("c(%s)", spec[, 4L]))) else seq(n_cl)
        cl_spec <- expand_list(c(list(k = n_cl), as.list(n0[-1L])), col = col, border = bd, which = w)
        # complete text
        if( length(cl_spec$which) == 1L && !is.null(cl_spec$text) ) cl_spec$text <- rep(cl_spec$text, length.out = cl_spec$k)
        if( length(bad_i <- which(cl_spec$which > n_cl)) ){
            stop(gettextf("all indexes in `@` or `which` must be between 1 and %d [bad indexes: %s]", n_cl, str_out(w[bad_i], 5, total = length(bad_i) > 5)))
        }
    }
    
    n0 <- n
    n <- abs(n)
    
    
	# exit early if n <=1: nothing to do
	if( n <= 1 ) return(x)
		
	# add node digest ids to x
	x <- dendrapply(x, function(n){
			attr(n, 'id') <- digest(attributes(n))
			n
	})

    # cut x in n groups
	# find the height where to cut
    if( use_dendextend && !is_NA(cl_spec$col) ){
        if( !is.null(cl_spec$col) ){
            # complete colors to avoid warnings
            if( length(cl_spec$which) == 1L ) cl_spec$col <- rep(cl_spec$col, length.out = cl_spec$k)
            
            x <- x %>% set("branches_k_color", k = n, value = cl_spec$col)
        }else x <- x %>% set("branches_k_color", k = n)
    }
	h <- cutheight(x, n)
	cfx <- cut(x, h)
	# get the ids of the upper nodes
	ids <- sapply(cfx$lower, function(sub) attr(sub, 'id'))
    
    # highlight the upper branches with dot lines
    if( use_dendextend ){
        dts <- c(lty=2, lwd=2)
        cl_spec$col <- unlist(sapply(cfx$lower, function(sub) attr(sub, 'edgePar')[['col']]))
        
    }else{
        dts <- c(lty=2, lwd=1.2, col = 8)
    }
    
    a <- dendrapply(x, function(node){
		a <- attributes(node) 
        if( n0 > 0 && (a$id %in% ids || (!is.leaf(node) && any(c(attr(node[[1]], 'id'), attr(node[[2]], 'id')) %in% ids))) # dash upper tree
            || (n0 < 0 && (!a$id %in% ids && !(!is.leaf(node) && any(c(attr(node[[1]], 'id'), attr(node[[2]], 'id')) %in% ids)))) ){# dash lower cluster trees
        
                if( use_dendextend ){
                    attr(node, 'edgePar') <- c(dts, attr(node, 'edgePar')['col'])
                } else attr(node, 'edgePar') <- dts
                
        }
#        else if( !is.leaf(node) && use_dendextend && attr(node[[1]], 'id') %in% ids && attr(node[[2]], 'id') %in% ids && nzchar(spec[, 3L]) ){
#            print(ids)
#            print(attr(node[[1]], 'height'))
#            print(attr(node, 'height'))
#            node[[1]] %>% raise.dendrogram(attr(node[[1]], 'height') - attr(node[[1]], 'height'))
#            node[[2]] <- node[[2]] %>% raise.dendrogram(- attr(node[[2]], 'height'))
#        }
        node
	})
    
    # store cluster specs
    if( use_dendextend ){
        attr(a, 'cluster.spec') <- cl_spec 
    }

    
    a
}

# internal class definition for 
as_treedef <- function(x, ...){
	res <- if( is(x, 'hclust') )
				list(dendrogram=as.dendrogram(x), dist.method=x$dist.method, method=x$method)
			else list(dendrogram=x, ...)
	class(res) <- "aheatmap_treedef"
	res
}
rev.aheatmap_treedef <- function(x){
	x$dendrogram <- rev(x$dendrogram)
	x
}

is_treedef <- function(x) is(x, 'aheatmap_treedef')

isLogical <- function(x) isTRUE(x) || identical(x, FALSE)


# Convert an index vector usable on the subset data into one usable on the 
# original data
subset2orginal_idx <- function(idx, subset){
	if( is.null(subset) || is.null(idx) ) idx
	else{
		res <- subset[idx]
		attr(res, 'subset') <- idx
		res
	}
}

#' Cluster Matrix Rows in Annotated Heatmaps
#'
#' @param mat original input matrix that has already been appropriately subset in 
#' the caller function (\code{aheatmap})
#' @param param clustering specifications
#' @param distfun Default distance method/function
#' @param hclustfun Default clustering (linkage) method/function
#' @param reorderfun Default reordering function
#' @param na.rm Logical that specifies if NA values should be removed
#' @param subset index (integer) vector specifying the subset indexes used to 
#' subset mat. This is required to be able to return the original indexes. 
#' 
#' @keywords internal
cluster_mat = function(mat, param, distfun, hclustfun, reorderfun, na.rm=TRUE, subset=NULL, verbose = FALSE){
	
	# do nothing if an hclust object is passed	
	parg <- deparse(substitute(param))
		
	Rowv <- 
	if( is(param, 'hclust') || is(param, 'dendrogram') ){ # hclust or dendrograms are honoured
		res <- as_treedef(param)
		
		# subset if requested: convert into an index vector
		# the actuval subsetting is done by first case (index vector)
		if( !is.null(subset) ){
			warning("Could not directly subset dendrogram/hclust object `", parg 
					,"`: using subset of the dendrogram's order instead.")
			# use dendrogram order instead of dendrogram itself
			param <- order.dendrogram(res$dendrogram)
			
		}else # EXIT: return treedef
			return(res)
	}else if( is(param, 'silhouette') ){ # use silhouette order 
        si <- sortSilhouette(param)
        param <- attr(si, 'iOrd')
    }

	# index vectors are honoured
	if( is.integer(param) && length(param) > 1 ){
		
		# subset if requested: reorder the subset indexes as in param
		if( !is.null(subset) )
			param <- order(match(subset, param))			
		
		param
	}else{ # will compute dendrogram (NB: mat was already subset before calling cluster_mat)
		
        use.cutdendro <- function(x){
            is.integer(x) || (isString(x) && grepl("^#-?[0-9]+", x)) || (is.list(x) && use.cutdendro(x[[1L]]))
        }
        
		param <- 
		if( use.cutdendro(param) ){
            param
		}else if( is.null(param) || isLogical(param) ) # use default reordering by rowMeans
			rowMeans(mat, na.rm=na.rm)
		else if( is.numeric(param) ){ # numeric reordering weights
			# subset if necessary
			if( !is.null(subset) )
				param <- param[subset]
			param
		}else if( is.character(param) || is.list(param) ){
			
			if( length(param) == 0 )
				stop("aheatmap - Invalid empty character argument `", parg, "`.")
			
			# set default names if no names were provided
			if( is.null(names(param)) ){			
				if( length(param) > 3 ){
					warning("aheatmap - Only using the three first elements of `", parg, "` for distfun and hclustfun respectively.")
					param <- param[1:3]
				}			
				n.allowed <- c('distfun', 'hclustfun', 'reorderfun')
				names(param) <- head(n.allowed, length(param))
			}
			
			# use the distance passed in param 
			if( 'distfun' %in% names(param) ) distfun <- param[['distfun']]
			# use the clustering function passed in param
			if( 'hclustfun' %in% names(param) ) hclustfun <- param[['hclustfun']]
			# use the reordering function passed in param
			if( 'reorderfun' %in% names(param) ) reorderfun <- param[['reorderfun']]
			
			TRUE
		}else		
			stop("aheatmap - Invalid value for argument `", parg, "`. See ?aheatmap.")
		
		# compute distances
		d <- if( isString(distfun) ){
			distfun <- distfun[1]
			
            corr.methods <- c("pearson", "kendall", "spearman")
			av <- c("correlation", corr.methods, "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")
			i <- pmatch(distfun, av)
			if( is_NA(i) )			
				stop("aheatmap - Invalid dissimilarity method, must be one of: ", str_out(av, Inf))
			
			distfun <- av[i]
            if(distfun == "correlation") distfun <- 'pearson'
			if(distfun %in% corr.methods){ # distance from correlation matrix 
                if( verbose ) message("Using distance method: correlation (", distfun, ')')
                d <- dist(1 - cor(t(mat), method = distfun))
                attr(d, 'method') <- distfun
                d 
            }else{
                if( verbose ) message("Using distance method: ", distfun)
                dist(mat, method = distfun)
            }
	
		}else if( is(distfun, "dist") ){
			if( verbose ) message("Using dist object: ", distfun)
			distfun
		}else if( is.function(distfun) ){
            if( verbose ) message("Using custom dist function")
			distfun(mat)
		}else
			stop("aheatmap - Invalid dissimilarity function: must be a character string, an object of class 'dist', or a function")
	
		# do hierarchical clustering 
		hc <- if( is.character(hclustfun) ){
			
			av <- c('ward.D', 'ward.D2', 'ward', 'single', 'complete', 'average', 'mcquitty', 'median', 'centroid')
			i <- pmatch(hclustfun, av)
			if( is.na(i) )
				stop("aheatmap - Invalid clustering method, must be one of: ", paste("'", av, "'", sep='', collapse=', '))
			
			hclustfun <- av[i]
			if( verbose ) message("Using clustering method: ", hclustfun)
			hclust(d, method=hclustfun)
			
		}else if( is.function(hclustfun) )
			hclustfun(d)
		else
			stop("aheatmap - Invalid clustering function: must be a character string or a function")
	
		#convert into a dendrogram
		dh <- as.dendrogram(hc)
	
		# customize the dendrogram plot: highlight clusters
		if( use.cutdendro(param) )
			dh <- cutdendro(dh, param)						
		else if( is.numeric(param) && length(param)==nrow(mat) ) # reorder the dendrogram if necessary
			dh <- reorderfun(dh, param)
		
		# wrap up into a aheatmap_treedef object
		as_treedef(dh, dist.method=hc$dist.method, method=hc$method)
	}
}

#scale_rows = function(x){
#	m = apply(x, 1, mean)
#	s = apply(x, 1, sd)
#	return((x - m) / s)
#}

scale_mat = function(x, scale, na.rm=TRUE){
	
	av <- c("none", "row", "column", 'r1', 'c1')
	i <- pmatch(scale, av)	
	if( is_NA(i) )
		stop("scale argument shoud take values: 'none', 'row' or 'column'")
	scale <- av[i]
		
	switch(scale, none = x
		, row = {
			x <- sweep(x, 1L, rowMeans(x, na.rm = na.rm), check.margin = FALSE)
			sx <- apply(x, 1L, sd, na.rm = na.rm)
			sweep(x, 1L, sx, "/", check.margin = FALSE)
		}
		, column = {
			x <- sweep(x, 2L, colMeans(x, na.rm = na.rm), check.margin = FALSE)
			sx <- apply(x, 2L, sd, na.rm = na.rm)
			sweep(x, 2L, sx, "/", check.margin = FALSE)
		}
		, r1 = sweep(x, 1L, rowSums(x, na.rm = na.rm), '/', check.margin = FALSE)
		, c1 = sweep(x, 2L, colSums(x, na.rm = na.rm), '/', check.margin = FALSE)
		)	
}

.Rd.seed <- new.env()

round.pretty <- function(x, min=2){
	
	if( is.null(x) ) return(NULL)		
	n <- 0
	y <- round(sort(x), n)
	if( all(diff(y)==0) ) return( round(x, min) ) 
	while( any(diff(y)==0) ){
		n <- n+1
		y <- round(sort(x), n)
	}	
	dec <- max(min,n)
	round(x, dec)
}

generate_annotation_colours = function(annotation, annotation_colors, seed=TRUE){
	if( is_NA(annotation_colors) ){
		annotation_colors = list()
	}
	# use names from annotations if necessary/possible
	if( length(annotation_colors) > 0L 
		&& length(annotation_colors) <=  length(annotation) 
		&& is.null(names(annotation_colors)) ){
		names(annotation_colors) <- head(names(annotation), length(annotation_colors))
	}
	
	count = 0
	annotationLevels <- list()	
	anames <- names(annotation)
	sapply(seq_along(annotation), function(i){
		a <- annotation[[i]]
		if( class(annotation[[i]]) %in% c("character", "factor")){			
			# convert to character vector
			a <- if( is.factor(a) ) levels(a) else unique(a)
			count <<- count + nlevels(a)
			# merge if possible
			if( !is.null(anames) && anames[i]!='' )	
				annotationLevels[[anames[i]]] <<- unique(c(annotationLevels[[anames[i]]], a))
			else 
				annotationLevels <<- c(annotationLevels, list(a))			
		}else
			annotationLevels <<- c(annotationLevels, annotation[i])
	})
	annotation <- annotationLevels
	#str(annotationLevels)

	factor_colors = hcl(h = seq(1, 360, length.out = max(count+1,20)), 100, 70)	
		
	# get random seeds to restore/update on exit
	rs <- RNGseed()
	on.exit({
		# update local random seed on exit
		.Rd.seed$.Random.seed <- getRNG()
		# restore global random seed
		RNGseed(rs)
	})
	# restore local random seed if it exists 
	if( !is.null(.Rd.seed$.Random.seed) )
		setRNG(.Rd.seed$.Random.seed)
	# set seed and restore on exit
	if( isTRUE(seed) ){
		# reset .Random.seed to a dummy RNG in case the current kind is user-supplied:
		# we do not want to draw even once from the current RNG
		setRNG(c(401L, 0L, 0L))
		set.seed(12345, 'default', 'default')
	}
	factor_colors <- sample(factor_colors)
#	pal(factor_colors); stop("sasa")
	
	res_colors <- list()
	for(i in 1:length(annotation)){
		ann <- annotation[[i]]
		aname <- names(annotation)[i]
		# skip already generated colors
		acol_def <- res_colors[[aname]]
		if( !is.null(acol_def) ) next;
		acol <- annotation_colors[[aname]]
		if( is.null(acol) ){
			res_colors[[aname]] <-
			if( class(annotation[[i]]) %in% c("character", "factor")){
				lev <- ann
				ind = 1:length(lev)
				acol <- setNames(factor_colors[ind], lev)
				factor_colors = factor_colors[-ind]
				# conserve NA value								
				acol[which(is.na(names(acol)))] <- NA
				acol
			}
			else{
				h = round(runif(1) * 360)
				rg <- range(ann, na.rm=TRUE)
				if( rg[1] == rg[2] ) rg <- sort(c(0, rg[1]))
				setNames(rev(sequential_hcl(2, h, l = c(50, 95))), round.pretty(rg))
			}
		
		}else{
						
			acol <- 
			if( length(acol) == 1 && grepl("^\\$", acol) ) # copy colors from other columns if the spec starts with '$'
				annotation_colors[[substr(acol, 2, nchar(acol))]]
			else if( !is.numeric(ann) ){				
				local({ #do this locally so that it does not affect `ann`
					
					# subset to the levels for which no colour has already been defined
					lev <- ann
					# subset to the levels for which no colour has already been defined
#					idx <- which(!lev %in% names(acol_def) & !is.na(lev))
#					lev <- lev[idx]
#					#idx <- idx + length(acol_def)
#					
#					if( length(lev) == 0L ) acol_def # nothing to add
#					else
					{
						# convert to a palette of the number of levels if necessary
						nl <- length(lev)
						acol <- ccPalette(acol, nl)
						if( is.null(names(acol)) )
							names(acol) <- lev
						c(acol_def, acol)
					}
				})
			}else{
				acol <- ccPalette(acol)
				if( is.null(names(acol)) )
					names(acol) <- round.pretty(seq(min(ann, na.rm=TRUE), max(ann, na.rm=TRUE), length.out=length(acol)))
			
				acol
			}
			
			# update the colors if necessary
			if( !is.null(acol) )
				res_colors[[aname]] <- acol
		}
		
		# store type information
		attr(res_colors[[aname]], 'afactor') <- !is.numeric(ann) 		
				
	}	
	
	# return ordered colors as the annotations
	res_colors[names(annotation)[!duplicated(names(annotation))]]	
}

# Create row/column names
generate_dimnames <- function(x, n, ref){
	if( is_NA(x) ) NULL
	else if( length(x) == n ) x
	else if( identical(x, 1) || identical(x, 1L) ) 1L:n
	else if( isString(x) ){
		regexp <- "^/(.+)/([0-9]+)?$"
		if( grepl(regexp, x) ){
			x <- str_match(x, regexp)
			p <- x[1,2]
			n <- if( x[1, 3] != '' ) as.numeric(x[1, 2]) else 2L
			s <- str_match(ref, p)[, n]
			ifelse(is.na(s), ref, s)
		}
		else paste(x, 1L:n, sep='')
		#print(str_match_all(x, "^/(([^%]*)(%[in])?)+/$"))
	}
	else stop("aheatmap - Invalid row/column label. Possible values are:"
			, " NA, a vector of correct length, value 1 (or 1L) or single character string.")
}


.make_annotation <- function(x, ord=NULL){
	# convert into a data.frame if necessary
	if( !is.data.frame(x) ){
		x <- if( is(x, 'ExpressionSet') ) Biobase::pData(x)
				else if( is.factor(x) || is.character(x) ) data.frame(Factor=x)
				else if( is.numeric(x) ) data.frame(Variable=x)
				else
					stop("aheatmap - Invalid annotation argument `", substitute(x), "`: must be a data.frame, a factor or a numeric vector")		
	}
	# reorder if necessary
	if( !is.null(ord) )
		x <- x[ord, , drop = F]
	
	# return modifed object 
	x
}

renderAnnotations <- function(annCol, annRow, annotation_colors, verbose=getOption('verbose')){
	
	# concatenate both col and row annotation
	annotation <- list()
	if( is_NA(annotation_colors) ) annotation_colors <- list() 
	nc <- length(annCol)
	nr <- length(annRow)
	flag <- function(x, f){ if( missing(f) ) attr(x, 'flag') else{ attr(x, 'flag') <- f; x} }
	if( !is_NA(annCol) ) annotation <- c(annotation, sapply(as.list(annCol), flag, 'col', simplify=FALSE))
	if( !is_NA(annRow) ) annotation <- c(annotation, sapply(as.list(annRow), flag, 'row', simplify=FALSE))
		
	if( length(annotation) == 0 ) return( list(annCol=NA, annRow=NA, colors=NA) )
	
	# generate the missing name
	n <- names(annotation)
	xnames <- paste('X', 1:length(annotation), sep='')
	if( is.null(n) ) names(annotation) <- xnames
	else names(annotation)[n==''] <- xnames[n==''] 
		
	# preprocess the annotation color links
	if( !is.null(cnames <- names(annotation_colors)) ){
		m <- str_match(cnames, "^@([^{]+)\\{([^}]+)\\}")		
		apply(m, 1L, function(x){
			# skip unmatched names
			if( is_NA(x[1]) ) return()
			
			acol <- annotation_colors[[x[1]]]
			# rename both annotation and annotation_colors if necessary
			if( x[2] != x[3] ){
				annotation[[x[3]]] <<- annotation[[x[2]]]				
				annotation[[x[2]]] <<- NULL
				if( !is_NA(acol) )
					annotation_colors[[x[3]]] <<- acol
				annotation_colors[[x[1]]] <<- NULL
			}
		})
	}
	
#	message("### ANNOTATION ###"); print(annotation)
#	message("### ANNOTATION COLORS ###"); print(annotation_colors)
	
	if( verbose ) message("Generate column annotation colours")
	annotation_colors  <- generate_annotation_colours(annotation, annotation_colors)
	if( verbose > 2 ){
		message("### Annotation colors ###")
		print(annotation_colors)
		message("#########################")
	}
	
	# bind each annotation with its respective color and regroup into column and row annotation
	res <- list()
	lapply(seq_along(annotation), function(i){
		aname <- names(annotation)[i]
		acol <- annotation_colors[[aname]]		
		if( is.null(acol) )
			stop("aheatmap - No color was defined for annotation '", aname, "'.")
		attr(annotation[[i]], 'color') <- acol
				
		# put into the right annotation list
		if( flag(annotation[[i]]) == 'col' ) res$annCol <<- c(res$annCol, annotation[i])
		else res$annRow <<- c(res$annRow, annotation[i])
	})
	
	res$annCol <- if( !is.null(res$annCol) ) convert_annotations(res$annCol, annotation_colors) else NA
	res$annRow <- if( !is.null(res$annRow) ) convert_annotations(res$annRow, annotation_colors) else NA
	res$colors <- annotation_colors
	
	# return result list
	res
	
}

# set/get special annotation handlers
specialAnnotation <- local({
	.empty <- list(list(), list())
	.cache <- .empty
	function(margin, name, fun, clear=FALSE){
	
		if( isTRUE(clear) ){
			if( nargs() > 1L )
				stop("Invalid call: no other argument can be passed when `clear=TRUE`")
			.cache <<- .empty
			return()
		}
		
		if( missing(name) && missing(fun) ){
			return(.cache[[margin]])
		}else if( is.list(name) ){
			.cache[[margin]] <<- c(.cache[[margin]], name)
		}else if( missing(fun) ){
			return(.cache[[margin]][[name]])
		}else{
			.cache[[margin]][[name]] <<- fun
		}
	}
})

# Converts Subset Specification into Indexes
subset_index <- function(x, margin, subset){
	
	# if null then do nothing
	if( is.null(subset) ) return( NULL )
	
	# get dimension
	n <- dim(x)[margin]
	dn <- dimnames(x)[[margin]]
	dt <- if( margin == 1L ) "rows" else "columns"
	
	so <- deparse(substitute(subset))
	if( length(subset) == 0 )
		stop("Invalid empty subset object `", so, "`")
	
	subIdx <- 
	if( is.logical(subset) ){			
		if( length(subset) != n ){
			if( n %% length(subset) == 0 )
				subset <- rep(subset, n / length(subset))
			else
				stop("Invalid length for logical subset argument `", so, "`: number of ", dt, " ["
						, n, "] is not a multiple of subset length [",length(subset),"].")
		}
		
		# convert into indexes
		which(subset)					
	}
	else if( is.integer(subset) || is.character(subset) ){
		if( length(subset) > n )
			stop("Invalid too long integer/character subset argument `", so
				, "`: length must not exceed the number of ", dt, " [", n, "].")
	
		if( anyDuplicated(subset) )
			warning("Duplicated index or name in subset argument `", so, "`.")
		
		# for character argument: match against dimname to convert into indexes
		if( is.character(subset) ){
			if( is.null(dn) )
				stop("Could not subset the ", dt, " with a character subset argument `", so, "`: no "
						, if( margin == 1L ) "rownames" else "colnames"
						, " are available.")
			msubset <- match(subset, dn)
			nas <- is.na(msubset)
			if( any(nas) ){
				warning("Mismatch in character subset argument `", so
						,"`: Could not find ", sum(nas), " out of ", length(subset), " names ("
						, paste("'", head(subset[nas], 5), "'", sep='', collapse=', ')
						, if( sum(nas) > 5 ) ", ... ", ").")
				msubset <- msubset[!nas]
			}
			subset <- msubset
		}
		subset
	}else
		stop("Invalid subset argument `", so, "`: should be a logical, integer or character vector.")
	
	# return the indexes sorted
	sort(subIdx)
}

#' Annotated Heatmaps
#' 
#' The function \code{aheatmap} plots high-quality heatmaps, with a detailed legend 
#' and unlimited annotation tracks for both columns and rows. 
#' The annotations are coloured differently according to their type
#' (factor or numeric covariate).
#' Although it uses grid graphics, the generated plot is compatible with base 
#' layouts such as the ones defined with \code{'mfrow'} or \code{\link{layout}}, 
#' enabling the easy drawing of multiple heatmaps on a single a plot -- at last!.
#'
#' The development of this function started as a fork of the function 
#' \code{pheatmap} from the \pkg{pheatmap} package, and provides
#' several enhancements such as:
#' \itemize{
#' \item argument names match those used in the base function \code{\link{heatmap}}; 
#' \item unlimited number of annotation for \strong{both} columns and rows, 
#' with simplified and more flexible interface;
#' \item easy specification of clustering methods and colors;
#' \item return clustering data, as well as grid grob object.
#' }
#' 
#' Please read the associated vignette for more information and sample code.
#' 
#' @section PDF graphic devices: if plotting on a PDF graphic device -- started with \code{\link{pdf}}, 
#' one may get generate a first blank page, due to internals of standard functions from 
#' the \pkg{grid} package that are called by \code{aheatmap}.
#' The \pkg{NMF} package ships a custom patch that fixes this issue.
#' However, in order to comply with CRAN policies, the patch is \strong{not} applied by default 
#' and the user must explicitly be enabled it.
#' This can be achieved on runtime by either setting the NMF specific option 'grid.patch' 
#' via \code{nmf.options(grid.patch=TRUE)}, or on load time if the environment variable 
#' 'R_PACKAGE_NMF_GRID_PATCH' is defined and its value is something that is not equivalent 
#' to \code{FALSE} (i.e. not '', 'false' nor 0).
#' 
#' @param x numeric matrix of the values to be plotted.
#' An \code{\link[Biobase:ExpressionSet-class]{ExpressionSet}} objects can also 
#' be passed, in which case the expression values are plotted (\code{exprs(x)}).
#' 
#' @param color colour specification for the heatmap. Default to palette 
#' '-RdYlBu2:100', i.e. reversed palette 'RdYlBu2' (a slight modification of 
#' RColorBrewer's palette 'RdYlBu') with 100 colors.
#' Possible values are:
#' \itemize{
#' \item a character/integer vector of length greater than 1 that is directly used 
#' and assumed to contain valid R color specifications.
#' \item a single color/integer (between 0 and 8)/other numeric value 
#' that gives the dominant colors. Numeric values are converted into a pallete  
#' by \code{rev(sequential_hcl(2, h = x, l = c(50, 95)))}. Other values are 
#' concatenated with the grey colour '#F1F1F1'.   
#' \item one of RColorBrewer's palette name (see \code{\link[RColorBrewer]{display.brewer.all}})
#' , or one of 'RdYlBu2', 'rainbow', 'heat', 'topo', 'terrain', 'cm'.
#' }
#' When the coluor palette is specified with a single value, and is negative or 
#' preceded a minus ('-'), the reversed palette is used.
#' The number of breaks can also be specified after a colon (':'). For example, 
#' the default colour palette is specified as '-RdYlBu2:100'.
#' 
#' @param na.color Specifies the colour to use for \code{NA} values.
#' Setting to \code{NA} (default) produces uncoloured cells (white).
#' 
#' It can also be a list of 2 elements, with the first element specifying the color and 
#' the second a given value or a range of values (as a 2-length vector) to be forced to NA. 
#' 
#' @param breaks a sequence of numbers that covers the range of values in \code{x} and is one 
#' element longer than color vector. Used for mapping values to colors. Useful, if needed 
#' to map certain values to certain colors. If value is NA then the 
#' breaks are calculated automatically. If \code{breaks} is a single value, 
#' then the colour palette is centered on this value. 
#' 
#' @param border_color color of cell borders on heatmap, use NA if no border should be 
#' drawn.
#' 
#' @param cellwidth individual cell width in points. If left as NA, then the values 
#' depend on the size of plotting window.
#' 
#' @param cellheight individual cell height in points. If left as NA, 
#' then the values depend on the size of plotting window.
#' 
#' @param scale character indicating how the values should scaled in 
#' either the row direction or the column direction. Note that the scaling is 
#' performed after row/column clustering, so that it has no effect on the 
#' row/column ordering.
#' Possible values are: 
#' \itemize{
#' \item \code{"row"}: center and standardize each row separately to row Z-scores 
#' \item \code{"column"}: center and standardize each column separately to column Z-scores
#' \item \code{"r1"}: scale each row to sum up to one
#' \item \code{"c1"}: scale each column to sum up to one
#' \item \code{"none"}: no scaling
#' }
#' 
#' @param Rowv clustering specification(s) for the rows. It allows to specify 
#' the distance/clustering/ordering/display parameters to be used for the 
#' \emph{rows only}.
#' 
#' See section \emph{Row/column ordering and display} for details on all supported values.
#' 
#' @param Colv clustering specification(s) for the columns. It accepts the same 
#' values as argument \code{Rowv} (modulo the expected length for vector specifications),  
#' and allow specifying the distance/clustering/ordering/display parameters to 
#' be used for the \emph{columns only}.
#' 
#' \code{Colv} may also be set to \code{"Rowv"}, in which case the dendrogram 
#' or ordering specifications applied to the rows are also applied to the 
#' columns. Note that this is allowed only for square matrices, 
#' and that the row ordering is in this case by default reversed 
#' (\code{revC=TRUE}) to obtain the diagonal in the standard way 
#' (from top-left to bottom-right).
#' 
#' See section \emph{Row/column ordering and display} for details on all supported values.
#' 
#' @param revC a logical that specify if the \emph{row order} defined by 
#' \code{Rowv} should be reversed. This is mainly used to get the rows displayed 
#' from top to bottom, which is not the case by default. Its default value is 
#' computed at runtime, to suit common situations where natural ordering is a 
#' more sensible choice: no or fix ordering of the rows (\code{Rowv=NA} or an 
#' integer vector of indexes -- of length > 1), and when a symmetric ordering is 
#' requested -- so that the diagonal is shown as expected.
#' An argument in favor of the "odd" default display (bottom to top) is that the 
#' row dendrogram is plotted from bottom to top, and reversing its reorder may 
#' take a not too long but non negligeable time.
#' 
#' @param distfun default distance measure used in clustering rows and columns. 
#' Possible values are:
#' \itemize{
#' \item all the distance methods supported by \code{\link{dist}} 
#' (e.g. "euclidean" or "maximum").
#' 
#' \item all correlation methods supported by \code{\link{cor}}, 
#' such as \code{"pearson"} or \code{"spearman"}.
#' The pairwise distances between rows/columns are then computed as 
#' \code{d <- dist(1 - cor(..., method = distfun))}.
#' 
#' One may as well use the string "correlation" which is an alias for "pearson".
#' 
#' \item an object of class \code{dist} such as returned by \code{\link{dist}} or 
#' \code{\link{as.dist}}.
#' } 
#' 
#' @param hclustfun default clustering method used to cluster rows and columns.
#' Possible values are:
#' \itemize{
#' \item a method name (a character string)  supported by \code{\link{hclust}} 
#' (e.g. \code{'average'}).
#' \item an object of class \code{hclust} such as returned by \code{\link{hclust}}
#' \item a dendrogram
#' }
#' 
#' @param reorderfun default dendrogram reordering function, used to reorder the 
#' dendrogram, when either \code{Rowv} or \code{Colv} is a numeric weight vector, 
#' or provides or computes a dendrogram. It must take 2 parameters: a dendrogram, 
#' and a weight vector.
#' 
#' @param subsetRow Specification of subsetting the rows before drawing the 
#' heatmap. 
#' Possible values are:
#' \itemize{
#' \item an integer vector of length > 1 specifying the indexes of the rows to 
#' keep;
#' \item a character vector of length > 1 specyfing the names of the rows to keep.
#' These are the original rownames, not the names specified in \code{labRow}.
#' \item a logical vector of length > 1, whose elements are recycled if the 
#' vector has not as many elements as rows in \code{x}. 
#' }
#' Note that in the case \code{Rowv} is a dendrogram or hclust object, it is first  
#' converted into an ordering vector, and cannot be displayed -- and a warning is thrown. 
#' 
#' @param subsetCol Specification of subsetting the columns before drawing the 
#' heatmap. It accepts the similar values as \code{subsetRow}. See details above.
#' 
#' @param txt character matrix of the same size as \code{x}, that contains text to 
#' display in each cell. 
#' \code{NA} values are allowed and are not displayed.
#' See demo for an example. 
#' 
#' @param treeheight how much space (in points) should be used to display 
#' dendrograms. If specified as a single value, it is used for both dendrograms. 
#' A length-2 vector specifies separate values for the row and 
#' column dendrogram respectively.    
#' Default value: 50 points.
#' 
#' @param legend boolean value that determines if a colour ramp for the heatmap's
#' colour palette should be drawn or not.
#' Default is \code{TRUE}.
#' 
#' @param annCol specifications of column annotation tracks displayed as coloured 
#' rows on top of the heatmaps. The annotation tracks are drawn from bottom to top.
#' A single annotation track can be specified as a single vector; multiple tracks 
#' are specified as a list, a data frame, or an 
#' \code{\link[Biobase:ExpressionSet-class]{ExpressionSet}} object, in 
#' which case the phenotypic data is used (\code{pData(eset)}).
#' Character or integer vectors are converted and displayed as factors.
#' Unnamed tracks are internally renamed into \code{Xi}, with i being incremented for 
#' each unamed track, across both column and row annotation tracks.
#' For each track, if no corresponding colour is specified in argument 
#' \code{annColors}, a palette or a ramp is automatically computed and named 
#' after the track's name.
#' 
#' @param annRow specifications of row annotation tracks displayed as coloured 
#' columns on the left of the heatmaps. The annotation tracks are drawn from 
#' left to right. The same conversion, renaming and colouring rules as for argument 
#' \code{annCol} apply.
#' 
#' @param annColors list for specifying annotation track colors manually. It is 
#' possible to define the colors for only some of the annotations. Check examples for 
#' details.
#' 
#' @param annLegend boolean value specifying if the legend for the annotation tracks 
#' should be drawn or not.
#' Default is \code{TRUE}.
#' 
#' @param cexAnn scaling coefficent for the size of the annotation tracks.
#' Values > 1 (resp. < 1) will increase (resp. decrease) the size of each annotation 
#' track.
#' This applies to the height (resp. width) of the column (resp. row) annotation tracks.
#' Separate row and column sizes can be specified as a vector \code{c(row_size, col_size)}, 
#' where an NA value means using the default for the corresponding track.
#' 
#' @param labRow labels for the rows.
#' @param labCol labels for the columns. See description for argument \code{labRow} 
#' for a list of the possible values.
#' 
#' @param layout layout specification that indicates the relative position 
#' of the heatmap's components.
#' Two layouts can be defined: one horizontal, which relates to components associated to rows, 
#' and one vertical, which relates to components associated with columns.
#' Each layout is specified as a character strings, composed of characters 
#' that encode the order of each component: dendrogram (d), annotation tracks (a), 
#' data matrix (m), labels (l) and legend (L).
#' 
#' See \code{\link{aheatmap_layout}} for more details on layout specifications. 
#'  
#' @param fontsize base fontsize for the plot 
#' @param cexRow fontsize for the rownames, specified as a fraction of argument 
#' \code{fontsize}. 
#' @param cexCol fontsize for the colnames, specified as a fraction of argument 
#' \code{fontsize}.
#' 
#' @param main Main title as a character string or a grob.
#' @param sub Subtitle as a character string or a grob.
#' @param info (experimental) Extra information as a character vector or a grob.
#'  If \code{info=TRUE}, information about the clustering methods is displayed
#' at the bottom of the plot.
#' 
#' @param filename file path ending where to save the picture. Currently following 
#' formats are supported: png, pdf, tiff, bmp, jpeg. Even if the plot does not fit into 
#' the plotting window, the file size is calculated so that the plot would fit there, 
#' unless specified otherwise.
#' @param width manual option for determining the output file width in
#' @param height manual option for determining the output file height in inches.
#' 
#' @param verbose if \code{TRUE} then verbose messages are displayed and the 
#' borders of some viewports are highlighted. It is entended for debugging 
#' purposes.
#' 
#' @param gp graphical parameters for the text used in plot. Parameters passed to 
#' \code{\link{grid.text}}, see \code{\link{gpar}}. 
#' 
#' @author 
#' Original version of \code{pheatmap}: Raivo Kolde
#' 
#' Enhancement into \code{aheatmap}: Renaud Gaujoux
#' 
#' 
#' @section Row/column ordering and display: 
#' 
#' Possible values are:
#' \itemize{
#' \item \code{TRUE} or \code{NULL} (to be consistent with \code{\link{heatmap}}):
#' compute a dendrogram from hierarchical clustering using the distance and 
#' clustering methods \code{distfun} and \code{hclustfun}.
#' 
#' \item \code{NA}: disable any ordering. In this case, and if not otherwise 
#' specified with argument \code{revC=FALSE}, the heatmap shows the input matrix 
#' with the rows in their original order, with the first row on top to the last 
#' row at the bottom. Note that this differ from the behaviour or \code{\link{heatmap}}, 
#' but seemed to be a more sensible choice when vizualizing a matrix without 
#' reordering.
#' 
#' \item an integer vector of length the number of rows of the input matrix 
#' (\code{nrow(x)}), that specifies the row order. As in the case \code{Rowv=NA}, 
#' the ordered matrix is shown first row on top, last row at the bottom. 
#' 
#' \item a character vector or a list specifying values to use instead of arguments 
#' \code{distfun}, \code{hclustfun} and \code{reorderfun} when clustering the 
#' rows (see the respective argument descriptions for a list of accepted 
#' values). 
#' If \code{Rowv} has no names, then the first element is used for \code{distfun},  
#' the second (if present) is used for \code{hclustfun}, and the third 
#' (if present) is used for \code{reorderfun}.
#'
#' \item a numeric vector of weights, of length the number of rows of the input matrix, 
#' used to reorder the internally computed dendrogram \code{d} 
#' by \code{reorderfun(d, Rowv)}.
#' 
#' \item \code{FALSE}: the dendrogram \emph{is} computed using methods \code{distfun}, 
#' \code{hclustfun}, and \code{reorderfun} but is not shown.
#' 
#' \item a single integer that specifies how many subtrees (i.e. clusters) 
#' should be highlighted, e.g., \code{aheatmap(x, Rowv = 3L)}.
#' 
#' If positive, then the dendrogram's branches upstream each cluster are faded out 
#' using dashed lines.
#' If negative, then the dendrogram's branches within each cluster are faded out 
#' using dashed lines, keeping the root upstream branches as is.
#' 
#' \item a single double that specifies how much space is used by the computed 
#' dendrogram. That is that this value is used in place of \code{treeheight}.
#' 
#' \item a single character string starting with a \code{'#'} or a list with its 
#' first element as such a string, e.g., \code{aheatmap(x, Rowv = '#3')} or 
#' \code{aheatmap(x, Colv = list('#3', text = LETTERS[1:3]))}.
#' 
#' 
#' 
#' }
#' 
#' @examples
#' 
#' ## See the demo 'aheatmap' for more examples:
#' \dontrun{
#' demo('aheatmap')
#' }
#' 
#' # Generate random data
#' n <- 50; p <- 20
#' x <- abs(rmatrix(n, p, rnorm, mean=4, sd=1))
#' x[1:10, seq(1, 10, 2)] <- x[1:10, seq(1, 10, 2)] + 3
#' x[11:20, seq(2, 10, 2)] <- x[11:20, seq(2, 10, 2)] + 2
#' rownames(x) <- paste("ROW", 1:n)
#' colnames(x) <- paste("COL", 1:p)
#'
#' ## Default heatmap
#' aheatmap(x)
#' 
#' ## Distance methods
#' aheatmap(x, Rowv = "correlation")
#' aheatmap(x, Rowv = "man") # partially matched to 'manhattan'
#' aheatmap(x, Rowv = "man", Colv="binary")
#' 
#' # Generate column annotations
#' annotation = data.frame(Var1 = factor(1:p %% 2 == 0, labels = c("Class1", "Class2")), Var2 = 1:10)
#' aheatmap(x, annCol = annotation)
#' 
#' @demo Annotated heatmaps 
#' 
#' # Generate random data
#' n <- 50; p <- 20
#' x <- abs(rmatrix(n, p, rnorm, mean=4, sd=1))
#' x[1:10, seq(1, 10, 2)] <- x[1:10, seq(1, 10, 2)] + 3
#' x[11:20, seq(2, 10, 2)] <- x[11:20, seq(2, 10, 2)] + 2
#' rownames(x) <- paste("ROW", 1:n)
#' colnames(x) <- paste("COL", 1:p)
#' 
#' ## Scaling
#' aheatmap(x, scale = "row")
#' aheatmap(x, scale = "col") # partially matched to 'column'
#' aheatmap(x, scale = "r1") # each row sum up to 1 
#' aheatmap(x, scale = "c1") # each colum sum up to 1
#' 
#' ## Heatmap colors
#' aheatmap(x, color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
#' # color specification as an integer: use R basic colors
#' aheatmap(x, color = 1L)
#' # color specification as a negative integer: use reverse basic palette
#' aheatmap(x, color = -1L)
#' # color specification as a numeric: use HCL color
#' aheatmap(x, color = 1)
#' # color for NA values
#' y <- x
#' y[sample(length(y), p)] <- NA
#' aheatmap(y)
#' aheatmap(y, na.color = 'black')
#' 
#' # do not cluster the rows 
#' aheatmap(x, Rowv = NA)
#' # no heatmap legend
#' aheatmap(x, legend = FALSE)
#' # cell and font size 
#' aheatmap(x, cellwidth = 10, cellheight = 5)
#' 
#' # directly write into a file
#' aheatmap(x, cellwidth = 15, cellheight = 12, fontsize = 8, filename = "aheatmap.pdf")
#' unlink('aheatmap.pdf')
#'
#' # Generate column annotations
#' annotation = data.frame(Var1 = factor(1:p %% 2 == 0, labels = c("Class1", "Class2")), Var2 = 1:10)
#'
#' aheatmap(x, annCol = annotation)
#' aheatmap(x, annCol = annotation, annLegend = FALSE)
#'
#'
#' # Specify colors
#' Var1 = c("navy", "darkgreen")
#' names(Var1) = c("Class1", "Class2")
#' Var2 = c("lightgreen", "navy")
#'
#' ann_colors = list(Var1 = Var1, Var2 = Var2)
#'
#' aheatmap(x, annCol = annotation, annColors = ann_colors)
#' 
#' # Specifying clustering from distance matrix
#' drows = dist(x, method = "minkowski")
#' dcols = dist(t(x), method = "minkowski")
#' aheatmap(x, Rowv = drows, Colv = dcols)
#' 
#' # Display text in each cells
#' t <- outer(as.character(outer(letters, letters, paste0)), letters, paste0)[1:n, 1:p] 
#' aheatmap(x, txt = t)
#' # NA values are shown as empty cells
#' t.na <- t
#' t.na[sample(length(t.na), 500)] <- NA # half of the cells
#' aheatmap(x, txt = t.na)
#' 
#' @export
aheatmap = function(x
, color = '-RdYlBu2:100', na.color = NA
, breaks = NA, border_color=NA, cellwidth = NA, cellheight = NA, scale = "none"
, Rowv=TRUE, Colv=TRUE
, revC = identical(Colv, "Rowv") || is_NA(Rowv) || (is.integer(Rowv) && length(Rowv) > 1)
    || is(Rowv, 'silhouette')
, distfun = "euclidean", hclustfun = "complete", reorderfun = function(d,w) reorder(d,w)
, treeheight = 50
, legend = TRUE, annCol = NA, annRow = NA, annColors = NA, annLegend = TRUE, cexAnn = NA
, labRow = NULL, labCol = NULL
, subsetRow = NULL, subsetCol = NULL
, txt = NULL, layout = '.'
, fontsize=10, cexRow = min(0.2 + 1/log10(nr), 1.2), cexCol = min(0.2 + 1/log10(nc), 1.2)
, filename = NA, width = NA, height = NA
, main = NULL, sub = NULL, info = NULL
, verbose=getOption('verbose'), gp = gpar()){

	# set verbosity level
	ol <- lverbose(verbose)
	on.exit( lverbose(ol) )
	
	# convert ExpressionSet into 
	if( is(x, 'ExpressionSet') ){
		library(Biobase)
		if( isTRUE(annCol) ) annCol <- atrack(x)
		x <- Biobase::exprs(x)
	}

	# rename to old parameter name
	mat <- x
    if( !is.null(txt) ){
        if( !all(dim(mat), dim(x)) ){
            stop("Incompatible data and text dimensions: arguments x and txt must have the same size.")
        }
    }
	
	# init result list
	res <- list()
	
	# treeheight: use common or separate spec for rows and columns 
	if( length(treeheight) == 1 )
		treeheight <- c(treeheight, treeheight)
	treeheight_row <- treeheight[1]
	treeheight_col <- treeheight[2]
	
	## SUBSET: process subset argument for rows/columsn if requested.
	# this has to be done before relabelling and clustering
	# but the actual subsetting has to be done after relabelling and before 
	# clustering.
	# Here one convert a subset argument into an interger vector with the indexes
	if( !is.null(subsetRow) ){
		if( verbose ) message("Compute row subset indexes")
		subsetRow <- subset_index(mat, 1L, subsetRow)
	}
	if( !is.null(subsetCol) ){
		if( verbose ) message("Compute column subset indexes")
		subsetCol <- subset_index(mat, 2L, subsetCol)
	}
	
	## LABELS: set the row/column labels
	# label row numerically if no rownames
	if( is.null(labRow) && is.null(rownames(mat)) )
		labRow <- 1L
	if( !is.null(labRow) ){
		if( verbose ) message("Process labRow")
		rownames(mat)  <- generate_dimnames(labRow, nrow(mat), rownames(mat))	
	}
	# label columns numerically if no colnames
	if( is.null(labCol) && is.null(colnames(mat)) )
		labCol <- 1L
	if( !is.null(labCol) ){
		if( verbose ) message("Process labCol")
		colnames(mat) <- generate_dimnames(labCol, ncol(mat), colnames(mat))
	}
	
	## DO SUBSET
	if( !is.null(subsetRow) ){
		mat <- mat[subsetRow, ]
	}	
	if( !is.null(subsetCol) ){
		mat <- mat[, subsetCol]
	}
	
	## CLUSTERING
	# Do row clustering	
	tree_row <- if( !is_NA(Rowv) ){
		if( verbose ) message("Cluster rows")
		# single numeric Rowv means treeheight
		if( isReal(Rowv) ){
			# treeheight
			treeheight_row <- Rowv			
			# do cluster the rows
			Rowv <- TRUE
		}
		cluster_mat(mat, Rowv
						, distfun=distfun, hclustfun=hclustfun
						, reorderfun=reorderfun, subset=subsetRow
						, verbose = verbose)		
	}
	else NA
	
	# do not show the tree if Rowv=FALSE or not a tree
	if( identical(Rowv, FALSE) || !is_treedef(tree_row) )
		treeheight_row <- 0
	
	# Do col clustering
	tree_col <- if( !is_NA(Colv) ){
		if( identical(Colv,"Rowv") ){ # use row indexing if requested
			if( ncol(mat) != nrow(mat) )
				stop("aheatmap - Colv='Rowv' but cannot treat columns and rows symmetrically: input matrix is not square.")
			treeheight_col <- treeheight_row 
			tree_row
		}else{
			# single numeric Colv means treeheight
			if( isReal(Colv) ){
				# tree height
				treeheight_col <- Colv				
				# do cluster the columns
				Colv <- TRUE
			}
			if( verbose ) message("Cluster columns")
			cluster_mat(t(mat), Colv
							, distfun=distfun, hclustfun=hclustfun
							, reorderfun=reorderfun, subset=subsetCol
							, verbose = verbose)
		}		
	}
	else NA
		
	# do not show the tree if Colv=FALSE
	if( identical(Colv, FALSE) || !is_treedef(tree_col) )
		treeheight_col <- 0
	
	## ORDER THE DATA
	if( !is_NA(tree_row) ){
		
		# revert the row order if requested
		if( revC ){
			if( verbose ) message("Reverse row clustering")
			tree_row <- rev(tree_row)
		}
		
		# store the order and row tree if possible
		if( is_treedef(tree_row) ){
			res$Rowv <- tree_row$dendrogram
			res$rowInd <- order.dendrogram(tree_row$dendrogram)
			if( length(res$rowInd) != nrow(mat) )
				stop("aheatmap - row dendrogram ordering gave index of wrong length (", length(res$rowInd), ")")
		}
		else{			
			res$rowInd <- tree_row
			tree_row <- NA
		}
		
	}else if( revC ){ # revert the row order if requested
		res$rowInd <- nrow(mat):1L			
	}
	# possibly map the index to the original data index 
	res$rowInd <- subset2orginal_idx(res$rowInd, subsetRow)
	
	# order the rows if necessary
	if( !is.null(res$rowInd) ){
		# check validity of ordering
		if( !is.integer(res$rowInd) || length(res$rowInd) != nrow(mat) )
			stop("aheatmap - Invalid row ordering: should be an integer vector of length nrow(mat)=", nrow(mat))
		
		if( verbose ) message("Order rows")
		subInd <- attr(res$rowInd, 'subset')
        ri <- if( is.null(subInd) ) res$rowInd else subInd
		mat <- mat[ri, , drop=FALSE] # data
        if( !is.null(txt) ) txt <- txt[ri, , drop = FALSE] # text 
	}
	
	if( !is_NA(tree_col) ){		
		# store the column order and tree if possible
		if( is_treedef(tree_col) ){
			res$Colv <- tree_col$dendrogram
			res$colInd <- order.dendrogram(tree_col$dendrogram)
			if( length(res$colInd) != ncol(mat) )
				stop("aheatmap - column dendrogram ordering gave index of wrong length (", length(res$colInd), ")")
		}else{
			res$colInd <- tree_col
			tree_col <- NA
		}
	}
	# possibly map the index to the original data index 
	res$colInd <- subset2orginal_idx(res$colInd, subsetCol)
	
	# order the columns if necessary
	if( !is.null(res$colInd) ){
		# check validity of ordering
		if( !is.integer(res$colInd) || length(res$colInd) != ncol(mat) )
			stop("aheatmap - Invalid column ordering: should be an integer vector of length ncol(mat)=", ncol(mat))
		
		if( verbose ) message("Order columns")
		subInd <- attr(res$colInd, 'subset')
        ci <- if( is.null(subInd) ) res$colInd else subInd
		mat <- mat[, ci, drop=FALSE] # data
        if( !is.null(txt) ) txt <- txt[, ci, drop = FALSE] # text
	}
	
	# adding clustering info
	if( isTRUE(info) || is.character(info) ){
		
		if( verbose ) message("Compute info")
		if( !is.character(info) ) info <- NULL
		linfo <- NULL
		if( is_treedef(tree_row) && !is.null(tree_row$dist.method) )
			linfo <- paste("rows:", tree_row$dist.method, '/', tree_row$method)
		if( is_treedef(tree_col) && !is.null(tree_col$dist.method) )
			linfo <- paste(linfo, paste(" - cols:", tree_col$dist.method, '/', tree_col$method))
		info <- c(info, linfo)
	}
	
	# drop extra info except dendrograms for trees
	if( is_treedef(tree_col) )
		tree_col <- tree_col$dendrogram
	if( is_treedef(tree_row) )
		tree_row <- tree_row$dendrogram
	
	# Preprocess matrix
	if( verbose ) message("Scale matrix")
	mat = as.matrix(mat)
	mat = scale_mat(mat, scale)
	
	## Colors and scales
    # generate colour scale
    if( verbose ) message("Generate colour scale (palette + breaks)")
    if( is_NA(breaks) ) breaks <- NULL
	colour_scale <- ccRamp(color, breaks = breaks, data = as.vector(mat))
    # store into result list
    res$col <- colour_scale
    # assign as separate objects (legacy)
	breaks <- setNames(colour_scale, NULL)
    color <- names(colour_scale)
        
    if( isTRUE(legend) ){
		if( verbose ) message("Generate colour scale ticks")
		legend = grid.pretty(range(as.vector(breaks)))
	}
	else {
		legend = NA
	}
    
    # convert numeric values to colours
    if( verbose ) message("Map values to colours")
    if( is.list(na.color) && length(na.color) > 1L ){ # force specified values to NA
        na_range <- na.color[[2L]]
        if( length(na_range) == 1L ) mat[mat %in% na_range] <- NA
        else mat[mat >= na_range[1L] & mat <= na_range[2L] ] <- NA
    }
    mat = scale_colours(mat, col = color, breaks = breaks)
    if( !is_NA(na.color) ){ # use specified color for NA values
        mat[is.na(mat)] <- na.color[[1L]]
    }
    
	annotation_legend <- annLegend
	annotation_colors <- annColors
	
	# render annotation tracks for both rows and columns
	annCol_processed <- atrack(annCol, order=res$colInd, .SPECIAL=specialAnnotation(2L), .DATA=amargin(x,2L), .CACHE=annRow)
	annRow_processed <- atrack(annRow, order=res$rowInd, .SPECIAL=specialAnnotation(1L), .DATA=amargin(x,1L), .CACHE=annCol)
	specialAnnotation(clear=TRUE)
	annTracks <- renderAnnotations(annCol_processed, annRow_processed 
								, annotation_colors = annotation_colors
								, verbose=verbose)
	#
	
	# retrieve dimension for computing cexRow and cexCol (evaluated from the arguments)
	nr <- nrow(mat); nc <- ncol(mat)
	# Draw heatmap	
	res$vp <- heatmap_motor(mat, border_color = border_color, cellwidth = cellwidth, cellheight = cellheight
	, treeheight_col = treeheight_col, treeheight_row = treeheight_row, tree_col = tree_col, tree_row = tree_row
	, filename = filename, width = width, height = height, breaks = breaks, color = color, legend = legend
	, annTracks = annTracks, annotation_legend = annotation_legend, cexAnn = cexAnn
    , txt = txt
	, fontsize = fontsize, fontsize_row = cexRow * fontsize, fontsize_col = cexCol * fontsize
	, main = main, sub = sub, info = info, layout = layout
	, verbose = verbose
	, gp = gp)
	
	# return info about the plot
	invisible(res)
}

#' @import gridBase
grid.base.mix <- function(opar, trace = getOption('verbose')){
	
	if( !missing(opar) ){
		if( !is.null(opar) ){
			if( trace ) message("grid.base.mix - restore")			
			upViewport(2)
			par(opar)			
		}
		return(invisible())
	}
	
	if( trace ) message("grid.base.mix - init")
	if( trace ) grid.rect(gp=gpar(lwd=40, col="blue"))
	opar <- par(xpd=NA)
	if( trace ) grid.rect(gp=gpar(lwd=30, col="green"))
	if( trace ) message("grid.base.mix - plot.new")
	plot.new()
	if( trace ) grid.rect(gp=gpar(lwd=20, col="black"))
	vps <- baseViewports()	
	pushViewport(vps$inner)
	if( trace ) grid.rect(gp=gpar(lwd=10, col="red"))
	pushViewport(vps$figure)
#	if( trace ) grid.rect(gp=gpar(lwd=3, col="green"))
#	pushViewport(vps$plot)
#	upViewport(2)
#	if( trace ) grid.rect(gp=gpar(lwd=3, col="pink"))
#	pushViewport(viewport(x=unit(0.5, "npc"), y=unit(0, "npc"), width=unit(0.5, "npc"), height=unit(1, "npc"), just=c("left","bottom")))
	if( trace )	grid.rect(gp=gpar(lwd=3, col="yellow"))
	opar
}

if( FALSE ){
	testmix <- function(){
		opar <- mixplot.start(FALSE)
		profplot(curated$data$model, curated$fits[[1]])
		mixplot.add(TRUE)
		basismarker(curated$fits[[1]], curated$data$markers)
		mixplot.end()
		par(opar)
	}
	
	dd <- function(d, horizontal = TRUE, ...){
		grid.rect(gp = gpar(col = "blue", lwd = 2))
		opar <- par(plt = gridPLT(), new = TRUE)
		on.exit(par(opar))
		plot(d, horiz=horizontal, xaxs="i", yaxs="i", axes=FALSE, leaflab="none", ...)
	}
	
	toto <- function(new=FALSE){
		
		library(RGraphics)
		set.seed(123456)
		x <- matrix(runif(30*20), 30)
		x <- crossprod(x)
		d <- as.dendrogram(hclust(dist(x)))
		
		#grid.newpage()
		if( new ) plot.new()
		lo <- grid.layout(nrow=2, ncol=2)
		pushViewport(viewport(layout=lo))
		
		pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
		dd(d)
		upViewport()
		pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
		dd(d, FALSE)
		upViewport()
		pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
		grid.imageGrob(nrow(x), ncol(x), x)
		upViewport()	
		
		popViewport()
		stop("END toto")
		
	} 
	
	test <- function(){
		pdf('aheatmap.pdf')
		#try(v <- aheatmap(consensus(res), color='grey:100', Colv=2L, verbose=TRUE))
		try(v <- consensusmap(res, color='grey:100', Colv=2L, verbose=TRUE))
		dev.off()
	}
	
	test2 <- function(){
		op <- par(mfrow=c(1,2))
		on.exit(par(op))
		#try(v <- aheatmap(consensus(res), color='grey:100', Colv=2L, verbose=TRUE))
		try(v <- consensusmap(res, verbose=TRUE))
		try(v <- consensusmap(res, color='grey:100', Colv=2L, verbose=TRUE))	
	}
	
	testsw <- function(file=TRUE){
		if(file ){
			pdf('asweave.pdf', width=20, height=7)
			on.exit(dev.off())
		}
		opar <- par(mfrow=c(1,2))
# removing all automatic annotation tracks
		coefmap(res, tracks=NA, verbose=TRUE)
# customized plot	
		coefmap(res, Colv = 'euclidean', Rowv='max', verbose=TRUE)
#			, main = "Metagene contributions in each sample", labCol = NULL
#			, tracks = c(Metagene='basis'), annCol = list(Class=a, Index=c)
#			, annColors = list(Metagene='Set2')
#			, info = TRUE)
		par(opar)	
	}
	
	testvp <- function(file=TRUE){
		if(file ){
			pdf('avp.pdf', width=20, height=7)
			on.exit(dev.off())
		}
		
		plot.new()
		lo <- grid.layout(nrow=1, ncol=2)
		pushViewport(viewport(layout=lo))
		
		pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
		basismap(res, Colv='eucl', verbose=TRUE)
		upViewport()
		pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
		coefmap(res, tracks=NA, verbose=TRUE)
		upViewport()
		
		popViewport()	
		
	} 
	
}
