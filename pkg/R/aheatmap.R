#' @include atracks.R
#' @include colorcode.R
NULL

library(grid)
library(colorspace)
library(stringr)

lo <- function (rown, coln, nrow, ncol, cellheight = NA, cellwidth = NA
, treeheight_col, treeheight_row, legend, main = NULL, sub = NULL, info = NULL
, annTracks, annotation_legend
, fontsize, fontsize_row, fontsize_col, ...){

	annotation_colors <- annTracks$colors
	row_annotation <- annTracks$annRow
	annotation <- annTracks$annCol
	
	coln_height <- unit(10, "bigpts")
	if(!is.null(coln)){
		longest_coln = which.max(nchar(coln))
		gp = gpar(fontsize = fontsize_col, ...)
		coln_height <- coln_height +  unit(1.1, "grobheight", textGrob(coln[longest_coln], rot = 90, gp = gp))
	}

	rown_width <- rown_width_min <- unit(10, "bigpts")
	if(!is.null(rown)){
		longest_rown = which.max(nchar(rown))
		gp = gpar(fontsize = fontsize_row, ...)
		rown_width <- rown_width_min + unit(1.2, "grobwidth", textGrob(rown[longest_rown], gp = gp))
	}
	
	gp = list(fontsize = fontsize, ...)
	# Legend position
	if( !isNA(legend) ){
		longest_break = which.max(nchar(as.character(legend)))
		longest_break = unit(1.1, "grobwidth", textGrob(as.character(legend)[longest_break], gp = do.call(gpar, gp)))
		# minimum fixed width: plan for 2 decimals and a sign 
		min_lw = unit(1.1, "grobwidth", textGrob("-00.00", gp = do.call(gpar, gp)))
		longest_break = max(longest_break, min_lw)
		title_length = unit(1.1, "grobwidth", textGrob("Scale", gp = gpar(fontface = "bold", ...)))
		legend_width = unit(12, "bigpts") + longest_break * 1.2
		legend_width = max(title_length, legend_width)
	}
	else{
		legend_width = unit(0, "bigpts")
	}
	
	.annLegend.dim <- function(annotation, fontsize){
		# Width of the corresponding legend
		longest_ann <- unlist(lapply(annotation, names))
		longest_ann <- longest_ann[which.max(nchar(longest_ann))]
		annot_legend_width = unit(1, "grobwidth", textGrob(longest_ann, gp = gpar(fontsize=fontsize, ...))) + unit(10, "bigpts")
		
		# width of the legend title
		annot_legend_title <- names(annotation)[which.max(nchar(names(annotation)))]
		annot_legend_title_width = unit(1, "grobwidth", textGrob(annot_legend_title, gp = gpar(fontface = "bold", fontsize=fontsize, ...)))
		
		# total width 
		max(annot_legend_width, annot_legend_title_width) + unit(5, "bigpts")
	}
	
	# Column annotations
	if( !isNA(annotation) ){
		# Column annotation height		
		annot_height = unit(ncol(annotation) * (8 + 2) + 2, "bigpts")
	}
	else{
		annot_height = unit(0, "bigpts")
	}
	
	# add a viewport for the row annotations
	if ( !isNA(row_annotation) ) {
		# Row annotation width		
		row_annot_width = unit(ncol(row_annotation) * (8 + 2) + 2, "bigpts")
	}
	else {
		row_annot_width = unit(0, "bigpts")
	}
	
	# Width of the annotation legend
	annot_legend_width <- if( annotation_legend ) .annLegend.dim(annotation_colors, fontsize)
			else unit(0, "npc")

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
		matwidth = unit(1, "npc") - rown_width - legend_width - row_annot_width  - treeheight_row - annot_legend_width
	}
	else{
		matwidth = unit(cellwidth * ncol, "bigpts")
	}

	if(is.na(cellheight)){
		matheight = unit(1, "npc") - treeheight_col - annot_height - main_height - coln_height - sub_height - info_height
	
		# recompute the cell width depending on the automatic fontsize
		if( is.na(cellwidth) && !is.null(rown) ){
			cellheight <- convertHeight(unit(1, "grobheight", rectGrob(0,0, matwidth, matheight)), "bigpts", valueOnly = T) / nrow
			fontsize_row <- convertUnit(min(unit(fontsize_row, 'points'), unit(0.6*cellheight, 'bigpts')), 'points')
			
			rown_width <- rown_width_min + unit(1.2, "grobwidth", textGrob(rown[longest_rown], gp = gpar(fontsize=fontsize_row, ...)))
			matwidth <- unit(1, "npc") - rown_width - legend_width - row_annot_width  - treeheight_row - annot_legend_width
		}
	}
	else{
		matheight = unit(cellheight * nrow, "bigpts")
	}	
		
	# HACK: 
	# - use 6 instead of 5 column for the row_annotation
	# - take into account the associated legend's width
	# Produce layout()
	unique.name <- vplayout(NULL)
	lo <- grid.layout(nrow = 7, ncol = 6
			, widths = unit.c(treeheight_row, row_annot_width, matwidth, rown_width, legend_width, annot_legend_width)
			, heights = unit.c(main_height, treeheight_col,  annot_height, matheight, coln_height, sub_height, info_height))
	hvp <- viewport( name=paste('aheatmap', unique.name, sep='-'), layout = lo)
	pushViewport(hvp)
	
	#grid.show.layout(lo); stop('sas')
	# Get cell dimensions
	vplayout('mat')
	cellwidth = convertWidth(unit(1, "npc"), "bigpts", valueOnly = T) / ncol
	cellheight = convertHeight(unit(1, "npc"), "bigpts", valueOnly = T) / nrow
	upViewport()
	
	height <- as.numeric(convertHeight(sum(lo$height), "inches"))
	width <- as.numeric(convertWidth(sum(lo$width), "inches"))
	# Return minimal cell dimension in bigpts to decide if borders are drawn
	mindim = min(cellwidth, cellheight) 
	return( list(width=width, height=height, vp=hvp, mindim=mindim, cellwidth=cellwidth, cellheight=cellheight) )
}

draw_dendrogram = function(hc, horizontal = T){
	
#	.draw.dendrodram <- function(hc){
#		
#		# convert into an hclust if necessary
#		if( is(hc, 'dendrogram') ){
#			hca <- attr(hc, 'hclust')
#			hc <- if( !is.null(hca) ) hca else as.hclust(hc)
#		}
#		
#		h = hc$height / max(hc$height) / 1.05
#		m = hc$merge
#		o = hc$order
#		n = length(o)
#	
#		m[m > 0] = n + m[m > 0] 
#		m[m < 0] = abs(m[m < 0])
#	
#		dist = matrix(0, nrow = 2 * n - 1, ncol = 2, dimnames = list(NULL, c("x", "y"))) 
#		dist[1:n, 1] = 1 / n / 2 + (1 / n) * (match(1:n, o) - 1)
#	
#		for(i in 1:nrow(m)){
#			dist[n + i, 1] = (dist[m[i, 1], 1] + dist[m[i, 2], 1]) / 2
#			dist[n + i, 2] = h[i]
#		}
#	
#		draw_connection = function(x1, x2, y1, y2, y){
#			grid.lines(x = c(x1, x1), y = c(y1, y))
#			grid.lines(x = c(x2, x2), y = c(y2, y))
#			grid.lines(x = c(x1, x2), y = c(y, y))
#		}
#		
#		# create a rotating viewport for vertical dendrogram 
#		if(!horizontal){
#			gr = rectGrob()
#			pushViewport(viewport(height = unit(1, "grobwidth", gr), width = unit(1, "grobheight", gr), angle = 90))
#			on.exit(upViewport())
#		}
#		
#		for(i in 1:nrow(m)){
#			draw_connection(dist[m[i, 1], 1], dist[m[i, 2], 1], dist[m[i, 1], 2], dist[m[i, 2], 2], h[i])
#		}		
#				
#	}
	
	.draw.dendrodram <- function(hc, ...){
		suppressWarnings( opar <- par(plt = gridPLT(), new = TRUE) )
		on.exit(par(opar))
		if( getOption('verbose') ) grid.rect(gp = gpar(col = "blue", lwd = 2))
		if( !is(hc, 'dendrogram') )
			hc <- as.dendrogram(hc)
		plot(hc, horiz=!horizontal, xaxs="i", yaxs="i", axes=FALSE, leaflab="none", ...)
	}
	
		
	# create a margin viewport
	if(!horizontal)
		pushViewport( viewport(x=0,y=0,width=0.9,height=1,just=c("left", "bottom")) )
	else
		pushViewport( viewport(x=0,y=0.1,width=1,height=0.9,just=c("left", "bottom")) )
	on.exit(upViewport())
	
	.draw.dendrodram(hc)

}

# draw a matrix first row at bottom, last at top
draw_matrix = function(matrix, border_color){
	n = nrow(matrix)
	m = ncol(matrix)
	x = (1:m)/m - 1/2/m
	y = (1:n)/n - 1/2/n
	for(i in 1:m){
		grid.rect(x = x[i], y = y, width = 1/m, height = 1/n, gp = gpar(fill = matrix[,i], col = border_color))
	}
}

draw_colnames = function(coln, ...){
	m = length(coln)
	x = (1:m)/m - 1/2/m
	grid.text(coln, x = x, y = unit(1, 'npc') - unit(5, "bigpts"), vjust = 0.5, hjust = 0, rot = 270, gp = gpar(...))
}

# draw rownames first row at bottom, last on top
draw_rownames = function(rown, ...){
	n = length(rown)
	y = (1:n)/n - 1/2/n
	grid.text(rown, x = unit(5, "bigpts"), y = y, vjust = 0.5, hjust = 0, gp = gpar(...))	
}

draw_legend = function(color, breaks, legend, ...){
	height = min(unit(1, "npc"), unit(150, "bigpts"))
	pushViewport(viewport(x = 0, y = unit(1, "npc"), just = c(0, 1), height = height))
	legend_pos = (legend - min(breaks)) / (max(breaks) - min(breaks))
	breaks = (breaks - min(breaks)) / (max(breaks) - min(breaks))
	h = breaks[-1] - breaks[-length(breaks)]
	grid.rect(x = 0, y = breaks[-length(breaks)], width = unit(10, "bigpts"), height = h, hjust = 0, vjust = 0, gp = gpar(fill = color, col = "#FFFFFF00"))
	grid.text(legend, x = unit(12, "bigpts"), y = legend_pos, hjust = 0, gp = gpar(...))
	upViewport()
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

	colnames(x) <- names(annotation)
	return(x)
	#return(as.matrix(new))
}

draw_annotations = function(converted_annotations, border_color, horizontal=TRUE){
	n = ncol(converted_annotations)
	m = nrow(converted_annotations)
	if( horizontal ){
		x = (1:m)/m - 1/2/m
		y = cumsum(rep(8, n)) - 4 + cumsum(rep(2, n))
		for(i in 1:m){
			grid.rect(x = x[i], unit(y[1:n], "bigpts"), width = 1/m, height = unit(8, "bigpts"), gp = gpar(fill = converted_annotations[i, ], col = border_color))
		}
	}else{
		x = cumsum(rep(8, n)) - 4 + cumsum(rep(2, n))
		y = (1:m)/m - 1/2/m
		for (i in 1:m) {
			grid.rect(x = unit(x[1:n], "bigpts"), y=y[i], width = unit(8, "bigpts"), 
					height = 1/m, gp = gpar(fill = converted_annotations[i,]
					, col = border_color))
		}
	}
}

draw_annotation_legend = function(annotation_colors, border_color, ...){
	
	y = unit(1, "npc")
	
	text_height = convertHeight(unit(1, "grobheight", textGrob("FGH", gp = gpar(...))), "bigpts")	
	for(i in names(annotation_colors)){
		grid.text(i, x = 0, y = y, vjust = 1, hjust = 0, gp = gpar(fontface = "bold", ...))
		y = y - 1.5 * text_height
		#if(class(annotation[[i]]) %in% c("character", "factor")){
		acol <- annotation_colors[[i]]
		if( attr(acol, 'afactor') ){
			sapply(seq_along(acol), function(j){
				grid.rect(x = unit(0, "npc"), y = y, hjust = 0, vjust = 1, height = text_height, width = text_height, gp = gpar(col = border_color, fill = acol[j]))
				grid.text(names(acol)[j], x = text_height * 1.3, y = y, hjust = 0, vjust = 1, gp = gpar(...))
				y <<- y - 1.5 * text_height
			})
		}
		else{
			yy = y - 4 * text_height + seq(0, 1, 0.01) * 4 * text_height
			h = 4 * text_height * 0.02
			grid.rect(x = unit(0, "npc"), y = yy, hjust = 0, vjust = 1, height = h, width = text_height, gp = gpar(col = "#FFFFFF00", fill = ccRamp(acol, 100)))
			txt = c(tail(names(acol),1), head(names(acol))[1])
			yy = y - c(0, 3) * text_height
			grid.text(txt, x = text_height * 1.3, y = yy, hjust = 0, vjust = 1, gp = gpar(...))
			y = y - 4.5 * text_height
		}
		y = y - 1.5 * text_height
	}
}

tryViewport <- function(name){
	vpo <- current.viewport()
	vp <- try( seekViewport(name), silent=TRUE )
	
	res <- if( is(vp, 'try-error') ){
		seekViewport(vpo$name)
		NULL
	}else vp
	
	invisible(res)
}

vplayout <- function ()
{
	graphic.name <- NULL
	
	function(x, y, verbose = getOption('verbose') ){
		# initialize the graph name
		if( is.null(x) ){
			graphic.name <<- grid:::vpAutoName()
			return(graphic.name)
		}
		name <- NULL
		if( !is.numeric(x) ){
					
			name <- paste(graphic.name, x, sep='-')
			
			if( !missing(y) && is(y, 'viewport') ){
				y$name <- name
				return(pushViewport(y))			
			}
			if( verbose ) message("vp - lookup for ", name)
			if( !is.null(tryViewport(name)) )
				return()
			
			switch(x
				, main={x<-1; y<-3;}
				, ctree={x<-2; y<-3;}
				, cann={x<-3; y<-3;}
				, rtree={x<-4; y<-1;}
				, rann={x<-4; y<-2;}
				, mat={x<-4; y<-3;}
				, rnam={x<-4; y<-4;}
				, leg={x<-4; y<-5;}
				, aleg={x<-4; y<-6;}
				, cnam={x<-5; y<-3;}
				, sub={x<-6; y<-3;}
				, info={x<-7; y<-3;}
				, stop("aheatmap - invalid viewport name")
			)
		}
		if( verbose ) message("vp - create ", name)
		pushViewport(viewport(layout.pos.row = x, layout.pos.col = y, name=name))
	}	
}
vplayout <- vplayout()

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
			stop("File type should be: pdf, png, bmp, jpg, tiff")
	)
	
	args <- c(list(filename), list(...))	
	if( !missing(width) ){
		args$width <- as.numeric(width)
		args$height <- as.numeric(height)
		if( !ending %in% c('pdf','svg') ){
			args$units <- "in"
			args$res <- 300
		}
	}	
	do.call('f', args)	
}

#' Custom Access to Viewport Paths
#' 
#' This function was defined to substitute current.vpPath because this latter 
#' calls grid.Call("L_currentViewport") which calls 'L_gridDirty' before 
#' effectively calling 'L_currentViewport'. 
#' This starts a new page (not sure why), which is not ok when writing to a 
#' device that records each page (e.g. PDF or SVG)
#' 
current.vpPath2 <- function(){
	names <- NULL
	pvp <- .Call('L_currentViewport', PACKAGE='grid')
	if( is.null(pvp) ) return(NULL)
	while (!grid:::rootVP(pvp)) {
		names <- c(names, pvp$name)
		pvp <- pvp$parent
	}
	if (!is.null(names)) 
		grid:::vpPathFromVector(rev(names))
	else names	
}

heatmap_motor = function(matrix, border_color, cellwidth, cellheight
	, tree_col, tree_row, treeheight_col, treeheight_row
	, filename=NA, width=NA, height=NA
	, breaks, color, legend
	, annTracks, annotation_legend=TRUE
	, new=TRUE, fontsize, fontsize_row, fontsize_col
	, main=NULL, sub=NULL, info=NULL
	, verbose=getOption('verbose')
	, ...){

	annotation_colors <- annTracks$colors
	row_annotation <- annTracks$annRow
	annotation <- annTracks$annCol
	writeToFile <- !is.na(filename)
	
	# open graphic device (dimensions will be changed after computation of the correct height)
	if( writeToFile ){
		gfile(filename)
		on.exit(dev.off())
	}
	 
	# identify the plotting context: base or grid
	#NB: use custom function current.vpPath2 instead of official 
	# grid::current.vpPath as this one creates a new page when called 
	# on a fresh graphic device	
	vpp <- current.vpPath2()
	if( is.null(vpp) ){ # we are at the root viewport
		if( verbose ) message("Detected path: [ROOT]")
		mf <- par('mfrow')		
		#print(mf)
		# if in in mfrow/layout context: setup fake-ROOT viewports with gridBase
		# and do not call plot.new as it is called in grid.base.mix. 
		new <- if( !identical(mf, c(1L,1L)) ){ 
			if( verbose ) message("Detected mfrow: ", mf[1], " - ", mf[2])
			opar <- grid.base.mix(trace=verbose>1)
			on.exit( grid.base.mix(opar) )
			FALSE
		}		
		else new
	}else{
		if( verbose ) message("Detected path: ", vpp)
		# if new is not specified: change the default behaviour by not calling 
		# plot.new so that drawing occurs in the current viewport
		if( missing(new) ){
			if( verbose ) message("Force no new plot")
			new <- FALSE
		}
	}
	# reset device if necessary or requested
	if( new ){
		if( verbose ) message("Call: plot.new")
		#grid.newpage()
		plot.new()
	}	
	
	# define grob for main 
	mainGrob <- if( !is.null(main) && !is.grob(main) ) textGrob(main, gp = gpar(fontsize = 1.2 * fontsize, fontface="bold", ...))
	subGrob <- if( !is.null(sub) && !is.grob(sub) ) textGrob(sub, gp = gpar(fontsize = 0.8 * fontsize, ...))
	infoGrob <- if( !is.null(info) && !is.grob(info) ){
#		infotxt <- paste(strwrap(paste(info, collapse=" | "), width=20), collapse="\n")
		grobTree(gList(rectGrob(gp = gpar(fill = "grey80"))
			,textGrob(paste(info, collapse=" | "), x=unit(5, 'bigpts'), y=0.5, just='left', gp = gpar(fontsize = 0.8 * fontsize, ...))))
	}
	
	# Set layout
	glo = lo(coln = colnames(matrix), rown = rownames(matrix), nrow = nrow(matrix), ncol = ncol(matrix)
	, cellwidth = cellwidth, cellheight = cellheight
	, treeheight_col = treeheight_col, treeheight_row = treeheight_row
	, legend = legend
	, annTracks = annTracks, annotation_legend = annotation_legend
	, fontsize = fontsize, fontsize_row = fontsize_row, fontsize_col = fontsize_col
	, main = mainGrob, sub = subGrob, info = infoGrob, ...)
	
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
		if( verbose ) grid.rect(width=unit(glo$width, 'inches'), height=unit(glo$height, 'inches'), gp=gpar(col='blue'))
	}

	#grid.show.layout(glo$layout); return()
	mindim <- glo$mindim
	# Omit border color if cell size is too small 
	if(mindim < 3) border_color = NA

	# Draw tree for the columns
	if (!isNA(tree_col) &&  treeheight_col != 0){
		#vplayout(1, 2)
		vplayout('ctree')
		draw_dendrogram(tree_col, horizontal = T)
		upViewport()
	}

	# Draw tree for the rows
	if(!isNA(tree_row) && treeheight_row !=0){
		#vplayout(3, 1)
		vplayout('rtree')
		draw_dendrogram(tree_row, horizontal = F)
		upViewport()
	}

	# Draw matrix
	#vplayout(3, 2)
	vplayout('mat')
	draw_matrix(matrix, border_color)
	#grid.imageGrob(nrow(matrix), ncol(matrix), matrix, byrow=FALSE, gp=gpar(col=border_color))
	upViewport()

	# Draw colnames
	if(length(colnames(matrix)) != 0){
		#vplayout(4, 2)
		vplayout('cnam')
		fontsize_col <- convertUnit(min(unit(fontsize_col, 'points'), unit(0.6*glo$cellwidth, 'bigpts')), 'points')
		draw_colnames(colnames(matrix), fontsize = fontsize_col, ...)
		upViewport()
	}
	
	# Draw rownames
	if(length(rownames(matrix)) != 0){
		#vplayout(3, 3)
		vplayout('rnam')
		fontsize_row <- convertUnit(min(unit(fontsize_row, 'points'), unit(0.6*glo$cellheight, 'bigpts')), 'points')
		draw_rownames(rownames(matrix), fontsize = fontsize_row, ...)
		upViewport()
	}

	# Draw annotation tracks
	if( !isNA(annotation) ){
		#vplayout(2, 2)
		vplayout('cann')
		draw_annotations(annotation, border_color)
		upViewport()
	}	
	
	# add row annotations if necessary	
	if ( !isNA(row_annotation) ) {
		vplayout('rann')
		draw_annotations(row_annotation, border_color, horizontal=FALSE)
		upViewport()
	}
	
	# Draw annotation legend
	if( annotation_legend ){
		#vplayout(3, 5)
		vplayout('aleg')
		draw_annotation_legend(annotation_colors, border_color, fontsize = fontsize, ...)
		upViewport()
	}

	# Draw legend
	if(!isNA(legend)){
		#vplayout(3, 4)
		vplayout('leg')
		draw_legend(color, breaks, legend, fontsize = fontsize, ...)
		upViewport()
	}

	# Draw main
	if(!is.null(mainGrob)){
		vplayout('main')
		grid.draw(mainGrob)
		upViewport()
	}
	
	# Draw subtitle
	if(!is.null(subGrob)){
		vplayout('sub')
		grid.draw(subGrob)
		upViewport()
	}
	
	# Draw info
	if(!is.null(infoGrob)){
		vplayout('info')
		grid.draw(infoGrob)
		upViewport()
	}
		
	# return current vp tree
	ct <- current.vpTree()
	
	popViewport()
	
	ct
}

generate_breaks = function(x, n, center=NA){
	if( missing(center) || isNA(center) )
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
	return(col[as.numeric(cut(x, breaks = breaks, include.lowest = T))])
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

#' Fade Out the Upper Branches from a Dendrogram
#' 
#' @param x a dendrogram
#' @param n the number of groups
#' 
#' @keywords internal
cutdendro <- function(x, n){
	
	# exit early if n <=1: nothing to do
	if( n <= 1 ) return(x)
		
	# add node digest ids to x
	x <- dendrapply(x, function(n){
			attr(n, 'id') <- digest(attributes(n))
			n
	})

	# cut x in n groups
	# find the height where to cut
	h <- cutheight(x, n)
	cfx <- cut(x, h)
	# get the ids of the upper nodes
	ids <- sapply(cfx$lower, function(sub) attr(sub, 'id'))		
	
	# highlight the upper branches with dot lines
	dts <- c(lty=2, lwd=1.2, col=8)
	a <- dendrapply(x, function(node){
		a <- attributes(node)
		if( a$id %in% ids || (!is.leaf(node) && any(c(attr(node[[1]], 'id'), attr(node[[2]], 'id')) %in% ids)) )
			attr(node, 'edgePar') <- dts
		node
	})	
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
	}

	# index vectors are honoured
	if( is.integer(param) && length(param) > 1 ){
		
		# subset if requested: reorder the subset indexes as in param
		if( !is.null(subset) )
			param <- order(match(subset, param))			
		
		param
	}else{ # will compute dendrogram (NB: mat was already subset before calling cluster_mat)
		
		param <- 
		if( is.integer(param) )
			param
		else if( is.null(param) || isLogical(param) ) # use default reordering by rowMeans
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
		d <- if( is.character(distfun) ){
			distfun <- distfun[1]
			
			av <- c("correlation", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")
			i <- pmatch(distfun, av)
			if( isNA(i) )			
				stop("aheatmap - Invalid dissimilarity method, must be one of: ", paste("'", av, "'", sep='', collapse=', '))
			
			distfun <- av[i]
			if( verbose ) message("Using distance method: ", distfun)
			if(distfun == "correlation"){ d <- dist(1 - cor(t(mat))); attr(d, 'method') <- 'correlation'; d }
			else dist(mat, method = distfun)
	
		}else if( is(distfun, "dist") ){
			if( verbose ) message("Using dist object: ", distfun)
			distfun
		}else if( is.function(distfun) )
			distfun(mat)
		else
			stop("aheatmap - Invalid dissimilarity function: must be a character string, an object of class 'dist', or a function")
	
		# do hierarchical clustering 
		hc <- if( is.character(hclustfun) ){
			
			av <- c('ward', 'single', 'complete', 'average', 'mcquitty', 'median', 'centroid')
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
		if( is.integer(param) )
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
	if( isNA(i) )
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

library(stringr)

round.pretty <- function(x, min=2){
		
	if( is.null(x) ) return(NULL)		
	d <- max(x) - min(x)
	n <- 0
	#message("d0=", d)
	while(d<1){
		#message("d=", d)
		d <- d*10
		n <- n+1
	}	
	#message("n=", n)
	round(x, max(min,n))
}

generate_annotation_colours = function(annotation, annotation_colors, seed=TRUE){
	if( isNA(annotation_colors) ){
		annotation_colors = list()
	}
	
	count = 0
	annotation2 <- list()
	anames <- names(annotation)
	sapply(seq_along(annotation), function(i){
		a <- annotation[[i]]
		if( class(annotation[[i]]) %in% c("character", "factor")){			
			# convert to character vector
			a <- if( is.factor(a) ) levels(a) else unique(a)
			count <<- count + nlevels(a)
			# merge if possible
			if( !is.null(anames) && anames[i]!='' )	
				annotation2[[anames[i]]] <<- unique(c(annotation2[[anames[i]]], a))
			else 
				annotation2 <<- c(annotation2, list(a))			
		}else
			annotation2 <<- a
	})
	annotation <- annotation2
	#str(annotation2)

	factor_colors = hcl(h = seq(1, 360, length.out = max(count+1,20)), 100, 70)	
		
	# get random seeds to restore/update on exit
	rs <- RNGscope()
	on.exit({
		# update local random seed on exit
		.Rd.seed$.Random.seed <- getRNG()
		# restore global random seed
		RNGscope(rs)
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
				setNames(rev(sequential_hcl(2, h, l = c(50, 95))), round.pretty(range(ann, na.rm=TRUE)))
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
generate_dimnames <- function(x, n){
	if( isNA(x) ) NULL
	else if( length(x) == n ) x
	else if( identical(x, 1) || identical(x, 1L) ) 1L:n
	else if( is.character(x) && length(x) == 1 ){
		#print(str_match_all(x, "^/(([^%]*)(%[in])?)+/$"))
		paste(x, 1L:n, sep='')
	}
	else stop("aheatmap - Invalid row/column label. Possible values are: NA, a vector of correct length, value 1 (or 1L) or single character string.")
}


.make_annotation <- function(x, ord=NULL){
	# convert into a data.frame if necessary
	if( !is.data.frame(x) ){
		x <- if( is(x, 'ExpressionSet') ) pData(x)
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
	if( isNA(annotation_colors) ) annotation_colors <- list() 
	nc <- length(annCol)
	nr <- length(annRow)
	flag <- function(x, f){ if( missing(f) ) attr(x, 'flag') else{ attr(x, 'flag') <- f; x} }
	if( !isNA(annCol) ) annotation <- c(annotation, sapply(as.list(annCol), flag, 'col', simplify=FALSE))
	if( !isNA(annRow) ) annotation <- c(annotation, sapply(as.list(annRow), flag, 'row', simplify=FALSE))
		
	if( length(annotation) == 0 ) return( list(annCol=NA, annRow=NA, colors=NA) )
	
	# generate the missing names
	n <- names(annotation)
	xnames <- paste('X', 1:length(annotation), sep='')
	if( is.null(n) ) names(annotation) <- xnames
	else names(annotation)[n==''] <- xnames[n==''] 
		
	# preprocess the annotation color links
	if( !is.null(cnames <- names(annotation_colors)) ){
		m <- str_match(cnames, "^@([^{]+)\\{([^}]+)\\}")		
		apply(m, 1L, function(x){
			# skip unmatched names
			if( isNA(x[1]) ) return()
			
			acol <- annotation_colors[[x[1]]]
			# rename both annotation and annotation_colors if necessary
			if( x[2] != x[3] ){
				annotation[[x[3]]] <<- annotation[[x[2]]]				
				annotation[[x[2]]] <<- NULL
				if( !isNA(acol) )
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
#' \code{\link[pheatmap]{pheatmap}} from the package \code{pheatmap}, and provides
#' several enhancements such as:
#' \itemize{
#' \item it replicates the arguments of the base \code{\link{heatmap}} function 
#' \item annotation of both columns and rows, with simplified interface
#' \item easy specification of clustering methods and colors
#' \item return clustering data and grid grob object
#' }
#' 
#' Please read the associated vignette for more information and sample code.
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
#' \item one of RColorBrewer's palette name, or one of 'RdYlBu2', 'rainbow', 
#' 'heat', 'topo', 'terrain', 'cm'.
#' }
#' When the coluor palette is specified with a single value, and is negative or 
#' preceded a minus ('-'), the reversed palette is used.
#' The number of breaks can also be specified after a colon (':'). For example, 
#' the default colour palette is specified as '-RdYlBu2:100'.
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
#' from the computed dendrogram should have their root faded out.
#' This can be used to better highlight the different clusters.  
#' 
#' \item a single double that specifies how much space is used by the computed 
#' dendrogram. That is that this value is used in place of \code{treeheight}.
#' }
#' 
#' @param Colv clustering specification(s) for the columns. It accepts the same 
#' values as argument \code{Rowv} (modulo the expected length for vector specifications),  
#' and allow specifying the distance/clustering/ordering/display parameters to 
#' be used for the \emph{columns only}.
#' \code{Colv} may also be set to \code{"Rowv"}, in which case the dendrogram 
#' or ordering specifications applied to the rows are also applied to the 
#' columns. Note that this is allowed only for square input matrices, 
#' and that the row ordering is in this case by default reversed 
#' (\code{revC=TRUE}) to obtain the diagonal in the standard way 
#' (from top-left to bottom-right).
#' See argument \code{Rowv} for other possible values.
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
#' \item \code{"correlation"} and all the distances supported by \code{\link{dist}} 
#' (e.g. \code{"euclidean"}).
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
#' @param labRow labels for the rows.
#' @param labCol labels for the columns. See description for argument \code{labRow} 
#' for a list of the possible values.
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
#' @param \dots graphical parameters for the text used in plot. Parameters passed to 
#' \code{\link{grid.text}}, see \code{\link{gpar}}. 
#' 
#' @author 
#' Original version of \code{\link[pheatmap]{pheatmap}}: Raivo Kolde <rkolde@@gmail.com>
#' 
#' Enhancement into \code{aheatmap}: Renaud Gaujoux
#' 
#' @examples
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
#' # do not cluster the rows 
#' aheatmap(x, Rowv = NA)
#' # no heatmap legend
#' aheatmap(x, legend = FALSE)
#' # cell and font size 
#' aheatmap(x, cellwidth = 10, cellheight = 5)
#' aheatmap(x, cellwidth = 15, cellheight = 12, fontsize = 8, filename = "aheatmap.pdf")
#'
#' # Generate column annotations
#' annotation = data.frame(Var1 = factor(1:10 \%\% 2 == 0, labels = c("Class1", "Class2")), Var2 = 1:10)
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
#' @export
aheatmap = function(x
, color = '-RdYlBu2:100'
, breaks = NA, border_color=NA, cellwidth = NA, cellheight = NA, scale = "none"
, Rowv=TRUE, Colv=TRUE, revC = identical(Colv, "Rowv") || isNA(Rowv) || (is.integer(Rowv) && length(Rowv) > 1)
, distfun = "euclidean", hclustfun = "complete", reorderfun = function(d,w) reorder(d,w)
, treeheight = 50
, legend = TRUE, annCol = NA, annRow = NA, annColors = NA, annLegend = TRUE
, labRow = NULL, labCol = NULL
, subsetRow = NULL, subsetCol = NULL
, fontsize=10, cexRow = min(0.2 + 1/log10(nr), 1.2), cexCol = min(0.2 + 1/log10(nc), 1.2)
, filename = NA, width = NA, height = NA
, main = NULL, sub = NULL, info = NULL
, verbose=getOption('verbose'), ...){

	# set verbosity level
	ol <- lverbose(verbose)
	on.exit( lverbose(ol) )
	
	# convert ExpressionSet into 
	if( is(x, 'ExpressionSet') ){
		library(Biobase)
		x <- exprs(x)
	}

	# rename to old parameter name
	mat <- x
	
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
		rownames(mat)  <- generate_dimnames(labRow, nrow(mat))	
	}
	# label columns numerically if no colnames
	if( is.null(labCol) && is.null(colnames(mat)) )
		labCol <- 1L
	if( !is.null(labCol) ){
		if( verbose ) message("Process labCol")
		colnames(mat) <- generate_dimnames(labCol, ncol(mat))
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
	tree_row <- if( !isNA(Rowv) ){
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
	tree_col <- if( !isNA(Colv) ){
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
	if( !isNA(tree_row) ){
		
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
		mat <- mat[if( is.null(subInd) ) res$rowInd else subInd, , drop=FALSE]
	}
	
	if( !isNA(tree_col) ){		
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
		mat <- mat[, if( is.null(subInd) ) res$colInd else subInd, drop=FALSE]		
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
	# load named palette if necessary
	color <- ccRamp(color)
	
	# generate breaks if necessary
	if( isNA(breaks) || isNumber(breaks) ){
		if( verbose ) message("Generate breaks")
		# if a single number: center the breaks on this value
		cbreaks <- if( isNumber(breaks) ) breaks else NA
		breaks = generate_breaks(as.vector(mat), length(color), center=cbreaks)
	}
	if (legend) {
		if( verbose ) message("Generate data legend breaks")
		legend = grid.pretty(range(as.vector(breaks)))
	}
	else {
		legend = NA
	}
	mat = scale_colours(mat, col = color, breaks = breaks)
	
	annotation_legend <- annLegend
	annotation_colors <- annColors
	
	# render annotation tracks for both rows and columns
	annTracks <- renderAnnotations(atrack(annCol, order=res$colInd, .DATA=amargin(x,2L))
								, atrack(annRow, order=res$rowInd, .DATA=amargin(x,1L))
								, annotation_colors = annotation_colors, verbose=verbose)
	
	# retrieve dimension for computing cexRow and cexCol (evaluated from the arguments)
	nr <- nrow(mat); nc <- ncol(mat)
	# Draw heatmap	
	res$vp <- heatmap_motor(mat, border_color = border_color, cellwidth = cellwidth, cellheight = cellheight
	, treeheight_col = treeheight_col, treeheight_row = treeheight_row, tree_col = tree_col, tree_row = tree_row
	, filename = filename, width = width, height = height, breaks = breaks, color = color, legend = legend
	, annTracks = annTracks, annotation_legend = annotation_legend
	, fontsize = fontsize, fontsize_row = cexRow * fontsize, fontsize_col = cexCol * fontsize
	, main = main, sub = sub, info = info
	, verbose = verbose
	, ...)
	
	# return info about the plot
	invisible(res)
}

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
	library(gridBase)
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
