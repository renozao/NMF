#' @include rmatrix.R
#' @include nmf-package.R
NULL

library(pkgmaker)

#' Utility Function in the NMF Package
#' 
#' @name utils-NMF
#' @rdname utils
NULL

#' Internal verbosity option
#' @param val logical that sets the verbosity level.
#' @return the old verbose level   
#' @keywords internal
lverbose <- local({
			.val <- 0
			function(val){
				if( missing(val) ) return(.val)
				oval <- .val
				.val <<- val
				invisible(oval)
			}
		})
vmessage <- function(..., appendLF=TRUE) if( lverbose() ) cat(..., if(appendLF) "\n", sep='')


nmf_stop <- function(name, ...){
	stop("NMF::", name , ' - ', ..., call.=FALSE)
}
nmf_warning <- function(name, ...){
	warning("NMF::", name , ' - ', ..., call.=FALSE)
}

# or-NULL operator
'%||%' <- function(x, y) if( !is.null(x) ) x else y 

# cat object or class for nice cat/message
quick_str <- function(x) if( is.atomic(x) ) x else class(x)[1]

# remove all attributes from an object
rmAttributes <- function(x){
	attributes(x) <- NULL
	x
}

#' \code{str_args} formats the arguments of a function using \code{\link{args}}, 
#' but returns the output as a string.
#' 
#' @param x a function
#' @param exdent indentation for extra lines if the output takes more than one line.  
#' 
#' @export
#' @rdname utils
#' 
#' @examples 
#' 
#' args(library)
#' str_args(library)
#' 
str_args <- function(x, exdent=10L){
	s <- capture.output(print(args(x)))
	paste(str_trim(s[-length(s)]), collapse=str_c('\n', paste(rep(' ', exdent), collapse='')))
}

#' Simple Progress Bar
#' 
#' Creates a simple progress bar with title.
#' This function is identical to \code{utils::txtProgressBar} but allow adding 
#' a title to the progress bar, and can be shared by multiple processes, 
#' e.g., in multicore or multi-hosts computations.
#' 
#' @inheritParams utils::txtProgressBar
#' @param shared specification of a shared directory to use when the progress
#' bar is to be used by multiple processes.
#' 
#' @author R Core Team
#' @keywords internal
txtProgressBar <- function (min = 0, max = 1, initial = 0, char = "=", width = NA, 
    	title= if( style == 3 ) ' ', label, style = 1, file = ""
		, shared = NULL)
{
    if (!identical(file, "") && !(inherits(file, "connection") && 
        		isOpen(file))) 
        stop("'file' must be \"\" or an open connection object")
    if (!style %in% 1L:3L) 
        style <- 1
    .val <- initial
    .killed <- FALSE
    .nb <- 0L
    .pc <- -1L
    nw <- nchar(char, "w")
    if (is.na(width)) {
        width <- getOption("width")
        if (style == 3L) 
            width <- width - 10L
        width <- trunc(width/nw)
    }
    if (max <= min) 
        stop("must have max > min")
	
	# setup shared directory
	.shared <- NULL
	if( isTRUE(shared) ) shared <- tempdir()
	if( is.character(shared) ){
		.shared <- tempfile('progressbar_', tmpdir=shared[1L])
		dir.create(.shared)
	}
	#
	
	getval <- function(value){
		if( value >= max || value <= min 
			|| is.null(.shared) || !file.exists(.shared) ){
			value
		}else{
			cat('', file=file.path(.shared, paste('_', value, sep='')))
			length(list.files(.shared))
		}
	}
	
    up1 <- function(value) {
		if (!is.finite(value) || value < min || value > max) 
            return()
		
		# get actual value 
		value <- getval(value)
		
        .val <<- value
        nb <- round(width * (value - min)/(max - min))
        if (.nb < nb) {
            cat(paste(rep.int(char, nb - .nb), collapse = ""), 
                	file = file)
            flush.console()
        }
        else if (.nb > nb) {
            cat("\r", title, paste(rep.int(" ", .nb * nw), collapse = ""), 
                	"\r", title, paste(rep.int(char, nb), collapse = ""), 
                	sep = "", file = file)
            flush.console()
        }
        .nb <<- nb
    }
    up2 <- function(value) {
        if (!is.finite(value) || value < min || value > max) 
            return()
		
		# get actual value 
		value <- getval(value)
		
        .val <<- value
        nb <- round(width * (value - min)/(max - min))
        if (.nb <= nb) {
            cat("\r", title, paste(rep.int(char, nb), collapse = ""), 
                	sep = "", file = file)
            flush.console()
        }
        else {
            cat("\r", title, paste(rep.int(" ", .nb * nw), collapse = ""), 
                	"\r", paste(rep.int(char, nb), collapse = ""), 
                	sep = "", file = file)
            flush.console()
        }
        .nb <<- nb
    }
    up3 <- function(value) {
        if (!is.finite(value) || value < min || value > max) 
            return()
		
		# get actual value 
		value <- getval(value)
		
        .val <<- value
        nb <- round(width * (value - min)/(max - min))
        pc <- round(100 * (value - min)/(max - min))
        if (nb == .nb && pc == .pc) 
            return()
        cat(paste(c("\r",title," |", rep.int(" ", nw * width + 6)), collapse = ""), 
            	file = file)
        cat(paste(c("\r",title," |", rep.int(char, nb), rep.int(" ", 
            							nw * (width - nb)), sprintf("| %3d%%", pc)), collapse = ""), 
            	file = file)
        flush.console()
        .nb <<- nb
        .pc <<- pc
    }
    getVal <- function() .val
    kill <- function(cleanup=TRUE) if (!.killed) {
        	cat("\n", file = file)
        	flush.console()
			.killed <<- TRUE
			
			# do some cleanup
			if( cleanup ){
				# delete shared directory
				if( !is.null(.shared) && file.exists(.shared) ) 
					unlink(.shared, recursive=TRUE)
				#
			}
			invisible(TRUE)
    	}
    up <- switch(style, up1, up2, up3)
    up(initial)
    structure(list(getVal = getVal, up = up, kill = kill), class = "txtProgressBar")
}

###% apply a function to each entry in a matrix
matapply <- function(x, FUN, ...){
	res <- sapply(x, FUN, ...)
	matrix(res, nrow(x))
}

###% try to convert a character string into a numeric
toNumeric <- function(x){
	suppressWarnings( as.numeric(x) )
}

###% Tells one is running in Sweave
isSweave <- function() !is.null(sweaveLabel()) 
	
sweaveLabel <- function(){
	if ((n.parents <- length(sys.parents())) >= 3) {
		for (i in seq_len(n.parents) - 1) {
			if ("chunkopts" %in% ls(envir = sys.frame(i))) {
				chunkopts = get("chunkopts", envir = sys.frame(i))
				if (all(c("prefix.string", "label") %in% names(chunkopts))) {
					img.name = paste(chunkopts$prefix.string, chunkopts$label, 
							sep = "-")
					return(img.name)
					break
				}
			}
		}
	}
}

sweaveFile <- function(){
	label <- sweaveLabel()
	if( !is.null(label) )
		paste(label, '.pdf', sep='')
}

fixSweaveFigure <- function(filename){
	if( missing(filename) ){
		filename <- sweaveLabel()
		if( is.null(filename) ) return()
		filename <- paste(filename, '.pdf', sep='')
	}
	filepath <- normalizePath(filename)
	tf <- tempfile()
	system(paste("pdftk", filepath, "cat 2-end output", tf, "; mv -f", tf, filepath))
}

###% 'more' functionality to read data progressively
more <- function(x, step.size=10, width=20, header=FALSE, pattern=NULL){
	
	if( !(is.matrix(x) || is.data.frame(x) || is.vector(x) || is.list(x)) )
		stop("NMF::more - invalid argument 'x': only 'matrix', 'data.frame', 'vector' and 'list' objects are handled.")
	
	one.dim <- is.null(dim(x))
	single.char <- FALSE
	n <-
		if( is.character(x) && length(x) == 1 ){			
			cat("<character string:", nchar(x), ">\n")
			single.char <- TRUE
			nchar(x)
		}
		else if( one.dim ){
			cat("<", class(x),":", length(x), ">\n")
			
			# limit to matching terms if necessary
			if( !is.null(pattern) )
				x[grep(pattern, x)]
			
			length(x)
		}else{
			cat("<", class(x),":", nrow(x), "x", ncol(x), ">\n")
			head.init <- colnames(x)
			head.on <- TRUE
			
			# limit to matching terms if necessary
			if( !is.null(pattern) ){
				idx <- apply(x, 2, grep, pattern=pattern)
				print(idx)
				idx <- unique(if( is.list(idx) ) unlist(idx) else as.vector(idx))
				x <- x[idx,, drop=FALSE]
			}
			
			nrow(x)
		}	
		
	i <- 0
	while( i < n ){
		# reduce 'step.size' if necessary
		step.size <- min(step.size, n-i)
		
		what2show <- if( single.char )
			substr(x, i+1, i+step.size)
		else if( one.dim )			
			if( !is.na(width) ) sapply(x[seq(i+1, i+step.size)], function(s) substr(s, 1, width) ) else x[seq(i+1, i+step.size)]
		else{
			w <- x[seq(i+1, i+step.size), , drop=FALSE]
			if( !is.na(width) ){ 
				w <- apply(w, 2, 
					function(s){
						ns <- toNumeric(s)
						if( !is.na(ns[1]) ) # keep numerical value as is
							ns
						else # limit output if required
							substr(s, 1, width)
						
					}) 
				rownames(w) <- rownames(x)[seq(i+1, i+step.size)]
			} 
				
			
			# remove header if not required
			if( !header && head.on ){
				colnames(x) <- sapply(colnames(x), function(c) paste(rep(' ', nchar(c)), collapse=''))
				head.on <- FALSE
			}
			
			# return the content
			w
		}
		
		cat( show(what2show) )
		i <- i + step.size
		
		# early break if necessary
		if( i >= n )
			break
		# ask user what to to next
		ans <- scan(what='character', quiet=TRUE, n=1, multi.line=FALSE)
		
		# process user command if any (otherwise carry on)
		if( length(ans) > 0 ){		
			if( !is.na(s <- toNumeric(ans)) ) # change step size
				step.size <- s
			else if( !header && ans %in% c('h', 'head') ){
				colnames(x) <- head.init
				head.on <- TRUE
			}
			else if( ans %in% c('q', 'quit') ) # quit
				break
		}
	}
	invisible()
} 

#' Randomizing Data
#' 
#' \code{randomize} permutates independently the entries in each column 
#' of a matrix-like object, to produce random data that can be used 
#' in permutation tests or bootstrap analysis.
#' 
#' In the context of NMF, it may be used to generate random data, whose 
#' factorization serves as a reference for selecting a factorization rank,
#' that does not overfit the data.
#' 
#' @param x data to be permutated. It must be an object suitable to be 
#' passed to the function \code{\link{apply}}.
#' @param ... extra arguments passed to the function \code{\link{sample}}.
#' 
#' @return a matrix
#' 
#' @export
#' @examples
#' x <- matrix(1:32, 4, 8)
#' randomize(x)
#' randomize(x)
#' 
randomize <- function(x, ...){
	
	if( is(x, 'ExpressionSet') ) x <- Biobase::exprs(x)
		
	# resample the columns
	apply(x, 2, function(c, ...) sample(c, size=length(c), ...), ...)
	
}

###% Returns the rank-k truncated SVD approximation of x
tsvd <- function(x, r, ...){
	stopifnot( r > 0 && r <= min(dim(x)))
	s <- svd(x, nu=r, nv=r, ...)
	s$d <- s$d[1:r]
	
	# return results
	s
}

###% Subset a list leaving only the arguments from a given function 
.extract.args <- function(x, fun, ...){
	
	fdef <- if( is.character(fun) )	getFunction(fun, ...)
			else if( is.function(fun) ) fun
			else stop("invalid argument 'fun': expected function name or definition")
	
	if( length(x) == 0 ) return(x)	
	x.ind <- charmatch(if( is.list(x) ) names(x) else x, args <- formalArgs(fdef))
	x[!is.na(x.ind)]
}

###% Returns the version of the package
nmfInfo <- function(command){	
	pkg <- 'NMF'
	curWarn <- getOption("warn")
	on.exit(options(warn = curWarn), add = TRUE)
	options(warn = -1)
	desc <- packageDescription(pkg, fields="Version")
	if (is.na(desc)) 
		stop(paste("Package", pkg, "not found"))
	desc
}

###% Returns TRUE if running under Mac OS X + GUI
is.Mac <- function(check.gui=FALSE){
	is.mac <- (length(grep("darwin", R.version$platform)) > 0)
	# return TRUE is running on Mac (adn optionally through GUI)
	is.mac && (!check.gui || .Platform$GUI == 'AQUA')
}

###% Hash a function body (using digest)
#' @import digest
hash_function <- function(f){
	b <- body(f)
	attributes(b) <- NULL
	fdef <- paste(c(capture.output(args(f))[1], capture.output(print(b))), collapse="\n")
	# print(fdef)
	digest(b)
}


###% compare function with copy and with no copy
cmp.cp <- function(...){
	res <- nmf(..., copy=F)
	resc <- nmf(..., copy=T)
	cat("identical: ", identical(fit(res), fit(resc))
			, " - all.equal: ", all.equal(fit(res), fit(resc))
			, " - diff: ", all.equal(fit(res), fit(resc), tol=0)
			, "\n"
	)
	invisible(res)
} 

# return the internal pointer address 
C.ptr <- function(x, rec=FALSE)
{	
	attribs <- attributes(x)
	if( !rec || is.null(attribs) )
		.Call("ptr_address", x, PACKAGE='NMF')
	else
		c( C.ptr(x), sapply(attribs, C.ptr, rec=TRUE))
	
}

is.same <- function(x, y){
	C.ptr(x) == C.ptr(y)
}

is.eset <- function(x) is(x, 'ExpressionSet')

# clone an object
clone <- function(x){
	.Call('clone_object', x, PACKAGE='NMF')
}

# deep-clone an object
clone2 <- function(x){
	if( is.environment(x) ){
		y <- Biobase::copyEnv(x)
		eapply(ls(x, all.names=TRUE), 
			function(n){
				if( is.environment(x[[n]]) ){
					y[[n]] <<- clone(x[[n]])
					if( identical(parent.env(x[[n]]), x) )
						parent.env(y[[n]]) <<- y
				}
		})
	}else{
		y <- .Call('clone_object', x, PACKAGE='NMF')		
		if( isS4(x) ){ ## deep copy R object
			lapply(slotNames(class(y)), 
				function(n){					
					slot(y, n) <<- clone(slot(x, n)) 
			})
		}else if( is.list(x) ){ ## copy list or vector
			sapply(seq_along(x), 
				function(i){					
					y[[i]] <<- clone(x[[i]])					
			})
		}
	}
	
	y
}

#compute RSS with C function
.rss <- function(x, y)
{	
	.Call("Euclidean_rss", x, y, PACKAGE='NMF')
}

#compute KL divergence with C function
.KL <- function(x, y)
{	
	.Call("KL_divergence", x, y, PACKAGE='NMF')
}

#' Updating Objects In Place
#' 
#' These functions modify objects (mainly matrix objects) in place, i.e. they 
#' act directly on the C pointer.
#' Due to their side-effect, they are not meant to be called by the end-user.
#' 
#' \code{pmax.inplace} is a version of \code{\link{pmax}} that updates its first
#' argument.    
#' 
#' @param x an object to update in place.
#' @param lim lower threshold value
#' @param skip indexes to skip
#' 
#' @export
#' @rdname inplace
#' @keywords internal
pmax.inplace <- function(x, lim, skip=NULL){
	
	.Call('ptr_pmax', x, lim, as.integer(skip), PACKAGE='NMF')
	
}

# colMin
colMin <- function(x){
	.Call('colMin', x, PACKAGE='NMF')
}

# colMax
colMax <- function(x){
	.Call('colMax', x, PACKAGE='NMF')
}

#' \code{neq.constraints.inplace} apply unequality constraints in place.
#' 
#' @param constraints constraint specification.
#' @param ratio fixed ratio on which the constraint applies.
#' @param value fixed value to enforce.
#' @param copy a logical that indicates if \code{x} should be updated in place 
#' or not.
#' 
#' @export
#' @rdname inplace
neq.constraints.inplace <- function(x, constraints, ratio=NULL, value=NULL, copy=FALSE){
	
	# if requested: clone data as neq.constrains.inplace modify the input data in place
	if( copy )
		x <- clone(x)
	
	.Call('ptr_neq_constraints', x, constraints, ratio, value, PACKAGE='NMF')	
}

# Test if an external pointer is nil
# Taken from package bigmemory 
ptr_isnil <- function (address) 
{
	if (class(address) != "externalptr") 
		stop("address is not an externalptr.")
	.Call("ptr_isnil", address, PACKAGE='NMF')	
}


###% Draw the palette of colors
###% 
###% Taken from the examples of colorspace::rainbow_hcl
###% 
pal <- function(col, h=1, border = "light gray")
{
	n <- length(col)	
	plot(0, 0, type="n", xlim = c(0, 1), ylim = c(0, h), axes = FALSE, xlab = "", ylab = "")
	rect(0:(n-1)/n, 0, 1:n/n, h, col = col, border = border)
}
###% Draw the Palette of Colors as a Wheel
###% 
###% Taken from the examples of colorspace::rainbow_hcl
###% 
wheel <- function(col, radius = 1, ...)
	pie(rep(1, length(col)), col = col, radius = radius, ...)

# Define a S4 class to handle function slots given as either a function definition 
# or a character string that gives the function's name. 
setClassUnion('.functionSlot', c('character', 'function'))

# Define a S4 class to handle function slots given as either a function definition 
# or a character string that gives the function's name or NULL.
setClassUnion('.functionSlotNULL', c('character', 'function', 'NULL'))
.validFunctionSlot <- function(slot, allow.empty=FALSE, allow.null=TRUE){
	if( is.null(slot) ){
		if( !allow.null ) return('NULL value is not allowed')
		return(TRUE)
	}
	if( is.character(slot) ){
		if( !allow.empty && slot == '' ) return('character string cannot be empty')
		if( length(slot) != 1 ) return(paste('character string must be a single value [length =', length(slot), ']', sep=''))
	}			
	
	return(TRUE)
}

# Returns the number of cores to use
getNCores <- function(){
	rv <- paste(R.version$major, R.version$minor, sep='.') 
	if( utils::compareVersion(rv, "2.14.0") >=0 && require(parallel) )
		parallel:::detectCores()
	else{
		# extracted from parallel::detectCores() in R-2.14.0
		all.tests <- FALSE; logical <- FALSE;
		systems <- list(darwin = "/usr/sbin/sysctl -n hw.ncpu 2>/dev/null", 
				freebsd = "/sbin/sysctl -n hw.ncpu 2>/dev/null", linux = "grep processor /proc/cpuinfo 2>/dev/null | wc -l", 
				irix = c("hinv | grep Processors | sed 's: .*::'", "hinv | grep '^Processor '| wc -l"), 
				solaris = if (logical) "/usr/sbin/psrinfo -v | grep 'Status of.*processor' | wc -l" else "/usr/sbin/psrinfo -p")
		for (i in seq(systems)) if (all.tests || length(grep(paste("^", 
									names(systems)[i], sep = ""), R.version$os))) 
				for (cmd in systems[i]) {
					a <- gsub("^ +", "", system(cmd, TRUE)[1])
					if (length(grep("^[1-9]", a))) 
						return(as.integer(a))
				}
		NA_integer_		
	}
}

####% Utility function needed in heatmap.plus.2
#invalid <- function (x) 
#{
#	if (missing(x) || is.null(x) || length(x) == 0) 
#		return(TRUE)
#	if (is.list(x)) 
#		return(all(sapply(x, invalid)))
#	else if (is.vector(x)) 
#		return(all(is.na(x)))
#	else return(FALSE)
#}
#

####% Modification of the function heatmap.2 including a small part of function 
####% heatmap.plus to allow extra annotation rows
#heatmap.plus.2 <- 
#		function (x, Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE, 
#				distfun = dist, hclustfun = hclust, dendrogram = c("both", 
#						"row", "column", "none"), symm = FALSE, scale = c("none", 
#						"row", "column"), na.rm = TRUE, revC = identical(Colv, 
#						"Rowv"), add.expr, breaks, symbreaks = min(x < 0, na.rm = TRUE) || 
#						scale != "none", col = "heat.colors", colsep, rowsep, 
#				sepcolor = "white", sepwidth = c(0.05, 0.05), cellnote, notecex = 1, 
#				notecol = "cyan", na.color = par("bg"), trace = c("column", 
#						"row", "both", "none"), tracecol = "cyan", hline = median(breaks), 
#				vline = median(breaks), linecol = tracecol, margins = c(5, 
#						5), ColSideColors, RowSideColors, cexRow = 0.2 + 1/log10(nr), 
#				cexCol = 0.2 + 1/log10(nc), labRow = NULL, labCol = NULL, 
#				key = TRUE, keysize = 1.5, density.info = c("histogram", 
#						"density", "none"), denscol = tracecol, symkey = min(x < 
#								0, na.rm = TRUE) || symbreaks, densadj = 0.25, main = NULL, 
#				xlab = NULL, ylab = NULL, lmat = NULL, lhei = NULL, lwid = NULL, 
#				...) 
#{
#	scale01 <- function(x, low = min(x), high = max(x)) {
#		x <- (x - low)/(high - low)
#		x
#	}
#	retval <- list()
#	scale <- if (symm && missing(scale)) 
#				"none"
#			else match.arg(scale)
#	dendrogram <- match.arg(dendrogram)
#	trace <- match.arg(trace)
#	density.info <- match.arg(density.info)
#	if (length(col) == 1 && is.character(col)) 
#		col <- get(col, mode = "function")
#	if (!missing(breaks) && (scale != "none")) 
#		warning("Using scale=\"row\" or scale=\"column\" when breaks are", 
#				"specified can produce unpredictable results.", "Please consider using only one or the other.")
#	if (is.null(Rowv) || is.na(Rowv)) 
#		Rowv <- FALSE
#	if (is.null(Colv) || is.na(Colv)) 
#		Colv <- FALSE
#	else if (Colv == "Rowv" && !isTRUE(Rowv)) 
#		Colv <- FALSE
#	if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
#		stop("`x' must be a numeric matrix")
#	nr <- di[1]
#	nc <- di[2]
#	if (nr <= 1 || nc <= 1) 
#		stop("`x' must have at least 2 rows and 2 columns")
#	if (!is.numeric(margins) || length(margins) != 2) 
#		stop("`margins' must be a numeric vector of length 2")
#	if (missing(cellnote)) 
#		cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
#	if (!inherits(Rowv, "dendrogram")) {
#		if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in% 
#					c("both", "row"))) {
#			if (is.logical(Colv) && (Colv)) 
#				dendrogram <- "column"
#			else dedrogram <- "none"
#			warning("Discrepancy: Rowv is FALSE, while dendrogram is `", 
#					dendrogram, "'. Omitting row dendogram.")
#		}
#	}
#	if (!inherits(Colv, "dendrogram")) {
#		if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in% 
#					c("both", "column"))) {
#			if (is.logical(Rowv) && (Rowv)) 
#				dendrogram <- "row"
#			else dendrogram <- "none"
#			warning("Discrepancy: Colv is FALSE, while dendrogram is `", 
#					dendrogram, "'. Omitting column dendogram.")
#		}
#	}
#	if (inherits(Rowv, "dendrogram")) {
#		ddr <- Rowv
#		rowInd <- order.dendrogram(ddr)
#	}
#	else if (is.integer(Rowv)) {
#		hcr <- hclustfun(distfun(x))
#		ddr <- as.dendrogram(hcr)
#		ddr <- reorder(ddr, Rowv)
#		rowInd <- order.dendrogram(ddr)
#		if (nr != length(rowInd)) 
#			stop("row dendrogram ordering gave index of wrong length")
#	}
#	else if (isTRUE(Rowv)) {
#		Rowv <- rowMeans(x, na.rm = na.rm)
#		hcr <- hclustfun(distfun(x))
#		ddr <- as.dendrogram(hcr)
#		ddr <- reorder(ddr, Rowv)
#		rowInd <- order.dendrogram(ddr)
#		if (nr != length(rowInd)) 
#			stop("row dendrogram ordering gave index of wrong length")
#	}
#	else {
#		rowInd <- nr:1
#	}
#	if (inherits(Colv, "dendrogram")) {
#		ddc <- Colv
#		colInd <- order.dendrogram(ddc)
#	}
#	else if (identical(Colv, "Rowv")) {
#		if (nr != nc) 
#			stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
#		if (exists("ddr")) {
#			ddc <- ddr
#			colInd <- order.dendrogram(ddc)
#		}
#		else colInd <- rowInd
#	}
#	else if (is.integer(Colv)) {
#		hcc <- hclustfun(distfun(if (symm) 
#									x
#								else t(x)))
#		ddc <- as.dendrogram(hcc)
#		ddc <- reorder(ddc, Colv)
#		colInd <- order.dendrogram(ddc)
#		if (nc != length(colInd)) 
#			stop("column dendrogram ordering gave index of wrong length")
#	}
#	else if (isTRUE(Colv)) {
#		Colv <- colMeans(x, na.rm = na.rm)
#		hcc <- hclustfun(distfun(if (symm) 
#									x
#								else t(x)))
#		ddc <- as.dendrogram(hcc)
#		ddc <- reorder(ddc, Colv)
#		colInd <- order.dendrogram(ddc)
#		if (nc != length(colInd)) 
#			stop("column dendrogram ordering gave index of wrong length")
#	}
#	else {
#		colInd <- 1:nc
#	}
#	retval$rowInd <- rowInd
#	retval$colInd <- colInd
#	retval$call <- match.call()
#	x <- x[rowInd, colInd]
#	x.unscaled <- x
#	cellnote <- cellnote[rowInd, colInd]
#	if (is.null(labRow)) 
#		labRow <- if (is.null(rownames(x))) 
#					(1:nr)[rowInd]
#				else rownames(x)
#	else labRow <- labRow[rowInd]
#	if (is.null(labCol)) 
#		labCol <- if (is.null(colnames(x))) 
#					(1:nc)[colInd]
#				else colnames(x)
#	else labCol <- labCol[colInd]
#	if (scale == "row") {
#		retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
#		x <- sweep(x, 1, rm)
#		retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
#		x <- sweep(x, 1, sx, "/")
#	}
#	else if (scale == "column") {
#		retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
#		x <- sweep(x, 2, rm)
#		retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
#		x <- sweep(x, 2, sx, "/")
#	}
#	if (missing(breaks) || is.null(breaks) || length(breaks) < 
#			1) {
#		if (missing(col) || is.function(col)) 
#			breaks <- 16
#		else breaks <- length(col) + 1
#	}
#	if (length(breaks) == 1) {
#		if (!symbreaks) 
#			breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm), 
#					length = breaks)
#		else {
#			extreme <- max(abs(x), na.rm = TRUE)
#			breaks <- seq(-extreme, extreme, length = breaks)
#		}
#	}
#	nbr <- length(breaks)
#	ncol <- length(breaks) - 1
#	if (class(col) == "function") 
#		col <- col(ncol)
#	min.breaks <- min(breaks)
#	max.breaks <- max(breaks)
#	x[x < min.breaks] <- min.breaks
#	x[x > max.breaks] <- max.breaks
#	if (missing(lhei) || is.null(lhei)) 
#		lhei <- c(keysize, 4)
#	if (missing(lwid) || is.null(lwid)) 
#		lwid <- c(keysize, 4)
#	if (missing(lmat) || is.null(lmat)) {
#		lmat <- rbind(4:3, 2:1)
#		
#		# hack for adding extra annotations
#		if (!missing(ColSideColors)) {
#			if (!is.matrix(ColSideColors)) 
#				stop("'ColSideColors' must be a matrix")
#			
#			if (!is.character(ColSideColors) || dim(ColSideColors)[1] != 
#					nc) 
#				stop("'ColSideColors' must be a character vector/matrix with length/ncol = ncol(x)")
#			lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
#			lhei <- c(lhei[1], 0.2, lhei[2])
#		}
#		if (!missing(RowSideColors)) {
#			if (!is.character(RowSideColors) || length(RowSideColors) != 
#					nr) 
#				stop("'RowSideColors' must be a character vector of length nrow(x)")
#			lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 
#											1), 1), lmat[, 2] + 1)
#			lwid <- c(lwid[1], 0.2, lwid[2])
#		}
#		lmat[is.na(lmat)] <- 0
#	}
#	if (length(lhei) != nrow(lmat)) 
#		stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
#	if (length(lwid) != ncol(lmat)) 
#		stop("lwid must have length = ncol(lmat) =", ncol(lmat))
#	op <- par(no.readonly = TRUE)
#	on.exit(par(op))
#	layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
#	if (!missing(RowSideColors)) {
#		par(mar = c(margins[1], 0, 0, 0.5))
#		image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
#	}
#	if (!missing(ColSideColors)) {
#		par(mar = c(0.5, 0, 0, margins[2]))
#		csc = ColSideColors[colInd, ]
#		csc.colors = matrix()
#		csc.names = names(table(csc))
#		csc.i = 1
#		for (csc.name in csc.names) {
#			csc.colors[csc.i] = csc.name
#			csc[csc == csc.name] = csc.i
#			csc.i = csc.i + 1
#		}
#		csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
#		image(csc, col = as.vector(csc.colors), axes = FALSE)
#		if (length(colnames(ColSideColors)) > 0) {
#			axis(2, 0:(dim(csc)[2] - 1)/(dim(csc)[2] - 1), colnames(ColSideColors), 
#					las = 2, tick = FALSE, cex.axis= cexRow)
#		}
#	}
#	par(mar = c(margins[1], 0, 0, margins[2]))
#	if (!symm || scale != "none") {
#		x <- t(x)
#		cellnote <- t(cellnote)
#	}
#	if (revC) {
#		iy <- nr:1
#		if (exists("ddr")) 
#			ddr <- rev(ddr)
#		x <- x[, iy]
#		cellnote <- cellnote[, iy]
#	}
#	else iy <- 1:nr
#	image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
#					c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, 
#			breaks = breaks, ...)
#	retval$carpet <- x
#	if (exists("ddr")) 
#		retval$rowDendrogram <- ddr
#	if (exists("ddc")) 
#		retval$colDendrogram <- ddc
#	retval$breaks <- breaks
#	retval$col <- col
#	if (!invalid(na.color) & any(is.na(x))) {
#		mmat <- ifelse(is.na(x), 1, NA)
#		image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "", 
#				col = na.color, add = TRUE)
#	}
#	axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0, 
#			cex.axis = cexCol)
#	if (!is.null(xlab)) 
#		mtext(xlab, side = 1, line = margins[1] - 1.25)
#	axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0, 
#			cex.axis = cexRow)
#	if (!is.null(ylab)) 
#		mtext(ylab, side = 4, line = margins[2] - 1.25)
#	if (!missing(add.expr)) 
#		eval(substitute(add.expr))
#	if (!missing(colsep)) 
#		for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, 
#							length(csep)), xright = csep + 0.5 + sepwidth[1], 
#					ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, 
#					col = sepcolor, border = sepcolor)
#	if (!missing(rowsep)) 
#		for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 
#								1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 
#								1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, 
#					col = sepcolor, border = sepcolor)
#	min.scale <- min(breaks)
#	max.scale <- max(breaks)
#	x.scaled <- scale01(t(x), min.scale, max.scale)
#	if (trace %in% c("both", "column")) {
#		retval$vline <- vline
#		vline.vals <- scale01(vline, min.scale, max.scale)
#		for (i in colInd) {
#			if (!is.null(vline)) {
#				abline(v = i - 0.5 + vline.vals, col = linecol, 
#						lty = 2)
#			}
#			xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
#			xv <- c(xv[1], xv)
#			yv <- 1:length(xv) - 0.5
#			lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
#		}
#	}
#	if (trace %in% c("both", "row")) {
#		retval$hline <- hline
#		hline.vals <- scale01(hline, min.scale, max.scale)
#		for (i in rowInd) {
#			if (!is.null(hline)) {
#				abline(h = i + hline, col = linecol, lty = 2)
#			}
#			yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
#			yv <- rev(c(yv[1], yv))
#			xv <- length(yv):1 - 0.5
#			lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
#		}
#	}
#	if (!missing(cellnote)) 
#		text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote), 
#				col = notecol, cex = notecex)
#	par(mar = c(margins[1], 0, 0, 0))
#	if (dendrogram %in% c("both", "row")) {
#		plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
#	}
#	else plot.new()
#	par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
#	if (dendrogram %in% c("both", "column")) {
#		plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
#	}
#	else plot.new()
#	if (!is.null(main)) 
#		title(main, cex.main = 1.5 * op[["cex.main"]])
#	if (key) {
#		par(mar = c(5, 4, 2, 1), cex = 0.75)
#		tmpbreaks <- breaks
#		if (symkey) {
#			max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
#			min.raw <- -max.raw
#			tmpbreaks[1] <- -max(abs(x))
#			tmpbreaks[length(tmpbreaks)] <- max(abs(x))
#		}
#		else {
#			min.raw <- min(x, na.rm = TRUE)
#			max.raw <- max(x, na.rm = TRUE)
#		}
#		z <- seq(min.raw, max.raw, length = length(col))
#		image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks, 
#				xaxt = "n", yaxt = "n")
#		par(usr = c(0, 1, 0, 1))
#		lv <- pretty(breaks)
#		xv <- scale01(as.numeric(lv), min.raw, max.raw)
#		axis(1, at = xv, labels = lv)
#		if (scale == "row") 
#			mtext(side = 1, "Row Z-Score", line = 2)
#		else if (scale == "column") 
#			mtext(side = 1, "Column Z-Score", line = 2)
#		else mtext(side = 1, "Value", line = 2)
#		if (density.info == "density") {
#			dens <- density(x, adjust = densadj, na.rm = TRUE)
#			omit <- dens$x < min(breaks) | dens$x > max(breaks)
#			dens$x <- dens$x[-omit]
#			dens$y <- dens$y[-omit]
#			dens$x <- scale01(dens$x, min.raw, max.raw)
#			lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol, 
#					lwd = 1)
#			axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
#			title("Color Key\nand Density Plot")
#			par(cex = 0.5)
#			mtext(side = 2, "Density", line = 2)
#		}
#		else if (density.info == "histogram") {
#			h <- hist(x, plot = FALSE, breaks = breaks)
#			hx <- scale01(breaks, min.raw, max.raw)
#			hy <- c(h$counts, h$counts[length(h$counts)])
#			lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s", 
#					col = denscol)
#			axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
#			title("Color Key\nand Histogram")
#			par(cex = 0.5)
#			mtext(side = 2, "Count", line = 2)
#		}
#		else title("Color Key")
#	}
#	else plot.new()
#	retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)], 
#			high = retval$breaks[-1], color = retval$col)
#	invisible(retval)
#}

# Extracted from the psychometric package (0.1.0)
# Copyright Thomas D. Fletcher
# Under Gnu GPL2
CI.Rsqlm <- function (obj, level = 0.95) 
{
	l <- level
	rsq <- summary(obj)$r.squared
	k <- summary(obj)$df[1] - 1
	n <- obj$df + k + 1
	mat <- CI.Rsq(rsq, n, k, level = l)
	return(mat)
}
# Extracted from the psychometric package (0.1.0)
# Copyright Thomas D. Fletcher
# Under Gnu GPL2
CI.Rsq <- function (rsq, n, k, level = 0.95) 
{
	noma <- 1 - level
	sersq <- sqrt((4 * rsq * (1 - rsq)^2 * (n - k - 1)^2)/((n^2 - 
							1) * (n + 3)))
	zs <- -qnorm(noma/2)
	mez <- zs * sersq
	lcl <- rsq - mez
	ucl <- rsq + mez
	mat <- data.frame(Rsq = rsq, SErsq = sersq, LCL = lcl, UCL = ucl)
	return(mat)
}
