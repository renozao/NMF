# Grid related functions
# 
# Mainly functions that duplicate grid functions but do not create a new plot
# if none is present.
#
# Author: Renaud Gaujoux
# Creation: 04 Jun 2012
###############################################################################

#' Grid Internal Functions
#' 
#' These functions are generally redefinitions of their counter parts in the 
#' \code{grid} package, except that they do not call \code{'L_gridDirty'}.
#' This is to enable mixing base and grid graphics in \code{\link{aheatmap}}.
#' 
#' @param fnname native \code{grid} function name.  
#' 
#' @keywords internal
#' @rdname grid 
grid.Call.graphics <- function (fnname, ...) 
{
	engineDLon <- .Call(grid:::L_getEngineDLon)
	if (engineDLon) {
		.Call.graphics(grid:::L_gridDirty, PACKAGE = "grid")
		result <- .Call.graphics(fnname, ..., PACKAGE = "grid")
	}
	else {
		#.Call(L_gridDirty)
		result <- .Call(fnname, ..., PACKAGE = "grid")
	}
	result
}

#' @rdname grid
grid.Call <- function (fnname, ...) 
{
	#.Call(L_gridDirty)
	.Call(fnname, ..., PACKAGE = "grid")
}

#' @rdname grid
current.vpPath <- function(){
	names <- NULL
	pvp <- current.viewport()
	if( is.null(pvp) ) return(NULL)
	while ( !is.null(pvp) && !grid:::rootVP(pvp)) {
		names <- c(names, pvp$name)
		pvp <- pvp$parent
	}
	if (!is.null(names)) 
		grid:::vpPathFromVector(rev(names))
	else names	
}

current.vpPath2 <- current.vpPath

## #' @rdname grid
current.viewport <- function()
{
	cv <- grid.Call(grid:::L_currentViewport)
	if( !is.null(cv) ) grid:::vpFromPushedvp(cv)
}
## 
## 
## #' @rdname grid
## vpDepth <- function() 
## {
##     pvp <- grid.Call(grid:::L_currentViewport)
##     count <- 0
##     while (!is.null(pvp$parent)) {
##         pvp <- pvp$parent
##         count <- count + 1
##     }
##     count
## }
## 
## #' @inheritParams grid::upViewport
## #' @rdname grid
## upViewport <- function (n = 1, recording = TRUE) 
## {
##     if (n < 0) 
##         stop("Must navigate up at least one viewport")
##     message('upViewport:', n)
##     on.exit(message("out upViewport"))
##     if (n == 0) {
##         n <- vpDepth()
##         upPath <- current.vpPath()
##     }
##     if (n > 0) {
##         path <- current.vpPath()
##         str(path)
##         upPath <- path[(grid:::depth(path) - n + 1):grid:::depth(path)]
##         grid.Call.graphics(grid:::L_upviewport, as.integer(n))
##         if (recording) {
##             class(n) <- "up"
##             message("upViewport to record")
##             record_grid(n)
##             message("upViewport recorded")
##         }
##     }
##     invisible(upPath)
## }
## 
## 
## #' @param x object
## #' @rdname grid
## is.vpPath <- function(x){
##     is(x, 'vpPath')
## }
## 
## #' @rdname grid
## inc.display.list <- function () 
## {
##     display.list <- grid.Call(grid:::L_getDisplayList)
##     dl.index <- grid.Call(grid:::L_getDLindex)
##     dl.index <- dl.index + 1
##     n <- length(display.list)
##     if (dl.index > (n - 1)) {
##         temp <- display.list
##         display.list <- vector("list", n + 100L)
##         display.list[1L:n] <- temp
##     }
##     grid.Call(grid:::L_setDisplayList, display.list)
##     grid.Call(grid:::L_setDLindex, as.integer(dl.index))
## }
## 
## 
## #' @rdname grid
## record_grid <- function (x) 
## {
##     grid.Call(grid:::L_setDLelt, x)
##     inc.display.list()
## }
## 
## #' @inheritParams grid::downViewport 
## #' @method downViewport default
## #' @rdname grid
## downViewport <- function (name, strict = FALSE, recording = TRUE) 
## {
##     if( !is.vpPath(name) ){
##         name <- as.character(name)
##         name <- grid:::vpPathDirect(name)
##     }
## 
##     message('downViewport:', name)
##     on.exit(message("out downViewport"))
##     if (name$n == 1) 
##         result <- grid.Call.graphics(grid:::L_downviewport, name$name, 
##                 strict)
##     else result <- grid.Call.graphics(grid:::L_downvppath, name$path, 
##                 name$name, strict)
##     pvp <- grid.Call(grid:::L_currentViewport)
##     grid.Call.graphics(grid:::L_setGPar, pvp$gpar)
##     if (recording) {
##         message("downViewport to record")
##         attr(name, "depth") <- result
##         record_grid(name)
##         message("downViewport recorded")
##     }
##     invisible(result)
## }
## 
## #' @rdname grid
## seekViewport <- function (name, recording = TRUE) 
## {
##     message("vp - seek for ", name)
##     upViewport(0, recording = recording)
##     message("vp - on top for ", name)
##     downViewport(name, recording = recording)
## }


#' \code{tryViewport} tries to go down to a viewport in the current tree.
#' 
#' \code{tryViewport}It uses \code{\link[grid]{grid.ls}} and not 
#' \code{\link{seekViewport}} as the latter would reset the graphic device.  
#' 
#' @param name viewport name 
#' @param verbose toggle verbosity
#' 
#' @rdname grid 
tryViewport <- function(name, verbose=FALSE){
	
	if( verbose ) message("vp - lookup for ", name)
	l <- grid.ls(viewports=TRUE, grobs=FALSE, print=FALSE)
	if( name %in% l$name ){
		downViewport(name)
	}	
}
