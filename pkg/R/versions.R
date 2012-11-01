# Tracking/Updating S4 class versions
# 
# Author: Renaud Gaujoux
###############################################################################

#' @include utils.R
NULL

objectUpdater <- local({
	.REGISTRY <- list()
	
	function(x, version=NULL, fun=NULL, vfun=NULL, verbose=FALSE){
		
		if( missing(x) ) return( .REGISTRY )
		
		if( is.null(version) ){
			cl <- class(x)
			UPDATER <- .REGISTRY[[cl]]
			
			vmsg <- 'Class' 
			if( is.character(verbose) ){
				vmsg <- paste(verbose, ':', sep='')
				verbose <- TRUE
			} 
			if( verbose ) message("# ", vmsg, " '", cl, "' ... ", appendLF=FALSE)
			if( !isS4(x) ){
				if( verbose) message("NO")
				return(x)
			}
			
			# create new object from old slots
			newObject <- if( verbose ){
				message()
				updateObjectFromSlots(x, verbose=verbose>1) 
			}else suppressWarnings( updateObjectFromSlots(x, verbose=verbose>1) )
				
			
			if( is.null(UPDATER) ){
				if( verbose) message("AUTO")
				return(newObject)
			}
			
			# find object version
			v <- sapply(UPDATER, function(f) f$vfun(x))
			v <- which(v)
			if( !length(v) ){
				if( verbose) message("SKIP [version unknown]")
				return(newObject)
			}
			if( length(v) > 1L ){
				warning("Object matched multiple version of class '", cl
						, "' [", paste(names(UPDATER)[v], collapse=", "), "]")
				if( verbose) message("SKIP [multiple versions]")
				return(newObject)
			}else if( verbose) message("UPDATING [", appendLF=FALSE)
			
			for(n in names(UPDATER[v[1L]])){
				f <- UPDATER[[n]]
				if( verbose ) message(n, ' -> ', appendLF=FALSE)
				newObject <- f$fun(x, newObject)
			}
			if( verbose ) message("*]")
			# return updated object
			return(newObject)
		}
		
		stopifnot( is.character(x) )
		if( is.null(version) ){
			if( !is.null(fun) || !is.null(vfun) )
				stop("Argument `version` is required for defining updater functions for class `", x, "`")
			return(.REGISTRY[[x]])
		}
		
		if( is.null(.REGISTRY[[x]]) ) .REGISTRY[[x]] <<- list() 
		# check result is a function
		stopifnot(is.function(fun))
		stopifnot(is.function(vfun))
		if( !is.null(.REGISTRY[[x]][[version]]) )
			stop("An update for class '", x, "' version ", version, " is already defined")
		
		.REGISTRY[[x]][[version]] <<- list(vfun=vfun, fun=fun)
		# put updaters in order
		.REGISTRY[[x]] <<- .REGISTRY[[x]][orderVersion(names(.REGISTRY[[x]]))]
		invisible(.REGISTRY[[x]])
	}
})

# Taken from BiocGenerics 2.16
getObjectSlots <- function (object) 
{
	if (!is.object(object) || isVirtualClass(class(object))) 
		return(NULL)
	value <- attributes(object)
	value$class <- NULL
	if (is(object, "vector")) {
		.Data <- as.vector(object)
		attr(.Data, "class") <- NULL
		attrNames <- c("comment", "dim", "dimnames", "names", 
				"row.names", "tsp")
		for (nm in names(value)[names(value) %in% attrNames]) attr(.Data, 
					nm) <- value[[nm]]
		value <- value[!names(value) %in% attrNames]
		value$.Data <- .Data
	}
	value
}

# Taken from BiocGenerics 2.16
updateObjectFromSlots <- function (object, objclass = class(object), ..., verbose = FALSE) 
{
	updateObject <- nmfObject
	if (is(object, "environment")) {
		if (verbose) 
			message("returning original object of class 'environment'")
		return(object)
	}
	classSlots <- slotNames(objclass)
	if (is.null(classSlots)) {
		if (verbose) 
			message("definition of '", objclass, "' has no slots; ", 
					"returning original object")
		return(object)
	}
	errf <- function(...) {
		function(err) {
			if (verbose) 
				message(..., ":\n    ", conditionMessage(err), 
						"\n    trying next method...")
			NULL
		}
	}
	if (verbose) 
		message("updateObjectFromSlots(object = '", class(object), 
				"' class = '", objclass, "')")
	objectSlots <- getObjectSlots(object)
	nulls <- sapply(names(objectSlots), function(slt) is.null(slot(object, 
								slt)))
	objectSlots[nulls] <- NULL
	joint <- intersect(names(objectSlots), classSlots)
	toUpdate <- joint[joint != ".Data"]
	objectSlots[toUpdate] <- lapply(objectSlots[toUpdate], updateObject, 
			..., verbose = verbose)
	toDrop <- which(!names(objectSlots) %in% classSlots)
	if (length(toDrop) > 0L) {
		warning("dropping slot(s) ", paste(names(objectSlots)[toDrop], 
						collapse = ", "), " from object = '", class(object), 
				"'")
		objectSlots <- objectSlots[-toDrop]
	}
	res <- NULL
	if (is.null(res)) {
		if (verbose) 
			message("heuristic updateObjectFromSlots, method 1")
		res <- tryCatch({
					do.call(new, c(objclass, objectSlots[joint]))
				}, error = errf("'new(\"", objclass, "\", ...)' from slots failed"))
	}
	if (is.null(res)) {
		if (verbose) 
			message("heuristic updateObjectFromSlots, method 2")
		res <- tryCatch({
					obj <- do.call(new, list(objclass))
					for (slt in joint) slot(obj, slt) <- updateObject(objectSlots[[slt]], 
								..., verbose = verbose)
					obj
				}, error = errf("failed to add slots to 'new(\"", objclass, 
						"\", ...)'"))
	}
	if (is.null(res)) 
		stop("could not updateObject to class '", objclass, "'", 
				"\nconsider defining an 'updateObject' method for class '", 
				class(object), "'")
	res
}


#' Updating NMF Objects
#' 
#' This function serves to update an objects created with previous versions of the 
#' NMF package, which would otherwise be incompatible with the current version, 
#' due to changes in their S4 class definition.
#' 
#' This function makes use of heuristics to automatically update object slots, 
#' which have been borrowed from the BiocGenerics package, the function 
#' \code{updateObjectFromSlots} in particular. 
#'
#' @param object an R object created by the NMF package, e.g., an object of class 
#' \code{\linkS4class{NMF}} or \code{\linkS4class{NMFfit}}.
#' @param verbose logical to toggle verbose messages.
#' 
#' @export
#' 
nmfObject <- function(object, verbose=FALSE){
	
	objectUpdater(object, verbose=verbose)
	
}
