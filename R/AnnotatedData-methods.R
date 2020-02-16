# Project: NMF
# 
# Author: Renaud Gaujoux
# Created: Jan 5, 2016
###############################################################################

#' @include AnnotatedData-class.R
NULL

## AnnotationData ----
AnnotationData <- function(data, varMetadata = data.frame(row.names = colnames(data)), dimLabels = c('', 'varNames'), ...){
    new('AnnotationData', data = data, varMetadata = varMetadata, ..., dimLabels = dimLabels)
}

setMethod("colData", "AnnotationData", function(object) object@varMetadata)
setMethod("varData", "AnnotationData", function(object) colData(object))

labels.AnnotationData <- function(object, ...){
    object[['.labels']] %||% names(object)
}
setReplaceMethod("colData",
        signature=signature(
                object="AnnotationData",
                value="data.frame"),
        function(object, value) {
            idx <- match(names(value), names(object))
            varMetadata <- varData(object)[idx,,drop=FALSE]
            row.names(varMetadata) <- names(value)
            initialize(object, data=value, varMetadata=varMetadata)
        })
#' @keywords internal
setMethod("names", "AnnotationData", function(x) rownames(colData(x)))
#' @keywords internal
setReplaceMethod("names",
        signature("AnnotationData", "ANY"),
        function(x, value) 
        {
            object <- x
            if (!is.null(value) && (length(value) != dim(object@data)[[2]]))
                stop("number of new variable names (", length(value), ") ",
                        "should equal number of columns in AnnotationData (",
                        dim(object)[[2]], ")")
            if (!is.null(value)) {
                ## silently ignore attempts to set colnames to NULL
                colnames(object@data) <- value
                row.names(object@varMetadata) <- value
            }
            object
        })

setMethod("mainData", "AnnotationData", function(object) object@data)

setReplaceMethod("mainData",
        signature=signature(
                object="AnnotationData",
                value="data.frame"),
        function(object, value) {
            idx <- match(names(value), names(object))
            varMetadata <- varData(object)[idx,,drop=FALSE]
            row.names(varMetadata) <- names(value)
            initialize(object, data=value, varMetadata=varMetadata)
        })

setAs("AnnotationData", "data.frame", function(from) {
            mainData(from)
        })
#' @keywords internal
setMethod("dimnames", "AnnotationData", function(x) {
            dimnames(mainData(x))
        })
#' @keywords internal
setReplaceMethod("dimnames", "AnnotationData", function(x, value) {
            dimnames(mainData(x)) <- value
            names(colData(x)) <- value[[2]]
            x
        })
#' @keywords internal
setMethod("[",
        signature(x="AnnotationData"),
        function(x, i, j, ..., drop) {
            if (missing(drop)) drop = FALSE
            else if (drop)
                stop("'AnnotationData' does not support drop = TRUE")
            if(missing(j)) {
                mD <- x@varMetadata
                pD <- x@data[i,,drop = drop]
            } else {
                mD <- x@varMetadata[j,,drop = drop]
                if( missing( i ))
                    pD <- x@data[,j,drop = drop]
                else
                    pD <- x@data[i,j,drop = drop]
            }
            initialize(x, data=pD, varMetadata=mD)
        })

##setMethod("$", "AnnotatedDataFrame", function(x, name) `$`(pData(x), name))
#' @keywords internal
setMethod("$", "AnnotationData", function(x, name) {
            eval(substitute(colData(x)$NAME_ARG, list(NAME_ARG=name)))
        })
#' @keywords internal
setReplaceMethod("$", "AnnotationData", function(x, name, value) {
            x[[name]] <- value
            x
        })

#' @keywords internal
setMethod("[[", "AnnotationData", function(x, i, j, ...) colData(x)[[i]] )

#' @keywords internal
setReplaceMethod("[[",
        signature=signature(x="AnnotationData"),
        function(x, i, j, ..., value) {
            colData(x)[[i]] <- value
            for (metadata in names(list(...)))
                colData(x)[i, metadata] <- list(...)[[metadata]]
            x
        })


## AnnotatedData ----

setMethod("colData", "AnnotatedData", function(object) {
            slot(object, "phenoData")
        })
setReplaceMethod("colData", "AnnotatedData", function(object, value) {
             slot(object, "phenoData") <- value
             object
        })
setMethod("rowData", "AnnotatedData", function(object) {
            slot(object, "featureData")
        })
setReplaceMethod("rowData", "AnnotatedData", function(object, value) {
            slot(object, "featureData") <- value
            object
        })

#' @keywords internal
setMethod("dimnames", "AnnotatedData", function(x) {
            dimnames(mainData(x))
        })
#' @keywords internal
setReplaceMethod("dimnames", "AnnotatedData", function(x, value) {
            # apply to main data
            dimnames(mainData(x)) <- value
            # set indices on annotation data
            rowData(x)$.index <- value[[1L]]
            colData(x)$.index <- value[[2L]]
            x
        })

#' @keywords internal
setMethod("dim", "AnnotatedData", function(x) dim(mainData(x)))
#' @keywords internal
setMethod("[", "AnnotatedData", function(x, i, j, ..., drop = FALSE) {
            if (missing(drop))
                drop <- FALSE
            if (missing(i) && missing(j)) {
                if (!missing(...))
                    stop("specify rows or columns to subset; use '",
                            substitute(x), "$", names(list(...))[[1]],
                            "' to access colData variables")
                return(x)
            }
            if (!missing(j)) {
                colData(x) <- colData(x)[j,, ..., drop = drop]
            }
            if (!missing(i))
                rowData(x) <- rowData(x)[i,,..., drop=drop]
            ## assayData; implemented here to avoid function call
            mainData(x) <- {
                orig <- mainData(x)
                if (missing(i))                     # j must be present
                    orig[, j, ..., drop = drop]
                else {                              # j may or may not be present
                    if (missing(j)) orig[i,, ..., drop = drop]
                    else orig[i, j, ..., drop = drop]
                }
            }
            x
        })

## $ stops dispatching ?!
#' @keywords internal
setMethod("$", "AnnotatedData", function(x, name) {
                eval(substitute(colData(x)$NAME_ARG, list(NAME_ARG=name)))
        })

.DollarNames.AnnotatedData <- function(x, pattern)
    grep(pattern, names(colData(x)), value=TRUE)
#' @keywords internal
setReplaceMethod("$", "AnnotatedData", function(x, name, value) {
              colData(x)[[name]] = value
              x
        })
#' @keywords internal
setMethod("[[", "AnnotatedData", function(x, i, j, ...) colData(x)[[i]])
#' @keywords internal
setReplaceMethod("[[", "AnnotatedData",
        function(x, i, j, ..., value) {
            colData(x)[[i, ...]] <- value
            x
        })


