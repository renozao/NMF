# Layer for Bioconductor
# 
# - define methods with signature for use within Bioconductor
# - define alias methods for use in the context of microarray analysis (metagenes, metaprofiles, ...)
# 
# Author: Renaud Gaujoux \email{renaud@@cbio.uct.ac.za}
###############################################################################


.onLoad.nmf.bioc <- function(){
	
	bioc.loaded <- FALSE
	if( "Biobase" %in% rownames(utils::installed.packages()) ){

	# load Biobase package
	library(Biobase)

	#' Performs NMF on an ExpressionSet: the target matrix is the expression matrix \code{exprs(x)}.
	setMethod('nmf', signature(x='ExpressionSet', rank='ANY', method='ANY'), 
		function(x, rank, method, ...)
		{
			# replace missing values by NULL values for correct dispatch
			if( missing(method) ) method <- NULL
			if( missing(rank) ) rank <- NULL
			
			# apply NMF to the gene expression matrix			
			nmf(Biobase::exprs(x), rank, method, ...)
		}
	)
	
	#' Run the algorithm on the expression matrix of an \code{ExpressionSet} object.
	setMethod('run', signature(method='NMFStrategy', x='ExpressionSet', seed='ANY'),
		function(method, x, seed, ...){
			
			run(method, Biobase::exprs(x), seed, ...)
			
		}
	)
	
	#' Computes the distance between the target ExpressionSet and its NMF fit 
	setMethod('distance', signature(target='ExpressionSet', x='NMF'), 
			function(target, x, ...){
								
				# compute the distance between the expression matrix and the fitted NMF model
				distance(Biobase::exprs(target), x, ...)
			}
	)
	
	#' Annotate the genes specific to each cluster.
	#'
	#' This function uses the \code{annaffy} package to generate an HTML table from the probe identifiers.
#	if ( is.null(getGeneric("annotate")) ) setGeneric('annotate', function(x, annotation, ...) standardGeneric('annotate') )
#	setMethod('annotate', signature(x='factor', annotation='character'), 
#		function(x, annotation, filename='NMF genes', outdir='.', name='Cluster specific genes', ...)
#		{
#			library(annaffy)
#			anncols<-aaf.handler()[c(1:3, 6:13)]			
#			
#			# add html suffix to filename if necessary
#			if( length(grep("\\.html$", filename)) == 0 ) filename <- paste(filename, 'html', sep='.')
#			
#			# for each cluster annotate the genes set		
#			print(head(x))		
#			by(names(x), x, function(g){	
#						print(head(g))
#						if( length(g) == 0 ) return()
#						g <- as.character(g)
#						anntable <- aafTableAnn(g, annotation, anncols)
#						# generate HTML output
#						saveHTML(anntable, file.path(outdir,filename), title=paste(name, '[top', nrow(anntable),']'))
#					}, simplify=FALSE)
#			
#			# return nothing
#			invisible()
#		}
#	)
#	
#	setMethod('annotate', signature(x='NMF', annotation='character'), 
#		function(x, annotation, ...)
#		{
#			s <- extractFeatures(x)
#			class <- .predict.nmf(t(s))
#			annotate(class, annotation=annotation, ...)
#		}
#	)

	# define generic for the rows/columns names, using the Biobase definition
	setGeneric('featureNames', package='Biobase')
	setGeneric('featureNames<-', package='Biobase')
	setGeneric('sampleNames', package='Biobase')
	setGeneric('sampleNames<-', package='Biobase')
		
	## Assign BioConductor aliases
	#' number of metagenes
	nmeta <- nbasis
	#' get/set methods of basis matrix
	metagenes <- basis
	`metagenes<-` <- `basis<-`
	#' get/set methods of mixture coefficients matrix
	metaprofiles <- coef
	`metaprofiles<-` <- `coef<-`
	
	# Export layer-specific methods [only if one is loading a namespace]
	ns <- getLoadingNamespace(getenv=TRUE)
	if( !is.null(ns) ){		
		namespaceExport(ns, c("nmeta"
						,"metagenes"
						,"metagenes<-"
						,"metaprofiles"
						,"metaprofiles<-")
		)
	}
	
	
	# set result to TRUE
	bioc.loaded <- TRUE
}else{
	# define generic for the rows/columns names, following the Biobase definition (no '...')
	setGeneric('featureNames', function(object) standardGeneric('featureNames') )
	setGeneric('featureNames<-', function(object, value) standardGeneric('featureNames<-') )
	setGeneric('sampleNames', function(object) standardGeneric('sampleNames') )
	setGeneric('sampleNames<-', function(object, value) standardGeneric('sampleNames<-') )
}


#' Get/Set methods for rows/columns names of the basis and mixture matrices
setMethod('featureNames', 'NMF',
	function(object){
		rownames(object)
	}
)
setReplaceMethod('featureNames', 'NMF',
	function(object, value){
		rownames(object) <- value
		object
	}
)
#' For NMFfitX objects: returns the featureNames of the best fit 
#' There is no replace method for NMFfitX objects
setMethod('featureNames', 'NMFfitX',
	function(object){
		rownames(fit(object))
	}
)

setMethod('sampleNames', 'NMF',
	function(object){
		colnames(object)
	}
)
setReplaceMethod('sampleNames', 'NMF',
	function(object, value){
		colnames(object) <- value
		object
	}
)
#' For NMFfitX objects: returns the sampleNames of the best fit 
#' There is no replace method for NMFfitX objects
setMethod('sampleNames', 'NMFfitX',
	function(object){
		colnames(fit(object))
	}
)

	# return the result
	bioc.loaded
}

