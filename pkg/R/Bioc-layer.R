# Layer for Bioconductor
# 
# - define methods with signature for use within Bioconductor
# - define alias methods for use in the context of microarray analysis (metagenes, metaprofiles, ...)
# 
# Author: Renaud Gaujoux \email{renaud@@cbio.uct.ac.za}
###############################################################################


#if( isClass('ExpressionSet') ){
.init.nmf.bioc <- function(){
	
	bioc.loaded <- FALSE
	if( "Biobase" %in% rownames(utils::installed.packages()) ){

	# load Biobase package
	library(Biobase)
	
	#' BioConductor alias to nbasis
	if ( !isGeneric("nmeta")) setGeneric('nmeta', function(object, ...) standardGeneric('nmeta') )
	setMethod('nmeta', signature(object='NMF'),
		function(object){
			nbasis(object)
		}
	)
	
	#' BioConductor alias get/set methods of basis matrix
	if ( !isGeneric("metagenes")) setGeneric('metagenes', function(object, ...) standardGeneric('metagenes') )
	setMethod('metagenes', signature(object='NMF'),
		function(object, ...){
			basis(object, ...)
		}
	)
	if ( !isGeneric("metagenes<-") ) setGeneric('metagenes<-', function(object, value) standardGeneric('metagenes<-') )
	setReplaceMethod('metagenes', signature(object='NMF', value='matrix'), 
		function(object, value){ 
			basis(object) <- value
			object
		} 
	)
	
	#' BioConductor alias to get/set methods of mixture coefficients matrix
	if ( !isGeneric("metaprofiles")) setGeneric('metaprofiles', function(object, ...) standardGeneric('metaprofiles') )
	setMethod('metaprofiles', signature(object='NMF'),
		function(object, ...){
			coef(object, ...)
		}
	)
	if ( !isGeneric("metaprofiles<-") ) setGeneric('metaprofiles<-', function(object, value) standardGeneric('metaprofiles<-') )
	setReplaceMethod('metaprofiles', signature(object='NMF', value='matrix'), 
		function(object, value){ 
			coef(object) <- value
			object
		} 
	)

	#' Performs NMF on an ExpressionSet: the target matrix is the expression matrix \code{exprs(x)}.
	setMethod('nmf', signature(x='ExpressionSet', rank='ANY', method='ANY'), 
		function(x, rank, method, ...)
		{
			# apply NMF to the gene expression matrix	
			nmf(Biobase::exprs(x), rank, method, ...)
		}
	)
	
	#' Run the algorithm on the expression matrix of an \code{ExpressionSet} object.
	setMethod('run', signature(object='NMFStrategy', target='ExpressionSet', start='ANY'),
		function(object, target, start, ...){
			
			run(object, Biobase::exprs(target), start, ...)
			
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
	
	# Export layer-specific methods [only if one is loading a namespace] 
	is.loading <- try(info <- loadingNamespaceInfo(), silent=TRUE)
	if( !is(is.loading, 'try-error') ){
		ns <- .Internal(getRegisteredNamespace(as.name(info$pkgname)))
		if ( is.null(ns) )
			stop("Error in exporting NMF-BioConductor layer: cannot find the loading namespace's environment");
		
		# export the methods into the loading namespace
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
			rownames(basis(object))
		}
)
setReplaceMethod('featureNames', 'NMF',
		function(object, value){
			rownames(basis(object)) <- value
			return(object)
		}
)
#' For NMFSet objects: returns the featureNames of the best fit 
#' There is no replace method for NMFSet objects
setMethod('featureNames', 'NMFSet',
		function(object){
			featureNames(fit(object))
		}
)

setMethod('sampleNames', 'NMF',
	function(object){
		colnames(coef(object))
	}
)
setReplaceMethod('sampleNames', 'NMF',
	function(object, value){
		colnames(coef(object)) <- value
		return(object)
	}
)
#' For NMFSet objects: returns the sampleNames of the best fit 
#' There is no replace method for NMFSet objects
setMethod('sampleNames', 'NMFSet',
		function(object){
			sampleNames(fit(object))
		}
)

	# return the result
	return(bioc.loaded)
}

