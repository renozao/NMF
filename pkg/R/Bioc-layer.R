# Layer for Bioconductor
# 
# - define methods with signature for use within Bioconductor
# - define alias methods for use in the context of microarray analysis (metagenes, metaprofiles, ...)
# 
# Author: Renaud Gaujoux \email{renaud@@cbio.uct.ac.za}
###############################################################################

#' Specific NMF Layer for Bioconductor 
#' 
#' The package NMF provides an optional layer for working with common objects
#' and functions defined in the Bioconductor platform.
#' 
#' It provides:
#' \itemize{
#' \item computation functions that support \code{ExpressionSet} objects as
#' inputs.
#' \item aliases and methods for generic functions defined and widely used by
#' Bioconductor base packages.
#' \item specialised visualisation methods that adapt the titles and legend
#' using bioinformatics terminology.
#' \item functions to link the results with annotations, etc...
#' }
#' 
#' @rdname bioc
#' @name bioc-NMF
#' 
#' @aliases nmf,ExpressionSet,ANY,ANY-method
#' @aliases nmf,matrix,ExpressionSet,ANY-method
#' 
#' @aliases seed,ExpressionSet,ANY,ANY-method
#' 
#' @aliases run,NMFStrategy,ExpressionSet,ANY-method
#' 
#' @aliases nmfModel,ExpressionSet,ANY-method
#' @aliases nmfModel,ANY,ExpressionSet-method
#' 
#' @aliases rnmf,ANY,ExpressionSet-method
#' 
#' @aliases nneg,ExpressionSet-method
#' @aliases rposneg,ExpressionSet-method
#' 
#' @aliases .atrack,ExpressionSet-method
#' 
#' @aliases sampleNames,NMF-method
#' @aliases sampleNames,NMFfitX-method
#' @aliases featureNames,NMF-method
#' @aliases featureNames,NMFfitX-method
#' 
#' @aliases nmeta
#' @aliases metagenes metagenes<-
#' @aliases metaprofiles metaprofiles<-
NULL

# add extra package Biobase
packageExtra('install_missing', 'Biobase', pkg='Biobase')

.onLoad.nmf.bioc <- function(){
	
if( pkgmaker::require.quiet('Biobase') ){

	# load Biobase package
	library(Biobase)

	#' Performs NMF on an ExpressionSet: the target matrix is the expression matrix \code{exprs(x)}.
	#' @rdname bioc
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
	
	#' Fits an NMF model partially seeding the computation with a given 
	#' ExpressionSet object passed in \code{rank}.
	#' 
	#' This method provides a shortcut for \code{nmf(x, exprs(rank), method, ...)}. 
	#' 
	#' @examples
	#' # partially seed with an ExpressionSet (requires package Biobase)
	#' \dontrun{
	#' data(esGolub)
	#' nmf(esGolub, esGolub[,1:3])
	#' }
	#' 
	setMethod('nmf', signature(x='matrix', rank='ExpressionSet', method='ANY'),
		function(x, rank, method, ...){
			# replace missing values by NULL values for correct dispatch
			if( missing(method) ) method <- NULL
			
			nmf(x, exprs(rank), method, ...)
		}
	)
	
	
	#' Seeds an NMF model directly on an ExpressionSet object.
	#' This method provides a shortcut for \code{seed(exprs(x), model, method, ...)}. 
	#' 
	#' @examples
	#' # run on an ExpressionSet (requires package Biobase)
	#' \dontrun{
	#' data(esGolub)
	#' nmf(esGolub, 3)
	#' }
	#' 
	setMethod('seed', signature(x='ExpressionSet', model='ANY', method='ANY'), 
		function(x, model, method, ...)
		{
			# replace missing values by NULL values for correct dispatch
			if( missing(method) ) method <- NULL
			if( missing(model) ) model <- NULL
			
			# apply NMF to the gene expression matrix			
			seed(Biobase::exprs(x), model, method, ...)
		}
	)
	
	#' Runs an NMF algorithm on the expression matrix of an \code{ExpressionSet} object.
	setMethod('run', signature(object='NMFStrategy', y='ExpressionSet', x='ANY'),
		function(object, y, x, ...){
			
			run(object, Biobase::exprs(y), x, ...)
			
		}
	)
		
	###% Method 'nmfModel' for 'ExpressionSet' target objects: 
	###% -> use the expression matrix of 'target' as the target matrix
	setMethod('nmfModel', signature(rank='ANY', target='ExpressionSet'),
			function(rank, target, ...){
				if( missing(rank) ) rank <- NULL
				# call nmfModel on the expression matrix
				nmfModel(rank, exprs(target), ...)
			}	
	)
	setMethod('nmfModel', signature(rank='ExpressionSet', target='ANY'),
			function(rank, target, ...){
				if( missing(target) ) target <- NULL
				# call nmfModel on the expression matrix
				nmfModel(exprs(rank), target, ...)
			}	
	)	
	
	###% Method 'rnmf' for 'ExpressionSet' target objects: 
	###% -> use the expression matrix of 'target' as the target matrix
	###% 
	setMethod('rnmf', signature(x='ANY', target='ExpressionSet'), 
		function(x, target, ...){
			rnmf(x, exprs(target), ...)
		}
	)
	
	###% The method for an \code{ExpressionSet} object returns the data.frame that 
	###% contains the phenotypic data (i.e. \code{pData(object)})
	setMethod('.atrack', 'ExpressionSet', 
		function(object, data=NULL, ...){
			if( is.null(data) ) data <- t(exprs(object))
			.atrack(pData(object), data=data, ...)	
		}
	)
	
	#' Apply \code{nneg} to the expression matrix of an \code{\link{ExpressionSet}} 
	#' object (i.e. \code{exprs(object)}). 
	#' All extra arguments in \code{...} are passed to the method \code{nneg,matrix}.
	#' 
	#' @examples
	#' 
	#' E <- ExpressionSet(x)
	#' nnE <- nneg(e)
	#' exprs(nnE)
	#' 
	setMethod('nneg', 'ExpressionSet'
			, function(object, ...){
				exprs(object) <- nneg(exprs(object), ...)
				object
			}
	)
	
	#' Apply \code{rposneg} to the expression matrix of an \code{\link{ExpressionSet}} 
	#' object (i.e. \code{exprs(object)}). 
	#' 
	#' @examples
	#' 
	#' E <- ExpressionSet(x)
	#' nnE <- posneg(E)
	#' E2 <- rposneg(nnE)
	#' all.equal(E, E2)
	#' 
	setMethod('rposneg', 'ExpressionSet'
			, function(object, ...){
				exprs(object) <- rposneg(exprs(object), ...)
				object
			}
	)
	
	###% Annotate the genes specific to each cluster.
	###%
	###% This function uses the \code{annaffy} package to generate an HTML table from the probe identifiers.
#	setGeneric('annotate', function(x, annotation, ...) standardGeneric('annotate') )
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

	## Assign BioConductor aliases
	###% number of metagenes
	nmeta <- nbasis
	###% get/set methods of basis matrix
	metagenes <- basis
	`metagenes<-` <- `basis<-`
	###% get/set methods of mixture coefficients matrix
	metaprofiles <- coef
	`metaprofiles<-` <- `coef<-`
	
	###% Get/Set methods for rows/columns names of the basis and mixture matrices
	# using the Biobase definition standard generics
	setGeneric('featureNames', package='Biobase')
	setGeneric('featureNames<-', package='Biobase')	
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
	###% For NMFfitX objects: returns the featureNames of the best fit 
	###% There is no replace method for NMFfitX objects
	setMethod('featureNames', 'NMFfitX',
		function(object){
			rownames(fit(object))
		}
	)
	
	setGeneric('sampleNames', package='Biobase')
	setGeneric('sampleNames<-', package='Biobase')	
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
	###% For NMFfitX objects: returns the sampleNames of the best fit 
	###% There is no replace method for NMFfitX objects
	setMethod('sampleNames', 'NMFfitX',
		function(object){
			colnames(fit(object))
		}
	)

	# Export layer-specific methods [only if one is loading a namespace]
	ns <- pkgmaker::addNamespaceExport(c("nmeta"
						,"featureNames", "featureNames<-"
						,"sampleNames", "sampleNames<-"
						,"metagenes", "metagenes<-"
						,"metaprofiles", "metaprofiles<-"))
	
	# return TRUE
	TRUE
}

}

