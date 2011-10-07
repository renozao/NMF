# Layer for Bioconductor
# 
# - define methods with signature for use within Bioconductor
# - define alias methods for use in the context of microarray analysis (metagenes, metaprofiles, ...)
# 
# Author: Renaud Gaujoux \email{renaud@@cbio.uct.ac.za}
###############################################################################


.onLoad.nmf.bioc <- function(){
	
if( !"Biobase" %in% rownames(utils::installed.packages()) )
	FALSE
else{

	# load Biobase package
	library(Biobase)

	###% Performs NMF on an ExpressionSet: the target matrix is the expression matrix \code{exprs(x)}.
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
	
	###% Seeds an NMF model from an ExpressionSet: the target matrix is the expression matrix \code{exprs(x)}.
	setMethod('seed', signature(x='ExpressionSet', model='ANY', method='ANY'), 
		function(x, model, method, ...)
		{
			# replace missing values by NULL values for correct dispatch
			if( missing(method) ) method <- nmf.getOption('default.seed')
			
			# apply NMF to the gene expression matrix			
			seed(Biobase::exprs(x), model, method, ...)
		}
	)
	
	###% Run the algorithm on the expression matrix of an \code{ExpressionSet} object.
	setMethod('run', signature(method='NMFStrategy', x='ExpressionSet', seed='ANY'),
		function(method, x, seed, ...){
			
			run(method, Biobase::exprs(x), seed, ...)
			
		}
	)
	
	###% Computes the distance between the target ExpressionSet and its NMF fit 
	setMethod('distance', signature(target='ExpressionSet', x='NMF'), 
			function(target, x, ...){
								
				# compute the distance between the expression matrix and the fitted NMF model
				distance(Biobase::exprs(target), x, ...)
			}
	)
	
	###% Method 'nmfModel' for 'ExpressionSet' target objects: 
	###% -> use the expression matrix of 'target' as the target matrix
	setMethod('nmfModel', signature(rank='ANY', target='ExpressionSet'),
			function(rank, target, ...){
				# call nmfModel on the expression matrix
				nmfModel(rank, exprs(target), ...)
			}	
	)
	setMethod('nmfModel', signature(rank='ExpressionSet', target='ANY'),
			function(rank, target, ...){
				# call nmfModel on the expression matrix swapping the arguments
				nmfModel(target, exprs(rank), ...)
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
	setMethod('.atrack', 'ExpressionSet', function(object, ...) pData(object) )
	
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

	# define generic for the rows/columns names, using the Biobase definition
	setGeneric('featureNames', package='Biobase')
	setGeneric('featureNames<-', package='Biobase')
	setGeneric('sampleNames', package='Biobase')
	setGeneric('sampleNames<-', package='Biobase')
		
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
	ns <- getLoadingNamespace(getenv=TRUE)
	if( !is.null(ns) ){
		namespaceExport(ns, c("nmeta"
						,"featureNames"
						,"featureNames<-"
						,"sampleNames"
						,"sampleNames<-"
						,"metagenes"
						,"metagenes<-"
						,"metaprofiles"
						,"metaprofiles<-")
		)
	}
	
	# return TRUE
	TRUE
}

}

