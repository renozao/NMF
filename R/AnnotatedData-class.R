# Project: NMF
# 
# Simplified version of the Bioconductor ExpressionSet class
# 
# Author: Renaud Gaujoux
# Created: Jan 5, 2016
###############################################################################

.AnnotationData <- setClass("AnnotationData",
        representation(varMetadata = "data.frame",
                data = "data.frame",
                dimLabels = "character"),
        prototype = prototype(
                varMetadata = new( "data.frame" ),
                data = new( "data.frame" ),
                dimLabels=c("rowNames", "columnNames")))
           
setClass("AnnotatedData",
        representation(assayData = "matrix",
                phenoData = "AnnotationData",
                featureData = "AnnotationData"),
        prototype = prototype(
                phenoData = .AnnotationData(
                        dimLabels=c("sampleNames", "sampleColumns")),
                featureData = .AnnotationData(
                        dimLabels=c("featureNames", "featureColumns"))
        )
)

# generics
setGeneric('colData', function(object) standardGeneric('colData'))
setGeneric('colData<-', function(object, value) standardGeneric('colData<-'))
#
setGeneric('mainData', function(object) standardGeneric('mainData'))
setGeneric('mainData<-', function(object, value) standardGeneric('mainData<-'))
#
setGeneric('varData', function(object) standardGeneric('varData'))
