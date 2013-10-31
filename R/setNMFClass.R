# Factory method for NMF model classes
# 
# Author: Renaud Gaujoux
# Creation: 18 Jul 2012
###############################################################################

#' @include NMF-class.R
NULL

if( FALSE ){ #START_DEACTIVATE
	
## #' Factory Method for NMF Model Classes
## #' 
## #' Defines two S4 classes for representing NMF models: one to hold data from 
## #' the actual model, the other one to hold fitting data for model estimated with 
## #' the function \code{\link{nmf}}. 
## #' 
#setNMFClass <- function(Class, ..., where=topns(), contains='NMFstd', VERBOSE=TRUE){
#	
#	# add 'NMF' to contains if necessary
#	wNMF <- sapply(contains, isNMFclass)
#	if( !length(wNMF) || !any(wNMF) ){
#		contains <- c(contains, 'NMFstd')
#		parentNMFClass <- 'NMFstd'
#	}else{
#		parentNMFClass <- contains[which(wNMF)]
#	}
#	
#	# extract NMF prefix if present
#	Class <- sub('^NMF(.*)', "\\1", Class)
#	# define class names
#	NMFClass <- str_c('NMF', Class)
#	NMFfitClass <- str_c(NMFClass, '_fit')
#	if( VERBOSE ){
#		message("Defining NMF classes: ", NMFClass , "(", parentNMFClass , ") and "
#				, NMFfitClass, ' in ', str_ns(where), ' ... '
#				, appendLF=FALSE)
#	}
#	# 1. Create model class
#	setClass(NMFClass, ..., where=where, contains=contains)
#	
#	# 2. Create model fit class (in the same environment as the model class)
#	e <- packageEnv(getClassDef(NMFClass)@package)
#	setClass(NMFfitClass, where=e, contains=c('NMFfit', NMFClass))
#	
#	if( VERBOSE ) message('OK')
#	
#	# return the name of the two new classes
#	c(NMFClass, NMFfitClass)
#}  

}#END_DEACTIVATE
