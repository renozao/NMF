#' Unit Testing script for NMF package: virtual class NMFStrategy
#'
#' @author Renaud Gaujoux
#' @creation 14 Aug 2009

# make the internal functions/objects visible
if( isNamespaceLoaded('NMF') ){
	name <- NMF:::name
	`name<-` <- NMF:::`name<-`
	is.mixed <- NMF:::is.mixed
}

.TestSeed <- 123456

#' Unit test for constructor
test.constructor <- function(){

	# no arguments
	checkException( new('NMFStrategy'), 'Class is virtual')
	
	# define a sub-class
	setClass('A', contains='NMFStrategy')
	on.exit(removeClass('A'), add=TRUE)
	
	# with just name
		#checkException( new('A', name=character()), 'Error if slot name is a null character string')
		checkException( new('A', name=''), 'Error if slot name is an empty character string')
		
		checkTrue( validObject(new('A', name='toto')), 'No error if object with non empty name')
	
	# with objective function	
		checkException( new('A', name='toto', objective=4), 'Error if slot objective is NOT a character or function (numeric)')
		checkException( new('A', name='toto', objective=''), 'Error if slot objective is an empty character string')
		
		checkTrue( validObject(new('A', name='toto', objective='tata')), 'No error if slot objective is a non empty character string')
		checkTrue( validObject(new('A', name='toto', objective=function(){})), 'No error if slot objective is a function')
		
	# with model 
		checkException( new('A', name='toto', model=''), 'Error if slot model is an empty character string')
		checkException( new('A', name='toto', model='toto'), 'Error if slot model is not an sub-class of class NMF')
		#Now it is allowed to set the model to 'NMF'
		#checkException( new('A', name='toto', model='NMF'), 'Error if slot model is class NMF')
		
		checkTrue( validObject( new('A', name='toto', model='NMFstd')) , "Valid object if slot model is 'NMFstd'")
		checkTrue( validObject( new('A', name='toto', model='NMFns')), "Valid object if slot model is a subclass of class 'NMF'")
		checkTrue( validObject( new('A', name='toto', model=c('NMFns','NMFOffset'))), "Valid object if slot model is vector of subclasses of class 'NMF'")
		
	# with mixed set
		checkException( new('A', name='toto', mixed='toto'), 'Error if slot mixed is not logical')
		checkException( new('A', name='toto', mixed=c(TRUE, FALSE)), 'Error if slot mixed is not length 1')
		checkTrue( validObject( new('A', name='toto', mixed=TRUE) ), 'Valid object if slot mixed is TRUE' )
		checkTrue( validObject( new('A', name='toto', mixed=TRUE) ), 'Valid object if slot mixed is FALSE' )
}

check.slots.methods <- function(obj, title=''){

	checkEquals( name(obj), obj@name, paste(title, ": slot methods for 'name' give same result"))
	checkEquals( objective(obj), obj@objective, paste(title, ": slot methods for 'objective' give same result"))
	checkEquals( modelname(obj), obj@model, paste(title, ": slot methods for 'model' give same result"))
	checkEquals( is.mixed(obj), obj@mixed, paste(title, ": slot methods for 'mixed' give same result"))
	
}

#' Unit test for accessors
test.accessors <- function(){
	
	# define a sub-class
	setClass('A', contains='NMFStrategy')
	on.exit(removeClass('A'), add=TRUE)
	
	# create an object
	a <- new('A')
	
	# prototype is not valid
	checkException( validObject(a), 'Prototype object is not valid')
	check.slots.methods(a, 'Prototype object')
	
	# slot name
#	checkException( a@name <- '', "Method @<-: Error if setting slot 'name' to ''")
#	checkException( a@name <- character(), "Method @<-: Error if setting slot 'name' to character()")
#	checkException( a@name <- 4, "Method @<-: Error if setting slot 'name' to not character")
	
	checkException( name(a) <- '', "Method name<-: Error if setting slot 'name' to ''")
	checkException( name(a) <- character(), "Method name<-: Error if setting slot 'name' to character()")
	checkException( name(a) <- 4, "Method name<-: Error if setting slot 'name' to not character")
	
	name(a) <- 'toto'
	checkEquals( a@name, 'toto', "Method name<-: set slot 'name' correctly")
	check.slots.methods(a, 'Object after name set by name<-')
	
	# slot objective
	checkException( objective(a) <- '', "Method objective<-: Error if setting slot 'objective' to ''")
	checkException( objective(a) <- character(), "Method objective<-: Error if setting slot 'objective' to character()")
	checkException( objective(a) <- 4, "Method objective<-: Error if setting slot 'objective' to not character")
	
	objective(a) <- 'toto'
	checkEquals( a@objective, 'toto', "Method objective<-: set slot 'objective' correctly")
	check.slots.methods(a, 'Object after name set by objective<-')
	objective(a) <- function(...){}
	checkEquals( a@objective, function(...){}, "Method objective<-: set slot 'objective' correctly")
	check.slots.methods(a, 'Object after name set by objective<-')
}

checkClass <- function(x, cl, msg){
	checkEquals(class(x)[1L], cl, paste(msg, "[object of class '", class(x)[1L],"']"))
}

test.constructorMethod <- function(){
	
	checkClass(NMFStrategy('a', Update=function(i, y, x, ...){}), 'NMFStrategyIterative', 'With argument `Update`')
	checkClass(NMFStrategy('a', function(i, y, x, ...){}), 'NMFStrategyFunction', 'With method=function')
	checkClass(NMFStrategy('a', algorithm=function(i, y, x, ...){}), 'NMFStrategyFunction', 'With argument `algorithm`')
	
}
