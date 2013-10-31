#' Unit Testing script for NMF package: package specific options
#'
#' @author Renaud Gaujoux
#' @creation 25 July 2009

# make the internal functions/objects visible
if( isNamespaceLoaded('NMF') ){
	#.nmf.Options.Runtime <- NMF:::.nmf.Options.Runtime
}


test.nmf.options <-function(){
	
	# clean up removing created options
	OLD <- nmf.options()
	on.exit( nmf.options(OLD) )
	
	checkTrue( is.list(nmf.options()), 'Options are returned as a list' )
	checkTrue( is.list(nmf.options('error.track', 'debug')), 'Options are returned as a list' )
	checkTrue( is.list(nmf.options('error.track')), 'Single option is returned as a list')	
	checkEquals( nmf.options('toto'), list(toto=NULL), 'Unknow in nmf.options returns NULL correctly named')
	checkTrue( !is.null(nmf.options(toto=6)), 'Can add new option')
	nmf.options(tata=10)
	opt <- nmf.options()
	checkEquals( o <- nmf.options(toto=5, tata=9, titi=25), list(toto=6, tata=10, titi=NULL), 'Changing options return old values')	
	checkIdentical({nmf.options(o); nmf.options()}, opt, "Restoring options works" )
	
	nmf.options(toto=4)
	checkEquals( nmf.options(toto=NULL), list(toto=4), 'Removing an option returns old value')
	checkTrue( !is.element('toto', names(nmf.options())), 'Removing an option actually removes it from the option list')
		
	#checkException( nmf.options(debug=NULL), 'Removing built-in option is throws an error')
}

test.nmf.getOption <-function(){
	
	# clean up removing created options
	OLD <- nmf.options()
	on.exit( nmf.options(OLD) )
	
	checkTrue( is.null(nmf.getOption('toto')), 'Unknow in nmf.getOption returns NULL')
	nmf.options(toto=5)
	checkEquals( nmf.getOption('toto'), 5, 'nmf.getOption returns correct value')
	
	# clean up removing created options
	nmf.options(toto=NULL)
}
