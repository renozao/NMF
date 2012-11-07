
# conditional import of class ExpressionSet from package Biobase
if( "Biobase" %in% rownames(utils::installed.packages()) ){
	if( isClass('ExpressionSet', where=getNamespace('Biobase')) ){ 
		message('Class ExpressionSet loaded from package Biobase');
		library(Biobase)
	}
}

isNamespaceLoaded <- function(name){
	!is.null(.Internal(getRegisteredNamespace(as.name(name))))
}

# store this file directory 
assign('.NMFpackageDir', getwd(), env=.GlobalEnv)

# simulate library(): package loading and attaching
nmf.library <- function(debug=FALSE){
	wd <- setwd(.NMFpackageDir)
	on.exit(setwd(wd), add=TRUE)
	message("######## Initiate package")
	if( debug )	verb <- options(verbose=TRUE)	
	.onLoad()
	.onAttach()
	if( debug )	options(verb)
	invisible()
}

###% Returns the files used in the package
pkg.collate <- function(pkg){
	if( missing(pkg) )
		pkg <- file.path(.NMFpackageDir,'..')
	
	files <- packageDescription('', pkg, fields='Collate');
	files <- strsplit(files, '[\n ]');
	files <- files[[1]];
	sapply(files, function(x) file.path(pkg, 'R', x))	
}

# Re-source the package's files
reload <- function(mode=c('none', 'comments', 'debug'), verbose=FALSE){
		
	if( isNamespaceLoaded('NMF') ){
		message('No need to load package: package loaded')
		return()
	}
	
	# load the right sourcing function
	mode <- match.arg(mode)
	source.fun <- switch(mode, none=source, comments=sourceV, debug=sourceD)
	
	# load required dependencies
	dep <- packageDescription('', file.path(.NMFpackageDir,'..'), fields='Depends')
	dep <- sub(" *([a-z_0-9]+).*", "\\1", strsplit(dep, '[\n,]')[[1]], ignore=TRUE)
	dep <- dep[dep!='R']
	message('# Load dependencies: ', paste(dep, collapse=', '), ' ... ', appendLF=FALSE)
	sapply(dep, function(p){				
				library(p, character=TRUE)
	})
	message('OK')
	# load the files as described in the DESCRIPTION file
	files <- pkg.collate()
	invisible(sapply(files, function(x){
						if( verbose ) cat("Source file: '", x, "'\n")		
						source.fun(x, chdir=TRUE)
					})
	)
	
	# simulate loading library
	nmf.library()	
}

# when sourced this file sources the package's files
reload()
