## Base abstract class that defines a strategy.
#setClass('Strategy'
#	, representation(
#			name = 'character' # the strategy name
#			, package = 'character' # the package that defines the strategy
#	)
#	, contains = 'VIRTUAL'
#	, validity=function(object){
#		
#		# slot 'name' must be a non-empty character string
#		obj <- name(object)
#		if( length(obj) != 0 && (!is.character(obj) || length(obj)!=1 || obj=='') )
#			return("Slot 'name' must be a non-empty character string.")
#		
#	}
#)
#
#setMethod('initialize', 'Strategy',
#	function(.Object, ..., name, package){
#				
#		print(.Object)
#		.Object <- callNextMethod(.Object, ...)
#		print(.Object)
#		
#		# set the name of the strategy
#		if( !missing(name) )
#			.Object@name <- name
#		
#		# set the name to the loading package		
#		if( missing(package) ){
#			ns <- topenv()
#			package <- if( FALSE && isNamespace(ns) )
#				package <- getNamespaceInfo(ns, 'spec')$name
#			else
#				''
#		}
#		
#		.Object@package <- package
#		print(.Object)
#		
#		.Object
#	}
#)
#
## Accessor methods to slot \code{name}
#setGeneric('name', function(object, ...) standardGeneric('name'))
#setMethod('name', signature(object='Strategy'),
#	function(object){
#		slot(object, 'name')
#	}
#)
#setGeneric('name<-', function(object, ..., value) standardGeneric('name<-'))
#setReplaceMethod('name', signature(object='Strategy', value='character'),
#	function(object, value){
#		slot(object, 'name') <- value
#		validObject(object)
#		object
#	}
#)

###% Base abstract class that defines the interface for NMF algorithms.
###%
###% @aliases NMFStrategy-class
###% 
###% @slot name character string giving the name of the strategy
###%
###% @slot objective the objective function associated with the algorithm (Frobenius, Kullback-Leibler, etc...). 
###% It is either a character string as a key registered by \code{nmfRegisterDistance} or a function definition. 
###% In the latter case, the given function must have the following signature (x=matrix, y=matrix) and return a nonnegative real value.
###%
###% @slot model a character string giving either the (sub)class name of the NMF-class instance used and returned by the strategy, or a function name.
###% In the latter case, the given function must have the following signature \code{(v=matrix, r=integer, ...)}, where the \code{v} is the target
###% matrix to approximate and \code{r} is the rank of factorization to achieve.
###%
setClass('NMFStrategy'
	, representation(
				name = 'character' # name of the method (also key)
				, objective = '.functionSlot' # the objective function used to compute the error (defined by name or function)
				, model = 'character' # NMF model to use
				, mixed = 'logical' # can the input data be negative?
	)
	, prototype=prototype(name='', objective='euclidean', model='NMFstd', mixed=FALSE)
	, validity=function(object){
		
		# slot 'name' must be a non-empty character string
		obj <- name(object)
		if( !is.character(obj) || length(obj)!=1 || obj=='' )
			return("Slot 'name' must be a non-empty character string.")
			
		# slot 'objective' must either be a non-empty character string or a function
		obj <- objective(object)
		if( is.character(obj) && obj == '' )
			return("Slot 'objective' must either be a non-empty character string or a function definition.")
			
		# slot 'model' must be the name of a class that extends class 'NMF'
		obj <- modelname(object)
		if( !is.character(obj) )
			return("Slot 'model' must be a character vector")
		if( any(inv <- !sapply(obj, is.nmf)) )
			return(paste("Slot 'model' must contain only names of a class that extends class 'NMF' [failure on class(es) "
					, paste( paste("'", obj[inv], "'", sep=''), collapse=', ')  
					,"]"
					, sep=''))
		
		# slot 'mixed' must be a single logical		
		obj <- slot(object, 'mixed')
		if( length(obj) != 1 )
			return( paste("Slot 'mixed' must be a single logical [length=", length(obj), "]", sep='') )
	}
	, contains = 'VIRTUAL'
)

setMethod('show', 'NMFStrategy',
		function(object){			
			cat('<object of class: ', class(object), ">\n")
			cat("name:\t", name(object), "\n")
			svalue <- objective(object)
			svalue <- if( is.function(svalue) ) '<function>' else paste("'", svalue,"'", sep='')
			cat("objective:\t", svalue, "\n")
			cat("NMF model:\t", modelname(object), "\n")
			return(invisible())
		}
)

# Coerce method for 'NMFStrategy' objects into 'character': give the main name
setAs('NMFStrategy', 'character'
	, def = function(from) name(from)	
) 

###% Main interface to run the algorithm
setGeneric('run', function(method, x, seed, ...) standardGeneric('run'))
setMethod('run', signature(method='NMFStrategy', x='matrix', seed='NMFfit'),
	function(method, x, seed, ...){
		stop("NMFStrategy::run is a pure virtual method that should be overloaded in class '", class(method),"'.")
	}
)

###% Accessor methods to slot \code{name}
setGeneric('name', function(object, ...) standardGeneric('name'))
setMethod('name', signature(object='NMFStrategy'),
	function(object, all=FALSE){
		if( !all ) slot(object, 'name')[1] else slot(object, 'name')
	}
)
setGeneric('name<-', function(object, ..., value) standardGeneric('name<-'))
setReplaceMethod('name', signature(object='NMFStrategy', value='character'),
	function(object, value){
		slot(object, 'name') <- value
		validObject(object)
		object
	}
)

###% Accessor methods to slot \code{objective}
setMethod('objective', signature(object='NMFStrategy'),
	function(object, x, y, ...){
	
		obj.fun <- slot(object, 'objective')
		
		# when both x and y are missing then returns slot objective
		if( missing(x) && missing(y) ) return(obj.fun)
		
		# return the distance computed using the strategy's objective function
		if( !is.function(obj.fun) )
			distance(x, y, method=obj.fun, ...)
		else # directly compute the objective function
			obj.fun(x, y, ...)
		
	}
)
setReplaceMethod('objective', signature(object='NMFStrategy', value='character'),
	function(object, value){
		#TODO: test for the existence of objective method
		slot(object, 'objective') <- value
		validObject(object)
		object
	}
)
setReplaceMethod('objective', signature(object='NMFStrategy', value='function'),
	function(object, value){
		slot(object, 'objective') <- value
		validObject(object)
		object
	}
)

###% Accessor methods to slot \code{model}
setMethod('modelname', signature(object='NMFStrategy'),
	function(object){
		slot(object, 'model')
	}
)
setGeneric('modelname<-', function(object, ..., value) standardGeneric('modelname<-'))
setReplaceMethod('modelname', signature(object='NMFStrategy', value='character'),
	function(object, value){
		slot(object, 'model') <- value
		validObject(object)
		object
	}
)

###% Accessor methods to slot \code{mixed}
setGeneric('is.mixed', function(object, ...) standardGeneric('is.mixed'))
setMethod('is.mixed', signature(object='NMFStrategy'),
		function(object){
			return( slot(object, 'mixed') )
		}
)

###########################################################################
# REGISTRY METHODS FOR ALGORITHMS
###########################################################################

###% Register a new algorithm into the NMF registry.
###%
setGeneric('nmfRegisterAlgorithm', function(method, key, ...) standardGeneric('nmfRegisterAlgorithm') )

###% Define an alias
setMethodNMF <- nmfRegisterAlgorithm

setMethod('nmfRegisterAlgorithm', signature(method='NMFStrategy', key='character'),
		function(method, key, ...){
			
			original.name <- name(method)
			fname <- key[1]
			if( original.name != fname )
				message("Registering method '", fname,"' based on template '", original.name, "'")
			
			# reset slots if necessary
			dots <- list(...)
			if( any(ndots <- names(dots) != '') > 0 ){
				dots <- dots[ndots]
				# partially match the slots
				mslots <- charmatch(names(dots), slotNames(class(method)))
				
				# if necessary: create object using 'method' as a template and the new values for the slots 
				if( length(dot.slots <- dots[!is.na(mslots)]) > 0 ){
					method <- do.call(new, c(list(class(method), method), dot.slots))
					# valid the new strategy
					validObject(method)
				}
			}
			
			# force the name to match the registering key, and warn the user if necessary			
			if( !is.null(dots$name) && dots$name != fname )
				warning("Argument 'name' discarded: the NMF algorithm's name was forced to match its key ['",fname,"'].")							
			name(method) <- fname
			
			# add to the algorithm registry 
			nmfRegister(method, key, registry.name='algorithm', ...)
		}
)
setMethod('nmfRegisterAlgorithm', signature(method='NMFStrategy', key='missing'), 
		function(method, key, ...){
			
			# use the name as a key
			key <- name(method)
			
			# register the method
			nmfRegisterAlgorithm(method, key, ...)		
		}
)
setMethod('nmfRegisterAlgorithm', signature(method='function', key='character'), 
	function(method, key, overwrite=FALSE, save=FALSE, ...){
		
		# build the NMFStrategy
		strategy <- newNMFStrategy(method, key, ...)
		
		# register the method
		nmfRegisterAlgorithm(strategy, key, overwrite, save, ...)
	}
)

setMethod('nmfRegisterAlgorithm', signature(method='character', key='ANY'), 
		function(method, key, ...){
			
			# register/update the already registered method
			method <- nmfAlgorithm(method)
			if( missing(key) )
				key <- name(method)
			nmfRegisterAlgorithm(nmfAlgorithm(method), key, ...)
		}
)


###% Factory method to create NMFStrategy objects.
###%
###% Create predefined NMFStrategy objects that implement algorithms from different papers.
setGeneric('newNMFStrategy', function(method, key, ...) standardGeneric('newNMFStrategy') )
setMethod('newNMFStrategy', signature(method='function', key='character'), 
	function(method, key, ...){
			
		# build a NMFStrategyFunction object on the fly to wrap function 'method'
		strategy.params <- list('NMFStrategyFunction', name=key, algorithm=method)
		strategy <- do.call('new', c(strategy.params, list(...)))
		
		# valid the new strategy
		validObject(strategy)
		
		# register the method
		strategy
	}
)

###% tells if an NMFStrategy can fit a given NMF model
setGeneric('canFit', function(x, y, ...) standardGeneric('canFit') )

###% tells if an NMFStrategy can fit a given NMF model given by its name
setMethod('canFit', signature(x='NMFStrategy', y='character'),
	function(x, y, exact=FALSE){
		
		if( !exact ){
			
			# check for one model amongst all the models fittable by the strategy
			can <- if( length(mo <- modelname(x)) > 1 )
						sapply(mo, function(m) extends(y, m))
					else extends(y, mo)
			any(can)
			
		}else
			is.element(y, modelname(x))
	}
)
###% tells if an NMFStrategy can fit a given NMF model given by an instance of it
setMethod('canFit', signature(x='NMFStrategy', y='NMF'),
		function(x, y, ...){
			canFit(x, modelname(y), ...)
		}
)
###% tells if an NMFStrategy given by its name can fit a given NMF model
setMethod('canFit', signature(x='character', y='ANY'),
		function(x, y, ...){
			canFit(nmfAlgorithm(x), y, ...)
		}
)

###% try to select the correct NMF algorithm given the NMF model
selectMethodNMF <- function(model, load=FALSE, exact=FALSE, all=FALSE, warning=TRUE){
	
	# lookup for all the algorithms that can fit the given model
	#NB: if only one model needs to be selected then first look for an exact fit as 
	# this would need to be done with exact=FALSE and TRUE anyways
	w <- sapply(algo <- nmfAlgorithm(), canFit, model, exact= if(all) exact else TRUE)	
	algo <- algo[w]
	
	# if no suitable algorithm was found, and an exact match is not required 
	# then look for other potential non-exact algorithms
	if( !all && !exact && length(algo) == 0 ){
		w <- sapply(algo <- nmfAlgorithm(), canFit, model, exact=FALSE)
		algo <- algo[w]
	}
	
	# return NULL if no algorithm was found
	if( length(algo) == 0 )		
		return(NULL)
			
	# if all=FALSE then try to choose the default algorithm if present in the list, or the first one
	res <- if( !all && length(algo) > 1 ){
		
		if( is.element(default <- nmf.getOption('default.algorithm'), algo) )
			default
		else{
			if( warning ) 
				warning("Non default NMF algorithm '", algo[1], "' was selected against other possible algorithm(s): "
						, paste(paste("'", algo[-1], "'", sep=''), collapse=", "))
			algo[1]
		}
		
	}else # otherwise return all the algorithms
		algo
	
	# load the methods if required
	if( load ){
		if( length(res) > 1 ) sapply(res, nmfAlgorithm) else nmfAlgorithm(res)
	}
	else
		res	
}

###% Access to registered algorithms
nmfAlgorithm <- function(name=NULL, model, type=c('C', 'R'), ...){	
	
	# if one passes an NMFStrategy just returns it
	if( is(name, 'NMFStrategy') )
		return(name)
	# check for type filtering
	if( missing(name) && !missing(type)  ){
		type <- match.arg(type)
		algo <- nmfGet(name, registry.name='algorithm', all=TRUE)
		algo <- algo[grep(paste("^\\.", type, '#', sep=''), algo)]
		names(algo) <- sub(paste(".", type, '#', sep=''), '', algo, fixed=TRUE)
	}
	else algo <- nmfGet(name, registry.name='algorithm', ...)
	
	# if no model was specified then return the selected algorithm(s)
	if( missing(model) )
		algo
	else{ # lookup for an algorithm suitable for the given NMF model
		if( !is.character(model) || (length(model) != 1) || nchar(model) == 0 || !extends(model, 'NMF') )
			stop("argument 'model' must be the name of a class that extends class 'NMF'")
		
		# if a single algorithm was found then check if it can fit the given model before returning it
		if( length(algo) == 1 ){
			if( !canFit(algo, model) )
				stop("NMF algorithm '", as(algo, 'character'), "' cannot fit model '", model, "'.")
			
			return( algo )
		}
		
		# returns the algorithms that are defined to fit the given model
		selectMethodNMF(model, exact=TRUE)
		
	}
}

###% Returns TRUE if the algorithm is registered FALSE otherwise
existsNMFAlgorithm <- function(name, exact=TRUE){	
	
	!is.null( nmfGet(name, registry.name='algorithm', all=TRUE, error=FALSE, exact=exact) )
	
}
