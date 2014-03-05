#' Unit Testing script for NMF package: NMF interface for algorithms.
#'
#' @author Renaud Gaujoux
#' @creation 14 May 2009

library(rngtools)
checkIdenticalRNG <- checkRNG

# make the internal functions/objects visible
if( isNamespaceLoaded('NMF') ){
	seed <- NMF:::seed
#	nmfUnregister <- NMF:::nmfUnregister	
	name <- NMF:::name
	`name<-` <- NMF:::`name<-`
}

.testData <- function(n=20, r=3, m=10, ...){
	syntheticNMF(n, r, m, ...)
}

test.registry <- function(){

	checkNotNull <- function(x, ...) checkTrue(!is.null(x), ...)
	# register function as a method
	dummy.method <- function(){}
		
		# methods that don't exist
#		checkException(nmfUnregister('algo.tata'), 'Unregister a method without specifying the registry name')
#		checkTrue(nmfUnregister('algo.tata', 'algorithm'), 'Unregister a method that does not exist: should not generate an error')
		checkIdentical(removeNMFMethod('algo.tata'), FALSE, 'removeNMFMethod a method that does not exist: should not generate an error')
		checkException( nmfAlgorithm('algo.toto'), 'Try to access a method that does not exist: should generate an error')
		checkTrue(is.null(nmfAlgorithm('algo.toto', error=FALSE)), 'Try to access a method that does not exist with error=FALSE: should NOT generate an error and return NULL')
		
	# Registration of new methods		
		# force un-registration of 'dummy' on exit
		on.exit({removeNMFMethod('dummy')}, add=TRUE)
		checkNotNull(setNMFMethod('dummy', dummy.method), 'Register works on dummy -- empty -- method')
		checkException(setNMFMethod('dummy', dummy.method), 'Try to register an algorithm with an existing name')
		checkNotNull(setNMFMethod('dummy', dummy.method, overwrite=TRUE), 'Overwrite an existing algorithm ')		
		
	# Access to methods
		checkTrue( is(nmfAlgorithm('dummy'), 'NMFStrategyFunction'), 'Get method by exact match')
		checkTrue( is(nmfAlgorithm('dum'), 'NMFStrategyFunction'), 'Get method by partial match')
		checkEquals( name(nmfAlgorithm('dum')),  'dummy', "The method's full name is set in slot 'name'")
}


#' Utility function for \code{test.seed}: performs a set of test on a seeded object
check.seed <- function(title, obj, V, r, expect.class){
		
	checkTrue( isNMFfit(obj), paste(title, ": class returned is a valid NMF fit object") )
	checkTrue( is(fit(obj), expect.class), paste(title, ": default class returned is ", expect.class) )
	checkTrue( !is.empty.nmf(fit(obj)) , paste(title, ": Seeded object is not empty"))
	checkEquals( nbasis(obj), r , paste(title, ": Seeded object has correct rank"))
	checkEquals( nrow(obj), nrow(V) , paste(title, ": Seeded object has correct number of rows"), checkNames=FALSE)
	checkEquals( ncol(obj), ncol(V) , paste(title, ": Seeded object has correct number of columns"), checkNames=FALSE)
	
}

check.res <- function(title, obj, V, r, expect.class, algo=NULL, seed=NULL
					, rng=NULL, rngref=NULL){
	
	# check the same thing as in the seed
	check.seed(title, obj, V, r, expect.class)
	
	# check if some slots are correctly set
	if( is.null(algo) ) algo <- nmf.getOption('default.algorithm')
	checkEquals( algorithm(obj), algo, paste(title, ": Slot 'method' is correctly set"))
	
	if( is.null(seed) ) seed <- nmf.getOption('default.seed')
	checkEquals( seeding(obj), seed, paste(title, ": Slot 'seed' is correctly set"))
	
	# rng
	if( !is.null(rng) ){
				
		# check RNG
		checkTrue( !rng.equal(rng), paste(title, ": The RNG after the fit is different from the one used to sed the computation"))
		if( nrun(obj) == 1 ){
			checkTrue(rng.equal(obj, rng), paste(title, ": The fit's RNG seed is correctly set"))
		}else{			
			if( is.list(rng) ){				
				if( is(obj, 'NMFfitXn') )
					checkTrue( all(mapply(rng.equal, obj, rng))
							, paste(title, ": The RNGs used in the multi-run computation are from the correct sequence"))
			}			
		}
		
		# check RNG_1
		rng1 <- if( is.list(rng) ) rng[[1]] else rng
		if( !is.null(rngref) )
			checkTrue( !rng.equal(rng1, rngref), paste(title, ": The initial current RNG is different from the first RNG used in computation"))		
		checkTrue(rng1.equal(obj, rng1), paste(title, ": The first fit's RNG seed is correctly set"))
		checkIdenticalRNG( getRNG1(obj), rng1, paste(title, ": The first RNG used in the computation is given by getRNG1"))
	}
	# ref rng 
	if( !is.null(rngref) ){
		checkTrue( rng.equal(rngref), paste(title, ": The current RNG was not affected by the computation"))
		if( is.null(rng) )
			checkTrue( rng.equal(obj, rngref), paste(title, ": The fit's RNG is the same as the reference RNG"))
	}
}


#' Unit test for the interface function 'seed'
test.seed <- function(){
	
	# create a random target matrix
	r <- 3; V <- .testData(r=r)
	n <- nrow(V); m <- ncol(V);
		
	# test default call
	obj <- seed(V, r)
	check.seed('Call with rank', obj, V, r, 'NMFstd')
	
	# test call with numeric value
	obj <- seed(V, r, 123456)
	check.seed('Call with rank and numeric seed', obj, V, r, 'NMFstd')
	
	# test call with name and extra parameters
	obj <- seed(V, r, 'nndsvd', densify='average')
	check.seed('Call with name and extra parameters', obj, V, r, 'NMFstd')
	# test error when unused argument is used
	checkException(seed(V, r, 'random', toto=1), "Throw an error when: unused parameter is passed to seeding method")
	
	# test providing the class to instantiate
	class.in <- 'NMFOffset'
	obj <- seed(V, list(class.in, r))
	check.seed('Call with class', obj, V, r, class.in)
	
	# test with an empty initalization object of another class: the class should not change
	class.in <- 'NMFOffset'
	obj <- nmfModel(r, model=class.in)
	obj <- seed(V, obj)
	check.seed('Call with object', obj, V, r, class.in)	
	
	# test calls of methods
	checkException({obj <- seed(V, r, 'seed.toto')}, 'Error when calling with an undefined seeding method')
	checkTrue(inherits(seed(V, r, 'random'), 'NMFfit'), 'No error when calling with a defined seeding method')
	
}

#' Unit test for the interface function 'nmf': minimum default call
test.nmf.default <- function(){
	
	# set random seed
	set.seed(123456)	
	# create a random target matrix
	r <- 3; V <- .testData(r=r)
	
	# run nmf with no argument 
	check.res('Call with rank (only)'
			, nmf(V, r)
			, V, r, 'NMFstd')
}

#' Unit test for the interface function 'nmf': argument 'method'
test.nmf.method <- function(){
	
	# set random seed
	set.seed(123456)	
	# create a random target matrix
	r <- 3; V <- .testData(r=r)
	
	# check errors
	checkException( nmf(V, r, 'zzz'), "Throw an error when: inexistent algorithm name")
	checkException( nmf(V, r, toto=3), "Throw an error when: unused argument is passed to algorithm")
		
	## ARGUMENT: method
	# run nmf with only an algorithm name 
	check.res('Call with algorithm name'
			, nmf(V, r, 'nsNMF')
			, V, r, 'NMFns', 'nsNMF')
	
	# run nmf with only an algorithm name (partial match) 
	check.res('Call with algorithm name and partial match'
			, nmf(V, r, 'ns')
			, V, r, 'NMFns', 'nsNMF')
	
	old.rseed <- getRNG()	
	res <- nmf(V, r, list('ns', 'br', 'lee'))
	checkIdentical(names(res), c('nsNMF', 'brunet', 'lee'), "Argument list(): names are set correctly to the complete method names")
	checkTrue( all(sapply(res, function(x) identical(getRNG(x), getRNG(res[[1]])))),
			"Initial RNG settings are the same for each method")
	checkTrue( !identical(old.rseed, getRNG()), "RNG setting is different after the run" )
	new.rseed <- getRNG()
	setRNG(old.rseed)
	nmf(V, r, 'lee')
	checkIdentical( new.rseed, getRNG(), "RNG setting after the run is the same as if one has run only the last method" )
	
    # list of methods
	res <- nmf(V, r, list('ns', 'br', 'lee'), nrun=3)
	checkIdentical(names(res), c('nsNMF', 'brunet', 'lee'), "Argument list() + multiple run: names are set correctly to the complete method names")
	checkTrue( all(sapply(res, function(x) identical(getRNG1(x), getRNG1(res[[1]])))),
			"Argument list() + multiple runs: Initial RNG settings are the same for each method")
    
    ml <- list('ns', 'brunet', 'lee')
    checkException(nmf(V, r, ml, .parameters = 2:3), "Error if argument .parameters not a list")
    checkException(nmf(V, r, ml, .parameters = list(list(copy = TRUE))), "Error if argument .parameters has no names")
    checkException(nmf(V, r, ml, .parameters = list(br = list(), list(copy = TRUE))), "Error if argument .parameters has missing names")
    checkException(nmf(V, r, ml, .parameters = list(br = list(), brun = list(copy = TRUE))), "Error if argument .parameters has multiple matching names")
    checkWarning(nmf(V, r, ml, .parameters = list(br = list(aaa = 1))), TRUE, "Error if unused argument in selected method-specific parameters")
    checkWarning(nmf(V, r, ml, .parameters = list(br = list(), toto = list())), TRUE, "Warning if unused elements in .parameters")
    checkTrue(all(isNMFfit(res <- nmf(V, r, ml, seed = 123, .parameters = list(br = list(maxIter = 10), ns = list(maxIter = 2)))))
                , "List of methods working if called with correct .parameters")
    checkIdentical( niter(res$nsNMF), 2L, )
    checkIdentical( niter(res$brunet), 10L)
    res_lee <- nmf(V, r, 'lee', seed = 123)
    checkTrue( niter(res$lee) > 10L, "Method without method-specific parameter specification correctly runs without them: niter > 10L" )
    checkIdentical( niter(res$lee), niter(res$lee), "Method without method-specific parameter specification correctly runs without them: niter equal" )
    checkTrue( nmf.equal(res$lee, res_lee), "Method without method-specific parameter specification correctly runs without them: identical result" )
    
}

#' Unit test for multiple rank
test.nmf.multirank <- function(){
	
	# set random seed
	set.seed(123456)
	# create a random target matrix
	r <- 3; V <- .testData(r=r)
	
	old.rseed <- getRNG()
	ranks <- 2:4
	old.rseed <- getRNG()
	res <- nmf(V, ranks, nrun=1)
	checkTrue( is(res, 'NMF.rank'), "Result is a NMF.rank object")
	checkIdentical(names(res), c('measures', 'consensus', 'fit'), "result names are corrects")
	checkIdentical(names(res$fit), as.character(ranks), "names of fits are the ranks as character strings")
	fits <- res$fit
	checkTrue( all(sapply(fits, function(x) identical(getRNG(x), getRNG(fits[[1]])))),
			"Initial RNG settings are the same for each rank")
	checkTrue( !identical(old.rseed, getRNG()), "RNG setting is different after the run" )
	new.rseed <- getRNG()
	setRNG(old.rseed)
	nmf(V, tail(ranks, 1))
	checkIdentical( new.rseed, getRNG(), "RNG setting after the run is the same as if one has run only the last rank" )
	
	res <- nmf(V, ranks, nrun=3)
	checkTrue( is(res, 'NMF.rank'), "multiple runs: Result is a NMF.rank object")
	checkIdentical(names(res), c('measures', 'consensus', 'fit'), "multiple runs: result names are corrects")
	checkIdentical(names(res$fit), as.character(ranks), "multiple runs: names of fits are the ranks as character strings")
	fits <- res$fit
	checkTrue( all(sapply(fits, function(x) identical(getRNG1(x), getRNG1(fits[[1]])))),
			"multiple runs: Initial RNG settings are the same for each rank")
}

#' Unit test for fault tolerance of .Random.seed
test.nmf.seed.fault <- function(){
	
	# set random seed
	set.seed(123456)
	# create a random target matrix
	r <- 3; V <- .testData(r=r)
	
	# .Random.seed is not changed after an error in the run
	os <- .Random.seed
	try(res <- nmf(V, r, method=function(...){}))
	checkIdentical( os, .Random.seed, ".Random.seed is NOT changed after error in single run without seed")
	#
	os <- .Random.seed
	try(res <- nmf(V, r, nrun=3, method=function(...){}, .opt='-p'))
	checkIdentical(os, .Random.seed, ".Random.seed is NOT changed after error in multiple runs without seed (sapply)")
	#
	os <- .Random.seed
	try(res <- nmf(V, r, nrun=3, method=function(...){}, .opt='P2'))
	checkIdentical(os, .Random.seed, ".Random.seed is NOT changed after error in multiple runs without seed (foreach-MC)")
	#
	os <- .Random.seed
	try(res <- nmf(V, r, nrun=3, method=function(...){}, .opt='P1'))
	checkIdentical(os, .Random.seed, ".Random.seed is NOT changed after error in multiple runs without seed (foreach-SEQ)")
	
}

#' Unit test for the interface function 'nmf': argument 'seed'
test.nmf.seed.argument <- function(){
	
	# set random seed
	set.seed(123456)
	# create a random target matrix
	r <- 3; V <- .testData(r=r)
	
	## ARGUMENT: seed
	
	# check errors
	checkException( nmf(V, r, seed='zzz'), "Throw an error when: inexistent seeding method name")
	checkException( nmf(V, r, seed=matrix(1,2,2)), "Throw an error when: invalid object as seeding method")
	checkException( nmf(V, r, seed=list('zzz')), "Throw an error when: inexistent seeding method name (passed as list)")
	checkException( nmf(V, r, seed=list(method='zzz')), "Throw an error when: inexistent seeding method name (passed as named list)")
	checkException( nmf(V, r, seed=list(toto=1, method='random')), "Throw an error when: unused argument is passed to seeding method")
    checkException( nmf(V, r, seed=numeric()), "Throw an error when: seed argument is an empty numeric")
    
	checkException( nmf(V, r, seed=c(1,2)), "Throw an error when: seed argument is a numeric of invalid length (2)")
	checkException( nmf(V, r, seed=rep(5,5)), "Throw an error when: seed argument is an invalid numeric value for .Random.seed (7)")
	
	# run nmf with only a seeding method name 
	set.seed(123)
	rngRef <- getRNG()
	check.res('Call with only a NON random seeding method name'
			, nmf(V, r, seed='nndsvd')
			, V, r, 'NMFstd', seed='nndsvd', rngref=rngRef)
	
	# run nmf with only a seeding method name partially matched
	check.res('Call with only a NON random seeding method name partially match'
			, nmf(V, r, seed='nnd')
			, V, r, 'NMFstd', seed='nndsvd')
	
	# run nmf with 'random' seeding method
	set.seed(1234)
	rngRef <- getRNG()
	check.res('Call with "random" seeding method name'
			, res <- nmf(V, r, seed='random')
			, V, r, 'NMFstd', seed='random', rng=rngRef)
	
	# run nmf with numeric seed
	msg <- function(...) paste("Call with only a numerical seed:", ...)
	rngRef <- getRNG()
	s <- nextRNG(123456)	
	check.res(msg()
			, res <- nmf(V, r, seed=123456)
			, V, r, 'NMFstd', seed='random', rng=s, rngref=rngRef)
	
	# run nmf with 6-length numeric seed
	msg <- function(...) paste("Call with a 6-length numerical seed:", ...)
	runif(10)
	rngRef <- getRNG()
	nseed <- c(1,2,3,4,5,6)
	s <- RNGseq(1, nseed)	
	check.res(msg()
			, res <- nmf(V, r, seed=nseed)
			, V, r, 'NMFstd', seed='random', rng=s, rngref=rngRef)
	
#	# run nmf with rstream object
#	msg <- function(...) paste("Call with only a rstream object:", ...)
#	rngRef <- getRNG()
#	s <- new('rstream.mrg32k3a')
#	check.res(msg()
#			, res <- nmf(V, r, seed=s)
#			, V, r, 'NMFstd', seed='random', rng=s, rngref=rngRef)
#		
	# run multi-nmf with numeric seed
	msg <- function(...) paste("Multirun - parallel + numeric seed (keep all):", ...)
	runif(10)
	rngRef <- getRNG()
	sRNG <- RNGseq(3, seed=5698)
	check.res(msg()
			, res <- nmf(V, r, nrun=3, seed=5698, .opt='kP')
			, V, r, 'NMFstd', seed='random', rng=sRNG, rngref=rngRef)
	
	# run multi-nmf with 7-length numeric seed (keep one)
	msg <- function(...) paste("Multirun - parallel + list of seeds (keep all):", ...)
	runif(10)
	rngRef <- getRNG()
	check.res(msg()
			, res2 <- nmf(V, r, nrun=3, seed=sRNG, .opt='kP')
			, V, r, 'NMFstd', seed='random', rng=sRNG, rngref=rngRef)
	checkIdenticalRNG( res2, res, msg("The best fit's RNG is the same as when seeding with corresponding single numeric seed"))
	
	# run multi-nmf with numeric seed (keep one)
	msg <- function(...) paste("Multirun - parallel + numeric seed (keep best):", ...)
	runif(10)
	res2 <- nmf(V, r, nrun=3, seed=5698, .opt='P')
	checkIdenticalRNG( res2, res, msg("The best fit's RNG is the same as when keeping all the fits"))
	checkIdenticalRNG( getRNG1(res2), sRNG[[1]], msg("The first RNG used in the computation of the NMFfitX1 object is given by getRNG1"))
	checkTrue( rng1.equal(res2, sRNG[[1]]), msg("The first RNG used in the computation is correct"))
	
	# run multi-nmf with 7-length numeric seed (keep one)
	msg <- function(...) paste("Multirun - parallel + single 7-length numeric seed (keep best):", ...)
	runif(10)
	res2 <- nmf(V, r, nrun=3, seed=sRNG[[1]], .opt='P')
	checkIdenticalRNG( res2, res, msg("The best fit's RNG is the same as when keeping all the fits"))
	checkIdenticalRNG( getRNG1(res2), sRNG[[1]], msg("The first RNG used in the computation of the NMFfitX1 object is given by getRNG1"))
	checkTrue( rng1.equal(res2, sRNG[[1]]), msg("The first RNG used in the computation is correct"))
	
#	# run multi-nmf with rstream object
#	msg <- function(...) paste("Multirun - parallel + rstream seed:", ...)
#	rngRef <- getRNG()
#	sRNG <- new('rstream.mrg32k3a')
#	check.res(msg()
#			, res <- nmf(V, r, nrun=3, seed=sRNG, .opt='kP')
#			, V, r, 'NMFstd', seed='random')	
#	checkIdenticalRNG( res[[1]], sRNG, msg("The first RNG used in the computation is correct"))	
#	checkTrue( !rng.equal(sRNG), msg("The current RNG is different from the first one used to seed the computation"))	
#	checkTrue( rng.equal(rngRef), msg("The current RNG was not affected by the computation"))	
#	checkIdenticalRNG( getRNG1(res), sRNG, msg("The first RNG used in the computation of the NMFfitXn object is given by getRNG1"))
	# run multi-nmf with rstream seed (keep one)
	msg <- function(...) paste("Multirun - parallel + rstream seed (keep best):", ...)	
	runif(10)	
	res2 <- nmf(V, r, nrun=3, seed=sRNG[[1]], .opt='P')
	checkIdenticalRNG( res2, res, msg("The best fit's RNG is the same as when keeping all the fits"))
	checkIdenticalRNG( getRNG1(res2), sRNG[[1]], msg("The first RNG used in the computation of the NMFfitX1 object is given by getRNG1"))
		
    # Seeding with NMF object
	obj.s <- rnmf(r, V)
    rngRef <- getRNG()
	res <- nmf(V, obj.s)
	check.res('Call with rank = <NMF object>', res, V, r, 'NMFstd', 'brunet', 'NMF', rngref = rngRef)
    checkTrue( nmf.equal(res, nmf(V, obj.s)), 'Run with rank=<NMF object> is deterministic')
    res.s <- nmf(V, seed = obj.s)
    check.res('Call with seed = <NMF object>', res, V, r, 'NMFstd', 'brunet', 'NMF', rngref = rngRef)
    checkTrue( nmf.equal(res, res.s), 'Run with rank=<NMF object> returns identical result as with seed=<NMF object>')
	
	# run nmf with only a seeding method name and some extra parameters 
	check.res('Call with only a seeding method name and some extra parameters (element method first and not named)'
			, nmf(V, r, seed=list('nndsvd', densify='average'))
			, V, r, 'NMFstd', seed='nndsvd')
	
	check.res('Call with only a seeding method name and some extra parameters (element method second and named)'
			, nmf(V, r, seed=list(densify='average', method='nndsvd'))
			, V, r, 'NMFstd', seed='nndsvd')
	
	# run nmf with both algorithm and seeding method 
	check.res('Call with both algorithm and seeding method'
			, nmf(V, r, 'lee', seed='nndsvd')
			, V, r, 'NMFstd', 'lee', 'nndsvd')
}

test.nmf.seed.equivalent <- function(){
	
	# set random seed
	set.seed(123456)
	# create a random target matrix
	r <- 3; V <- .testData(r=r)
	
	# multiple run nmf with numeric seed is equivalent to set.seed before the call
	set.seed(1234)
	ss_SEQ <- nmf(V, r, nrun=3, .opt='k-p')
	#
	runif(10)
	SEQ <- nmf(V, r, nrun=3, seed=1234, .opt='k-p')	
	checkTrue( nmf.equal(ss_SEQ, SEQ, all=TRUE)
			, "Multiple run using sapply with a numeric seed is equivalent to set.seed + run (sapply)")
	
	# multiple run nmf with numeric seed is equivalent to set.seed before the call
	set.seed(1234)
	ss_PAR <- nmf(V, r, nrun=3, .opt='kP2')
	runif(10)
	PAR <- nmf(V, r, nrun=3, seed=1234, .opt='kP2')
	checkTrue( nmf.equal(ss_SEQ, PAR, all=TRUE)
			, "Multiple run using foreach with a numeric seed is equivalent to set.seed + run (sapply)")
	checkTrue( nmf.equal(ss_PAR, PAR, all=TRUE)
			, "Multiple run using foreach with a numeric seed is equivalent to set.seed + run (foreach)")
	PAR_SEQ <- nmf(V, r, nrun=3, seed=1234, .opt='kP', .pbackend='seq')
	checkTrue( nmf.equal(PAR_SEQ, PAR, all=TRUE)
			, "Multiple run using foreach with a numeric seed is equivalent to foreach sequential with numeric seed")
	set.seed(1234)
	ss_PAR_SEQ <- nmf(V, r, nrun=3, .opt='kP', .pbackend='seq')
	checkTrue( nmf.equal(ss_PAR_SEQ, PAR, all=TRUE)
			, "Multiple run using foreach with a numeric seed is equivalent to set.seed + foreach sequential")
	#and:
	set.seed(1234)
	ss_SEQ_noR <- nmf(V, r, nrun=3, .opt='k-pR')
	runif(10)
	SEQ_noR <- nmf(V, r, nrun=3, seed=1234, .opt='k-pR')	
	checkTrue( nmf.equal(ss_SEQ_noR, SEQ_noR, all=TRUE)
			, "Multiple run using sapply with a numeric seed WITHOUT option R is equivalent to set.seed + run (sapply) WITHOUT option R")
	checkTrue( !nmf.equal(SEQ_noR, SEQ, all=TRUE)
			, "Multiple run using sapply with a numeric seed WITHOUT option R is NOT equivalent to set.seed + run (sapply)")
	#

	# fits of multiple runs are not the same as the ones obtained from separate fits	
	set.seed(1234)
	ss_SEPA <- replicate(3, nmf(V,r))	
	checkTrue( !nmf.equal(ss_SEPA, ss_SEQ, all=TRUE)
			, "Fits of multiple runs with sapply and a seed are NOT the same as the ones obtained from set.seed + separate fits")
	# but:
	checkTrue( nmf.equal(ss_SEPA, ss_SEQ_noR, all=TRUE)
			, "Fits of multiple runs WITHOUT option R with sapply and a seed are the same as the ones obtained from set.seed + separate fits")
	#
	checkTrue( !nmf.equal(ss_SEPA, PAR, all=TRUE)
			, "Fits of multiple runs with foreach and a seed are NOT the same as the ones obtained from set.seed + separate fits")
	# and:
	PAR_noR <- nmf(V, r, nrun=3, seed=1234, .opt='kP2-R')
	checkTrue( nmf.equal(PAR_noR, PAR, all=TRUE), "Option -R has no effect on true parallel computations")
	
}

test.nmf.seed.repro <- function(){
	
	# set random seed
	set.seed(123456)
	# create a random target matrix
	r <- 3; V <- .testData(r=r)
	
	# check reproducibility
	res <- replicate(2, nmf(V, r, seed=123))
	checkTrue( nmf.equal(res[[1]], res[[2]]), "Results are reproducible: single run" )
	res <- replicate(2, nmf(V, r))
	checkTrue( !nmf.equal(res[[1]], res[[2]]), "Results are NOT the same if not seeded: single run" )
	#
	res <- replicate(2, nmf(V, r, nrun=3, seed=123, .opt='kP'), simplify=FALSE)
	checkTrue( nmf.equal(res[[1]], res[[2]], all=TRUE), "Results are reproducible: multiple run - Parallel" )
	res <- replicate(2, nmf(V, r, nrun=3, .opt='kP'), simplify=FALSE)
	checkTrue( !nmf.equal(res[[1]], res[[2]]), "Results are NOT the same if not seeded: multiple run - Parallel" )
	#
	res <- replicate(2, nmf(V, r, nrun=3, seed=123, .opt='k-p'), simplify=FALSE)
	checkTrue( nmf.equal(res[[1]], res[[2]], all=TRUE), "Results are reproducible: multiple run - sapply" )
	res <- replicate(2, nmf(V, r, nrun=3, .opt='k-p'), simplify=FALSE)
	checkTrue( !nmf.equal(res[[1]], res[[2]]), "Results are NOT the same if not seeded: multiple run - sapply" )
	#
	res <- list(nmf(V, r, nrun=3, seed=123, .opt='kP', .pbackend='seq')
			, nmf(V, r, nrun=3, seed=123, .opt='kP', .pbackend='mc'))
	checkTrue( nmf.equal(res[[1]], res[[2]], all=TRUE), "Identical results from seeded foreach MC and SEQ" )
	res <- list(nmf(V, r, nrun=3, .opt='kP', .pbackend='seq')
			, nmf(V, r, nrun=3, .opt='kP', .pbackend='mc'))
	checkTrue( !nmf.equal(res[[1]], res[[2]], all=TRUE), "NON-identical results from non-seeded foreach MC and SEQ" )
	
}

#' Unit test for the interface function 'nmf': argument 'model'
test.nmf.model <- function(){
	
	# set random seed
	set.seed(123456)	
	# create a random target matrix
	r <- 3; V <- .testData(r=r)
	
	## ARGUMENT: model
	# run nmf with empty argument'model' 
	check.res("Call with empty argument 'model'"
			, nmf(V, r, model=list())
			, V, r, 'NMFstd')
	
	# run nmf with bad types in 'model' 
	checkException(nmf(V, r, model=NA), "Error if there argument 'model' is of a bad type: NA")
	checkException(nmf(V, r, model='toto'), "Error if there argument 'model' is of a bad type: character")
	checkException(nmf(V, r, model=12), "Error if there argument 'model' is of a bad type: numeric")
	
	# run nmf with not named element in argument 'model' 
	checkException(nmf(V, r, model=list('toto')), "Error if there is a not named element in argument 'model'")
	
	# run nmf with bad slot name in argument 'model' 
	checkException(nmf(V, r, model=list(toto=5)), "Error if there is a bad slot name in argument 'model'")
			
	# run nmf specifying arguments for initialization in argument 'model'
	res <- nmf(V, r, 'nsNMF', model=list(theta=0.6))
	check.res("Call with argument 'model' to specify extra initialization parameter"
			, res
			, V, r, 'NMFns', 'nsNMF')
	checkEquals(fit(res)@theta, 0.6
			, "Call algo:nsNMF with theta in argument 'model': argument correctly passed to model")				
}

str_dim <- NMF:::str_dim

test.nmfModel.formula <- function(){
    
    set.seed(123456)
    r <- 3
    V <- .testData(r = r)
    w <- rmatrix(nrow(V), r)
    h <- rmatrix(r, ncol(V))
    cx <- runif(ncol(V))
    bx <- runif(nrow(V))
    
    .check <- function(res, dims, msg, cterm = NULL, bterm = NULL){
        .msg <- function(...) paste0(msg, ': ', ...)
        checkTrue(isNMFfit(res), .msg('Result is an NMFfit object'))
        checkEquals(dim(res), dims, .msg('Dimensions [', str_dim(res), '] are as expected [', str_dim(dims=dims), ']'))
        
        # check fixed terms don't change
        if( !is.null(cterm) ){
            if( is.null(dim(cterm)) ) cterm <- matrix(cterm, 1L)
            else if( is.data.frame(cterm) ) t(as.matrix(cterm))
            else if( !is.matrix(cterm) ) stop("Unexpected error: invalid data type [", class(cterm), ']')
            n <- nrow(cterm)
            ft <- coef(res)[tail(1:nbasis(res), n), , drop = FALSE]
            dimnames(ft) <- NULL
            checkIdentical(cterm, ft, "Fixed coef term don't change")
        }
        if( !is.null(bterm) ){
            if( is.null(dim(bterm)) ) bterm <- matrix(bterm, ncol = 1L)
            else if( is.data.frame(bterm) ) as.matrix(bterm)
            else if( !is.matrix(cterm) ) stop("Unexpected error: invalid data type [", class(bterm), ']')
            n <- ncol(bterm)
            ft <- basis(res)[, tail(1:nbasis(res), n), drop = FALSE]
            dimnames(ft) <- NULL
            checkIdentical(bterm, ft, "Fixed basis term don't change")
        }
    }
    
    # coef terms
    .check(nmf(V ~ cx), c(dim(V), 1L), cterm = cx, 'Single coef term')
    .check(nmf(V ~ h), c(dim(V), nrow(h)), cterm = h, 'Matrix coef term')
    .check(nmf(t(V) ~ t(w)), c(dim(t(V)), ncol(w)), cterm = t(w), 'Matrix coef term (transpose)')
    .check(nmf(V ~ data.frame(t(h))), c(dim(V), nrow(h)), cterm = h, 'Data frame coef term')
    # basis terms
#    .check(nmf(V ~ bx), c(dim(V), 1L), bterm = bx, 'Single basis term')
#    .check(nmf(V ~ w), c(dim(V), ncol(w)), bterm = w, 'Matrix basis term')
#    .check(nmf(V ~ data.frame(w)), c(dim(V), ncol(w)), bterm = w, 'Data frame basis term')
    
}

#' Unit test for the interface function 'nmf': argument '...'
test.nmf.dots <- function(){
	
	# set random seed
	set.seed(123456)	
	# create a random target matrix
	r <- 3; V <- .testData(r=r)
	
	## ARGUMENT: ...
	# run nmf with unused parameter in '...'
	checkException(nmf(V, r, toto=5), "Error if there is an unused parameter in '...'")
	
	# run nmf forcing using argument in '...' for algorithm
	checkException(nmf(V, r, 'nsNMF', model=list(), theta=0.6), "Forcing argument to go to algo: error if there is an unused parameter in '...'")
	
	# run nmf specifying arguments for initialization in argument '...'
	res <- nmf(V, r, 'nsNMF', theta=0.6)
	check.res("Call with argument '...' to specify extra initialization parameter"
			, res
			, V, r, 'NMFns', 'nsNMF')
	checkEquals(fit(res)@theta, 0.6
			, "Call algo:nsNMF with theta in argument '...': argument correctly passed to model")
	

}

test.nmf.callback <- function(){

	# create a random target matrix
	r <- 3; V <- .testData(r=r)
	
	# check that the result of the callback are stored in the result object
	cb <- function(object, i){ 1 }
	res <- nmf(V, r, nrun=3, .callback=cb, .opt='P')
	checkEquals(res$.callback, rep(1, 3), 'Result of callback is in: res$.callback (PAR)')
	res <- nmf(V, r, nrun=3, .callback=cb, .opt='-p')
	checkEquals(res$.callback, rep(1, 3), 'Result of callback is in: res$.callback (SEQ)')
	
	# check that callback can make use of the result of each run
	cb <- algorithm
	res <- nmf(V, r, nrun=3, .callback=cb, .opt='P')
	checkEquals(res$.callback, rep('brunet', 3), 'Result of callback can use the result object of each run (PAR)')	
	res <- nmf(V, r, nrun=3, .callback=cb, .opt='-p')
	checkEquals(res$.callback, rep('brunet', 3), 'Result of callback can use the result object of each run (SEQ)')
	
	# check that the callback is not used with option 'keep.all' 
	cb <- function(object, i){	stop() }
	checkWarning(res <- nmf(V, r, 'br', nrun=3, .callback=cb, .opt='Pk'), "discarding argument .*\\.callback")
	checkTrue( is(res, 'NMFfitXn'), "Callback function is not used with option 'keep.all=TRUE' (PAR)")
	checkWarning(res <- nmf(V, r, 'br', nrun=3, .callback=cb, .opt='k-p'), "discarding argument .*\\.callback")
	checkTrue( is(res, 'NMFfitXn'), "Callback function is not used with option 'keep.all=TRUE' (SEQ)")
	
	# check that an error in the callback function stops the computation with an error  
	cb <- function(object, i){	stop('BIG ERROR') }	
	checkTrue(isNMFfit(res <- nmf(V, r, 'br', nrun=3, .callback=cb, .opt='P'))
		, 'Error in callback function does not stop the copmutation (PAR)')
	checkTrue(is.list(res$.callback), 'res$.callback is a list when there is an error (PAR)')
	checkEquals(sapply(res$.callback, function(x) is(x, 'error')), rep(TRUE, 3), checkNames = FALSE
		, 'Error in callback function returns errors in res$.callback (PAR)')
	
	checkTrue(isNMFfit(res <- nmf(V, r, 'br', nrun=3, .callback=cb, .opt='-p')), 'Error in callback function does not stop the copmutation (SEQ)')
	checkTrue(is.list(res$.callback), 'res$.callback is a list when there is an error (SEQ)')
	checkEquals(sapply(res$.callback, function(x) is(x, 'error')), rep(TRUE, 3), checkNames = FALSE
		, 'Error in callback function returns errors in res$.callback (SEQ)')
	
	# simplification from list if no error 
	res <- nmf(V, r, 'br', nrun=3, .callback=summary, .opt='P')
	checkTrue(is.matrix(res$.callback), 'res$.callback is a list when there is NO error (PAR)')
	res <- nmf(V, r, 'br', nrun=3, .callback=summary, .opt='P-S')
	checkTrue(is.list(res$.callback), 'res$.callback is a list when there is NO error (PAR) and simplifyCB=FALSE')
	
	res <- nmf(V, r, 'br', nrun=3, .callback=summary, .opt='-p')
	checkTrue(is.matrix(res$.callback), 'res$.callback is a list when there is NO error (SEQ)')
	res <- nmf(V, r, 'br', nrun=3, .callback=summary, .opt='-pS')
	checkTrue(is.list(res$.callback), 'res$.callback is a list when there is NO error (SEQ) and simplifyCB=FALSE')
#	
	# no simplification from list if there is at least one error
	cb <- function(object, i){ if( i ==1 ) stop('BIG ERROR ', i); summary(object) }
	res <- nmf(V, r, 'br', nrun=3, .callback=cb, .opt='P')
	checkTrue(is.list(res$.callback), 'res$.callback is a list when there is at least one error (PAR)')
	checkEquals(sapply(res$.callback, function(x) is(x, 'error')), c(TRUE, FALSE, FALSE), checkNames = FALSE
			, 'Error in callback function returns errors mixed with values in res$.callback (PAR)')
	
	res <- nmf(V, r, 'br', nrun=3, .callback=cb, .opt='-p')
	checkTrue(is.list(res$.callback), 'res$.callback is a list when there is at least one error (SEQ)')
	checkEquals(sapply(res$.callback, function(x) is(x, 'error')), c(TRUE, FALSE, FALSE), checkNames = FALSE
			, 'Error in callback function returns errors mixed with values in res$.callback (SEQ)')
	
}

test.nmf.options <- function(){
	
	x <- rmatrix(20,10)
	
	.check <- function(msg, it, ...){
		
		.msg <- function(...) paste(msg, ':', ...)
		res <- nmf(x, 2, ...)
		t <- residuals(res, track=TRUE)
		checkTrue( !is.null(names(t)), .msg("Track has names"))
		checkTrue( 0 == names(t)[1], .msg("First value in track is for iteration 0"))
		t <- t[-1]
		lags <- head(diff(as.numeric(names(t))), length(t)-2)
		checkIdentical( lags, rep(it, length(t)-2), .msg("Track interval is correct"))
		
	}
	
	.test <- function(msg, ...){
		.check(paste(msg, '- single run'), ..., nrun=1)
		.check(paste(msg, '- multiple runs'), ..., nrun=3)
	}
	
	.test('Default call -> use option', nmf.getOption('track.interval'), .options='t')
	.test('Specified in .options="t5" -> use value from "t"', 5, .options='t5')
	nmf.options(track.interval=7)
	.test('Default call after changing option -> use new option', 7, .options='t')
	
	
}

test.nmf.custom <- function(){
	
	# create a random target matrix
	r <- 3; V <- .testData(r=r)
	
	# define a dummy nmf algorithm with an argument with the same name as a slot	
	my.algo <- function(x, seed, theta=0){
		seed$extra.param <- theta
		seed
	}
	
	# check if everything works fine plain
	res <- nmf(V, r, my.algo, name='dummy')
	check.res('No argument (only name)'
			, res
			, V, r, 'NMFstd', 'dummy')
	checkEquals(res$extra.param, 0, "No argument: NO argument is provided to the algorithm")
	
	# check if everything works fine if model is an empty list
	res <- nmf(V, r, my.algo, name='dummy', model=list())
	check.res("Argument 'model' an empty list"
			, res
			, V, r, 'NMFstd', 'dummy')
	checkEquals(res$extra.param, 0, "No argument: NO argument is provided to the algorithm")
	
	# with standard model: theta is not a model parameter => Error
	checkException(nmf(V, r, my.algo, name='dummy', model=list(theta=1))
			, "Error when passing non model parameter in argument 'model', with same name as an algo parameter")
	
	# with standard model: extra argument in '...' used in algo
	res <- nmf(V, r, my.algo, name='dummy', theta=10)
	check.res('No argument (only name)'
			, res
			, V, r, 'NMFstd', 'dummy')
	checkEquals(res$extra.param, 10
			, "NMFstd + Extra argument in '...': extra argument IS provided to the algorithm")
	
	# unsued extra argument in algorithm => Error
	checkException(nmf(V, r, my.algo, name='dummy', toto=1)
			, "Error if NMFstd + unused parameter in argument '...'")
	
	# run model nsNMF plain to get default value of theta
	res <- nmf(V, r, my.algo, model='NMFns')
	default.theta <- fit(res)@theta
	custom.theta <- 1 
	checkEquals(res$extra.param, 0, "NMFns + no args: extra argument is NOT provided to the algorithm")
	
	# with model nsNMF: the parameter should be used in the model if in argument 'model'	
	res <- nmf(V, r, my.algo, name='dummy', model=list('NMFns', theta=custom.theta))
	check.res("NMFns + With extra argument in argument 'model'"
			, res
			, V, r, 'NMFns', 'dummy')
	checkEquals(res$extra.param, 0, "NMFns + Argument in 'model': extra argument is NOT provided to the algorithm")
	checkEquals(fit(res)@theta, custom.theta, "NMFns + Argument in 'model': extra argument IS used in model")
		
	# with model nsNMF: the parameter should be used in the model if in 
	# argument '...' and 'model' is not specified	
	res <- nmf(V, r, my.algo, name='dummy', model='NMFns', theta=custom.theta)
	check.res("NMFns + With extra argument in argument '...' and 'model' a model name"
			, res
			, V, r, 'NMFns', 'dummy')
	checkEquals(res$extra.param, 0, "NMFns + Argument in '...', 'model' a model name: extra argument is NOT provided to the algorithm")
	checkEquals(fit(res)@theta, custom.theta, "NMFns + Argument in '...' and 'model' a model name: extra argument used in model")
	
	# with model nsNMF: the parameter should be used in the algorithm if in 
	# argument '...' and 'model' is a list
	res <- nmf(V, r, my.algo, name='dummy', model=list('NMFns'), theta=1)
	check.res("NMFns + With extra argument in argument '...' and 'model' a list with model name"
			, res
			, V, r, 'NMFns', 'dummy')
	checkEquals(res$extra.param, 1, "NMFns + Argument in '...' and 'model' is a list with model name: extra argument IS provided to the algorithm")
	checkEquals(fit(res)@theta, default.theta, "NMFns + Argument in '...' and 'model' is list with model name: extra argument NOT used in model")
	
	# with model nsNMF: conflicts in names resolved passing different values 
	# in arguments '...' and 'model'
	res <- nmf(V, r, my.algo, name='dummy', model=list('NMFns', theta=custom.theta), theta=1)
	check.res("NMFns + With extra argument in argument 'model'"
			, res
			, V, r, 'NMFns', 'dummy')
	checkEquals(res$extra.param, 1, "NMFns + Different values in argument in '...' and 'model': correct extra argument IS provided to the algorithm")
	checkEquals(fit(res)@theta, custom.theta, "NMFns + Different values in '...' and 'model': correct extra argument IS used in model")
	
	
	# TODO: run nmf with both algorithm and seeding method
	
	# test with negative input entries
	V.neg <- V
	V.neg[1,1] <- -1
	checkException( nmf(V.neg, r, my.algo, name='dummy'), 'Throw an error if some input entries are negative and algoritham is declared NOT MIXED')
	res <- nmf(V.neg, r, my.algo, name='dummy', mixed=TRUE)
	check.res('Call with dummy MIXED algorithm on input with negative entries'
			, res
			, V.neg, r, 'NMFstd', 'dummy')
	
	# test with negative output entries
	my.algo <- function(target, start, param1, param2){
		basis(start)[1,] <- -1 
		start
	}
	res <- nmf(V, r, my.algo, name='dummy')
	check.res('Call with dummy algorithm and MIXED output'
			, res
			, V, r, 'NMFstd', 'dummy')
		
}

#' Unit test for interface nmf: testing the passage of parameters
test.nmf.parameters <- function(){
	
	# define a dummy nmf algorithm
	my.algo <- function(target, start, param1, param2){
		start
	}
	
	# create a random target matrix
	r <- 3; V <- .testData(r=r)
	
	check.res('Call with custom algo, model NMF', 
			nmf(V, r, my.algo, name='dummy', model='NMFstd')
			, V, r, 'NMFstd', 'dummy', 'random')
	
	check.res('Call with custom algo, model NMFns', 
			nmf(V, r, my.algo, name='dummy', model='NMFns')
			, V, r, 'NMFns', 'dummy', 'random')
	
	res <- nmf(V, r, my.algo, name='dummy', model=list('NMFns', theta=0.3))
	check.res('Call with custom algo, model NMFns, theta in argument model: TEST RES',
			res
			, V, r, 'NMFns', 'dummy', 'random')
	checkEquals(fit(res)@theta, 0.3
			, 'Call with custom algo, model NMFns, theta in argument model: argument correctly passed to model')
	
	res <- nmf(V, r, my.algo, name='dummy', model='NMFns', theta=0.3)
	check.res('Call with custom algo, model NMFns, theta in extra argument: TEST RES',
			res
			, V, r, 'NMFns', 'dummy', 'random')
	checkEquals(fit(res)@theta, 0.3
			, 'Call with custom algo, model NMFns, theta in extra argument: argument correctly passed to model')
	
	res <- nmf(V, r, my.algo, name='dummy', model='NMFstd', param1=0.6)	
	check.res('Call with custom algo, model NMFns, plus an extra argument: TEST RES',
			res
			, V, r, 'NMFstd', 'dummy', 'random')
	checkEquals(res@parameters, list(param1=0.6)
			, 'Call with custom algo, model NMFns, plus an extra argument: argument is passed correctly to algorithm')
	
	# redefine a dummy nmf algorithm
	my.algo2 <- function(target, start, theta){
		start
	}	
	res <- nmf(V, r, my.algo2, name='dummy', model=list('NMFns', theta=0.3), theta=0.6)
	check.res('Call with custom algo, model NMFns, theta in argument model AND extra argument: TEST RES',
			res
			, V, r, 'NMFns', 'dummy', 'random')
	checkEquals(fit(res)@theta, 0.3
		, 'Call with custom algo, model NMFns, theta in argument model AND extra argument: argument model passed to model')	
	checkEquals(res@parameters, list(theta=0.6)
			, 'Call with custom algo, model NMFns, theta in argument model AND extra argument: extra argument passed to algorithm')
	
	
	# test seeding 
#	# define a dummy seeding method
#	my.seed <- function(model, target, param.s1, param.s2){
#		rnmf(model, target)
#	}

	#res <- nmf(V, r, my.algo2, name='dummy', seed=list(), model=list('NMFns', theta=0.3), theta=0.6)

}

#' Unit test for the interface function: compare
test.compare <- function(){
	
	# create a random target matrix
	r <- 3; V <- .testData(r=r)
	m <- ncol(V)
	
	# init a list for the results
	res <- list()
	
	# compute NMF using different algorithms
	res$brunet <- nmf(V, r, 'brunet')
	res$ns <- nmf(V, r, 'ns')
	res$offset <- nmf(V, r, 'off')
	res$snmfr <- nmf(V, r, 'snmf/r')
	res$snmfl <- nmf(V, r, 'snmf/l')
	
	classes <- as.factor(sample(seq(r), m, replace=TRUE))
	
	# compare the results with a list argument
	checkTrue( is.data.frame(compare(res, class=classes)), "Result of method 'compare' (list) is a data.frame" )
	
	# compare the results with sorting
	checkTrue( is.data.frame(compare(res, class=classes, sort.by='method'))
			, "Result of method 'compare' (list) is a data.frame" )

	# try with multiple runs
	res$brunet <- nmf(V, r, 'brunet', nrun=10)
	checkTrue( is.data.frame(compare(res, class=classes)), "Result of method 'compare' (list with multiple runs) is a data.frame" )
	
	# try with multiple methods
	res <- nmf(V, r, list('brunet', 'lee', 'offset'))
	checkTrue( is.data.frame(compare(res, class=classes)), "Result of method 'compare' (list with multiple methods) is a data.frame" )
	
	# try with multiple runs multiple methods
	res <- nmf(V, r, list('brunet', 'lee', 'ns'), nrun=3)
	checkTrue( is.data.frame(compare(res, class=classes))
			, "Result of method 'compare' (list with multiple runs + multiple methods) is a data.frame" )
			
}


test.summary <- function(){
	
	# create a random target matrix
	r <- 3; V <- .testData(r=r)
	
	# test on a single run
	checkTrue( is.numeric(summary(nmf(V, r))), 'Single run: result is numeric')
	
	# test on a multiple run (no keep)
	checkTrue( is.numeric(summary(nmf(V, r, nrun=3))), 'Multiple run (no keep): result is numeric')
	
	# test on a multiple run with keep
	checkTrue( is.numeric(summary(nmf(V, r, nrun=3, .options='k'))), 'Multiple run with keep: result is numeric')
		
}

test.parallel <- function(){
		
	# create a random target matrix
	r <- 3; V <- .testData(r=r)
	# 
	# run seeded standard sequential run
	ref <- nmf(V,r, nrun=3, .opt='k-p', seed=123456)
		
	# identical results with .pbackend='seq' and standard sequential
	res <- nmf(V,r, nrun=3, .opt='kP', seed=123456, .pbackend='seq')
	checkTrue( nmf.equal(res, ref), "Identical results with seeded parallel .pbackend='seq' and standard sequential '-p'" )
	
	# identical results with .pbackend='par' and standard sequential
	res <- nmf(V,r, nrun=3, .opt='kP', seed=123456, .pbackend='par')
	checkTrue( nmf.equal(res, ref), "Identical results with seeded parallel .pbackend='par' and standard sequential '-p'" )
	
	# identical results with .pbackend='mc' and standard sequential
	res <- nmf(V,r, nrun=3, .opt='kP', seed=123456, .pbackend='mc')
	checkTrue( nmf.equal(res, ref), "Identical results with seeded parallel .pbackend='mc' and standard sequential '-p'" )
	
	# identical results with .pbackend='seq' and NA
	res <- nmf(V,r, nrun=3, .pbackend=NA, seed=123456)
	checkTrue( nmf.equal(res, ref), "Identical results with seeded .pbackend=NA and standard sequential" )
	
	# identical results with .pbackend=cluster
	cl <- makeCluster(2)
	on.exit( stopCluster(cl), add=TRUE)
	res <- nmf(V,r, nrun=3, .pbackend=cl, seed=123456)
	checkTrue( nmf.equal(res, ref), "Identical results with seeded .pbackend=cl and standard sequential" )
	
	# identical results with .pbackend=NULL and registered backend
	registerDoParallel(cl)
	on.exit( registerDoSEQ(), add=TRUE)
	res <- nmf(V,r, nrun=3, .pbackend=NULL, seed=123456)
	checkTrue( nmf.equal(res, ref), "Identical results with seeded .pbackend=NULL + registered backend and standard sequential" )
	
}

#' Unit test for the stopping criterium
test.nmf.stop <- function(){
	
	# create a random target matrix
	r <- 3; V <- .testData(r=r)
	#
	
    if( is.null(maxIter <- nmf.getOption('maxIter')) ) maxIter <- 2000L
    maxIter <- maxIter + 100L
	checkIdentical( niter(nmf(V, r, .stop=37L)), 37L, "Integer stopping criterium: fixed number of iterations")
	checkIdentical( niter(nmf(V, r, .stop=maxIter)), maxIter
		, "Integer stopping criterium greater than default maxIter: fixed number of iterations is honoured")
	checkIdentical( niter(nmf(V, r, .stop=200L, maxIter=67L)), 67L
		, "Integer stopping criterium greater than provided maxIter: maxIter is honoured")
	checkTrue( niter(nmf(V, r, .stop=37)) != 37, "Numeric stopping criterium: stationarity threshold")
	checkTrue( niter(nmf(V, r, .stop='nmf.stop.stationary')) != 37, "stopping criterium 'stationary' as full name is passed correctly")
	checkTrue( niter(nmf(V, r, .stop='stationary')) != 37, "stopping criterium 'stationary' as short name is passed correctly")
	
}