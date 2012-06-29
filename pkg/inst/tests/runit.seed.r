#' Unit Testing script for NMF package: seeding methods.
#'
#' @author Renaud Gaujoux
#' @creation 17 Jul 2009

# make the internal functions/objects visible
if( isNamespaceLoaded('NMF') ){
	seed <- NMF:::seed
}

RNGseed <- getRNG

checkIdenticalRNG <- function(a, b, ...){	
	checkTrue(rng.equal(a, b), ...)
}

check.setRNG <- function(value, target, title){
	
	runif(10)
	old <- getRNG()
	checkIdenticalRNG(setRNG(value), old, paste(title, ": correctly returns old value of RNG"))	
	checkIdentical(.Random.seed, target, paste(title, ": correctly sets the seed"))
	
	if( is.numeric(value) ){
		runif(100)
		old <- getRNG()
		checkIdenticalRNG(setRNG(value, test=FALSE), old, paste(title, ", test=FALSE: correctly returns old value of RNG"))	
		checkIdentical(.Random.seed, target, paste(title, ", test=FALSE: correctly sets the seed"))
		
		runif(1000)
		old <- getRNG()
		checkIdentical(setRNG(value, test=TRUE)@seed, target, paste(title, ", test=TRUE: correctly returns new value of RNG"))
		checkIdentical(.Random.seed, old@seed, paste(title, ", test=TRUE: correctly DO NOT set RNG"))
	}
	
}

#' Unit test for setRNG
test.setRNG <- function(){
	DEACTIVATED("Algorithm 'setRNG' was changed a lot.")		
	checkTrue(is(getRNG(), 'rstream'), "No arguments: returns an rstream oject")
	checkException(setRNG(c(1,2)), "Error if argument is of wrong length")
	checkException(setRNG(numeric()), "Error if argument is of wrong length")
	
	set.seed(123456)
	seed123456 <- .Random.seed
	rng123456 <- getRNG()
	
	# Single numeric argument
	check.setRNG(123456, seed123456, 'Single numeric')
		
	# .Random.seed-length numeric argument
	check.setRNG(rng123456, seed123456, '.Random.seed-length numeric')	
	checkException( setRNG(seed123456), 'Error called with a numeric too long')
	
}

test.getRNG <- function(){
	
	set.seed(123456)
	seed123456 <- .Random.seed
	rng123456 <- getRNG()
	
	checkIdentical(seed123456, rng123456, "No arguments: returns .Random.seed")
	
	# check test calls of getRNG
	runif(10)
	oseed <- .Random.seed
	rngtest <- getRNG(123456)
	checkIdentical(seed123456, rngtest, "Single numeric argument: returns .Random.seed as it would be after setting the seed")
	checkIdentical(.Random.seed, oseed, "Single numeric argument: does not change .Random.seed")
	
	oseed <- .Random.seed
	rngtest <- getRNG(.Random.seed)
	checkIdentical(oseed, rngtest, "Numeric vector argument: returns its argument unchanged")
	checkIdentical(.Random.seed, oseed, "Numeric vector argument: does not change .Random.seed")
	
}

test.setRNG2 <- function(){
	
	on.exit( RNGrecovery() )
	
	set.seed(123456)
	refseed <- .Random.seed
	
	runif(10)
	setRNG(123456)
	checkIdentical(refseed, .Random.seed, "Single numeric: sets current RNG with seed")
	RNGrecovery()
		
	# setting kind with a character string
	runif(100)
	RNGkind('Mar')
	refseed <- .Random.seed
	RNGrecovery()
	runif(100)
	setRNG('Mar')
	checkIdentical(refseed, .Random.seed, "Single character: change RNG kind")
	RNGrecovery()
	
	# setting kind with a character string
	runif(100)
	RNGkind('Mar', 'Ahrens')
	refseed <- .Random.seed
	RNGrecovery()
	runif(100)
	setRNG('Mar', 'Ahrens')
	checkIdentical(refseed, .Random.seed, "Two character strings: change RNG kind and normal kind")
	RNGrecovery()
	
	# setting kind
	runif(100)
	set.seed(123456, kind='Mar')
	refseed <- .Random.seed
	RNGrecovery()
	runif(100)
	setRNG(123456, kind='Mar')
	checkIdentical(refseed, .Random.seed, "Single numeric + kind: change RNG kind + set seed")
	RNGrecovery()
	
	# setting Nkind
	runif(100)
	set.seed(123456, normal.kind='Ahrens')
	refseed <- .Random.seed
	RNGrecovery()
	runif(100)
	setRNG(123456, normal.kind='Ahrens')
	checkIdentical(refseed, .Random.seed, "Single numeric + normal.kind: change RNG normal kind + set seed")
	RNGrecovery()
	
	# setting kind and Nkind
	runif(100)
	set.seed(123456, kind='Mar', normal.kind='Ahrens')
	refseed <- .Random.seed
	RNGrecovery()
	runif(100)
	setRNG(123456, kind='Mar', normal.kind='Ahrens')
	checkIdentical(refseed, .Random.seed, "Single numeric + kind + normal.kind: change RNG all kinds + set seed")
	RNGrecovery()
	
	# with seed length > 1
	runif(100)
	refseed <- as.integer(c(201, 0, 0))
	runif(100)
	setRNG(refseed)
	checkIdentical(refseed, .Random.seed, "numeric vector: directly set seed")
	checkIdentical(RNGkind(), c("Marsaglia-Multicarry", "Box-Muller"), "after numeric vector: RNGkind is correct")
	RNGrecovery()
	refseed <- .Random.seed
	checkException( setRNG(c(906, 1, 1)), "numeric vector: throws an error if invalid value for .Random.seed")
	checkIdentical( .Random.seed, refseed, ".Random.seed is not changed in case of an error in setRNG")
	
}

#' Unit test for seeding method: none
test.none <- function(){	
	
	# create a random target matrix
	n <- 50; r <- 3; m <- 20
	V <- syntheticNMF(n, r, m, noise=TRUE) 
	
	# seed with the matrix
	obj <- seed(V, r, 'none')
	checkTrue( is(obj, 'NMFfit') , "Seeded object is an instance of class 'NMFfit'")
	checkTrue( is(fit(obj), 'NMFstd') , "Seeded model is an instance of class 'NMFstd'")
	checkTrue( is.empty.nmf(obj), 'Should not initialize the NMF object')
	checkEquals( nbasis(obj), r, 'Seeded object have the correct rank')
	
	# seed with empty object
	obj.init <- nmfModel(r)
	obj <- seed(V, obj.init, 'none')
	checkTrue( identical(fit(obj), obj.init), 'Empty object: seeded object is identical to the initial one')
	
	# seed with dummy object
	obj.init <- nmfModel(r, model='NMFstd', W=matrix(seq(n*r), n, r), H=matrix(seq(r*m), r, m))
	obj <- seed(V, obj.init, 'none')
	checkTrue( identical(fit(obj), obj.init), 'Dummy object: seeded object is identical to the initial one')
	
}


#' Utility function for \code{test.seed}: performs a set of test on a seeded object
check.seed <- function(title, obj, V, r, seeding.meth, expect.class, exact.class=TRUE){
	
	checkTrue( inherits(obj, 'NMFfit'), paste(title, ": result class inherits from 'NMFfit'") )
	checkTrue( inherits(fit(obj), 'NMF'), paste(title, ": model class inherits from 'NMF'") )
	
	if( exact.class )
		checkTrue( is(fit(obj), expect.class), paste(title, ": default class returned is '", expect.class, "'") )
	else
		checkTrue( inherits(fit(obj), expect.class), paste(title, ": default class returned inherits from '", expect.class, "'") )
	
	checkTrue( !is.empty.nmf(obj) , paste(title, ": Seeded object is not empty"))
	checkEquals( nbasis(obj), r , paste(title, ": Seeded object has correct rank"))
	checkEquals( nrow(obj), nrow(V) , paste(title, ": Seeded object has correct number of rows"), checkNames=FALSE)
	checkEquals( ncol(obj), ncol(V) , paste(title, ": Seeded object has correct number of columns"), checkNames=FALSE)
	checkEquals( seeding(obj), seeding.meth, "Seeding method's name is correctly set")
	
}

#' Unit test for compatibility of algorithms and seeding methods
test.zzz.all <- function(){
	
	set.seed(123)
	# create a random target matrix
	n <- 50; r <- 3; m <- 20
	V <- syntheticNMF(n, r, m, noise=TRUE)
	
	# list the available algorithms
	algorithms <- nmfAlgorithm()
	algorithms <- algorithms[!algorithms %in% c('ls-nmf', 'pe-nmf')]
	# list the available seeding methods
	seed.methods <- nmfSeed()
	seed.methods <- seed.methods[which(seed.methods != 'none')]
	
	test_algo <- function(name.algo, ...){
		sapply(seed.methods,
				function(name.seed, ...){
					message("\n###########\n# ", name.algo, " + ", name.seed, "\n#############")
					err <- try(obj <- nmf(V, r, name.algo, seed=name.seed, ...))
					checkTrue( !is(err, 'try-error'), paste('Run OK - Algo:', name.algo, '+ Seed:', name.seed, if( is(err, 'try-error') ) paste('[Error: ', err, ']') else NULL) )
					check.seed(paste('Algo:', name.algo, '+ Seed:', name.seed), obj, V, r, name.seed, 'NMF', exact.class=FALSE)							
				}
		, ...)
	}
	
	sapply(algorithms, test_algo)
	test_algo('pe-nmf', alpha=1, beta=0.1)
	test_algo('ls-nmf', weight=rmatrix(V))
}


#' Utility check function: checks the range of value in a NMF object
check.range <- function(title, nmf.fit, max){
	obj <- fit(nmf.fit)
	checkTrue( all(basis(obj) <= max & basis(obj) >= 0), paste(title, ': All entries of W are between 0 and', max) )
	checkTrue( all(coef(obj) <= max & coef(obj) >= 0), paste(title, ': All entries of H are between 0 and', max) )
}

#' Unit test for seeding method: random
test.random <- function(){
	
	.seedTest <- 123456
	# create dummy matrix	
	V.na <- matrix(NA, 50, 20)
	r <- 3
		
	# check the range of the generated values
	obj <- seed(V.na, r, 'random')
	check.range('NA matrix', obj, 1)
	
	max <- 0.05
	obj <- seed(matrix(max, 50, 20), r, 'random')
	check.range(paste('Matrix', max), obj, max)
		
	# seed with the matrix
	set.seed(.seedTest)
	obj <- seed(V.na, r, 'random')
	check.seed('With matrix', obj, V.na, r, 'random', 'NMFstd')
	
	# test randomness
	obj2 <- seed(V.na, r, 'random')
	checkTrue( !identical(obj2, obj), 'Seeded objects are different if seed has not been fixed before seeding')
	# test reproducibility
	set.seed(.seedTest)
	obj.bis <- seed(V.na, r, 'random')
	checkTrue( nmf.equal(obj.bis, obj), 'Seeded NMF models are identical if seed has been reset to the same value before seeding')
	
	# seed with object
	set.seed(.seedTest)
	nmfOff <- nmfModel(r, model='NMFOffset') 
	obj <- seed(V.na, nmfOff, 'random')
	max <- 1
	check.seed('With NMFOffset object', obj, V.na, r, 'random', 'NMFOffset')
	check.range(paste('Object NMFOffset', max), obj, max)
	checkTrue( all(offset(fit(obj)) <= max & offset(fit(obj)) >= 0), paste('Object NMFOffset: All entries of Offset are between 0 and', max) )
	
	# seed with numeric value
	res <- seed(V.na, r, .seedTest)
	# manually reset name for seeding method
	set.seed(.seedTest)
	checkTrue( nmf.equal(seed(V.na, r, 'random'), res), "Seeded NMF models are identical when setting random generator seed and call method 'seed' with the same numerical seed (except for name of seeding method)")
	
}

check.seed.change <- function(msg, expr, base){
				
	bs <- .Random.seed
	e <- parent.frame()
	eval(expr, env=e)
	
	if( base ) checkTrue( any(bs != .Random.seed), paste(msg, ": .Random.seed IS changed")) 
	else checkIdentical(bs, .Random.seed, paste(msg, ": .Random.seed is NOT changed"))
				
}

#' Unit tests for checking the impact of seeding nmf computation on the .Random.seed 
test.seed.effect <- function(){
	# set random seed
	set.seed(123456)
	# create a random target matrix
	n <- 50; r <- 3; m <- 20
	V <- syntheticNMF(n, r, m, noise=TRUE)
	
	
	
	# Single runs
	check.seed.change("After single run without seed", nmf(V, r), TRUE)	
	check.seed.change("After single run without seed (NO-REPRO has no effect)", nmf(V, r), TRUE)
	check.seed.change("After single run with seed", nmf(V, r, seed=123), FALSE)
	#	
	# Multiple runs: NO seed
	check.seed.change("After multiple runs without seed (sapply)", nmf(V, r, nrun=3, .opt='-p'), TRUE)
	check.seed.change("NO-REPRO: After multiple runs without seed (sapply)", nmf(V, r, nrun=3, .opt='-pR'), TRUE)
	check.seed.change("After multiple runs without seed (foreach-MC)", nmf(V, r, nrun=3, .opt='P', .pbackend='par'), TRUE)
	check.seed.change("NO-REPRO: After multiple runs without seed (foreach-MC)", nmf(V, r, nrun=3, .opt='P-R', .pbackend='par'), TRUE)
	check.seed.change("After multiple runs without seed (foreach-SEQ)", nmf(V, r, nrun=3, .opt='P', .pbackend='seq'), TRUE)
	check.seed.change("NO-REPRO: After multiple runs without seed (foreach-SEQ)", nmf(V, r, nrun=3, .opt='P-R', .pbackend='seq'), TRUE)
	
	# Multiple runs: WITH numeric seed
	check.seed.change("After multiple runs with seed (sapply)", nmf(V, r, nrun=3, seed=1234, .opt='-p'), FALSE)
	check.seed.change("NO-REPRO: After multiple runs with seed (sapply)", nmf(V, r, nrun=3, seed=1234, .opt='-pR'), FALSE)
	check.seed.change("After multiple runs with seed (foreach-MC)", nmf(V, r, nrun=3, seed=1234, .opt='P', .pback='par'), FALSE)
	check.seed.change("NO-REPRO: After multiple runs with seed (foreach-MC)", nmf(V, r, nrun=3, seed=1234, .opt='P-R', .pback='par'), FALSE)
	check.seed.change("After multiple runs with seed (foreach-SEQ)", nmf(V, r, nrun=3, seed=1234, .opt='P', .pback='seq'), FALSE)
	check.seed.change("NO-REPRO: After multiple runs with seed (foreach-SEQ)", nmf(V, r, nrun=3, seed=1234, .opt='P-R', .pback='seq'), FALSE)
	
	# Multiple runs: WITH NA seed
#	check.seed.change("After multiple runs with NA seed (sapply)", nmf(V, r, nrun=3, seed=NA, .opt='-p'), FALSE)
#	check.seed.change("NO-REPRO: After multiple runs with NA seed (sapply)", nmf(V, r, nrun=3, seed=NA, .opt='-pR'), TRUE)
#	check.seed.change("After multiple runs with NA seed (foreach-MC)", nmf(V, r, nrun=3, seed=NA, .opt='P', .pback='par'), FALSE)
#	check.seed.change("NO-REPRO: After multiple runs with NA seed (foreach-MC)", nmf(V, r, nrun=3, seed=NA, .opt='P-R', .pback='par'), FALSE)
#	check.seed.change("After multiple runs with NA seed (foreach-SEQ)", nmf(V, r, nrun=3, seed=NA, .opt='P', .pback='seq'), FALSE, TRUE)
#	check.seed.change("NO-REPRO: After multiple runs with NA seed (foreach-SEQ)", nmf(V, r, nrun=3, seed=NA, .opt='P-R', .pback='seq'), TRUE)
}

#' test the restoration of the random seed
test.restore <- function(){
	
	DEACTIVATED("The option 'restore.seed' is deprecated. Related tests are now in test.seed.effect")
	
	# create a random target matrix
	n <- 50; r <- 3; m <- 20
	V <- syntheticNMF(n, r, m, noise=TRUE)
	
	# default call no seed
	os <- .Random.seed
	nmf(V, r)
	checkTrue( !(all.equal(os, .Random.seed) == TRUE ), "call with no seed: seed is correctly NOT restored")
	
	# default call
	os <- .Random.seed
	nmf(V, r, seed=1)
	checkIdentical(os, .Random.seed, "Default behaviour is to restore the seed: seed is correctly restored")
	
	# force restore
	os <- .Random.seed
	nmf(V, r, .opt='r', seed=12)
	checkIdentical(os, .Random.seed, "force seed restoration with 'r': seed correctly restored")
	os <- .Random.seed
	nmf(V, r, .opt=list(restore.seed=TRUE), seed=123)
	checkIdentical(os, .Random.seed, "force seed restoration with 'restore.seed=TRUE': seed correctly restored")
	
	# do not restore
	os <- .Random.seed
	nmf(V, r, .opt='-r', seed=1234)
	checkTrue( !(all.equal(os, .Random.seed) == TRUE), "Disable seed restoration with '-r': seed correctly NOT restored")
	os <- .Random.seed
	nmf(V, r, .opt=list(restore.seed=FALSE), seed=12345)
	checkTrue( !(all.equal(os, .Random.seed) == TRUE), "force seed restoration with 'restore.seed=FALSE': seed correctly NOT restored")
	
}

#' Unit test for seeding method: Non-Negative Double SVD (nndsvd)
test.nndsvd <- function(){
	
	.seedTest <- 123456
	set.seed(.seedTest)
	
	# create a random target matrix
	n <- 50; r <- 3; m <- 20
	V <- syntheticNMF(n, r, m, noise=TRUE)
		
	# perform NMF with seed 'nndsvd'
	check.seed.change('seeding with "nndsvd"', obj <- seed(V, r, 'nndsvd'), FALSE)
	check.seed('With matrix', obj, V, r, 'nndsvd', 'NMF')
	# redo: should be the same (one needs identical=FALSE because svd() may return slightly different values) 
	obj.bis <- seed(V, r, 'nndsvd')
	checkTrue( nmf.equal(obj.bis, obj, identical=FALSE), 'Seeded NMF models are identical for every run')
	
}
