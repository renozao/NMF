# Unit tests for the package doRepro
# 
# Author: renaud gaujoux
# Creation: 30 Jun 2011
###############################################################################

test.xapply <- function(){
	
	DEACTIVATED("Development is not finished yet")
	
	.seed <- 1:6
	n <- 3
	
	checkIdentical(xapply(1:n, .seed, function(i) i), 1:n, "Main argument is correctly passed")
	checkIdentical(xapply(1:n, .seed, function(i, b){b}, b='a'), rep('a',n), "Other arguments are correctly passed (1)")
	checkIdentical(xapply(1:n, .seed, function(i, b, c){c}, c='a'), rep('a',n), "Other arguments are correctly passed (2)")
	
	rngs <- sapply(RNGseq(n, .seed), RNGdigest)
	checkIdentical(xapply(1:n, .seed, function(i){ RNGdigest() }), rngs, "RNG are correctly set")
	
	# check that the stream seed is restored
	rngs <- sapply(RNGseq(n-1, .seed), RNGdigest)
	res <- xapply(1:n, .seed, function(i){ RNGdigest() })
	rngs <- cbind(rngs, RNGdigest())
	checkIdentical(res, rngs, "RNG are correctly set")
	
	# results are reproducible
	orng <- getRNG()
	res <- xapply(1:n, .seed, function(i) runif(i) )
	checkTrue( rng.equal(orng), "RNG is restored")
	checkIdentical(res, xapply(1:n, .seed, function(i) runif(i) ), "Results are reproducible")
	checkIdentical(sapply(res, length), 1:n, "Test results have correct dimension")
	
}

test.reproduce <- function(){
	
	DEACTIVATED("Development is not finished yet")
	
	.seed <- 1:6
	n <- 3
	p <- 5
	
	rngs <- sapply(RNGseq(3, .seed), RNGdigest)
	checkIdentical(reproduce(n, .seed, RNGdigest()), rngs, "RNG are correctly set")
	
	# results are reproducible
	orng <- getRNG()
	res <- reproduce(n, .seed, runif(p) )
	checkTrue( rng.equal(orng), "RNG is restored")
	checkIdentical(res, reproduce(n, .seed, runif(p) ), "Results are reproducible")
	checkIdentical(dim(res), as.integer(c(p,n)), "Test results have correct dimension")
}

