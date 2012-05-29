# Unit tests for the rstream extension
# 
# Author: renaud
# Creation: 13 Jul 2011
###############################################################################

# make the internal functions/objects visible
if( isNamespaceLoaded('NMF') ){
	rngtools.get.provider <- NMF:::rngtools.get.provider
	rstream.provider <- NMF:::rstream.provider
	rstream.user.provider <- NMF:::rstream.user.provider
	RNGwrapper.Libpath <- NMF:::RNGwrapper.Libpath
	RNGwrapper.Libname <- NMF:::RNGwrapper.Libname
	is.RNGlib <- NMF:::is.RNGlib
	isRNGwrapper <- NMF:::isRNGwrapper
}

test.wrap <- function(){
	DEACTIVATED("Algorithm 'RNGwrap' was removed.")	
	# reset things	
	try(detach('package:rlecuyer', unload=TRUE))
	RNGrecovery()
	getRNG()	
	n <- 10
	
	# setup a user-supplied RNG
	a <- new('rstream.mrg32k3a', seed=rep(123,6), force=TRUE)
	r5 <- r(a, n)	
	rstream.reset(a)
	setRNG(a)
	checkIdentical( RNGlib(), 'rstream', "After activating rstream with RNGkind: RNG library is 'rstream'")
	checkIdentical(r5, runif(n), "Correctly set rstream.mrg32k3a object")
	setRNG(a)
	
	# load rlecuyer
	library(rlecuyer)
	checkIdentical( RNGlib(), 'rlecuyer', "After loading package rlecuyer: RNG library is 'rlecuyer'")
	
	.lec.CreateStream('a')
	.lec.CurrentStream('a')
	checkIdentical( RNGlib(), 'rlecuyer', "After activating rlecuyer with RNGkind: RNG library is still 'rlecuyer'")
	rlec <- runif(n)
	checkTrue( !all(r5 == rlec), "Lecuyer seed is not the same as rstream.mrg32k3a")
	
	.lec.CurrentStream('a')			
	checkTrue( !is.RNGlib(RNGwrapper.Libname()), "RNG wrapper library is not loaded")
	checkIdentical( rngtools.get.provider(), '', "rngtools.get.provider returns ''")
	checkIdentical( rlec, runif(10), "The values of random draws are correct")
	
	# reset things	
	try(detach('package:rlecuyer', unload=TRUE))
	RNGrecovery()
	getRNG()
	
	# redo with a RNGwrapper loaded before rlecuyer
	rstream.RNG(new('rstream.runif'))
	checkTrue( rngtools.get.provider() != 'rlecuyer', "current rstream.runif: RNG user provider is not 'rlecuyer' at start")
	
	library(randtoolbox)
	setRNG(a)	
	checkTrue( is.RNGlib(RNGwrapper.Libname()), "RNG provider not directly activable: RNG wrapper library is loaded")
	checkTrue( isRNGwrapper(), "RNG provider not directly activable: RNG wrapper library is activated")
	checkTrue( rngtools.get.provider() == 'rstream', "RNG provider not directly activable: RNG wrapper points to correct RNG library")
	
	library(rlecuyer)
	.lec.CreateStream('a')
	.lec.CurrentStream('a')
	checkTrue( rngtools.get.provider() == 'rlecuyer', "Auto-detect change: RNG wrapper points to correct RNG library")	
	
	# reset things	
	try(detach('package:rlecuyer', unload=TRUE))
	RNGrecovery()
	getRNG()
	
}

test.RNGdigest <- function(){
	
	DEACTIVATED("Algorithm 'RNGdigest' was downgraded.")
	s <- new('rstream.runif')
	checkTrue( is.character(RNGdigest(s)), "runif: RNGdigest of a stream is a character string")
	checkTrue( is.character(RNGdigest(s)), "runif: RNGdigest (call 2) of a stream is a character string")
	
	s <- new('rstream.mrg32k3a')
	checkTrue( is.character(RNGdigest(s)), "mrg32k3a: RNGdigest of a stream is a character string")
	checkTrue( is.character(RNGdigest(s)), "mrg32k3a: RNGdigest (call 2) of a stream is a character string")
	
	checkTrue( is.character(RNGdigest()), "current: RNGdigest of a stream is a character string")
	checkTrue( is.character(RNGdigest()), "current: RNGdigest (call 2) of a stream is a character string")
	
}

test.all <- function(){
	
	DEACTIVATED("Algorithm RNGtools functions were removed.")
	rseed <- .Random.seed
	rng <- getRNG()
	checkIdentical(rstream.provider(rng), 'base', "At start current provider is 'base'")
	checkIdentical(rstream.user.provider(rng), 'rstream', "At start current user provider is 'rstream'")
	checkTrue(is(rng, 'rstream.runif'), "At start RNG is a rstream.runif")
	checkIdentical(rng@seed, rseed, "At start current RNG seed is correct")
	checkIdentical(rng@xstate$state, rseed, "At start current RNG state is correct")
	
	rseed <- .Random.seed
	rng <- new('rstream.user')
	checkIdentical(rstream.provider(rng), 'rstream', "At start a rstream.user object is provided by 'rstream'")
	checkIdentical(rstream.user.provider(rng), 'rstream', "At start a rstream.user object is user-provided by 'rstream'")
	checkIdentical(.Random.seed, rseed, "Creating a rstream.user(rstream) object does not change .Random.seed")
	
	try(detach('package:rlecuyer', unload=TRUE))
	library(rlecuyer)
	rng <- getRNG()
	checkIdentical(rstream.provider(rng), 'base', "After loading rlecuyer current provider is 'base'")
	checkIdentical(rstream.user.provider(rng), 'rlecuyer', "After loading rlecuyer current user provider is 'rlecuyer'")
	checkTrue(is(rng, 'rstream.runif'), "After loading rlecuyer RNG is still a rstream.runif")
	checkIdentical(rng@seed, rseed, "After loading rlecuyer current RNG seed is still correct")
	checkIdentical(rng@xstate$state, rseed, "After loading rlecuyer current RNG state is correct")
	try(detach('package:rlecuyer', unload=TRUE))
	
} 