# Unit testing for multiple NMF runs
# 
# Author: Renaud Gaujoux
###############################################################################

# make the internal functions visible
if( isNamespaceLoaded('NMF') ){
	join <- NMF:::NMFfitX
}

.testData <- function(n=20, r=3, m=10, ...){
	syntheticNMF(n, r, m, ...)
}

check.result <- function(obj, V, r, size, msg=NULL){
	checkTrue( is(obj, 'NMFfitX'), paste(msg, ' -> result is an object of class NMFfitX', sep=''))
	checkEquals( nrun(obj), size, paste(msg, ' -> number of run is correctly set', sep=''))
}

test.join.singleRuns <- function(){
	
	# set random seed
	set.seed(123456)
	
	# create a random target matrix
	r <- 3; V <- .testData(r=r)
	
	# init a list for the results
	resL <- list()
	resL$brunet <- nmf(V, r)
	resL$brunet2 <- nmf(V, r)
	resL$brunet3 <- nmf(V, r)
		
	# simple join 
	res <- join(resL)
	check.result(res, V, r, length(resL), 'Simple join')
	
	res <- join(resL, runtime.all={tt<-runif(5)})
	check.result(res, V, r, length(resL), 'Simple join + set time')
	checkTrue(all.equal(tt, as.numeric(runtime.all(res))), 'Simple join + set time: time is correctly set')
	
	# merging join 
	res <- join(resL, .merge=TRUE)
	check.result(res, V, r, length(resL), 'Merging join')
	
}

test.join.multipleRuns <- function(){
	
	# set random seed
	set.seed(123456)
	
	# create a random target matrix
	r <- 3; V <- .testData(r=r)
		
	# init a list for the results
	resL <- list()
	nruns <- c(2,3,4)
	resL$brunet <- nmf(V, r, nrun=nruns[1])
	resL$brunet2 <- nmf(V, r, nrun=nruns[2])
	resL$brunet3 <- nmf(V, r, nrun=nruns[3])
	res <- join(resL)
	check.result(res, V, r, sum(nruns), 'Join multiple runs')

}


test.join.multipleAndSingleRunsMethods <- function(){
	
	# set random seed
	set.seed(123456)
	
	# create a random target matrix
	r <- 3; V <- .testData(r=r)
		
	# init a list for the results
	resL <- list()
	nruns <- c(1,3,4,1)
	resL$a <- nmf(V, r)
	resL$b <- nmf(V, r, nrun=nruns[2])
	resL$c <- nmf(V, r, nrun=nruns[3])
	resL$d <- nmf(V, r)
	res <- join(resL)
	check.result(res, V, r, sum(nruns), 'Join multiple runs + single runs')
	
}

test.multipleruns <- function(){
	
	# set random seed
	set.seed(123456)
	
	# create a random target matrix
	r <- 3; V <- .testData(r=r)
	
	# multiple runs
	res <- nmf(V, r, nrun=5)
	check.result(res, V, r, 5, 'Multiple runs')
	
}

test.interface <- function(){
	
	# set random seed
	set.seed(123456)
	
	# create a random target matrix
	r <- 3; V <- .testData(r=r)
	rownames(V) <- seq(nrow(V))
	colnames(V) <- seq(ncol(V))
	
	# perform multiple runs
	res <- nmf(V, r, nrun=5)
	
	# clusters
	checkTrue( is.factor(predict(res)), 'Clusters are computed without error')
	checkEquals( length(predict(res)), ncol(V), 'Clusters for samples have correct length')
	checkEquals( length(predict(res, 'features')), nrow(V), 'Clusters for features have correct length')
	checkTrue( is.list(predict(res, prob=TRUE)), 'Clusters are computed without error')
	
	# featureNames
	checkTrue( is.character(featureNames(res)), 'Feature names are accessible without error')
	
	# sampleNames
	checkTrue( is.character(sampleNames(res)), 'Sample names are accessible without error')
	
}

#' Unit test for the function 'fit'
test.fit <- function(){
	
	# set random seed
	set.seed(123456)
	
	# create a random target matrix
	n <- 20; r <- 3; m <- 30
	V <- rmatrix(n,m)
		
	# check result
	check.result <- function(res){
		cl <- class(res)
		
		checkTrue( is( fit(res), 'NMF' ),  paste(cl, "- fit: returns an object of class 'NMF'"))
		checkTrue( !isNMFfit( fit(res) ),  paste(cl, '- fit: does not return the complete fit'))
		checkTrue( is( minfit(res), 'NMFfit' ),  paste(cl, "- minfit: returns a 'NMFfit' object"))
	}
	
	# perform multiple runs not keeping
	check.result( nmf(V, r) )
	# perform multiple runs not keeping
	check.result( nmf(V, r, nrun=3, .opt='-k') )
	# perform multiple runs keeping
	check.result( nmf(V, r, nrun=3, .opt='k') )
}