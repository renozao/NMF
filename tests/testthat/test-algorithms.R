#' Unit Testing script for NMF package: NMF algorithms.
#'
#' @author Renaud Gaujoux
#' @creation 22 April 2009
# Converted from RUnit: 16 Feb 2020

library(rngtools)

# RUnit ports
DEACTIVATED <- function(msg = ""){
    skip(msg)
}

##
    
.TestSeed <- 123456

.testData <- function(n=20, r=3, m=10, ...){
	syntheticNMF(n, r, m, ...)
}

checkNMFPlot <- function(V, res, prefix=''){
	
	if( isCHECK() ) return()
	
	# check heatmaps of the target matrix, the metaprofiles and the metagenes
	checkPlot( aheatmap(V), paste(prefix, ': Target'))
	checkPlot( coefmap(res), paste(prefix, ': Metaprofiles'))
	checkPlot( basismap(res), paste(prefix, ': Metagenes'))
}

#' Tests an algorithm on all data: synthetic, Golub, Cell-type
checkData <- function(meth, ...){
	
	message("\n###############\n#", meth, "\n###############")
	# On synthetic data	
		
		#1. with noise		
		set.seed(.TestSeed)
		r <- 3; V <- .testData(r=r)
		n <- nrow(V); m <- ncol(V);
		
		res <- nmf(V, r, meth, ..., rng=1234)
		
		# check the consistency of the result
		expect_true( validObject(res), 'Returned object is valid')
		expect_equal(nrow(res), n, info = 'W has correct number of rows')
		expect_equal(nbasis(res), r, info = 'Result has correct rank')
		expect_equal(ncol(res), m, info = 'H has correct number of columns')
		expect_equal(algorithm(res), meth, info = "Algorithm's name is correctly set")
		
		wr <- nmfWrapper(meth)
		expect_true( isNMFfit(resWrapper <- wr(V, r, ..., rng=1234)), 'Wrapper call works')
		expect_true( nmf.equal(res, resWrapper), 'Wrapper call gives identical result')
		
		# check heatmaps
		checkNMFPlot(V, res, 'Synthetic [noise]')

		#2. with noise + offset
		n.offset <- 15
		o <- c(rep(1, n.offset), rep(0, n-n.offset))
		set.seed(.TestSeed)
		V <- syntheticNMF(n, r, m, offset=o, noise=TRUE)
		res <- nmf(V, r, meth, ...)		
		# check heatmaps
		checkNMFPlot(V, res, 'Synthetic [noise + offset]')
	
	# On Golub data
		data(esGolub)
		eset <- esGolub[1:50,]
		# check the heatmap of the target expression matrix		
		set.seed(.TestSeed)
		res <- nmf(eset, r, meth, ...)
		# check heatmaps
		checkNMFPlot(exprs(eset), res, 'Golub')		
	
	return(TRUE)
	# On Cell-type data
		data(CellType)
		eset <- esCellType[1:100,]
		V <- exp(exprs(eset))
		# check the heatmap of the target expression matrix		
		set.seed(.TestSeed)
		res <- nmf(V, r, meth, ...)
		# check heatmaps
		checkNMFPlot(V, res, 'Cell-type')

}

##' Tests multiple runs of NMF, using Brunet algorithm Golub data.
#atest.zzz.runs <- function(){
#	
#	# define the number of runs
#	N <- 3
#	r <- 3
#	
#	# load data 
#	data(esGolub)
#	eset <- esGolub[1:50, ]
#	
#	# run nmf N times	
#	set.seed(.TestSeed)
#	res <- nmf(eset, r, 'brunet', nrun=N)
#	
#	# check the consensus matrix
#	checkPlot( basicHM(connectivity(res)), 'Consensus matrix')
#	
#}

#' Unit test for identical results if NMF algorithms
check.algos <- function(algo1, algo2, identical=FALSE){
	
	r <- 3
	data(esGolub)
	eset <- esGolub[1:50,]
	res1 <- nmf(eset, r, algo1, seed=.TestSeed)
	res2 <- nmf(eset, r, algo2, seed=.TestSeed)
	expect_true( nmf.equal(res2, res2, identical=identical)
			, paste("Results are the same for '", algo1, "' and '", algo2, "'", sep=''))
	
}

#' Unit test for C and R versions of algorithms
check.cversion <- function(algo){
		
	check.algos(paste('.R#', algo, sep=''), algo)
	
}


test_that("test.brunet", {
    checkData("brunet")
})

test_that("test.cversions.brunet", {
    check.cversion("brunet")
})

test_that("test.cversions.lee", {
    check.cversion("lee")
})

test_that("test.cversions.lnmf", {
    DEACTIVATED("Algorithm 'lnmf' is not fully working.")
    check.cversion("lnmf")
})

test_that("test.cversions.ns", {
    check.cversion("ns")
})

test_that("test.cversions.offset", {
    check.cversion("offset")
})

test_that("test.frobenius", {
    checkData("Frobenius")
})

test_that("test.KL", {
    checkData("KL")
})

test_that("test.lee", {
    checkData("lee")
})

test_that("test.lnmf", {
    DEACTIVATED("Algorithm 'lnmf' is not fully working.")
    checkData("lnmf")
})

test_that("test.ns", {
    checkData("nsNMF", theta = 0.6)
})

test_that("test.offset", {
    checkData("offset")
})

test_that("test.port_brunet", {
    if (!require(RcppOctave)) {
        DEACTIVATED("Package RcppOctave not available.")
    }
    o_source(file.path(packagePath("m-files", package = "NMF"), 
        "brunet.m"))
    setRNG("default", "default")
    set.seed(1234)
    x <- rmatrix(100, 20)
    x0 <- rnmf(3, x)
    o <- .CallOctave("brunet", x, 3, FALSE, basis(x0), coef(x0))
    ofit <- nmfModel(o$W, o$H)
    expect_false(nmf.equal(ofit, x0), info = "MATLAB version returned something different than the seed model")
    o2 <- .CallOctave("brunet", x, 3, FALSE, basis(x0), coef(x0))
    ofit2 <- nmfModel(o2$W, o2$H)
    expect_true(nmf.equal(ofit, ofit2), info = "MATLAB version really uses the seed model")
    o_rm("brunet")
    tol <- 10^-14
    res <- nmf(x, 3, ".R#brunet", seed = x0, maxIter = 2000L)
    expect_equal(o$niter, niter(res), info = "Pure R and MATLAB use same number of iterations")
    expect_true(nmf.equal(ofit, res, tolerance = tol), info = paste("Pure R port and MATLAB results are identical at tolerance ", 
        tol))
    res <- nmf(x, 3, "brunet", seed = x0, copy = TRUE, maxIter = 2000L)
    expect_equal(o$niter, niter(res), info = "C version without copy and MATLAB use same number of iterations")
    expect_true(isTRUE(nmf.equal(ofit, res, tolerance = tol)), 
        info = paste("C version with copy and MATLAB results are identical at tolerance ", 
            tol))
    res <- nmf(x, 3, "brunet", seed = x0, copy = FALSE, maxIter = 2000L)
    expect_equal(o$niter, niter(res), info = "C version without copy and MATLAB use same number of iterations")
    expect_true(isTRUE(nmf.equal(ofit, res, tolerance = tol)), 
        info = paste("C version without copy and MATLAB results are identical at tolerance ", 
            tol))
    check.algos("brunet", ".M#brunet", identical = TRUE)
})

test_that("test.snmf", {
    checkData("snmf/r")
    checkData("snmf/l")
    set.seed(.TestSeed)
    V <- .testData()
    expect_error(nmf(V, 3, "snmf/r", beta = 0), info = "beta: 0 is not a valid value")
    expect_error(nmf(V, 3, "snmf/r", beta = -1), info = "beta: <0 is not a valid value")
    expect_error(nmf(V, 3, "snmf/r", bi_conv = -1), info = "bi_conv: 1-length vector is not a valid value")
    expect_error(nmf(V, 3, "snmf/r", bi_conv = c(-1, 10, 3)), 
        info = "bi_conv: > 2-length vector is not a valid value")
    expect_error(nmf(V, 3, "snmf/r", bi_conv = c(-1, 10)), info = "wminchange: <0 is not a valid value")
    expect_error(nmf(V, 3, "snmf/r", bi_conv = c(0, -1)), info = "iconv: <0 is not a valid value")
    expect_error(nmf(V, 3, "snmf/r", eps_conv = 0), info = "eps_conv: 0 is not a valid value")
    expect_error(nmf(V, 3, "snmf/r", eps_conv = -1), info = "eps_conv: <0 is not a valid value")
    res <- nmf(V, 3, "snmf/r", eps_conv = 1)
    expect_equal(1, res@parameters$eps_conv, info = "eps_conv: argument is passed to algorithm")
    res <- nmf(V, 3, "snmf/r", eta = 1)
    expect_equal(1, res@parameters$eta, info = "eta: argument is passed to algorithm")
    res <- nmf(V, 3, "snmf/r", beta = 0.05)
    expect_equal(0.05, res@parameters$beta, info = "beta: argument is passed to algorithm")
    res <- nmf(V, 3, "snmf/r", bi_conv = c(1, 10))
    expect_equal(c(1, 10), res@parameters$bi_conv, info = "bi_conv: argument is passed to algorithm")
})

