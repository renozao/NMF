#' Unit Testing script for NMF package: NMF algorithms.
#'
#' @author Renaud Gaujoux
#' @creation 22 April 2009

.TestSeed <- 123456

.testData <- function(n=20, r=3, m=10, ...){
	syntheticNMF(n, r, m, ...)
}

checkNMFPlot <- function(V, res, prefix=''){
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
		
		res <- nmf(V, r, meth, ...)
		
		# check the consistency of the result
		checkTrue( validObject(res), 'Returned object is valid')
		checkEquals(nrow(res), n, 'W has correct number of rows', checkNames=FALSE)
		checkEquals(nbasis(res), r, 'Result has correct rank', checkNames=FALSE)
		checkEquals(ncol(res), m, 'H has correct number of columns', checkNames=FALSE)
		checkEquals(algorithm(res), meth, "Algorithm's name is correctly set", checkNames=FALSE)
		
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
	
	return()
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

#' Tests Brunet algorithm
test.brunet <- function(){
	checkData('brunet')
}

#' Tests Lee algorithm
test.lee <- function(){
	checkData('lee')
}

#' Test NMF with offset
test.offset <- function(){
	checkData('offset')
}

#' Test Nonsmooth NMF
test.ns <- function(){
	checkData('nsNMF', theta=0.6)	
}

#' Test Sparse NMF (i.e. NNLS)
test.snmf <- function(){
	checkData('snmf/r')#, version='R')
	checkData('snmf/l')#, version='L')
	
	# check errors due to parameter checks
	set.seed(.TestSeed)
	V <- .testData()
	# beta
	checkException(nmf(V, 3, 'snmf/r', beta=0), 'beta: 0 is not a valid value')
	checkException(nmf(V, 3, 'snmf/r', beta=-1), 'beta: <0 is not a valid value')
	# bi_conv
	checkException(nmf(V, 3, 'snmf/r', bi_conv=-1), 'bi_conv: 1-length vector is not a valid value')
	checkException(nmf(V, 3, 'snmf/r', bi_conv=c(-1, 10, 3)), 'bi_conv: > 2-length vector is not a valid value')
	checkException(nmf(V, 3, 'snmf/r', bi_conv=c(-1, 10)), 'wminchange: <0 is not a valid value')
	checkException(nmf(V, 3, 'snmf/r', bi_conv=c(0, -1)), 'iconv: <0 is not a valid value')
	# eps_conv
	checkException(nmf(V, 3, 'snmf/r', eps_conv=0), 'eps_conv: 0 is not a valid value')
	checkException(nmf(V, 3, 'snmf/r', eps_conv=-1), 'eps_conv: <0 is not a valid value')
	
	# check passage of parameters
	# eps_conv
	res <- nmf(V, 3, 'snmf/r', eps_conv=1)
	checkEquals(res@parameters$eps_conv, 1, 'eps_conv: argument is passed to algorithm')
	# eta
	res <- nmf(V, 3, 'snmf/r', eta=1)
	checkEquals(res@parameters$eta, 1, 'eta: argument is passed to algorithm')
	# beta
	res <- nmf(V, 3, 'snmf/r', beta=0.05)
	checkEquals(res@parameters$beta, 0.05, 'beta: argument is passed to algorithm')
	# bi_conv
	res <- nmf(V, 3, 'snmf/r', bi_conv=c(1, 10))
	checkEquals(res@parameters$bi_conv, c(1, 10), 'bi_conv: argument is passed to algorithm')
}

#' Test Local NMF
test.lnmf <- function(){
	DEACTIVATED("Algorithm 'lnmf' was removed, for the moment.")
	checkData('lnmf')
}

#' Tests multiple runs of NMF, using Brunet algorithm Golub data.
atest.zzz.runs <- function(){
	
	# define the number of runs
	N <- 3
	r <- 3
	
	# load data 
	data(esGolub)
	eset <- esGolub[1:50, ]
	
	# run nmf N times	
	set.seed(.TestSeed)
	res <- nmf(eset, r, 'brunet', nrun=N)
	
	# check the consensus matrix
	checkPlot( basicHM(connectivity(res)), 'Consensus matrix')
	
}

#' Unit test for C and R versions of the algorithm: Brunet
check.cversion <- function(algo){
		
	r <- 3
	data(esGolub)
	eset <- esGolub[1:50,]
	resR <- nmf(eset, r, paste('.R#', algo, sep=''), seed=.TestSeed)
	res <- nmf(eset, r, algo, seed=.TestSeed)
	checkTrue( nmf.equal(res, resR, identical=FALSE), 'Results are the same for C and R version' )
	
}

#' Unit test for C and R versions of the algorithm: Brunet
test.cversions.brunet <- function(){
	
	check.cversion('brunet')
	
}

#' Unit test for C and R versions of the algorithm: Lee
test.cversions.lee <- function(){

	check.cversion('lee')
	
}

#' Unit test for C and R versions of the algorithm: nsNMF
test.cversions.ns <- function(){
	
	check.cversion('ns')
	
}

#' Unit test for C and R versions of the algorithm: Offset
test.cversions.offset <- function(){
	
	check.cversion('offset')
	
}

#' Unit test for C and R versions of the algorithm: Lee
test.cversions.lnmf <- function(){
	
	DEACTIVATED("Algorithm 'lnmf' was removed, for the moment.")
	check.cversion('lnmf')
	
}

#' Unit test for the port of `brunet`
test.port_brunet <- function(){
	
	# load RcppOctave if possible
	if( !require(RcppOctave) ) return()
	
	# source
	o_source(file.path(packagePath('inst/matlab', package='NMF'), 'brunet.m'))
	
	# define input data
	set.seed(1234)
	x <- rmatrix(100,20)
	x0 <- rnmf(3, x)
	
	# run MATLAB code: brunet(v,r,verbose, w, h)
	o <- .CallOctave('brunet', x, 3, FALSE, basis(x0), coef(x0))
	ofit <- nmfModel(o$w, o$h)
	checkTrue( !nmf.equal(ofit, x0), "MATLAB version returned something different than the seed model")
	o2 <- .CallOctave('brunet', x, 3, FALSE, basis(x0), coef(x0))
	ofit2 <- nmfModel(o2$w, o2$h)
	checkTrue( nmf.equal(ofit, ofit2), "MATLAB version really uses the seed model")
	
	# run R port
	tol <- 10^-14
	checkTrue(nmf.equal(ofit, res <- nmf(x, 3, '.R#brunet', seed=x0), tol=tol)
			, paste("Pure R port and MATLAB results are identical at ", tol))
	checkEquals(niter(res), o$niter, "Pure R and MATLAB use same number of iterations")
	
	# C version with copy	
	checkTrue(isTRUE(nmf.equal(ofit, res <- nmf(x, 3, 'brunet', seed=x0, copy=TRUE), tolerance=10^-14))
			, paste("C version with copy and MATLAB results are identical"))
	checkEquals(niter(res), o$niter, "C version without copy and MATLAB use same number of iterations")
	# C version without copy
	checkTrue(isTRUE(nmf.equal(ofit, res <- nmf(x, 3, 'brunet', seed=x0, copy=FALSE), tolerance=10^-14))
			, paste("C version without copy and MATLAB results are identical"))
	checkEquals(niter(res), o$niter, "C version without copy and MATLAB use same number of iterations")
	
}
