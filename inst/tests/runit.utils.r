#' Unit Testing script for NMF package: NMF utility functions.
#'
#' @author Renaud Gaujoux
#' @creation 10 Aug 2010


#' Unit test for rmatrix: random matrix generation
test.rmatrix <- function(){
	
	n <- 100; p <- 20
	A <- matrix(1, n, p)
	
	# square matrix if y is missing
	set.seed(123456); M <- matrix(runif(n*n), n, n)	
	set.seed(123456); checkIdentical(M, rmatrix(n), "Square matrix if 'y' is missing")
	set.seed(123456); checkIdentical(M, rmatrix(matrix(NA, nrow(M), ncol(M))), "Correct if 'x' is a matrix")
	
	# from NMF model
	model <- rnmf(3, A)
	set.seed(123456); M <- fitted(model) + matrix(runif(n*p), n, p)
	set.seed(123456); checkIdentical(M, rmatrix(model), "Correct if 'x' is an NMF model")
	set.seed(123456); M <- fitted(model) + matrix(rnorm(n*p), n, p)
	set.seed(123456); checkIdentical(M, rmatrix(model, dist=rnorm), "dist is passed correctly if 'x' is an NMF model")
	
	# default dist is uniform
	set.seed(123456); M <- matrix(runif(n*p), n, p)
	set.seed(123456); checkIdentical(M, rmatrix(n, p), "Default correctly to 'runif'")
	set.seed(123456); checkIdentical(M, rmatrix(A), "Default correctly to 'runif' (arg: matrix)")
	
	# argument byrow is correctly passed
	set.seed(123456); M <- matrix(runif(n*p), n, p, byrow=TRUE)
	set.seed(123456); checkIdentical(M, rmatrix(n, p, byrow=TRUE), "argument byrow is correctly passed")
	set.seed(123456); checkIdentical(M, rmatrix(A, byrow=TRUE), "argument byrow is correctly passed (arg: matrix)")
	
	# argument dimnames is correctly passed
	dims <-list(rep('a',n), rep('b',p))
	set.seed(123456); M <- matrix(runif(n*p), n, p, dimnames=dims)
	set.seed(123456); checkIdentical(M, rmatrix(n, p, dimnames=dims), "argument dimnames is correctly passed")
	set.seed(123456); checkIdentical(M, rmatrix(A, dimnames=dims), "argument dimnames is correctly passed (arg: matrix)")
	
	# can pass distribution function 
	set.seed(123456); M <- matrix(rnorm(n*p), n, p)
	set.seed(123456); checkIdentical(M, rmatrix(n, p, dist=rnorm), "argument dist is correctly passed")
	set.seed(123456); checkIdentical(M, rmatrix(A, dist=rnorm), "argument dist is correctly passed (arg: matrix)")
	
	# can pass distribution functions as third argument 
	set.seed(123456); M <- matrix(rnorm(n*p), n, p)
	set.seed(123456); checkIdentical(M, rmatrix(n, p, rnorm), "argument dist is the third argument")
	set.seed(123456); checkIdentical(M, rmatrix(A, rnorm), "argument dist is the second argument (arg: matrix)")
	
	# can pass extra arguments to distribution function 
	set.seed(123456); M <- matrix(rnorm(n*p, 20), n, p)
	set.seed(123456); checkIdentical(M, rmatrix(n, p, rnorm, mean=20), "extra arguments are passed to the distribution function")
	set.seed(123456); checkIdentical(M, rmatrix(A, rnorm, mean=20), "extra arguments are passed to the distribution function (arg: matrix)")
}

#test.ptr_neq_constraints <- function(){
#	
#	.do_constrain <- function(...){
#		
#	}
#	
#	.check <- function(c, msg){
#		.msg <- function(...) paste(msg, ':', ...)
#		x <- rmatrix(20,3)
#		y <- NMF:::neq.constraints.inplace(x, copy=TRUE)		
#		checkIdentical(max.col(y[1:9,]), c(rep(1,3), rep(2,3), rep(3,3)), .msg("Max are ok"))	
#		checkIdentical(y[-(1:9,], , .msg("Non constrained rows are identical"))
#	}
#	
#	#.check(list(1:3,4:6,7:9), "")
#	
#	
#	
#}


test.nmfWrapper <- function(){
	
	.msg <- NULL
	msg <- function(...) paste(.msg, ': ', ..., sep='')
	
	f <- nmfWrapper('lee')
	x <- rmatrix(20, 10)
	checkTrue( isNMFfit(res <- f(x, 3)), msg('result is an NMFfit object') )
	checkIdentical(nbasis(res), 3L, msg('result was computed using the correct rank') )
	checkIdentical(algorithm(res), 'lee', msg('result was computed using the correct algorithm') )
	
	
	.msg <- 'with default maxIter and seed value'
	f <- nmfWrapper('nsNMF', maxIter=3, seed='nndsvd')
	checkTrue( isNMFfit(res <- f(x, 2)), msg('result is an NMFfit object' ))
	checkIdentical(nbasis(res), 2L, msg('result was computed using the correct rank') )
	checkIdentical(algorithm(res), 'nsNMF', msg('result was computed using the correct algorithm') )
	checkIdentical(niter(res), 3L, msg('result was computed using the correct number of iterations') )
	checkIdentical(seeding(res), 'nndsvd', msg('result was computed using the correct seed') )
	# overwrite default in call
	.msg <- 'overwriting defaults in call'
	checkTrue( isNMFfit(res <- f(x, 4, seed='random')), msg('result is an NMFfit object' ))
	checkIdentical(nbasis(res), 4L, msg('result was computed using the correct rank') )
	checkIdentical(algorithm(res), 'nsNMF', msg('result was computed using the correct algorithm') )
	checkIdentical(niter(res), 3L, msg('result was computed using the correct number of iterations') )
	checkIdentical(seeding(res), 'random', msg('result was computed using the correct seed') )
	# pass method as well
	.msg <- 'overwriting defaults in call + try overwrite method'
	checkWarning(res <- f(x, 4, method='lee', seed='random'), 'Discarding fixed arguments.*', msg('a warning is thrown'))
	checkTrue( isNMFfit(res), msg('result is an NMFfit object' ))
	checkIdentical(algorithm(res), 'nsNMF', msg('result was still computed using the correct algorithm defined in nmfWrapper') )
	
}