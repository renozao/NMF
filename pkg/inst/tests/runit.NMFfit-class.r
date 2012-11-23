#' Unit Testing script for NMF package: NMFfit objects
#'
#' @author Renaud Gaujoux
#' @creation 22 April 2009

# make the internal functions/objects visible
if( isNamespaceLoaded('NMF') ){
	NMFfitX <- NMF:::NMFfitX
}

.TestSeed <- 123456

.testData <- function(n=20, r=3, m=10, ...){
	syntheticNMF(n, r, m, ...)
}


#' Unit test for the number of iterations
test.niter <- function(){
	
	# set random seed
	set.seed(.TestSeed)
	# generate random target matrix
	r <- 3; V <- .testData(r=r)
	
	# fit an iterative model
	res <- nmf(V,r)
	checkTrue(!is.null(niter(res)), "The number of iterations is set by the default -- iterative -- algorithm")
	# fit with SNMF/R(L)
	res <- nmf(V,r, method='snmf/r')
	checkTrue(!is.null(niter(res)), "The number of iterations is set by the SNMF/R algorithms")
	res <- nmf(V,r, method='snmf/l')
	checkTrue(!is.null(niter(res)), "The number of iterations is set by the SNMF/R algorithms")
	
	# fix number of iterations
	res <- nmf(V, r, .stop=function(strat, i, target, data, ...) if(i>=10) TRUE else FALSE)
	checkEquals(niter(res), 10, "The number of iterations is correctly set in the case of a fixed number of iterations .stop function")
	
}

test.isNMFfit <- function(){
	
	# set random seed
	set.seed(.TestSeed)
	# generate random target matrix
	r <- 3; V <- .testData(r=r)
	
	# single run
	res <- nmf(V, 2)
	checkTrue(isNMFfit(res), "isNMFfit returns TRUE on the result of a single run")
	
	# multiple runs - keeping single fit
	resm <- nmf(V, 2, nrun=3)	
	checkTrue(isNMFfit(resm), "isNMFfit returns TRUE on the result of a multiple runs - keep best")
	
	# multiple runs - keeping all fits
	resM <- nmf(V, 2, nrun=3, .opt='k') 
	checkTrue(isNMFfit(resM), "isNMFfit returns TRUE on the result of a multiple runs - keep best")
	
	# with a list of results
	checkEquals(isNMFfit(list(res, resm, resM)), rep(TRUE, 3), "isNMFfit returns [TRUE TRUE TRUE] on a list of 3 NMF results")
	checkEquals(isNMFfit(list(res, resm, resM, 'not a result')), c(rep(TRUE, 3), FALSE), "isNMFfit returns [TRUE TRUE TRUE FALSE] on a list of 3 NMF results + 1 not result")
	checkEquals(isNMFfit(list(res, resm, resM), recursive=FALSE), FALSE, "isNMFfit returns FALSE on a list of 3 NMF results when 'recursive=FALSE'")
	
}

#' Unit test for function nmf.equal
test.nmf.equal <- function(){
	
	check.nmf.equal <- function(type=c('NMF', 'NMFfit', 'NMFfitX1', 'NMFfitXn')){
		
		n <- 100; r <- 3; m <- 20
		
		# create an NMF model
		set.seed(123)
		V <- rmatrix(n, m)
		resM <- nmf(V, 3)
		resM <- NMFfitX(list(resM))
		
		## utility functions
		create.type <- function(type, obj){
			a <- if( type=='NMF' )
				fit(obj)
			else if( type=='NMFfit' )
				minfit(obj)
			else if( type=='NMFfitX1' )
				NMFfitX(obj, .merge=TRUE)
			else
				obj
	
			if( type != 'NMF' ){
				#print(class(a))
				stopifnot( class(a) == type )
			}
			a
		}
		
		add.diff <- function(obj, addon){		
			basis(fit(obj[[1]])) <- basis(fit(obj[[1]])) + rmatrix(basis(fit(obj[[1]])), max=addon)
			obj
		}
		##
		
		a <- create.type(type, resM)
		sapply(c('NMF', 'NMFfit', 'NMFfitX1', 'NMFfitXn'), function(type2){
		
				b <- create.type(type2, resM)
				
				type.pair <- paste(type, "x", type2)
				# on same object
				checkTrue( nmf.equal(a, a), paste(type.pair, "- Default: returns TRUE on same object"))
				checkTrue( nmf.equal(a, a, identical=TRUE), paste(type.pair, "- With identical=TRUE: returns TRUE on same object"))	
				checkTrue( nmf.equal(a, a, identical=FALSE), paste(type.pair, "- With identical=FALSE: returns TRUE on same object"))
				checkTrue( nmf.equal(a, a, identical=FALSE, tol=0), paste(type.pair, "- With identical=FALSE, tol=0: returns TRUE on same object"))
				checkTrue( nmf.equal(a, a, tol=0), paste(type.pair, "- With only argument tol=0: returns TRUE on same object"))
				
				# on almost same object	
				b <- add.diff(resM, .Machine$double.eps ^ 0.6)
				b <- create.type(type2, b)			
				
				checkTrue( !nmf.equal(a, b), paste(type.pair, "- Default: returns FALSE on almost same object"))
				checkTrue( !nmf.equal(a, b, identical=TRUE), paste(type.pair, "- With identical=TRUE: returns FALSE on almost same object"))
				checkTrue( nmf.equal(a, b, identical=FALSE), paste(type.pair, "- With identical=FALSE: returns TRUE on almost same object"))
				checkTrue( nmf.equal(a, b, identical=FALSE, tol=10^-4)
						, paste(type.pair, "- With identical=FALSE, tol > difference: returns TRUE on almost same object"))
				checkTrue( nmf.equal(a, b, tol=10^-4)
						, paste(type.pair, "- With only argument tol > difference: returns TRUE on almost same object"))
				checkTrue( !isTRUE(nmf.equal(a, b, identical=FALSE, tolerance= .Machine$double.eps * 2))
						, paste(type.pair, "- With identical=FALSE, tol < difference: returns FALSE on almost same object"))	
				checkTrue( !isTRUE(nmf.equal(a, b, tolerance= .Machine$double.eps * 2))
						, paste(type.pair, "- With only argument tol < difference: returns FALSE on almost same object"))
				
				# on very different object
				b <- add.diff(resM, 10)
				b <- create.type(type2, b)
				
				checkTrue( !nmf.equal(a, b), paste(type.pair, "- Default: returns FALSE on very different object"))
				checkTrue( !nmf.equal(a, b, identical=TRUE), paste(type.pair, "- With identical=TRUE: returns FALSE on very different object"))
				checkTrue( !isTRUE(nmf.equal(a, b, identical=FALSE)), paste(type.pair, "- With identical=FALSE: returns FALSE on very different object"))
				checkTrue( nmf.equal(a, b, identical=FALSE, tol=11)
						, paste(type.pair, "- With identical=FALSE, tol > difference: returns TRUE on very different object"))
				checkTrue( nmf.equal(a, b, tol=11)
						, paste(type.pair, "- With only argument tol > difference: returns TRUE on very different object"))
				checkTrue( !isTRUE(nmf.equal(a, b, identical=FALSE, tol=0.5))
						, paste(type.pair, "- With identical=FALSE, tol < difference: returns FALSE on very different object"))
				checkTrue( !isTRUE(nmf.equal(a, b, tol=0.5))
						, paste(type.pair, "- With only argument tol < difference: returns FALSE on very different object"))
		})	
	
	}
	
	
	sapply(c('NMF', 'NMFfit', 'NMFfitX1', 'NMFfitXn'), check.nmf.equal)	

}

test.deviance <- function(){
	
	
	
}