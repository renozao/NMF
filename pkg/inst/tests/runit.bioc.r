# Unit tests for the Bioconductor layer
# 
# Author: Renaud Gaujoux
# Created: Mar 6, 2013
###############################################################################


# check extension of Bioc access methods
# => this also serves to check that they are correctly exported
test.access <- function(){
	
	x <- nmfModel(3, 20, 10)	
	.check <- function(fun, val, newval){
		msg <- function(...) paste(fun, ': ', ..., sep='')
		f <- match.fun(fun)
		checkIdentical(f(x), val, msg('on fresh object returns NULL'))
		
		if( isNumber(newval) ) newval <- paste('aaa', 1:newval)
		res <- try(eval(parse(text=paste(fun, '(x) <- newval', sep=''))))
		checkTrue(!is(res, 'try-error'), msg('setting value works'))
		checkIdentical(f(x), newval, msg('new value is correct'))
		res <- try(eval(parse(text=paste(fun, '(x) <- val', sep=''))))
		checkTrue(!is(res, 'try-error'), msg('resetting value works'))
		checkIdentical(f(x), val, msg('reset value is correct'))
	}
	
	checkIdentical(nbasis(x), nmeta(x), 'nmeta is defined')
	.check('featureNames', NULL, 20)
	.check('sampleNames', NULL, 10)
	.check('basisnames', NULL, 3)
	.check('metaprofiles', coef(x), rmatrix(coef(x)))	
	.check('metagenes', basis(x), rmatrix(basis(x)))
	
}