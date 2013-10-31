#' Unit Testing script for NMF package: NMF models and utilities.
#'
#' @author Renaud Gaujoux
#' @creation 22 April 2009

# make the internal functions/objects visible
if( isNamespaceLoaded('NMF') ){
	.predict.nmf <- NMF:::.predict.nmf
	is.same <- NMF:::is.same
}

.TestSeed <- 123456

checkOnlyNAs <- function(x, msg){
	checkTrue( all(is.na(basis(x))) , paste(msg, ': Only NAs in W'))	
	checkTrue( all(is.na(coef(x))), paste(msg, ': Only NAs in H'))	
}

basicHM <- aheatmap

checkPlot <- function(...){
	if( isCHECK() ) return()
	pkgmaker::checkPlot(...)
}

#' checks the validity of a NMF object: dimensions
check.object <- function(obj, n, m, r, title='Check NMF object', class='NMFstd'){
	
	# convert dimensions to integers
	n <- as.integer(n)
	m <- as.integer(m)
	r <- as.integer(r)
	
	# check class
	checkTrue(is(obj, class), paste(title, ': object is of class "', class, "'"));
	# check virtual interface accessors
	checkTrue(is.same(obj@W, basis(obj)), msg('method `basis` returns slot W'))
	checkTrue(is.same(obj@H, coef(obj)), msg('method `coef` returns slot H'))
	
	checkIdentical(nrow(basis(obj)), n, paste(title, ': number of rows in basis matrix is correct'))
	checkIdentical(nrow(obj), n, paste(title, ': nrow returns correct value'))
	checkIdentical(ncol(basis(obj)), r, paste(title, ': number of columns in basis matrix is correct'))
	checkIdentical(nbasis(obj), r, paste(title, ': nbasis returns correct value'))
	checkIdentical(nrow(coef(obj)), r, paste(title, ': number of rows in coef matrix correct'))
	checkIdentical(ncol(coef(obj)), m, paste(title, ': number of columns in coef matrix is correct'))
	checkIdentical(ncol(obj), m, paste(title, ': ncol returns correct value'))
}


check.NMF.class <- function(class, ...){
	
	set.seed(.TestSeed)
	
	msg <- function(...) paste('New -', class,':', ...) 
	
	# base constructor: should build an empty model
	a <- new(class, ...)
	
	# check if there is only NAs
	checkOnlyNAs(a, msg('with no parameters'))	
	check.object(a, 0, 0, 0, msg('with no parameters'), class)
	
	# define some dimensions to use as template
	n <- 50; m <- 10; r <- 3;
	w <- rmatrix(n, r); h <- rmatrix(r, m)
	
	# base constructor with one of the two matrices: exceptions	
	checkException(new(class, W=w, ...), msg('error if only basis matrix'))
	checkException(new(class, H=h, ...), msg('error if only coef matrix'))	
	
	# base constructor with two matrices		
	a <- new(class, W=w, H=h, ...)
	
	# check the dimensions of the internal matrices
	check.object(a, n, m, r, msg('with two matrices - '), class)
	checkIdentical(basis(a), w, msg('with two matrices - basis matrix is correctly set')); 
	checkIdentical(coef(a), h, msg('with two matrices - coef matrix is correctly set'));
	
	# check error with wrong dimension
	checkException(a <- new(class, W=matrix(1,n,r+1), H=matrix(1,r,m), ...), msg('Error if incompatible dimensions (W)'))
	checkException(a <- new(class, W=matrix(1,n,r), H=matrix(1,r+1,m), ...), msg('Error if incompatible dimensions (H)'))	
}

#' Tests the class 'NMF'
test.class.NMFstd <- function(){
	check.NMF.class('NMFstd')	
}

#' Tests the class 'NMFns'
test.class.NMFns <- function(){
	# with no parameter
	check.NMF.class('NMFns')
	
	# with a parameter
	t <- 0.8
	check.NMF.class('NMFns', theta=t)
	
	# test specific stuff 
	a <- new('NMFns', theta=t)
	checkIdentical(a@theta, t, 'Slot theta is correctly set')
	checkException(a <- new('NMFns', theta=-1), 'Negative value of theta throws an exception')
	checkException(a <- new('NMFns', theta=1.2), 'Value of theta > 1 throws an exception')
	
	# test equivalence with standard model when theta=0
	set.seed(.TestSeed)
	n <- 50; m <- 10; r <- 3	
	W <- rmatrix(n, r)
	H <- rmatrix(r, m)	
	a.ns <- new('NMFns', W=W, H=H, theta=0)
	a.std <- new('NMFstd', W=W, H=H)
	checkIdentical(fitted(a.ns), fitted(a.std), 'Values fitted correspond to standard model if theta=0')
	
	# check method smoothing
	a.ns <- new('NMFns', W=W, H=H, theta=0.4)
	s <- smoothing(a.ns)
	checkEquals(dim(s), c(r, r), 'Smoothing matrix: dimension are correct')
	checkTrue( all(s >= 0), 'Smoothing matrix: all entries are nonnegative')
	checkEquals( rowSums(s), rep(1, r), 'Smoothing matrix: sum of rows are ones')
	checkEquals( colSums(s), rep(1, r), 'Smoothing matrix: sum of columns are ones')
	checkEquals( fitted(a.ns), basis(a.ns) %*% s %*% coef(a.ns), 'Fitted values are correct (product of basis, smoothing and coef).')
	
	#TODO: put this in file test.algorithms.R 
#	n <- 100; m <- 20; r <- 3			
#	set.seed(.TestSeed)
#	V <- syntheticNMF(n, r, m, noise=TRUE) # add noise for stability
#	set.seed(.TestSeed)
#	res.std <- fit(nmf(V, r, 'brunet', seed='nndsvd'))		
#	set.seed(.TestSeed)
#	res.ns <- fit(nmf(V, r, 'ns', seed='nndsvd', theta=0))
	
	#checkEqualsNumeric(basis(res.ns), basis(res.std), 'W from nsNMF(theta=0) are the same as the ones from the standard model')
	#checkEqualsNumeric(coef(res.ns), coef(res.std), 'H from nsNMF(theta=0) are the same as the ones from the standard model')
}
	
test.nmfModel <- function(){
	
	set.seed(.TestSeed)
	
	# define some dimensions to use as template
	n <- as.integer(25); m <- as.integer(10); r <- as.integer(3);
	
	check.empty.model <- function(x, n, m, r, msg, class='NMFstd'){
		check.object(x, n, m, r, msg, class)
		checkOnlyNAs(x, msg)
	}
	
	# check basic errors
	checkException(nmfModel(numeric()), 'Error if negative rank')
	checkException(nmfModel(c(1,2)), 'Error if rank of length != 1')
	checkException(nmfModel(r, -1), 'Error if negative target dimension')
	checkException(nmfModel(r, 1:3), 'Error if target dimension of length > 2')
	checkException(nmfModel(r, 1, -1), 'Error if target ncol negative')
	checkException(nmfModel(r, 1, 1:2), 'Error if target ncol of length > 1')
	checkException(nmfModel(r, 1, matrix(1,2,2)), 'Error if target ncol not vector')
	
	# constructor of empty model 		
		check.empty.model(nmfModel(), 0, 0, 0, 'Constructor with no arguments returns empty model')
	# constructor of empty model with ncol specified 		
		check.empty.model(nmfModel(ncol=5), 0, 5, 0, 'Constructor with only with ncol specified')
	# constructor with missing target
		check.empty.model(nmfModel(r), 0, 0, r, 'Constructor with missing target')
	# constructor using dimension vector 
		check.empty.model(nmfModel(r, c(n,m)), n, m, r, 'Constructor with dimensions as a vector')				
	# constructor using separate dimensions 
		check.empty.model(nmfModel(r, n, m), n, m, r, 'Constructor with separate dimensions')
	# constructor with single dimension
		check.empty.model(nmfModel(r, n), n, n, r, 'Constructor with single dimension (nrow)')
	# constructor with single dimension (ncol)
		check.empty.model(nmfModel(r, ncol=m), 0, m, r, 'Constructor with single dimension (ncol)')
	#TODO: check exceptions if passing negative numbers
		
	# constructor using a target matrix
		msg <- function(...) paste('Constructor with target matrix -', ...)
		V <- rmatrix(n, m) 
		check.empty.model(nmfModel(r, V), n, m, r, msg('second argument')) 
		check.empty.model(nmfModel(V, r), n, m, r, msg('first argument'))
		check.empty.model(nmfModel(V), n, m, 0, msg('single argument'))
		dimnames(V) <- list(rows=letters[1:n], cols=letters[1:m])				
		checkIdentical(dimnames(nmfModel(V)), c(dimnames(V), list(NULL)), msg('dimnames are correctly set'))
		checkIdentical(dimnames(nmfModel(V, use.names=FALSE)), NULL, msg('dimnames are not used if use.names=FALSE'))
		
		
	w <- rmatrix(n, r); h <- rmatrix(r, m)
	# constructor supplying matrices W and H
		msg <- function(...) paste('Constructor with target rank and both basis and coef matrices -', ...)	
		a <- nmfModel(r, W=w, H=h)
		check.object(a, n, m, r, msg())
		checkIdentical(basis(a), w, msg('basis matrix is correctly set'))
		checkIdentical(coef(a), h, msg('coef matrix is correctly set'))	
		checkException(nmfModel(r, c(n+1,m), W=w, H=h), msg('error if bad number of rows in basis matrix'))
		checkException(nmfModel(r, c(n,m), W=w[,-r], H=h), msg('error if bad number of columns in basis matrix'))
		checkException(nmfModel(r, c(n,m), W=w, H=h[-r,]), msg('error if bad number of rows in coef matrix'))
		checkException(nmfModel(r, c(n,m+1), W=w, H=h), msg('error if bad number of columns in coef matrix'))
		
		# reducing dimensions
		rmsg <- function(...) msg('reduce rank -', ...)
		a <- nmfModel(r-1, W=w, H=h)
		check.object(a, n, m, r-1, rmsg('dimensions are OK'))
		checkIdentical(basis(a), w[,-r], rmsg("entries for basis are OK"))
		checkIdentical(coef(a), h[-r,], rmsg("entries for coef are OK"))
		
		rmsg <- function(...) msg('reduce nrow -', ...)
		a <- nmfModel(r, n-1, W=w, H=h)
		check.object(a, n-1, m, r, rmsg('dimensions are OK'))
		checkIdentical(basis(a), w[-n,], rmsg("entries for basis are OK"))
		checkIdentical(coef(a), h, rmsg("entries for coef are OK"))
		
		rmsg <- function(...) msg('reduce ncol -', ...)
		a <- nmfModel(r, ncol=m-1, W=w, H=h)
		check.object(a, n, m-1, r, rmsg('dimensions are OK'))
		checkIdentical(basis(a), w, rmsg("entries for basis are OK"))
		checkIdentical(coef(a), h[,-m], rmsg("entries for coef are OK"))
		
	# constructor supplying only matrices W and H
		msg <- function(...) paste('Constructor with only basis and coef matrices (named arguments) -', ...)	
		a <- nmfModel(W=w, H=h)
		check.object(a, n, m, r, msg())
		checkIdentical(basis(a), w, msg('basis matrix is correctly set'))
		checkIdentical(coef(a), h, msg('coef matrix is correctly set'))
		checkException(nmfModel(W=matrix(1, n, r+1), H=matrix(1, r, m)), msg('error if incompatible dimensions (1)'))
		checkException(nmfModel(W=matrix(1, n, r), H=matrix(1, r+1, m)), msg('error if incompatible dimensions (2)'))		
	# constructor supplying only matrices W and H and not naming the arguments
		msg <- function(...) paste('Constructor with only basis and coef matrices (unamed argument) -', ...)	
		a <- nmfModel(w, h)
		check.object(a, n, m, r, msg())
		checkIdentical(basis(a), w, msg('basis matrix is correctly set'))
		checkIdentical(coef(a), h, msg('coef matrix is correctly set'))
		checkException(nmfModel(matrix(1, n, r+1), matrix(1, r, m)), msg('error if incompatible dimensions (1)'))
		checkException(nmfModel(matrix(1, n, r), matrix(1, r+1, m)), msg('error if incompatible dimensions (2)'))
		
	# constructor supplying only W and target rank
		msg <- function(...) paste('Constructor with target rank and basis matrix only -', ...)
		a <- nmfModel(r, W=w)
		check.object(a, n, 0, r, msg())
		checkIdentical(basis(a), w, msg('basis matrix is correctly set'))
		checkTrue(all(is.na(coef(a))), msg('only NA in coef matrix'))
		checkException(nmfModel(r+1, W=w), msg('error if smaller number of columns in basis matrix'))
		
		# reducing dimensions
		rmsg <- function(...) msg('reduce rank -', ...)
		a <- nmfModel(r-1, W=w)
		check.object(a, n, 0, r-1, rmsg('dimensions are OK'))
		checkIdentical(basis(a), w[,-r], rmsg("entries for basis are OK"))
		checkTrue(all(is.na(coef(a))), rmsg('only NA in coef matrix'))
		checkException(nmfModel(r-1, W=w, force.dim=FALSE), msg('error if greater number of columns in basis matrix and force.dim=FALSE'))
				
	# constructor supplying only W
		msg <- function(...) paste('Constructor with basis matrix only -', ...)
		a <- nmfModel(W=w)
		check.object(a, n, 0, r, msg())
		checkIdentical(basis(a), w, msg('basis matrix is correctly set'))
		checkTrue(all(is.na(coef(a))), msg('only NA in coef matrix'))	
		
	# constructor supplying only H and target rank
		msg <- function(...) paste('Constructor with target rank and coef matrix only -', ...)
		a <- nmfModel(r, H=h)
		check.object(a, 0, m, r, msg())
		checkIdentical(coef(a), h, msg('coef matrix is correctly set'))
		checkTrue(all(is.na(basis(a))), msg('only NA in basis matrix'))
		checkException(nmfModel(r+1, H=h), msg('error if smaller number of rows in coef matrix'))
		
		# reducing dimensions
		rmsg <- function(...) msg('reduce rank -', ...)
		a <- nmfModel(r-1, H=h)
		check.object(a, 0, m, r-1, rmsg('dimensions are OK'))
		checkIdentical(coef(a), h[-r,], rmsg('coef matrix is correctly set'))
		checkTrue(all(is.na(basis(a))), rmsg('only NA in basis matrix'))
		checkException(nmfModel(r-1, H=h, force.dim=FALSE), msg('error if greater number of rows in coef matrix and force.dim=FALSE'))
		
	# constructor supplying only H
		msg <- function(...) paste('Constructor with coef matrix only -', ...)
		a <- nmfModel(H=h)
		check.object(a, 0, m, r, msg())
		checkIdentical(coef(a), h, msg('coef matrix is correctly set'))
		checkTrue(all(is.na(basis(a))), msg('only NA in basis matrix'))
						
	# constructor supplying W and target dimensions
		msg <- function(...) paste('Constructor with basis matrix and both target dimensions -', ...)
		a <- nmfModel(r, n, m, W=w)
		check.object(a, n, m, r, msg())
		checkIdentical(basis(a), w, rmsg("entries for basis are OK"))
		checkTrue(all(is.na(coef(a))), rmsg('only NA in coef matrix'))
		
		rmsg <- function(...) msg('reduce nrow -', ...)
		a <- nmfModel(r, n-1, m, W=w)
		check.object(a, n-1, m, r, rmsg('dimensions are OK'))
		checkIdentical(basis(a), w[-n,], rmsg("entries for basis are OK"))
		checkTrue(all(is.na(coef(a))), rmsg('only NA in coef matrix'))		
		
		checkException(nmfModel(r, n+1, m, W=w), msg('error if smaller number of rows in basis matrix'))
		checkException(nmfModel(r, n-1, W=w, force.dim=FALSE), msg('error if greater number of rows in basis matrix and force.dim=FALSE'))
		checkException(nmfModel(r+1, n, m, W=w), msg('error if smaller number of columns in basis matrix'))		
		
	# constructor supplying H and target dimensions		
		msg <- function(...) paste('Constructor with coef matrix and both target dimensions -', ...)
		a <- nmfModel(r, n, m, H=h)
		check.object(a, n, m, r, msg())
		checkTrue(all(is.na(basis(a))), msg('only NA in basis matrix'))
		checkIdentical(coef(a), h, msg('coef matrix is correctly set'))
		
		rmsg <- function(...) msg('reduce ncol -', ...)
		a <- nmfModel(r, n , m-1, H=h)
		check.object(a, n, m-1, r, rmsg('dimensions are OK'))
		checkIdentical(coef(a), h[,-m], rmsg('coef matrix is correctly set'))
		checkTrue(all(is.na(basis(a))), rmsg('only NA in basis matrix'))
		
		checkException(nmfModel(r+1, n, m, H=h), msg('error if smaller number of rows in coef matrix'))
		checkException(nmfModel(r, n, m-1, H=h, force.dim=FALSE), msg('error if greater number of columns in coef matrix and force.dim=FALSE'))
		checkException(nmfModel(r, n, m+1, H=h), msg('error if smaller number of columns in coef matrix'))
	
	# check basisnames passage		

		check.model.dimnames <- function(x, dn, title){
			msg <- function(...) paste(title, '-', ...)
			checkIdentical(dimnames(x), dn, msg("dimnames are correct"))
			checkIdentical(colnames(basis(x)), dn[[3]], msg("colnames of basis matrix are correct"))
			checkIdentical(rownames(coef(x)), dn[[3]], msg("rownames of coef matrix are correct"))
		}

		dn <- letters[1:nrow(h)]
		h2 <- h; rownames(h2) <- dn; 
		a <- nmfModel(r, H=h2)
		check.model.dimnames(a, list(NULL, NULL, dn), "Basis names are passed from input coef matrix")
		
		w2 <- w; colnames(w2) <- dn; 
		a <- nmfModel(r, W=w2)
		check.model.dimnames(a, list(NULL, NULL, dn), "Basis names are passed from input basis matrix")
		
		w2 <- w; colnames(w2) <- dn;
		h2 <- h; rownames(h2) <- dn; 
		a <- nmfModel(W=w2, H=h2)
		check.model.dimnames(a, list(NULL, NULL, dn), "Basis names are used unchanged if equal in input basis and coef matrices")
		
		msg <- function(...) paste("Basis names from input basis matrix are used to order the components - ", ...)	
		w2 <- w; colnames(w2) <- dn;
		h2 <- h; rownames(h2) <- rev(dn);
		a <- nmfModel(W=w2, H=h2)
		check.model.dimnames(a, list(NULL, NULL, dn), msg("rownames of input basis are enforced"))
		checkIdentical(basis(a), w2, msg("basis unchanged"))
		checkIdentical(coef(a), h2[nrow(h):1,], msg("coef entries are reordered"))
		
		msg <- function(...) paste("Basis names from input basis matrix are NOT used to order the components if argument order.basis=FALSE - ", ...)
		w2 <- w; colnames(w2) <- dn;
		h2 <- h; rownames(h2) <- rev(dn);
		a <- nmfModel(W=w2, H=h2, order.basis=FALSE)
		check.model.dimnames(a, list(NULL, NULL, dn), msg("rownames of input basis are enforced"))
		checkIdentical(basis(a), w2, msg("basis unchanged"))
		checkEquals(coef(a), h2, msg("coef entries are not ordered"), check.attributes=FALSE)
		
		msg <- function(...) paste("Basis names from input basis matrix are forced onto to the coef matrix - ", ...)	
		w2 <- w; colnames(w2) <- dn;
		h2 <- h; rownames(h2) <- paste(letters[1:nrow(h)], 2);
		a <- nmfModel(W=w2, H=h2)
		check.model.dimnames(a, list(NULL, NULL, dn), msg("rownames of input basis are enforced"))
		checkIdentical(basis(a), w2, msg("basis is unchanged"))
		checkEquals(coef(a), h2, msg("coef entries are correct"), check.attributes=FALSE)
}

test.dimensions <- function(){
	# define some dimensions to use as template
	n <- as.integer(50); m <- as.integer(10); r <- as.integer(3);
	
	# create an NMF object
	a <- nmfModel(r, n, m)
	
	# check the dimensions	
	checkIdentical(dim(a), c(n, m, r), "Function 'dim' is OK")
	
	# check the row number
	checkIdentical(nrow(a), n, "Function 'nrow' is OK")
	
	# check the column number
	checkIdentical(ncol(a), m, "Function 'ncol' is OK")
	
	# check the number of basis
	checkIdentical(nbasis(a), r, "Function 'nbasis' is OK")

}

check.dimnames <- function(x, dn, msg){
	
	checkIdentical(dimnames(x), dn, paste(msg, '-',"dimnames returns correct value"))
	checkIdentical(dimnames(x)[c(1,3)], dimnames(basis(x)), paste(msg, '-', "dimnames returns value consistent with basis"))
	checkIdentical(dimnames(x)[c(3,2)], dimnames(coef(x)), paste(msg, '-', "dimnames returns value consistent with coef"))
	
	checkIdentical(rownames(x), dn[[1]], paste(msg, '-', "rownames returns correct value"))
	checkIdentical(rownames(basis(x)), rownames(x), paste(msg, '-', "rownames returns same value as rownames(basis)"))
	
	checkIdentical(colnames(x), dn[[2]], paste(msg, '-', "colnames returns correct value"))
	checkIdentical(colnames(coef(x)), colnames(x), paste(msg, '-', "colnames returns same value as colnames(basis)"))
	
	checkIdentical(basisnames(x), dn[[3]], paste(msg, '-', "basisnames returns correct value"))
	checkIdentical(colnames(basis(x)), basisnames(x), paste(msg, '-', "basisnames returns same value as colnames(basis)"))
	checkIdentical(rownames(coef(x)), basisnames(x), paste(msg, '-', "basisnames returns same value as rownames(coef)"))
	
}

test.dimnames <- function(){
	
	# set random seed
	set.seed(.TestSeed)	
	# define some dimensions to use as template
	n <- 20; m <- 10; r <- 3;
	# generate matrices
	w <- rmatrix(n, r); h <- rmatrix(r, m)
	
	M <- nmfModel(r, n, m)
	
	# check errors	
	a <- M
	checkException({dimnames(a) <- 1:n}, 'set to vector')
	checkException({dimnames(a) <- list(seq(n-1))}, 'Error of wrong dimension (nrow-1)')
	checkException({dimnames(a) <- list(seq(n+1))}, 'Error of wrong dimension (nrow+1)')
	
	# check with no elements
	msg <- function(...) paste('Dimnames with 2 elements -', ...)
	a <- M
	check.dimnames(a, NULL, 'No dimnames => NULL')
	check.dimnames({dimnames(a) <- NULL; a}, NULL, msg('set to NULL'))
	check.dimnames({dimnames(a) <- list(); a}, NULL, msg('set to list()'))	
		
	# check with one elements
	msg <- function(...) paste('Dimnames with 1 element -', ...)
	a <- M
	dn <- list(letters[1:nrow(a)])
	dn.name <- setNames(dn, 'rows')	
	check.dimnames({a <- M; dimnames(a) <- dn; a}, c(dn, list(NULL, NULL)), msg('Set dimnames'))
	check.dimnames({a <- M; dimnames(a) <- dn.name; a}, c(dn.name, list(NULL, NULL)), msg('Set with names'))
		
	check.dimnames({a <- M; rownames(a) <- dn[[1]]; a}, c(dn, list(NULL, NULL)), msg('Set rownames'))
	check.dimnames({a <- M; colnames(a) <- letters[1:ncol(a)]; a}, list(NULL, letters[1:ncol(a)], NULL), msg('Set colnames'))
	check.dimnames({a <- M; basisnames(a) <- letters[1:nbasis(a)]; a}, list(NULL, NULL, letters[1:nbasis(a)]), msg('Set basisnames'))	
	check.dimnames({dimnames(a) <- NULL; a}, NULL, msg('Reset to NULL'))
		
	# check with two elements
	msg <- function(...) paste('Dimnames with 2 elements -', ...)
	a <- M
	dn <- list(letters[1:nrow(a)], letters[seq(nrow(a)+1, nrow(a)+ncol(a))])
	dn.name <- setNames(dn, c('rows', 'cols'))	
	check.dimnames({dimnames(a) <- dn; a}, c(dn, list(NULL)), msg('Set dimnames'))
	check.dimnames({dimnames(a) <- NULL; a}, NULL, msg('Reset to NULL'))
	check.dimnames({dimnames(a) <- dn.name; a}, c(dn.name, list(NULL)), msg('Set with names'))
	check.dimnames({dimnames(a) <- NULL; a}, NULL, msg('Reset to NULL (2)'))
		
	# check with three elements
	msg <- function(...) paste('Dimnames with 3 elements -', ...)
	a <- M
	dn <- list(letters[1:nrow(a)], letters[seq(nrow(a)+1, nrow(a)+ncol(a))], letters[seq(nrow(a)+ncol(a)+1, nrow(a)+ncol(a)+nbasis(a))])
	dn.name <- setNames(dn, c('rows', 'cols', 'basis'))	
	check.dimnames({dimnames(a) <- dn; a}, dn, msg('Set dimnames'))	
	check.dimnames({dimnames(a) <- NULL; a}, NULL, msg('Reset to NULL'))
	check.dimnames({dimnames(a) <- dn.name; a}, dn.name, msg('Set with names'))
	check.dimnames({dimnames(a) <- NULL; a}, NULL, msg('Reset to NULL (2)'))	
}

#' Unit test for accessor method: \code{basis}
test.basis <- function(){
	
	# define some dimensions to use as template
	n <- 50; m <- 10; r <- 3;	
	
	# create an empty object
	a <- nmfModel(r);
	
	# get factor W
	checkTrue(is.same(basis(a), a@W), "Method 'basis' correctly returns slot W");
	
	# set factor W
	W.ext <- matrix(1, n, r)
	basis(a) <- W.ext
	checkIdentical(a@W, W.ext, "Method 'basis<-' correctly sets slot W");
	#checkException(basis(a) <- matrix(1, n, r+1), "Error if setting basis with a matrix of incompatible dimensions");
}

#' Unit test for accessor method: \code{coef}
test.coef <- function(){
	
	# define some dimensions to use as template
	n <- 50; m <- 10; r <- 3;	
	
	# create an empty object
	a <- nmfModel(r);
	
	# get factor H
	checkTrue(is.same(coef(a), a@H), "Method 'coef' correctly returns slot H");
	
	# set factor H
	ext <- matrix(1, r, m)
	coef(a) <- ext
	checkIdentical(a@H, ext, "Method 'coef<-' correctly sets slot H");
	#checkException(coef(a) <- matrix(1, r+1, m), "Error if setting coef with a matrix of incompatible dimensions");
}

test.NMF.rnmf <- function(){
	
	# define some dimensions to use as template
	n <- 50; m <- 10; r <- 3;	
	
	# create an empty NMF object
	a <- nmfModel(r)
	
	# check parameter validity checks
	set.seed(.TestSeed)
	checkException(rnmf(a, as.numeric(NA)), 'Error thrown if single NA value')
	checkException(rnmf(a, as.numeric(c(1,NA))), 'Error thrown if some NA value')
	checkException(rnmf(a, c(1,2,3)), 'Error thrown if some length greater than 2')
	checkException(rnmf(a, numeric()), 'Error thrown if some length is 0')
	
	# create a random NMF of given dimension
	set.seed(.TestSeed)
	a <- rnmf(r, n, m)
	checkPlot(basicHM(basis(a)), 'Random NMF basis (target dimension)')
	checkPlot(basicHM(coef(a)), 'Random NMF coef (target dimension)')
	# check the NMF dimensions	
	checkEquals(nrow(a), n, 'number of rows matches target dimension')
	checkEquals(ncol(a), m, 'number of columns matches target dimension')
	checkEquals(nbasis(a), r, 'number of basis matches target dimension')
	# check equivalence of calls
	set.seed(.TestSeed)
	checkIdentical(a, rnmf(r, c(n, m)), "calling with numeric target of length 2 is equivalent to separate dimensions")
	set.seed(.TestSeed)
	checkIdentical(a, rnmf(r, n, m), "calling with numeric target of length 2 is equivalent to separate dimensions")	
	
	# create a random NMF based on a template NMF	
	a <- nmfModel(r, n, m)
	set.seed(.TestSeed)
	b <- rnmf(a)
	checkPlot(basicHM(basis(b)), 'Random NMF basis (target NMF)')
	checkPlot(basicHM(coef(b)), 'Random NMF coef (target NMF)')
	# check the NMF dimensions	
	checkEquals(nrow(a), nrow(b), 'number of rows matches NMF target')
	checkEquals(ncol(a), ncol(b), 'number of columns matches NMF target')
	checkEquals(nbasis(a), nbasis(b), 'number of basis matches NMF model')	
		
	# create a random NMF based on a target ExpressionSet
	set.seed(.TestSeed)
	max.entry <- 100 
	V <- new('ExpressionSet', exprs=rmatrix(n, m, min=0, max=max.entry))	
	a <- nmfModel(r, n, m)
	check.object(a <- rnmf(a, V), n, m, r, 'NMFobject + target ExpressionSet')
	check.dimnames(a, c(dimnames(exprs(V)), list(NULL)), 'NMFobject + target ExpressionSet')	
	check.object(a <- rnmf(r, V), n, m, r, 'rank + target ExpressionSet')
	check.dimnames(a, c(dimnames(exprs(V)), list(NULL)), 'rank + target ExpressionSet')
	
	# create a random NMF based on a target matrix
	set.seed(.TestSeed)	
	max.entry <- 100 
	V <- rmatrix(n, m, min=0, max=max.entry)	
	a <- nmfModel(r, n, m)	
	a <- rnmf(a, V)
	checkPlot(basicHM(basis(a)), 'Random NMF basis (target matrix)')
	checkPlot(basicHM(coef(a)), 'Random NMF coef (target matrix)')
	# check the NMF dimensions	
	checkEquals(nrow(a), n, 'number of rows matches matrix target')
	checkEquals(ncol(a), m, 'number of columns matches matrix target')
	checkEquals(nbasis(a), r, 'matrix target: number of basis matches NMF model')
	# check maximum value
	msg <- function(...) paste('Set max in target matrix -', ...)
	.check_max <- function(x, max.entry){
		checkTrue( max(basis(a)) <= max.entry, msg("Basis maximum entry is OK (<=)"))
		checkTrue( max(basis(a)) >= max.entry/2, msg("Basis maximum entry is OK (>=)"))
		checkTrue( max(coef(a)) <= max.entry, msg("Coef maximum entry is OK (<=)"))
		checkTrue( max(coef(a)) >= max.entry/2, msg("Coef maximum entry is OK (>=)"))
	}
	.check_max(a, max.entry)
	# check dimnames
	dn <- list(rows=letters[1:nrow(V)], cols=letters[seq(25-ncol(V)+1,25)])
	dimnames(V) <- dn	
	check.dimnames(rnmf(a, V), c(dn, list(NULL)), "rnmf on target matrix with dimnames: dimnames are passed")
	check.dimnames(rnmf(a, V, use.dimnames=FALSE), NULL, "rnmf on target matrix with dimnames and use.dimnames=FALSE: dimnames are not passed")
	
	# setting maximum on NMFOffset model
	msg <- function(...) paste('Set max by argument -', ...)
	set.seed(.TestSeed)
	a <- rnmf(r, V, model='NMFOffset')
	.check_max(a, max.entry)
	
	# setting maximum entry by argument
	msg <- function(...) paste('Set max by argument -', ...)
	set.seed(.TestSeed)
	max.entry <- 5
	a <- rnmf(r, n, m, dist=list(max=max.entry))
	.check_max(a, max.entry)
	
	
}

#' Tests synthetic data generation
test.syntheticNMF <- function(){
	
	# define sizes	
	n <- 100; m <- 20; r <- 3
	
	# standard
	set.seed(.TestSeed)
	d <- syntheticNMF(n, r, m)	
	checkPlot(basicHM(d), 'Synthetic data plain')
	
	# with offset
	set.seed(.TestSeed)
	n.offset <- 15
	o <- c(rep(1, n.offset), rep(0, n-n.offset))
	d <- syntheticNMF(n, r, m, offset=o)	
	# the function should throw an error when called with an offset of wrong length
	checkException(syntheticNMF(n, r, m, offset=o[-1]))
	checkPlot(basicHM(d), 'Synthetic data offset')
	
	# with noise
	set.seed(.TestSeed)
	d <- syntheticNMF(n, r, m, noise=TRUE)
	checkPlot(basicHM(d), 'Synthetic data with noise')
	
	# with offset and noise
	set.seed(.TestSeed)
	d <- syntheticNMF(n, r, m, offset=o, noise=TRUE)
	checkPlot(basicHM(d), 'Synthetic data with offset and noise')
}

# Sparseness
test.sparseness <- function(){
	
	# local function to check bounds of sparseness
	checkBounds <- function(s){
		checkTrue(s >= 0, 'greater than 0')
		checkTrue(s <= 1, 'lower than 1')	 
	}
	
	# test with a perfectly non-sparse vector (constant vector)
		x <- rep(1, 100)
		s <- sparseness(x)
		
		# check values
		checkBounds(s)
		checkEqualsNumeric(s, 0, 'should be 0')
		
	# test with a perfectly sparse vector
		x <- c(1, rep(0, 100))
		s <- sparseness(x)
		
		# check values
		checkBounds(s)
		checkEqualsNumeric(s, 1)
		
	# define sizes	
	n <- 100; m <- 20; r <- 3
		
	# test with a random matrix (not necessarly positive)
		set.seed(.TestSeed)		
		V <- matrix(rnorm(n*m), n, m)
		s <- sparseness(V)
		#check values
		checkBounds(s)
		
	# test with a random NMF object
		a <- nmfModel(r, V)
		set.seed(.TestSeed)
		a <- rnmf(a, V)
		s <- sparseness(a)
		# check return vector
		checkTrue( length(s) == 2, "Method 'sparseness' returns a 2-length vector" )
		#check values
		checkBounds(s[1])
		checkBounds(s[2])
}

# Clusters
test.predict <- function(){

	# define some dimensions 
	n <- 100; m <- 20; r <- 3

	.msg <- NULL
	mess <- function(...) paste(.msg, ':', ...)
	
	# test on a known matrix
		.msg <- 'Artificial matrix'
		V <- matrix(c( c(1, rep(0,n-1)), c(0, 1, rep(0,n-2)), rep(c(0,0,1, rep(0, n-3)), m-2) ) , n, m)
		# compute clusters
		res <- .predict.nmf(t(V))
		# check the result
		checkEquals(res, as.factor(c(1, 2, rep(3, m-2))), mess('Return known clusters'))

	# test on a random matrix
		.msg <- 'Random matrix'
		set.seed(.TestSeed)		
		V <- matrix(sapply(sample(n), function(i){x <- rep(0, n); x[i] <- 1; x}), n, m)
			
		# compute clusters
		res <- .predict.nmf(V)
		# check the result
		checkTrue( is.factor(res), mess('Result is a factor'))
		checkTrue( length(res) == nrow(V), mess('Result is the right size' ))
		checkTrue( nlevels(res) == ncol(V), mess('Result has the right number of levels' ))		
	
	# test on a random NMF	
		.msg <- 'NMF model'
		a <- nmfModel(r, V)
		set.seed(.TestSeed)
		a <- rnmf(a)
	
		res <- predict(a, 'samples')
		# check the result
		checkTrue( is.factor(res), mess('Result is a factor'))
		checkTrue( length(res) == ncol(a), mess('Result has right size' ))
		checkTrue( nlevels(res) == nbasis(a), mess('Result has right number of levels' ))
	
	# factor
	.msg <- 'Factor'
		res <- predict(a, 'features')
		# check the result
		checkTrue( is.factor(res), mess('Result is a factor'))
		checkTrue( length(res) == nrow(a), mess('Result has right size' ))
		checkTrue( nlevels(res) == nbasis(a), mess('Result has right number of levels' ))

}

# Purity
test.purity <- function(){

	# local function to check bounds of purity
	checkBounds <- function(x){
		checkTrue(x >= 0, 'greater than 0')
		checkTrue(x <= 1, 'lower than 1')	 
	}
	
	# generate factor
	x <- as.factor(c(rep(1,5), rep(2,10), rep(3,15)))
	
	# compute perfect clustering
	p <- purity(x, x)
	checkBounds(p)
	checkEqualsNumeric(p, 1)
	
	# compute random clustering
	set.seed(.TestSeed)
	p <- purity(as.factor(sample(x)), x)
	checkBounds(p)
}

# Entropy
test.entropy <- function(){
	
	# local function to check bounds of entropy
	checkBounds <- function(x){
		checkTrue(x >= 0, 'greater than 0')
		checkTrue(x <= 1, 'lower than 1')	 
	}
	
	# generate factor
	x <- as.factor(c(rep(1,5), rep(2,10), rep(3,15)))
	
	# comppute perfect entropy
	e <- entropy(x, x)
	checkBounds(e)
	checkEqualsNumeric(e, 0)
	
	# compute random clustering
	set.seed(.TestSeed)
	e <- entropy(as.factor(sample(x)), x)
	checkBounds(e)
	
}

#' Tests the distance method
test.deviance <- function(){

	# define some dimensions and matrices
	n <- 10; m <- 5;
	set.seed(.TestSeed)
	y <- rmatrix(n, m) # random matrix 1	
	x <- rnmf(nmfModel(3, y))	
		
	# check default is no computation
	checkEquals( deviance(x, y), as.numeric(NA), 'With no method: NA')
	# check error if undefined method	
	checkException( deviance(x, y, 'toto'), 'Error if undefined method')
	
	# check if works when passing a function name
		norm1 <- function(x, y){ sum(abs(fitted(x)-y)) } # local function
		# generate a random function name
		funname <- paste('test', paste(sample(1:9, 20, replace=TRUE), collapse=''), sep='.')
		assign(funname, norm1, envir=.GlobalEnv)
		on.exit(rm(list=funname, envir=.GlobalEnv), add=TRUE) 
		#TODO: make it work without assigning in .GloablEnv
		checkTrue( deviance(x, y, norm1) > 0, 'Works with a user-defined function in .GlobalEnv')
		checkTrue( deviance(x, y, funname) > 0, 'Works with name of a user-defined function .GlobalEnv')	
	
	# check euclidean distance
		meth <- 'euclidean'
		#checkEqualsNumeric( deviance(x, x, meth), 0, 'Euclidean: separation')
		checkTrue( deviance(x, y, meth) > 0, 'Euclidean: positive')
		checkEqualsNumeric( deviance(x, y, meth), sum((fitted(x) - y)^2)/2, 'Euclidean: OK')
	
	# check Kullback-Leibler divergence
		meth <- 'KL'
		#checkEqualsNumeric( deviance(x, x, meth), 0, 'KL: separation')
		checkTrue( deviance(x, y, meth) > 0, 'KL: positive')
		#checkTrue( deviance(x, y, meth) != deviance(y, x, meth), 'KL: not symetric')
		# check if not infinite when there is some zero entries in the target matrix
		z <- y; z[1,] <- 0
		checkTrue( deviance(x, z, meth) != Inf, 'Ok if some zeros in the target' )
		# check if infinite when there is some zero entries in the second term
		z <- x; basis(z)[1,] <- 0
		checkIdentical( deviance(z, y, meth), Inf, 'Infinite if some zeros in the estimate' )	
}


#' Tests the connectivity method
test.connectivity <- function(){
	# define some dimensions to use as template
	n <- 50; m <- 10; r <- 3;		
	
	# build random NMF
	set.seed(.TestSeed)
	a <- nmfModel(r, c(n,m))
	a <- rnmf(a)
	con <- connectivity(a)
	
	# check properties of the connectivity matrix
	checkTrue( is.matrix(con), 'The result is a matrix')
	checkTrue( all(con %in% c(0,1)), 'All entries are 0 or 1')
	checkTrue( all(t(con) == con) , 'The connectivity matrix is symmetric')	
	
}

#' test subsetting function
test.subset <- function(){
	
	# create a random NMF object
	n <- 30; r <- 5; p <- 20
	a <- nmfModel(r, n, p)
	a <- rnmf(a)
	
	# fake subset
	checkTrue(identical(a[], a), "subset [] is OK (identical)")
	checkTrue(identical(a[,], a), "subset [,] is OK (identical)")
	checkTrue(identical(a[TRUE,], a), "subset [TRUE,] is OK (identical)")
	checkTrue(identical(a[,TRUE], a), "subset [,TRUE] is OK (identical)")
	checkTrue(identical(a[TRUE,TRUE], a), "subset [TRUE,TRUE] is OK (identical)")
	checkTrue(identical(a[,,TRUE], a), "subset [,,TRUE] is OK (identical)")
	checkTrue(identical(a[TRUE,TRUE,TRUE], a), "subset [TRUE,TRUE,TRUE] is OK (identical)")
	# with NULL
	checkTrue(identical(a[NULL,], a[0,]), "subset [NULL,] is OK")
	checkTrue(identical(a[,NULL], a[,0]), "subset [,NULL] is OK")
	checkTrue(identical(a[NULL,NULL], a[0,0]), "subset [NULL,NULL] is OK")
	checkTrue(identical(a[NULL,NULL,NULL], a[0,0,0]), "subset [NULL,NULL,NULL] is OK")
	checkException(a[,,], "Error when called with [,,]")
		
	# subset on features
	checkEquals(dim(a[1,]), c(1, p, r), "subset 1 feature is OK (dim)")
	checkEquals(basis(a[5,]), basis(a)[5,, drop=FALSE], "subset 1 feature is OK (basis)")
	checkEquals(coef(a[5,]), coef(a), "subset 1 feature is OK (coef)")
	checkEquals(a[5,,drop=TRUE], basis(a)[5,, drop=TRUE], "subset 1 feature dropping is OK (return basis)")
	
	checkEquals(dim(a[1:10,]), c(10, p, r), "subset more than 1 feature is OK (dim)")
	checkEquals(basis(a[1:10,]), basis(a)[1:10,, drop=FALSE], "subset more than 1 feature is OK (basis)")	
	checkEquals(coef(a[1:10,]), coef(a), "subset more than 1 feature is OK (coef)")
	checkEquals(a[1:10,,drop=TRUE], basis(a)[1:10,, drop=TRUE], "subset more than 1 feature dropping is OK (return basis)")
	
	# subset on samples
	checkEquals(dim(a[,1]), c(n, 1, r), "subset 1 sample is OK (dim)")
	checkEquals(coef(a[,5]), coef(a)[,5, drop=FALSE], "subset 1 sample is OK (coef)")	
	checkEquals(basis(a[,1]), basis(a), "subset 1 sample is OK (basis)")
	checkEquals(a[,5,drop=TRUE], coef(a)[,5, drop=TRUE], "subset 1 sample dropping is OK (return coef)")
	
	checkEquals(dim(a[,1:10]), c(n, 10, r), "subset more then 1 sample is OK (dim)")
	checkEquals(coef(a[,1:10]), coef(a)[,1:10, drop=FALSE], "subset more than 1 sample is OK (coef)")	
	checkEquals(basis(a[,1:10]), basis(a), "subset more than 1 sample is OK (basis)")	
	checkEquals(a[,1:10,drop=TRUE], coef(a)[,1:10, drop=TRUE], "subset more than 1 sample dropping is OK (return coef)")
	
	# subset on basis	
	checkEquals(dim(a[,,1]), c(n, p, 1), "subset 1 basis is OK (dim)")
	checkEquals(coef(a[,,3]), coef(a)[3,, drop=FALSE], "subset 1 basis is OK (coef)")	
	checkEquals(basis(a[,,3]), basis(a)[,3, drop=FALSE], "subset 1 basis is OK (basis)")
	checkTrue(identical(a[,,3,drop=TRUE], a[,,3,drop=TRUE]), "subset 1 basis dropping is OK (do nothing)")
	
	checkEquals(dim(a[,,2:4]), c(n, p, 3), "subset more than 1 basis is OK (dim)")
	checkEquals(coef(a[,,2:4]), coef(a)[2:4, , drop=FALSE], "subset more than 1 basis is OK (coef)")	
	checkEquals(basis(a[,,2:4]), basis(a)[, 2:4, drop=FALSE], "subset more than 1 basis is OK (basis)")
	checkTrue(identical(a[,,2:4,drop=TRUE], a[,,2:4,drop=FALSE]), "subset more than 1 basis dropping is OK (do nothing)")
	
	checkEquals(dim(a[NULL]), c(n, p, 0), "subset basis NULL is OK (dim)")
    checkIdentical(a[2], basis(a)[,2], "subset with single index + drop missing returns single basis as vector")
	checkIdentical(a[2, drop = FALSE], a[,,2], "subset with single index with drop=FALSE returns the complete NMF object")
	checkIdentical(a[2, drop=TRUE], basis(a)[,2, drop=TRUE], "subset with single index + drop=TRUE returns single basis as vector")
	
	checkTrue(is.nmf(a[2:3]), "subset with single vector index returns NMF object")
	checkEquals(basis(a[2:3]), basis(a)[,2:3], "subset with single vector index returns subset of NMF object")
	checkEquals(a[2:3, drop=TRUE], basis(a)[,2:3, drop=TRUE], "subset with single vector index + dropping returns correct matrix if length > 1")
	checkIdentical(a[2:3, drop=FALSE], a[,,2:3], "subset with single vector index + NOT dropping returns correct matrix if length > 1")
	
	# subset on both features and samples
	checkEquals(dim(a[1,1]), c(1, 1, r), "subset 1 feature x 1 sample is OK (dim)")
	checkEquals(coef(a[3,5]), coef(a)[,5, drop=FALSE], "subset 1 feature x 1 sample is OK (coef)")
	checkEquals(basis(a[3,5]), basis(a)[3,, drop=FALSE], "subset 1 feature x 1 sample is OK (basis)")
	
	checkEquals(dim(a[10:19,5:11]), c(10, 7, r), "subset more than 1 feature x sample is OK (dim)")
	checkEquals(coef(a[10:19,5:11]), coef(a)[,5:11], "subset more than 1 feature x sample is OK (coef)")
	checkEquals(basis(a[10:19,5:11]), basis(a)[10:19,], "subset more than 1 feature x sample is OK (basis)")
}


#' Tests for get/set misc elements in NMF models
test.misc <- function(){
	
	x <- nmfModel()
	m <- slot(x, 'misc')
	checkTrue(is.list(m) && length(m)==0L, 'On empty model misc is an empty list')
	checkEquals(misc(x), m, 'On empty model misc() returns an empty list')
	x$a <- 3
	checkEquals(slot(x, 'misc'), list(a=3), "Setting misc with $ works")
	checkEquals(misc(x), list(a=3), "Getting misc with misc() works")
	checkEquals(x$a, 3, "Getting misc with $ works")
	
	checkEquals(misc(list()), NULL, 'On empty list misc is NULL')
	checkEquals(misc(list(misc=4)), 4, 'On list with a `misc` element, misc() returns the element')
	checkEquals(misc(1), NULL, 'On non list object misc() is NULL')
	a <- 1
	attr(a, 'misc') <- 2
	checkEquals(misc(a), 2, 'On non list object with `misc` attribute, misc() returns the attribute')
	
}