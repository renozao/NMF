#' Unit Testing script for NMF package: NMF distances.
#'
#' @author Renaud Gaujoux
#' @creation 08 April 2010

.TestSeed <- 123456

# make the internal function visible
if( isNamespaceLoaded('NMF') ){
	.rss <- NMF:::.rss
	.KL <- NMF:::.KL	
}



check.storage.mode <- function(R.fun, C.fun, name, ...){
		
	n <- 5000; p <- 50;
	
	# create random matrices in storage mode double 
	x <- rmatrix(n, p, min=1, max=100)
	y <- rmatrix(n, p, min=1, max=100)
	
	# create random matrices in storage mode integer 
	xi <- x; storage.mode(xi) <- 'integer'
	yi <- y; storage.mode(yi) <- 'integer'
	
	# check result for all combinations
	checkEquals( R.fun(x, y, ...), C.fun(x, y), paste(name, "- Version double-double: OK" ))
	checkEquals( R.fun(x, yi), C.fun(x, yi), paste(name, "- Version double-integer: OK" ))
	checkEquals( R.fun(xi, y), C.fun(xi, y), paste(name, "- Version integer-double: OK" ))
	checkEquals( R.fun(xi, yi), C.fun(xi, yi), paste(name, "- Version integer-integer: OK" ))
	
}

test.rss <- function(){
	
	set.seed(.TestSeed)	
	
	# create R version for RSS
	R_rss <- function(x, y) sum((x-y)^2)	
		
	# check the storage mode
	check.storage.mode(R_rss, .rss, 'RSS')	
}

test.KL <- function(){
	
	set.seed(.TestSeed)	
	
	# create R version for RSS
	R_kl <- function(x, y) sum( ifelse(x==0, y, x * log(x/y) - x + y) );
		
	# check the storage mode
	check.storage.mode(R_kl, .KL, 'KL divergence')
			
	# create random matrices
	n <- 5000; p <- 50;
	x <- rmatrix(n, p, min=1, max=100)
	y <- rmatrix(n, p, min=1, max=100)
	
	# check result with NA values
	z <- x; z[1,1] <- NA 
	checkTrue( is.na(.KL(z, y)), "Return NA if there is a NA in x")
	checkTrue( is.na(.KL(x, z)), "Return NA if there is a NA in y")
	
	# check result with NaN values
	z <- x; z[1,1] <- NaN 
	checkTrue( is.na(.KL(z, y)), "Return NA if there is a NaN in x")
	checkTrue( is.na(.KL(x, z)), "Return NA if there is a NaN in y")
	
	# check result with 0 values in y
	z <- y; z[1,1] <- 0 
	checkEquals( .KL(x, z), Inf, "Return Inf if there is a 0 in y")	
	
}