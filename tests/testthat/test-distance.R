#' Unit Testing script for NMF package: NMF distances.
#'
#' @author Renaud Gaujoux
#' @creation 08 April 2010
# Converted from RUnit: 16 Feb 2020

.TestSeed <- 123456

check.storage.mode <- function(R.fun, C.fun, name, ...){
		
	n <- 5000; p <- 50;
	
	# create random matrices in storage mode double 
	x <- rmatrix(n, p, min=1, max=100)
	y <- rmatrix(n, p, min=1, max=100)
	
	# create random matrices in storage mode integer 
	xi <- x; storage.mode(xi) <- 'integer'
	yi <- y; storage.mode(yi) <- 'integer'
	
	# check result for all combinations
	expect_equal( R.fun(x, y, ...), C.fun(x, y), info = paste(name, "- Version double-double: OK" ))
	expect_equal( R.fun(x, yi), C.fun(x, yi), info = paste(name, "- Version double-integer: OK" ))
	expect_equal( R.fun(xi, y), C.fun(xi, y), info = paste(name, "- Version integer-double: OK" ))
	expect_equal( R.fun(xi, yi), C.fun(xi, yi), info = paste(name, "- Version integer-integer: OK" ))
	
}

test_that("test.KL", {
    set.seed(.TestSeed)
    R_kl <- function(x, y) sum(ifelse(x == 0, y, x * log(x/y) - 
        x + y))
    check.storage.mode(R_kl, .KL, "KL divergence")
    n <- 5000
    p <- 50
    x <- rmatrix(n, p, min = 1, max = 100)
    y <- rmatrix(n, p, min = 1, max = 100)
    z <- x
    z[1, 1] <- NA
    expect_true(is.na(.KL(z, y)), info = "Return NA if there is a NA in x")
    expect_true(is.na(.KL(x, z)), info = "Return NA if there is a NA in y")
    z <- x
    z[1, 1] <- NaN
    expect_true(is.na(.KL(z, y)), info = "Return NA if there is a NaN in x")
    expect_true(is.na(.KL(x, z)), info = "Return NA if there is a NaN in y")
    z <- y
    z[1, 1] <- 0
    expect_equal(Inf, .KL(x, z), info = "Return Inf if there is a 0 in y")
})

test_that("test.rss", {
    set.seed(.TestSeed)
    R_rss <- function(x, y) sum((x - y)^2)
    check.storage.mode(R_rss, .rss, "RSS")
})

